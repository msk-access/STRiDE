import pickle
from pathlib import Path
from typing import Any, Optional, Union

import joblib
import numpy as np
import pandas as pd

from stride.core.base import BasePredictor
from stride.core.dataset import process_single_sample
from stride.utils.torch_patches import apply_cpu_patches, force_cpu, setup_tabpfn_shims

DEFAULT_TABPFN_DIR = Path(__file__).parent



class TabPFNPredictor(BasePredictor):
    def __init__(
        self,
        model_path: Optional[Union[str, Path]] = None,
        imputer_path: Optional[Union[str, Path]] = None,
        columns_path: Optional[Union[str, Path]] = None,
        variant: Optional[str] = None,
        threshold: Optional[float] = None,
        device: str = "cpu",
    ):
        self.variant = variant
        self.device = device
        self.default_threshold = threshold if threshold is not None else 0.660658

        # Resolve variant directory if variant is specified
        variant_dir = DEFAULT_TABPFN_DIR
        if variant:
            cand_variant_dir = DEFAULT_TABPFN_DIR / variant
            if cand_variant_dir.is_dir():
                variant_dir = cand_variant_dir

        if model_path is None:
            # Check for best model pipeline first, then finetuned joblib
            if (variant_dir / "tabpfn_best_model_pipeline.joblib").exists():
                model_path = variant_dir / "tabpfn_best_model_pipeline.joblib"
            else:
                model_path = variant_dir / "tabpfn_finetuned.joblib"

        if imputer_path is None:
            imputer_path = variant_dir / "imputer.joblib"
        if columns_path is None:
            if (variant_dir / "feature_columns_best_combo.txt").exists():
                columns_path = variant_dir / "feature_columns_best_combo.txt"
            else:
                columns_path = variant_dir / "feature_columns.txt"

        self.model_path = Path(model_path)
        self.imputer_path = Path(imputer_path)
        self.columns_path = Path(columns_path)
        self.threshold = self.default_threshold

        if self.device == "cpu":
            apply_cpu_patches()

        self._load_components()

    def _load_components(self):
        setup_tabpfn_shims()
        # Check if model_path points to an all-in-one pipeline artifact dictionary

        if self.model_path.exists():
            try:
                raw = joblib.load(self.model_path)
            except Exception:
                with open(self.model_path, "rb") as f:
                    raw = pickle.load(f)

            if isinstance(raw, dict) and "model" in raw:
                self.model = raw["model"]
                self.imputer = raw.get("imputer", None)
                self.selector = raw.get("selector", None)
                self.expected_columns = raw.get("best_columns", raw.get("combo_cols", raw.get("feature_columns", [])))
                if "threshold" in raw:
                    self.threshold = float(raw["threshold"])
                self._apply_device_fixes()
                return

            self.model = raw
        else:
            raise FileNotFoundError(f"Model path does not exist: {self.model_path}")

        # If not an all-in-one dictionary, load individual components
        # 1) Load expected feature columns
        if not self.columns_path.exists():
            raise FileNotFoundError(f"Columns path does not exist: {self.columns_path}")
        with open(self.columns_path) as f:
            self.expected_columns = [line.strip() for line in f if line.strip()]

        # 2) Load SimpleImputer
        if not self.imputer_path.exists():
            raise FileNotFoundError(f"Imputer path does not exist: {self.imputer_path}")
        try:
            self.imputer = joblib.load(self.imputer_path)
        except Exception:
            with open(self.imputer_path, "rb") as f:
                self.imputer = pickle.load(f)

        self._apply_device_fixes()

    def _apply_device_fixes(self):
        # Handle CPU-specific relocations
        if self.device == "cpu" and self.model is not None:
            force_cpu(self.model)
            if hasattr(self.model, "device"):
                self.model.device = "cpu"
            if hasattr(self.model, "model_") and hasattr(self.model.model_, "to"):
                self.model.model_ = self.model.model_.to("cpu")
            if hasattr(self.model, "model") and hasattr(self.model.model, "to"):
                self.model.model = self.model.model.to("cpu")

    def predict_sample(self, tsv_path: Union[str, Path]) -> dict[str, Any]:
        """Evaluate a single sample TSV file."""
        tsv_path = Path(tsv_path)
        if not tsv_path.exists():
            raise FileNotFoundError(f"Sample TSV not found: {tsv_path}")

        # Format features to wide vector matching expected columns schema
        wide_vector = process_single_sample(tsv_path, self.expected_columns)
        X_sample = wide_vector.values.reshape(1, -1)

        # Impute missing features
        if self.imputer is not None:
            X_sample = self.imputer.transform(X_sample).astype(np.float32, copy=False)
        else:
            X_sample = np.nan_to_num(X_sample, nan=0.0).astype(np.float32, copy=False)

        # Apply feature selector if present
        if hasattr(self, "selector") and self.selector is not None:
            X_sample = self.selector.transform(X_sample).astype(np.float32, copy=False)

        # Predict probability
        prob = float(self.model.predict_proba(X_sample)[0, 1])
        pred_class = 1 if prob >= self.threshold else 0
        pred_label = "MSI" if pred_class == 1 else "MSS"

        return {
            "sample_id": tsv_path.stem,
            "file_path": str(tsv_path.resolve()),
            "msi_status": pred_label,
            "msi_score": prob,
            "p_msi": prob,
            "threshold": self.threshold,
            "y_pred": pred_class,
            "prediction": pred_label,
        }


    def get_baseline_matrix(self, n_background: int = 30) -> np.ndarray:
        """Returns the background baseline matrix for Shapley attribution."""
        if hasattr(self.imputer, "statistics_") and self.imputer.statistics_ is not None:
            bg_vec = self.imputer.statistics_
        else:
            bg_vec = np.zeros(len(self.expected_columns), dtype=np.float32)
        return np.tile(bg_vec, (max(5, min(n_background, 30)), 1)).astype(np.float32)

    def explain_sample(
        self,
        tsv_path: Union[str, Path],
        budget: int = 128,
        n_background: int = 30,
        top_n: int = 15,
    ) -> dict[str, Any]:
        """
        Evaluate a single sample TSV file and compute ShapIQ / Shapley site attributions.
        """
        from stride.core.explainability import (
            aggregate_features_to_sites,
            build_waterfall_figure,
            compute_sample_shapley_values,
            extract_positive_probs,
        )

        tsv_path = Path(tsv_path)
        if not tsv_path.exists():
            raise FileNotFoundError(f"Sample TSV not found: {tsv_path}")

        wide_vector = process_single_sample(tsv_path, self.expected_columns)
        X_sample = wide_vector.values.reshape(1, -1)

        if self.imputer is not None:
            X_imputed = self.imputer.transform(X_sample).astype(np.float32, copy=False)

        else:
            X_imputed = np.nan_to_num(X_sample, nan=0.0).astype(np.float32, copy=False)

        class ModelWithSelector:
            def __init__(self, model, selector):
                self.model = model
                self.selector = selector

            def predict_proba(self, X):
                X_sel = self.selector.transform(X) if self.selector is not None else X
                return self.model.predict_proba(X_sel)

        effective_model = ModelWithSelector(self.model, getattr(self, "selector", None))
        prob = float(extract_positive_probs(effective_model.predict_proba(X_imputed))[0])
        pred_class = 1 if prob >= self.threshold else 0
        pred_label = "MSI" if pred_class == 1 else "MSS"

        # Shapley computation
        bg_matrix = self.get_baseline_matrix(n_background=n_background)
        base_prob = float(extract_positive_probs(effective_model.predict_proba(bg_matrix)).mean())

        phi_features = compute_sample_shapley_values(
            model=effective_model,
            x_sample=X_imputed[0],
            background_matrix=bg_matrix,
            feature_names=self.expected_columns,
            budget=budget,
        )


        site_attributions = aggregate_features_to_sites(
            feature_names=self.expected_columns,
            x_sample=X_imputed[0],
            phi_features=phi_features,
        )

        fig = build_waterfall_figure(
            sample_id=tsv_path.stem,
            p_msi=prob,
            threshold=self.threshold,
            base_prob=base_prob,
            site_attributions=site_attributions,
            top_n=top_n,
        )

        top_site = site_attributions[0]["site_id"] if site_attributions else "N/A"
        top_phi = site_attributions[0]["phi"] if site_attributions else 0.0

        return {
            "sample_id": tsv_path.stem,
            "file_path": str(tsv_path.resolve()),
            "p_msi": prob,
            "y_pred": pred_class,
            "prediction": pred_label,
            "threshold": self.threshold,
            "base_prob": base_prob,
            "site_attributions": site_attributions,
            "top_driver_site": top_site,
            "top_driver_phi": top_phi,
            "waterfall_fig": fig,
        }

    def predict_batch(
        self, tsv_paths: list[Union[str, Path]], output_tsv: Optional[Union[str, Path]] = None
    ) -> pd.DataFrame:
        """Evaluate multiple sample TSV files."""
        results = []
        for path in tsv_paths:
            path = Path(path)
            try:
                res = self.predict_sample(path)
                results.append(res)
            except Exception as e:
                results.append(
                    {
                        "sample_id": path.stem,
                        "file_path": str(path.resolve()),
                        "p_msi": np.nan,
                        "y_pred": -1,
                        "prediction": f"ERROR: {str(e)}",
                    }
                )

        df = pd.DataFrame(results)

        if output_tsv is not None:
            output_tsv = Path(output_tsv)
            output_tsv.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_tsv, sep="\t", index=False)
            print(f"Saved evaluation results to: {output_tsv}")

        return df

    def predict_batch_from_dir(
        self,
        input_dir: Union[str, Path],
        file_ext: str = "tsv",
        output_tsv: Optional[Union[str, Path]] = None,
    ) -> pd.DataFrame:
        """Scan directory for TSVs and evaluate them in batch."""
        input_dir = Path(input_dir)
        if not input_dir.is_dir():
            raise FileNotFoundError(f"Input directory not found: {input_dir}")

        tsv_files = sorted(input_dir.glob(f"*.{file_ext}"))
        print(f"Found {len(tsv_files)} files to process in {input_dir}")
        return self.predict_batch(tsv_files, output_tsv=output_tsv)

