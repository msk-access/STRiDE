import pickle
from pathlib import Path
from typing import Any, Optional, Union

import joblib
import numpy as np
import pandas as pd

from stride.core.base import BasePredictor
from stride.core.dataset import process_single_sample
from stride.utils.torch_patches import apply_cpu_patches, force_cpu

DEFAULT_TABPFN_DIR = Path(__file__).parent


class TabPFNPredictor(BasePredictor):
    def __init__(
        self,
        model_path: Optional[Union[str, Path]] = None,
        imputer_path: Optional[Union[str, Path]] = None,
        columns_path: Optional[Union[str, Path]] = None,
        threshold: float = 0.660658,
        device: str = "cpu",
    ):
        if model_path is None:
            model_path = DEFAULT_TABPFN_DIR / "tabpfn_finetuned.joblib"
        if imputer_path is None:
            imputer_path = DEFAULT_TABPFN_DIR / "imputer.joblib"
        if columns_path is None:
            columns_path = DEFAULT_TABPFN_DIR / "feature_columns_best_combo.txt"

        self.model_path = Path(model_path)
        self.imputer_path = Path(imputer_path)
        self.columns_path = Path(columns_path)
        self.threshold = threshold
        self.device = device

        if self.device == "cpu":
            apply_cpu_patches()

        self._load_components()

    def _load_components(self):
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

        # 3) Load TabPFN classifier model
        if not self.model_path.exists():
            raise FileNotFoundError(f"Model path does not exist: {self.model_path}")
        try:
            self.model = joblib.load(self.model_path)
        except Exception:
            with open(self.model_path, "rb") as f:
                self.model = pickle.load(f)

        # 4) Handle CPU-specific relocations
        if self.device == "cpu":
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
        X_imputed = self.imputer.transform(X_sample).astype(np.float32, copy=False)

        # Predict probability
        prob = float(self.model.predict_proba(X_imputed)[0, 1])
        pred_class = 1 if prob >= self.threshold else 0
        pred_label = "MSI" if pred_class == 1 else "MSS"

        return {
            "sample_id": tsv_path.stem,
            "file_path": str(tsv_path.resolve()),
            "p_msi": prob,
            "y_pred": pred_class,
            "prediction": pred_label,
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
