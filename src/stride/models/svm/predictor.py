import json
import pickle
from pathlib import Path
from typing import Any, Optional, Union

import joblib
import numpy as np
import pandas as pd

from stride.core.base import BasePredictor
from stride.core.dataset import process_single_sample

DEFAULT_SVM_DIR = Path(__file__).parent


class SVMPredictor(BasePredictor):
    def __init__(
        self,
        model_path: Optional[Union[str, Path]] = None,
        columns_path: Optional[Union[str, Path]] = None,
        sign_map_path: Optional[Union[str, Path]] = None,
        threshold: float = 0.0,
    ):
        if model_path is None:
            model_path = DEFAULT_SVM_DIR / "svm_incremental.joblib"
        if columns_path is None:
            columns_path = DEFAULT_SVM_DIR / "feature_columns_best_combo.txt"

        self.model_path = Path(model_path)
        self.columns_path = Path(columns_path)
        self.sign_map_path = Path(sign_map_path) if sign_map_path else None
        self.threshold = threshold

        self._load_components()

    def _load_components(self):
        # 1) Load expected feature columns
        if not self.columns_path.exists():
            raise FileNotFoundError(f"Columns path does not exist: {self.columns_path}")
        with open(self.columns_path) as f:
            self.expected_columns = [line.strip() for line in f if line.strip()]

        # 2) Load sign map if present
        self.sign_map = None
        if self.sign_map_path and self.sign_map_path.exists():
            with open(self.sign_map_path) as f:
                self.sign_map = json.load(f)

        # 3) Load SVM model
        if not self.model_path.exists():
            raise FileNotFoundError(f"Model path does not exist: {self.model_path}")
        try:
            self.model = joblib.load(self.model_path)
        except Exception:
            with open(self.model_path, "rb") as f:
                self.model = pickle.load(f)

    def _get_scores(self, X: np.ndarray) -> float:
        try:
            return float(self.model.decision_function(X)[0])
        except Exception:
            try:
                proba = self.model.predict_proba(X)
                return float(proba[0, 1])
            except Exception:
                return float(self.model.predict(X)[0])

    def predict_sample(self, tsv_path: Union[str, Path]) -> dict[str, Any]:
        """Evaluate a single sample TSV file."""
        tsv_path = Path(tsv_path)
        if not tsv_path.exists():
            raise FileNotFoundError(f"Sample TSV not found: {tsv_path}")

        # Format features to wide vector matching expected columns schema, applying sign map if present
        wide_vector = process_single_sample(tsv_path, self.expected_columns, sign_map=self.sign_map)
        X_sample = wide_vector.values.reshape(1, -1)

        # Impute missing values (represented as nan) with 0.0 (since SVM matrices are zero-filled)
        X_filled = np.nan_to_num(X_sample, nan=0.0)

        # Predict score
        score = self._get_scores(X_filled)
        pred_class = 1 if score >= self.threshold else 0
        pred_label = "MSI" if pred_class == 1 else "MSS"

        return {
            "sample_id": tsv_path.stem,
            "file_path": str(tsv_path.resolve()),
            "score": score,
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
                        "score": np.nan,
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
