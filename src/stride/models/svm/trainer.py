import json
import joblib
import pickle
import copy
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
from typing import List, Optional, Tuple, Dict, Any

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import SGDClassifier
from sklearn.utils.class_weight import compute_class_weight
from sklearn.utils import shuffle as sk_shuffle
from sklearn.model_selection import StratifiedKFold

from stride.core.base import BaseTrainer
from stride.core.dataset import build_matrix, load_sfs_sign_map
from stride.core.metrics import (
    set_seed,
    release_runtime_memory,
    confusion_counts,
    metrics_from_counts,
    safe_roc_auc,
    safe_pr_auc,
    compute_metrics
)
from stride.core.qc import select_decision_threshold


class SVMTrainer(BaseTrainer):
    def __init__(
        self,
        cv_folds: int = 5,
        cv_seed: int = 42,
        incremental_epochs: int = 5,
        incremental_batch_size: int = 256,
        incremental_weight: float = 1.0,
        update_scaler: bool = False,
        min_spec: float = 0.95,
        rank_by: str = "sensitivity",
        sfs_tsv: Optional[str] = None,
        sfs_rank1_only: bool = True
    ):
        self.cv_folds = cv_folds
        self.cv_seed = cv_seed
        self.incremental_epochs = incremental_epochs
        self.incremental_batch_size = incremental_batch_size
        self.incremental_weight = incremental_weight
        self.update_scaler = update_scaler
        self.min_spec = min_spec
        self.rank_by = rank_by
        self.sfs_tsv = sfs_tsv
        self.sfs_rank1_only = sfs_rank1_only

    def _train_base_on_impact(self, X: pd.DataFrame, y: np.ndarray, random_state: int = 42):
        classes = np.array([0, 1])
        class_weights = compute_class_weight(class_weight="balanced", classes=classes, y=y)
        wmap = dict(zip(classes, class_weights))
        w = np.array([wmap[i] for i in y])

        scaler = StandardScaler()
        Xs = scaler.fit_transform(X)

        clf = SGDClassifier(
            loss="hinge",
            max_iter=5000,
            tol=1e-4,
            random_state=random_state,
            average=True,
        )
        clf.partial_fit(Xs, y, classes=classes, sample_weight=w)
        pipe = make_pipeline(scaler, clf)
        return pipe, classes

    def _incremental_update_on_access(
        self,
        model,
        X_B: pd.DataFrame,
        y_B: np.ndarray,
        classes: np.ndarray,
        epochs: int = 5,
        batch_size: int = 256,
        update_scaler: bool = False,
        incremental_weight: float = 1.0,
    ):
        scaler = model.named_steps["standardscaler"]
        clf = model.named_steps["sgdclassifier"]

        class_weights = compute_class_weight(class_weight="balanced", classes=classes, y=y_B)
        wmap = dict(zip(classes, class_weights))
        base_w = np.array([wmap[i] for i in y_B], dtype=float) * float(incremental_weight)

        idx = np.arange(len(y_B))
        for ep in range(int(epochs)):
            idx = sk_shuffle(idx, random_state=self.cv_seed + ep)
            for start in range(0, len(idx), batch_size):
                sl = idx[start : start + batch_size]
                Xb = X_B.iloc[sl].values
                yb = y_B[sl]
                wb = base_w[sl]
                if update_scaler:
                    scaler.partial_fit(Xb)
                Xb_scaled = scaler.transform(Xb)
                clf.partial_fit(Xb_scaled, yb, classes=classes, sample_weight=wb)
        return model

    def _get_scores(self, model, X: pd.DataFrame) -> np.ndarray:
        try:
            return model.decision_function(X)
        except Exception:
            try:
                proba = model.predict_proba(X)
                return proba[:, 1]
            except Exception:
                return model.predict(X).astype(float)

    def train_cv_refit(
        self,
        impact_records: pd.DataFrame,
        access_records: pd.DataFrame,
        test_records: pd.DataFrame,
        metric_pool: List[str],
        outdir: Path,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Train and evaluate the SVM incremental model pipeline.
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        set_seed(self.cv_seed)

        sign_map = None
        if self.sfs_tsv:
            print(f"Loading SFS sign map from {self.sfs_tsv}...")
            sign_map = load_sfs_sign_map(self.sfs_tsv, rank1_only=self.sfs_rank1_only)

        print("Building datasets...")
        X_impact, y_impact, expected_sites = build_matrix(
            impact_records,
            metric_pool,
            expected_sites=None,
            max_sites=0,
            seed=self.cv_seed,
            sign_map=sign_map
        )
        X_access, y_access, _ = build_matrix(
            access_records,
            metric_pool,
            expected_sites=expected_sites,
            max_sites=0,
            seed=self.cv_seed,
            sign_map=sign_map
        )
        X_test, y_test, _ = build_matrix(
            test_records,
            metric_pool,
            expected_sites=expected_sites,
            max_sites=0,
            seed=self.cv_seed,
            sign_map=sign_map
        )

        # Retain common columns
        common_cols = X_impact.columns.intersection(X_access.columns).intersection(X_test.columns)
        if len(common_cols) == 0:
            raise RuntimeError("No common columns across datasets.")
        X_impact = X_impact[common_cols].fillna(0.0)
        X_access = X_access[common_cols].fillna(0.0)
        X_test = X_test[common_cols].fillna(0.0)
        feature_columns = list(common_cols)

        print(f"Matrix shape: impact={X_impact.shape}, access={X_access.shape}, test={X_test.shape}")

        # Train BASE model on Impact once
        print("Training base model on Impact cohort...")
        base_model, classes = self._train_base_on_impact(X_impact, y_impact, random_state=self.cv_seed)

        # Cross Validation Setup on Access
        cv = StratifiedKFold(n_splits=self.cv_folds, shuffle=True, random_state=self.cv_seed)
        oof_scores = np.zeros(len(X_access))

        print(f"Starting {self.cv_folds}-fold CV on Access...")
        for fold, (tr_idx, val_idx) in enumerate(cv.split(X_access, y_access), start=1):
            X_tr = X_access.iloc[tr_idx]
            y_tr = y_access[tr_idx]
            X_va = X_access.iloc[val_idx]
            
            # Copy base weights to ensure starting from same fixed Impact pretrained model
            fold_model = copy.deepcopy(base_model)
            fold_model = self._incremental_update_on_access(
                fold_model,
                X_tr,
                y_tr,
                classes,
                epochs=self.incremental_epochs,
                batch_size=self.incremental_batch_size,
                update_scaler=self.update_scaler,
                incremental_weight=self.incremental_weight
            )

            # Predict scores on validation fold
            oof_scores[val_idx] = self._get_scores(fold_model, X_va)
            release_runtime_memory()

        # Select decision threshold based on Access OOF scores to meet min_spec constraint
        best_threshold, best_sens, best_spec = select_decision_threshold(
            oof_scores, y_access, min_spec=self.min_spec, rank_by=self.rank_by
        )
        print(f"Selected threshold: {best_threshold:.4f} (OOF Sens: {best_sens:.4f}, Spec: {best_spec:.4f})")

        # Refit final model on all data (Impact + Access)
        print("Refitting final incremental model on all Access data...")
        final_model = copy.deepcopy(base_model)
        final_model = self._incremental_update_on_access(
            final_model,
            X_access,
            y_access,
            classes,
            epochs=self.incremental_epochs,
            batch_size=self.incremental_batch_size,
            update_scaler=self.update_scaler,
            incremental_weight=self.incremental_weight
        )

        # Evaluate on Test data
        p_test = self._get_scores(final_model, X_test)
        yhat_test = (p_test >= best_threshold).astype(int)

        test_tp, test_fp, test_tn, test_fn = confusion_counts(y_test, yhat_test)
        test_met = metrics_from_counts(test_tp, test_fp, test_tn, test_fn)
        test_roc = safe_roc_auc(y_test, p_test)
        test_pr = safe_pr_auc(y_test, p_test)

        print(f"Test Score: ROC-AUC={test_roc:.4f}, Spec={test_met['specificity']:.4f}, Sens={test_met['sensitivity']:.4f}")

        # Save all trained model components
        print(f"Saving outputs to {outdir}...")
        
        # Save base model & final incremental model
        joblib.dump(base_model, outdir / "svm_base.joblib")
        joblib.dump(final_model, outdir / "svm_incremental.joblib")
        
        # Save feature column lists and site IDs
        (outdir / "feature_columns_best_combo.txt").write_text("\n".join(feature_columns) + "\n")
        (outdir / "site_ids.txt").write_text("\n".join(expected_sites) + "\n")
        (outdir / "metric_pool.txt").write_text("\n".join(metric_pool) + "\n")

        # Save sign map JSON if applicable
        if sign_map:
            with open(outdir / "sign_map.json", "w") as f:
                json.dump(sign_map, f, indent=2)

        # Generate run manifest & artifacts payload
        artifacts = {
            "created": datetime.now().isoformat(timespec="seconds"),
            "epochs": self.incremental_epochs,
            "best_threshold": float(best_threshold),
            "cv_metrics": {
                "oof_roc_auc": float(safe_roc_auc(y_access, oof_scores)),
                "oof_pr_auc": float(safe_pr_auc(y_access, oof_scores)),
            },
            "test_metrics": {
                "roc_auc": float(test_roc),
                "pr_auc": float(test_pr),
                "tp": int(test_tp),
                "fp": int(test_fp),
                "tn": int(test_tn),
                "fn": int(test_fn),
                **test_met
            },
            "expected_sites_count": len(expected_sites),
            "feature_columns_count": len(feature_columns)
        }

        with open(outdir / "artifacts.json", "w") as f:
            json.dump(artifacts, f, indent=2)

        # Save OOF predictions and test predictions as TSV files
        oof_df = access_records.copy()
        oof_df["score"] = oof_scores
        oof_df["y_pred"] = (oof_scores >= best_threshold).astype(int)
        oof_df.to_csv(outdir / "access_oof_predictions.tsv", sep="\t", index=False)

        test_df = test_records.copy()
        test_df["score"] = p_test
        test_df["y_pred"] = yhat_test
        test_df.to_csv(outdir / "test_predictions.tsv", sep="\t", index=False)

        print("SVM Training run completed successfully.")
        return artifacts
