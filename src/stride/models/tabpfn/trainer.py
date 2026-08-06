import json
import joblib
import pickle
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
from typing import List, Optional, Tuple, Dict, Any

from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold, StratifiedShuffleSplit

# Robust import of FinetunedTabPFNClassifier
FinetunedTabPFNClassifier = None
try:
    from tabpfn import FinetunedTabPFNClassifier
except ImportError:
    try:
        from tabpfn.finetuning.finetuned_classifier import FinetunedTabPFNClassifier
    except ImportError:
        try:
            from tabpfn.finetuning import FinetunedTabPFNClassifier
        except ImportError:
            pass

from stride.core.base import BaseTrainer
from stride.core.dataset import build_matrix
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

class TabPFNTrainer(BaseTrainer):
    def __init__(
        self,
        cv_folds: int = 5,
        cv_repeats: int = 1,
        cv_seed: int = 42,
        min_spec: float = 0.95,
        device: str = "cpu",
        refit_val_frac: float = 0.15,
        max_sites: int = 3000,
        rank_by: str = "sensitivity"
    ):
        self.cv_folds = cv_folds
        self.cv_repeats = cv_repeats
        self.cv_seed = cv_seed
        self.min_spec = min_spec
        self.device = device
        self.refit_val_frac = refit_val_frac
        self.max_sites = max_sites
        self.rank_by = rank_by

        if FinetunedTabPFNClassifier is None:
            raise ImportError(
                "FinetunedTabPFNClassifier could not be imported. Please ensure tabpfn is installed in the current environment."
            )

    def train_cv_refit(
        self,
        impact_records: pd.DataFrame,
        access_records: pd.DataFrame,
        test_records: pd.DataFrame,
        metric_pool: List[str],
        outdir: Path,
        epochs: int = 30,
        learning_rate: float = 1e-5,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Train and evaluate the TabPFN pipeline.
        1) Build matrices for Impact, Access, and Test data.
        2) Fit SimpleImputer and transform features.
        3) Run Cross-Validation on Access data (with Impact as fixed training).
        4) Pick decision threshold based on Access OOF predictions.
        5) Refit on combined (Impact + Access) training data.
        6) Evaluate on Test data.
        7) Save all trained model components and metadata to outdir.
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        set_seed(self.cv_seed)

        print("Building datasets...")
        X_impact, y_impact, expected_sites = build_matrix(
            impact_records,
            metric_pool,
            expected_sites=None,
            max_sites=self.max_sites,
            seed=self.cv_seed
        )
        X_access, y_access, _ = build_matrix(
            access_records,
            metric_pool,
            expected_sites=expected_sites,
            max_sites=0,
            seed=self.cv_seed
        )
        X_test, y_test, _ = build_matrix(
            test_records,
            metric_pool,
            expected_sites=expected_sites,
            max_sites=0,
            seed=self.cv_seed
        )

        # Retain common columns
        common_cols = X_impact.columns.intersection(X_access.columns).intersection(X_test.columns)
        if len(common_cols) == 0:
            raise RuntimeError("No common columns across datasets.")
        X_impact = X_impact[common_cols]
        X_access = X_access[common_cols]
        X_test = X_test[common_cols]
        feature_columns = list(common_cols)

        print(f"Matrix shape: impact={X_impact.shape}, access={X_access.shape}, test={X_test.shape}")
        
        # Fit Imputer
        imputer = SimpleImputer(strategy="median")
        X_impact_imp = imputer.fit_transform(X_impact).astype(np.float32)
        X_access_imp = imputer.transform(X_access).astype(np.float32)
        X_test_imp = imputer.transform(X_test).astype(np.float32)

        # Cross Validation Setup
        if self.cv_repeats > 1:
            cv = RepeatedStratifiedKFold(n_splits=self.cv_folds, n_repeats=self.cv_repeats, random_state=self.cv_seed)
        else:
            cv = StratifiedKFold(n_splits=self.cv_folds, shuffle=True, random_state=self.cv_seed)

        oof_probs = np.zeros(len(X_access_imp))
        print(f"Starting {self.cv_folds}-fold CV (repeats={self.cv_repeats})...")
        split_count = 0
        for tr_idx, val_idx in cv.split(X_access_imp, y_access):
            split_count += 1
            print(f"Fold {split_count} / {self.cv_folds * self.cv_repeats}")
            
            # Combine Impact with Access train fold
            X_tr_fold = np.concatenate([X_impact_imp, X_access_imp[tr_idx]], axis=0)
            y_tr_fold = np.concatenate([y_impact, y_access[tr_idx]], axis=0)
            
            X_va_fold = X_access_imp[val_idx]
            y_va_fold = y_access[val_idx]

            # Fit TabPFN Fold Model
            fold_model = FinetunedTabPFNClassifier(
                device=self.device,
                epochs=epochs,
                learning_rate=learning_rate
            )
            fold_model.fit(X_tr_fold, y_tr_fold, X_val=X_va_fold, y_val=y_va_fold)
            
            # Predict probabilities
            probs_va = fold_model.predict_proba(X_va_fold)[:, 1]
            if self.cv_repeats == 1:
                oof_probs[val_idx] = probs_va
            else:
                oof_probs[val_idx] += probs_va / self.cv_repeats

            release_runtime_memory(self.device)

        # Select decision threshold to satisfy min_spec
        best_threshold, best_sens, best_spec = select_decision_threshold(
            oof_probs, y_access, min_spec=self.min_spec, rank_by=self.rank_by
        )
        print(f"Selected threshold: {best_threshold:.4f} (OOF Sens: {best_sens:.4f}, Spec: {best_spec:.4f})")

        # 5) Refit final model on all data (Impact + Access)
        print("Refitting final model on all training data...")
        X_train_full = np.concatenate([X_impact_imp, X_access_imp], axis=0)
        y_train_full = np.concatenate([y_impact, y_access], axis=0)

        refit_X = X_train_full
        refit_y = y_train_full
        refit_val_X = X_train_full
        refit_val_y = y_train_full

        if 0.0 < self.refit_val_frac < 1.0:
            try:
                sss = StratifiedShuffleSplit(n_splits=1, test_size=self.refit_val_frac, random_state=self.cv_seed)
                tr_idx, va_idx = next(sss.split(X_train_full, y_train_full))
                refit_X = X_train_full[tr_idx]
                refit_y = y_train_full[tr_idx]
                refit_val_X = X_train_full[va_idx]
                refit_val_y = y_train_full[va_idx]
            except ValueError:
                print("WARNING: Could not create stratified refit split; using full train.")

        final_model = FinetunedTabPFNClassifier(
            device=self.device,
            epochs=epochs,
            learning_rate=learning_rate
        )
        
        final_model_dir = outdir / "final_model"
        final_model_dir.mkdir(parents=True, exist_ok=True)
        final_model.fit(refit_X, refit_y, X_val=refit_val_X, y_val=refit_val_y, output_dir=final_model_dir)

        # 6) Evaluate on Test data
        p_test = final_model.predict_proba(X_test_imp)[:, 1]
        yhat_test = (p_test >= best_threshold).astype(int)

        test_tp, test_fp, test_tn, test_fn = confusion_counts(y_test, yhat_test)
        test_met = metrics_from_counts(test_tp, test_fp, test_tn, test_fn)
        test_roc = safe_roc_auc(y_test, p_test)
        test_pr = safe_pr_auc(y_test, p_test)

        print(f"Test Score: ROC-AUC={test_roc:.4f}, Spec={test_met['specificity']:.4f}, Sens={test_met['sensitivity']:.4f}")

        # 7) Save all trained model components
        print(f"Saving outputs to {outdir}...")
        
        # Save imputer
        joblib.dump(imputer, outdir / "imputer.joblib")
        with open(outdir / "imputer.pkl", "wb") as f:
            pickle.dump(imputer, f)

        # Save feature column lists and site IDs
        (outdir / "feature_columns_best_combo.txt").write_text("\n".join(feature_columns) + "\n")
        (outdir / "site_ids.txt").write_text("\n".join(expected_sites) + "\n")
        (outdir / "metric_pool.txt").write_text("\n".join(metric_pool) + "\n")

        # Save TabPFN model joblib
        joblib.dump(final_model, outdir / "tabpfn_finetuned.joblib")
        with open(outdir / "tabpfn_finetuned.pkl", "wb") as f:
            pickle.dump(final_model, f)

        # Save model state dict if possible
        if hasattr(final_model, "model_") and hasattr(final_model.model_, "state_dict"):
            try:
                import torch
                torch.save(final_model.model_.state_dict(), outdir / "tabpfn_finetuned_state_dict.pt")
            except Exception as e:
                print(f"WARNING: Could not save torch state_dict: {e}")

        # Generate run manifest & artifacts payload
        artifacts = {
            "created": datetime.now().isoformat(timespec="seconds"),
            "epochs": epochs,
            "learning_rate": learning_rate,
            "best_threshold": float(best_threshold),
            "cv_metrics": {
                "oof_roc_auc": float(safe_roc_auc(y_access, oof_probs)),
                "oof_pr_auc": float(safe_pr_auc(y_access, oof_probs)),
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
        oof_df["p_msi"] = oof_probs
        oof_df["y_pred"] = (oof_probs >= best_threshold).astype(int)
        oof_df.to_csv(outdir / "access_oof_predictions.tsv", sep="\t", index=False)

        test_df = test_records.copy()
        test_df["p_msi"] = p_test
        test_df["y_pred"] = yhat_test
        test_df.to_csv(outdir / "test_predictions.tsv", sep="\t", index=False)

        print("Training run completed successfully.")
        return artifacts
