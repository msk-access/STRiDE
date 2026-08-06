import random
import gc
import numpy as np
from typing import Dict, Tuple, Optional, Any

# Robust imports of sklearn metrics
try:
    from sklearn.metrics import (
        average_precision_score, 
        roc_auc_score, 
        accuracy_score, 
        precision_score, 
        recall_score, 
        balanced_accuracy_score,
        f1_score, 
        matthews_corrcoef,
        confusion_matrix
    )
except ImportError:
    average_precision_score = None
    roc_auc_score = None
    accuracy_score = None
    precision_score = None
    recall_score = None
    balanced_accuracy_score = None
    f1_score = None
    matthews_corrcoef = None
    confusion_matrix = None

# Robust import of torch
try:
    import torch
except ImportError:
    torch = None


def confusion_counts(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[int, int, int, int]:
    """Calculate counts of TP, FP, TN, FN."""
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tp, fp, tn, fn


def metrics_from_counts(tp: int, fp: int, tn: int, fn: int) -> Dict[str, float]:
    """Derive metric scores from TP, FP, TN, FN counts."""
    sens = tp / (tp + fn + 1e-12)
    spec = tn / (tn + fp + 1e-12)
    prec = tp / (tp + fp + 1e-12)
    acc = (tp + tn) / max(tp + fp + tn + fn, 1)
    f1 = 2 * tp / (2 * tp + fp + fn + 1e-12)
    bal = (sens + spec) / 2.0
    return {
        "sensitivity": float(sens),
        "specificity": float(spec),
        "precision": float(prec),
        "accuracy": float(acc),
        "f1_score": float(f1),
        "balanced_accuracy": float(bal),
    }


def safe_roc_auc(y_true: np.ndarray, p: np.ndarray) -> float:
    if roc_auc_score is None:
        return float("nan")
    try:
        if len(np.unique(y_true)) > 1:
            return float(roc_auc_score(y_true, p))
        return float("nan")
    except Exception:
        return float("nan")


def safe_pr_auc(y_true: np.ndarray, p: np.ndarray) -> float:
    if average_precision_score is None:
        return float("nan")
    try:
        if len(np.unique(y_true)) > 1:
            return float(average_precision_score(y_true, p))
        return float("nan")
    except Exception:
        return float("nan")


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray, scores: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """
    Comprehensive performance metrics dictionary combining scikit-learn standard metrics 
    and detailed confusion counts.
    """
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    
    tp, fp, tn, fn = confusion_counts(y_true, y_pred)
    derived = metrics_from_counts(tp, fp, tn, fn)
    
    out = {
        "accuracy": float(accuracy_score(y_true, y_pred)) if accuracy_score else derived["accuracy"],
        "precision": float(precision_score(y_true, y_pred, zero_division=0)) if precision_score else derived["precision"],
        "recall_sensitivity": float(recall_score(y_true, y_pred, zero_division=0)) if recall_score else derived["sensitivity"],
        "specificity": float(recall_score(y_true, y_pred, pos_label=0, zero_division=0)) if recall_score else derived["specificity"],
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)) if balanced_accuracy_score else derived["balanced_accuracy"],
        "f1_score": float(f1_score(y_true, y_pred, zero_division=0)) if f1_score else derived["f1_score"],
    }
    
    if matthews_corrcoef:
        try:
            out["mcc"] = float(matthews_corrcoef(y_true, y_pred))
        except Exception:
            out["mcc"] = float("nan")
    else:
        out["mcc"] = float("nan")
        
    if scores is not None:
        out["auc"] = safe_roc_auc(y_true, scores)
        out["pr_auc"] = safe_pr_auc(y_true, scores)
    else:
        out["auc"] = float("nan")
        out["pr_auc"] = float("nan")
        
    out.update({
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "n_pos": int((y_true == 1).sum()),
        "n_neg": int((y_true == 0).sum()),
        "n_total": len(y_true)
    })
    
    return out


def set_seed(seed: int, deterministic: bool = False) -> None:
    """Set standard random seeds across random, numpy, and torch."""
    random.seed(seed)
    np.random.seed(seed)
    if torch is not None:
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        if deterministic:
            try:
                torch.backends.cudnn.deterministic = True
                torch.backends.cudnn.benchmark = False
            except Exception:
                pass


def release_runtime_memory(device: Optional[str] = None) -> None:
    """Explicitly call garbage collection and empty CUDA cache if PyTorch is available."""
    gc.collect()
    if torch is not None and device == "cuda" and torch.cuda.is_available():
        try:
            torch.cuda.empty_cache()
        except Exception:
            pass
