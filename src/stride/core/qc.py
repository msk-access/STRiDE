from typing import Any, Optional

import numpy as np

from stride.core.metrics import compute_metrics, confusion_counts, metrics_from_counts


def select_decision_threshold(
    oof_scores: np.ndarray,
    y_true: np.ndarray,
    min_spec: float = 0.95,
    rank_by: str = "sensitivity",
    grid_steps: int = 200,
) -> tuple[float, float, float]:
    """
    Shared threshold scanning calculation:
    Scans candidate decision thresholds on cross-validation out-of-fold scores
    to identify the threshold meeting the minimum specificity requirement (min_spec).

    Returns:
        Tuple[float, float, float]: (best_threshold, sensitivity_at_threshold, specificity_at_threshold)
    """
    oof_scores = np.asarray(oof_scores, dtype=float)
    y_true = np.asarray(y_true, dtype=int)

    if len(oof_scores) == 0:
        return 0.5, 0.0, 0.0

    min_score = float(np.min(oof_scores))
    max_score = float(np.max(oof_scores))

    if min_score == max_score:
        threshold_grid = np.array([min_score])
    else:
        threshold_grid = np.linspace(min_score - 0.05, max_score + 0.05, grid_steps)

    passed_thresholds = []

    for thr in threshold_grid:
        y_pred = (oof_scores >= thr).astype(int)
        tp, fp, tn, fn = confusion_counts(y_true, y_pred)
        met = metrics_from_counts(tp, fp, tn, fn)
        if met["specificity"] >= min_spec:
            rank_val = met.get(rank_by, met["sensitivity"])
            passed_thresholds.append((thr, rank_val, met["specificity"]))

    if passed_thresholds:
        passed_thresholds.sort(key=lambda x: x[1], reverse=True)
        best_threshold, best_sens, best_spec = passed_thresholds[0]
        return float(best_threshold), float(best_sens), float(best_spec)

    # Fallback: select threshold that maximizes specificity
    best_thr_fallback = float(threshold_grid[0])
    best_spec_fallback = -1.0
    for thr in threshold_grid:
        y_pred = (oof_scores >= thr).astype(int)
        tp, fp, tn, fn = confusion_counts(y_true, y_pred)
        spec = metrics_from_counts(tp, fp, tn, fn)["specificity"]
        if spec > best_spec_fallback:
            best_spec_fallback = spec
            best_thr_fallback = float(thr)

    y_pred_fb = (oof_scores >= best_thr_fallback).astype(int)
    tp, fp, tn, fn = confusion_counts(y_true, y_pred_fb)
    sens_fallback = metrics_from_counts(tp, fp, tn, fn)["sensitivity"]
    return float(best_thr_fallback), float(sens_fallback), float(best_spec_fallback)


def evaluate_model_performance(
    y_true: np.ndarray, y_pred: np.ndarray, scores: Optional[np.ndarray] = None
) -> dict[str, Any]:
    """
    Perform comprehensive quality control and metric evaluation.
    """
    return compute_metrics(y_true, y_pred, scores=scores)
