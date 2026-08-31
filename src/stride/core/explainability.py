"""
STRiDE Model Explainability & ShapIQ / Shapley Site-Level Attribution Engine.

Computes game-theoretic Shapley attributions for TabPFN and other MSI predictors,
aggregates multi-metric feature attributions to microsatellite loci, and generates
publication-grade interactive Waterfall charts and ranked driver tables.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Library availability flags
SHAPIQ_AVAILABLE = False
try:
    from tabpfn_extensions.interpretability.shapiq import get_tabpfn_imputation_explainer
    SHAPIQ_AVAILABLE = True
except ImportError:
    try:
        import shapiq  # noqa: F401
        SHAPIQ_AVAILABLE = True
    except ImportError:
        SHAPIQ_AVAILABLE = False

SHAP_AVAILABLE = False
try:
    import shap
    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False


def extract_positive_probs(raw_probs: Any) -> np.ndarray:
    """Helper to extract positive (class 1) probabilities from model output."""
    arr = np.asarray(raw_probs)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2:
        if arr.shape[1] >= 2:
            return arr[:, 1]
        return arr[:, 0]
    return arr.ravel()


def compute_sample_shapley_values(
    model: Any,
    x_sample: np.ndarray,
    background_matrix: np.ndarray,
    feature_names: List[str],
    budget: int = 128,
) -> np.ndarray:
    """
    Computes Shapley values (phi) for a single sample across selected features.
    
    Execution Strategy:
    1. PriorLabs ShapIQ extension (if available)
    2. Standard SHAP SamplingExplainer (if available)
    3. Fast vectorized marginal permutation sampling fallback (zero extra dependencies)
    """
    n_features = len(feature_names)
    x_2d = x_sample.reshape(1, -1)

    # 1. PriorLabs ShapIQ extension
    if SHAPIQ_AVAILABLE:
        try:
            from tabpfn_extensions.interpretability.shapiq import get_tabpfn_imputation_explainer
            explainer = get_tabpfn_imputation_explainer(
                model=model,
                data=background_matrix,
                max_order=1,
                index="SV",
                budget=budget,
            )
            interaction_values = explainer.explain(x_2d)
            if hasattr(interaction_values, "get_values"):
                phi = interaction_values.get_values(order=1)
                if len(phi) == n_features:
                    return np.asarray(phi, dtype=np.float32)
        except Exception as err:
            logger.debug("PriorLabs ShapIQ failed, falling back to SHAP/sampling: %s", err)

    # 2. Standard SHAP Explainer
    if SHAP_AVAILABLE:
        try:
            def predict_fn(X_eval):
                raw = model.predict_proba(X_eval)
                return extract_positive_probs(raw)

            explainer = shap.SamplingExplainer(predict_fn, background_matrix)
            shap_vals = explainer.shap_values(x_2d, nsamples=budget)
            if isinstance(shap_vals, list):
                shap_vals = shap_vals[-1]
            return np.asarray(shap_vals[0], dtype=np.float32)
        except Exception as err:
            logger.debug("SHAP SamplingExplainer failed, falling back to sampling: %s", err)

    # 3. Vectorized Marginal Permutation Sampling Fallback
    np.random.seed(42)
    m_samples = max(32, min(budget, 128))
    phi = np.zeros(n_features, dtype=np.float32)
    bg_mean = background_matrix.mean(axis=0)

    # Predict baseline and sample full probability
    p_full = float(extract_positive_probs(model.predict_proba(x_2d))[0])
    p_base = float(extract_positive_probs(model.predict_proba(bg_mean.reshape(1, -1)))[0])

    # Sample random feature subsets
    for _ in range(m_samples):
        perm = np.random.permutation(n_features)
        split_point = np.random.randint(1, n_features)
        feat_idx = perm[0]

        # Hybrid masks
        mask_with = np.copy(bg_mean)
        mask_without = np.copy(bg_mean)

        active_indices = perm[:split_point]
        mask_with[active_indices] = x_sample[active_indices]
        mask_without[active_indices] = x_sample[active_indices]
        mask_without[feat_idx] = bg_mean[feat_idx]

        batch = np.vstack([mask_with, mask_without])
        probs = extract_positive_probs(model.predict_proba(batch))
        contrib = float(probs[0] - probs[1])
        phi[feat_idx] += contrib

    # Normalize contributions so sum approximates (p_full - p_base)
    sum_phi = np.sum(np.abs(phi))
    if sum_phi > 1e-7:
        scale = (p_full - p_base) / np.sum(phi) if np.abs(np.sum(phi)) > 1e-7 else 1.0
        phi = phi * np.clip(scale, -2.0, 2.0)

    return phi


def aggregate_features_to_sites(
    feature_names: List[str],
    x_sample: np.ndarray,
    phi_features: np.ndarray,
) -> List[Dict[str, Any]]:
    """
    Aggregates multi-metric feature attributions into per-microsatellite locus Shapley scores.
    """
    site_phi_map: Dict[str, float] = {}
    site_metrics_map: Dict[str, Dict[str, float]] = {}

    for f_idx, feat_col in enumerate(feature_names):
        phi_val = float(phi_features[f_idx])
        if "__" in feat_col:
            sid, metric = feat_col.split("__", 1)
        else:
            sid, metric = feat_col, "value"

        site_phi_map[sid] = site_phi_map.get(sid, 0.0) + phi_val
        if sid not in site_metrics_map:
            site_metrics_map[sid] = {}
        site_metrics_map[sid][metric] = float(x_sample[f_idx])

    site_attributions = []
    for sid, s_phi in site_phi_map.items():
        entry: Dict[str, Any] = {
            "site_id": sid,
            "phi": s_phi,
        }
        # Parse locus coordinates if standard format chr:pos:repeat:unit
        parts = sid.split(":")
        if len(parts) >= 2:
            entry["chrom"] = parts[0]
            entry["pos"] = parts[1]
        if len(parts) >= 4:
            entry["repeat_info"] = f"{parts[2]}x{parts[3]}"

        entry.update(site_metrics_map.get(sid, {}))
        site_attributions.append(entry)

    # Sort descending by absolute attribution
    site_attributions.sort(key=lambda x: abs(x["phi"]), reverse=True)
    return site_attributions


def build_waterfall_figure(
    sample_id: str,
    p_msi: float,
    threshold: float,
    base_prob: float,
    site_attributions: List[Dict[str, Any]],
    top_n: int = 15,
    true_label: Optional[Union[int, str]] = None,
) -> Any:
    """
    Builds an interactive Plotly horizontal Waterfall chart showing top contributing loci.
    """
    import plotly.graph_objects as go

    if not site_attributions:
        return go.Figure()

    # Top N by absolute magnitude
    sorted_sites = sorted(site_attributions, key=lambda x: abs(x["phi"]), reverse=True)[:top_n]
    sorted_sites.reverse()  # Top contributor at the top

    labels = []
    values = []
    colors = []
    hover_texts = []

    for s in sorted_sites:
        sid = s["site_id"]
        parts = sid.split(":")
        if len(parts) >= 4:
            short_sid = f"chr{parts[0]}:{parts[1]} ({parts[2]}x{parts[3]})"
        elif len(parts) >= 2:
            short_sid = f"chr{parts[0]}:{parts[1]}"
        else:
            short_sid = sid

        feat_details = []
        for k, v in s.items():
            if k not in ("site_id", "phi", "chrom", "pos", "repeat_info") and isinstance(v, (int, float)):
                feat_details.append(f"{k}: {v:.3f}")
        feat_str = f" | {', '.join(feat_details)}" if feat_details else ""

        labels.append(short_sid)
        val = s["phi"]
        values.append(val)

        # Color: Coral/Red for MSI (+), Steel Blue for MSS (-)
        colors.append("#d9534f" if val >= 0 else "#337ab7")
        direction_str = "Pushes toward MSI" if val >= 0 else "Pushes toward MSS"
        hover_texts.append(f"<b>{short_sid}</b><br>Attribution (\u03c6): <b>{val:+.4f}</b> ({direction_str}){feat_str}")

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=labels,
            x=values,
            orientation="h",
            marker=dict(
                color=colors,
                line=dict(color="#1a1d29", width=1),
            ),
            text=[f"{v:+.4f}" for v in values],
            textposition="outside",
            hovertext=hover_texts,
            hoverinfo="text",
            showlegend=False,
        )
    )

    pred_str = "MSI" if p_msi >= threshold else "MSS"
    title_text = (
        f"<b>MSI Driver Attribution (ShapIQ) — {sample_id}</b><br>"
        f"<span style='font-size:12px; color:#9aa0a6;'>"
        f"Prediction: <b>{pred_str}</b> (p\u2098\u209b\u1d62 = {p_msi:.4f}, Cutoff = {threshold:.4f}, Baseline = {base_prob:.4f})"
        f"</span>"
    )

    fig.update_layout(
        title=title_text,
        xaxis_title="<b>Shapley Attribution (\u03c6, Contribution to MSI Probability)</b>",
        yaxis_title="",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="Inter, sans-serif", color="#e8eaed", size=12),
        margin=dict(l=150, r=60, t=70, b=50),
        xaxis=dict(
            gridcolor="#30333D",
            zerolinecolor="#e8eaed",
            zerolinewidth=1.5,
        ),
        yaxis=dict(
            gridcolor="#30333D",
        ),
        height=max(450, top_n * 32),
    )


    return fig


def export_driver_tsv(
    site_attributions: List[Dict[str, Any]],
    out_tsv: Union[str, Path],
) -> Path:
    """Exports ranked site attributions to a TSV file."""
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(site_attributions)
    df.to_csv(out_tsv, sep="\t", index=False)
    return out_tsv
