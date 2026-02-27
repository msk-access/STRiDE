"""
QC Report Generator for STRiDE (v2.1).

Generates a premium interactive HTML quality-control dashboard
from STRiDE feature TSVs.  Designed to help clinical reviewers
differentiate biological MSI from sequencing / alignment artifacts.

Design System
─────────────
- Dark theme with CVD-safe colour encoding
- Tab-based progressive disclosure (Dashboard → Explorer → Table)
- All distance metrics (L1, L2, Wasserstein) surfaced
- Tabulator.js for searchable / sortable data table

Dependencies: plotly ≥ 5.15  (Tabulator loaded from CDN)
"""

from __future__ import annotations

import json
import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

if TYPE_CHECKING:
    import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# ── Design Tokens ──────────────────────────────────────────────────────────
BG_PAGE = "#0f1117"
BG_CARD = "#1a1d29"
BG_CARD_BORDER = "#2a2d3a"
TEXT_PRIMARY = "#e8eaed"
TEXT_SECONDARY = "#9aa0a6"
ACCENT = "#4fc3f7"

CLR_NORMAL = "#5B8FF9"  # Steel Blue
CLR_TUMOR = "#F6BD16"  # Warm Amber
CLR_GOOD = "#5AD8A6"  # Sage Green
CLR_WARN = "#FF8C42"  # Warm Orange (CVD-safe separation from green)
CLR_SIG = "#9270CA"  # Lavender
CLR_GRID = "#30333D"  # Charcoal

# ── Quality Thresholds (single source of truth) ────────────────────────────
QC_THRESHOLDS = {
    "mapq": 20,
    "baseq": 20,
    "coverage": 200,
    "insert_delta": 50,
}

PLOTLY_TEMPLATE = {
    "layout": {
        "paper_bgcolor": BG_CARD,
        "plot_bgcolor": BG_CARD,
        "font": {"family": "Inter, -apple-system, sans-serif", "color": TEXT_PRIMARY, "size": 13},
        "title_font": {"size": 15, "color": TEXT_PRIMARY},
        "xaxis": {"gridcolor": CLR_GRID, "zerolinecolor": CLR_GRID},
        "yaxis": {"gridcolor": CLR_GRID, "zerolinecolor": CLR_GRID},
        "margin": {"l": 50, "r": 20, "t": 50, "b": 40},
        "hoverlabel": {
            "bgcolor": BG_CARD,
            "bordercolor": BG_CARD_BORDER,
            "font_color": TEXT_PRIMARY,
        },
    }
}


# ── Public API ─────────────────────────────────────────────────────────────
def is_qc_available() -> bool:
    """Check if optional QC dependencies are installed."""
    return PLOTLY_AVAILABLE


def parse_feature_tsv(feature_tsv: str) -> pd.DataFrame:
    """Read feature TSV and parse space-separated array columns."""
    df = pd.read_csv(feature_tsv, sep="\t")

    array_cols = [
        "normal_freqs",
        "tumor_freqs",
        "normal_norm_freqs",
        "tumor_norm_freqs",
    ]
    for col in array_cols:
        if col in df.columns:
            df[col] = (
                df[col]
                .fillna("")
                .apply(
                    lambda x: (
                        np.array([float(v) for v in str(x).split()])
                        if str(x).strip()
                        else np.array([])
                    )
                )
            )

    df["locus"] = (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(str)
        + " ("  # type: ignore[operator]
        + df["repeat_count"].astype(str)
        + "x"  # type: ignore[operator]
        + df["repeat_unit"].astype(str)  # type: ignore[operator]
        + ")"  # type: ignore[operator]
    )
    return df


def generate_report(
    feature_tsv: str,
    output_path: str,
    format: str = "html",
    prediction_result: dict | None = None,
    top_n: int = 99999,
) -> None:
    """Main entrypoint — generates the interactive HTML QC report."""
    if not is_qc_available():
        logger.error("QC dependencies missing.  pip install 'stride[qc]'")
        return

    df = parse_feature_tsv(feature_tsv)
    if len(df) == 0:
        logger.warning("Feature TSV %s is empty.", feature_tsv)
        return

    generate_html_report(df, output_path, prediction_result, top_n)


# ── Dashboard Cards ────────────────────────────────────────────────────────
def _apply_theme(fig: go.Figure) -> go.Figure:
    """Apply the dark STRiDE theme to a figure."""
    fig.update_layout(**PLOTLY_TEMPLATE["layout"])
    return fig


def _waterfall_colour_gradient(values: pd.Series, base_rgb: tuple) -> list:
    """Generate a colour gradient from charcoal (low) to base_rgb (high)."""
    mx = values.max() or 1
    return [
        f"rgba({int(48 + (base_rgb[0]-48)*(v/mx))}, "
        f"{int(51 + (base_rgb[1]-51)*(v/mx))}, "
        f"{int(61 + (base_rgb[2]-61)*(v/mx))}, 0.85)"
        for v in values
    ]


def _card_waterfalls(df: pd.DataFrame) -> go.Figure:
    """Card 1 — Triple waterfall: L1, L2, Wasserstein side-by-side."""
    metrics = [
        ("l1_distance", "L1 Distance", (79, 195, 247)),  # Cyan
        ("l2_distance", "L2 Distance", (146, 112, 202)),  # Lavender
        ("wasserstein_distance", "Wasserstein", (90, 216, 166)),  # Green
    ]

    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=[m[1] for m in metrics],
        horizontal_spacing=0.06,
    )

    for col_idx, (col, title, rgb) in enumerate(metrics, 1):
        ds = df.sort_values(col, ascending=True).reset_index(drop=True)
        colours = _waterfall_colour_gradient(ds[col], rgb)

        fig.add_trace(
            go.Bar(
                y=ds["locus"],
                x=ds[col],
                orientation="h",
                marker_color=colours,
                hovertemplate=("<b>%{y}</b><br>" f"{title}: %{{x:.4f}}<extra></extra>"),
                showlegend=False,
            ),
            row=1,
            col=col_idx,
        )
        # Show labels for top-5 highest-value loci only
        top_n_loci = min(5, len(ds))
        tick_indices = list(range(len(ds) - top_n_loci, len(ds)))
        fig.update_yaxes(
            showticklabels=True,
            tickvals=[ds["locus"].iloc[i] for i in tick_indices],
            ticktext=[ds["locus"].iloc[i] for i in tick_indices],
            tickfont={"size": 9, "color": TEXT_SECONDARY},
            title="",
            row=1,
            col=col_idx,
        )
        fig.update_xaxes(title_text=title, row=1, col=col_idx)

    fig.update_layout(height=500)
    _apply_theme(fig)
    # Ensure gridlines on ALL subplots (PLOTLY_TEMPLATE only sets axis1)
    for i in range(1, 4):
        fig.update_xaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID, row=1, col=i)
        fig.update_yaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID, row=1, col=i)
    return fig


def _card_distance_correlation(df: pd.DataFrame) -> go.Figure:
    """Card 2 — L1 vs L2 scatter, coloured by Wasserstein."""
    fig = go.Figure(
        go.Scatter(
            x=df["l1_distance"],
            y=df["l2_distance"],
            mode="markers",
            marker={
                "size": 8,
                "color": df["wasserstein_distance"],
                "colorscale": "Viridis",
                "colorbar": {"title": "Wass.", "len": 0.6, "thickness": 12},
                "opacity": 0.8,
                "line": {"width": 0.5, "color": BG_CARD_BORDER},
            },
            hovertext=df["locus"],
            hovertemplate=(
                "<b>%{hovertext}</b><br>"
                "L1: %{x:.3f}<br>L2: %{y:.3f}<br>"
                "Wass: %{marker.color:.4f}<extra></extra>"
            ),
        )
    )
    mx = max(df["l1_distance"].max(), df["l2_distance"].max()) * 1.05
    fig.add_shape(
        type="line",
        x0=0,
        y0=0,
        x1=mx,
        y1=mx,
        line={"dash": "dot", "color": TEXT_SECONDARY, "width": 1},
    )
    fig.update_layout(
        title="Distance Correlation (L1 vs L2)",
        xaxis_title="L1 Distance",
        yaxis_title="L2 Distance",
        height=500,
    )
    return _apply_theme(fig)


def _card_volcanoes(df: pd.DataFrame) -> go.Figure:
    """Card 3 — Triple volcano: L1/L2/Wasserstein vs –log10(p).

    Sites are coloured by a composite quality flag:
      * MapQ < 20, BaseQ < 20, or Coverage < 200 ⇒ flagged.
    """
    neg_log_p = -np.log10(df["p_value"].clip(lower=1e-300))

    # Build composite quality flag (uses shared QC_THRESHOLDS)
    flag_mapq = df["tumor_mapq_mean"] < QC_THRESHOLDS["mapq"]
    flag_bq = (
        (df["tumor_bq_mean"] < QC_THRESHOLDS["baseq"])
        if "tumor_bq_mean" in df.columns
        else pd.Series(False, index=df.index)
    )
    flag_cov = df["tumor_total_coverage"] < QC_THRESHOLDS["coverage"]
    flagged = flag_mapq | flag_bq | flag_cov

    # Build per-site flag description for hover
    flag_reasons: list[str] = []
    for i in df.index:
        reasons = []
        if flag_mapq[i]:
            reasons.append(f"MapQ={df.loc[i, 'tumor_mapq_mean']:.0f}")
        if flag_bq[i]:
            reasons.append(f"BaseQ={df.loc[i, 'tumor_bq_mean']:.1f}")
        if flag_cov[i]:
            reasons.append(f"Cov={df.loc[i, 'tumor_total_coverage']:.0f}")
        flag_reasons.append(", ".join(reasons) if reasons else "")
    flag_arr = np.array(flag_reasons)

    metrics = [
        ("l1_distance", "L1 Distance"),
        ("l2_distance", "L2 Distance"),
        ("wasserstein_distance", "Wasserstein"),
    ]

    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=[m[1] for m in metrics],
        horizontal_spacing=0.06,
    )

    for col_idx, (col, title) in enumerate(metrics, 1):
        # Good sites
        fig.add_trace(
            go.Scatter(
                x=df[~flagged][col],
                y=neg_log_p[~flagged],
                mode="markers",
                marker={"size": 5, "color": CLR_GOOD, "opacity": 0.7},
                name="Pass QC",
                showlegend=(col_idx == 1),
                legendgroup="good",
                hovertext=df[~flagged]["locus"],
                hovertemplate=f"<b>%{{hovertext}}</b><br>{title}: %{{x:.4f}}<br>–log₁₀(p): %{{y:.1f}}<extra></extra>",
            ),
            row=1,
            col=col_idx,
        )
        # Flagged sites
        fig.add_trace(
            go.Scatter(
                x=df[flagged][col],
                y=neg_log_p[flagged],
                mode="markers",
                marker={"size": 7, "color": CLR_WARN, "symbol": "x", "opacity": 0.9},
                name="QC Flag",
                showlegend=(col_idx == 1),
                legendgroup="flag",
                hovertext=df[flagged]["locus"],
                customdata=flag_arr[flagged],
                hovertemplate=(
                    f"<b>%{{hovertext}}</b><br>{title}: %{{x:.4f}}<br>"
                    "–log₁₀(p): %{y:.1f}<br>"
                    "⚠ %{customdata}<extra></extra>"
                ),
            ),
            row=1,
            col=col_idx,
        )
        fig.update_xaxes(title_text=title, row=1, col=col_idx)

    # Add p=0.05 threshold line to each subplot
    threshold = -np.log10(0.05)
    for i in range(1, 4):
        fig.add_hline(y=threshold, line_dash="dash", line_color=TEXT_SECONDARY, row=1, col=i)

    fig.update_yaxes(title_text="–log₁₀(p-value)", row=1, col=1)
    fig.update_layout(
        height=500,
        showlegend=True,
        legend={"x": 0.01, "y": 0.99, "bgcolor": "rgba(0,0,0,0)"},
    )
    _apply_theme(fig)
    for i in range(1, 4):
        fig.update_xaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID, row=1, col=i)
        fig.update_yaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID, row=1, col=i)
    return fig


def _card_entropy(df: pd.DataFrame) -> go.Figure:
    """Card 4 — Entropy scatter (Tumor vs Normal)."""
    fig = go.Figure(
        go.Scatter(
            x=df["normal_entropy"],
            y=df["tumor_entropy"],
            mode="markers",
            marker={
                "size": 7,
                "color": df["entropy_diff"],
                "colorscale": [[0, CLR_NORMAL], [0.5, TEXT_SECONDARY], [1, CLR_TUMOR]],
                "colorbar": {"title": "Δ Entropy", "len": 0.6, "thickness": 12},
                "opacity": 0.8,
            },
            hovertext=df["locus"],
            hovertemplate=(
                "<b>%{hovertext}</b><br>"
                "Normal ent: %{x:.2f}<br>Tumor ent: %{y:.2f}<extra></extra>"
            ),
        )
    )
    mx = max(df["normal_entropy"].max(), df["tumor_entropy"].max()) * 1.05
    fig.add_shape(
        type="line",
        x0=0,
        y0=0,
        x1=mx,
        y1=mx,
        line={"dash": "dot", "color": TEXT_SECONDARY, "width": 1},
    )
    fig.update_layout(
        title="Entropy Scatter",
        xaxis_title="Normal Entropy",
        yaxis_title="Tumor Entropy",
        height=500,
    )
    return _apply_theme(fig)


def _card_quality_metrics(df: pd.DataFrame) -> list[go.Figure]:
    """Card 5 — Three separate box+strip charts for MapQ, BaseQ, Coverage."""
    figs: list[go.Figure] = []

    def _box_strip(
        normal_vals: pd.Series,
        tumor_vals: pd.Series,
        title: str,
        y_label: str,
    ) -> go.Figure:
        """Build one box+strip figure for Normal vs Tumor."""
        fig = go.Figure()

        for vals, label, clr in [
            (normal_vals, "Normal", CLR_NORMAL),
            (tumor_vals, "Tumor", CLR_TUMOR),
        ]:
            fig.add_trace(
                go.Box(
                    y=vals,
                    name=label,
                    marker_color=clr,
                    line_color=clr,
                    boxmean="sd",
                    opacity=0.7,
                    hoverinfo="y+name",
                )
            )
            fig.add_trace(
                go.Scatter(
                    y=vals,
                    x=[label] * len(vals),
                    mode="markers",
                    marker={"size": 4, "color": clr, "opacity": 0.35},
                    showlegend=False,
                    hovertext=df["locus"],
                    hovertemplate=(
                        "<b>%{hovertext}</b><br>" f"{y_label}: %{{y:.2f}}<extra></extra>"
                    ),
                )
            )

        fig.update_layout(
            title=title,
            yaxis_title=y_label,
            height=400,
            showlegend=False,
            boxmode="group",
        )
        _apply_theme(fig)
        fig.update_xaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID)
        fig.update_yaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID)
        return fig

    # 1 ─ MapQ
    figs.append(
        _box_strip(
            df["normal_mapq_mean"],
            df["tumor_mapq_mean"],
            "Mapping Quality (MapQ)",
            "Mean MapQ",
        )
    )

    # 2 ─ Base Quality
    if "normal_bq_mean" in df.columns:
        figs.append(
            _box_strip(
                df["normal_bq_mean"],
                df["tumor_bq_mean"],
                "Base Quality",
                "Mean Base Quality",
            )
        )

    # 3 ─ Coverage (log₁₀) — custom build with real coverage on hover
    fig_cov = go.Figure()
    for raw_vals, label, clr in [
        (df["normal_total_coverage"], "Normal", CLR_NORMAL),
        (df["tumor_total_coverage"], "Tumor", CLR_TUMOR),
    ]:
        log_vals = np.log10(raw_vals.clip(lower=1))
        fig_cov.add_trace(
            go.Box(
                y=log_vals,
                name=label,
                marker_color=clr,
                line_color=clr,
                boxmean="sd",
                opacity=0.7,
                hoverinfo="y+name",
            )
        )
        fig_cov.add_trace(
            go.Scatter(
                y=log_vals,
                x=[label] * len(log_vals),
                mode="markers",
                marker={"size": 4, "color": clr, "opacity": 0.35},
                showlegend=False,
                hovertext=df["locus"],
                customdata=raw_vals.values,
                hovertemplate=(
                    "<b>%{hovertext}</b><br>"
                    "log₁₀(Cov): %{y:.2f}<br>"
                    "Coverage: %{customdata:,.0f}<extra></extra>"
                ),
            )
        )
    fig_cov.update_layout(
        title="Read Depth (log₁₀)",
        yaxis_title="log₁₀(Coverage)",
        height=400,
        showlegend=False,
        boxmode="group",
    )
    _apply_theme(fig_cov)
    fig_cov.update_xaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID)
    fig_cov.update_yaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID)
    figs.append(fig_cov)

    return figs


def _card_distance_histograms(df: pd.DataFrame) -> go.Figure:
    """Card 6 — Overlapping L1 / L2 / Wasserstein distributions."""
    fig = go.Figure()
    y_offsets = [-20, -40, -60]
    for idx, (col, name, clr) in enumerate(
        [
            ("l1_distance", "L1", ACCENT),
            ("l2_distance", "L2", CLR_SIG),
            ("wasserstein_distance", "Wasserstein", CLR_GOOD),
        ]
    ):
        fig.add_trace(
            go.Histogram(
                x=df[col],
                name=name,
                marker_color=clr,
                opacity=0.55,
                nbinsx=30,
            )
        )
        fig.add_vline(
            x=df[col].mean(),
            line_dash="dash",
            line_color=clr,
            annotation_text=f"{name} μ={df[col].mean():.3f}",
            annotation_font_color=clr,
            annotation_yshift=y_offsets[idx],
        )

    fig.update_layout(
        title="Distance Distributions",
        xaxis_title="Distance Value",
        yaxis_title="Count",
        barmode="overlay",
        height=500,
        legend={"x": 0.7, "y": 0.95, "bgcolor": "rgba(0,0,0,0)"},
    )
    return _apply_theme(fig)


def _card_insert_size(df: pd.DataFrame) -> go.Figure:
    """Card 7 — Insert size distributions (Normal vs Tumor): All, Ref, Alt."""
    panels = [
        ("normal_insert_mean_all", "tumor_insert_mean_all", "All Reads"),
        ("normal_insert_mean_ref", "tumor_insert_mean_ref", "Ref Reads"),
        ("normal_insert_mean_alt", "tumor_insert_mean_alt", "Alt Reads"),
    ]
    # Check columns exist
    available = [(n, t, lbl) for n, t, lbl in panels if n in df.columns and t in df.columns]
    if not available:
        return None  # type: ignore[return-value]

    fig = make_subplots(
        rows=1,
        cols=len(available),
        subplot_titles=[lbl for _, _, lbl in available],
        horizontal_spacing=0.06,
    )

    for col_idx, (n_col, t_col, _label) in enumerate(available, 1):
        for vals, name, clr in [
            (df[n_col], "Normal", CLR_NORMAL),
            (df[t_col], "Tumor", CLR_TUMOR),
        ]:
            fig.add_trace(
                go.Histogram(
                    x=vals,
                    name=name,
                    marker_color=clr,
                    opacity=0.55,
                    nbinsx=30,
                    showlegend=(col_idx == 1),
                    legendgroup=name,
                ),
                row=1,
                col=col_idx,
            )
            fig.add_vline(
                x=vals.mean(),
                line_dash="dash",
                line_color=clr,
                row=1,
                col=col_idx,
            )
        fig.update_xaxes(title_text="Insert Size (bp)", row=1, col=col_idx)

    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_layout(
        title="Insert Size Distribution (Normal vs Tumor)",
        barmode="overlay",
        height=450,
        showlegend=True,
        legend={"x": 0.01, "y": 0.99, "bgcolor": "rgba(0,0,0,0)"},
    )
    _apply_theme(fig)
    for i in range(1, len(available) + 1):
        fig.update_xaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID, row=1, col=i)
        fig.update_yaxes(gridcolor=CLR_GRID, zerolinecolor=CLR_GRID, row=1, col=i)
    return fig


# ── Quality Badge ──────────────────────────────────────────────────────────
def _quality_badge(row: pd.Series) -> str:
    """Return a styled HTML badge span based on quality metrics.

    Uses shared QC_THRESHOLDS — same criteria as volcano plot flags.
    Binary: ⚠ REVIEW (any threshold breached) or ✓ PASS.
    """
    mapq = row["tumor_mapq_mean"]
    bq = row.get("tumor_bq_mean", 99)
    cov = row["tumor_total_coverage"]
    insert_delta = abs(row.get("tumor_insert_mean_alt", 0) - row.get("tumor_insert_mean_ref", 0))

    flags: list[str] = []
    if mapq < QC_THRESHOLDS["mapq"]:
        flags.append(f"MapQ={mapq:.0f}")
    if bq < QC_THRESHOLDS["baseq"]:
        flags.append(f"BaseQ={bq:.1f}")
    if cov < QC_THRESHOLDS["coverage"]:
        flags.append(f"Cov={cov:.0f}")
    if insert_delta > QC_THRESHOLDS["insert_delta"]:
        flags.append(f"InsΔ={insert_delta:.0f}")

    if flags:
        reason_title = ", ".join(flags)
        return (
            f"<span style='background:rgba(232,104,74,0.25);color:{CLR_WARN};"
            "padding:2px 10px;border-radius:8px;font-size:11px;font-weight:600;"
            f"border:1px solid {CLR_WARN}' title='{reason_title}'>⚠ REVIEW</span>"
        )
    return (
        f"<span style='background:rgba(90,216,166,0.2);color:{CLR_GOOD};"
        "padding:2px 10px;border-radius:8px;font-size:11px;font-weight:600;"
        f"border:1px solid {CLR_GOOD}'>✓ PASS</span>"
    )


# ── Site Explorer ──────────────────────────────────────────────────────────
def _build_site_explorer(df: pd.DataFrame) -> go.Figure:
    """Build a single explorer figure for ALL sites, sorted by chrom + position.

    Each site gets 3 traces: Normal bars, Tumor bars, Ref line (vertical).
    All metric values are embedded in each button's args[2] so the JS
    metrics card can display them without needing multiple figures.
    """
    df_s = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    n_sites = len(df_s)
    traces_per_site = 3  # Normal, Tumor, Ref line

    # Compute global max nonzero repeat length for tight x-axis
    global_xmax = 0
    for _, _row in df_s.iterrows():
        for col in ("tumor_norm_freqs", "normal_norm_freqs"):
            nz = np.nonzero(_row[col])[0]
            if len(nz):
                global_xmax = max(global_xmax, int(nz[-1]) + 1)
    x_range = [0, global_xmax + 3]  # small padding beyond last data

    fig = go.Figure()
    buttons = []

    # Add 3 always-visible legend-only traces so legend persists across site switches
    fig.add_trace(
        go.Bar(
            x=[None],
            y=[None],
            name="Normal",
            marker_color=CLR_NORMAL,
            opacity=0.75,
            showlegend=True,
            legendgroup="Normal",
            visible=True,
        )
    )
    fig.add_trace(
        go.Bar(
            x=[None],
            y=[None],
            name="Tumor",
            marker_color=CLR_TUMOR,
            opacity=0.85,
            showlegend=True,
            legendgroup="Tumor",
            visible=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            line={"color": CLR_GOOD, "width": 2, "dash": "dash"},
            name="Ref",
            showlegend=True,
            legendgroup="Ref",
            visible=True,
        )
    )
    n_legend_traces = 3

    for i, row in df_s.iterrows():
        t_norm = row["tumor_norm_freqs"]
        n_norm = row["normal_norm_freqs"]
        t_raw = row["tumor_freqs"]
        n_raw = row["normal_freqs"]
        # Convert to Python lists so Plotly serialises as JSON arrays,
        # not binary blobs — required for JS toggle (.slice()).
        x_list = list(range(1, len(t_norm) + 1))
        n_norm_l = [round(float(v), 6) for v in n_norm]
        t_norm_l = [round(float(v), 6) for v in t_norm]
        n_raw_l = [int(v) for v in n_raw]
        t_raw_l = [int(v) for v in t_raw]
        ref_count = int(row["repeat_count"])

        visible = i == 0

        # Trace 1: Normal bars
        fig.add_trace(
            go.Bar(
                x=x_list,
                y=n_norm_l,
                name="Normal",
                marker_color=CLR_NORMAL,
                opacity=0.75,
                visible=visible,
                customdata=n_raw_l,
                hovertemplate="Normal<br>Repeat %{x}<br>Freq: %{y:.3f}<br>Reads: %{customdata}<extra></extra>",
                showlegend=False,
                legendgroup="Normal",
            )
        )
        # Trace 2: Tumor bars
        fig.add_trace(
            go.Bar(
                x=x_list,
                y=t_norm_l,
                name="Tumor",
                marker_color=CLR_TUMOR,
                opacity=0.85,
                visible=visible,
                customdata=t_raw_l,
                hovertemplate="Tumor<br>Repeat %{x}<br>Freq: %{y:.3f}<br>Reads: %{customdata}<extra></extra>",
                showlegend=False,
                legendgroup="Tumor",
            )
        )
        # Trace 3: Ref line (vertical scatter at repeat_count)
        max_y = (
            max(
                float(n_norm.max()) if len(n_norm) else 0,
                float(t_norm.max()) if len(t_norm) else 0,
                0.01,
            )
            * 1.1
        )
        fig.add_trace(
            go.Scatter(
                x=[ref_count, ref_count],
                y=[0, max_y],
                mode="lines",
                line={"color": CLR_GOOD, "width": 2, "dash": "dash"},
                name="Ref",
                visible=visible,
                showlegend=False,
                legendgroup="Ref",
                hovertemplate=f"Reference: {ref_count}x<extra></extra>",
            )
        )

        # Build visibility array: first 3 legend traces always True
        vis = [True] * n_legend_traces + [False] * (n_sites * traces_per_site)
        base = n_legend_traces + traces_per_site * int(i)  # type: ignore[call-overload]
        vis[base] = True
        vis[base + 1] = True
        vis[base + 2] = True

        buttons.append(
            {
                "label": row["locus"],
                "method": "update",
                "args": [
                    {"visible": vis},
                    {"title": ""},
                    {
                        "_locus": row["locus"],
                        "_l1": float(row["l1_distance"]),
                        "_l2": float(row["l2_distance"]),
                        "_wass": float(row["wasserstein_distance"]),
                        "_pval": float(row["p_value"]),
                        "_entropy": float(row["entropy_diff"]),
                        "_t_cov": int(row["tumor_total_coverage"]),
                        "_n_cov": int(row["normal_total_coverage"]),
                        "_t_mapq": float(row["tumor_mapq_mean"]),
                        "_t_bq": float(row.get("tumor_bq_mean", 0)),
                        "_badge": _quality_badge(row),
                    },
                ],
            }
        )

    fig.update_layout(
        updatemenus=[
            {
                "active": 0,
                "buttons": buttons,
                "visible": False,
                "x": 0.0,
                "xanchor": "left",
                "y": 1.22,
                "yanchor": "top",
                "bgcolor": BG_CARD_BORDER,
                "font": {"color": TEXT_PRIMARY, "size": 11},
                "pad": {"r": 10, "t": 10},
            }
        ],
        barmode="group",
        bargap=0.15,
        bargroupgap=0.05,
        xaxis_title="Repeat Length",
        xaxis_range=x_range,
        yaxis_title="Normalized Frequency",
        title="",
        height=520,
        legend={"x": 0.85, "y": 1.0, "bgcolor": "rgba(0,0,0,0)"},
    )

    return _apply_theme(fig)


# ── Data Table (Tabulator.js) ──────────────────────────────────────────────
def _build_data_table_json(df: pd.DataFrame) -> str:
    """Build JSON data array for Tabulator.js table."""
    ds = df.sort_values("l1_distance", ascending=False).reset_index(drop=True)
    ds["rank"] = range(1, len(ds) + 1)

    rows = []
    for _, r in ds.iterrows():
        row = {
            "rank": int(r["rank"]),
            "chrom": str(r["chrom"]),
            "start": int(r["start"]),
            "unit": str(r["repeat_unit"]),
            "count": int(r["repeat_count"]),
            "l1": round(float(r["l1_distance"]), 4),
            "l2": round(float(r["l2_distance"]), 4),
            "wass": round(float(r["wasserstein_distance"]), 5),
            "pvalue": float(r["p_value"]),
            "entropy_d": round(float(r["entropy_diff"]), 3),
            "t_cov": int(r["tumor_total_coverage"]),
            "n_cov": int(r["normal_total_coverage"]),
            "t_mapq": round(float(r["tumor_mapq_mean"]), 1),
        }
        if "tumor_bq_mean" in r.index:
            row["t_bq"] = round(float(r["tumor_bq_mean"]), 1)
        rows.append(row)

    return json.dumps(rows)


# ── HTML Assembly ──────────────────────────────────────────────────────────
_CSS = f"""
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400&display=swap');

/* ── Dual Palette (CSS Custom Properties) ──────── */
.theme-dark {{
    --bg-page: {BG_PAGE};
    --bg-card: {BG_CARD};
    --bg-card-border: {BG_CARD_BORDER};
    --text-primary: {TEXT_PRIMARY};
    --text-secondary: {TEXT_SECONDARY};
    --accent: {ACCENT};
    --grid: {CLR_GRID};
    --header-bg: #252836;
    --row-even: #1f2233;
    --row-hover: #252a40;
    --hero-grad-end: #1e2235;
    --shadow: rgba(0,0,0,0.3);
    --filter-input-bg: {BG_PAGE};
    --pill-active-bg: rgba(79,195,247,0.15);
}}
.theme-light {{
    --bg-page: #f0f2f5;
    --bg-card: #ffffff;
    --bg-card-border: #d8dce3;
    --text-primary: #1a1d29;
    --text-secondary: #5f6368;
    --accent: #1a73e8;
    --grid: #e8ebed;
    --header-bg: #f1f3f5;
    --row-even: #f8f9fa;
    --row-hover: #e8f0fe;
    --hero-grad-end: #e3e7ef;
    --shadow: rgba(0,0,0,0.06);
    --filter-input-bg: #f0f2f5;
    --pill-active-bg: rgba(26,115,232,0.10);
}}
@media print {{
    html {{ --bg-page:#fff!important; --bg-card:#fff!important; --bg-card-border:#ccc!important;
           --text-primary:#000!important; --text-secondary:#555!important; --grid:#ddd!important;
           --header-bg:#f5f5f5!important; --row-even:#fafafa!important; --shadow:none!important; }}
    .tab-bar {{ display: none !important; }}
    .tab-content {{ display: block !important; page-break-before: always; }}
    .tab-content:first-of-type {{ page-break-before: avoid; }}
    .card {{ page-break-inside: avoid; box-shadow: none !important; }}
    .card-grid {{ display: block !important; }}
    .card-grid > * {{ margin-bottom: 16px; }}
    .explorer-wrap, .table-wrap {{ page-break-inside: avoid; }}
    .tabulator-footer {{ display: none !important; }}
    #theme-toggle {{ display: none !important; }}
    body {{ padding: 12px; }}
}}

/* ── Smooth Transitions ────────────────────────── */
body, .card, .hero, .tab-btn, .table-wrap, .explorer-wrap,
.tabulator, .hero-metric .lbl {{
    transition: background-color 0.2s ease, color 0.2s ease, border-color 0.2s ease;
}}

body {{
    background: var(--bg-page);
    color: var(--text-primary);
    font-family: Inter, -apple-system, sans-serif;
    padding: 24px;
}}
.wrapper {{ max-width: 1500px; margin: 0 auto; }}

/* ── Hero Banner ───────────────────────────────── */
.hero {{
    background: linear-gradient(135deg, var(--bg-card) 0%, var(--hero-grad-end) 100%);
    border: 1px solid var(--bg-card-border);
    border-radius: 16px;
    padding: 28px 36px;
    margin-bottom: 24px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: 16px;
}}
.hero-title {{
    font-size: 22px; font-weight: 700;
    color: var(--accent);
    letter-spacing: -0.3px;
}}
.hero-subtitle {{ color: var(--text-secondary); font-size: 13px; margin-top: 4px; }}
.hero-badge {{
    display: inline-block;
    padding: 6px 18px;
    border-radius: 20px;
    font-weight: 600;
    font-size: 14px;
    letter-spacing: 0.5px;
}}
.badge-msi {{ background: rgba(230,100,70,0.25); color: {CLR_WARN}; border: 1px solid {CLR_WARN}; }}
.badge-mss {{ background: rgba(90,216,166,0.2); color: {CLR_GOOD}; border: 1px solid {CLR_GOOD}; }}
.badge-unknown {{ background: rgba(154,160,166,0.15); color: var(--text-secondary); border: 1px solid var(--text-secondary); }}
.hero-metrics {{
    display: flex; gap: 28px; flex-wrap: wrap; align-items: center;
}}
.hero-metric {{
    text-align: center;
}}
.hero-metric .val {{
    font-size: 24px; font-weight: 700;
    font-family: 'JetBrains Mono', monospace;
}}
.hero-metric .lbl {{
    font-size: 11px; color: var(--text-secondary); text-transform: uppercase; letter-spacing: 0.8px;
}}
/* ── Theme Toggle ──────────────────────────────── */
#theme-toggle {{
    background: var(--bg-card-border); border: none;
    width: 40px; height: 40px; border-radius: 50%;
    font-size: 20px; cursor: pointer;
    display: flex; align-items: center; justify-content: center;
    transition: background 0.2s, transform 0.2s;
}}
#theme-toggle:hover {{ transform: scale(1.1); }}

/* ── Tab Bar ───────────────────────────────────── */
.tab-bar {{
    display: flex; gap: 0;
    margin-bottom: 24px;
    border-bottom: 2px solid var(--bg-card-border);
}}
.tab-btn {{
    padding: 10px 24px;
    background: none; border: none;
    color: var(--text-secondary);
    font-family: Inter, sans-serif;
    font-size: 14px; font-weight: 500;
    cursor: pointer;
    border-bottom: 2px solid transparent;
    margin-bottom: -2px;
    transition: color 0.2s, border-color 0.2s;
}}
.tab-btn:hover {{ color: var(--text-primary); }}
.tab-btn.active {{
    color: var(--accent);
    border-bottom-color: var(--accent);
}}
.tab-content {{ display: none; }}
.tab-content.active {{ display: block; }}

/* ── Card Grid ─────────────────────────────────── */
.card-grid {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 20px;
}}
@media (max-width: 1100px) {{
    .card-grid {{ grid-template-columns: repeat(2, 1fr); }}
}}
.card {{
    background: var(--bg-card);
    border: 1px solid var(--bg-card-border);
    border-radius: 12px;
    overflow: hidden;
    box-shadow: 0 2px 8px var(--shadow);
}}
.card-full {{
    grid-column: 1 / -1;
}}

/* ── Section Labels ───────────────────────────── */
.section-label {{
    grid-column: 1 / -1;
    padding: 16px 4px 4px;
}}
.section-label h3 {{
    font-size: 15px; font-weight: 600; color: var(--accent);
    margin: 0 0 2px; letter-spacing: -0.2px;
}}
.section-label p {{
    font-size: 12px; color: var(--text-secondary); margin: 0;
}}

/* ── Explorer ──────────────────────────────────── */
.explorer-wrap {{
    background: var(--bg-card); border: 1px solid var(--bg-card-border);
    border-radius: 12px; padding: 16px; box-shadow: 0 2px 8px var(--shadow);
}}
.explorer-hint {{
    color: var(--text-secondary); font-size: 13px; margin-bottom: 12px;
}}

/* ── Locus Metrics Card ───────────────────────── */
.locus-metrics-card {{
    background: var(--bg-card);
    border: 1px solid var(--bg-card-border);
    border-radius: 10px;
    padding: 14px 20px;
    margin-bottom: 12px;
}}
.lm-header {{
    display: flex; align-items: center; gap: 12px;
    margin-bottom: 10px;
}}
.lm-locus {{
    font-size: 15px; font-weight: 700;
    color: var(--text-primary);
    font-family: 'JetBrains Mono', monospace;
}}
.lm-grid {{
    display: flex; flex-wrap: wrap; gap: 8px 24px;
}}
.lm-item {{
    display: flex; flex-direction: column; min-width: 70px;
}}
.lm-label {{
    font-size: 10px; font-weight: 600; letter-spacing: 0.5px;
    text-transform: uppercase;
    color: var(--text-secondary);
}}
.lm-val {{
    font-size: 14px; font-weight: 600;
    color: var(--text-primary);
    font-family: 'JetBrains Mono', monospace;
}}

/* ── Searchable Combobox ─────────────────────── */
.locus-combobox {{
    position: relative;
    display: inline-block;
    width: 420px;
    margin-bottom: 12px;
}}
.locus-combobox input {{
    width: 100%; padding: 9px 36px 9px 14px;
    border-radius: 8px;
    border: 1px solid var(--bg-card-border);
    background: var(--filter-input-bg, var(--bg-card));
    color: var(--text-primary);
    font-family: Inter, 'SF Mono', monospace;
    font-size: 13px;
    outline: none;
    transition: border-color 0.2s;
}}
.locus-combobox input:focus {{
    border-color: var(--accent);
    box-shadow: 0 0 0 2px rgba(64,169,255,0.15);
}}
.cb-arrow {{
    position: absolute; right: 12px; top: 50%;
    transform: translateY(-50%);
    pointer-events: none;
    color: var(--text-secondary);
    font-size: 14px;
}}
.cb-list {{
    display: none;
    position: absolute; z-index: 100;
    top: calc(100% + 4px); left: 0; right: 0;
    max-height: 340px; overflow-y: auto;
    background: var(--bg-card);
    border: 1px solid var(--bg-card-border);
    border-radius: 8px;
    box-shadow: 0 8px 24px var(--shadow);
    list-style: none; margin: 0; padding: 4px 0;
}}
.cb-list.open {{ display: block; }}
.cb-list li {{
    padding: 7px 14px;
    font-size: 12px; font-family: Inter, 'SF Mono', monospace;
    color: var(--text-primary);
    cursor: pointer;
    transition: background 0.12s;
    white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
}}
.cb-list li:hover,
.cb-list li.kb-focus {{
    background: var(--pill-active-bg, rgba(64,169,255,0.08));
}}
.cb-list li mark {{
    background: transparent;
    color: var(--accent); font-weight: 700;
}}
.cb-list li.selected {{
    color: var(--accent); font-weight: 600;
}}
.cb-no-match {{
    padding: 12px 14px;
    font-size: 12px; color: var(--text-secondary);
    font-style: italic;
}}
.locus-nav {{
    display: flex; align-items: center; justify-content: center;
    gap: 10px; margin-top: 8px;
}}
.locus-nav-btn {{
    width: 32px; height: 28px;
    background: transparent;
    border: 1px solid var(--bg-card-border);
    border-radius: 6px;
    color: var(--accent); font-size: 18px; line-height: 1;
    cursor: pointer;
    transition: background 0.15s, opacity 0.15s;
}}
.locus-nav-btn:hover:not(:disabled) {{
    background: var(--pill-active-bg, rgba(64,169,255,0.08));
}}
.locus-nav-btn:disabled {{
    opacity: 0.3; cursor: not-allowed;
}}
.locus-counter {{
    font-size: 12px; color: var(--text-secondary);
    font-variant-numeric: tabular-nums;
    min-width: 70px; text-align: center;
}}
.view-toggle {{
    display: flex; justify-content: center; gap: 0; margin-top: 8px;
}}
.vt-btn {{
    padding: 5px 14px; font-size: 11px;
    border: 1px solid var(--bg-card-border);
    background: transparent; color: var(--text-secondary);
    cursor: pointer; transition: all 0.15s;
    font-family: Inter, 'SF Mono', monospace;
}}
.vt-btn:first-child {{ border-radius: 6px 0 0 6px; }}
.vt-btn:last-child  {{ border-radius: 0 6px 6px 0; border-left: none; }}
.vt-btn.active {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
.vt-btn:hover:not(.active) {{ background: var(--pill-active-bg, rgba(64,169,255,0.08)); }}

/* ── Table ─────────────────────────────────────── */
.table-wrap {{
    background: var(--bg-card); border: 1px solid var(--bg-card-border);
    border-radius: 12px; padding: 16px; box-shadow: 0 2px 8px var(--shadow);
}}
.table-title {{
    font-size: 16px; font-weight: 600; color: var(--text-primary);
    margin-bottom: 12px;
}}
.btn-primary {{
    background: var(--accent); color: #fff; border: none;
    padding: 8px 18px; border-radius: 8px; cursor: pointer;
    font-size: 13px; font-weight: 600; font-family: Inter, sans-serif;
    transition: filter 0.2s, transform 0.15s;
}}
.btn-primary:hover {{ filter: brightness(1.15); transform: translateY(-1px); }}
.btn-primary:active {{ transform: translateY(0); }}

/* ── Tabulator Theme Overrides ─────────────────── */
.tabulator {{
    background-color: var(--bg-card) !important;
    border: 1px solid var(--bg-card-border) !important;
    border-radius: 8px;
    font-family: 'JetBrains Mono', monospace;
    font-size: 12px;
}}
.tabulator .tabulator-header {{
    background-color: var(--header-bg) !important;
    border-bottom: 1px solid var(--bg-card-border) !important;
}}
.tabulator .tabulator-header .tabulator-col {{
    background-color: var(--header-bg) !important;
    border-right: 1px solid var(--bg-card-border) !important;
    color: var(--text-primary) !important;
}}
.tabulator .tabulator-header .tabulator-col .tabulator-col-content {{
    padding: 8px 10px;
}}
.tabulator .tabulator-header .tabulator-col.tabulator-sortable .tabulator-col-title {{
    color: var(--text-primary) !important;
}}
.tabulator .tabulator-tableHolder .tabulator-table .tabulator-row {{
    background-color: var(--bg-card) !important;
    color: var(--text-primary) !important;
    border-bottom: 1px solid var(--grid) !important;
}}
.tabulator .tabulator-tableHolder .tabulator-table .tabulator-row:nth-child(even) {{
    background-color: var(--row-even) !important;
}}
.tabulator .tabulator-tableHolder .tabulator-table .tabulator-row:hover {{
    background-color: var(--row-hover) !important;
}}
.tabulator .tabulator-tableHolder .tabulator-table .tabulator-row .tabulator-cell {{
    border-right: 1px solid var(--grid) !important;
    padding: 6px 10px;
}}
.tabulator .tabulator-footer {{
    background-color: var(--header-bg) !important;
    border-top: 1px solid var(--bg-card-border) !important;
    color: var(--text-secondary) !important;
}}
.tabulator .tabulator-footer .tabulator-page {{
    color: var(--text-secondary) !important;
    border: 1px solid var(--bg-card-border) !important;
    background: transparent !important;
}}
.tabulator .tabulator-footer .tabulator-page.active {{
    color: var(--accent) !important;
    border-color: var(--accent) !important;
}}
.tabulator .tabulator-header .tabulator-col .tabulator-header-filter input {{
    background: var(--filter-input-bg) !important;
    border: 1px solid var(--bg-card-border) !important;
    color: var(--text-primary) !important;
    border-radius: 4px;
    padding: 4px 8px;
    font-family: Inter, sans-serif;
    font-size: 11px;
}}
.tabulator .tabulator-header .tabulator-col .tabulator-header-filter input::placeholder {{
    color: var(--text-secondary) !important;
}}
/* Header text wrapping */
.tabulator .tabulator-header .tabulator-col .tabulator-col-content .tabulator-col-title {{
    white-space: normal !important;
    word-wrap: break-word !important;
    line-height: 1.3;
    overflow: visible !important;
    text-overflow: clip !important;
}}
"""

_JS = """
// ── Theme Toggle ──────────────────────────────────
function strideToggleTheme() {
    var html = document.documentElement;
    var isDark = html.classList.contains('theme-dark');
    var next = isDark ? 'light' : 'dark';
    html.className = 'theme-' + next;
    localStorage.setItem('stride-theme', next);

    // Update toggle button
    var btn = document.getElementById('theme-toggle');
    if (btn) {
        btn.innerHTML = next === 'dark' ? '&#x1F319;' : '&#x2600;&#xFE0F;';
        btn.setAttribute('aria-checked', String(next === 'light'));
    }

    // Swap Tabulator stylesheet
    var link = document.getElementById('tabulator-theme');
    if (link) {
        link.href = next === 'dark'
            ? 'https://unpkg.com/tabulator-tables@6/dist/css/tabulator_midnight.min.css'
            : 'https://unpkg.com/tabulator-tables@6/dist/css/tabulator_simple.min.css';
    }

    // Redraw Tabulator after stylesheet swap
    setTimeout(function(){
        if (window._strideTable) window._strideTable.redraw(true);
    }, 100);

    // Re-theme all Plotly charts via relayout
    var s = getComputedStyle(document.documentElement);
    var bgCard    = s.getPropertyValue('--bg-card').trim();
    var txtPrim   = s.getPropertyValue('--text-primary').trim();
    var txtSec    = s.getPropertyValue('--text-secondary').trim();
    var gridClr   = s.getPropertyValue('--grid').trim();
    var hoverBg   = s.getPropertyValue('--bg-card').trim();
    var hoverBord = s.getPropertyValue('--bg-card-border').trim();

    var plotUpdate = {
        paper_bgcolor: bgCard,
        plot_bgcolor:  bgCard,
        'font.color':  txtPrim,
        'hoverlabel.bgcolor':     hoverBg,
        'hoverlabel.bordercolor': hoverBord,
        'hoverlabel.font.color':  txtPrim,
    };
    // Cover all subplot axes (waterfalls=3, volcanoes=3, insert=3)
    ['','2','3','4','5','6'].forEach(function(n){
        plotUpdate['xaxis'+n+'.gridcolor'] = gridClr;
        plotUpdate['yaxis'+n+'.gridcolor'] = gridClr;
        plotUpdate['xaxis'+n+'.zerolinecolor'] = gridClr;
        plotUpdate['yaxis'+n+'.zerolinecolor'] = gridClr;
    });
    document.querySelectorAll('.js-plotly-plot').forEach(function(el) {
        Plotly.relayout(el, plotUpdate);
    });
}

document.addEventListener('DOMContentLoaded', function() {
    // ── Set toggle icon on load ───────────────────────
    var isDark = document.documentElement.classList.contains('theme-dark');
    var btn = document.getElementById('theme-toggle');
    if (btn) {
        btn.innerHTML = isDark ? '&#x1F319;' : '&#x2600;&#xFE0F;';
        btn.setAttribute('aria-checked', String(!isDark));
        btn.addEventListener('click', strideToggleTheme);
    }

    // ── Tab navigation ─────────────────────────────
    document.querySelectorAll('.tab-btn').forEach(function(btn) {
        btn.addEventListener('click', function() {
            document.querySelectorAll('.tab-content').forEach(function(c) {
                c.classList.remove('active');
            });
            document.querySelectorAll('.tab-btn').forEach(function(b) {
                b.classList.remove('active');
                b.setAttribute('aria-selected', 'false');
            });
            document.getElementById(btn.dataset.target).classList.add('active');
            btn.classList.add('active');
            btn.setAttribute('aria-selected', 'true');
            window.dispatchEvent(new Event('resize'));
        });
    });

    // ── Searchable Combobox (Unified Site Explorer) ─
    var cbInput = document.getElementById('locus-search');
    var cbList  = document.getElementById('locus-listbox');
    var cbWrap  = cbInput ? cbInput.closest('.locus-combobox') : null;
    var kbIdx   = -1;
    var currentIdx = 0;

    function getPlot() {
        var el = document.querySelector('#explorer-chart .js-plotly-plot');
        if (!el || !el.layout || !el.layout.updatemenus) return null;
        return el;
    }

    function getButtons() {
        var plot = getPlot();
        return plot ? (plot.layout.updatemenus[0].buttons || []) : [];
    }

    function getMeta(btn) {
        return (btn.args && btn.args[2]) || {};
    }

    function highlightMatch(text, query) {
        if (!query) return text;
        var idx = text.toLowerCase().indexOf(query.toLowerCase());
        if (idx === -1) return text;
        return text.slice(0, idx) + '<mark>' + text.slice(idx, idx + query.length) + '</mark>' + text.slice(idx + query.length);
    }

    function updateMetricsCard(meta) {
        document.getElementById('lm-locus').textContent = meta._locus || '—';
        document.getElementById('lm-badge').innerHTML = meta._badge || '';
        document.getElementById('lm-l1').textContent = (meta._l1 != null) ? meta._l1.toFixed(3) : '—';
        document.getElementById('lm-l2').textContent = (meta._l2 != null) ? meta._l2.toFixed(3) : '—';
        document.getElementById('lm-wass').textContent = (meta._wass != null) ? meta._wass.toFixed(4) : '—';
        document.getElementById('lm-pval').textContent = (meta._pval != null) ? meta._pval.toExponential(2) : '—';
        document.getElementById('lm-entropy').textContent = (meta._entropy != null) ? meta._entropy.toFixed(3) : '—';
        document.getElementById('lm-tcov').textContent = (meta._t_cov != null) ? meta._t_cov.toLocaleString() : '—';
        document.getElementById('lm-ncov').textContent = (meta._n_cov != null) ? meta._n_cov.toLocaleString() : '—';
        document.getElementById('lm-mapq').textContent = (meta._t_mapq != null) ? meta._t_mapq.toFixed(1) : '—';
        document.getElementById('lm-bq').textContent = (meta._t_bq != null) ? meta._t_bq.toFixed(1) : '—';
    }

    function updateCounter() {
        var el = document.getElementById('locus-counter');
        if (el) el.textContent = (currentIdx + 1) + ' / ' + getButtons().length;
        var prevBtn = document.getElementById('locus-prev');
        var nextBtn = document.getElementById('locus-next');
        if (prevBtn) prevBtn.disabled = (currentIdx === 0);
        if (nextBtn) nextBtn.disabled = (currentIdx >= getButtons().length - 1);
    }

    function renderList(query) {
        if (!cbList) return;
        var btns = getButtons();
        var html = '', count = 0;
        var q = (query || '').toLowerCase();
        for (var i = 0; i < btns.length; i++) {
            var locus = getMeta(btns[i])._locus || btns[i].label;
            if (q && locus.toLowerCase().indexOf(q) === -1) continue;
            var cls = (i === currentIdx) ? ' class="selected"' : '';
            html += '<li role="option" data-idx="' + i + '"' + cls + '>' + highlightMatch(locus, q) + '</li>';
            count++;
        }
        if (count === 0) {
            html = '<li class="cb-no-match">No loci match &ldquo;' + (query||'') + '&rdquo;</li>';
        }
        cbList.innerHTML = html;
        kbIdx = -1;
        cbList.classList.add('open');
        if (cbWrap) cbWrap.setAttribute('aria-expanded', 'true');
        // Scroll selected item into view
        var sel = cbList.querySelector('li.selected');
        if (sel) setTimeout(function(){ sel.scrollIntoView({block: 'nearest'}); }, 0);
    }

    var isRaw = false;

    function applyViewMode() {
        var plot = getPlot(); if (!plot) return;
        var indices = [];
        for (var t = 3; t < plot.data.length; t++) {
            if (plot.data[t].visible === true && plot.data[t].customdata) indices.push(t);
        }
        indices.forEach(function(idx) {
            var d = plot.data[idx];
            var newY = d.customdata.slice(), newC = d.y.slice();
            Plotly.restyle(plot, { y:[newY], customdata:[newC] }, [idx]);
        });
        var yLabel = isRaw ? 'Raw Read Count' : 'Normalized Frequency';
        Plotly.relayout(plot, {'yaxis.title.text': yLabel});
        indices.forEach(function(idx) {
            var label = plot.data[idx].legendgroup;
            var ht = isRaw
                ? label+'<br>Repeat %{x}<br>Reads: %{y}<br>Freq: %{customdata:.3f}<extra></extra>'
                : label+'<br>Repeat %{x}<br>Freq: %{y:.3f}<br>Reads: %{customdata}<extra></extra>';
            Plotly.restyle(plot, {hovertemplate: ht}, [idx]);
        });
        document.getElementById('toggle-norm').classList.toggle('active', !isRaw);
        document.getElementById('toggle-raw').classList.toggle('active', isRaw);
    }

    function selectByIndex(btnIdx) {
        var plot = getPlot();
        var btns = getButtons();
        if (!plot || btnIdx < 0 || btnIdx >= btns.length) return;
        var btn = btns[btnIdx];
        Plotly.update(plot, btn.args[0], btn.args[1]);
        currentIdx = btnIdx;
        var meta = getMeta(btn);
        updateMetricsCard(meta);
        updateCounter();
        if (cbInput) cbInput.value = meta._locus || btn.label;
        closeList();
        if (isRaw) applyViewMode();
    }

    function closeList() {
        if (cbList) cbList.classList.remove('open');
        if (cbWrap) cbWrap.setAttribute('aria-expanded', 'false');
        kbIdx = -1;
    }

    function updateKbFocus() {
        if (!cbList) return;
        var items = cbList.querySelectorAll('li[data-idx]');
        items.forEach(function(li, i) {
            li.classList.toggle('kb-focus', i === kbIdx);
        });
        if (kbIdx >= 0 && kbIdx < items.length) {
            items[kbIdx].scrollIntoView({block: 'nearest'});
        }
    }

    // Populate metrics card and counter with initial locus on load
    (function() {
        var btns = getButtons();
        if (btns.length > 0) {
            updateMetricsCard(getMeta(btns[0]));
            updateCounter();
        }
    })();

    // Prev / Next navigation
    document.getElementById('locus-prev').addEventListener('click', function() {
        if (currentIdx > 0) selectByIndex(currentIdx - 1);
    });
    document.getElementById('locus-next').addEventListener('click', function() {
        var btns = getButtons();
        if (currentIdx < btns.length - 1) selectByIndex(currentIdx + 1);
    });

    // Raw / Normalized toggle
    document.getElementById('toggle-norm').addEventListener('click', function() {
        if (!isRaw) return;
        isRaw = false; applyViewMode();
    });
    document.getElementById('toggle-raw').addEventListener('click', function() {
        if (isRaw) return;
        isRaw = true; applyViewMode();
    });

    if (cbInput) {
        cbInput.addEventListener('focus', function() {
            // If input still shows the selected locus, show full unfiltered list
            var btns = getButtons();
            if (currentIdx >= 0 && currentIdx < btns.length) {
                var selLocus = getMeta(btns[currentIdx])._locus || btns[currentIdx].label;
                if (this.value === selLocus) { renderList(''); return; }
            }
            renderList(this.value);
        });
        cbInput.addEventListener('input', function() { renderList(this.value); });
        cbInput.addEventListener('keydown', function(e) {
            var items = cbList ? cbList.querySelectorAll('li[data-idx]') : [];
            if (e.key === 'ArrowDown') {
                e.preventDefault();
                if (!cbList.classList.contains('open')) { renderList(this.value); return; }
                kbIdx = Math.min(kbIdx + 1, items.length - 1);
                updateKbFocus();
            } else if (e.key === 'ArrowUp') {
                e.preventDefault();
                kbIdx = Math.max(kbIdx - 1, 0);
                updateKbFocus();
            } else if (e.key === 'Enter') {
                e.preventDefault();
                if (kbIdx >= 0 && kbIdx < items.length) {
                    selectByIndex(parseInt(items[kbIdx].dataset.idx));
                }
            } else if (e.key === 'Escape') {
                closeList();
                this.blur();
            }
        });
        document.addEventListener('click', function(e) {
            if (cbWrap && !cbWrap.contains(e.target)) closeList();
        });
    }
    if (cbList) {
        cbList.addEventListener('click', function(e) {
            var li = e.target.closest('li[data-idx]');
            if (li) selectByIndex(parseInt(li.dataset.idx));
        });
    }
});
"""


def generate_html_report(
    df: pd.DataFrame,
    output_path: str,
    prediction_result: dict | None = None,
    top_n: int = 99999,
) -> None:
    """Assemble the full interactive HTML report."""

    # ── Hero data ──────────────────────────────────────────────────────────
    n_sites = len(df)
    msi_status = "UNKNOWN"
    msi_score = 0.0
    if prediction_result:
        msi_status = prediction_result.get("msi_status", "UNKNOWN")
        msi_score = prediction_result.get("msi_score", 0.0)
    threshold = prediction_result.get("threshold", 0.50) if prediction_result else 0.50

    badge_cls = "badge-unknown"
    if msi_status == "MSI":
        badge_cls = "badge-msi"
    elif msi_status in ("NA", "MSS"):
        badge_cls = "badge-mss"

    hero_html = f"""
    <div class="hero">
        <div>
            <div class="hero-title">STRiDE MSI Quality Control Report</div>
            <div class="hero-subtitle">Interactive QC Dashboard  ·  {n_sites} microsatellite loci analysed</div>
        </div>
        <div class="hero-metrics">
            <div class="hero-metric">
                <div class="val"><span class="hero-badge {badge_cls}">{msi_status}</span></div>
                <div class="lbl">Prediction</div>
            </div>
            <div class="hero-metric">
                <div class="val" style="color:{ACCENT}">{msi_score:.3f}</div>
                <div class="lbl">MSI Score</div>
            </div>
            <div class="hero-metric">
                <div class="val">{n_sites}</div>
                <div class="lbl">Sites</div>
            </div>
            <div class="hero-metric">
                <div class="val" style="font-size:16px;color:{TEXT_SECONDARY}">{threshold:.2f}</div>
                <div class="lbl">Threshold</div>
            </div>
            <button id="theme-toggle" role="switch" aria-label="Toggle light/dark mode"
                    aria-checked="false" title="Toggle theme">&#x1F319;</button>
        </div>
    </div>
    """

    fig_waterfalls = _card_waterfalls(df)
    fig_dist_corr = _card_distance_correlation(df)
    fig_dist_hist = _card_distance_histograms(df)
    fig_volcanoes = _card_volcanoes(df)
    fig_entropy = _card_entropy(df)
    fig_quality_list = _card_quality_metrics(df)
    fig_insert_size = _card_insert_size(df)

    # Build single site explorer figure (sorted by chrom+position)
    fig_explorer = _build_site_explorer(df)

    # Data table JSON
    table_json = _build_data_table_json(df)

    _plotly_cfg = {"displayModeBar": "hover"}
    # Convert to HTML divs
    # Volcano is rendered first in the dashboard layout, so it loads PlotlyJS
    pjs = "cdn"
    volcano_html = f'<div class="card card-full">{fig_volcanoes.to_html(full_html=False, include_plotlyjs=pjs, config=_plotly_cfg)}</div>'
    # Waterfalls card (full width) — Plotly already loaded
    waterfall_html = f'<div class="card card-full">{fig_waterfalls.to_html(full_html=False, include_plotlyjs=False, config=_plotly_cfg)}</div>'

    # Half-width cards: dist corr, dist hist, entropy
    other_cards = []
    for fig in [fig_dist_corr, fig_dist_hist, fig_entropy]:
        other_cards.append(
            f'<div class="card">{fig.to_html(full_html=False, include_plotlyjs=False, config=_plotly_cfg)}</div>'
        )
    # Quality metric cards (each is a separate card)
    quality_cards = []
    for fig in fig_quality_list:
        quality_cards.append(
            f'<div class="card">{fig.to_html(full_html=False, include_plotlyjs=False, config=_plotly_cfg)}</div>'
        )
    # Insert size card (full width, optional)
    insert_html = ""
    if fig_insert_size is not None:
        insert_html = f'<div class="card card-full">{fig_insert_size.to_html(full_html=False, include_plotlyjs=False, config=_plotly_cfg)}</div>'

    explorer_html = fig_explorer.to_html(
        full_html=False, include_plotlyjs=False, config=_plotly_cfg
    )

    # Compute max values for progress bar formatters
    max_l1 = float(df["l1_distance"].max()) if len(df) else 1.0
    max_l2 = float(df["l2_distance"].max()) if len(df) else 1.0
    max_wass = float(df["wasserstein_distance"].max()) if len(df) else 1.0

    # Tabulator column definitions
    has_bq = "tumor_bq_mean" in df.columns
    tabulator_cols = f"""
    [
        {{title:"Rank", field:"rank", sorter:"number", width:70, hozAlign:"center",
         headerFilter:"number", frozen:true}},
        {{title:"Chr", field:"chrom", sorter:"string", width:70, hozAlign:"center",
         headerFilter:"input"}},
        {{title:"Start", field:"start", sorter:"number", width:120, hozAlign:"right",
         headerFilter:"input",
         formatter:function(cell){{return cell.getValue().toLocaleString();}}}},
        {{title:"Unit", field:"unit", sorter:"string", width:70, hozAlign:"center",
         headerFilter:"input"}},
        {{title:"Count", field:"count", sorter:"number", width:75, hozAlign:"center",
         headerFilter:"number"}},
        {{title:"L1", field:"l1", sorter:"number", width:140, hozAlign:"left",
         formatter:"progress", formatterParams:{{min:0, max:{max_l1:.4f},
           color:["{CLR_GOOD}", "{CLR_WARN}"], legend:true,
           legendColor:"#fff", legendAlign:"right"}},
         headerFilter:"number", tooltip:function(e,cell){{return "L1: "+cell.getValue().toFixed(4);}}}},
        {{title:"L2", field:"l2", sorter:"number", width:140, hozAlign:"left",
         formatter:"progress", formatterParams:{{min:0, max:{max_l2:.4f},
           color:["{CLR_GOOD}", "{CLR_WARN}"], legend:true,
           legendColor:"#fff", legendAlign:"right"}},
         headerFilter:"number", tooltip:function(e,cell){{return "L2: "+cell.getValue().toFixed(4);}}}},
        {{title:"Wass.", field:"wass", sorter:"number", width:140, hozAlign:"left",
         formatter:"progress", formatterParams:{{min:0, max:{max_wass:.6f},
           color:["{CLR_GOOD}", "{CLR_WARN}"], legend:true,
           legendColor:"#fff", legendAlign:"right"}},
         headerFilter:"number", tooltip:function(e,cell){{return "Wass: "+cell.getValue().toFixed(5);}}}},
        {{title:"p-value", field:"pvalue", sorter:"number", width:100, hozAlign:"right",
         formatter:function(cell){{return cell.getValue().toExponential(2);}},
         headerFilter:"number"}},
        {{title:"Ent. Δ", field:"entropy_d", sorter:"number", width:85, hozAlign:"right",
         formatter:function(cell){{var v=cell.getValue(); return v.toFixed(3);}},
         headerFilter:"number"}},
        {{title:"T_Cov", field:"t_cov", sorter:"number", width:100, hozAlign:"right",
         formatter:function(cell){{return cell.getValue().toLocaleString();}},
         headerFilter:"number"}},
        {{title:"N_Cov", field:"n_cov", sorter:"number", width:100, hozAlign:"right",
         formatter:function(cell){{return cell.getValue().toLocaleString();}},
         headerFilter:"number"}},
        {{title:"T_MapQ", field:"t_mapq", sorter:"number", width:85, hozAlign:"right",
         formatter:function(cell){{
           var v=cell.getValue();
           if(v<{QC_THRESHOLDS['mapq']}){{cell.getElement().style.color="{CLR_WARN}";cell.getElement().style.fontWeight="600";}}
           return v.toFixed(1);
         }},
         headerFilter:"number"}},
    """
    if has_bq:
        tabulator_cols += """    {title:"T_BQ", field:"t_bq", sorter:"number", width:80, hozAlign:"right",
         formatter:function(cell){return cell.getValue().toFixed(1);},
         headerFilter:"number"},
    """
    tabulator_cols += "]"

    tabulator_init_js = f"""
    var tableData = {table_json};
    window._strideTable = new Tabulator("#stride-table", {{
        data: tableData,
        layout: "fitDataFill",
        height: "calc(100vh - 260px)",
        maxHeight: "800px",
        pagination: "local",
        paginationSize: 50,
        paginationSizeSelector: [25, 50, 100, true],
        movableColumns: true,
        resizableColumns: "header",
        resizableRows: false,
        placeholder: "No Data Available",
        initialSort: [{{column:"rank", dir:"asc"}}],
        rowFormatter: function(row) {{
            var data = row.getData();
            if (data.l1 > {max_l1 * 0.7:.4f}) {{
                row.getElement().style.borderLeft = "3px solid {CLR_WARN}";
            }}
        }},
        columns: {tabulator_cols}
    }});

    // CSV Download button
    document.getElementById("download-csv-btn").addEventListener("click", function(){{
        window._strideTable.download("csv", "stride_qc_loci.csv");
    }});
    """

    # ── Assemble ───────────────────────────────────────────────────────────
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>STRiDE QC Report</title>
    <script>
    (function(){{
        var t = localStorage.getItem('stride-theme');
        if (!t) t = window.matchMedia && matchMedia('(prefers-color-scheme: light)').matches ? 'light' : 'dark';
        document.documentElement.className = 'theme-' + t;
        var sheet = t === 'dark' ? 'tabulator_midnight' : 'tabulator_simple';
        document.write('<link id="tabulator-theme" rel="stylesheet" '
          + 'href="https://unpkg.com/tabulator-tables@6/dist/css/' + sheet + '.min.css">');
    }})();
    </script>
    <style>{_CSS}</style>
</head>
<body>
<noscript>
  <div style="padding:40px;text-align:center;color:#e8eaed;background:#1a1d29;border-radius:12px;margin:24px;">
    <h2>JavaScript Required</h2>
    <p>This interactive QC report requires JavaScript to display charts and tables.</p>
  </div>
</noscript>
<div class="wrapper">
    {hero_html}

    <div class="tab-bar" role="tablist">
        <button class="tab-btn active" data-target="tab-dash"
                role="tab" aria-selected="true">Dashboard</button>
        <button class="tab-btn" data-target="tab-explorer"
                role="tab" aria-selected="false">Site Explorer ({n_sites} sites)</button>
        <button class="tab-btn" data-target="tab-table"
                role="tab" aria-selected="false">Data Table</button>
    </div>

    <div id="tab-dash" class="tab-content active" role="tabpanel">
        <div class="card-grid">
            <div class="section-label">
                <h3>Distance Significance</h3>
                <p>Sites above the dashed line (p &lt; 0.05) show significant distributional shift. Orange × markers flag low MapQ, BaseQ, or coverage.</p>
            </div>
            {volcano_html}
            <div class="section-label">
                <h3>Distance Agreement, Distribution &amp; Entropy</h3>
                <p>Correlation, histogram, and entropy views of the three distance metrics. Points near the diagonal indicate agreement; entropy shift highlights repeat-length instability.</p>
            </div>
            {other_cards[0]}
            {other_cards[1]}
            {other_cards[2]}
            <div class="section-label">
                <h3>Site Ranking</h3>
                <p>All loci ranked by distance; gradient intensity reflects magnitude.</p>
            </div>
            {waterfall_html}
            <div class="section-label">
                <h3>Sequencing Quality</h3>
                <p>MapQ, BaseQ, coverage, and insert-size distributions help distinguish biological signal from artifacts.</p>
            </div>
            {''.join(quality_cards)}
            {insert_html}
        </div>
    </div>

    <div id="tab-explorer" class="tab-content" role="tabpanel">
        <div class="explorer-wrap">
            <div class="explorer-hint">Search or browse loci to view repeat-length distributions. All distance metrics and quality stats appear below.</div>
            <div class="locus-combobox" role="combobox"
                 aria-expanded="false" aria-haspopup="listbox"
                 aria-owns="locus-listbox">
                <input id="locus-search" type="text"
                       placeholder="Search locus (e.g. chr4:55…)"
                       autocomplete="off" aria-autocomplete="list"
                       aria-controls="locus-listbox">
                <span class="cb-arrow">&#9662;</span>
                <ul id="locus-listbox" role="listbox" class="cb-list"></ul>
            </div>
            <div class="locus-nav">
                <button id="locus-prev" class="locus-nav-btn" title="Previous locus" disabled>&#8249;</button>
                <span id="locus-counter" class="locus-counter">1 / —</span>
                <button id="locus-next" class="locus-nav-btn" title="Next locus">&#8250;</button>
            </div>
            <div class="view-toggle">
                <button id="toggle-norm" class="vt-btn active">Normalized</button>
                <button id="toggle-raw" class="vt-btn">Raw Counts</button>
            </div>
            <div id="locus-metrics" class="locus-metrics-card">
                <div class="lm-header">
                    <span id="lm-locus" class="lm-locus">—</span>
                    <span id="lm-badge"></span>
                </div>
                <div class="lm-grid">
                    <div class="lm-item"><span class="lm-label">L1</span><span id="lm-l1" class="lm-val">—</span></div>
                    <div class="lm-item"><span class="lm-label">L2</span><span id="lm-l2" class="lm-val">—</span></div>
                    <div class="lm-item"><span class="lm-label">Wasserstein</span><span id="lm-wass" class="lm-val">—</span></div>
                    <div class="lm-item"><span class="lm-label">p-value</span><span id="lm-pval" class="lm-val">—</span></div>
                    <div class="lm-item"><span class="lm-label">Entropy Δ</span><span id="lm-entropy" class="lm-val">—</span></div>
                    <div class="lm-item"><span class="lm-label">T Cov</span><span id="lm-tcov" class="lm-val">—</span></div>
                    <div class="lm-item"><span class="lm-label">N Cov</span><span id="lm-ncov" class="lm-val">—</span></div>
                    <div class="lm-item"><span class="lm-label">T MapQ</span><span id="lm-mapq" class="lm-val">—</span></div>
                    <div class="lm-item"><span class="lm-label">T BQ</span><span id="lm-bq" class="lm-val">—</span></div>
                </div>
            </div>
            <div id="explorer-chart">
                {explorer_html}
            </div>
        </div>
    </div>

    <div id="tab-table" class="tab-content" role="tabpanel">
        <div class="table-wrap">
            <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:12px;">
                <div class="table-title" style="margin-bottom:0;">All Loci — Multi-Distance Ranking</div>
                <button id="download-csv-btn" class="btn-primary">⬇ Download CSV</button>
            </div>
            <div id="stride-table"></div>
        </div>
    </div>
</div>
<script src="https://unpkg.com/tabulator-tables@6/dist/js/tabulator.min.js"></script>
<script>{_JS}</script>
<script>
document.addEventListener('DOMContentLoaded', function() {{
    {tabulator_init_js}
}});
</script>
</body>
</html>"""

    with open(output_path, "w") as f:
        f.write(html)
    logger.info("Generated HTML QC report at %s", output_path)
