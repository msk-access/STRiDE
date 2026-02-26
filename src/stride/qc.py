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

import json
import logging
from typing import Optional

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

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
CLR_WARN = "#E8684A"  # Coral
CLR_SIG = "#9270CA"  # Lavender
CLR_GRID = "#30333D"  # Charcoal

PLOTLY_TEMPLATE = {
    "layout": {
        "paper_bgcolor": BG_CARD,
        "plot_bgcolor": BG_CARD,
        "font": {"family": "Inter, -apple-system, sans-serif", "color": TEXT_PRIMARY, "size": 12},
        "title_font": {"size": 14, "color": TEXT_PRIMARY},
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
                    lambda x: np.array([float(v) for v in str(x).split()])
                    if str(x).strip()
                    else np.array([])
                )
            )

    df["locus"] = (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(str)
        + " ("
        + df["repeat_count"].astype(str)
        + "x"
        + df["repeat_unit"]
        + ")"
    )
    return df


def generate_report(
    feature_tsv: str,
    output_path: str,
    format: str = "html",
    prediction_result: Optional[dict] = None,
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
        ("l1_distance", "L1 Distance", (79, 195, 247)),    # Cyan
        ("l2_distance", "L2 Distance", (146, 112, 202)),    # Lavender
        ("wasserstein_distance", "Wasserstein", (90, 216, 166)),  # Green
    ]

    fig = make_subplots(
        rows=1, cols=3,
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
                hovertemplate=(
                    "<b>%{y}</b><br>"
                    f"{title}: %{{x:.4f}}<extra></extra>"
                ),
                showlegend=False,
            ),
            row=1, col=col_idx,
        )
        fig.update_yaxes(showticklabels=False, title="", row=1, col=col_idx)
        fig.update_xaxes(title_text=title, row=1, col=col_idx)

    fig.update_layout(height=500)
    return _apply_theme(fig)


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
    fig.add_shape(type="line", x0=0, y0=0, x1=mx, y1=mx,
                  line={"dash": "dot", "color": TEXT_SECONDARY, "width": 1})
    fig.update_layout(
        title="Distance Correlation (L1 vs L2)",
        xaxis_title="L1 Distance",
        yaxis_title="L2 Distance",
        height=500,
    )
    return _apply_theme(fig)


def _card_volcanoes(df: pd.DataFrame) -> go.Figure:
    """Card 3 — Triple volcano: L1/L2/Wasserstein vs –log10(p)."""
    neg_log_p = -np.log10(df["p_value"].clip(lower=1e-300))
    low_q = df["tumor_mapq_mean"] < 40

    metrics = [
        ("l1_distance", "L1 Distance"),
        ("l2_distance", "L2 Distance"),
        ("wasserstein_distance", "Wasserstein"),
    ]

    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=[m[1] for m in metrics],
        horizontal_spacing=0.06,
    )

    for col_idx, (col, title) in enumerate(metrics, 1):
        fig.add_trace(go.Scatter(
            x=df[~low_q][col], y=neg_log_p[~low_q],
            mode="markers",
            marker={"size": 5, "color": CLR_GOOD, "opacity": 0.7},
            name="MapQ ≥ 40",
            showlegend=(col_idx == 1),
            legendgroup="good",
            hovertext=df[~low_q]["locus"],
            hovertemplate=f"<b>%{{hovertext}}</b><br>{title}: %{{x:.4f}}<br>–log₁₀(p): %{{y:.1f}}<extra></extra>",
        ), row=1, col=col_idx)
        fig.add_trace(go.Scatter(
            x=df[low_q][col], y=neg_log_p[low_q],
            mode="markers",
            marker={"size": 7, "color": CLR_WARN, "symbol": "x", "opacity": 0.9},
            name="MapQ < 40",
            showlegend=(col_idx == 1),
            legendgroup="low",
            hovertext=df[low_q]["locus"],
            hovertemplate=f"<b>%{{hovertext}}</b><br>{title}: %{{x:.4f}}<br>–log₁₀(p): %{{y:.1f}}<extra></extra>",
        ), row=1, col=col_idx)
        fig.update_xaxes(title_text=title, row=1, col=col_idx)

    # Add p=0.05 threshold line to each subplot
    threshold = -np.log10(0.05)
    for i in range(1, 4):
        fig.add_hline(y=threshold, line_dash="dash", line_color=TEXT_SECONDARY,
                      row=1, col=i)

    fig.update_yaxes(title_text="–log₁₀(p-value)", row=1, col=1)
    fig.update_layout(
        height=500,
        legend={"x": 0.01, "y": 0.99, "bgcolor": "rgba(0,0,0,0)"},
    )
    return _apply_theme(fig)


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
    fig.add_shape(type="line", x0=0, y0=0, x1=mx, y1=mx,
                  line={"dash": "dot", "color": TEXT_SECONDARY, "width": 1})
    fig.update_layout(
        title="Entropy Scatter",
        xaxis_title="Normal Entropy",
        yaxis_title="Tumor Entropy",
        height=500,
    )
    return _apply_theme(fig)


def _card_quality_violins(df: pd.DataFrame) -> go.Figure:
    """Card 5 — Split violins for MapQ, BaseQ, Coverage."""
    fig = make_subplots(rows=1, cols=3, subplot_titles=("MapQ", "Base Quality", "log₁₀(Coverage)"))

    fig.add_trace(go.Violin(y=df["normal_mapq_mean"], name="Normal", side="negative",
                            line_color=CLR_NORMAL, meanline_visible=True, showlegend=False), row=1, col=1)
    fig.add_trace(go.Violin(y=df["tumor_mapq_mean"], name="Tumor", side="positive",
                            line_color=CLR_TUMOR, meanline_visible=True, showlegend=False), row=1, col=1)
    if "normal_bq_mean" in df.columns:
        fig.add_trace(go.Violin(y=df["normal_bq_mean"], name="Normal", side="negative",
                                line_color=CLR_NORMAL, meanline_visible=True, showlegend=False), row=1, col=2)
        fig.add_trace(go.Violin(y=df["tumor_bq_mean"], name="Tumor", side="positive",
                                line_color=CLR_TUMOR, meanline_visible=True, showlegend=False), row=1, col=2)
    fig.add_trace(go.Violin(y=np.log10(df["normal_total_coverage"].clip(lower=1)), name="Normal", side="negative",
                            line_color=CLR_NORMAL, meanline_visible=True, showlegend=False), row=1, col=3)
    fig.add_trace(go.Violin(y=np.log10(df["tumor_total_coverage"].clip(lower=1)), name="Tumor", side="positive",
                            line_color=CLR_TUMOR, meanline_visible=True, showlegend=False), row=1, col=3)

    fig.update_layout(
        title="Quality Metrics (Normal vs Tumor)",
        height=500,
        violingap=0, violinmode="overlay",
    )
    return _apply_theme(fig)


def _card_distance_histograms(df: pd.DataFrame) -> go.Figure:
    """Card 6 — Overlapping L1 / L2 / Wasserstein distributions."""
    fig = go.Figure()
    for col, name, clr in [
        ("l1_distance", "L1", ACCENT),
        ("l2_distance", "L2", CLR_SIG),
        ("wasserstein_distance", "Wasserstein", CLR_GOOD),
    ]:
        fig.add_trace(go.Histogram(
            x=df[col], name=name, marker_color=clr, opacity=0.55, nbinsx=30,
        ))
        fig.add_vline(x=df[col].mean(), line_dash="dash", line_color=clr,
                      annotation_text=f"{name} μ={df[col].mean():.3f}",
                      annotation_font_color=clr)

    fig.update_layout(
        title="Distance Distributions",
        xaxis_title="Distance Value",
        yaxis_title="Count",
        barmode="overlay",
        height=500,
        legend={"x": 0.7, "y": 0.95, "bgcolor": "rgba(0,0,0,0)"},
    )
    return _apply_theme(fig)


# ── Quality Badge ──────────────────────────────────────────────────────────
def _quality_badge(row: pd.Series) -> str:
    """Return a styled HTML badge span based on quality metrics."""
    mapq = row["tumor_mapq_mean"]
    insert_delta = abs(
        row.get("tumor_insert_mean_alt", 0) - row.get("tumor_insert_mean_ref", 0)
    )
    if mapq < 30 or insert_delta > 50:
        return (
            f"<span style='background:rgba(232,104,74,0.25);color:{CLR_WARN};"
            "padding:2px 10px;border-radius:8px;font-size:11px;font-weight:600;"
            f"border:1px solid {CLR_WARN}'>⚠ FAIL</span>"
        )
    elif mapq < 40 or insert_delta > 30:
        return (
            "<span style='background:rgba(246,189,22,0.2);color:#F6BD16;"
            "padding:2px 10px;border-radius:8px;font-size:11px;font-weight:600;"
            "border:1px solid #F6BD16'>⚠ REVIEW</span>"
        )
    else:
        return (
            f"<span style='background:rgba(90,216,166,0.2);color:{CLR_GOOD};"
            "padding:2px 10px;border-radius:8px;font-size:11px;font-weight:600;"
            f"border:1px solid {CLR_GOOD}'>✓ PASS</span>"
        )


# ── Site Explorer ──────────────────────────────────────────────────────────
def _build_site_explorer_for_metric(
    df: pd.DataFrame, metric_col: str, metric_label: str
) -> go.Figure:
    """Build an explorer figure with dropdown for ALL sites, sorted by metric_col.

    Each site gets 3 traces: Normal bars, Tumor bars, Ref line (vertical).
    The ref line correctly moves per-site because each site's traces are
    toggled via visibility.
    """
    df_s = df.sort_values(metric_col, ascending=False).reset_index(drop=True)
    n_sites = len(df_s)
    traces_per_site = 3  # Normal, Tumor, Ref line

    fig = go.Figure()
    buttons = []

    # Add 3 always-visible legend-only traces so legend persists across site switches
    fig.add_trace(go.Bar(
        x=[None], y=[None], name="Normal", marker_color=CLR_NORMAL,
        opacity=0.75, showlegend=True, legendgroup="Normal", visible=True,
    ))
    fig.add_trace(go.Bar(
        x=[None], y=[None], name="Tumor", marker_color=CLR_TUMOR,
        opacity=0.85, showlegend=True, legendgroup="Tumor", visible=True,
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="lines",
        line={"color": CLR_GOOD, "width": 2, "dash": "dash"},
        name="Ref", showlegend=True, legendgroup="Ref", visible=True,
    ))
    n_legend_traces = 3

    for i, row in df_s.iterrows():
        t_norm = row["tumor_norm_freqs"]
        n_norm = row["normal_norm_freqs"]
        t_raw = row["tumor_freqs"]
        n_raw = row["normal_freqs"]
        x = np.arange(1, len(t_norm) + 1)
        ref_count = int(row["repeat_count"])

        visible = i == 0

        # Trace 1: Normal bars
        fig.add_trace(go.Bar(
            x=x, y=n_norm, name="Normal", marker_color=CLR_NORMAL, opacity=0.75,
            visible=visible, customdata=n_raw,
            hovertemplate="Repeat %{x}<br>Freq: %{y:.3f}<br>Reads: %{customdata}<extra></extra>",
            showlegend=False,
            legendgroup="Normal",
        ))
        # Trace 2: Tumor bars
        fig.add_trace(go.Bar(
            x=x, y=t_norm, name="Tumor", marker_color=CLR_TUMOR, opacity=0.85,
            visible=visible, customdata=t_raw,
            hovertemplate="Repeat %{x}<br>Freq: %{y:.3f}<br>Reads: %{customdata}<extra></extra>",
            showlegend=False,
            legendgroup="Tumor",
        ))
        # Trace 3: Ref line (vertical scatter at repeat_count)
        max_y = max(float(n_norm.max()) if len(n_norm) else 0,
                    float(t_norm.max()) if len(t_norm) else 0, 0.01) * 1.1
        fig.add_trace(go.Scatter(
            x=[ref_count, ref_count], y=[0, max_y],
            mode="lines",
            line={"color": CLR_GOOD, "width": 2, "dash": "dash"},
            name="Ref",
            visible=visible,
            showlegend=False,
            legendgroup="Ref",
            hovertemplate=f"Reference: {ref_count}x<extra></extra>",
        ))

        # Build visibility array: first 3 legend traces always True
        vis = [True] * n_legend_traces + [False] * (n_sites * traces_per_site)
        vis[n_legend_traces + traces_per_site * i] = True
        vis[n_legend_traces + traces_per_site * i + 1] = True
        vis[n_legend_traces + traces_per_site * i + 2] = True

        badge = _quality_badge(row)
        metric_val = row[metric_col]
        fmt = ".4f" if metric_col == "wasserstein_distance" else ".3f"

        title_html = (
            f"<b>{row['locus']}</b>  {badge}<br>"
            f"<span style='font-size:11px;color:{TEXT_SECONDARY}'>"
            f"L1: {row['l1_distance']:.3f}  ·  L2: {row['l2_distance']:.3f}  ·  "
            f"Wass: {row['wasserstein_distance']:.4f}  ·  p: {row['p_value']:.2e}  ·  "
            f"Entropy Δ: {row['entropy_diff']:.3f}<br>"
            f"T_Cov: {int(row['tumor_total_coverage']):,}  ·  N_Cov: {int(row['normal_total_coverage']):,}  ·  "
            f"T_MapQ: {row['tumor_mapq_mean']:.1f}  ·  T_BQ: {row.get('tumor_bq_mean', 0):.1f}"
            f"</span>"
        )

        buttons.append({
            "label": f"Rank {i+1} · {row['locus']} — {metric_label}: {metric_val:{fmt}}",
            "method": "update",
            "args": [{"visible": vis}, {"title": title_html}],
        })

    fig.update_layout(
        updatemenus=[{
            "active": 0, "buttons": buttons,
            "x": 0.0, "xanchor": "left", "y": 1.22, "yanchor": "top",
            "bgcolor": BG_CARD_BORDER, "font": {"color": TEXT_PRIMARY, "size": 11},
            "pad": {"r": 10, "t": 10},
        }],
        barmode="group",
        xaxis_title="Repeat Length",
        yaxis_title="Normalized Frequency",
        height=520,
        legend={"x": 0.85, "y": 1.0, "bgcolor": "rgba(0,0,0,0)"},
    )

    # Set initial title
    if len(df_s) > 0:
        r = df_s.iloc[0]
        badge = _quality_badge(r)
        fig.update_layout(title=(
            f"<b>{r['locus']}</b>  {badge}<br>"
            f"<span style='font-size:11px;color:{TEXT_SECONDARY}'>"
            f"L1: {r['l1_distance']:.3f}  ·  L2: {r['l2_distance']:.3f}  ·  "
            f"Wass: {r['wasserstein_distance']:.4f}  ·  p: {r['p_value']:.2e}  ·  "
            f"Entropy Δ: {r['entropy_diff']:.3f}<br>"
            f"T_Cov: {int(r['tumor_total_coverage']):,}  ·  N_Cov: {int(r['normal_total_coverage']):,}  ·  "
            f"T_MapQ: {r['tumor_mapq_mean']:.1f}  ·  T_BQ: {r.get('tumor_bq_mean', 0):.1f}"
            f"</span>"
        ))

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
            "pvalue": f"{r['p_value']:.2e}",
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
body {{
    background: {BG_PAGE};
    color: {TEXT_PRIMARY};
    font-family: Inter, -apple-system, sans-serif;
    padding: 24px;
}}
.wrapper {{ max-width: 1500px; margin: 0 auto; }}

/* ── Hero Banner ───────────────────────────────── */
.hero {{
    background: linear-gradient(135deg, {BG_CARD} 0%, #1e2235 100%);
    border: 1px solid {BG_CARD_BORDER};
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
    color: {ACCENT};
    letter-spacing: -0.3px;
}}
.hero-subtitle {{ color: {TEXT_SECONDARY}; font-size: 13px; margin-top: 4px; }}
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
.badge-unknown {{ background: rgba(154,160,166,0.15); color: {TEXT_SECONDARY}; border: 1px solid {TEXT_SECONDARY}; }}
.hero-metrics {{
    display: flex; gap: 28px; flex-wrap: wrap;
}}
.hero-metric {{
    text-align: center;
}}
.hero-metric .val {{
    font-size: 24px; font-weight: 700;
    font-family: 'JetBrains Mono', monospace;
}}
.hero-metric .lbl {{
    font-size: 11px; color: {TEXT_SECONDARY}; text-transform: uppercase; letter-spacing: 0.8px;
}}

/* ── Tab Bar ───────────────────────────────────── */
.tab-bar {{
    display: flex; gap: 0;
    margin-bottom: 24px;
    border-bottom: 2px solid {BG_CARD_BORDER};
}}
.tab-btn {{
    padding: 10px 24px;
    background: none; border: none;
    color: {TEXT_SECONDARY};
    font-family: Inter, sans-serif;
    font-size: 14px; font-weight: 500;
    cursor: pointer;
    border-bottom: 2px solid transparent;
    margin-bottom: -2px;
    transition: color 0.2s, border-color 0.2s;
}}
.tab-btn:hover {{ color: {TEXT_PRIMARY}; }}
.tab-btn.active {{
    color: {ACCENT};
    border-bottom-color: {ACCENT};
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
    background: {BG_CARD};
    border: 1px solid {BG_CARD_BORDER};
    border-radius: 12px;
    overflow: hidden;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}}
.card-full {{
    grid-column: 1 / -1;
}}

/* ── Explorer ──────────────────────────────────── */
.explorer-wrap {{
    background: {BG_CARD}; border: 1px solid {BG_CARD_BORDER};
    border-radius: 12px; padding: 16px; box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}}
.explorer-hint {{
    color: {TEXT_SECONDARY}; font-size: 13px; margin-bottom: 12px;
}}

/* ── Metric Pills ─────────────────────────────── */
.metric-pills {{
    display: flex; gap: 8px; margin-bottom: 16px;
}}
.metric-pill {{
    padding: 6px 20px;
    border-radius: 20px;
    border: 1px solid {BG_CARD_BORDER};
    background: transparent;
    color: {TEXT_SECONDARY};
    font-family: Inter, sans-serif;
    font-size: 13px; font-weight: 500;
    cursor: pointer;
    transition: all 0.2s;
}}
.metric-pill:hover {{
    color: {TEXT_PRIMARY};
    border-color: {TEXT_SECONDARY};
}}
.metric-pill.active {{
    background: rgba(79, 195, 247, 0.15);
    color: {ACCENT};
    border-color: {ACCENT};
}}
.explorer-panel {{ display: none; }}
.explorer-panel.active {{ display: block; }}

/* ── Table ─────────────────────────────────────── */
.table-wrap {{
    background: {BG_CARD}; border: 1px solid {BG_CARD_BORDER};
    border-radius: 12px; padding: 16px; box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}}
.table-title {{
    font-size: 16px; font-weight: 600; color: {TEXT_PRIMARY};
    margin-bottom: 12px;
}}

/* Tabulator dark theme overrides */
.tabulator {{
    background-color: {BG_CARD} !important;
    border: 1px solid {BG_CARD_BORDER} !important;
    border-radius: 8px;
    font-family: 'JetBrains Mono', monospace;
    font-size: 12px;
}}
.tabulator .tabulator-header {{
    background-color: #252836 !important;
    border-bottom: 1px solid {BG_CARD_BORDER} !important;
}}
.tabulator .tabulator-header .tabulator-col {{
    background-color: #252836 !important;
    border-right: 1px solid {BG_CARD_BORDER} !important;
    color: {TEXT_PRIMARY} !important;
}}
.tabulator .tabulator-header .tabulator-col .tabulator-col-content {{
    padding: 8px 10px;
}}
.tabulator .tabulator-header .tabulator-col.tabulator-sortable .tabulator-col-title {{
    color: {TEXT_PRIMARY} !important;
}}
.tabulator .tabulator-tableHolder .tabulator-table .tabulator-row {{
    background-color: {BG_CARD} !important;
    color: {TEXT_PRIMARY} !important;
    border-bottom: 1px solid {CLR_GRID} !important;
}}
.tabulator .tabulator-tableHolder .tabulator-table .tabulator-row:nth-child(even) {{
    background-color: #1f2233 !important;
}}
.tabulator .tabulator-tableHolder .tabulator-table .tabulator-row:hover {{
    background-color: #252a40 !important;
}}
.tabulator .tabulator-tableHolder .tabulator-table .tabulator-row .tabulator-cell {{
    border-right: 1px solid {CLR_GRID} !important;
    padding: 6px 10px;
}}
.tabulator .tabulator-footer {{
    background-color: #252836 !important;
    border-top: 1px solid {BG_CARD_BORDER} !important;
    color: {TEXT_SECONDARY} !important;
}}
.tabulator .tabulator-footer .tabulator-page {{
    color: {TEXT_SECONDARY} !important;
    border: 1px solid {BG_CARD_BORDER} !important;
    background: transparent !important;
}}
.tabulator .tabulator-footer .tabulator-page.active {{
    color: {ACCENT} !important;
    border-color: {ACCENT} !important;
}}
.tabulator .tabulator-header .tabulator-col .tabulator-header-filter input {{
    background: {BG_PAGE} !important;
    border: 1px solid {BG_CARD_BORDER} !important;
    color: {TEXT_PRIMARY} !important;
    border-radius: 4px;
    padding: 4px 8px;
    font-family: Inter, sans-serif;
    font-size: 11px;
}}
.tabulator .tabulator-header .tabulator-col .tabulator-header-filter input::placeholder {{
    color: {TEXT_SECONDARY} !important;
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
document.addEventListener('DOMContentLoaded', function() {
    // ── Tab navigation ─────────────────────────────
    document.querySelectorAll('.tab-btn').forEach(function(btn) {
        btn.addEventListener('click', function() {
            document.querySelectorAll('.tab-content').forEach(function(c) {
                c.classList.remove('active');
            });
            document.querySelectorAll('.tab-btn').forEach(function(b) {
                b.classList.remove('active');
            });
            document.getElementById(btn.dataset.target).classList.add('active');
            btn.classList.add('active');
            window.dispatchEvent(new Event('resize'));
        });
    });

    // ── Metric pill switching (Site Explorer) ──────
    document.querySelectorAll('.metric-pill').forEach(function(pill) {
        pill.addEventListener('click', function() {
            document.querySelectorAll('.metric-pill').forEach(function(p) {
                p.classList.remove('active');
            });
            document.querySelectorAll('.explorer-panel').forEach(function(ep) {
                ep.classList.remove('active');
            });
            pill.classList.add('active');
            document.getElementById(pill.dataset.panel).classList.add('active');
            window.dispatchEvent(new Event('resize'));
        });
    });
});
"""


def generate_html_report(
    df: pd.DataFrame,
    output_path: str,
    prediction_result: Optional[dict] = None,
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
        </div>
    </div>
    """

    # ── Build figures ──────────────────────────────────────────────────────
    fig_waterfalls = _card_waterfalls(df)
    fig_dist_corr = _card_distance_correlation(df)
    fig_dist_hist = _card_distance_histograms(df)
    fig_volcanoes = _card_volcanoes(df)
    fig_entropy = _card_entropy(df)
    fig_violins = _card_quality_violins(df)

    # Build 3 site explorer figures (one per metric)
    fig_explorer_l1 = _build_site_explorer_for_metric(df, "l1_distance", "L1")
    fig_explorer_l2 = _build_site_explorer_for_metric(df, "l2_distance", "L2")
    fig_explorer_wass = _build_site_explorer_for_metric(df, "wasserstein_distance", "Wass")

    # Data table JSON
    table_json = _build_data_table_json(df)

    # Convert to HTML divs
    pjs = "cdn"  # first figure loads Plotly JS from CDN
    # Waterfalls card (full width)
    waterfall_html = f'<div class="card card-full">{fig_waterfalls.to_html(full_html=False, include_plotlyjs=pjs)}</div>'
    # Volcanoes card (full width)
    volcano_html = f'<div class="card card-full">{fig_volcanoes.to_html(full_html=False, include_plotlyjs=False)}</div>'

    # Order: dist corr, dist hist, entropy, violins
    other_cards = []
    for fig in [fig_dist_corr, fig_dist_hist, fig_entropy, fig_violins]:
        other_cards.append(
            f'<div class="card">{fig.to_html(full_html=False, include_plotlyjs=False)}</div>'
        )

    explorer_l1_html = fig_explorer_l1.to_html(full_html=False, include_plotlyjs=False)
    explorer_l2_html = fig_explorer_l2.to_html(full_html=False, include_plotlyjs=False)
    explorer_wass_html = fig_explorer_wass.to_html(full_html=False, include_plotlyjs=False)

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
        {{title:"p-value", field:"pvalue", sorter:"string", width:100, hozAlign:"right",
         headerFilter:"input"}},
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
           if(v<40){{cell.getElement().style.color="{CLR_WARN}";cell.getElement().style.fontWeight="600";}}
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
    var table = new Tabulator("#stride-table", {{
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
        table.download("csv", "stride_qc_loci.csv");
    }});
    """

    # ── Assemble ───────────────────────────────────────────────────────────
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>STRiDE QC Report</title>
    <link href="https://unpkg.com/tabulator-tables@6/dist/css/tabulator_midnight.min.css" rel="stylesheet">
    <style>{_CSS}</style>
</head>
<body>
<div class="wrapper">
    {hero_html}

    <div class="tab-bar">
        <button class="tab-btn active" data-target="tab-dash">Dashboard</button>
        <button class="tab-btn" data-target="tab-explorer">Site Explorer ({n_sites} sites)</button>
        <button class="tab-btn" data-target="tab-table">Data Table</button>
    </div>

    <div id="tab-dash" class="tab-content active">
        <div class="card-grid">
            {waterfall_html}
            {''.join(other_cards[:2])}
            {volcano_html}
            {''.join(other_cards[2:])}
        </div>
    </div>

    <div id="tab-explorer" class="tab-content">
        <div class="explorer-wrap">
            <div class="explorer-hint">Select a locus from the dropdown to visualise its repeat-length distribution shift. Use the metric pills to re-rank sites.</div>
            <div class="metric-pills">
                <button class="metric-pill active" data-panel="panel-l1">L1 Distance</button>
                <button class="metric-pill" data-panel="panel-l2">L2 Distance</button>
                <button class="metric-pill" data-panel="panel-wass">Wasserstein</button>
            </div>
            <div id="panel-l1" class="explorer-panel active">
                {explorer_l1_html}
            </div>
            <div id="panel-l2" class="explorer-panel">
                {explorer_l2_html}
            </div>
            <div id="panel-wass" class="explorer-panel">
                {explorer_wass_html}
            </div>
        </div>
    </div>

    <div id="tab-table" class="tab-content">
        <div class="table-wrap">
            <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:12px;">
                <div class="table-title" style="margin-bottom:0;">All Loci — Multi-Distance Ranking</div>
                <button id="download-csv-btn" style="
                    background:{ACCENT};color:#fff;border:none;padding:8px 18px;
                    border-radius:8px;cursor:pointer;font-size:13px;font-weight:600;
                    font-family:Inter,sans-serif;
                ">⬇ Download CSV</button>
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
