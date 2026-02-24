"""STRiDE CLI — Typer-based command-line interface.

Provides three sub-commands:

* ``stride features``  — generate MSI feature TSVs from tumor/normal BAMs
* ``stride predict``   — predict MSI/MSS from pre-computed feature TSVs
* ``stride run``       — end-to-end pipeline (features → prediction)

All parameters are explicit ``--options``; no positional arguments.
"""

import logging
from typing import Optional

import typer

from . import __version__
from .models import get_default_model_path, get_default_sites_path
from .utils import setup_logging

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# App / version callback
# ---------------------------------------------------------------------------

app = typer.Typer(
    name="stride",
    help="STRiDE — MSI prediction pipeline for MSK-ACCESS.",
    add_completion=False,
    no_args_is_help=True,
)


def _version_callback(value: bool) -> None:
    """Print the version string and exit."""
    if value:
        typer.echo(f"stride {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    """STRiDE — Microsatellite Instability prediction for MSK-ACCESS."""


# ---------------------------------------------------------------------------
# stride features
# ---------------------------------------------------------------------------

@app.command()
def features(
    tumor_bam: str = typer.Option(..., "--tumor-bam", help="Path to tumor BAM file."),
    normal_bam: str = typer.Option(..., "--normal-bam", help="Path to normal BAM file."),
    out_dir: str = typer.Option("out", "--out-dir", help="Base output directory."),
    site_list: Optional[str] = typer.Option(
        None, "--site-list", help="MSI site list TSV. Default: bundled 170-site list."
    ),
    output_tsv: Optional[str] = typer.Option(
        None, "--output-tsv", help="Explicit output TSV path (overrides --out-dir)."
    ),
    sample_id: Optional[str] = typer.Option(
        None, "--sample-id", help="Sample ID for output filename."
    ),
    min_coverage: int = typer.Option(20, "--min-coverage", help="Minimum read coverage per site."),
    max_repeat_bins: int = typer.Option(
        100, "--max-repeat-bins", help="Maximum repeat-count bins."
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable debug logging."),
) -> None:
    """Generate MSI feature TSV from paired tumor/normal BAMs."""
    import os

    from .pipeline import run_feature_generation
    from .utils import safe_name, strip_ext

    setup_logging(verbose=verbose)

    # Resolve site list — use bundled default if not provided
    resolved_sites = site_list or get_default_sites_path()

    # Determine output path
    if output_tsv:
        out_tsv = output_tsv
    else:
        sid = sample_id.strip() if sample_id else strip_ext(tumor_bam)
        features_dir = os.path.join(out_dir, "features")
        os.makedirs(features_dir, exist_ok=True)
        out_tsv = os.path.join(features_dir, f"msi_features_{safe_name(sid)}.tsv")

    run_feature_generation(
        sites_file=resolved_sites,
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        out_features_tsv=out_tsv,
        min_coverage=min_coverage,
        max_repeat_bins=max_repeat_bins,
    )

    logger.info("Feature TSV written to: %s", out_tsv)


# ---------------------------------------------------------------------------
# stride predict
# ---------------------------------------------------------------------------

@app.command()
def predict(
    model_joblib: Optional[str] = typer.Option(
        None, "--model-joblib", help="Trained model .joblib. Default: bundled SGD model."
    ),
    features_dir: Optional[str] = typer.Option(
        None, "--features-dir", help="Directory containing feature TSVs (recursive search)."
    ),
    feature_files: Optional[list[str]] = typer.Option(
        None, "--feature-files", help="One or more feature TSV files."
    ),
    list_file: Optional[str] = typer.Option(
        None, "--list-file", help="Text file with one feature TSV path per line."
    ),
    out_dir: str = typer.Option(..., "--out-dir", help="Output directory for prediction files."),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable debug logging."),
) -> None:
    """Predict MSI/MSS from pre-computed feature TSVs."""
    from rich.console import Console
    from rich.table import Table

    from .predictor import (
        gather_samples_from_inputs,
        load_model,
        get_expected_features,
        build_matrix,
        get_scores,
        write_one_output_per_sample,
    )

    setup_logging(verbose=verbose)
    console = Console()

    # Resolve model — use bundled default if not provided
    resolved_model = model_joblib or get_default_model_path()

    # Delegate file gathering to predictor.gather_samples_from_inputs()
    try:
        sample_ids, feature_bag = gather_samples_from_inputs(
            samples_dir=features_dir,
            sample_files=feature_files,
            list_file=list_file,
        )
    except FileNotFoundError:
        logger.error("No feature TSV files found. Provide --features-dir, --feature-files, or --list-file.")
        raise typer.Exit(code=1)

    logger.info("Found %d sample(s)", len(sample_ids))

    # Load model and predict
    import numpy as np
    import pandas as pd

    model = load_model(resolved_model)
    expected = get_expected_features(model)
    X = build_matrix(sample_ids, feature_bag, expected)
    y_pred = np.asarray(model.predict(X)).astype(int)
    scores = np.asarray(get_scores(model, X), dtype=float)

    df_preds = pd.DataFrame({
        "Sample_ID": sample_ids,
        "MSI_class_predicted": y_pred,
        "msi_score": np.round(scores, 6),
    })
    paths = write_one_output_per_sample(df_preds, out_dir)

    # Display a Rich summary table
    table = Table(title="MSI Predictions")
    table.add_column("Sample", style="cyan", no_wrap=True)
    table.add_column("Prediction", style="bold")
    table.add_column("Score", justify="right")

    for _, row in df_preds.iterrows():
        pred_label = "MSI-H" if int(row["MSI_class_predicted"]) == 1 else "MSS"
        style = "red bold" if pred_label == "MSI-H" else "green"
        table.add_row(str(row["Sample_ID"]), f"[{style}]{pred_label}[/{style}]", f"{row['msi_score']:.6f}")

    console.print(table)
    logger.info("Wrote %d prediction file(s) to %s", len(paths), out_dir)


# ---------------------------------------------------------------------------
# stride run
# ---------------------------------------------------------------------------

@app.command()
def run(
    tumor_bam: Optional[str] = typer.Option(
        None, "--tumor-bam", help="Tumor BAM (single-sample mode)."
    ),
    normal_bam: Optional[str] = typer.Option(
        None, "--normal-bam", help="Normal BAM (single-sample mode)."
    ),
    out_dir: str = typer.Option(..., "--out-dir", help="Output directory."),
    site_list: Optional[str] = typer.Option(
        None, "--site-list", help="MSI site list TSV. Default: bundled 170-site list."
    ),
    model_joblib: Optional[str] = typer.Option(
        None, "--model-joblib", help="Trained model .joblib. Default: bundled SGD model."
    ),
    samples_list: Optional[str] = typer.Option(
        None,
        "--samples-list",
        help="CSV/TSV with sample_id, tumor_bam, normal_bam columns (batch mode).",
    ),
    sample_id: Optional[str] = typer.Option(
        None, "--sample-id", help="Sample ID (single-sample mode)."
    ),
    min_coverage: int = typer.Option(20, "--min-coverage", help="Minimum read coverage per site."),
    max_repeat_bins: int = typer.Option(
        100, "--max-repeat-bins", help="Maximum repeat-count bins."
    ),
    delete_features: bool = typer.Option(
        False, "--delete-features", help="Delete feature TSVs after prediction."
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable debug logging."),
) -> None:
    """End-to-end MSI pipeline: BAMs → features → prediction."""
    from rich.console import Console
    from rich.table import Table

    from .pipeline import run_end_to_end_batch, run_end_to_end_single

    setup_logging(verbose=verbose)
    console = Console()

    # Resolve defaults
    resolved_sites = site_list or get_default_sites_path()
    resolved_model = model_joblib or get_default_model_path()

    # --- Batch mode ---
    if samples_list:
        logger.info("Running batch mode from: %s", samples_list)
        results = run_end_to_end_batch(
            samples_list_path=samples_list,
            sites_file=resolved_sites,
            model_joblib=resolved_model,
            out_dir=out_dir,
            min_coverage=min_coverage,
            max_repeat_bins=max_repeat_bins,
            keep_features=not delete_features,
        )

        # Summary table
        table = Table(title=f"Batch Complete — {len(results)} samples")
        table.add_column("Sample", style="cyan", no_wrap=True)
        table.add_column("Prediction File")
        for r in results:
            table.add_row(r["sample_id"], r["prediction_txt"])
        console.print(table)
        return

    # --- Single-sample mode ---
    if not tumor_bam or not normal_bam:
        logger.error(
            "Provide --samples-list (batch) OR both --tumor-bam and --normal-bam (single)."
        )
        raise typer.Exit(code=1)

    logger.info("Running single-sample mode")
    res = run_end_to_end_single(
        sites_file=resolved_sites,
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        model_joblib=resolved_model,
        out_dir=out_dir,
        sample_id=sample_id,
        min_coverage=min_coverage,
        max_repeat_bins=max_repeat_bins,
        keep_features=not delete_features,
    )
    logger.info("Sample: %s", res["sample_id"])
    logger.info("Prediction: %s", res["prediction_txt"])
