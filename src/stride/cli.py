"""STRiDE CLI — Typer-based command-line interface.

Provides three sub-commands:

* ``stride features``  — generate MSI feature TSVs from tumor/normal BAMs
* ``stride predict``   — predict MSI status from pre-computed feature TSVs
* ``stride run``       — end-to-end pipeline (features → prediction)
* ``stride qc``        — generate interactive QC reports from feature TSVs

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
    model: str = typer.Option(
        "svm",
        "--model",
        help="Model architecture: 'svm', 'tabpfn', 'tabpfn_access_only', or 'tabpfn_access_impact'.",
    ),
    model_joblib: Optional[str] = typer.Option(
        None,
        "--model-joblib",
        "--model-path",
        help="Path to custom trained model file (overrides default).",
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
    matched_norm_sample_barcode: Optional[str] = typer.Option(
        None,
        "--matched-norm-sample-barcode",
        help="Matched normal sample barcode (MAF-standard). Populates output column.",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable debug logging."),
) -> None:
    """Predict MSI status from pre-computed feature TSVs using SVM or TabPFN."""
    from rich.console import Console
    from rich.table import Table

    from .models import get_predictor
    from .predictor import (
        gather_samples_from_inputs,
        write_one_output_per_sample,
    )

    setup_logging(verbose=verbose)
    console = Console()

    # Delegate file gathering
    try:
        sample_ids, feature_bag = gather_samples_from_inputs(
            samples_dir=features_dir,
            sample_files=feature_files,
            list_file=list_file,
        )
    except FileNotFoundError as err:
        logger.error(
            "No feature TSV files found. Provide --features-dir, --feature-files, or --list-file."
        )
        raise typer.Exit(code=1) from err

    logger.info("Found %d sample(s) for prediction using model='%s'", len(sample_ids), model)

    # Initialize model predictor via factory
    predictor_inst = get_predictor(method=model, model_path=model_joblib)

    # Run batch evaluation
    results_df = predictor_inst.predict_batch(feature_bag)
    paths = write_one_output_per_sample(results_df, out_dir)

    # Display a Rich summary table
    table = Table(title=f"MSI Predictions (Model: {model.upper()})")
    table.add_column("Tumor Sample", style="cyan", no_wrap=True)
    table.add_column("Prediction", style="bold")
    table.add_column("Score", justify="right")

    for _, row in results_df.iterrows():
        pred_label = str(row.get("prediction", "MSS"))
        score_val = row.get("score", 0.0)
        style = "red bold" if pred_label == "MSI" else "green"
        table.add_row(
            str(row.get("sample_id", "")),
            f"[{style}]{pred_label}[/{style}]",
            f"{score_val:.6f}",
        )

    console.print(table)
    logger.info("Wrote %d prediction file(s) to %s", len(paths), out_dir)


# ---------------------------------------------------------------------------
# stride run
# ---------------------------------------------------------------------------


@app.command()
def run(
    model: str = typer.Option(
        "svm",
        "--model",
        help="Model architecture: 'svm', 'tabpfn', 'tabpfn_access_only', or 'tabpfn_access_impact'.",
    ),
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
        None,
        "--model-joblib",
        "--model-path",
        help="Path to custom trained model file (overrides default).",
    ),
    samples_list: Optional[str] = typer.Option(
        None,
        "--samples-list",
        help="CSV/TSV with sample_id, tumor_bam, normal_bam columns (batch mode).",
    ),
    sample_id: Optional[str] = typer.Option(
        None, "--sample-id", help="Sample ID / Tumor_Sample_Barcode (single-sample mode)."
    ),
    matched_norm_sample_barcode: Optional[str] = typer.Option(
        None,
        "--matched-norm-sample-barcode",
        help="Matched normal sample barcode (MAF-standard, single-sample mode). "
        "Defaults to normal BAM basename.",
    ),
    min_coverage: int = typer.Option(20, "--min-coverage", help="Minimum read coverage per site."),
    max_repeat_bins: int = typer.Option(
        100, "--max-repeat-bins", help="Maximum repeat-count bins."
    ),
    delete_features: bool = typer.Option(
        False, "--delete-features", help="Delete feature TSVs after prediction."
    ),
    generate_qc: bool = typer.Option(
        False, "--generate-qc", help="Generate an interactive QC report after prediction."
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
    resolved_model = model_joblib or (get_default_model_path() if model == "svm" else None)

    # --- Batch mode ---
    if samples_list:
        logger.info("Running batch mode from: %s (model: %s)", samples_list, model)
        results = run_end_to_end_batch(
            samples_list_path=samples_list,
            sites_file=resolved_sites,
            model_joblib=resolved_model,
            out_dir=out_dir,
            model_method=model,
            min_coverage=min_coverage,
            max_repeat_bins=max_repeat_bins,
            keep_features=not delete_features,
            generate_qc=generate_qc,
        )

        # Summary table
        table = Table(title=f"Batch Complete — {len(results)} samples")
        table.add_column("Sample", style="cyan", no_wrap=True)
        table.add_column("Prediction File")
        if generate_qc:
            table.add_column("QC Report")
            for r in results:
                table.add_row(r["sample_id"], r["prediction_txt"], r.get("qc_report", ""))
        else:
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

    logger.info("Running single-sample mode (model: %s)", model)
    res = run_end_to_end_single(
        sites_file=resolved_sites,
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        model_joblib=resolved_model,
        out_dir=out_dir,
        model_method=model,
        sample_id=sample_id,
        matched_norm_sample_barcode=matched_norm_sample_barcode,
        min_coverage=min_coverage,
        max_repeat_bins=max_repeat_bins,
        keep_features=not delete_features,
        generate_qc=generate_qc,
    )
    logger.info("Sample: %s", res["sample_id"])
    logger.info("Prediction: %s", res["prediction_txt"])
    if generate_qc:
        logger.info("QC Report: %s", res.get("qc_report"))


# ---------------------------------------------------------------------------
# stride qc
# ---------------------------------------------------------------------------


@app.command()
def qc(
    feature_tsv: str = typer.Option(
        ..., "--feature-tsv", help="Feature TSV from stride run or features."
    ),
    prediction: Optional[str] = typer.Option(
        None, "--prediction", help="Prediction output TSV to extract MSI status."
    ),
    model: Optional[str] = typer.Option(
        "tabpfn", "--model", help="Model architecture for attribution: 'tabpfn' or 'svm'."
    ),
    model_joblib: Optional[str] = typer.Option(
        None, "--model-joblib", "--model-path", help="Path to custom trained model file."
    ),
    output: str = typer.Option(
        "qc_report.html", "--output", help="Output path for the HTML report."
    ),
    explain: bool = typer.Option(
        True, "--explain/--no-explain", help="Compute ShapIQ site attributions and embed Waterfall chart."
    ),
    shapiq_budget: int = typer.Option(
        128, "--shapiq-budget", help="Evaluation budget for Shapley value sampling."
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable debug logging."),
) -> None:
    """Generate an interactive HTML QC report for a sample with ShapIQ explainability."""
    import pandas as pd
    from pathlib import Path

    from .models import get_predictor
    from .qc import generate_report, is_qc_available

    setup_logging(verbose=verbose)

    if not is_qc_available():
        logger.error("QC dependencies missing. Install using: pip install '.[qc]'")
        raise typer.Exit(code=1)

    pred_info = None
    if prediction:
        try:
            df_pred = pd.read_csv(prediction, sep="\t")
            if not df_pred.empty:
                pred_info = {
                    "msi_status": df_pred.iloc[0].get("MSI_class_predicted", df_pred.iloc[0].get("prediction", "UNKNOWN")),
                    "msi_score": float(df_pred.iloc[0].get("msi_score", df_pred.iloc[0].get("p_msi", 0.0))),
                }
        except Exception as e:
            logger.warning(f"Failed to parse prediction file: {e}")

    att_info = None
    if explain and model:
        try:
            predictor_inst = get_predictor(method=model, model_path=model_joblib)
            if hasattr(predictor_inst, "explain_sample"):
                logger.info("Computing ShapIQ locus attributions for %s...", feature_tsv)
                att_info = predictor_inst.explain_sample(feature_tsv, budget=shapiq_budget)
                if pred_info is None:
                    pred_info = {
                        "msi_status": att_info.get("prediction", "UNKNOWN"),
                        "msi_score": float(att_info.get("p_msi", 0.0)),
                    }
                # Save accompanying driver TSV in same directory as output report
                out_p = Path(output)
                driver_tsv = out_p.parent / f"{out_p.stem.replace('_qc', '')}_drivers.tsv"
                from stride.core.explainability import export_driver_tsv
                export_driver_tsv(att_info["site_attributions"], driver_tsv)
                logger.info("Saved driver loci attributions to %s", driver_tsv)
        except Exception as e:
            logger.warning("Could not compute model explainability: %s", e)

    generate_report(
        feature_tsv=feature_tsv,
        output_path=output,
        prediction_result=pred_info,
        attribution_result=att_info,
    )



# ---------------------------------------------------------------------------
# stride train
# ---------------------------------------------------------------------------


@app.command()
def train(
    method: str = typer.Option("svm", "--method", help="Model training method: 'svm' or 'tabpfn'."),
    access_msi_dir: str = typer.Option(
        ..., "--access-msi-dir", help="Directory with ACCESS MSI feature TSVs."
    ),
    access_mss_dir: str = typer.Option(
        ..., "--access-mss-dir", help="Directory with ACCESS MSS feature TSVs."
    ),
    impact_msi_dir: Optional[str] = typer.Option(
        None, "--impact-msi-dir", help="Directory with IMPACT MSI feature TSVs."
    ),
    impact_mss_dir: Optional[str] = typer.Option(
        None, "--impact-mss-dir", help="Directory with IMPACT MSS feature TSVs."
    ),
    out_dir: str = typer.Option(
        "trained_model", "--out-dir", help="Output directory for trained model."
    ),
    min_spec: float = typer.Option(0.95, "--min-spec", help="Minimum specificity constraint."),
    cv_folds: int = typer.Option(5, "--cv-folds", help="Number of cross-validation folds."),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable debug logging."),
) -> None:
    """Train an MSI prediction model (SVM or TabPFN) from feature TSV cohorts."""
    from pathlib import Path

    from .core.dataset import collect_binary_records
    from .models import get_trainer

    setup_logging(verbose=verbose)
    logger.info("Initializing STRiDe training workflow (method=%s)...", method)

    acc_msi = [Path(access_msi_dir)]
    acc_mss = [Path(access_mss_dir)]
    access_records = collect_binary_records(acc_msi, acc_mss, cohort="access")

    imp_records = None
    if impact_msi_dir and impact_mss_dir:
        imp_msi = [Path(impact_msi_dir)]
        imp_mss = [Path(impact_mss_dir)]
        imp_records = collect_binary_records(imp_msi, imp_mss, cohort="impact")

    trainer = get_trainer(method=method, cv_folds=cv_folds, min_spec=min_spec)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    trainer.train_cv_refit(
        impact_records=imp_records if imp_records is not None else access_records,
        access_records=access_records,
        test_records=access_records,
        metric_pool=["norm_l1", "norm_l2", "norm_wasserstein", "entropy_diff"],
        outdir=out_path,
    )

    logger.info("Training complete. Artifacts saved to: %s", out_path)
