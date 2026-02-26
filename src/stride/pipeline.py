"""Pipeline orchestration for STRiDE.

Provides three run modes:

* :func:`run_feature_generation` — single-sample feature extraction
* :func:`run_end_to_end_single` — single-sample features → prediction
* :func:`run_end_to_end_batch`  — batch features → prediction

All functions use the context-manager form of :class:`MSIProfileGenerator`
to ensure BAM file handles are released promptly.
"""

import logging
import os
from typing import Optional

from .feature_generator import MSIProfileGenerator
from .predictor import predict_from_feature_tsvs, write_one_output_per_sample
from .qc import generate_report, is_qc_available
from .utils import read_samples_list, safe_name, strip_ext

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Feature generation (single sample)
# ---------------------------------------------------------------------------


def run_feature_generation(
    sites_file: str,
    tumor_bam: str,
    normal_bam: str,
    out_features_tsv: str,
    min_coverage: int = 20,
    max_repeat_bins: int = 100,
) -> str:
    """Generate a per-site feature TSV for one tumor/normal pair.

    Returns the path to the written TSV.
    """
    os.makedirs(os.path.dirname(out_features_tsv) or ".", exist_ok=True)

    with MSIProfileGenerator(
        sites_file=sites_file,
        tumor_bam_path=tumor_bam,
        normal_bam_path=normal_bam,
        min_coverage=min_coverage,
        max_repeat_bins=max_repeat_bins,
    ) as gen:
        gen.run(out_features_tsv)

    return out_features_tsv


# ---------------------------------------------------------------------------
# End-to-end: single sample
# ---------------------------------------------------------------------------


def run_end_to_end_single(
    sites_file: str,
    tumor_bam: str,
    normal_bam: str,
    model_joblib: str,
    out_dir: str,
    sample_id: Optional[str] = None,
    matched_norm_sample_barcode: Optional[str] = None,
    min_coverage: int = 20,
    max_repeat_bins: int = 100,
    keep_features: bool = True,
    generate_qc: bool = False,
) -> dict[str, str]:
    """Run the full pipeline for a single sample.

    Parameters
    ----------
    matched_norm_sample_barcode : str | None
        Explicit matched-normal barcode for the output.  When *None*,
        derived automatically from *normal_bam* via :func:`strip_ext`.

    Returns a dict with keys ``sample_id``, ``features_tsv``, and
    ``prediction_txt``.
    """
    os.makedirs(out_dir, exist_ok=True)

    sid = sample_id.strip() if sample_id else strip_ext(tumor_bam)
    normal_bc = matched_norm_sample_barcode or strip_ext(normal_bam)
    logger.info("Processing sample: %s (matched normal: %s)", sid, normal_bc)

    features_dir = os.path.join(out_dir, "features")
    preds_dir = os.path.join(out_dir, "predictions")
    os.makedirs(features_dir, exist_ok=True)
    os.makedirs(preds_dir, exist_ok=True)

    feat_tsv = os.path.join(features_dir, f"msi_features_{safe_name(sid)}.tsv")

    run_feature_generation(
        sites_file,
        tumor_bam,
        normal_bam,
        feat_tsv,
        min_coverage=min_coverage,
        max_repeat_bins=max_repeat_bins,
    )

    df_preds = predict_from_feature_tsvs(model_joblib, [feat_tsv], normal_barcodes=[normal_bc])
    out_paths = write_one_output_per_sample(df_preds, preds_dir)

    # Optional QC Generation
    qc_path = ""
    if generate_qc:
        if not is_qc_available():
            logger.warning("QC requested but dependencies are missing. Skipping QC generation.")
        else:
            qc_dir = os.path.join(out_dir, "qc")
            os.makedirs(qc_dir, exist_ok=True)
            qc_path = os.path.join(qc_dir, f"{safe_name(sid)}_qc.html")

            # Extract basic prediction info
            pred_info = None
            if len(df_preds) > 0:
                pred_info = {
                    "msi_status": df_preds.iloc[0]["MSI_class_predicted"],
                    "msi_score": df_preds.iloc[0]["msi_score"]
                }

            logger.info("Generating QC report: %s", qc_path)
            generate_report(feat_tsv, qc_path, prediction_result=pred_info)

    if not keep_features:
        try:
            os.remove(feat_tsv)
            logger.debug("Removed intermediate feature TSV: %s", feat_tsv)
        except OSError:
            pass

    return {
        "sample_id": sid,
        "features_tsv": feat_tsv,
        "prediction_txt": out_paths[0],
        "qc_report": qc_path
    }


# ---------------------------------------------------------------------------
# End-to-end: batch
# ---------------------------------------------------------------------------


def run_end_to_end_batch(
    samples_list_path: str,
    sites_file: str,
    model_joblib: str,
    out_dir: str,
    min_coverage: int = 20,
    max_repeat_bins: int = 100,
    keep_features: bool = True,
    generate_qc: bool = False,
) -> list[dict[str, str]]:
    """Run the full pipeline for every sample in a manifest file.

    Returns a list of dicts (one per sample), each containing
    ``sample_id``, ``features_tsv``, and ``prediction_txt``.
    """
    os.makedirs(out_dir, exist_ok=True)

    samples: list[dict[str, str]] = read_samples_list(samples_list_path)
    logger.info("Loaded %d samples from %s", len(samples), samples_list_path)

    features_dir = os.path.join(out_dir, "features")
    preds_dir = os.path.join(out_dir, "predictions")
    os.makedirs(features_dir, exist_ok=True)
    os.makedirs(preds_dir, exist_ok=True)

    # 1) Generate feature TSV per sample
    feature_tsvs: list[str] = []
    sample_ids: list[str] = []
    normal_barcodes: list[str] = []
    for i, s in enumerate(samples, 1):
        sid = s["sample_id"]
        logger.info("[%d/%d] Generating features for %s", i, len(samples), sid)

        feat_tsv = os.path.join(features_dir, f"msi_features_{safe_name(sid)}.tsv")
        run_feature_generation(
            sites_file,
            s["tumor_bam"],
            s["normal_bam"],
            feat_tsv,
            min_coverage=min_coverage,
            max_repeat_bins=max_repeat_bins,
        )
        feature_tsvs.append(feat_tsv)
        sample_ids.append(sid)
        # Use explicit barcode from manifest if available, otherwise derive
        normal_barcodes.append(s.get("matched_norm_sample_barcode") or strip_ext(s["normal_bam"]))

    # 2) Predict all (batch) then write one output per sample
    logger.info("Running batch prediction for %d samples", len(sample_ids))
    df_preds = predict_from_feature_tsvs(
        model_joblib, feature_tsvs, normal_barcodes=normal_barcodes
    )
    out_paths = write_one_output_per_sample(df_preds, preds_dir)

    # 3) Optional QC Generation
    qc_paths: list[str] = [""] * len(sample_ids)
    if generate_qc:
        if not is_qc_available():
            logger.warning("QC requested but dependencies are missing. Skipping QC generation.")
        else:
            qc_dir = os.path.join(out_dir, "qc")
            os.makedirs(qc_dir, exist_ok=True)
            logger.info("Generating QC reports for %d samples", len(sample_ids))
            for idx, (sid, feat) in enumerate(zip(sample_ids, feature_tsvs)):
                qc_path = os.path.join(qc_dir, f"{safe_name(sid)}_qc.html")

                pred_info = None
                # Match the row in df_preds
                row_match = df_preds[df_preds["Tumor_Sample_Barcode"] == sid]
                if not row_match.empty:
                    pred_info = {
                        "msi_status": row_match.iloc[0]["MSI_class_predicted"],
                        "msi_score": row_match.iloc[0]["msi_score"]
                    }

                generate_report(feat, qc_path, prediction_result=pred_info)
                qc_paths[idx] = qc_path

    # 4) Optional cleanup of intermediate feature TSVs
    if not keep_features:
        for fp in feature_tsvs:
            try:
                os.remove(fp)
            except OSError:
                pass
        logger.debug("Removed %d intermediate feature TSV(s)", len(feature_tsvs))

    # 5) Map outputs back to sample_id
    results: list[dict[str, str]] = []
    pred_map = {os.path.basename(p).replace("_msi.txt", ""): p for p in out_paths}
    for i, sid in enumerate(sample_ids):
        results.append(
            {
                "sample_id": sid,
                "features_tsv": feature_tsvs[i],
                "prediction_txt": pred_map.get(safe_name(sid), ""),
                "qc_report": qc_paths[i],
            }
        )

    logger.info("Batch complete — %d samples processed", len(results))
    return results
