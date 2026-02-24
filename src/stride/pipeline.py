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
    min_coverage: int = 20,
    max_repeat_bins: int = 100,
    keep_features: bool = True,
) -> dict[str, str]:
    """Run the full pipeline for a single sample.

    Returns a dict with keys ``sample_id``, ``features_tsv``, and
    ``prediction_txt``.
    """
    os.makedirs(out_dir, exist_ok=True)

    sid = sample_id.strip() if sample_id else strip_ext(tumor_bam)
    logger.info("Processing sample: %s", sid)

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

    df_preds = predict_from_feature_tsvs(model_joblib, [feat_tsv])
    out_paths = write_one_output_per_sample(df_preds, preds_dir)

    if not keep_features:
        try:
            os.remove(feat_tsv)
            logger.debug("Removed intermediate feature TSV: %s", feat_tsv)
        except OSError:
            pass

    return {"sample_id": sid, "features_tsv": feat_tsv, "prediction_txt": out_paths[0]}


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

    # 2) Predict all (batch) then write one output per sample
    logger.info("Running batch prediction for %d samples", len(sample_ids))
    df_preds = predict_from_feature_tsvs(model_joblib, feature_tsvs)
    out_paths = write_one_output_per_sample(df_preds, preds_dir)

    # 3) Optional cleanup of intermediate feature TSVs
    if not keep_features:
        for fp in feature_tsvs:
            try:
                os.remove(fp)
            except OSError:
                pass
        logger.debug("Removed %d intermediate feature TSV(s)", len(feature_tsvs))

    # 4) Map outputs back to sample_id
    results: list[dict[str, str]] = []
    pred_map = {os.path.basename(p).replace("_msi.txt", ""): p for p in out_paths}
    for sid, feat in zip(sample_ids, feature_tsvs):
        results.append(
            {
                "sample_id": sid,
                "features_tsv": feat,
                "prediction_txt": pred_map.get(safe_name(sid), ""),
            }
        )

    logger.info("Batch complete — %d samples processed", len(results))
    return results
