import os
from typing import List, Dict, Optional

from .feature_generator import MSIProfileGenerator
from .predictor import predict_from_feature_tsvs, write_one_output_per_sample
from .utils import safe_name, read_samples_list, strip_ext

def run_feature_generation(
    sites_file: str,
    tumor_bam: str,
    normal_bam: str,
    out_features_tsv: str,
    min_coverage: int = 20,
    max_repeat_bins: int = 100
) -> str:
    os.makedirs(os.path.dirname(out_features_tsv) or ".", exist_ok=True)
    gen = MSIProfileGenerator(
        sites_file=sites_file,
        tumor_bam_path=tumor_bam,
        normal_bam_path=normal_bam,
        min_coverage=min_coverage,
        max_repeat_bins=max_repeat_bins
    )
    gen.run(out_features_tsv)
    return out_features_tsv

def run_end_to_end_single(
    sites_file: str,
    tumor_bam: str,
    normal_bam: str,
    model_joblib: str,
    out_dir: str,
    sample_id: Optional[str] = None,
    min_coverage: int = 20,
    max_repeat_bins: int = 100,
    keep_features: bool = True
) -> Dict[str, str]:
    os.makedirs(out_dir, exist_ok=True)

    sid = sample_id.strip() if sample_id else strip_ext(tumor_bam)

    features_dir = os.path.join(out_dir, "features")
    preds_dir = os.path.join(out_dir, "predictions")
    os.makedirs(features_dir, exist_ok=True)
    os.makedirs(preds_dir, exist_ok=True)

    feat_tsv = os.path.join(features_dir, f"msi_features_{safe_name(sid)}.tsv")

    run_feature_generation(
        sites_file, tumor_bam, normal_bam, feat_tsv,
        min_coverage=min_coverage, max_repeat_bins=max_repeat_bins
    )

    df_preds = predict_from_feature_tsvs(model_joblib, [feat_tsv])
    out_paths = write_one_output_per_sample(df_preds, preds_dir)

    if not keep_features:
        try:
            os.remove(feat_tsv)
        except OSError:
            pass

    return {"sample_id": sid, "features_tsv": feat_tsv, "prediction_txt": out_paths[0]}

def run_end_to_end_batch(
    samples_list_path: str,
    sites_file: str,
    model_joblib: str,
    out_dir: str,
    min_coverage: int = 20,
    max_repeat_bins: int = 100,
    keep_features: bool = True
) -> List[Dict[str, str]]:
    os.makedirs(out_dir, exist_ok=True)

    samples: List[Dict[str, str]] = read_samples_list(samples_list_path)

    features_dir = os.path.join(out_dir, "features")
    preds_dir = os.path.join(out_dir, "predictions")
    os.makedirs(features_dir, exist_ok=True)
    os.makedirs(preds_dir, exist_ok=True)

    # 1) generate feature TSV per sample
    feature_tsvs = []
    sample_ids = []
    for s in samples:
        sid = s["sample_id"]
        tumor_bam = s["tumor_bam"]
        normal_bam = s["normal_bam"]

        feat_tsv = os.path.join(features_dir, f"msi_features_{safe_name(sid)}.tsv")
        run_feature_generation(
            sites_file, tumor_bam, normal_bam, feat_tsv,
            min_coverage=min_coverage, max_repeat_bins=max_repeat_bins
        )
        feature_tsvs.append(feat_tsv)
        sample_ids.append(sid)

    # 2) predict all (batch) but write ONE output per sample
    df_preds = predict_from_feature_tsvs(model_joblib, feature_tsvs)
    out_paths = write_one_output_per_sample(df_preds, preds_dir)

    # 3) optional cleanup
    if not keep_features:
        for fp in feature_tsvs:
            try:
                os.remove(fp)
            except OSError:
                pass

    # 4) map outputs back to sample_id
    results = []
    pred_map = {os.path.basename(p).replace("_msi.txt", ""): p for p in out_paths}
    for sid, feat in zip(sample_ids, feature_tsvs):
        results.append({
            "sample_id": sid,
            "features_tsv": feat,
            "prediction_txt": pred_map.get(safe_name(sid), "")
        })
    return results
