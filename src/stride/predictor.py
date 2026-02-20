import os
import glob
import joblib
import numpy as np
import pandas as pd
from typing import Optional, Dict, List, Any, Tuple

from .utils import safe_name

META_IGNORE_COLS = {"chrom", "chr", "start", "end", "ref", "alt", "site_id", "label", "sample_id"}

def extract_all_features_from_tsv(file_path: str) -> Optional[Dict[str, float]]:
    try:
        if os.path.getsize(file_path) == 0:
            print(f"Skipping empty file: {file_path}")
            return None
        df = pd.read_csv(file_path, sep="\t")
        if "start" not in df.columns or ("chrom" not in df.columns and "chr" not in df.columns):
            print(f"Skipping malformed file (missing chrom/start): {file_path}")
            return None
        chrom_col = "chrom" if "chrom" in df.columns else "chr"
        df["site_id"] = df[chrom_col].astype(str) + "_" + df["start"].astype(str)

        metric_cols = [c for c in df.columns if c not in META_IGNORE_COLS]
        feats: Dict[str, float] = {}
        for _, row in df.iterrows():
            site = row["site_id"]
            for m in metric_cols:
                feats[f"{site}_{m}"] = row.get(m, np.nan)
        return feats
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def gather_samples_from_inputs(
    samples_dir: Optional[str],
    sample_files: Optional[List[str]],
    list_file: Optional[str]
) -> Tuple[List[str], Dict[str, Dict[str, float]]]:
    tsvs: List[str] = []

    if samples_dir:
        tsvs.extend(glob.glob(os.path.join(samples_dir, "**", "*.tsv"), recursive=True))

    if sample_files:
        tsvs.extend(sample_files)

    if list_file:
        with open(list_file, "r") as f:
            for ln in f:
                p = ln.strip()
                if p:
                    tsvs.append(p)

    seen = set()
    tsvs2 = []
    for p in tsvs:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            tsvs2.append(ap)
    tsvs = tsvs2

    if not tsvs:
        raise FileNotFoundError("No .tsv sample files found.")

    feature_bag: Dict[str, Dict[str, float]] = {}
    sample_ids: List[str] = []

    for fp in tsvs:
        fn = os.path.basename(fp)
        sample_id = os.path.splitext(fn)[0]
        if sample_id.startswith("msi_features_"):
            sample_id = sample_id[len("msi_features_"):]
        sample_id = sample_id.strip()

        feats = extract_all_features_from_tsv(fp)
        if feats is None:
            continue

        feature_bag[sample_id] = feats
        sample_ids.append(sample_id)

    sample_ids = list(dict.fromkeys(sample_ids))
    if not sample_ids:
        raise RuntimeError("No valid sample TSVs could be parsed.")
    return sample_ids, feature_bag

def unwrap_model(obj: Any):
    from sklearn.pipeline import Pipeline

    if hasattr(obj, "named_steps"):
        return obj

    if isinstance(obj, dict):
        for k in ("pipeline", "model", "estimator"):
            val = obj.get(k)
            if hasattr(val, "named_steps"):
                return val
        if "scaler" in obj and "clf" in obj:
            try:
                return Pipeline([("standardscaler", obj["scaler"]), ("sgdclassifier", obj["clf"])])
            except Exception:
                return None
    return None

def get_expected_features(model) -> np.ndarray:
    scaler = None
    if hasattr(model, "named_steps"):
        scaler = model.named_steps.get("standardscaler")

    if scaler is not None:
        if hasattr(scaler, "feature_names_in_"):
            return scaler.feature_names_in_
        if hasattr(scaler, "get_feature_names_out"):
            return scaler.get_feature_names_out()

    for part in ("sgdclassifier", "svc", "svm", "estimator", "linearsvc"):
        if hasattr(model, "named_steps") and part in model.named_steps:
            est = model.named_steps[part]
            if hasattr(est, "feature_names_in_"):
                return est.feature_names_in_

    if hasattr(model, "feature_names_in_"):
        return model.feature_names_in_

    raise ValueError("Cannot determine expected feature names from the model/scaler.")

def build_matrix(sample_ids: List[str], feature_bag: Dict[str, Dict[str, float]], expected_features: np.ndarray) -> pd.DataFrame:
    rows = []
    for sid in sample_ids:
        feats = feature_bag.get(sid, {})
        row = {k: feats.get(k, 0.0) for k in expected_features}
        rows.append(row)
    return pd.DataFrame(rows, index=sample_ids)

def get_scores(model, X: pd.DataFrame) -> np.ndarray:
    try:
        scores = model.decision_function(X)
        scores = np.asarray(scores)
        if scores.ndim > 1:
            scores = scores[:, -1]
        return scores
    except Exception:
        pass

    try:
        proba = getattr(model, "predict_proba", None)
        if proba is not None:
            p = np.asarray(proba(X))
            if p.ndim == 2 and p.shape[1] > 1:
                return p[:, 1]
    except Exception:
        pass

    return np.full((X.shape[0],), np.nan, dtype=float)

def load_model(model_joblib: str):
    raw = joblib.load(model_joblib)
    model = unwrap_model(raw)
    if model is None:
        raise ValueError("Unsupported model object: not a Pipeline (or dict containing one).")
    return model

def predict_from_feature_tsvs(model_joblib: str, feature_tsvs: List[str]) -> pd.DataFrame:
    model = load_model(model_joblib)
    expected = get_expected_features(model)

    sample_ids, feature_bag = gather_samples_from_inputs(None, feature_tsvs, None)
    X = build_matrix(sample_ids, feature_bag, expected)
    y_pred = np.asarray(model.predict(X)).astype(int)
    scores = np.asarray(get_scores(model, X), dtype=float)

    return pd.DataFrame({
        "Sample_ID": sample_ids,
        "MSI_class_predicted": y_pred,
        "msi_score": np.round(scores, 6)
    })

def write_one_output_per_sample(df_preds: pd.DataFrame, out_dir: str) -> List[str]:
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for _, row in df_preds.iterrows():
        sid = str(row["Sample_ID"])
        out_path = os.path.join(out_dir, f"{safe_name(sid)}_msi.txt")
        pd.DataFrame([row]).to_csv(out_path, sep="\t", index=False)
        written.append(out_path)
    return written
