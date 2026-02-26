import glob
import logging
import os
from typing import Any, Optional

import joblib
import numpy as np
import pandas as pd

from .utils import safe_name

# Module logger — configured by the CLI layer via utils.setup_logging()
logger = logging.getLogger(__name__)

META_IGNORE_COLS = {"chrom", "chr", "start", "end", "ref", "alt", "site_id", "label", "sample_id"}


def extract_all_features_from_tsv(file_path: str) -> Optional[dict[str, float]]:
    try:
        if os.path.getsize(file_path) == 0:
            logger.warning("Skipping empty file: %s", file_path)
            return None
        df = pd.read_csv(file_path, sep="\t")
        if "start" not in df.columns or ("chrom" not in df.columns and "chr" not in df.columns):
            logger.warning("Skipping malformed file (missing chrom/start): %s", file_path)
            return None
        chrom_col = "chrom" if "chrom" in df.columns else "chr"
        df["site_id"] = df[chrom_col].astype(str) + "_" + df["start"].astype(str)

        metric_cols = [c for c in df.columns if c not in META_IGNORE_COLS]
        feats: dict[str, float] = {}
        for _, row in df.iterrows():
            site = row["site_id"]
            for m in metric_cols:
                feats[f"{site}_{m}"] = row.get(m, np.nan)
        return feats
    except Exception as e:
        logger.error("Error reading file %s: %s", file_path, e)
        return None


def gather_samples_from_inputs(
    samples_dir: Optional[str], sample_files: Optional[list[str]], list_file: Optional[str]
) -> tuple[list[str], dict[str, dict[str, float]]]:
    tsvs: list[str] = []

    if samples_dir:
        tsvs.extend(glob.glob(os.path.join(samples_dir, "**", "*.tsv"), recursive=True))

    if sample_files:
        tsvs.extend(sample_files)

    if list_file:
        with open(list_file) as f:
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

    feature_bag: dict[str, dict[str, float]] = {}
    sample_ids: list[str] = []

    for fp in tsvs:
        fn = os.path.basename(fp)
        sample_id = os.path.splitext(fn)[0]
        if sample_id.startswith("msi_features_"):
            sample_id = sample_id[len("msi_features_") :]
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
            return np.asarray(scaler.feature_names_in_)  # type: ignore[no-any-return]
        if hasattr(scaler, "get_feature_names_out"):
            return np.asarray(scaler.get_feature_names_out())  # type: ignore[no-any-return]

    for part in ("sgdclassifier", "svc", "svm", "estimator", "linearsvc"):
        if hasattr(model, "named_steps") and part in model.named_steps:
            est = model.named_steps[part]
            if hasattr(est, "feature_names_in_"):
                return np.asarray(est.feature_names_in_)  # type: ignore[no-any-return]

    if hasattr(model, "feature_names_in_"):
        return np.asarray(model.feature_names_in_)  # type: ignore[no-any-return]

    raise ValueError("Cannot determine expected feature names from the model/scaler.")


def build_matrix(
    sample_ids: list[str], feature_bag: dict[str, dict[str, float]], expected_features: np.ndarray
) -> pd.DataFrame:
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
        return scores  # type: ignore[no-any-return]
    except Exception:
        pass

    try:
        proba = getattr(model, "predict_proba", None)
        if proba is not None:
            p = np.asarray(proba(X))
            if p.ndim == 2 and p.shape[1] > 1:
                return np.asarray(p[:, 1])  # type: ignore[no-any-return]
    except Exception:
        pass

    return np.full((X.shape[0],), np.nan, dtype=float)


def load_model(model_joblib: str):
    raw = joblib.load(model_joblib)
    model = unwrap_model(raw)
    if model is None:
        raise ValueError("Unsupported model object: not a Pipeline (or dict containing one).")
    return model


def map_msi_label(prediction: int) -> str:
    """Map an integer model prediction to a MAF-style MSI status label.

    Parameters
    ----------
    prediction : int
        Binary model output (``1`` = MSI-H positive, ``0`` = not MSI).

    Returns
    -------
    str
        ``"MSI"`` when *prediction* is ``1``, ``"NA"`` otherwise.
    """
    return "MSI" if prediction == 1 else "NA"


def predict_from_feature_tsvs(
    model_joblib: str,
    feature_tsvs: list[str],
    normal_barcodes: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Load a model and predict MSI status from feature TSVs.

    Parameters
    ----------
    model_joblib : str
        Path to the serialised scikit-learn model (``.joblib``).
    feature_tsvs : list[str]
        Paths to per-sample feature TSV files.
    normal_barcodes : list[str] | None
        Optional matched-normal sample barcodes, one per sample.  When
        *None*, the ``Matched_Norm_Sample_Barcode`` column is filled with
        empty strings.

    Returns
    -------
    pd.DataFrame
        MAF-aligned prediction table with columns
        ``Tumor_Sample_Barcode``, ``Matched_Norm_Sample_Barcode``,
        ``MSI_class_predicted`` (``"MSI"`` / ``"NA"``), and ``msi_score``.
    """
    model = load_model(model_joblib)
    expected = get_expected_features(model)

    sample_ids, feature_bag = gather_samples_from_inputs(None, feature_tsvs, None)
    X = build_matrix(sample_ids, feature_bag, expected)
    y_pred = np.asarray(model.predict(X)).astype(int)
    scores = np.asarray(get_scores(model, X), dtype=float)

    msi_labels = [map_msi_label(p) for p in y_pred]
    norm_bc = normal_barcodes if normal_barcodes else [""] * len(sample_ids)

    logger.info(
        "Prediction complete: %d sample(s), %d MSI, %d NA",
        len(sample_ids),
        msi_labels.count("MSI"),
        msi_labels.count("NA"),
    )

    return pd.DataFrame(
        {
            "Tumor_Sample_Barcode": sample_ids,
            "Matched_Norm_Sample_Barcode": norm_bc,
            "MSI_class_predicted": msi_labels,
            "msi_score": np.round(scores, 6),
        }
    )


def write_one_output_per_sample(df_preds: pd.DataFrame, out_dir: str) -> list[str]:
    """Write one tab-delimited prediction file per sample.

    Each file is named ``<tumor_barcode>_msi.txt`` and contains the full
    prediction row (MAF-aligned columns).

    Parameters
    ----------
    df_preds : pd.DataFrame
        Prediction table produced by :func:`predict_from_feature_tsvs`.
    out_dir : str
        Directory to write the per-sample files into.

    Returns
    -------
    list[str]
        Paths to the written output files.
    """
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for _, row in df_preds.iterrows():
        sid = str(row["Tumor_Sample_Barcode"])
        out_path = os.path.join(out_dir, f"{safe_name(sid)}_msi.txt")
        pd.DataFrame([row]).to_csv(out_path, sep="\t", index=False)
        logger.debug("Wrote prediction: %s", out_path)
        written.append(out_path)
    return written
