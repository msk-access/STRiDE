import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

ID_COLS_REQUIRED = ["chrom", "start", "repeat_unit", "repeat_count"]


def make_site_id(df: pd.DataFrame) -> pd.Series:
    """Build unique site IDs matching the format used in training: chrom:start:repeat_unit:repeat_count."""
    return (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(str)
        + ":"
        + df["repeat_unit"].astype(str)
        + ":"
        + df["repeat_count"].astype(str)
    )


def sample_tsv_to_wide(
    tsv_path: Path,
    features: list[str],
    expected_sites: Optional[list[str]] = None,
    sign_map: Optional[dict[str, float]] = None,
) -> tuple[pd.Series, list[str]]:
    """Convert a single sample TSV to a flattened wide Series, applying sign orientation if provided."""
    df = pd.read_csv(tsv_path, sep="\t")

    for col in ID_COLS_REQUIRED:
        if col not in df.columns:
            raise ValueError(f"{tsv_path.name} missing required column: {col}")

    base = pd.DataFrame({"site_id": make_site_id(df)})

    if base["site_id"].duplicated().any():
        dup_n = int(base["site_id"].duplicated().sum())
        raise ValueError(f"{tsv_path.name} has {dup_n} duplicate site_ids.")

    # Expose site_id feature values, applying sign map if present
    for col in features:
        if col in df.columns:
            val_series = pd.to_numeric(df[col], errors="coerce")
            if sign_map:
                # Apply sign map for specific site_id_feature key
                for idx, row_site in enumerate(base["site_id"]):
                    key = f"{row_site}_{col}"
                    if key in sign_map:
                        val_series.iloc[idx] = val_series.iloc[idx] * sign_map[key]
            base[col] = val_series
        else:
            base[col] = np.nan

    if expected_sites is not None:
        base = base.set_index("site_id").reindex(expected_sites).reset_index()

    base = base.sort_values("site_id")
    site_ids = base["site_id"].tolist()
    wide_df = base.set_index("site_id")[features]

    try:
        wide_s = wide_df.stack(future_stack=True)
    except TypeError:
        wide_s = wide_df.stack(dropna=False)

    wide_s.index = [f"{sid}__{feat}" for sid, feat in wide_s.index]
    return wide_s, site_ids


def process_single_sample(
    tsv_path: Path, expected_columns: list[str], sign_map: Optional[dict[str, float]] = None
) -> pd.Series:
    """
    Load a single sample TSV, extract features, apply optional sign map,
    and pivot to exactly match expected_columns schema.
    """
    df = pd.read_csv(tsv_path, sep="\t")

    # Check required ID columns
    for col in ID_COLS_REQUIRED:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' missing in {tsv_path.name}")

    # Parse expected columns: format is {site_id}__{feature_name}
    site_to_features = {}
    for col in expected_columns:
        parts = col.split("__")
        if len(parts) == 2:
            site_id, feat = parts
            if site_id not in site_to_features:
                site_to_features[site_id] = []
            site_to_features[site_id].append(feat)

    expected_sites = list(site_to_features.keys())
    features = list({feat for feats in site_to_features.values() for feat in feats})

    base = pd.DataFrame({"site_id": make_site_id(df)})

    if base["site_id"].duplicated().any():
        dup_count = base["site_id"].duplicated().sum()
        raise ValueError(f"Found {dup_count} duplicate site IDs in {tsv_path.name}")

    for feat in features:
        if feat in df.columns:
            val_series = pd.to_numeric(df[feat], errors="coerce")
            if sign_map:
                for idx, row_site in enumerate(base["site_id"]):
                    key = f"{row_site}_{feat}"
                    if key in sign_map:
                        val_series.iloc[idx] = val_series.iloc[idx] * sign_map[key]
            base[feat] = val_series
        else:
            base[feat] = np.nan

    base = base.set_index("site_id").reindex(expected_sites).reset_index()
    base = base.sort_values("site_id")

    wide_df = base.set_index("site_id")[features]
    try:
        wide_s = wide_df.stack(future_stack=True)
    except TypeError:
        wide_s = wide_df.stack(dropna=False)

    wide_s.index = [f"{sid}__{feat}" for sid, feat in wide_s.index]
    wide_s = wide_s.reindex(expected_columns)
    return wide_s


def build_matrix(
    records: pd.DataFrame,
    features: list[str],
    expected_sites: Optional[list[str]] = None,
    max_sites: int = 0,
    seed: int = 42,
    sign_map: Optional[dict[str, float]] = None,
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """Stack sample vectors from records into wide feature matrix X and target y."""
    if records.empty:
        raise ValueError("No records were provided to build_matrix().")

    if expected_sites is None:
        first_path = Path(records.iloc[0]["path"])
        _, expected_sites = sample_tsv_to_wide(
            first_path, features, expected_sites=None, sign_map=sign_map
        )

        if max_sites and len(expected_sites) > max_sites:
            rng = np.random.default_rng(seed)
            expected_sites = sorted(
                rng.choice(expected_sites, size=max_sites, replace=False).tolist()
            )

    rows = []
    labels = []
    sample_ids = []

    for rec in records.itertuples(index=False):
        wide_s, _ = sample_tsv_to_wide(
            Path(rec.path), features=features, expected_sites=expected_sites, sign_map=sign_map
        )
        wide_s.name = rec.sample_id
        rows.append(wide_s)
        labels.append(int(rec.label))
        sample_ids.append(str(rec.sample_id))

    X = pd.DataFrame(rows)
    X.index = sample_ids
    X.index.name = "sample_id"
    X = X.reindex(sorted(X.columns), axis=1)
    y = np.asarray(labels, dtype=int)
    return X, y, expected_sites


def collect_binary_records(
    msi_roots: list[Path], mss_roots: list[Path], cohort: str, file_ext: str = "tsv"
) -> pd.DataFrame:
    """Scan directories for samples and construct metadata DataFrame with labels."""
    msi_files: list[Path] = []
    mss_files: list[Path] = []

    for root in msi_roots:
        if root.is_dir():
            msi_files.extend(sorted(root.rglob(f"*.{file_ext}")))
    for root in mss_roots:
        if root.is_dir():
            mss_files.extend(sorted(root.rglob(f"*.{file_ext}")))

    if not msi_files and not mss_files:
        return empty_binary_records()

    rows = []
    for path in msi_files:
        rows.append({"sample_id": path.stem, "label": 1, "cohort": cohort, "path": str(path)})
    for path in mss_files:
        rows.append({"sample_id": path.stem, "label": 0, "cohort": cohort, "path": str(path)})

    df = pd.DataFrame(rows)
    df = df.sort_values("sample_id").reset_index(drop=True)
    return df


def empty_binary_records() -> pd.DataFrame:
    """Return an empty DataFrame matching the records schema."""
    return pd.DataFrame(columns=["sample_id", "label", "cohort", "path"])


def expand_roots(roots_csv: str) -> list[Path]:
    """Helper to split a comma-separated list of roots and return Path objects."""
    return [Path(r.strip()) for r in (roots_csv or "").split(",") if r.strip()]


def collect_records_from_roots(
    msi_roots_csv: str, mss_roots_csv: str, cohort: str, file_ext: str = "tsv"
) -> pd.DataFrame:
    """Collect sample records from comma-separated list of root directories."""
    msi_roots = expand_roots(msi_roots_csv)
    mss_roots = expand_roots(mss_roots_csv)
    return collect_binary_records(msi_roots, mss_roots, cohort, file_ext)


def load_sfs_sign_map(sfs_path: str, rank1_only: bool = True) -> dict[str, float]:
    """Load SFS directions mapping to determine sign multipliers."""
    if not sfs_path or not os.path.exists(sfs_path):
        return {}
    df = pd.read_csv(sfs_path, sep="\t")
    if "best_metric" in df.columns and "metric" not in df.columns:
        df = df.rename(columns={"best_metric": "metric"})
    if rank1_only and "rank" in df.columns:
        df = df[df["rank"] == 1].copy()
    required = {"site_id", "metric", "direction"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"SFS TSV missing columns: {missing}")
    sign = np.where(df["direction"].astype(str).str.contains("higher→MSI"), 1.0, -1.0)
    cols = df["site_id"].astype(str) + "_" + df["metric"].astype(str)
    return dict(zip(cols, sign))
