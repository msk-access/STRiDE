import os
import pandas as pd
from typing import List, Dict, Optional

def safe_name(s: str) -> str:
    return "".join([c if c.isalnum() or c in "._-+=" else "_" for c in s])

def strip_ext(path: str) -> str:
    base = os.path.basename(path)
    for ext in (".bam", ".cram", ".sam"):
        if base.endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]

def read_samples_list(list_path: str) -> List[Dict[str, str]]:
    """
    Read a samples list file containing at least:
      - sample_id
      - tumor_bam
      - normal_bam

    Accepts TSV/CSV (auto-detected). Header names are flexible:
      sample_id: sample_id, sample, id
      tumor_bam: tumor_bam, tumor, tumor_path
      normal_bam: normal_bam, normal, normal_path
    """
    df = pd.read_csv(list_path, sep=None, engine="python")
    df.columns = [c.strip().lower() for c in df.columns]

    sample_cols = [c for c in df.columns if c in {"sample_id", "sample", "id"}]
    tumor_cols  = [c for c in df.columns if c in {"tumor_bam", "tumor", "tumor_path"}]
    normal_cols = [c for c in df.columns if c in {"normal_bam", "normal", "normal_path"}]

    if not sample_cols or not tumor_cols or not normal_cols:
        raise ValueError(
            f"Samples list must contain columns for sample_id/tumor_bam/normal_bam. "
            f"Found columns: {df.columns.tolist()}"
        )

    s_col, t_col, n_col = sample_cols[0], tumor_cols[0], normal_cols[0]

    out = []
    for _, row in df.iterrows():
        sid = str(row[s_col]).strip()
        tb  = str(row[t_col]).strip()
        nb  = str(row[n_col]).strip()
        if not sid or sid.lower() == "nan":
            continue
        out.append({"sample_id": sid, "tumor_bam": tb, "normal_bam": nb})
    return out
