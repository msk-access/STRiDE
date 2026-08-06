"""Shared utility functions for STRiDE.

Provides file-handling helpers, sample list parsing, and the centralized
Rich-powered logging setup used by the CLI layer.
"""

import logging
import os

import pandas as pd
from rich.logging import RichHandler

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------


def setup_logging(verbose: bool = False) -> None:
    """Configure the root logger with a Rich handler.

    Call once from the CLI entry point so that every module that uses
    ``logging.getLogger(__name__)`` inherits Rich formatting automatically.

    Parameters
    ----------
    verbose : bool
        If *True*, set level to DEBUG and show file paths in log output.
        Otherwise default to INFO.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, show_path=verbose)],
        force=True,  # override any prior basicConfig calls
    )


# ---------------------------------------------------------------------------
# Filename / path helpers
# ---------------------------------------------------------------------------


def safe_name(s: str) -> str:
    """Sanitise a string for use in filenames.

    Replaces any character that is not alphanumeric or in ``._-+=`` with
    an underscore.
    """
    return "".join(c if c.isalnum() or c in "._-+=" else "_" for c in s)


def strip_ext(path: str) -> str:
    """Return the basename of *path* without a recognised alignment extension.

    Recognised extensions: ``.bam``, ``.cram``, ``.sam``.  Falls back to
    :func:`os.path.splitext` for anything else.
    """
    base = os.path.basename(path)
    for ext in (".bam", ".cram", ".sam"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return os.path.splitext(base)[0]


# ---------------------------------------------------------------------------
# Sample-list I/O
# ---------------------------------------------------------------------------


def read_samples_list(list_path: str) -> list[dict[str, str]]:
    """Read a CSV/TSV sample manifest.

    The file must contain columns mappable to ``sample_id``, ``tumor_bam``
    and ``normal_bam``.  Header names are flexible — common aliases such as
    ``sample``, ``tumor``, ``normal_path``, etc. are accepted.

    An optional ``matched_norm_sample_barcode`` column (aliases:
    ``normal_barcode``, ``normal_sample_barcode``, ``matched_normal``)
    provides an explicit matched-normal identifier for the output.  When
    absent, the pipeline layer derives it from ``normal_bam`` via
    :func:`strip_ext`.

    Parameters
    ----------
    list_path : str
        Path to the tab- or comma-delimited manifest file.

    Returns
    -------
    list[dict]
        Each dict has keys ``sample_id``, ``tumor_bam``, ``normal_bam``,
        and ``matched_norm_sample_barcode`` (may be ``""``).

    Raises
    ------
    ValueError
        If required columns cannot be identified.
    """
    df = pd.read_csv(list_path, sep=None, engine="python")
    df.columns = [c.strip().lower() for c in df.columns]

    sample_cols = [c for c in df.columns if c in {"sample_id", "sample", "id"}]
    tumor_cols = [c for c in df.columns if c in {"tumor_bam", "tumor", "tumor_path"}]
    normal_cols = [c for c in df.columns if c in {"normal_bam", "normal", "normal_path"}]

    if not sample_cols or not tumor_cols or not normal_cols:
        raise ValueError(
            f"Samples list must contain columns for sample_id/tumor_bam/normal_bam. "
            f"Found columns: {df.columns.tolist()}"
        )

    s_col, t_col, n_col = sample_cols[0], tumor_cols[0], normal_cols[0]

    # Optional: matched-normal barcode column
    norm_bc_aliases = {
        "matched_norm_sample_barcode",
        "normal_barcode",
        "normal_sample_barcode",
        "matched_normal",
    }
    norm_bc_cols = [c for c in df.columns if c in norm_bc_aliases]
    norm_bc_col = norm_bc_cols[0] if norm_bc_cols else None

    out: list[dict[str, str]] = []
    for _, row in df.iterrows():
        sid = str(row[s_col]).strip()
        tb = str(row[t_col]).strip()
        nb = str(row[n_col]).strip()
        if not sid or sid.lower() == "nan":
            continue
        # Barcode defaults to "" when column is absent or value is NaN
        norm_bc = ""
        if norm_bc_col is not None:
            raw = str(row[norm_bc_col]).strip()
            if raw and raw.lower() != "nan":
                norm_bc = raw
        out.append(
            {
                "sample_id": sid,
                "tumor_bam": tb,
                "normal_bam": nb,
                "matched_norm_sample_barcode": norm_bc,
            }
        )
    return out
