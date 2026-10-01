"""TabPFN Model Registry and Manifest Resolver.

Reads `models.tsv` to discover bundled models, calibrated decision boundaries,
and feature sets. Supports resolution by short alias, combo name, cohort, or path.
"""

from __future__ import annotations

import csv
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

TABPFN_DIR = Path(__file__).parent
MANIFEST_PATH = TABPFN_DIR / "models.tsv"


@dataclass
class TabPFNModelInfo:
    model_id: str
    cohort: str
    rank: str
    combo_name: str
    threshold: float
    features: list[str]
    artifact_path: Path
    description: str

    @property
    def exists(self) -> bool:
        return self.artifact_path.is_symlink() or self.artifact_path.exists()


def load_manifest(manifest_path: Optional[Path] = None) -> list[TabPFNModelInfo]:
    """Load all registered models from models.tsv."""
    path = manifest_path or MANIFEST_PATH
    if not path.exists():
        logger.warning("TabPFN models manifest not found at: %s", path)
        return []

    models: list[TabPFNModelInfo] = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            artifact_rel = row.get("artifact_relpath", "").strip()
            artifact_path = TABPFN_DIR / artifact_rel
            feats = [x.strip() for x in row.get("features", "").split(",") if x.strip()]
            thr = float(row.get("threshold", 0.5))

            info = TabPFNModelInfo(
                model_id=row.get("model_id", "").strip(),
                cohort=row.get("cohort", "").strip(),
                rank=row.get("rank", "").strip(),
                combo_name=row.get("combo_name", "").strip(),
                threshold=thr,
                features=feats,
                artifact_path=artifact_path,
                description=row.get("description", "").strip(),
            )
            models.append(info)
    return models


def parse_threshold_from_filename(filename_or_path: str | Path) -> Optional[float]:
    """Extract trailing decision boundary float from filename.
    
    e.g. SPECIAL_k2__ed__te_0.6975.joblib -> 0.6975
    """
    name = Path(filename_or_path).name
    match = re.search(r"_([0-9]+\.[0-9]+)(?:\.joblib)?$", name)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None


def resolve_tabpfn_model(
    model_identifier: Optional[str] = None,
    custom_path: Optional[str | Path] = None,
    cohort: Optional[str] = None,
    manifest_path: Optional[Path] = None,
) -> tuple[Path, float, Optional[TabPFNModelInfo]]:
    """Resolve model file path and calibrated decision boundary threshold.

    Parameters
    ----------
    model_identifier : str | None
        Can be:
        - short alias: 'ao_top1', 'ao_top2', 'ai_top1', 'ai_top2', etc.
        - combo name: 'SPECIAL_k2__ed__te', 'SPECIAL_k3__wd__te__ne', etc.
        - cohort indicator: 'access_only', 'access_impact', 'tabpfn_access_only', 'tabpfn_access_impact'
        - filename: 'SPECIAL_k2__ed__te_0.6975.joblib'
    custom_path : str | Path | None
        Explicit file path to a .joblib model.
    cohort : str | None
        Optional cohort filter ('access_only' or 'access_impact') to resolve ambiguous combo names.

    Returns
    -------
    tuple[Path, float, TabPFNModelInfo | None]
        (joblib_path, calibrated_threshold, model_info_or_None)
    """
    manifest = load_manifest(manifest_path)

    # 1. If explicit custom path provided
    if custom_path:
        p = Path(custom_path)
        # Check if matches any known manifest model
        for m in manifest:
            if m.artifact_path.name == p.name or str(m.artifact_path) == str(p):
                return m.artifact_path, m.threshold, m
        # Fallback to regex from filename
        thr = parse_threshold_from_filename(p)
        default_thr = thr if thr is not None else 0.5
        return p, default_thr, None

    clean_id = (model_identifier or "").strip().lower()

    # Normalize method strings
    if clean_id.startswith("tabpfn_"):
        clean_id = clean_id[len("tabpfn_") :]
    if not clean_id or clean_id == "tabpfn" or clean_id == "default":
        # Global default: Access-Only Rank 1 (ed + te)
        clean_id = "ao_top1"
    elif clean_id in ("access_only", "ao"):
        clean_id = "ao_top1"
    elif clean_id in ("access_impact", "ai", "access_plus_impact"):
        clean_id = "ai_top1"

    # Match by model_id exact
    for m in manifest:
        if m.model_id.lower() == clean_id:
            return m.artifact_path, m.threshold, m

    # Match by combo name with optional cohort filter
    norm_cohort = (cohort or "").lower().replace("-", "_")
    for m in manifest:
        if m.combo_name.lower() == clean_id:
            if norm_cohort and norm_cohort not in m.cohort.lower():
                continue
            return m.artifact_path, m.threshold, m

    # Match by joblib filename
    for m in manifest:
        if m.artifact_path.name.lower() == clean_id or m.artifact_path.stem.lower() == clean_id:
            return m.artifact_path, m.threshold, m

    # If clean_id looks like a direct file path on disk
    if model_identifier and ("/" in model_identifier or model_identifier.endswith(".joblib")):
        cand_path = Path(model_identifier)
        if cand_path.is_symlink() or cand_path.exists():
            thr = parse_threshold_from_filename(cand_path)
            return cand_path, thr if thr is not None else 0.5, None

    # Fallback default: Access-Only Rank 1
    if manifest:
        top_model = manifest[0]
        logger.warning(
            "Could not resolve model identifier '%s'. Falling back to default %s (%s).",
            model_identifier,
            top_model.model_id,
            top_model.description,
        )
        return top_model.artifact_path, top_model.threshold, top_model

    # Ultimate fallback if manifest is missing
    fallback_joblib = TABPFN_DIR / "tabpfn_finetuned.joblib"
    return fallback_joblib, 0.660658, None


if __name__ == "__main__":
    manifest = load_manifest()
    print(f"Loaded {len(manifest)} models from manifest:")
    for m in manifest:
        exists_mark = "OK" if m.exists else "MISSING"
        print(f"  [{m.model_id:<7}] {m.cohort:<13} #{m.rank} | {m.combo_name:<32} | Thr: {m.threshold:.4f} | [{exists_mark}]")
