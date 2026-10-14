"""TabPFN Model Registry and Manifest Resolver.

Reads `models.tsv` to look up the calibrated decision boundary and feature set of
each known TabPFN preset. Preset model files are not bundled with STRiDE: the
caller always supplies the model file path, and the manifest is matched against
that file's name to find its calibrated threshold.
"""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

TABPFN_DIR = Path(__file__).parent
MANIFEST_PATH = TABPFN_DIR / "models.tsv"

# The one model file that does ship inside the package, with its calibrated threshold.
BUNDLED_MODEL_NAME = "tabpfn_finetuned.joblib"
BUNDLED_MODEL_THRESHOLD = 0.660658


@dataclass
class TabPFNModelInfo:
    model_id: str
    cohort: str
    rank: str
    combo_name: str
    threshold: float
    features: list[str]
    filename: str
    description: str


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
            feats = [x.strip() for x in row.get("features", "").split(",") if x.strip()]
            thr = float(row.get("threshold", 0.5))

            info = TabPFNModelInfo(
                model_id=row.get("model_id", "").strip(),
                cohort=row.get("cohort", "").strip(),
                rank=row.get("rank", "").strip(),
                combo_name=row.get("combo_name", "").strip(),
                threshold=thr,
                features=feats,
                filename=row.get("filename", "").strip(),
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


def find_preset(
    model_identifier: Optional[str],
    manifest: list[TabPFNModelInfo],
    cohort: Optional[str] = None,
) -> Optional[TabPFNModelInfo]:
    """Look up a preset by alias, combo name, or filename.

    Accepts short aliases ('ao_top1', 'ai_top2'), combo names ('SPECIAL_k2__ed__te'),
    cohort indicators ('access_only', 'tabpfn_access_impact') and model filenames.
    An empty identifier means the cohort's top preset ('ai_top1' for access_impact,
    otherwise 'ao_top1').
    """
    clean_id = (model_identifier or "").strip().lower()

    # Normalize method strings
    if clean_id.startswith("tabpfn_"):
        clean_id = clean_id[len("tabpfn_") :]
    if not clean_id or clean_id == "tabpfn" or clean_id == "default":
        clean_id = "ai_top1" if "impact" in (cohort or "").lower() else "ao_top1"
    elif clean_id in ("access_only", "ao"):
        clean_id = "ao_top1"
    elif clean_id in ("access_impact", "ai", "access_plus_impact"):
        clean_id = "ai_top1"

    # Match by model_id exact
    for m in manifest:
        if m.model_id.lower() == clean_id:
            return m

    # Match by combo name with optional cohort filter
    norm_cohort = (cohort or "").lower().replace("-", "_")
    for m in manifest:
        if m.combo_name.lower() == clean_id:
            if norm_cohort and norm_cohort not in m.cohort.lower():
                continue
            return m

    # Match by joblib filename
    clean_name = Path(clean_id).name
    for m in manifest:
        if m.filename.lower() in (clean_name, f"{clean_name}.joblib"):
            return m

    return None


def resolve_tabpfn_model(
    model_identifier: Optional[str] = None,
    custom_path: Optional[str | Path] = None,
    cohort: Optional[str] = None,
    manifest_path: Optional[Path] = None,
) -> tuple[Path, Optional[float], Optional[TabPFNModelInfo]]:
    """Resolve the model file path and its calibrated decision boundary threshold.

    Parameters
    ----------
    model_identifier : str | None
        Preset alias, combo name, cohort indicator or filename (see
        :func:`find_preset`). Only used to name the expected file when
        *custom_path* is missing.
    custom_path : str | Path | None
        Path to the model .joblib file. Required.
    cohort : str | None
        Optional cohort filter ('access_only' or 'access_impact') to resolve ambiguous combo names.

    Returns
    -------
    tuple[Path, float | None, TabPFNModelInfo | None]
        (joblib_path, calibrated_threshold_or_None, model_info_or_None).
        The threshold comes from the manifest entry matching the file name, else
        the ``_<float>.joblib`` filename suffix, else None.

    Raises
    ------
    FileNotFoundError
        If *custom_path* is not given.
    """
    manifest = load_manifest(manifest_path)

    if not custom_path:
        preset = find_preset(model_identifier, manifest, cohort)
        if preset is not None:
            hint = f" Preset '{preset.model_id}' expects a file named {preset.filename}."
        else:
            known = ", ".join(m.model_id for m in manifest)
            hint = f" Unknown preset '{model_identifier}' (known presets: {known})."
        raise FileNotFoundError(
            "No TabPFN model file given. Preset model files are not bundled with STRiDE; "
            "pass the model path with --model-joblib." + hint
        )

    p = Path(custom_path)
    match = next((m for m in manifest if m.filename == p.name), None)

    if model_identifier:
        requested = find_preset(model_identifier, manifest, cohort)
        if requested is not None and requested is not match:
            logger.warning(
                "--tabpfn-model '%s' does not match model file %s; using the model file.",
                model_identifier,
                p.name,
            )

    if match is not None:
        return p, match.threshold, match
    if p.name == BUNDLED_MODEL_NAME:
        return p, BUNDLED_MODEL_THRESHOLD, None
    return p, parse_threshold_from_filename(p), None


if __name__ == "__main__":
    manifest = load_manifest()
    print(f"Loaded {len(manifest)} models from manifest:")
    for m in manifest:
        print(f"  [{m.model_id:<7}] {m.cohort:<13} #{m.rank} | {m.combo_name:<32} | Thr: {m.threshold:.4f} | {m.filename}")
