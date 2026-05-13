"""IRI-calibrated reference density for plasma screening ansatz."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


def load_iri_trajectory_profiles(project_root: Path) -> dict[str, Any]:
    iri_file = project_root / "results" / "step033_iri_trajectory_profiles.json"
    if not iri_file.exists():
        raise FileNotFoundError(
            "IRI trajectory profiles are required (results/step033_iri_trajectory_profiles.json)"
        )
    with open(iri_file, encoding="utf-8") as handle:
        return json.load(handle)


def mission_peak_electron_density_cm3(
    profiles: dict[str, Any],
    mission_key: str,
) -> float | None:
    profile = profiles.get(mission_key)
    if not isinstance(profile, dict):
        return None
    densities = np.asarray(profile.get("trajectory", {}).get("iri_ne_cm3", []), dtype=float)
    if densities.size == 0:
        return None
    return float(np.max(densities))


def iri_reference_electron_density_cm3(project_root: Path) -> float:
    """Median peak IRI n_e along analyzed flyby trajectories (cm^-3)."""
    profiles = load_iri_trajectory_profiles(project_root)
    peaks: list[float] = []
    for data in profiles.values():
        densities = np.asarray(data["trajectory"]["iri_ne_cm3"], dtype=float)
        if densities.size:
            peaks.append(float(np.max(densities)))
    if not peaks:
        raise RuntimeError("No IRI electron densities available for n_ref calibration")
    return float(np.median(peaks))
