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


def plasma_screening_sensitivity(
    plasma_density_cm3: float,
    n_ref: float | None = None,
) -> dict[str, float]:
    """
    Compute plasma screening factor under multiple functional forms.

    Tests sensitivity of the plasma attenuation to the choice of
    phenomenological ansatz. The three forms span plausible physical
    behaviours for scalar-field screening in a plasma:

    - exponential: S = exp(-n_e / n_ref)  [current default]
    - power_law:   S = (1 + n_e/n_ref)^(-gamma) with gamma=0.3
    - linear:      S = max(0, 1 - n_e/n_ref)  [most aggressive]

    Parameters
    ----------
    plasma_density_cm3 : float
        Electron density at perigee (cm^-3).
    n_ref : float or None
        Reference density. If None, uses a fixed 10^4 cm^-3 scale.

    Returns
    -------
    dict mapping functional form name to screening factor in [0, 1].
    """
    if n_ref is None:
        n_ref = 1e4

    ne = max(plasma_density_cm3, 0.0)
    ratio = ne / n_ref

    return {
        'exponential': float(np.exp(-ratio)),
        'power_law': float((1.0 + ratio) ** (-0.3)),
        'linear': float(max(0.0, 1.0 - ratio)),
        'n_e_cm3': float(ne),
        'n_ref_cm3': float(n_ref),
        'ratio': float(ratio),
    }


def plasma_envelope_spread(
    per_flyby_densities: dict[str, float],
    n_ref: float | None = None,
) -> dict[str, Any]:
    """
    Quantify how much the choice of plasma screening ansatz affects
    the spread in geometry envelope factors across flybys.

    For each flyby, computes the envelope contribution under all three
    functional forms. Reports the coefficient of variation (std/mean)
    of the envelope factors under each form. If the CV is similar
    across forms, the conclusions are robust to the plasma model choice.

    Parameters
    ----------
    per_flyby_densities : dict
        Mapping of flyby name -> peak electron density (cm^-3).
    n_ref : float or None
        Reference density.

    Returns
    -------
    dict with per-form CV and per-flyby breakdown.
    """
    if n_ref is None:
        n_ref = 1e4

    forms = ['exponential', 'power_law', 'linear']
    per_form: dict[str, list[float]] = {f: [] for f in forms}
    per_flyby: dict[str, dict[str, float]] = {}

    for name, ne in per_flyby_densities.items():
        sens = plasma_screening_sensitivity(ne, n_ref)
        per_flyby[name] = {f: sens[f] for f in forms}
        for f in forms:
            per_form[f].append(sens[f])

    result: dict[str, Any] = {
        'n_ref_cm3': n_ref,
        'n_flybys': len(per_flyby_densities),
        'per_flyby': per_flyby,
    }
    for f in forms:
        vals = np.array(per_form[f])
        result[f] = {
            'min': float(np.min(vals)),
            'max': float(np.max(vals)),
            'median': float(np.median(vals)),
            'cv': float(np.std(vals) / np.mean(vals)) if np.mean(vals) > 0 else float('inf'),
            'spread_ratio': float(np.max(vals) / np.min(vals)) if np.min(vals) > 0 else float('inf'),
        }

    return result
