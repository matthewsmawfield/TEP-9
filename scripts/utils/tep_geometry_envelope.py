"""Deterministic geometry envelope for TEP flyby predictions (path-resolved; Step 039 uses Step 008 pooled scale)."""

from __future__ import annotations

import math
from typing import Any, Dict, Optional

import numpy as np

from scripts.utils.physics import (
    DISFORMAL_COUPLING_STRENGTH,
    DISFORMAL_VELOCITY_THRESHOLD_KM_S,
    J2_EARTH,
    J3_EARTH,
    J4_EARTH,
    R_EARTH,
)

# |Ξ| bands for regime labels (manuscript-aligned; thresholds are documented outputs).
XI_REGIME_LOW: float = 0.03
XI_REGIME_HIGH: float = 0.25

# Near-symmetric trajectories require stronger multipole cancellation than the
# leading-order cos(delta_in) - cos(delta_out) factor alone.
ASYMMETRY_COHERENCE_THRESHOLD = 0.12

# Nominal geometry-envelope heuristics (Step 007 / Step 041 stress tests).
# Each carries a nominal ±50% systematic band in Monte Carlo sensitivity sweeps;
# they are not least-squares free parameters in Step 008 (k_fit = 1 for TEP restricted).
GEOMETRY_ENVELOPE_HEURISTIC_KEYS = (
    "inclination_scale",
    "j2_cos2_amp",
    "j2_alt_decay_km",
    "plasma_density_scale_cm3",
    "plasma_density_exponent",
    "velocity_screening_exponent",
    "asymmetry_coherence_threshold",
)


def default_geometry_envelope_heuristics() -> Dict[str, float]:
    """Return a fresh copy of nominal envelope heuristics (Step 007 defaults)."""
    return {
        "inclination_scale": 0.15,
        "j2_cos2_amp": 0.00054,
        "j2_alt_decay_km": 2000.0,
        "plasma_density_scale_cm3": 5000.0,
        "plasma_density_exponent": -0.3,
        "velocity_screening_exponent": 4.0,
        "asymmetry_coherence_threshold": float(ASYMMETRY_COHERENCE_THRESHOLD),
    }


def geometry_modulation_factors_from_heuristics(
    altitude_km: float,
    latitude_deg: float,
    velocity_km_s: float,
    plasma_density_cm3: float,
    heuristics: Dict[str, float],
    *,
    v_threshold_km_s: float = DISFORMAL_VELOCITY_THRESHOLD_KM_S,
) -> Dict[str, float]:
    """
    Deterministic modulation factors for the geometry envelope core (no plasma S_ansatz).

    ``heuristics`` must contain all keys in ``GEOMETRY_ENVELOPE_HEURISTIC_KEYS``.
    """
    missing = [k for k in GEOMETRY_ENVELOPE_HEURISTIC_KEYS if k not in heuristics]
    if missing:
        raise ValueError(
            f"geometry envelope heuristics missing keys: {missing}; "
            f"required: {GEOMETRY_ENVELOPE_HEURISTIC_KEYS}"
        )
    h = heuristics
    f_inclination = 1.0 + h["inclination_scale"] * abs(np.sin(np.radians(latitude_deg)))
    f_j2 = (1.0 - h["j2_cos2_amp"] * np.cos(np.radians(latitude_deg)) ** 2) * np.exp(
        -altitude_km / h["j2_alt_decay_km"]
    )
    f_plasma = (
        1.0 + plasma_density_cm3 / h["plasma_density_scale_cm3"]
    ) ** h["plasma_density_exponent"]
    exp_v = h["velocity_screening_exponent"]
    f_velocity = (
        (v_threshold_km_s / velocity_km_s) ** exp_v
        if velocity_km_s > v_threshold_km_s
        else 1.0
    )
    f_total_core = f_inclination * f_j2 * f_plasma * f_velocity
    return {
        "f_inclination": float(f_inclination),
        "f_j2": float(f_j2),
        "f_plasma": float(f_plasma),
        "f_velocity": float(f_velocity),
        "f_total_core": float(f_total_core),
    }


def zonal_harmonic_bracket(latitude_deg: float, r_m: float) -> float:
    """J2/J3/J4 zonal bracket scaled by (R_Earth / r)^2."""
    sin_lat = math.sin(math.radians(latitude_deg))
    r_ratio_sq = (R_EARTH / r_m) ** 2
    p4 = (35.0 * sin_lat**4 - 30.0 * sin_lat**2 + 3.0) / 8.0
    return (J2_EARTH + J3_EARTH * sin_lat + J4_EARTH * p4) * r_ratio_sq


def j2_only_bracket(latitude_deg: float, r_m: float) -> float:
    """Legacy J2-only bracket retained for component decomposition."""
    sin_lat = math.sin(math.radians(latitude_deg))
    r_ratio_sq = (R_EARTH / r_m) ** 2
    return J2_EARTH * r_ratio_sq


def near_symmetry_cancellation_factor(
    cos_asymmetry: float,
    coherence_threshold: float = ASYMMETRY_COHERENCE_THRESHOLD,
) -> float:
    """
    Quadratic suppression for trajectories with |cos asymmetry| below the
    coherence threshold. Unity for asymmetric detections.
    """
    if coherence_threshold <= 0:
        raise ValueError("coherence_threshold must be positive")
    asym = abs(cos_asymmetry)
    if asym >= coherence_threshold:
        return 1.0
    return (asym / coherence_threshold) ** 3


def disformal_envelope_factor(
    v_sc_m_s: float,
    cos_asymmetry: float,
    v_cmb_frame_kms: Optional[float] = None,
    sc_cmb_cos_theta: Optional[float] = None,
) -> float:
    """
  Disformal response with asymmetry-gated amplitude. The conformal-gradient
  term already carries cos asymmetry; only the disformal excess is gated.
  """
    v_th = DISFORMAL_VELOCITY_THRESHOLD_KM_S
    alpha_b = DISFORMAL_COUPLING_STRENGTH
    asym = abs(cos_asymmetry)

    if v_cmb_frame_kms is not None and sc_cmb_cos_theta is not None and sc_cmb_cos_theta > 0:
        v_km_s = v_cmb_frame_kms
    else:
        v_km_s = v_sc_m_s / 1e3

    if v_km_s > v_th and cos_asymmetry < 0:
        raw = -1.0 * abs(1.0 + alpha_b * (v_km_s / v_th))
    else:
        raw = 1.0 + alpha_b * (v_km_s / v_th) * math.copysign(1.0, cos_asymmetry or 0.0)

    return 1.0 + (raw - 1.0) * asym


def compute_disformal_transition_xi(
    *,
    velocity_km_s: float,
    cos_asymmetry: float,
    field_gradient_ratio: float,
    v_trans_km_s: float = DISFORMAL_VELOCITY_THRESHOLD_KM_S,
    xi_low: float = XI_REGIME_LOW,
    xi_high: float = XI_REGIME_HIGH,
) -> Dict[str, Any]:
    """
    Disformal transition scalar Ξ (velocity-activated, manuscript definition).

    Ξ = (v / v_trans)² × |asym| × (|∇φ| / |∇φ_⊕|) × sgn(asym)

    ``field_gradient_ratio`` must be the ratio of field-gradient magnitudes at
    perigee to the reference surface value (|∇φ(r_p)| / |∇φ(R_⊕)|), both from the
    same screening model instance.
    """
    if v_trans_km_s <= 0:
        raise ValueError("v_trans_km_s must be positive")
    if field_gradient_ratio < 0:
        raise ValueError("field_gradient_ratio must be non-negative")

    v_ratio = float(velocity_km_s) / float(v_trans_km_s)
    asym = float(cos_asymmetry)
    sgn = math.copysign(1.0, asym) if asym != 0.0 else 0.0
    xi = (v_ratio**2) * abs(asym) * float(field_gradient_ratio) * sgn
    xi_abs = abs(xi)

    if xi_abs < xi_low:
        regime = "conformal_dominated"
    elif xi_abs < xi_high:
        regime = "mixed"
    else:
        regime = "disformal_dominated"

    return {
        "xi": float(xi),
        "xi_abs": float(xi_abs),
        "regime_classification": regime,
        "v_trans_km_s": float(v_trans_km_s),
        "v_over_v_trans": float(v_ratio),
        "cos_dec_asymmetry": float(asym),
        "field_gradient_ratio": float(field_gradient_ratio),
        "definition": "(v/v_trans)^2 * |asym| * (|grad_phi|/|grad_phi_earth|) * sgn(asym)",
    }


def compose_geometry_envelope(
    *,
    altitude_km: float,
    latitude_deg: float,
    velocity_km_s: float,
    cos_asymmetry: float,
    plasma_density_cm3: float,
    modulation: Dict[str, float],
    plasma_screening_factor: float,
    plasma_sign_factor: float,
    asymmetry_coherence_threshold: float = ASYMMETRY_COHERENCE_THRESHOLD,
) -> Dict[str, Any]:
    """Aggregate deterministic geometry, plasma, and symmetry factors."""
    asymmetry_factor = near_symmetry_cancellation_factor(
        cos_asymmetry, coherence_threshold=asymmetry_coherence_threshold
    )
    envelope = (
        modulation["f_inclination"]
        * modulation["f_j2"]
        * modulation["f_plasma"]
        * modulation["f_velocity"]
        * asymmetry_factor
    )
    return {
        "f_inclination": modulation["f_inclination"],
        "f_j2": modulation["f_j2"],
        "f_plasma_core": modulation["f_plasma"],
        "f_velocity": modulation["f_velocity"],
        "f_geometry_core": modulation["f_total_core"],
        "plasma_density_cm3": plasma_density_cm3,
        "plasma_screening_factor": plasma_screening_factor,
        "plasma_sign_factor": plasma_sign_factor,
        "asymmetry_cancellation_factor": asymmetry_factor,
        "geometry_envelope": envelope,
    }
