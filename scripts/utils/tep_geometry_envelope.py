"""Deterministic geometry envelope for universal-beta TEP flyby predictions."""

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

# Near-symmetric trajectories require stronger multipole cancellation than the
# leading-order cos(delta_in) - cos(delta_out) factor alone.
ASYMMETRY_COHERENCE_THRESHOLD = 0.12


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


def near_symmetry_cancellation_factor(cos_asymmetry: float) -> float:
    """
  Quadratic suppression for trajectories with |cos asymmetry| below the
  coherence threshold. Unity for asymmetric detections.
  """
    asym = abs(cos_asymmetry)
    if asym >= ASYMMETRY_COHERENCE_THRESHOLD:
        return 1.0
    return (asym / ASYMMETRY_COHERENCE_THRESHOLD) ** 3


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
) -> Dict[str, Any]:
    """Aggregate deterministic geometry, plasma, and symmetry factors."""
    asymmetry_factor = near_symmetry_cancellation_factor(cos_asymmetry)
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
