"""Shared flyby selection rules for ensemble β fitting and model comparison."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

from scripts.utils.physics import validate_screened_coupling

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
_PIPELINE_CONFIG_PATH = PROJECT_ROOT / "config" / "pipeline_config.json"


def strict_sign_gate_from_config() -> bool:
    """
    When True, ensemble rows require Δv_obs and Δv_TEP(β_ref) to share a sign.

    Default is False: S/N-qualified rows remain in the likelihood layer; the
    sign-gated subset is reported separately by Step 008.
    """
    with open(_PIPELINE_CONFIG_PATH, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    tep = cfg["parameters"]["analysis"]["tep_physics"]
    return bool(tep.get("strict_sign_gate", False))


# Frozen policy text for JSON audit exports (Step 026 / provenance).
ENSEMBLE_GATE_POLICY = (
    "Likelihood rows require a published (Δv_obs, σ) pair with no imputation and "
    "S/N = |Δv_obs|/σ > 2. Sign agreement between Δv_obs and the TEP reference "
    "prediction at β_ref is enforced only when parameters.analysis.tep_physics."
    "strict_sign_gate is true (default: false). Exclusion reasons are explicit "
    "(missing_observation, snr_below_threshold, sign_mismatch, missing_prediction) "
    "rather than silent omission."
)

SNR_THRESHOLD: float = 2.0


def flyby_snr(dv_obs: Optional[float], dv_unc: Optional[float]) -> float:
    if dv_obs is None or dv_unc is None or dv_unc <= 0:
        return 0.0
    return abs(dv_obs) / dv_unc


def flyby_sign_product(dv_tep: Optional[float], dv_obs: Optional[float]) -> Optional[float]:
    if dv_tep is None or dv_obs is None:
        return None
    return dv_tep * dv_obs


def flyby_ensemble_exclusion_reason(
    prediction: Dict[str, Any],
    *,
    enforce_sign_agreement: Optional[bool] = None,
) -> Optional[str]:
    observed = prediction.get("observed", {})
    dv_obs = observed.get("dv_obs_mm_s")
    dv_unc = observed.get("sigma_mm_s")
    if dv_obs is None or dv_unc is None:
        return "missing_observation"

    snr = flyby_snr(dv_obs, dv_unc)
    if snr < SNR_THRESHOLD:
        return "snr_below_threshold"

    if enforce_sign_agreement is None:
        enforce_sign_agreement = strict_sign_gate_from_config()

    dv_tep = prediction.get("tep_predictions", {}).get("dv_tep_mm_s")
    sign_product = flyby_sign_product(dv_tep, dv_obs)
    if sign_product is None:
        return "missing_prediction"
    if enforce_sign_agreement and sign_product < 0:
        return "sign_mismatch"
    return None


def flyby_ensemble_eligible(prediction: Dict[str, Any]) -> bool:
    return flyby_ensemble_exclusion_reason(prediction) is None


def load_step008_ensemble_summary(step008: Dict[str, Any]) -> Dict[str, float]:
    overall = step008.get("overall_analysis")
    if not overall:
        raise ValueError("step008_fitting_results.json is missing overall_analysis")

    beta_stats = overall.get("beta_statistics")
    beta_eff_stats = overall.get("beta_eff_statistics")
    if not beta_stats or beta_stats.get("weighted_mean") is None:
        raise ValueError("step008_fitting_results.json is missing beta_statistics.weighted_mean")
    if not beta_eff_stats or beta_eff_stats.get("weighted_mean") is None:
        raise ValueError("step008_fitting_results.json is missing beta_eff_statistics.weighted_mean")

    uncertainty = overall.get("recommended_uncertainty")
    if uncertainty is None:
        uncertainty = beta_stats.get("weighted_uncertainty")
    if uncertainty is None:
        raise ValueError("step008_fitting_results.json is missing recommended_uncertainty")

    return {
        "beta_weighted": float(beta_stats["weighted_mean"]),
        "beta_uncertainty": float(uncertainty),
        "beta_eff_weighted": float(beta_eff_stats["weighted_mean"]),
        "n_fits": int(overall.get("n_fits", 0)),
    }


def load_validated_step008_ensemble_summary(step008: Dict[str, Any]) -> Dict[str, float]:
    summary = load_step008_ensemble_summary(step008)
    validate_screened_coupling(summary["beta_weighted"], summary["beta_eff_weighted"])
    return summary
