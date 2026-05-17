#!/usr/bin/env python3
"""
Step 039b: Juno Falsification Pathway
=====================================

Computes the β₀ thresholds at which the TEP model's Juno prediction
would drop below observational significance levels.

The Juno flyby (2013-10-09, altitude 817.4 km) is the sole raw-tension
case: the TEP model predicts a positive Δv at universal β while the
published bound is 0.00 ± 0.02 mm/s.

This script computes:
  1. β thresholds for 1σ, 2σ, and 3σ consistency with the null
  2. The tension metric at the current fitted β
  3. The β interval where TEP would be falsified by the Juno null

Output: results/step039b_juno_falsification_bounds.json
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.physics import R_EARTH
from scripts.utils.step_logger import StepLogger


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def compute_beta_threshold(dv_ref: float, beta_ref: float, target_dv: float) -> float:
    """
    Compute the β at which dv_tep = target_dv, given the reference.
    Scaling law: dv ∝ β^(3/4).
    """
    if dv_ref <= 0 or target_dv <= 0:
        return float("inf")
    return beta_ref * (target_dv / dv_ref) ** (4.0 / 3.0)


def main():
    logger = StepLogger("step_039b_juno_falsification_pathway", PROJECT_ROOT)
    logger.header("STEP 039B: JUNO FALSIFICATION PATHWAY")

    results_dir = PROJECT_ROOT / "results"
    out_path = results_dir / "step039b_juno_falsification_bounds.json"

    # ------------------------------------------------------------------
    # Load prerequisite results
    # ------------------------------------------------------------------
    step007 = load_json(results_dir / "step007_tep_predictions.json")
    step008 = load_json(results_dir / "step008_fitting_results.json")

    juno_pred = step007["predictions"]["Juno"]
    juno_tep = juno_pred["tep_predictions"]
    dv_tep_ref = juno_tep["dv_tep_mm_s"]  # at beta_ref = 1e-4
    beta_ref = juno_tep["beta_reference"]
    altitude_km = juno_pred["perigee"]["altitude_km"]
    r_perigee_m = R_EARTH + altitude_km * 1e3

    # Observed bound
    dv_obs = juno_pred["observed"]["dv_obs_mm_s"]  # 0.0
    sigma_obs = juno_pred["observed"]["sigma_mm_s"]  # 0.02

    # Current fitted β
    beta_stats = step008["overall_analysis"]["beta_statistics"]
    beta_weighted = beta_stats["weighted_mean"]
    beta_uncertainty = beta_stats.get("random_effects_uncertainty",
                                      beta_stats["weighted_uncertainty"])

    # ------------------------------------------------------------------
    # Compute falsification thresholds
    # ------------------------------------------------------------------
    scale_exponent = 0.75

    # At current β
    dv_current = dv_tep_ref * (beta_weighted / beta_ref) ** scale_exponent
    tension_sigma_current = dv_current / sigma_obs if sigma_obs > 0 else float("inf")

    # β thresholds for consistency with null at 1σ, 2σ, 3σ
    thresholds = {}
    for nsigma in [1, 2, 3]:
        target_dv = nsigma * sigma_obs
        beta_thresh = compute_beta_threshold(dv_tep_ref, beta_ref, target_dv)
        dv_at_thresh = dv_tep_ref * (beta_thresh / beta_ref) ** scale_exponent

        # Sign-gated threshold (conservative: require prediction < nsigma)
        falsified = beta_weighted > beta_thresh

        thresholds[f"{nsigma}sigma"] = {
            "target_dv_mm_s": round(target_dv, 4),
            "beta_threshold": round(beta_thresh, 6),
            "dv_at_threshold_mm_s": round(dv_at_thresh, 4),
            "current_beta_exceeds_threshold": bool(falsified),
            "tension_at_current_beta_sigma": round(tension_sigma_current, 2)
            if falsified else None,
        }

    # ------------------------------------------------------------------
    # β interval where TEP is falsified by Juno null
    # ------------------------------------------------------------------
    beta_falsified_upper = thresholds["3sigma"]["beta_threshold"]
    beta_falsified_fraction = beta_weighted / beta_falsified_upper

    # If current β is above the 3σ threshold, TEP is in tension
    tep_falsified_by_juno = beta_weighted > beta_falsified_upper

    # ------------------------------------------------------------------
    # Compute what β would be needed to explain the *upper bound* as TEP
    # (i.e., treat 0.02 mm/s as a 1σ TEP detection)
    # ------------------------------------------------------------------
    beta_explain_1sigma = compute_beta_threshold(dv_tep_ref, beta_ref, sigma_obs)
    beta_explain_2sigma = compute_beta_threshold(dv_tep_ref, beta_ref, 2 * sigma_obs)

    output = {
        "spacecraft": "Juno",
        "flyby_date": "2013-10-09",
        "altitude_km": altitude_km,
        "perigee_radius_m": round(r_perigee_m, 1),
        "observed": {
            "dv_obs_mm_s": dv_obs,
            "sigma_mm_s": sigma_obs,
        },
        "reference_prediction": {
            "dv_tep_ref_mm_s": round(dv_tep_ref, 6),
            "beta_ref": beta_ref,
        },
        "current_fitted_beta": {
            "beta_weighted": round(beta_weighted, 6),
            "beta_uncertainty": round(beta_uncertainty, 6),
            "dv_prediction_mm_s": round(dv_current, 4),
            "tension_sigma": round(tension_sigma_current, 2),
        },
        "falsification_thresholds": thresholds,
        "falsification_assessment": {
            "tep_falsified_by_juno_at_3sigma": bool(tep_falsified_by_juno),
            "beta_falsified_interval": {
                "upper_bound": round(beta_falsified_upper, 6),
                "current_beta_fraction_of_threshold": round(beta_falsified_fraction, 2),
            },
            "beta_needed_to_explain_observation": {
                "at_1sigma_mm_s": round(beta_explain_1sigma, 6),
                "at_2sigma_mm_s": round(beta_explain_2sigma, 6),
            },
            "interpretation": (
                "At the current fitted β, the TEP model predicts a Juno anomaly "
                f"{tension_sigma_current:.1f}σ above the published null. "
                f"To drop below 3σ consistency, β would need to be < {beta_falsified_upper:.4e}. "
                f"The current fitted β ({beta_weighted:.4e}) is "
                f"{beta_falsified_fraction:.1f}x the 3σ threshold."
            ),
        },
        "notes": [
            "Scaling law: dv_tep ∝ β^(3/4) (from Step 008 fitting).",
            "Falsification requires independent minimal-OD reanalysis with TEP-inclusive force modeling.",
            "OD survival factors (F_OD) are withheld until mission OD configuration data are available.",
            "These bounds assume the published null (0.00 ± 0.02 mm/s) is correct and systematic-free.",
        ],
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    logger.info(f"✓ Saved Juno falsification bounds to: {out_path}")
    logger.info(f"  Current β = {beta_weighted:.4e} → Juno prediction = {dv_current:.3f} mm/s")
    logger.info(f"  Tension = {tension_sigma_current:.1f}σ vs published null")
    logger.info(f"  3σ falsification threshold: β < {beta_falsified_upper:.4e}")
    logger.info(f"  TEP falsified by Juno at 3σ: {tep_falsified_by_juno}")

    logger.success("STEP 039B COMPLETE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
