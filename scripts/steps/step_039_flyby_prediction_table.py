#!/usr/bin/env python3
"""
Step 039: Per-Flyby Post-OD Prediction and Classification Table

Computes a rigorous classification table for all 12 flybys using:
  1. Step 007 raw TEP predictions (at reference beta = 1e-4)
  2. Step 008 inverse-variance pooled beta (sign-gated restricted tier), with
     random-effects scatter used for prediction uncertainty when available
  3. Step 021 OD survival factors (F_OD), emitted only when mission-specific
     OD configuration data provide defensible values
  4. Algorithmic raw-layer classification logic with falsifiability criterion

Output: results/step039_flyby_prediction_table.json
"""

import datetime
import json
import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _detection_flags(
    observed_mm_s,
    sigma_mm_s,
    raw_mm_s,
    post_od_mm_s=None,
    raw_sigma_mm_s=None,
    *,
    uncertainty_aware_prediction: bool = False,
):
    if sigma_mm_s is None or sigma_mm_s <= 0:
        return None
    if observed_mm_s is None or raw_mm_s is None:
        return None

    threshold = 3.0 * sigma_mm_s
    raw_threshold = threshold
    if uncertainty_aware_prediction:
        raw_sigma = float(raw_sigma_mm_s or 0.0)
        raw_threshold = 3.0 * float(np.sqrt(float(sigma_mm_s) ** 2 + raw_sigma**2))
    return {
        "threshold_mm_s": threshold,
        "raw_prediction_threshold_mm_s": raw_threshold,
        "obs_detected": abs(observed_mm_s) > threshold,
        "raw_predicted": abs(raw_mm_s) > raw_threshold,
        "post_predicted": (
            abs(post_od_mm_s) > threshold if post_od_mm_s is not None else None
        ),
    }


def classify_raw_pooled_beta(observed_mm_s, sigma_mm_s, raw_mm_s):
    """
    Raw-layer classification at fixed Step 008 pooled amplitude (no OD survival applied).

    Classifications:
      true_positive      : anomaly observed AND raw prediction exceeds threshold
      true_null          : null observed AND raw prediction below threshold
      raw_tension        : null observed BUT raw prediction exceeds threshold
      raw_surplus        : anomaly observed BUT raw prediction below threshold
    """
    flags = _detection_flags(observed_mm_s, sigma_mm_s, raw_mm_s)
    if flags is None:
        return "uncertainty_unavailable"
    if flags["obs_detected"] and flags["raw_predicted"]:
        return "true_positive"
    if not flags["obs_detected"] and not flags["raw_predicted"]:
        return "true_null"
    if not flags["obs_detected"] and flags["raw_predicted"]:
        return "raw_tension"
    if flags["obs_detected"] and not flags["raw_predicted"]:
        return "raw_surplus"
    return "unclassified"


def classify_raw_pooled_beta_uncertainty_aware(observed_mm_s, sigma_mm_s, raw_mm_s, raw_sigma_mm_s):
    """
    Raw-layer classification that treats the Step 008 amplitude scatter as
    prediction uncertainty. This is the falsification/stress-test label.
    """
    flags = _detection_flags(
        observed_mm_s,
        sigma_mm_s,
        raw_mm_s,
        raw_sigma_mm_s=raw_sigma_mm_s,
        uncertainty_aware_prediction=True,
    )
    if flags is None:
        return "uncertainty_unavailable"
    if flags["obs_detected"] and flags["raw_predicted"]:
        return "true_positive"
    if not flags["obs_detected"] and not flags["raw_predicted"]:
        return "true_null"
    if not flags["obs_detected"] and flags["raw_predicted"]:
        return "raw_tension"
    if flags["obs_detected"] and not flags["raw_predicted"]:
        return "observed_detection_prediction_uncertain"
    return "unclassified"


def classify_flyby(observed_mm_s, sigma_mm_s, raw_mm_s, post_od_mm_s):
    """
    Post-OD classification with 3-sigma threshold.

    Classifications:
      true_positive      : anomaly observed AND post-OD predicts anomaly
      true_null          : null observed AND post-OD predicts null
      suppressed_positive: null observed BUT raw predicts anomaly AND post-OD predicts null
      false_positive     : anomaly observed BUT post-OD predicts null
      false_negative     : null observed BUT post-OD predicts anomaly
    """
    flags = _detection_flags(observed_mm_s, sigma_mm_s, raw_mm_s, post_od_mm_s)
    if flags is None:
        return "uncertainty_unavailable"
    if post_od_mm_s is None or flags["post_predicted"] is None:
        return "unclassified"

    if flags["obs_detected"] and flags["post_predicted"]:
        return "true_positive"
    if not flags["obs_detected"] and not flags["post_predicted"]:
        return "true_null"
    if not flags["obs_detected"] and flags["raw_predicted"] and not flags["post_predicted"]:
        return "suppressed_positive"
    if flags["obs_detected"] and not flags["post_predicted"]:
        return "false_positive"
    if not flags["obs_detected"] and flags["post_predicted"]:
        return "false_negative"
    return "unclassified"


def compute_full_catalog_raw_likelihood(rows):
    """
    Compute a raw-layer full-catalog stress-test likelihood.

    This is deliberately not a replacement for mission-specific OD reanalysis.
    It compares the null model against the Step 008 pooled-amplitude raw TEP prediction
    for rows with published observations and explicit raw TEP predictions.
    Uncertainty combines the published observational uncertainty and the
    propagated pooled-amplitude prediction uncertainty.
    """
    included = []
    for row in rows:
        observed = row.get("observed_dv_mm_s")
        sigma_obs = row.get("observed_sigma_mm_s")
        raw = row.get("raw_tep_prediction_mm_s")
        sigma_raw = row.get("raw_tep_uncertainty_mm_s")
        if observed is None or sigma_obs is None or raw is None:
            continue
        sigma_raw = sigma_raw or 0.0
        sigma_total = float(np.sqrt(float(sigma_obs) ** 2 + float(sigma_raw) ** 2))
        if sigma_total <= 0:
            continue
        included.append(
            {
                "flyby": row["flyby"],
                "observed": float(observed),
                "raw_prediction": float(raw),
                "sigma_total": sigma_total,
                "tep_residual": float(observed - raw),
                "null_residual": float(observed),
            }
        )

    def _loglike(residual, sigma):
        return -0.5 * (residual / sigma) ** 2 - 0.5 * np.log(2 * np.pi * sigma**2)

    log_l_null = sum(_loglike(row["null_residual"], row["sigma_total"]) for row in included)
    log_l_tep = sum(_loglike(row["tep_residual"], row["sigma_total"]) for row in included)
    chi2_null = sum((row["null_residual"] / row["sigma_total"]) ** 2 for row in included)
    chi2_tep = sum((row["tep_residual"] / row["sigma_total"]) ** 2 for row in included)

    return {
        "description": (
            "Raw fixed-amplitude (Step 008 pooled beta) full-catalog stress test over published rows with "
            "explicit raw TEP predictions; excludes geometry-unavailable Rosetta 2009 "
            "and no-public-report flybys."
        ),
        "n_rows": len(included),
        "included_flybys": [row["flyby"] for row in included],
        "log_likelihood_null": float(log_l_null),
        "log_likelihood_tep_raw": float(log_l_tep),
        "delta_log_likelihood_tep_minus_null": float(log_l_tep - log_l_null),
        "chi2_null": float(chi2_null),
        "chi2_tep_raw": float(chi2_tep),
        "delta_chi2_null_minus_tep": float(chi2_null - chi2_tep),
        "notes": [
            "This is a raw-layer stress test, not a post-OD mission likelihood.",
            "Juno and Cassini remain explicit stress cases in this likelihood.",
            "Prediction uncertainty is propagated from Step 008; random-effects amplitude scatter is used when available."
        ],
        "per_flyby": included,
    }


def main():
    logger = StepLogger("step_039_flyby_prediction_table", PROJECT_ROOT)
    logger.header("STEP 039: PER-FLYBY POST-OD PREDICTION AND CLASSIFICATION TABLE")

    results_dir = PROJECT_ROOT / "results"
    out_path = results_dir / "step039_flyby_prediction_table.json"

    # ------------------------------------------------------------------
    # Load prerequisite results
    # ------------------------------------------------------------------
    step007 = load_json(results_dir / "step007_tep_predictions.json")
    step008 = load_json(results_dir / "step008_fitting_results.json")
    step021 = load_json(results_dir / "step021_od_simulation_validation.json")
    catalog = load_json(results_dir / "step003_archival_flyby_catalog.json")

    predictions = step007.get("predictions", {})
    fits = step008.get("individual_fits", {})

    beta_stats = step008["overall_analysis"]["beta_statistics"]
    beta_weighted = beta_stats["weighted_mean"]
    beta_formal_uncertainty = beta_stats["weighted_uncertainty"]
    beta_uncertainty = beta_stats.get("random_effects_uncertainty", beta_formal_uncertainty)
    beta_uncertainty_model = (
        "random_effects_uncertainty"
        if beta_stats.get("random_effects_uncertainty") is not None
        else "weighted_uncertainty"
    )
    beta_initial = 1e-4

    # Scaling exponent from Step 008 (3/4 power law)
    scale_factor = (beta_weighted / beta_initial) ** 0.75
    rel_unc_beta = beta_uncertainty / beta_weighted
    rel_unc_raw = 0.75 * rel_unc_beta  # from power-law error propagation

    logger.info(
        f"Step 008 pooled beta = {beta_weighted:.6e}; prediction uncertainty "
        f"uses {beta_uncertainty_model} = {beta_uncertainty:.6e}"
    )
    logger.info(f"Scale factor (beta^(3/4)) = {scale_factor:.4f}")
    logger.info(f"Relative uncertainty on raw prediction = {rel_unc_raw:.4f}")

    # ------------------------------------------------------------------
    # Build lookup maps
    # ------------------------------------------------------------------
    catalog_by_name = {entry["mission_name"]: entry for entry in catalog["flybys"]}
    fod_results = step021.get("results", {})

    def get_fod(mission_name: str):
        """Retrieve F_OD from step021 when a real value exists."""
        entry = fod_results.get(mission_name)
        if isinstance(entry, dict) and entry.get("f_od_estimate") is not None:
            return entry
        return None

    def infer_data_class(entry: dict) -> str:
        """Classify catalog provenance for manuscript table regeneration."""
        observed = entry.get("published_anomaly_mm_s")
        if observed is None:
            return "No public anomaly report"
        sigma = entry.get("published_anomaly_uncertainty_mm_s")
        if sigma is None:
            sigma = entry.get("tracking_precision_mm_s")
        if observed == 0 or (sigma is not None and abs(observed) < 2 * sigma):
            return "Published null/bound"
        return "Published anomaly"

    # ------------------------------------------------------------------
    # Ordered list of the 12 flybys as they appear in the manuscript
    # ------------------------------------------------------------------
    flyby_order = [
        ("NEAR", "NEAR"),
        ("Galileo 1990", "Galileo_1990"),
        ("Rosetta 2005", "Rosetta_2005"),
        ("Cassini", "Cassini"),
        ("Galileo 1992", "Galileo_1992"),
        ("MESSENGER", "MESSENGER_2005"),
        ("Rosetta 2009", "Rosetta_2009"),
        ("Juno", "Juno"),
        ("Rosetta 2007", "Rosetta_2007"),
        ("Stardust", "Stardust"),
        ("OSIRIS-REx", "OSIRIS-REx"),
        ("BepiColombo", "BepiColombo"),
    ]

    rows = []
    classification_counts = {
        "true_positive": 0,
        "true_null": 0,
        "suppressed_positive": 0,
        "false_positive": 0,
        "false_negative": 0,
    }
    raw_classification_counts = {
        "true_positive": 0,
        "true_null": 0,
        "raw_tension": 0,
        "raw_surplus": 0,
    }
    raw_uncertainty_aware_counts = {
        "true_positive": 0,
        "true_null": 0,
        "raw_tension": 0,
        "raw_surplus": 0,
        "observed_detection_prediction_uncertain": 0,
    }

    for display_name, key in flyby_order:
        cat = catalog_by_name.get(key, {})
        pred = predictions.get(key, {})
        fit = fits.get(key, {})

        observed = cat.get("published_anomaly_mm_s")
        sigma = cat.get("published_anomaly_uncertainty_mm_s")
        if sigma is None:
            sigma = cat.get("tracking_precision_mm_s")

        altitude_km = cat.get("perigee_altitude_km")
        flyby_date = cat.get("flyby_date", "")
        data_class = infer_data_class(cat)

        if pred:
            dv_tep_ref = pred["tep_predictions"]["dv_tep_mm_s"]
            dv_grad = pred["tep_predictions"]["dv_grad_mm_s"]
            dv_disf = pred["tep_predictions"]["dv_disf_mm_s"]
            cos_asym = pred.get("geometry", {}).get("cos_dec_asymmetry")

            # Scale reference prediction to Step 008 pooled amplitude (beta^(3/4) law)
            dv_raw = dv_tep_ref * scale_factor
            sigma_raw = abs(dv_raw) * rel_unc_raw

            raw_classification = classify_raw_pooled_beta(observed, sigma, dv_raw)
            raw_classification_uncertainty_aware = classify_raw_pooled_beta_uncertainty_aware(
                observed, sigma, dv_raw, sigma_raw
            )
            pooled_beta_residual = (
                round(observed - dv_raw, 3)
                if observed is not None and dv_raw is not None
                else None
            )

            # Apply F_OD only when mission OD configuration produced a value
            fod_data = get_fod(key)
            if fod_data is None:
                dv_post_od = None
                sigma_post_od = None
                f_od = None
                sigma_fod = None
                classification = "f_od_unavailable"
                notes = (
                    "Post-OD classification withheld until mission OD configuration "
                    "yields F_OD. Raw pooled-amplitude classification is authoritative."
                )
            else:
                f_od = fod_data["f_od_estimate"]
                sigma_fod = fod_data["f_od_uncertainty"]
                dv_post_od = dv_raw * f_od
                rel_unc_fod = sigma_fod / f_od if f_od > 0 else 1.0
                sigma_post_od = abs(dv_post_od) * np.sqrt(rel_unc_raw**2 + rel_unc_fod**2)
                classification = classify_flyby(observed, sigma, dv_raw, dv_post_od)
                notes = ""

            row = {
                "flyby": display_name,
                "data_class": data_class,
                "altitude_km": altitude_km,
                "cos_asymmetry": cos_asym,
                "observed_dv_mm_s": observed,
                "observed_sigma_mm_s": sigma,
                "raw_tep_prediction_mm_s": round(dv_raw, 3) if dv_raw is not None else None,
                "raw_tep_uncertainty_mm_s": round(sigma_raw, 3) if sigma_raw is not None else None,
                "pooled_beta_residual_mm_s": pooled_beta_residual,
                "raw_classification": raw_classification,
                "raw_classification_uncertainty_aware": raw_classification_uncertainty_aware,
                "f_od": f_od,
                "f_od_uncertainty": sigma_fod,
                "post_od_prediction_mm_s": round(dv_post_od, 3) if dv_post_od is not None else None,
                "post_od_uncertainty_mm_s": round(sigma_post_od, 3) if sigma_post_od is not None else None,
                "classification": classification,
                "tep_geometry_available": True,
                "notes": notes,
            }
        else:
            # No TEP prediction available (insufficient geometry data)
            # For very high altitude flybys (>5000 km), classify as true null
            # because field gradient is negligible regardless of trajectory
            if altitude_km and altitude_km > 5000:
                classification = "true_null"
                notes = "High altitude; negligible TEP signal expected. Insufficient declination data for explicit prediction."
            elif observed == 0:
                classification = "insufficient_geometry_published_null"
                notes = "Published null/bound; insufficient geometry data for explicit pooled-amplitude TEP prediction."
            else:
                classification = "insufficient_data"
                notes = "Insufficient geometry data for TEP prediction."

            row = {
                "flyby": display_name,
                "data_class": data_class,
                "altitude_km": altitude_km,
                "cos_asymmetry": None,
                "observed_dv_mm_s": observed,
                "observed_sigma_mm_s": sigma,
                "raw_tep_prediction_mm_s": None,
                "raw_tep_uncertainty_mm_s": None,
                "f_od": None,
                "f_od_uncertainty": None,
                "post_od_prediction_mm_s": None,
                "post_od_uncertainty_mm_s": None,
                "classification": classification,
                "tep_geometry_available": False,
                "notes": notes,
            }

        if classification in classification_counts:
            classification_counts[classification] += 1
        elif classification not in (
            "f_od_unavailable",
            "uncertainty_unavailable",
            "unclassified",
            "insufficient_data",
            "insufficient_geometry_published_null",
        ):
            classification_counts[classification] = classification_counts.get(classification, 0) + 1

        raw_classification = row.get("raw_classification")
        if raw_classification in raw_classification_counts:
            raw_classification_counts[raw_classification] += 1
        raw_unc = row.get("raw_classification_uncertainty_aware")
        if raw_unc in raw_uncertainty_aware_counts:
            raw_uncertainty_aware_counts[raw_unc] += 1

        rows.append(row)
        logger.info(
            f"{display_name:20s} | obs={observed!s:>8} | raw={row.get('raw_tep_prediction_mm_s')!s:>8} | "
            f"raw_cls={raw_classification!s:>14} | raw_unc={raw_unc!s:>14} | FOD={row.get('f_od')!s:>6} | "
            f"post={row.get('post_od_prediction_mm_s')!s:>8} | {classification}"
        )

    # ------------------------------------------------------------------
    # Build output
    # ------------------------------------------------------------------
    output = {
        "metadata": {
            "step": "039_flyby_prediction_table",
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "beta_pooled_step008": beta_weighted,
            "beta_uncertainty": beta_uncertainty,
            "beta_formal_weighted_uncertainty": beta_formal_uncertainty,
            "beta_prediction_uncertainty_model": beta_uncertainty_model,
            "scale_factor": scale_factor,
            "rel_unc_raw": rel_unc_raw,
            "n_flybys": len(rows),
            "classification_summary": classification_counts,
            "raw_classification_summary": raw_classification_counts,
            "raw_uncertainty_aware_classification_summary": raw_uncertainty_aware_counts,
            "f_od_policy": (
                "Post-OD columns are emitted only when Step 021 supplies mission-specific "
                "F_OD from real OD configuration data. Raw pooled-amplitude classification "
                "uses the Step 008 inverse-variance mean scaled prediction without OD survival factors."
            ),
            "geometry_envelope": "step007 v5.4 deterministic multipole/plasma/symmetry envelope",
        },
        "falsifiability_criterion": {
            "description": (
                "Uncertainty-aware raw pooled-amplitude tension cases (observed null where "
                "the scaled TEP prediction exceeds the 3-sigma combined observation+prediction "
                "threshold) define model stress tests independent of OD survival factors. "
                "The deterministic raw_classification column is retained only as a fixed-amplitude "
                "diagnostic. Post-OD false negatives are counted only when mission-specific F_OD is available."
            ),
            "raw_tensions_found": raw_uncertainty_aware_counts["raw_tension"],
            "fixed_amplitude_raw_tensions_found": raw_classification_counts["raw_tension"],
            "false_negatives_found": classification_counts["false_negative"],
            "threshold_sigma": 3.0,
        },
        "full_catalog_raw_likelihood": compute_full_catalog_raw_likelihood(rows),
        "rows": rows,
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    logger.info(f"")
    logger.info(f"Raw pooled-amplitude classification summary:")
    for cls, count in raw_classification_counts.items():
        logger.info(f"  {cls:25s}: {count}")
    logger.info(f"")
    logger.info(f"Uncertainty-aware raw classification summary:")
    for cls, count in raw_uncertainty_aware_counts.items():
        logger.info(f"  {cls:25s}: {count}")
    logger.info(f"")
    logger.info(f"Post-OD classification summary:")
    for cls, count in classification_counts.items():
        logger.info(f"  {cls:25s}: {count}")
    logger.info(f"")
    logger.info(f"Output written to: {out_path}")
    logger.log_step_summary(len(rows), "SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
