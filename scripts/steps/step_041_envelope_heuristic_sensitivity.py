#!/usr/bin/env python3
"""
Step 041: Monte Carlo sensitivity of Step 007 geometry-envelope heuristics.

Independent uniform ±50% multiplicative stress on each nominal coefficient
(``GEOMETRY_ENVELOPE_HEURISTIC_KEYS``), recomputing catalog predictions and
Step 008-style fits.  Emits leverage summaries and BIC stability fractions.

Monte Carlo throughput (many-core hosts such as M4 Pro):
- ``TEP_EFA_STEP041_MC_WORKERS`` (default ``min(10, max(2, cpu_count - 2))``): run draws in
  parallel processes, each with an independent ``numpy.random.SeedSequence([42, i])`` stream.
  Set to ``1`` to force the legacy single-threaded generator (slowest, but matches the
  original draw sequence bit-for-bit).
"""

from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.steps.step_007_tep_model import compute_tep_predictions_for_catalog
from scripts.steps.step_008_fitting import (
    analyze_fit_quality,
    bayesian_model_comparison,
    convert_to_native_types,
    fit_beta_to_observation,
)
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep_geometry_envelope import (
    GEOMETRY_ENVELOPE_HEURISTIC_KEYS,
    default_geometry_envelope_heuristics,
)

# Populated by ``ProcessPoolExecutor`` initializer in parallel Monte Carlo mode.
_MC_PLASMA = None


def _mc_pool_init() -> None:
    global _MC_PLASMA
    from scripts.steps.step_017_plasma_modulation import PlasmaModulationModel

    _MC_PLASMA = PlasmaModulationModel()


def _mc_pool_one(sample_index: int) -> dict:
    """One Monte Carlo draw (module-level for pickling under ``spawn``)."""
    if _MC_PLASMA is None:
        raise RuntimeError("Plasma model not initialized in worker")
    rng = np.random.default_rng(np.random.SeedSequence([42, int(sample_index)]))
    h = _random_heuristics(rng)
    pred = compute_tep_predictions_for_catalog(
        geometry_envelope_heuristics=h,
        plasma_model=_MC_PLASMA,
    )
    fits = _build_fits(pred)
    mc = bayesian_model_comparison(fits)
    row: dict = {"sample_index": int(sample_index), "heuristics": h}
    if mc.get("status") != "insufficient_data":
        row["best_model_bic"] = mc["model_comparison"]["best_model_bic"]
        row["best_model_bic_pessimistic_tep_vs_null"] = mc["model_comparison"].get(
            "best_model_bic_pessimistic_tep_vs_null"
        )
        q = analyze_fit_quality(fits)
        row["beta_weighted"] = float(q["beta_statistics"]["weighted_mean"])
    else:
        row["best_model_bic"] = None
        row["beta_weighted"] = None
        row["best_model_bic_pessimistic_tep_vs_null"] = None
    return row


def _build_fits(predictions: dict) -> dict:
    fits = {}
    for name, pred in predictions.items():
        fit_result = fit_beta_to_observation(pred, logger=None)
        fits[name] = {
            "spacecraft": pred["spacecraft"],
            "perigee": pred["perigee"],
            "observed": pred["observed"],
            "tep_predictions": pred["tep_predictions"],
            "cos_dec_asymmetry": pred.get("geometry", {}).get("cos_dec_asymmetry", 0.0),
            "fit": fit_result,
        }
    return fits


def _random_heuristics(rng: np.random.Generator) -> dict:
    base = default_geometry_envelope_heuristics()
    return {k: float(base[k] * rng.uniform(0.5, 1.5)) for k in GEOMETRY_ENVELOPE_HEURISTIC_KEYS}


def _one_at_a_time_ranges(plasma_model) -> dict:
    nominal = default_geometry_envelope_heuristics()
    out: dict = {}
    for key in GEOMETRY_ENVELOPE_HEURISTIC_KEYS:
        betas = []
        for factor in (0.5, 1.5):
            h = {**nominal, key: float(nominal[key] * factor)}
            pred = compute_tep_predictions_for_catalog(
                geometry_envelope_heuristics=h,
                plasma_model=plasma_model,
            )
            fits = _build_fits(pred)
            q = analyze_fit_quality(fits)
            if q.get("status") == "no_fits":
                betas.append(None)
            else:
                betas.append(q["beta_statistics"]["weighted_mean"])
        valid = [b for b in betas if b is not None]
        out[key] = {
            "beta_weighted_at_low": betas[0],
            "beta_weighted_at_high": betas[1],
            "beta_range": float(max(valid) - min(valid)) if len(valid) == 2 else None,
        }
    order = sorted(
        out.keys(),
        key=lambda k: (-(out[k]["beta_range"] or -1.0)),
    )
    return {"per_key": out, "leverage_rank_by_beta_range": order}


def main() -> int:
    logger = StepLogger("step_041_envelope_heuristic_sensitivity", PROJECT_ROOT)
    t0 = time.time()
    logger.header("STEP 041: GEOMETRY ENVELOPE HEURISTIC SENSITIVITY")

    cat_path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
    if not cat_path.exists():
        logger.error(f"Missing {cat_path}")
        logger.log_step_summary(time.time() - t0, "FAILED")
        return 1

    rng = np.random.default_rng(42)
    n_samples = 400

    from scripts.steps.step_017_plasma_modulation import PlasmaModulationModel

    shared_plasma = PlasmaModulationModel()

    nominal = default_geometry_envelope_heuristics()
    baseline_pred = compute_tep_predictions_for_catalog(plasma_model=shared_plasma)
    baseline_fits = _build_fits(baseline_pred)
    base_quality = analyze_fit_quality(baseline_fits)
    base_mc = bayesian_model_comparison(baseline_fits)

    oa = _one_at_a_time_ranges(shared_plasma)

    cpu_n = os.cpu_count() or 8
    default_mc_workers = max(2, min(10, cpu_n - 2))
    mc_workers = max(1, int(os.environ.get("TEP_EFA_STEP041_MC_WORKERS", str(default_mc_workers))))
    logger.info(f"Monte Carlo parallel workers: {mc_workers} (cpu_count={cpu_n})")

    records: list[dict] = []
    bic_wins = 0
    bic_pess_wins = 0
    usable = 0
    beta_samples: list[float] = []

    if mc_workers <= 1:
        for i in range(n_samples):
            h = _random_heuristics(rng)
            pred = compute_tep_predictions_for_catalog(
                geometry_envelope_heuristics=h,
                plasma_model=shared_plasma,
            )
            fits = _build_fits(pred)
            mc = bayesian_model_comparison(fits)
            row = {"sample_index": i, "heuristics": h}
            if mc.get("status") != "insufficient_data":
                usable += 1
                bmw = mc["model_comparison"]["best_model_bic"]
                row["best_model_bic"] = bmw
                row["best_model_bic_pessimistic_tep_vs_null"] = mc["model_comparison"].get(
                    "best_model_bic_pessimistic_tep_vs_null"
                )
                q = analyze_fit_quality(fits)
                bw = q["beta_statistics"]["weighted_mean"]
                row["beta_weighted"] = float(bw)
                beta_samples.append(float(bw))
                if bmw == "TEP":
                    bic_wins += 1
                if mc["model_comparison"].get("best_model_bic_pessimistic_tep_vs_null") == "TEP":
                    bic_pess_wins += 1
            else:
                row["best_model_bic"] = None
                row["beta_weighted"] = None
            records.append(row)
    else:
        with ProcessPoolExecutor(max_workers=mc_workers, initializer=_mc_pool_init) as pool:
            records = list(pool.map(_mc_pool_one, range(n_samples)))
        for row in records:
            if row.get("best_model_bic") is not None:
                usable += 1
                beta_samples.append(float(row["beta_weighted"]))
                if row["best_model_bic"] == "TEP":
                    bic_wins += 1
                if row.get("best_model_bic_pessimistic_tep_vs_null") == "TEP":
                    bic_pess_wins += 1

    beta_arr = np.array(beta_samples, dtype=float)
    payload = {
        "metadata": {
            "step": "041",
            "n_monte_carlo": n_samples,
            "random_seed": 42,
            "monte_carlo_parallel_workers": mc_workers,
            "monte_carlo_rng_scheme": (
                "single_generator_legacy" if mc_workers <= 1 else "seed_sequence_per_draw_[42,i]"
            ),
            "monte_carlo_usable_samples": usable,
            "nominal_geometry_envelope_heuristics": nominal,
        },
        "baseline": {
            "analyze_fit_quality_status": base_quality.get("status"),
            "n_fits": base_quality.get("n_fits"),
            "beta_weighted": base_quality.get("beta_statistics", {}).get("weighted_mean"),
            "model_comparison": base_mc,
        },
        "monte_carlo_summary": {
            "fraction_best_bic_tep": bic_wins / n_samples if n_samples else 0.0,
            "fraction_best_bic_tep_pessimistic_vs_null": bic_pess_wins / n_samples
            if n_samples
            else 0.0,
            "beta_weighted_std_monte_carlo": float(np.std(beta_arr)) if len(beta_arr) else None,
            "beta_weighted_min": float(np.min(beta_arr)) if len(beta_arr) else None,
            "beta_weighted_max": float(np.max(beta_arr)) if len(beta_arr) else None,
        },
        "one_at_a_time_leverage": oa,
        "monte_carlo_records": records,
    }

    out_path = PROJECT_ROOT / "results" / "step041_envelope_heuristic_sensitivity.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(convert_to_native_types(payload), f, indent=2)

    logger.success(f"Wrote {out_path}")
    logger.add_output_file(out_path, "Envelope heuristic sensitivity")
    logger.log_step_summary(time.time() - t0, "SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
