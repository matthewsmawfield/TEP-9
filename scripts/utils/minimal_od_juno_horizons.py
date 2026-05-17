"""
Public-data minimal dynamics fit for Juno 2013 Earth flyby (no proprietary TRK).

Uses only repository-visible **JPL Horizons** ephemeris
(``data/raw/jpl_horizons/Juno_2013/Juno_2013_trajectory.json``) plus
``results/step038_3d_state_vectors.json`` for a perigee Cartesian anchor.

Observables (Horizons ``CENTER='Geocentric'`` OBS table):
  * ``range_m``  → geocentric distance |r|
  * ``velocity_m_s`` → geocentric range rate (``deldot``)

Two batch adjustments are reported:
  * **velocity_only** — match ``deldot`` alone (weaker geometry).
  * **range_and_velocity** — stacked residuals (stronger use of public data).

Dynamics: ``OrbitalMechanics3D`` from Step 012 (pure Kepler, point-mass Earth).
Estimator: ``scipy.optimize.least_squares`` (TRF) with per-state scaling.

This is **not** DSN two-way Doppler or MONTE-class OD. Step 030 merges the
returned structure as ``horizons_public_ephemeris_batch`` when inputs exist
unless ``TEP_030_SKIP_HORIZONS_PUBLIC_OD=1``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import least_squares

from scripts.steps.step_012_od_filter_simulation import OrbitalMechanics3D, PhysicalConstants


class GeocentricRangeRateObservable:
    """h_rdot = (v · r) / |r| (Horizons ``deldot`` for geocentric OBS)."""

    def __init__(self, noise_sigma: float):
        self.noise_sigma = float(noise_sigma)

    def compute_range_rate(self, sc_state: np.ndarray, t: float) -> float:
        del t
        r = sc_state[:3]
        v = sc_state[3:6]
        rm = float(np.linalg.norm(r))
        if rm < 1e4:
            return 0.0
        return float(np.dot(v, r) / rm)


_GEO = GeocentricRangeRateObservable(1.0)


def _rdot_series_from_states(states: np.ndarray) -> np.ndarray:
    return np.array(
        [_GEO.compute_range_rate(states[i], 0.0) for i in range(len(states))],
        dtype=float,
    )


def _propagate_arc(
    propagator: OrbitalMechanics3D, state0: np.ndarray, t_rel: np.ndarray
) -> np.ndarray:
    return propagator.propagate(state0, t_rel)


def _parse_ts(s: str) -> datetime:
    return datetime.strptime(s.strip(), "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _state_at_perigee_from_038(project_root: Path) -> Tuple[np.ndarray, int, str]:
    path = project_root / "results" / "step038_3d_state_vectors.json"
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {path}; run Step 038 before Horizons public OD block."
        )
    data = _load_json(path)
    juno = data.get("Juno")
    if not isinstance(juno, dict):
        raise KeyError("step038_3d_state_vectors.json has no 'Juno' block.")
    idx = int(juno["perigee_index"])
    rx = float(juno["rx_km"]) * 1000.0
    ry = float(juno["ry_km"]) * 1000.0
    rz = float(juno["rz_km"]) * 1000.0
    vx = float(juno["vx_km_s"]) * 1000.0
    vy = float(juno["vy_km_s"]) * 1000.0
    vz = float(juno["vz_km_s"]) * 1000.0
    state = np.array([rx, ry, rz, vx, vy, vz], dtype=float)
    dt_utc = str(juno.get("datetime_utc", ""))
    return state, idx, dt_utc


def _propagate_two_body_to_time(
    propagator: OrbitalMechanics3D,
    state_at_anchor: np.ndarray,
    t_anchor_s: float,
    t_target_s: float,
) -> np.ndarray:
    dt = float(t_target_s - t_anchor_s)
    if dt == 0.0:
        return state_at_anchor.copy()
    n = max(20, min(2000, int(abs(dt) / 30.0) + 2))
    t_arr = np.linspace(0.0, dt, n)
    return propagator.propagate(state_at_anchor, t_arr)[-1]


@dataclass(frozen=True)
class HorizonsMinimalODConfig:
    half_width_minutes: int = 180
    max_nfev: int = 500
    nominal_sigma_rdot_m_s: float = 50.0
    """Floor for range sigma (m); actual sigma uses max(floor, prefit std)."""
    nominal_sigma_range_floor_m: float = 5000.0


def _slice_arc(
    project_root: Path, cfg: HorizonsMinimalODConfig
) -> Tuple[Dict[str, Any], Path, np.ndarray, np.ndarray, Optional[np.ndarray], np.ndarray, float]:
    traj_path = project_root / "data" / "raw" / "jpl_horizons" / "Juno_2013" / "Juno_2013_trajectory.json"
    if not traj_path.is_file():
        raise FileNotFoundError(f"Missing Horizons trajectory: {traj_path}")

    traj = _load_json(traj_path)
    ts: List[str] = traj["timestamp"]
    v_los: List[float] = traj["velocity_m_s"]
    n = len(ts)
    if n != len(v_los):
        raise ValueError("Juno_2013_trajectory.json timestamp / velocity_m_s length mismatch.")

    range_m = traj.get("range_m")
    if range_m is not None:
        if len(range_m) != n:
            raise ValueError("range_m length mismatch vs timestamp in Juno_2013_trajectory.json")
        meas_r = np.array([float(range_m[i]) for i in range(n)], dtype=float)
    else:
        meas_r = None

    perigee_state, perigee_idx, perigee_dt = _state_at_perigee_from_038(project_root)
    lo = max(0, perigee_idx - cfg.half_width_minutes)
    hi = min(n - 1, perigee_idx + cfg.half_width_minutes)
    if hi - lo < 30:
        raise RuntimeError(
            f"Horizons arc slice too short after windowing: indices {lo}–{hi} (need ≥30 points)."
        )

    t_perigee = _parse_ts(ts[perigee_idx]).timestamp()
    t0 = _parse_ts(ts[lo]).timestamp()
    t_rel = np.array(
        [_parse_ts(ts[i]).timestamp() - t0 for i in range(lo, hi + 1)],
        dtype=float,
    )
    meas_v = np.array([v_los[i] for i in range(lo, hi + 1)], dtype=float)
    meas_range = meas_r[lo : hi + 1] if meas_r is not None else None

    meta = {
        "horizons_trajectory_file": str(traj_path.relative_to(project_root)),
        "step038_anchor": {
            "perigee_index": perigee_idx,
            "perigee_datetime_utc": perigee_dt,
            "arc_index_lo": lo,
            "arc_index_hi": hi,
            "half_width_minutes": cfg.half_width_minutes,
        },
    }
    return meta, traj_path, t_rel, meas_v, meas_range, perigee_state, float(t_perigee - t0)


def _trf_solve(
    residual_fn,
    x0: np.ndarray,
    x_scale: np.ndarray,
    max_nfev: int,
) -> Tuple[np.ndarray, Any]:
    lsq = least_squares(
        residual_fn,
        x0,
        method="trf",
        x_scale=x_scale,
        max_nfev=int(max_nfev),
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    return np.asarray(lsq.x, dtype=float), lsq


def _pack_velocity_only_result(
    lsq: Any,
    x_hat: np.ndarray,
    prefit: np.ndarray,
    post: np.ndarray,
    x_scale: np.ndarray,
    cfg: HorizonsMinimalODConfig,
) -> Dict[str, Any]:
    rms_pre = float(np.sqrt(np.mean(prefit**2)))
    rms_post = float(np.sqrt(np.mean(post**2)))
    return {
        "status": "computed" if lsq.success else "completed_non_success",
        "scipy_least_squares": {
            "success": bool(lsq.success),
            "message": str(lsq.message),
            "nfev": int(lsq.nfev),
            "cost": float(lsq.cost) if lsq.cost is not None else float("nan"),
            "method": "trf",
            "x_scale_m_m_s": x_scale.tolist(),
            "residuals_scaled_by_sigma_rdot_m_s": float(cfg.nominal_sigma_rdot_m_s),
        },
        "batch_least_squares": {
            "prefit_rms_rdot_residual_m_s": rms_pre,
            "prefit_rms_rdot_residual_mm_s": rms_pre * 1000.0,
            "post_fit_rms_rdot_residual_m_s": rms_post,
            "post_fit_rms_rdot_residual_mm_s": rms_post * 1000.0,
            "prefit_rms_m_s": rms_pre,
            "prefit_rms_mm_s": rms_pre * 1000.0,
            "post_fit_rms_m_s": rms_post,
            "post_fit_rms_mm_s": rms_post * 1000.0,
            "state_estimate_m_m_s": x_hat.tolist(),
            "residual_definition": (
                "RMS of (pure-Kepler geocentric deldot model − Horizons velocity_m_s); "
                "prefit/post_fit_rms_m_s duplicate prefit/post_fit_rms_rdot_residual_m_s for legacy readers."
            ),
        },
        "residuals_summary": {
            "n_obs": int(post.size),
            "mean_m_s": float(np.mean(post)),
            "median_abs_m_s": float(np.median(np.abs(post))),
        },
    }


def run_horizons_public_ephemeris_batch(
    project_root: Path,
    cfg: Optional[HorizonsMinimalODConfig] = None,
) -> Dict[str, Any]:
    """
    Best-effort public-data batch fits: Horizons ``deldot`` and, when ``range_m``
    is present, joint ``|r|`` + ``deldot`` vs pure Kepler.
    """
    cfg = cfg or HorizonsMinimalODConfig()
    meta, traj_path, t_rel, meas_v, meas_range, perigee_state, t_anchor_rel = _slice_arc(
        project_root, cfg
    )

    propagator = OrbitalMechanics3D(tep_model=None)
    state0_guess = _propagate_two_body_to_time(
        propagator, perigee_state, t_anchor_rel, t_rel[0]
    )
    x_scale = np.array([5e5, 5e5, 5e5, 5.0, 5.0, 5.0], dtype=float)

    # --- velocity-only ---
    def res_vel(x: np.ndarray) -> np.ndarray:
        st = _propagate_arc(propagator, x, t_rel)
        mod = _rdot_series_from_states(st)
        return (mod - meas_v) / float(cfg.nominal_sigma_rdot_m_s)

    st_pre = _propagate_arc(propagator, state0_guess, t_rel)
    pre_v = _rdot_series_from_states(st_pre) - meas_v
    x_v, lsq_v = _trf_solve(res_vel, state0_guess.copy(), x_scale, cfg.max_nfev)
    st_post_v = _propagate_arc(propagator, x_v, t_rel)
    post_v = _rdot_series_from_states(st_post_v) - meas_v

    velocity_block = _pack_velocity_only_result(
        lsq_v, x_v, pre_v, post_v, x_scale, cfg
    )
    velocity_block["interpretation"] = (
        "Match Horizons geocentric ``deldot`` only. Weakest public-data path; "
        "included for comparison with the range+velocity joint fit."
    )

    out: Dict[str, Any] = {
        "data_sources": {
            "horizons_trajectory": str(traj_path.relative_to(project_root)),
            "step038": "results/step038_3d_state_vectors.json (Juno block)",
        },
        "physics": {
            "mu_earth_m3_s2": PhysicalConstants.MU_EARTH,
            "dynamics": "OrbitalMechanics3D pure Kepler (Step 012)",
        },
        "arc": meta,
        "velocity_only": velocity_block,
    }

    # --- range + velocity (when range available) ---
    if meas_range is None:
        out["range_and_velocity"] = {
            "status": "skipped",
            "reason": "Juno_2013_trajectory.json has no range_m column",
        }
        out["status"] = velocity_block["status"]
        return out

    sigma_r = max(
        float(cfg.nominal_sigma_range_floor_m),
        float(np.std(np.linalg.norm(st_pre[:, :3], axis=1) - meas_range)),
    )
    sigma_v = float(cfg.nominal_sigma_rdot_m_s)

    def res_joint(x: np.ndarray) -> np.ndarray:
        st = _propagate_arc(propagator, x, t_rel)
        r_mod = np.linalg.norm(st[:, :3], axis=1)
        v_mod = _rdot_series_from_states(st)
        return np.concatenate(
            [(r_mod - meas_range) / sigma_r, (v_mod - meas_v) / sigma_v]
        )

    pre_r = np.linalg.norm(st_pre[:, :3], axis=1) - meas_range
    pre_joint = np.concatenate([pre_r / sigma_r, pre_v / sigma_v])
    cost_pre = float(np.dot(pre_joint, pre_joint))

    x_j, lsq_j = _trf_solve(res_joint, state0_guess.copy(), x_scale, cfg.max_nfev)
    st_j = _propagate_arc(propagator, x_j, t_rel)
    post_r = np.linalg.norm(st_j[:, :3], axis=1) - meas_range
    post_v2 = _rdot_series_from_states(st_j) - meas_v

    out["range_and_velocity"] = {
        "status": "computed" if lsq_j.success else "completed_non_success",
        "interpretation": (
            "Joint fit to Horizons geocentric ``range_m`` and ``velocity_m_s`` under "
            "pure Kepler. Uses all public ephemeris columns in the Step 004 JSON; "
            "still not TRK Doppler or empirical-acceleration OD."
        ),
        "sigmas_used_m": {
            "range_m": float(sigma_r),
            "rdot_m_s": float(sigma_v),
            "range_sigma_from": "max(floor, prefit std of |r_kepler|-|r_Horizons|)",
        },
        "scipy_least_squares": {
            "success": bool(lsq_j.success),
            "message": str(lsq_j.message),
            "nfev": int(lsq_j.nfev),
            "cost": float(lsq_j.cost) if lsq_j.cost is not None else float("nan"),
            "method": "trf",
            "prefit_cost_l2": cost_pre,
        },
        "batch_least_squares": {
            "state_estimate_m_m_s": x_j.tolist(),
            "post_fit_rms_range_m": float(np.sqrt(np.mean(post_r**2))),
            "post_fit_rms_rdot_m_s": float(np.sqrt(np.mean(post_v2**2))),
            "post_fit_rms_rdot_mm_s": float(np.sqrt(np.mean(post_v2**2)) * 1000.0),
            "prefit_rms_range_m": float(np.sqrt(np.mean(pre_r**2))),
            "prefit_rms_rdot_m_s": float(np.sqrt(np.mean(pre_v**2))),
        },
    }
    out["status"] = (
        "computed"
        if velocity_block["status"] == "computed"
        and out["range_and_velocity"]["status"] == "computed"
        else "partial"
    )
    return out


def run_horizons_los_minimal_od_batch(
    project_root: Path,
    cfg: Optional[HorizonsMinimalODConfig] = None,
) -> Dict[str, Any]:
    """Backward-compatible name: full public report (velocity + joint when possible)."""
    return run_horizons_public_ephemeris_batch(project_root, cfg)
