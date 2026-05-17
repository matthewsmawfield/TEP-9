"""
Step 042: Time-Resolved Cosmographic Geometry Along Flyby Arcs

Re-evaluates the Step 040 CMB / solar scalar modulation construction at every
JPL Horizons epoch on each flyby arc (from step038_trajectory_series.json),
using the same Earth-barycentric + geocentric-SC velocity composition as Step 040.

Each epoch also records Earth barycentric position in ICRS (AU) from Astropy and,
when Step 038 supplies geocentric SC Cartesian position (km), the projection of
that position onto the CMB dipole apex direction (J2000 equatorial), for mapping
Earth and spacecraft trajectories relative to the dipole frame.

Along-arc time derivatives (np.gradient vs UTC) of selected CMB-related scalars
are summarized under ``arc_temporal_derivative_geometry`` (geometry-only; high |r|
can follow from coupled kinematics).

Published flyby anomalies are often single scalars per mission. Optional per-mission
residual time series may be supplied under ``data/time_resolved_flyby_residuals/<Mission>.json``
(see ``schema.json``); when present and valid, Step 042 aligns them to Horizons epochs
and reports Pearson correlations. When there are enough consecutive matched minutes,
``sequential_derivative_correlations`` compares d(residual)/dt to d(geometry)/dt
(Doppler-like structure along the arc; requires time-varying ``value_mm_s``).
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

try:
    from scripts.utils.step_logger import StepLogger
except ImportError:

    class StepLogger:
        def __init__(self, *args, **kwargs):
            pass

        def header(self, s):
            print(f"\n=== {s} ===")

        def section(self, s):
            print(f"\n--- {s} ---")

        def subsection(self, s):
            print(f"\n--- {s} ---")

        def info(self, s):
            print(s)

        def warning(self, s):
            print(f"WARNING: {s}")

        def log_step_summary(self, *args, **kwargs):
            pass


def _load_step040_geometry_module():
    path = PROJECT_ROOT / "scripts" / "steps" / "step_040_cosmographic_shear.py"
    spec = importlib.util.spec_from_file_location("tep_step040_geometry", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load Step 040 geometry module")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _utc_iso_to_datetime(iso: str) -> datetime:
    if iso.endswith("Z"):
        return datetime.fromisoformat(iso.replace("Z", "+00:00"))
    return datetime.fromisoformat(iso)


def _pearson_along_arc(x: List[float], y: List[float]) -> Tuple[Optional[float], Optional[float], int]:
    xc: List[float] = []
    yc: List[float] = []
    for a, b in zip(x, y):
        if math.isfinite(a) and math.isfinite(b):
            xc.append(a)
            yc.append(b)
    n = len(xc)
    if n < 3:
        return None, None, n
    r, p = stats.pearsonr(xc, yc)
    return float(r), float(p), n


def _both_aligned_flag(sc_cmb_cos_theta: float, earth_cmb_proj_kms: float) -> float:
    earth_proj = earth_cmb_proj_kms / 30.0
    return 1.0 if (sc_cmb_cos_theta > 0.0 and earth_proj > 0.0) else 0.0


def _to_json_serializable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _to_json_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_json_serializable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, float)):
        return float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, bool) or obj is None:
        return obj
    if isinstance(obj, (str,)):
        return obj
    return obj


def _compute_arc_epochs(
    traj: Dict[str, Any],
    g40: Any,
    logger: StepLogger,
) -> List[Dict[str, Any]]:
    n = int(traj["n_points"])
    hours = [float(h) for h in traj["hours_from_perigee"]]
    range_km = [float(r) for r in traj["range_km"]]
    altitude_km = [float(a) for a in traj["altitude_km"]]
    utc_iso = list(traj["utc_iso"])

    vx_eq = [float(v) for v in traj["vx_km_s"]]
    vy_eq = [float(v) for v in traj["vy_km_s"]]
    vz_eq = [float(v) for v in traj["vz_km_s"]]
    vx_ec = [float(v) for v in traj["vx_ecl_km_s"]]
    vy_ec = [float(v) for v in traj["vy_ecl_km_s"]]
    vz_ec = [float(v) for v in traj["vz_ecl_km_s"]]

    n_cmb = g40.unit_vector_equatorial(g40.CMB_DIPOLE_RA_DEG, g40.CMB_DIPOLE_DEC_DEG)
    sc_rx = traj.get("sc_rx_km")
    sc_ry = traj.get("sc_ry_km")
    sc_rz = traj.get("sc_rz_km")
    has_sc_pos = (
        isinstance(sc_rx, list)
        and isinstance(sc_ry, list)
        and isinstance(sc_rz, list)
        and len(sc_rx) == n
        and len(sc_ry) == n
        and len(sc_rz) == n
    )

    epochs: List[Dict[str, Any]] = []
    for i in range(n):
        dt = _utc_iso_to_datetime(utc_iso[i])
        ephem = g40.earth_sun_ephemeris_state(dt, logger, quiet=True)
        r_au = float(ephem["r_au"])
        lon_deg = float(ephem["ecliptic_longitude_deg"])
        v_earth_ecl = np.asarray(ephem["earth_helio_vel_kms_ecliptic"], dtype=float)
        v_earth_equ = np.asarray(ephem["earth_barycentric_vel_kms_icrs"], dtype=float)

        v_sc_equ = np.array([vx_eq[i], vy_eq[i], vz_eq[i]], dtype=float)
        v_sc_ecl = np.array([vx_ec[i], vy_ec[i], vz_ec[i]], dtype=float)

        solar_mod = g40.compute_solar_scalar_modulation(r_au, lon_deg, v_earth_ecl, v_sc_ecl)
        cmb_mod = g40.compute_cmb_dipole_modulation(v_earth_equ, v_sc_equ)

        sc_cos = float(cmb_mod["sc_cmb_cos_theta"])
        earth_cmb = float(cmb_mod["earth_orbital_cmb_proj_kms"])
        enh = float(cmb_mod["cmb_disformal_enhancement"])

        pos_e_icrs_au = np.asarray(ephem["earth_barycentric_pos_au_icrs"], dtype=float).reshape(3)
        cmb_map = {
            "earth_barycentric_pos_icrs_au_x": float(pos_e_icrs_au[0]),
            "earth_barycentric_pos_icrs_au_y": float(pos_e_icrs_au[1]),
            "earth_barycentric_pos_icrs_au_z": float(pos_e_icrs_au[2]),
            "earth_barycentric_cmb_parallel_au": float(np.dot(pos_e_icrs_au, n_cmb)),
        }
        if has_sc_pos:
            r_sc = np.array(
                [float(sc_rx[i]), float(sc_ry[i]), float(sc_rz[i])],
                dtype=float,
            )
            cmb_map["sc_geocentric_cmb_parallel_km"] = float(np.dot(r_sc, n_cmb))

        row = {
            "utc_iso": utc_iso[i],
            "hours_from_perigee": hours[i],
            "range_km": range_km[i],
            "altitude_km": altitude_km[i],
            "heliocentric_distance_au": r_au,
            "ecliptic_longitude_deg": lon_deg,
            "earth_sun_ephemeris": ephem["ephemeris"],
            **cmb_map,
            **solar_mod,
            **cmb_mod,
            "both_aligned_flag": _both_aligned_flag(sc_cos, earth_cmb),
        }
        row = _to_json_serializable(row)
        epochs.append(row)

    return epochs


def _arc_temporal_derivative_geometry(epochs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Along-arc time derivatives via np.gradient against UTC seconds from arc start.
    Pairs are geometry-only (Doppler-like structure in CMB projections vs. range-driven factors).
    """
    if len(epochs) < 5:
        return []
    t0 = _utc_iso_to_datetime(epochs[0]["utc_iso"])
    t_sec = np.array(
        [(_utc_iso_to_datetime(e["utc_iso"]) - t0).total_seconds() for e in epochs],
        dtype=float,
    )
    if np.any(np.diff(t_sec) <= 0):
        return []

    specs: List[Tuple[str, str, str]] = [
        (
            "total_velocity_cmb_proj_kms",
            "cmb_disformal_enhancement",
            "dt_total_cmb_proj_vs_dt_cmb_disformal_enhancement",
        ),
        (
            "sc_velocity_cmb_proj_kms",
            "cmb_modulation_factor",
            "dt_sc_cmb_proj_vs_dt_cmb_modulation_factor",
        ),
        (
            "total_velocity_cmb_proj_kms",
            "range_km",
            "dt_total_cmb_proj_vs_dt_range_km",
        ),
    ]
    if all("sc_geocentric_cmb_parallel_km" in e for e in epochs):
        specs.append(
            (
                "sc_velocity_cmb_proj_kms",
                "sc_geocentric_cmb_parallel_km",
                "dt_sc_cmb_proj_vs_dt_sc_geocentric_cmb_parallel_km",
            )
        )
    if all("earth_barycentric_cmb_parallel_au" in e for e in epochs):
        specs.append(
            (
                "earth_orbital_cmb_proj_kms",
                "earth_barycentric_cmb_parallel_au",
                "dt_earth_cmb_proj_vs_dt_earth_barycentric_cmb_parallel_au",
            )
        )

    out: List[Dict[str, Any]] = []
    for ka, kb, label in specs:
        try:
            ya = np.array([float(e[ka]) for e in epochs], dtype=float)
            yb = np.array([float(e[kb]) for e in epochs], dtype=float)
        except (KeyError, TypeError, ValueError):
            continue
        if not np.all(np.isfinite(ya)) or not np.all(np.isfinite(yb)):
            continue
        ga = np.gradient(ya, t_sec)
        gb = np.gradient(yb, t_sec)
        rp, pp, nn = _pearson_along_arc(ga.tolist(), gb.tolist())
        out.append({"test": label, "key_a": ka, "key_b": kb, "pearson_r": rp, "p_value": pp, "n": nn})
    return out


def _residual_sequential_derivative_correlations(
    idxs: List[int],
    y_res: List[float],
    epochs: List[Dict[str, Any]],
) -> Tuple[int, List[Dict[str, Any]]]:
    """
    For matched Horizons epochs in chronological order, use pairs where geometry
    indices differ by exactly 1 (consecutive 1-minute samples). Relates d(residual)/dt
    to d(geometry)/dt along the real trajectory (requires time-varying sidecar values).
    """
    keys: List[Tuple[str, str]] = [
        ("cmb_disformal_enhancement", "d_residual_dt_vs_dt_cmb_disformal_enhancement"),
        ("total_velocity_cmb_proj_kms", "d_residual_dt_vs_dt_total_velocity_cmb_proj_kms"),
        ("sc_velocity_cmb_proj_kms", "d_residual_dt_vs_dt_sc_velocity_cmb_proj_kms"),
        ("sc_geocentric_cmb_parallel_km", "d_residual_dt_vs_dt_sc_geocentric_cmb_parallel_km"),
        ("cmb_modulation_factor", "d_residual_dt_vs_dt_cmb_modulation_factor"),
    ]
    out: List[Dict[str, Any]] = []
    drv_res: List[float] = []
    geom_series: Dict[str, List[float]] = {g: [] for g, _ in keys}

    for k in range(len(idxs) - 1):
        j0, j1 = idxs[k], idxs[k + 1]
        if j1 != j0 + 1:
            continue
        t_a = _utc_iso_to_datetime(epochs[j0]["utc_iso"])
        t_b = _utc_iso_to_datetime(epochs[j1]["utc_iso"])
        dt = (t_b - t_a).total_seconds()
        if dt <= 0:
            continue
        drv_res.append((float(y_res[k + 1]) - float(y_res[k])) / dt)
        for gkey, _ in keys:
            v0 = epochs[j0].get(gkey)
            v1 = epochs[j1].get(gkey)
            if v0 is None or v1 is None:
                geom_series[gkey].append(float("nan"))
            else:
                fv0, fv1 = float(v0), float(v1)
                if not (math.isfinite(fv0) and math.isfinite(fv1)):
                    geom_series[gkey].append(float("nan"))
                else:
                    geom_series[gkey].append((fv1 - fv0) / dt)

    n_pairs = len(drv_res)
    if n_pairs < 5:
        return n_pairs, []

    for gkey, label in keys:
        ys = geom_series[gkey]
        xs: List[float] = []
        yc: List[float] = []
        for a, b in zip(drv_res, ys):
            if math.isfinite(a) and math.isfinite(b):
                xs.append(a)
                yc.append(b)
        rp, pp, nn = _pearson_along_arc(xs, yc)
        out.append(
            {
                "test": label,
                "geometry_key": gkey,
                "pearson_r": rp,
                "p_value": pp,
                "n": nn,
            }
        )
    return n_pairs, out


def _shape_metrics(epochs: List[Dict[str, Any]], perigee_index: int) -> Dict[str, Any]:
    hours = [float(e["hours_from_perigee"]) for e in epochs]
    enh = [float(e["cmb_disformal_enhancement"]) for e in epochs]
    sc_cos = [float(e["sc_cmb_cos_theta"]) for e in epochs]
    cmb_f = [float(e["cmb_modulation_factor"]) for e in epochs]
    rng = [float(e["range_km"]) for e in epochs]
    alt = [float(e["altitude_km"]) for e in epochs]
    ba = [float(e["both_aligned_flag"]) for e in epochs]

    idx_max_enh = int(np.argmax(enh))
    idx_max_abs_cos = int(np.argmax(np.abs(sc_cos)))
    idx_min_range = int(np.argmin(rng))

    inv_range = [1.0 / r if r > 1.0 else float("nan") for r in rng]
    r_inv, p_inv, n_inv = _pearson_along_arc(enh, inv_range)
    r_h, p_h, n_h = _pearson_along_arc(sc_cos, hours)
    r_alt, p_alt, n_alt = _pearson_along_arc(sc_cos, alt)

    flips = 0
    for i in range(1, len(ba)):
        if ba[i] != ba[i - 1]:
            flips += 1

    max_enh = max(enh) if enh else float("nan")
    enh_pg = float(enh[perigee_index]) if perigee_index < len(enh) else float("nan")
    frac_pg = float(enh_pg / max_enh) if max_enh > 0 and math.isfinite(max_enh) else float("nan")

    return {
        "perigee_index": perigee_index,
        "hours_at_perigee": float(hours[perigee_index]) if hours else float("nan"),
        "index_max_cmb_disformal_enhancement": idx_max_enh,
        "hours_peak_cmb_disformal_minus_perigee": float(hours[idx_max_enh]),
        "index_max_abs_sc_cmb_cos_theta": idx_max_abs_cos,
        "hours_peak_abs_sc_cmb_cos_minus_perigee": float(hours[idx_max_abs_cos]),
        "index_min_geocentric_range_km": idx_min_range,
        "hours_min_range_minus_perigee": float(hours[idx_min_range]),
        "cmb_disformal_enhancement_at_perigee": enh_pg,
        "max_cmb_disformal_enhancement": float(max_enh),
        "fraction_disformal_at_perigee_vs_max": frac_pg,
        "mean_both_aligned_flag": float(np.mean(ba)) if ba else float("nan"),
        "count_both_aligned_flag_transitions": flips,
        "temporal_pearson_disformal_vs_inv_range_km": r_inv,
        "temporal_p_value_disformal_vs_inv_range_km": p_inv,
        "temporal_n_disformal_vs_inv_range_km": n_inv,
        "temporal_pearson_sc_cmb_cos_vs_hours_from_perigee": r_h,
        "temporal_p_value_sc_cmb_cos_vs_hours_from_perigee": p_h,
        "temporal_n_sc_cmb_cos_vs_hours": n_h,
        "temporal_pearson_sc_cmb_cos_vs_altitude_km": r_alt,
        "temporal_p_value_sc_cmb_cos_vs_altitude_km": p_alt,
        "temporal_n_sc_cmb_cos_vs_altitude_km": n_alt,
    }


def _optional_residual_geometry_correlations(
    mission: str,
    epochs: List[Dict[str, Any]],
    logger: StepLogger,
) -> Dict[str, Any]:
    """
    If data/time_resolved_flyby_residuals/<mission>.json exists and validates,
    align residual samples to Horizons arc epochs and correlate with geometry.

    Raises RuntimeError on schema / mission mismatch (no silent skip for corrupt files).
    """
    from scripts.utils.flyby_time_series_residuals import (
        residuals_dir,
        load_validated_residual_file,
        build_residual_timeline,
        align_geometry_epochs_to_residuals,
    )

    path = residuals_dir(PROJECT_ROOT) / f"{mission}.json"
    if not path.is_file():
        return {
            "status": "no_residual_file",
            "expected_path": str(path.relative_to(PROJECT_ROOT)),
            "schema": "data/time_resolved_flyby_residuals/schema.json",
        }

    data = load_validated_residual_file(path)
    if str(data.get("mission")) != str(mission):
        raise RuntimeError(
            f"{path}: field mission={data.get('mission')!r} must match filename key {mission!r}"
        )

    max_dt = float(os.environ.get("TEP_042_RESIDUAL_MAX_DT_S", "7200"))
    timeline = build_residual_timeline(data)
    utc_geo = [str(e["utc_iso"]) for e in epochs]
    idxs, y_res, delta_ts, skipped = align_geometry_epochs_to_residuals(
        utc_geo, timeline, max_dt
    )

    n_match = len(y_res)
    if n_match < 5:
        return {
            "status": "insufficient_time_matches",
            "n_geometry_epochs": len(epochs),
            "n_residual_points": len(timeline),
            "n_matched": n_match,
            "n_epochs_skipped_no_neighbor": skipped,
            "max_delta_t_s": max_dt,
            "quantity": data.get("quantity"),
            "reference": data.get("reference"),
        }

    keys_xy = [
        ("cmb_disformal_enhancement", "residual_vs_cmb_disformal_enhancement"),
        ("sc_cmb_cos_theta", "residual_vs_sc_cmb_cos_theta"),
        ("cmb_modulation_factor", "residual_vs_cmb_modulation_factor"),
        ("both_aligned_flag", "residual_vs_both_aligned_flag"),
    ]
    correlations: List[Dict[str, Any]] = []
    for gkey, label in keys_xy:
        xs: List[float] = []
        ys: List[float] = []
        for k in range(n_match):
            j = idxs[k]
            v = epochs[j].get(gkey)
            yv = y_res[k]
            if v is not None and math.isfinite(float(v)) and math.isfinite(float(yv)):
                xs.append(float(v))
                ys.append(float(yv))
        rp, pp, nn = _pearson_along_arc(xs, ys)
        correlations.append(
            {
                "test": label,
                "geometry_key": gkey,
                "pearson_r": rp,
                "p_value": pp,
                "n": nn,
            }
        )

    mean_dt = float(np.mean(np.abs(delta_ts))) if delta_ts else float("nan")

    n_seq, deriv_corr = _residual_sequential_derivative_correlations(idxs, y_res, epochs)

    logger.info(f"  Residual sidecar: matched {n_match} epochs (max |Δt| ≤ {max_dt:.0f} s)")
    if n_seq >= 5:
        logger.info(f"  Sequential d/dt pairs (consecutive Horizons minutes): {n_seq}")
    elif n_seq > 0:
        logger.info(
            f"  Sequential d/dt pairs: {n_seq} (< 5; derivative Pearson tests skipped)"
        )

    return {
        "status": "computed",
        "residual_file": str(path.relative_to(PROJECT_ROOT)),
        "quantity": data.get("quantity"),
        "quantity_description": data.get("quantity_description"),
        "reference": data.get("reference"),
        "max_delta_t_s": max_dt,
        "n_matched": n_match,
        "mean_abs_delta_t_geometry_minus_residual_s": mean_dt,
        "n_epochs_skipped_no_neighbor": skipped,
        "n_residual_points": len(timeline),
        "correlations": correlations,
        "n_sequential_ddt_pairs": n_seq,
        "sequential_derivative_correlations": deriv_corr,
    }


def _load_catalog_published_anomalies(catalog_path: Path) -> Dict[str, Optional[float]]:
    if not catalog_path.is_file():
        raise RuntimeError(
            f"Missing {catalog_path.name}; run Step 003 before Step 042 for pooled epoch export."
        )
    with open(catalog_path, encoding="utf-8") as f:
        c = json.load(f)
    out: Dict[str, Optional[float]] = {}
    for fb in c.get("flybys", []):
        if not fb.get("usable_for_analysis", False):
            continue
        name = str(fb.get("mission_name", ""))
        if name:
            out[name] = fb.get("published_anomaly_mm_s")
    return out


def _build_pooled_epoch_rows(
    out_missions: Dict[str, Any],
    published_by_mission: Dict[str, Optional[float]],
) -> List[Dict[str, Any]]:
    keys_geom = (
        "utc_iso",
        "hours_from_perigee",
        "range_km",
        "altitude_km",
        "sc_geocentric_cmb_parallel_km",
        "earth_barycentric_pos_icrs_au_x",
        "earth_barycentric_pos_icrs_au_y",
        "earth_barycentric_pos_icrs_au_z",
        "earth_barycentric_cmb_parallel_au",
        "sc_cmb_cos_theta",
        "cmb_disformal_enhancement",
        "cmb_modulation_factor",
        "earth_orbital_cmb_proj_kms",
        "both_aligned_flag",
        "heliocentric_distance_au",
    )
    pooled: List[Dict[str, Any]] = []
    for mission, block in sorted(out_missions.items()):
        if not isinstance(block, dict):
            continue
        pub = published_by_mission.get(mission)
        for e in block.get("epoch_table", []):
            if not isinstance(e, dict):
                continue
            row: Dict[str, Any] = {"mission": mission}
            for k in keys_geom:
                row[k] = e.get(k)
            row["published_anomaly_mm_s_mission_aggregate"] = pub
            pooled.append(row)
    return pooled


def main() -> int:
    logger = StepLogger("step_042_time_resolved_cosmography", PROJECT_ROOT)
    logger.header("STEP 042: TIME-RESOLVED COSMOGRAPHIC GEOMETRY ALONG FLYBY ARCS")

    traj_path = PROJECT_ROOT / "results" / "step038_trajectory_series.json"
    if not traj_path.exists():
        raise RuntimeError(
            f"Missing {traj_path.name}; run scripts/steps/step_038_extract_3d_vectors.py first."
        )

    with open(traj_path) as f:
        traj_root = json.load(f)

    missions = traj_root.get("missions")
    if not isinstance(missions, dict) or not missions:
        raise RuntimeError(f"{traj_path.name} has no missions; re-run Step 038.")

    g40 = _load_step040_geometry_module()
    logger.info(
        "Using Earth/Sun ephemeris from Astropy "
        f"(TEP_SOLAR_SYSTEM_EPHEMERIS={os.environ.get('TEP_SOLAR_SYSTEM_EPHEMERIS', 'de440')!r}) "
        "and the same CMB/solar modulation functions as Step 040."
    )

    out_missions: Dict[str, Any] = {}

    for mission, traj in sorted(missions.items()):
        if not isinstance(traj, dict):
            continue
        logger.subsection(mission)
        n_pts = int(traj.get("n_points", 0))
        if n_pts < 3:
            logger.warning(f"  Skipping {mission}: n_points={n_pts} (< 3)")
            continue

        logger.info(f"  Arc epochs: {n_pts}")
        epochs = _compute_arc_epochs(traj, g40, logger)

        perigee_index = int(traj.get("perigee_index", 0))
        metrics = _shape_metrics(epochs, perigee_index)

        out_missions[mission] = {
            "n_epochs": len(epochs),
            "perigee_index": perigee_index,
            "epoch_table": epochs,
            "shape_metrics": metrics,
            "arc_temporal_derivative_geometry": _arc_temporal_derivative_geometry(epochs),
            "time_aligned_residual_analysis": _optional_residual_geometry_correlations(
                mission, epochs, logger
            ),
        }
        logger.info(
            f"  Peak CMB disformal enhancement vs perigee: "
            f"{metrics['hours_peak_cmb_disformal_minus_perigee']:+.2f} h"
        )
        logger.info(
            f"  Temporal r(disformal, 1/range): "
            f"{metrics['temporal_pearson_disformal_vs_inv_range_km']}"
        )

    if not out_missions:
        raise RuntimeError("No mission arcs processed; check Step 038 trajectory output.")

    deltas_disformal = [
        float(m["shape_metrics"]["hours_peak_cmb_disformal_minus_perigee"])
        for m in out_missions.values()
        if isinstance(m.get("shape_metrics"), dict)
    ]
    deltas_abs_cos = [
        float(m["shape_metrics"]["hours_peak_abs_sc_cmb_cos_minus_perigee"])
        for m in out_missions.values()
        if isinstance(m.get("shape_metrics"), dict)
    ]

    n_with_computed_residual_corr = sum(
        1
        for m in out_missions.values()
        if isinstance(m.get("time_aligned_residual_analysis"), dict)
        and m["time_aligned_residual_analysis"].get("status") == "computed"
    )

    aggregate = {
        "n_missions": len(out_missions),
        "n_missions_with_computed_residual_correlations": n_with_computed_residual_corr,
        "mean_hours_peak_cmb_disformal_minus_perigee": float(np.mean(deltas_disformal)),
        "std_hours_peak_cmb_disformal_minus_perigee": float(np.std(deltas_disformal, ddof=1))
        if len(deltas_disformal) > 1
        else 0.0,
        "mean_hours_peak_abs_sc_cmb_cos_minus_perigee": float(np.mean(deltas_abs_cos)),
        "std_hours_peak_abs_sc_cmb_cos_minus_perigee": float(np.std(deltas_abs_cos, ddof=1))
        if len(deltas_abs_cos) > 1
        else 0.0,
    }

    output = {
        "metadata": {
            "step": "042",
            "description": (
                "Time-resolved cosmographic proxies along each Horizons flyby arc, "
                "including Earth barycentric ICRS position (AU) and SC geocentric "
                "position projected onto the CMB dipole direction when Step 038 "
                "includes Cartesian arc positions. Optional residual JSON under "
                "data/time_resolved_flyby_residuals/ (see schema.json): nearest-time "
                "alignment within TEP_042_RESIDUAL_MAX_DT_S (default 7200 s) yields "
                "Pearson tests; consecutive-minute d/dt pairs add sequential_derivative_correlations. "
                "Missing residual files do not fail the step."
            ),
            "inputs": {
                "trajectory_series": str(traj_path.relative_to(PROJECT_ROOT)),
                "geometry_module": "scripts/steps/step_040_cosmographic_shear.py",
                "pooled_epochs_default": "results/step042_pooled_geometry_epochs.json (disable with TEP_042_WRITE_POOLED_EPOCHS=0)",
                "archival_catalog": "results/step003_archival_flyby_catalog.json",
            },
            "earth_sun_ephemeris": os.environ.get("TEP_SOLAR_SYSTEM_EPHEMERIS", "builtin").strip() or "builtin",
            "earth_sun_ephemeris_note": "Per-epoch Earth/Sun state via Astropy (same entry points as Step 040).",
            "interpretation_notes": [
                "Per-epoch published flyby Delta-v(t) is not in the archival pipeline; epoch_table is geometry-only unless a residual sidecar supplies measurements.",
                "temporal_pearson_disformal_vs_inv_range_km is strongly kinematic: cmb_disformal_enhancement "
                "scales with (|v_total|/|v_sc|)^2 while range varies along the hyperbola; high |r| is not, by itself, "
                "evidence of cosmographic coupling.",
                "arc_temporal_derivative_geometry compares time derivatives of geometry scalars along the Horizons cadence; interpret alongside kinematic coupling.",
                "earth_barycentric_cmb_parallel_au is r_Earth·n_CMB in ICRS (AU); its time derivative aligns with earth_orbital_cmb_proj_kms up to frame conventions.",
                "sc_geocentric_cmb_parallel_km is r_sc·n_CMB in geocentric equatorial km; its time derivative is not identical to sc_velocity_cmb_proj_kms unless motion is purely radial in the dipole direction.",
                "Residual correlations use whatever quantity is documented in each sidecar file; "
                "pairwise-Doppler proxies (e.g. Step 030 export) are not OD-scale post-fit residuals unless "
                "the sidecar explicitly supplies such a series.",
                "sequential_derivative_correlations requires both varying value_mm_s and enough consecutive matched Horizons minutes; repeating one aggregate at every epoch yields zero residual slope.",
            ],
        },
        "aggregate_cross_mission": aggregate,
        "missions": out_missions,
    }

    catalog_path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
    published_by_mission = _load_catalog_published_anomalies(catalog_path)
    pooled_rows = _build_pooled_epoch_rows(out_missions, published_by_mission)

    write_pooled = os.environ.get("TEP_042_WRITE_POOLED_EPOCHS", "1").strip().lower() not in (
        "0",
        "false",
        "no",
    )
    if write_pooled:
        pooled_path = PROJECT_ROOT / "results" / "step042_pooled_geometry_epochs.json"
        pooled_payload = {
            "metadata": {
                "step": "042",
                "companion_primary": "results/step042_time_resolved_cosmography.json",
                "n_rows": len(pooled_rows),
                "n_missions": len(out_missions),
                "catalog_source": str(catalog_path.relative_to(PROJECT_ROOT)),
                "published_anomaly_mm_s_mission_aggregate_note": (
                    "One catalog scalar per mission repeated on every arc epoch for that mission: "
                    "use for cross-mission geometry vs aggregate-anomaly exploration only; not a "
                    "time-resolved Delta-v along the trajectory. Per-epoch anomaly-like series belong "
                    "in data/time_resolved_flyby_residuals/<Mission>.json (see schema.json)."
                ),
            },
            "epochs": pooled_rows,
        }
        with open(pooled_path, "w", encoding="utf-8") as f:
            json.dump(pooled_payload, f, indent=2)
        logger.info(f"Pooled epoch table: {pooled_path.name} ({len(pooled_rows)} rows)")

    out_file = PROJECT_ROOT / "results" / "step042_time_resolved_cosmography.json"
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2)

    logger.info(f"\nResults saved to: {out_file}")
    logger.log_step_summary(len(out_missions), "SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
