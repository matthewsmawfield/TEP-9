#!/usr/bin/env python3
"""
Step 033: Fetch Continuous IRI Electron Density Along Trajectories

This step fetches IRI electron density data for continuous points along
the full historical trajectories, replacing the Chapman layer approximation
with real empirical ionospheric data.

For each flyby trajectory, IRI data is fetched using Celestrak F10.7 values
and sub-satellite latitude/longitude reconstructed from cached JPL Horizons
RA/Dec responses. Missing IRI or space-weather inputs raise errors instead
of substituting defaults.

Performance (Apple Silicon / multi-core)
---------------------------------------
PyIRI is invoked once per retained trajectory sample. Dense Horizons arcs
(~4k points per mission) make naive all-point sampling extremely slow.

Environment overrides (defaults chosen for a good speed / coverage tradeoff):

- ``TEP_EFA_IRI_TRAJECTORY_STRIDE`` (default ``6``): retain every Nth trajectory
  point along the arc, plus a dense band around perigee.
  Set to ``1`` to approximate the legacy all-point behavior (very slow).
- ``TEP_EFA_IRI_PERIGEE_HALFWIDTH`` (default ``200``): number of trajectory indices
  on each side of perigee that are always retained (union with strided samples).
- ``TEP_EFA_IRI_MISSION_WORKERS`` (default ``min(12, cpu_count)``): process missions
  in parallel processes. Set ``1`` to force serial execution.

Step 020 interpolates ``iri_ne_cm3`` vs altitude; subsampling along time is
acceptable provided altitude coverage remains representative (stride + perigee
window + existing perigee altitude ladder preserve the profile shape).
"""

from __future__ import annotations

import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.celestrak_space_weather import lookup_space_weather
from scripts.steps.step_038_extract_3d_vectors import parse_raw_response
from scripts.utils.step_logger import StepLogger

# Check if PyIRI is available
IRI_AVAILABLE = False
try:
    from PyIRI.main_library import IRI_density_1day
    import PyIRI
    IRI_AVAILABLE = True
except ImportError:
    pass


def parse_datetime(dt_str: str) -> datetime:
    """Parse datetime string from JPL Horizons format."""
    formats = [
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d %H:%M:%S.%f",
    ]

    for fmt in formats:
        try:
            return datetime.strptime(dt_str, fmt)
        except ValueError:
            continue

    raise ValueError(f"Cannot parse datetime: {dt_str}")


def geodetic_lat_lon_deg(ra_deg: float, dec_deg: float) -> Tuple[float, float]:
    """Sub-satellite geodetic latitude/longitude from JPL Horizons RA/Dec."""
    lat_deg = float(dec_deg)
    lon_deg = float(ra_deg)
    if lon_deg > 180.0:
        lon_deg -= 360.0
    return lat_deg, lon_deg


def horizons_geometry_by_timestamp(mission: str) -> Dict[str, Tuple[float, float]]:
    """Map trajectory timestamps to geodetic coordinates from cached Horizons data."""
    raw_path = PROJECT_ROOT / "data" / "raw" / "jpl_horizons" / mission / f"{mission}_raw_response.txt"
    if not raw_path.exists():
        raise FileNotFoundError(f"Missing JPL Horizons raw response for {mission}: {raw_path}")

    parsed = parse_raw_response(raw_path)
    if not parsed or not parsed.get("timestamps"):
        raise RuntimeError(f"Could not parse JPL Horizons geometry for {mission}")

    geometry: Dict[str, Tuple[float, float]] = {}
    for dt, ra_deg, dec_deg in zip(parsed["timestamps"], parsed["ra_deg"], parsed["dec_deg"]):
        geometry[dt.strftime("%Y-%m-%d %H:%M:%S")] = geodetic_lat_lon_deg(ra_deg, dec_deg)
    return geometry


def fetch_iri_data(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: int,
    altitude_km: float,
    latitude_deg: float,
    longitude_deg: float,
    f107_sf: float,
) -> float:
    """Fetch IRI electron density data for a specific location and time."""
    if not IRI_AVAILABLE:
        raise RuntimeError("PyIRI library not available. Install with: pip install PyIRI")

    ahr = np.array([hour + minute / 60.0])
    alon = np.array([longitude_deg])
    alat = np.array([latitude_deg])
    aalt = np.array([altitude_km])
    coeff_dir = PyIRI.coeff_dir
    _, _, _, _, _, _, edp = IRI_density_1day(
        year, month, day, ahr, alon, alat, aalt, f107_sf, coeff_dir, 0
    )
    ne = edp[0, 0, 0] / 1e6
    if not np.isfinite(ne) or ne < 0:
        raise RuntimeError(
            f"IRI returned invalid electron density at {year:04d}-{month:02d}-{day:02d} "
            f"{hour:02d}:{minute:02d}:{second:02d} UT"
        )
    return float(ne)


def _chapman_ne_cm3(altitude_km: float, f107_actual: float) -> float:
    h_max = 300.0
    n_max = 2e5 * (f107_actual / 150.0)
    if altitude_km <= h_max:
        scale_height = 50.0
        z = (altitude_km - h_max) / scale_height
        return float(n_max * math.exp(0.5 * (1.0 - z - math.exp(-z))))
    return float(n_max * (h_max / altitude_km) ** 4.5)


def _subsample_indices(n_points: int, perigee_index: int, stride: int, halfwidth: int) -> list[int]:
    """Indices along the trajectory arc for which PyIRI will be evaluated."""
    if n_points <= 0:
        return []
    stride = max(1, int(stride))
    halfwidth = max(0, int(halfwidth))
    coarse = set(range(0, n_points, stride))
    lo = max(0, perigee_index - halfwidth)
    hi = min(n_points - 1, perigee_index + halfwidth)
    dense = set(range(lo, hi + 1))
    coarse.update({0, n_points - 1, perigee_index})
    return sorted(coarse | dense)


def _mission_worker(payload: Tuple[str, str, int, int]) -> Tuple[str, dict]:
    """Process one mission (picklable entry point for ``ProcessPoolExecutor``)."""
    mission, traj_relpath, stride, halfwidth = payload
    if not IRI_AVAILABLE:
        raise RuntimeError("PyIRI not available in worker")

    traj_path = PROJECT_ROOT / traj_relpath
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory file not found: {traj_path}")

    with open(traj_path, "r", encoding="utf-8") as f:
        traj_data = json.load(f)

    timestamps = traj_data["timestamp"]
    range_m = np.asarray(traj_data["range_m"], dtype=np.float64)
    altitude_km_all = (range_m - 6371000.0) / 1000.0
    n_points = len(timestamps)
    perigee_index = int(np.argmin(altitude_km_all))
    perigee_altitude_km = float(altitude_km_all[perigee_index])

    geometry_by_timestamp = horizons_geometry_by_timestamp(mission)
    f107_by_date: Dict[str, float] = {}

    def f107_for_date(date_key: str) -> float:
        if date_key not in f107_by_date:
            f107_by_date[date_key] = float(lookup_space_weather(date_key)["f10_7"])
        return f107_by_date[date_key]

    indices = _subsample_indices(n_points, perigee_index, stride, halfwidth)

    iri_profile: list[float] = []
    chapman_profile: list[float] = []
    altitude_profile: list[float] = []
    timestamp_profile: list[str] = []
    latitude_profile: list[float] = []
    longitude_profile: list[float] = []

    for i in indices:
        timestamp_str = timestamps[i]
        altitude_km = float(altitude_km_all[i])
        if timestamp_str not in geometry_by_timestamp:
            raise RuntimeError(f"No JPL Horizons geometry for {mission} at {timestamp_str}")
        lat_deg, lon_deg = geometry_by_timestamp[timestamp_str]
        dt = parse_datetime(timestamp_str)
        date_key = dt.strftime("%Y-%m-%d")
        f107_actual = f107_for_date(date_key)

        iri_ne = fetch_iri_data(
            year=dt.year,
            month=dt.month,
            day=dt.day,
            hour=dt.hour,
            minute=dt.minute,
            second=dt.second,
            altitude_km=altitude_km,
            latitude_deg=lat_deg,
            longitude_deg=lon_deg,
            f107_sf=f107_actual,
        )
        chapman_ne = _chapman_ne_cm3(altitude_km, f107_actual)

        iri_profile.append(iri_ne)
        chapman_profile.append(chapman_ne)
        altitude_profile.append(altitude_km)
        timestamp_profile.append(timestamp_str)
        latitude_profile.append(lat_deg)
        longitude_profile.append(lon_deg)

    perigee_timestamp = timestamps[perigee_index]
    perigee_lat, perigee_lon = geometry_by_timestamp[perigee_timestamp]
    perigee_dt = parse_datetime(perigee_timestamp)
    perigee_sw = lookup_space_weather(perigee_dt.strftime("%Y-%m-%d"))
    f107_perigee = float(perigee_sw["f10_7"])

    min_traj_altitude = float(np.min(altitude_km_all))
    additional_altitudes = np.arange(perigee_altitude_km, min_traj_altitude, 100.0)

    for alt in additional_altitudes:
        if any(abs(float(alt) - existing) < 50.0 for existing in altitude_profile):
            continue
        alt_f = float(alt)
        iri_ne = fetch_iri_data(
            year=perigee_dt.year,
            month=perigee_dt.month,
            day=perigee_dt.day,
            hour=perigee_dt.hour,
            minute=perigee_dt.minute,
            second=perigee_dt.second,
            altitude_km=alt_f,
            latitude_deg=perigee_lat,
            longitude_deg=perigee_lon,
            f107_sf=f107_perigee,
        )
        chapman_ne = _chapman_ne_cm3(alt_f, f107_perigee)
        iri_profile.append(iri_ne)
        chapman_profile.append(chapman_ne)
        altitude_profile.append(alt_f)
        timestamp_profile.append(perigee_timestamp)
        latitude_profile.append(perigee_lat)
        longitude_profile.append(perigee_lon)

    sorted_indices = np.argsort(np.asarray(altitude_profile, dtype=float))
    altitude_profile = [altitude_profile[i] for i in sorted_indices]
    iri_profile = [iri_profile[i] for i in sorted_indices]
    chapman_profile = [chapman_profile[i] for i in sorted_indices]
    timestamp_profile = [timestamp_profile[i] for i in sorted_indices]
    latitude_profile = [latitude_profile[i] for i in sorted_indices]
    longitude_profile = [longitude_profile[i] for i in sorted_indices]

    alt_arr = np.asarray(altitude_profile, dtype=float)
    pidx = int(np.argmin(alt_arr))
    iri_perigee = float(iri_profile[pidx])
    chapman_perigee = float(chapman_profile[pidx])
    alt_perigee = float(altitude_profile[pidx])

    result = {
        "trajectory": {
            "timestamp": timestamp_profile,
            "altitude_km": altitude_profile,
            "latitude_deg": latitude_profile,
            "longitude_deg": longitude_profile,
            "iri_ne_cm3": iri_profile,
            "chapman_ne_cm3": chapman_profile,
        },
        "coordinates": {
            "perigee_timestamp": perigee_timestamp,
            "perigee_altitude_km": perigee_altitude_km,
            "perigee_latitude_deg": perigee_lat,
            "perigee_longitude_deg": perigee_lon,
        },
        "f107_sfu": f107_perigee,
        "f10_7_source": perigee_sw["data_source"],
        "f10_7_field": perigee_sw["f10_7_field"],
        "geometry_source": "JPL_Horizons_RA_Dec",
        "iri_mean_cm3": float(np.mean(iri_profile)),
        "chapman_mean_cm3": float(np.mean(chapman_profile)),
        "iri_perigee_cm3": iri_perigee,
        "chapman_perigee_cm3": chapman_perigee,
        "iri_perigee_altitude_km": alt_perigee,
        "iri_max_cm3": float(np.max(iri_profile)),
        "chapman_max_cm3": float(np.max(chapman_profile)),
        "altitude_range_km": [float(np.min(alt_arr)), float(np.max(alt_arr))],
        "n_points": len(altitude_profile),
        "sampling_metadata": {
            "trajectory_stride": stride,
            "perigee_index_full_arc": perigee_index,
            "perigee_dense_halfwidth_indices": halfwidth,
            "n_horizon_samples_retained": len(indices),
            "n_horizon_samples_total": n_points,
        },
    }
    return mission, result


def main() -> int:
    """Execute IRI trajectory profile fetching."""
    logger = StepLogger("step_033_iri_trajectory_profile", PROJECT_ROOT)
    logger.header("STEP 033: CONTINUOUS IRI TRAJECTORY PROFILES")

    if not IRI_AVAILABLE:
        logger.error("PyIRI library not available. Install with: pip install PyIRI")
        return 1

    primary_flybys = {
        "NEAR_1998": "data/raw/jpl_horizons/NEAR_1998/NEAR_1998_trajectory.json",
        "Galileo_1990": "data/raw/jpl_horizons/Galileo_1990/Galileo_1990_trajectory.json",
        "Galileo_1992": "data/raw/jpl_horizons/Galileo_1992/Galileo_1992_trajectory.json",
        "Cassini_1999": "data/raw/jpl_horizons/Cassini_1999/Cassini_1999_trajectory.json",
        "Rosetta_2005": "data/raw/jpl_horizons/Rosetta_2005/Rosetta_2005_trajectory.json",
        "Rosetta_2007": "data/raw/jpl_horizons/Rosetta_2007/Rosetta_2007_trajectory.json",
        "Rosetta_2009": "data/raw/jpl_horizons/Rosetta_2009/Rosetta_2009_trajectory.json",
        "MESSENGER_2005": "data/raw/jpl_horizons/MESSENGER_2005/MESSENGER_2005_trajectory.json",
        "Juno_2013": "data/raw/jpl_horizons/Juno_2013/Juno_2013_trajectory.json",
        "Stardust_2001": "data/raw/jpl_horizons/Stardust_2001/Stardust_2001_trajectory.json",
        "OSIRIS-REx_2017": "data/raw/jpl_horizons/OSIRIS-REx_2017/OSIRIS-REx_2017_trajectory.json",
        "BepiColombo_2020": "data/raw/jpl_horizons/BepiColombo_2020/BepiColombo_2020_trajectory.json",
    }

    stride = max(1, int(os.environ.get("TEP_EFA_IRI_TRAJECTORY_STRIDE", "6")))
    halfwidth = max(0, int(os.environ.get("TEP_EFA_IRI_PERIGEE_HALFWIDTH", "200")))
    cpu_n = os.cpu_count() or 8
    default_workers = max(1, min(len(primary_flybys), cpu_n))
    workers = max(1, min(len(primary_flybys), int(os.environ.get("TEP_EFA_IRI_MISSION_WORKERS", str(default_workers)))))

    logger.info(
        f"IRI sampling: stride={stride}, perigee_halfwidth={halfwidth}, "
        f"mission_parallel_workers={workers} (cpu_count={cpu_n})"
    )

    payloads = [(m, p, stride, halfwidth) for m, p in primary_flybys.items()]
    results: Dict[str, dict] = {}

    if workers <= 1:
        for pl in payloads:
            mission, payload = _mission_worker(pl)
            results[mission] = payload
            logger.subsection(f"Processed {mission} ({payload['n_points']} altitude samples)")
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_mission_worker, pl): pl[0] for pl in payloads}
            for fut in as_completed(futures):
                expected = futures[fut]
                mission_name, payload = fut.result()
                if mission_name != expected:
                    raise RuntimeError(f"Mission key mismatch: {mission_name} vs {expected}")
                results[mission_name] = payload
                logger.subsection(f"Processed {mission_name} ({payload['n_points']} altitude samples)")

    ordered = {m: results[m] for m in primary_flybys}

    output_file = PROJECT_ROOT / "results" / "step033_iri_trajectory_profiles.json"
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(ordered, f, indent=2)

    logger.success(f"IRI trajectory profiles saved to: {output_file}")
    logger.success("STEP 033 COMPLETE")

    return 0


if __name__ == "__main__":
    sys.exit(main())
