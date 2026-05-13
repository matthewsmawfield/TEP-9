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
"""

import json
import math
import numpy as np
from pathlib import Path
from typing import Dict, Tuple
import sys
from datetime import datetime
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.celestrak_space_weather import lookup_space_weather
from scripts.steps.step_040a_extract_3d_vectors import parse_raw_response
from scripts.utils.step_logger import StepLogger

# Check if PyIRI is available
IRI_AVAILABLE = False
try:
    from PyIRI.main_library import IRI_density_1day
    import PyIRI
    IRI_AVAILABLE = True
except ImportError:
    pass


def parse_datetime(dt_str):
    """Parse datetime string from JPL Horizons format."""
    formats = [
        '%Y-%m-%d %H:%M:%S',
        '%Y-%m-%dT%H:%M:%S',
        '%Y-%m-%d %H:%M:%S.%f',
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

    geometry = {}
    for dt, ra_deg, dec_deg in zip(parsed["timestamps"], parsed["ra_deg"], parsed["dec_deg"]):
        geometry[dt.strftime("%Y-%m-%d %H:%M:%S")] = geodetic_lat_lon_deg(ra_deg, dec_deg)
    return geometry


def fetch_iri_data(year, month, day, hour, minute, second,
                  altitude_km, latitude_deg, longitude_deg, f107_sf: float) -> float:
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


def main():
    """Execute IRI trajectory profile fetching."""
    logger = StepLogger("step_033_iri_trajectory_profile", PROJECT_ROOT)
    logger.header("STEP 033: CONTINUOUS IRI TRAJECTORY PROFILES")
    
    if not IRI_AVAILABLE:
        logger.error("PyIRI library not available. Install with: pip install PyIRI")
        return 1

    primary_flybys = {
        'NEAR_1998': 'data/raw/jpl_horizons/NEAR_1998/NEAR_1998_trajectory.json',
        'Galileo_1990': 'data/raw/jpl_horizons/Galileo_1990/Galileo_1990_trajectory.json',
        'Galileo_1992': 'data/raw/jpl_horizons/Galileo_1992/Galileo_1992_trajectory.json',
        'Cassini_1999': 'data/raw/jpl_horizons/Cassini_1999/Cassini_1999_trajectory.json',
        'Rosetta_2005': 'data/raw/jpl_horizons/Rosetta_2005/Rosetta_2005_trajectory.json',
        'Rosetta_2007': 'data/raw/jpl_horizons/Rosetta_2007/Rosetta_2007_trajectory.json',
        'Rosetta_2009': 'data/raw/jpl_horizons/Rosetta_2009/Rosetta_2009_trajectory.json',
        'MESSENGER_2005': 'data/raw/jpl_horizons/MESSENGER_2005/MESSENGER_2005_trajectory.json',
        'Juno_2013': 'data/raw/jpl_horizons/Juno_2013/Juno_2013_trajectory.json',
        'Stardust_2001': 'data/raw/jpl_horizons/Stardust_2001/Stardust_2001_trajectory.json',
    }

    results = {}

    for mission, traj_file in primary_flybys.items():
        logger.subsection(f"Processing {mission}")

        traj_path = PROJECT_ROOT / traj_file
        if not traj_path.exists():
            raise FileNotFoundError(f"Trajectory file not found: {traj_path}")

        with open(traj_path, 'r', encoding='utf-8') as f:
            traj_data = json.load(f)

        geometry_by_timestamp = horizons_geometry_by_timestamp(mission)

        iri_profile = []
        chapman_profile = []
        altitude_profile = []
        timestamp_profile = []
        latitude_profile = []
        longitude_profile = []

        total_points = len(traj_data['timestamp'])
        perigee_index = 0
        perigee_altitude_km = math.inf
        perigee_timestamp = None
        perigee_lat = None
        perigee_lon = None

        for i, (timestamp_str, range_m) in enumerate(zip(traj_data['timestamp'], traj_data['range_m'])):
            altitude_km = (range_m - 6371000.0) / 1000.0
            if altitude_km < perigee_altitude_km:
                perigee_altitude_km = altitude_km
                perigee_index = i

            if timestamp_str not in geometry_by_timestamp:
                raise RuntimeError(f"No JPL Horizons geometry for {mission} at {timestamp_str}")

            lat_deg, lon_deg = geometry_by_timestamp[timestamp_str]
            dt = parse_datetime(timestamp_str)
            f107_actual = lookup_space_weather(dt.strftime("%Y-%m-%d"))["f10_7"]

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

            h_max = 300
            n_max = 2e5 * (f107_actual / 150.0)
            if altitude_km <= h_max:
                scale_height = 50
                z = (altitude_km - h_max) / scale_height
                chapman_ne = n_max * math.exp(0.5 * (1 - z - math.exp(-z)))
            else:
                chapman_ne = n_max * (h_max / altitude_km) ** 4.5

            iri_profile.append(iri_ne)
            chapman_profile.append(chapman_ne)
            altitude_profile.append(altitude_km)
            timestamp_profile.append(timestamp_str)
            latitude_profile.append(lat_deg)
            longitude_profile.append(lon_deg)

            if (i + 1) % 50 == 0:
                logger.info(f"Progress: {i + 1}/{total_points} points processed")
            time.sleep(0.01)

        perigee_timestamp = timestamp_profile[perigee_index]
        perigee_lat = latitude_profile[perigee_index]
        perigee_lon = longitude_profile[perigee_index]
        perigee_dt = parse_datetime(perigee_timestamp)
        perigee_sw = lookup_space_weather(perigee_dt.strftime("%Y-%m-%d"))
        f107_perigee = perigee_sw["f10_7"]

        min_traj_altitude = min(altitude_profile)
        additional_altitudes = np.arange(perigee_altitude_km, min_traj_altitude, 100)
        logger.info(
            f"Adding {len(additional_altitudes)} additional altitude points from "
            f"{perigee_altitude_km:.1f} km to {min_traj_altitude:.1f} km"
        )

        for alt in additional_altitudes:
            if any(abs(alt - existing) < 50 for existing in altitude_profile):
                continue

            iri_ne = fetch_iri_data(
                year=perigee_dt.year,
                month=perigee_dt.month,
                day=perigee_dt.day,
                hour=perigee_dt.hour,
                minute=perigee_dt.minute,
                second=perigee_dt.second,
                altitude_km=float(alt),
                latitude_deg=perigee_lat,
                longitude_deg=perigee_lon,
                f107_sf=f107_perigee,
            )

            n_max = 2e5 * (f107_perigee / 150.0)
            if alt <= 300:
                scale_height = 50
                z = (alt - 300) / scale_height
                chapman_ne = n_max * math.exp(0.5 * (1 - z - math.exp(-z)))
            else:
                chapman_ne = n_max * (300 / alt) ** 4.5

            iri_profile.append(iri_ne)
            chapman_profile.append(chapman_ne)
            altitude_profile.append(float(alt))
            timestamp_profile.append(perigee_timestamp)
            latitude_profile.append(perigee_lat)
            longitude_profile.append(perigee_lon)
            time.sleep(0.005)

        sorted_indices = np.argsort(altitude_profile)
        altitude_profile = [altitude_profile[i] for i in sorted_indices]
        iri_profile = [iri_profile[i] for i in sorted_indices]
        chapman_profile = [chapman_profile[i] for i in sorted_indices]
        timestamp_profile = [timestamp_profile[i] for i in sorted_indices]
        latitude_profile = [latitude_profile[i] for i in sorted_indices]
        longitude_profile = [longitude_profile[i] for i in sorted_indices]

        logger.info(f"Total altitude points after adding perigee coverage: {len(altitude_profile)}")
        logger.info(f"Altitude range: {min(altitude_profile):.1f} km to {max(altitude_profile):.1f} km")

        results[mission] = {
            'trajectory': {
                'timestamp': timestamp_profile,
                'altitude_km': altitude_profile,
                'latitude_deg': latitude_profile,
                'longitude_deg': longitude_profile,
                'iri_ne_cm3': iri_profile,
                'chapman_ne_cm3': chapman_profile,
            },
            'coordinates': {
                'perigee_timestamp': perigee_timestamp,
                'perigee_altitude_km': perigee_altitude_km,
                'perigee_latitude_deg': perigee_lat,
                'perigee_longitude_deg': perigee_lon,
            },
            'f107_sfu': f107_perigee,
            'f10_7_source': perigee_sw["data_source"],
            'f10_7_field': perigee_sw["f10_7_field"],
            'geometry_source': 'JPL_Horizons_RA_Dec',
            'iri_mean_cm3': float(np.mean(iri_profile)),
            'chapman_mean_cm3': float(np.mean(chapman_profile)),
            'iri_max_cm3': float(np.max(iri_profile)),
            'chapman_max_cm3': float(np.max(chapman_profile)),
            'altitude_range_km': [float(min(altitude_profile)), float(max(altitude_profile))],
            'n_points': len(altitude_profile),
        }

        logger.info(f"Perigee location: Lat {perigee_lat:.2f}°, Lon {perigee_lon:.2f}°")
        logger.info(f"Historical F10.7: {f107_perigee:.1f} sfu ({perigee_sw['data_source']})")
        logger.info(f"IRI mean density: {results[mission]['iri_mean_cm3']:.2e} cm^-3")
        logger.info(f"Chapman mean density: {results[mission]['chapman_mean_cm3']:.2e} cm^-3")
    
    # Save results
    output_file = PROJECT_ROOT / 'results' / 'step033_iri_trajectory_profiles.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.success(f"IRI trajectory profiles saved to: {output_file}")
    logger.success("STEP 033 COMPLETE")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
