#!/usr/bin/env python3
"""
Step 033: Fetch Continuous IRI Electron Density Along Trajectories

This step fetches IRI electron density data for continuous points along
the full historical trajectories, replacing the Chapman layer approximation
with real empirical ionospheric data.

For each flyby trajectory (301 points at 1-minute intervals), IRI data
is fetched using historical F10.7 solar flux values and approximate
perigee coordinates.

Author: TEP-EFA Pipeline
Date: 2026-05-10
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
import sys
from datetime import datetime
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

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


def fetch_iri_data(year, month, day, hour, minute, second, 
                  altitude_km, latitude_deg, longitude_deg, f107_override=None):
    """
    Fetch IRI electron density data for specific location and time.
    
    Returns: electron density in electrons/cm^3
    """
    if not IRI_AVAILABLE:
        return None
    
    try:
        # PyIRI IRI_density_1day requires:
        # - year, month, day: integers
        # - ahr: array of UT hours
        # - alon: array of longitudes
        # - alat: array of latitudes
        # - aalt: array of altitudes
        # - f107: F10.7 solar flux
        # - coeff_dir: coefficient directory
        # - ccir_or_ursi: coefficient model choice
        
        # Set up arrays for single point query
        ahr = np.array([hour + minute/60.0])
        alon = np.array([longitude_deg])
        alat = np.array([latitude_deg])
        aalt = np.array([altitude_km])
        
        # Use provided F10.7 value or default
        f107 = f107_override if f107_override is not None else 100.0
        
        # Get coefficient directory from PyIRI
        coeff_dir = PyIRI.coeff_dir
        
        # Use CCIR coefficients (0 = CCIR, 1 = URSI)
        ccir_or_ursi = 0
        
        # Call IRI_density_1day
        # Returns: F2, F1, E, Es, sun, mag, EDP (7 values)
        F2, F1, E, Es, sun, mag, edp = IRI_density_1day(
            year, month, day, ahr, alon, alat, aalt, f107, coeff_dir, ccir_or_ursi
        )
        
        # Extract electron density from EDP array
        # EDP shape: (altitudes, longitudes, latitudes, parameters)
        # We want electron density at our single point
        ne = edp[0, 0, 0]  # electrons/m^3
        ne = ne / 1e6  # Convert from m^-3 to cm^-3
        
        return ne
    except Exception as e:
        print(f"Error fetching IRI data: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Execute IRI trajectory profile fetching."""
    logger = StepLogger("step_033_iri_trajectory_profile", PROJECT_ROOT)
    logger.header("STEP 033: CONTINUOUS IRI TRAJECTORY PROFILES")
    
    if not IRI_AVAILABLE:
        logger.error("PyIRI library not available. Install with: pip install PyIRI")
        return 1
    
    # Historical F10.7 solar flux data from step_016_space_weather.py
    historical_f107 = {
        'Galileo_1990': 230.4,  # Solar maximum
        'NEAR_1998': 96.9,      # Moderate
        'Cassini_1999': 130.7,   # Moderate
        'Rosetta_2005': 78.9     # Solar minimum
    }
    
    # Approximate perigee coordinates from step_001c_download_ionospheric.py
    flyby_coordinates = {
        'NEAR_1998': {'lat': -16.0, 'lon': 120.0, 'perigee_altitude_km': 567.9, 'perigee_date': '1998-01-23'},
        'Galileo_1990': {'lat': 15.0, 'lon': -80.0, 'perigee_altitude_km': 972.3, 'perigee_date': '1990-12-08'},
        'Cassini_1999': {'lat': -20.0, 'lon': 130.0, 'perigee_altitude_km': 1197.3, 'perigee_date': '1999-08-18'},
        'Rosetta_2005': {'lat': 5.0, 'lon': 20.0, 'perigee_altitude_km': 1968.7, 'perigee_date': '2005-03-04'}
    }
    
    # Primary flybys with their trajectory files
    primary_flybys = {
        'NEAR_1998': 'data/raw/jpl_horizons/NEAR_1998/NEAR_1998_trajectory.json',
        'Galileo_1990': 'data/raw/jpl_horizons/Galileo_1990/Galileo_1990_trajectory.json',
        'Cassini_1999': 'data/raw/jpl_horizons/Cassini_1999/Cassini_1999_trajectory.json',
        'Rosetta_2005': 'data/raw/jpl_horizons/Rosetta_2005/Rosetta_2005_trajectory.json'
    }
    
    results = {}
    
    for mission, traj_file in primary_flybys.items():
        logger.subsection(f"Processing {mission}")
        
        traj_path = PROJECT_ROOT / traj_file
        if not traj_path.exists():
            logger.warning(f"Trajectory file not found: {traj_path}")
            continue
        
        # Load trajectory data
        with open(traj_path, 'r') as f:
            traj_data = json.load(f)
        
        # Get coordinates and F10.7 for this mission
        coords = flyby_coordinates.get(mission, {'lat': 0.0, 'lon': 0.0})
        f107_actual = historical_f107.get(mission, 100.0)
        
        logger.info(f"Trajectory points: {len(traj_data['timestamp'])}")
        logger.info(f"Location: Lat {coords['lat']:.1f}°, Lon {coords['lon']:.1f}°")
        logger.info(f"Historical F10.7: {f107_actual:.1f} sfu")
        
        # Fetch IRI data for each trajectory point
        iri_profile = []
        chapman_profile = []
        altitude_profile = []
        timestamp_profile = []
        
        total_points = len(traj_data['timestamp'])
        for i, (timestamp_str, range_m) in enumerate(zip(traj_data['timestamp'], traj_data['range_m'])):
            # Calculate altitude
            R_EARTH = 6371000  # meters
            altitude_km = (range_m - R_EARTH) / 1000.0
            
            # Parse datetime
            dt = parse_datetime(timestamp_str)
            
            # Fetch IRI data
            iri_ne = fetch_iri_data(
                year=dt.year,
                month=dt.month,
                day=dt.day,
                hour=dt.hour,
                minute=dt.minute,
                second=dt.second,
                altitude_km=altitude_km,
                latitude_deg=coords['lat'],
                longitude_deg=coords['lon'],
                f107_override=f107_actual
            )
            
            # Compute Chapman layer density for comparison
            h_max = 300  # km
            n_max = 2e5 * (f107_actual / 150.0)  # Scale with F10.7
            if altitude_km <= h_max:
                scale_height = 50
                z = (altitude_km - h_max) / scale_height
                chapman_ne = n_max * np.exp(0.5 * (1 - z - np.exp(-z)))
            else:
                alpha_topside = 4.5
                chapman_ne = n_max * (h_max / altitude_km) ** alpha_topside
            
            iri_profile.append(iri_ne if iri_ne is not None else 0.0)
            chapman_profile.append(chapman_ne)
            altitude_profile.append(altitude_km)
            timestamp_profile.append(timestamp_str)
            
            # Progress update every 50 points
            if (i + 1) % 50 == 0:
                logger.info(f"Progress: {i + 1}/{total_points} points processed")
            
            # Small delay to avoid overwhelming the IRI library
            time.sleep(0.01)
        
        # Add additional altitude points from perigee up to trajectory minimum
        # This ensures IRI data coverage for perigee region where plasma matters most
        min_traj_altitude = min(altitude_profile) if altitude_profile else 60000
        perigee_alt = coords.get('perigee_altitude_km', 500)
        perigee_date_str = coords.get('perigee_date', '2000-01-01')
        perigee_dt = parse_datetime(perigee_date_str + ' 12:00:00')
        
        # Generate altitude points from perigee to min trajectory altitude
        # Use 100 km spacing for efficient coverage
        additional_altitudes = np.arange(perigee_alt, min_traj_altitude, 100)
        logger.info(f"Adding {len(additional_altitudes)} additional altitude points from {perigee_alt:.1f} km to {min_traj_altitude:.1f} km")
        
        for alt in additional_altitudes:
            # Skip if already close to existing point
            if any(abs(alt - existing) < 50 for existing in altitude_profile):
                continue
            
            # Fetch IRI data at this altitude using perigee time and location
            iri_ne = fetch_iri_data(
                year=perigee_dt.year,
                month=perigee_dt.month,
                day=perigee_dt.day,
                hour=perigee_dt.hour,
                minute=perigee_dt.minute,
                second=perigee_dt.second,
                altitude_km=alt,
                latitude_deg=coords['lat'],
                longitude_deg=coords['lon'],
                f107_override=f107_actual
            )
            
            # Compute Chapman layer density
            h_max = 300  # km
            n_max = 2e5 * (f107_actual / 150.0)
            if alt <= h_max:
                scale_height = 50
                z = (alt - h_max) / scale_height
                chapman_ne = n_max * np.exp(0.5 * (1 - z - np.exp(-z)))
            else:
                alpha_topside = 4.5
                chapman_ne = n_max * (h_max / alt) ** alpha_topside
            
            # Add to profiles
            iri_profile.append(iri_ne if iri_ne is not None else 0.0)
            chapman_profile.append(chapman_ne)
            altitude_profile.append(alt)
            timestamp_profile.append(perigee_date_str + ' 12:00:00')
            
            # Small delay
            time.sleep(0.005)
        
        # Sort all profiles by altitude for proper interpolation
        sorted_indices = np.argsort(altitude_profile)
        altitude_profile = [altitude_profile[i] for i in sorted_indices]
        iri_profile = [iri_profile[i] for i in sorted_indices]
        chapman_profile = [chapman_profile[i] for i in sorted_indices]
        timestamp_profile = [timestamp_profile[i] for i in sorted_indices]
        
        logger.info(f"Total altitude points after adding perigee coverage: {len(altitude_profile)}")
        logger.info(f"Altitude range: {min(altitude_profile):.1f} km to {max(altitude_profile):.1f} km")
        
        results[mission] = {
            'trajectory': {
                'timestamp': timestamp_profile,
                'altitude_km': altitude_profile,
                'iri_ne_cm3': iri_profile,
                'chapman_ne_cm3': chapman_profile
            },
            'coordinates': coords,
            'f107_sfu': f107_actual,
            'iri_mean_cm3': float(np.mean(iri_profile)),
            'chapman_mean_cm3': float(np.mean(chapman_profile)),
            'iri_max_cm3': float(np.max(iri_profile)),
            'chapman_max_cm3': float(np.max(chapman_profile)),
            'altitude_range_km': [float(min(altitude_profile)), float(max(altitude_profile))],
            'n_points': len(altitude_profile)
        }
        
        logger.info(f"IRI mean density: {results[mission]['iri_mean_cm3']:.2e} cm^-3")
        logger.info(f"Chapman mean density: {results[mission]['chapman_mean_cm3']:.2e} cm^-3")
    
    # Save results
    output_file = PROJECT_ROOT / 'results' / 'step033_iri_trajectory_profiles.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.success(f"IRI trajectory profiles saved to: {output_file}")
    logger.success("STEP 027 COMPLETE")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
