"""
Step 026: IRI Plasma Validation for Flyby Dates

This step fetches historical International Reference Ionosphere (IRI) data
for the specific flyby dates (1990, 1998, 1999, 2005) to validate the
Chapman layer approximation used in plasma modulation.

This turns the theoretical plasma modulation into an observational fact.
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger

try:
    from PyIRI.main_library import IRI_density_1day
    IRI_AVAILABLE = True
except ImportError as e:
    IRI_AVAILABLE = False
    print(f"PyIRI library not available: {e}")
    print("Install with: pip install pyiri")


def find_perigee_point(trajectory_data):
    """
    Find the perigee point (minimum range) in trajectory data.
    
    Returns: dict with datetime, altitude_km, latitude_deg, longitude_deg
    """
    # JPL Horizons data has range_m, not altitude_km
    # Compute altitude from range (range - Earth radius)
    R_earth = 6371000  # Earth radius in meters
    
    ranges = np.array(trajectory_data['range_m'])
    min_idx = np.argmin(ranges)
    
    altitude_m = ranges[min_idx] - R_earth
    altitude_km = altitude_m / 1000.0
    
    # JPL Horizons doesn't provide lat/lon in the basic ephemeris
    # We'll estimate from the timestamp and use a representative location
    # For Earth flybys, we can use the sub-solar point or a standard location
    # Since the flyby is near Earth, we'll use 0° lat, 0° lon as a reasonable approximation
    # for ionospheric density at the spacecraft altitude
    
    return {
        'datetime': trajectory_data['timestamp'][min_idx],
        'altitude_km': altitude_km,
        'latitude_deg': 0.0,  # Equatorial approximation
        'longitude_deg': 0.0  # Greenwich meridian approximation
    }


def parse_datetime(dt_str):
    """Parse datetime string from trajectory data."""
    formats = [
        '%Y-%m-%d %H:%M:%S',
        '%Y-%b-%d %H:%M:%S.%f',
        'A.D. %Y-%b-%d %H:%M:%S.%f',
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
        import PyIRI
        coeff_dir = PyIRI.coeff_dir
        
        # Use CCIR coefficients (0 = CCIR, 1 = URSI)
        ccir_or_ursi = 0
        
        # Call IRI_density_1day
        # Returns: F2, F1, E, Es, sun, mag, EDP (7 values)
        F2, F1, E, Es, sun, mag, edp = IRI_density_1day(
            year, month, day, ahr, alon, alat, aalt, f107, coeff_dir, ccir_or_ursi
        )
        
        # Extract electron density from EDP (Electron Density Profile)
        # EDP shape is [N_T, N_V, N_G] = [time, altitude, location]
        # We have 1 time, 1 altitude, 1 location, so we extract the single value
        ne = edp[0, 0, 0]  # electrons/m^3
        
        # Convert from m^-3 to cm^-3
        ne = ne / 1e6
        
        return ne
    except Exception as e:
        print(f"Error fetching IRI data: {e}")
        import traceback
        traceback.print_exc()
        return None


def chapman_layer_density(altitude_km, h_max=300, H_scale=50, n_max=1e6):
    """
    Compute Chapman layer electron density.
    
    Parameters:
    - altitude_km: altitude in km
    - h_max: peak altitude (300 km)
    - H_scale: scale height (50 km)
    - n_max: peak density (1e6 cm^-3 for solar maximum)
    
    Returns: electron density in cm^-3
    """
    if altitude_km < 0:
        altitude_km = 0
    
    delta_h = altitude_km - h_max
    exponent = 0.5 * (1 - delta_h / H_scale - np.exp(-delta_h / H_scale))
    n_iono = n_max * np.exp(exponent)
    
    return n_iono


def main():
    """Execute IRI validation for all primary flybys."""
    logger = StepLogger("step_026_iri_validation", PROJECT_ROOT)
    logger.header("STEP 026: IRI PLASMA VALIDATION")
    
    if not IRI_AVAILABLE:
        logger.error("IRI2016 library not available. Install with: pip install iri2016")
        return
    
    # JPL Horizons doesn't provide lat/lon, but we have approximate perigee coordinates
    # from step_001c_download_ionospheric.py for more accurate IRI calculations
    flyby_coordinates = {
        'NEAR_1998': {'lat': -16.0, 'lon': 120.0},
        'Galileo_1990': {'lat': 15.0, 'lon': -80.0},
        'Cassini_1999': {'lat': -20.0, 'lon': 130.0},
        'Rosetta_2005': {'lat': 5.0, 'lon': 20.0}
    }
    
    # Historical F10.7 solar flux data from step_016_space_weather.py
    # (downloaded from Celestrak SW-All.csv)
    historical_f107 = {
        'Galileo_1990': 230.4,  # Solar maximum
        'NEAR_1998': 96.9,      # Moderate
        'Cassini_1999': 130.7,   # Moderate
        'Rosetta_2005': 78.9     # Solar minimum
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
        
        # Find perigee point
        perigee = find_perigee_point(traj_data)
        logger.info(f"Perigee: {perigee['datetime']}, Altitude: {perigee['altitude_km']:.1f} km")
        
        # Get actual coordinates and F10.7 for this mission
        coords = flyby_coordinates.get(mission, {'lat': 0.0, 'lon': 0.0})
        f107_actual = historical_f107.get(mission, 100.0)
        
        logger.info(f"Location: Lat {coords['lat']:.1f}°, Lon {coords['lon']:.1f}°")
        logger.info(f"Historical F10.7: {f107_actual:.1f} sfu")
        
        # Parse datetime
        dt = parse_datetime(perigee['datetime'])
        
        # Fetch IRI data with actual coordinates and F10.7
        logger.info("Fetching IRI electron density...")
        iri_ne = fetch_iri_data(
            year=dt.year,
            month=dt.month,
            day=dt.day,
            hour=dt.hour,
            minute=dt.minute,
            second=dt.second,
            altitude_km=perigee['altitude_km'],
            latitude_deg=coords['lat'],
            longitude_deg=coords['lon'],
            f107_override=f107_actual
        )
        
        # Compute Chapman layer density
        chapman_ne = chapman_layer_density(perigee['altitude_km'])
        
        logger.info(f"IRI electron density: {iri_ne:.2e} cm^-3" if iri_ne else "IRI data unavailable")
        logger.info(f"Chapman layer density: {chapman_ne:.2e} cm^-3")
        
        if iri_ne is not None:
            ratio = iri_ne / chapman_ne if chapman_ne > 0 else 0
            logger.info(f"IRI/Chapman ratio: {ratio:.2f}")
        
        results[mission] = {
            'perigee': perigee,
            'datetime_str': str(dt),
            'coordinates': coords,
            'f107_sfu': f107_actual,
            'iri_ne_cm3': float(iri_ne) if iri_ne is not None else None,
            'chapman_ne_cm3': float(chapman_ne),
            'ratio': float(iri_ne / chapman_ne) if (iri_ne is not None and chapman_ne > 0) else None
        }
    
    # Save results
    output_file = PROJECT_ROOT / 'results' / 'step026_iri_validation.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Results saved to: {output_file}")
    logger.success("STEP 026 COMPLETE")


if __name__ == '__main__':
    main()
