#!/usr/bin/env python3
"""
Convert Rosetta SPICE Kernels to JSON Trajectory Format

This script reads the BSP files and extracts state vectors around the flyby epochs,
converting them to the same JSON format used by the pipeline.
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import sys
import urllib.request
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger

# Check for spiceypy availability
try:
    import spiceypy as spice
    HAS_SPICE = True
except ImportError:
    HAS_SPICE = False

KERNEL_DIR = PROJECT_ROOT / 'data' / 'raw' / 'spice_kernels'
OUTPUT_DIR = PROJECT_ROOT / 'data' / 'raw' / 'flyby_trajectories'

# Generic SPICE kernels needed
LEAPSECONDS_KERNEL = 'naif0012.tls'
NAIF_GENERIC_URL = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/'

# Rosetta flyby information
ROSETTA_FLYBYS = {
    'Rosetta_2005': {
        'kernel': 'ORER_______________00031.BSP',
        'flyby_date': '2005-03-04',
        'perigee_time': '2005-03-04T22:10:00',  # Actual perigee from SPICE
        'jpl_id': '-85',
        'observed_anomaly_mm_s': 1.82,
        'anomaly_uncertainty_mm_s': 0.05,
        'days_window': 1,  # ±1 day around perigee
    },
    'Rosetta_2007': {
        'kernel': 'ORFR_______________00067.BSP',
        'flyby_date': '2007-11-13',
        'perigee_time': None,  # Use noon as default
        'jpl_id': '-85',
        'observed_anomaly_mm_s': 0.02,
        'anomaly_uncertainty_mm_s': 0.05,
        'days_window': 2,
    }
}


def download_leapseconds_kernel(logger: StepLogger) -> bool:
    """Download the generic leapseconds kernel if not present."""
    kernel_path = KERNEL_DIR / LEAPSECONDS_KERNEL
    
    if kernel_path.exists():
        logger.success(f"Leapseconds kernel already present: {kernel_path}")
        return True
    
    url = f"{NAIF_GENERIC_URL}{LEAPSECONDS_KERNEL}"
    logger.progress("Downloading leapseconds kernel...")
    logger.debug(f"URL: {url}")
    
    try:
        urllib.request.urlretrieve(url, kernel_path)
        file_size = kernel_path.stat().st_size
        logger.success(f"Downloaded: {file_size:,} bytes")
        logger.add_output_file(kernel_path, "SPICE leapseconds kernel")
        return True
    except Exception as e:
        logger.error(f"Download failed: {e}")
        return False


def load_kernel(kernel_path: Path) -> bool:
    """Load a SPICE kernel file."""
    if not HAS_SPICE:
        return False
    
    try:
        spice.furnsh(str(kernel_path))
        return True
    except Exception as e:
        print(f"  ✗ Error loading kernel: {e}")
        return False


def extract_trajectory(flyby_date: str, perigee_time: str, days_window: int, 
                       target: str = 'ROSETTA', observer: str = 'EARTH') -> dict:
    """
    Extract spacecraft trajectory from SPICE kernel.
    
    Parameters
    ----------
    flyby_date : str
        Date string (YYYY-MM-DD) for naming
    perigee_time : str or None
        Specific perigee time (YYYY-MM-DDTHH:MM:SS) or None to use noon
    days_window : int
        Days before/after perigee to extract
    
    Returns trajectory data in the format matching existing pipeline JSON.
    """
    if not HAS_SPICE:
        return None
    
    try:
        # Use specified perigee time or default to noon on flyby date
        if perigee_time:
            perigee_dt = datetime.strptime(perigee_time, '%Y-%m-%dT%H:%M:%S')
        else:
            perigee_dt = datetime.strptime(flyby_date, '%Y-%m-%d').replace(hour=12, minute=0)
        
        start_dt = perigee_dt - timedelta(days=days_window)
        end_dt = perigee_dt + timedelta(days=days_window)
        
        # Generate time points (30-minute intervals like existing data)
        time_points = []
        current = start_dt
        while current <= end_dt:
            time_points.append(current)
            current += timedelta(minutes=30)
        
        # Extract state vectors
        ephemeris = []
        ranges = []
        
        for t in time_points:
            # Convert to ephemeris time
            et = spice.str2et(t.strftime('%Y-%m-%dT%H:%M:%S'))
            
            # Get state vector (position and velocity)
            state, lt = spice.spkezr(target, et, 'J2000', 'NONE', observer)
            
            # state is a 6-element array: [x, y, z, vx, vy, vz] in km and km/s
            x, y, z = state[0], state[1], state[2]
            vx, vy, vz = state[3], state[4], state[5]
            
            # Compute range
            r = np.sqrt(x**2 + y**2 + z**2)
            ranges.append(r)
            
            # Format datetime string (matching existing format)
            dt_str = f"A.D. {t.strftime('%Y-%b-%d %H:%M:%S')}.0000"
            
            ephemeris.append({
                'datetime': dt_str,
                'x_km': float(x),
                'y_km': float(y),
                'z_km': float(z),
                'vx_km_s': float(vx),
                'vy_km_s': float(vy),
                'vz_km_s': float(vz),
                'range_km': float(r)
            })
        
        # Find perigee (minimum range)
        min_idx = np.argmin(ranges)
        perigee = ephemeris[min_idx]
        
        # Compute velocity magnitude at perigee
        v_mag = np.sqrt(
            perigee['vx_km_s']**2 + 
            perigee['vy_km_s']**2 + 
            perigee['vz_km_s']**2
        )
        
        trajectory = {
            'spacecraft': None,  # Will be set by caller
            'jpl_id': None,  # Will be set by caller
            'flyby_date': flyby_date,
            'n_points': len(ephemeris),
            'perigee': {
                'datetime': perigee['datetime'],
                'range_km': perigee['range_km'],
                'altitude_km': perigee['range_km'] - 6371.0,  # Subtract Earth radius
                'velocity_km_s': float(v_mag)
            },
            'ephemeris': ephemeris,
            'source': 'ESA SPICE Service / ESOC Flight Dynamics',
            'kernel_file': None,  # Will be set by caller
            'extraction_date': datetime.now().isoformat()
        }
        
        return trajectory
        
    except Exception as e:
        print(f"  ✗ Error extracting trajectory: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Convert SPICE kernels to JSON trajectory files."""
    logger = StepLogger("step_001b_spice_to_json", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("ROSETTA SPICE KERNEL TO JSON CONVERTER")
    
    if not HAS_SPICE:
        logger.error("spiceypy is required but not installed")
        logger.info("Install with: pip install spiceypy")
        logger.log_step_summary(0, "FAILED")
        return 1
    
    logger.info(f"Kernel directory: {KERNEL_DIR}")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    
    logger.section("Configuration")
    logger.parameter("KERNEL_DIR", KERNEL_DIR)
    logger.parameter("OUTPUT_DIR", OUTPUT_DIR)
    logger.parameter("HAS_SPICE", HAS_SPICE)
    logger.parameter("Number of flybys", len(ROSETTA_FLYBYS))
    
    # Create directories
    KERNEL_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Created directories")
    
    # Clear any existing kernels
    spice.kclear()
    logger.debug("Cleared existing SPICE kernels")
    
    # Download and load leapseconds kernel
    logger.section("LOADING GENERIC KERNELS")
    if not download_leapseconds_kernel(logger):
        logger.error("Failed to obtain leapseconds kernel")
        logger.log_step_summary(time.time() - start_time, "FAILED")
        return 1
    
    lsk_path = KERNEL_DIR / LEAPSECONDS_KERNEL
    logger.progress("Loading leapseconds kernel...")
    if not load_kernel(lsk_path):
        logger.error("Failed to load leapseconds kernel")
        logger.log_step_summary(time.time() - start_time, "FAILED")
        return 1
    logger.success(f"Loaded: {LEAPSECONDS_KERNEL}")
    
    success_count = 0
    
    for spacecraft_name, info in ROSETTA_FLYBYS.items():
        logger.subheader(f"Processing: {spacecraft_name}")
        logger.info(f"Kernel: {info['kernel']}")
        logger.info(f"Flyby date: {info['flyby_date']}")
        
        kernel_path = KERNEL_DIR / info['kernel']
        
        # Check if kernel exists
        if not kernel_path.exists():
            logger.error(f"Kernel not found: {kernel_path}")
            logger.info("Run 'step_001a_download_spice.py' first to download kernels")
            continue
        
        # Load kernel
        logger.progress("Loading kernel...")
        if not load_kernel(kernel_path):
            continue
        
        # Extract trajectory
        logger.progress(f"Extracting trajectory (±{info['days_window']} days)...")
        if info.get('perigee_time'):
            logger.info(f"Centered on perigee: {info['perigee_time']}")
        trajectory = extract_trajectory(
            info['flyby_date'],
            info.get('perigee_time'),
            info['days_window']
        )
        
        if trajectory is None:
            spice.unload(str(kernel_path))
            continue
        
        # Add metadata
        trajectory['spacecraft'] = spacecraft_name
        trajectory['jpl_id'] = info['jpl_id']
        trajectory['kernel_file'] = info['kernel']
        
        # Add observed anomaly data
        trajectory['observed_anomaly'] = {
            'dv_mm_s': info['observed_anomaly_mm_s'],
            'uncertainty_mm_s': info['anomaly_uncertainty_mm_s'],
            'reference': info.get('anomaly_reference', 'Literature value'),
            'data_source': 'peer_reviewed_literature'
        }
        trajectory['observed_anomaly_mm_s'] = info['observed_anomaly_mm_s']
        trajectory['anomaly_uncertainty_mm_s'] = info['anomaly_uncertainty_mm_s']
        
        # Save to JSON
        output_file = OUTPUT_DIR / f"{spacecraft_name}_trajectory.json"
        logger.progress(f"Saving to: {output_file}")
        
        with open(output_file, 'w') as f:
            json.dump(trajectory, f, indent=2)
        
        file_size = output_file.stat().st_size
        logger.success(f"Saved: {file_size:,} bytes")
        logger.info(f"Perigee altitude: {trajectory['perigee']['altitude_km']:.1f} km")
        logger.info(f"Perigee velocity: {trajectory['perigee']['velocity_km_s']:.2f} km/s")
        logger.info(f"Data points: {trajectory['n_points']}")
        logger.add_output_file(output_file, f"{spacecraft_name} trajectory JSON")
        
        success_count += 1
        
        # Unload kernel
        spice.unload(str(kernel_path))
    
    # Clear all kernels
    spice.kclear()
    
    # Summary
    duration = time.time() - start_time
    logger.section("CONVERSION SUMMARY")
    logger.info(f"Successfully converted: {success_count}/{len(ROSETTA_FLYBYS)} flybys")
    
    if success_count == len(ROSETTA_FLYBYS):
        logger.success("All SPICE kernels converted to JSON!")
        logger.info("Next steps:")
        logger.info("  1. Update 'step_001_data_ingestion.py' to include Rosetta 2005/2007")
        logger.info("  2. Update 'archival_flyby_catalog.json' processed_ephemeris_available = true")
        logger.info("  3. Re-run the pipeline with expanded dataset")
        logger.log_step_summary(duration, "SUCCESS")
        return 0
    else:
        logger.error(f"Only {success_count} of {len(ROSETTA_FLYBYS)} flybys converted")
        logger.log_step_summary(duration, "PARTIAL")
        return 1


if __name__ == "__main__":
    sys.exit(main())
