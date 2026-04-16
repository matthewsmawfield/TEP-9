#!/usr/bin/env python3
"""
Download Rosetta 2005/2007 SPICE Kernels from ESA SPICE Service

This script downloads the BSP trajectory files for Rosetta's Earth flybys.
"""

import urllib.request
import urllib.error
from pathlib import Path
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger

KERNEL_DIR = PROJECT_ROOT / 'data' / 'raw' / 'spice_kernels'

# ESA SPICE Service FTP URL
ESA_SPICE_URL = "https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/spk/"

# Kernel files needed
ROSETTA_KERNELS = {
    'Rosetta_2005': {
        'filename': 'ORER_______________00031.BSP',
        'description': '1st Earth flyby (March 4, 2005)',
        'flyby_date': '2005-03-04',
        'altitude_km': 1956,
        'anomaly_mm_s': 1.82,
    },
    'Rosetta_2007': {
        'filename': 'ORFR_______________00067.BSP',
        'description': '2nd Earth flyby (November 13, 2007)',
        'flyby_date': '2007-11-13',
        'altitude_km': 5322,
        'anomaly_mm_s': 0.02,
    }
}


def download_kernel(filename: str, url_base: str, output_dir: Path, logger: StepLogger) -> bool:
    """
    Download a SPICE kernel file.
    
    Returns True if successful, False otherwise.
    """
    url = f"{url_base}{filename}"
    output_path = output_dir / filename
    
    logger.progress(f"Downloading: {filename}")
    logger.debug(f"  URL: {url}")
    logger.debug(f"  Output: {output_path}")
    
    start_time = time.time()
    
    try:
        urllib.request.urlretrieve(url, output_path)
        file_size = output_path.stat().st_size
        duration = time.time() - start_time
        logger.success(f"Downloaded: {filename} ({file_size:,} bytes in {duration:.2f}s)")
        logger.add_output_file(output_path, "SPICE BSP kernel file")
        return True
    except urllib.error.HTTPError as e:
        logger.error(f"HTTP Error {e.code}: {e.reason}")
        logger.debug(f"  Failed URL: {url}")
        return False
    except Exception as e:
        logger.error(f"Download failed: {e}")
        logger.debug(f"  Exception type: {type(e).__name__}")
        return False


def main():
    """Download Rosetta SPICE kernels."""
    logger = StepLogger("step_001a_download_spice", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("ROSETTA 2005/2007 SPICE KERNEL DOWNLOAD")
    logger.info(f"Source: ESA SPICE Service")
    logger.info(f"URL: {ESA_SPICE_URL}")
    logger.info(f"Output directory: {KERNEL_DIR}")
    
    logger.section("Configuration")
    logger.parameter("ESA_SPICE_URL", ESA_SPICE_URL)
    logger.parameter("KERNEL_DIR", KERNEL_DIR)
    logger.parameter("Number of kernels", len(ROSETTA_KERNELS))
    
    # Create output directory
    KERNEL_DIR.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Created output directory: {KERNEL_DIR}")
    
    # Download kernels
    success_count = 0
    for mission, info in ROSETTA_KERNELS.items():
        logger.subheader(f"Mission: {mission}")
        logger.info(f"Description: {info['description']}")
        logger.info(f"Flyby date: {info['flyby_date']}")
        logger.info(f"Perigee altitude: {info['altitude_km']} km")
        logger.info(f"Observed anomaly: {info['anomaly_mm_s']} mm/s")
        
        logger.section(f"Downloading {mission}")
        logger.parameter("Filename", info['filename'])
        
        if download_kernel(info['filename'], ESA_SPICE_URL, KERNEL_DIR, logger):
            success_count += 1
        else:
            logger.warning(f"Failed to download {mission} kernel")
    
    # Summary
    duration = time.time() - start_time
    logger.section("DOWNLOAD SUMMARY")
    logger.info(f"Successfully downloaded: {success_count}/{len(ROSETTA_KERNELS)} kernels")
    
    if success_count == len(ROSETTA_KERNELS):
        logger.success("All kernels downloaded successfully!")
        logger.info("Next step: Run 'step_001b_spice_to_json.py' to convert BSP to JSON")
        logger.log_step_summary(duration, "SUCCESS")
        return 0
    else:
        logger.error(f"Only {success_count} of {len(ROSETTA_KERNELS)} kernels downloaded")
        logger.info("Troubleshooting:")
        logger.info("  1. Check internet connection")
        logger.info("  2. Verify ESA SPICE Service is accessible")
        logger.info("  3. Try alternative NASA NAIF URL:")
        logger.info("     https://naif.jpl.nasa.gov/pub/naif/ROSETTA/kernels/spk/former_versions/")
        logger.log_step_summary(duration, "PARTIAL")
        return 1


if __name__ == "__main__":
    sys.exit(main())
