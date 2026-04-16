#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 003: DSN Data Request Workflow

Generates actual data request packages for NASA/ESA archives.
This creates the necessary documentation to request raw Level-1 DSN tracking data.
"""

import json
from pathlib import Path
from datetime import datetime, timezone
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


def generate_nasa_dsn_request(mission: str, archive_info: dict) -> dict:
    """Generate NASA DSN data request package."""
    
    request = {
        'request_metadata': {
            'request_type': 'DSN Tracking Data',
            'mission': mission,
            'request_date': datetime.now(timezone.utc).isoformat(),
            'requestor': 'TEP Research Collaboration',
            'scientific_justification': 'Independent verification of Earth flyby anomaly',
            'data_level_requested': 'Level 1 (Processed Doppler Observables)',
            'urgency': 'Standard (2-4 weeks processing)'
        },
        'mission_details': {
            'spacecraft_name': mission.split('_')[0],
            'flyby_date': archive_info.get('flyby_date', 'See catalog'),
            'jpl_horizons_id': archive_info.get('jpl_id', 'See catalog'),
            'tracking_stations': archive_info.get('stations', []),
            'frequency_bands': archive_info.get('bands', [])
        },
        'data_specifications': {
            'observation_types': ['2-way Doppler', '3-way Doppler (if available)'],
            'time_coverage': '±48 hours around perigee passage',
            'required_precision': '0.1 mm/s velocity accuracy'
        },
        'archive_destination': {
            'primary_archive': 'NASA NAIF',
            'access_method': 'PDS download',
            'estimated_volume_gb': archive_info.get('size_gb', 'Unknown')
        },
        'submission_details': {
            'url': 'https://naif.jpl.nasa.gov/naif/data_archives.html',
            'contact': 'naif@jpl.nasa.gov',
            'processing_time': '2-4 weeks'
        }
    }
    
    return request


def main():
    """Generate DSN data request packages."""
    logger = StepLogger("step_003_dsn_data_ingestion", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("STEP 003: DSN Data Request Workflow")
    logger.info("Generating actual data request packages for NASA/ESA archives")
    
    logger.section("Configuration")
    missions = {
        'NEAR_1998': {'agency': 'NASA', 'size_gb': 2.5, 'stations': ['DSS-43']},
        'Galileo_1990': {'agency': 'NASA', 'size_gb': 1.8, 'stations': ['DSS-43']},
        'Cassini_1999': {'agency': 'NASA', 'size_gb': 3.2, 'stations': ['DSS-43']},
        'Juno_2013': {'agency': 'NASA', 'size_gb': 5.5, 'stations': ['DSS-25']},
        'Rosetta_2005': {'agency': 'ESA', 'size_gb': 2.1, 'stations': ['NNO']}
    }
    logger.parameter("Number of missions", len(missions))
    logger.parameter("Total data volume", f"{sum(m['size_gb'] for m in missions.values()):.1f} GB")
    
    output_dir = PROJECT_ROOT / 'data' / 'processed' / 'dsn_requests'
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Output directory: {output_dir}")
    
    generated = {}
    
    for mission, info in missions.items():
        logger.subheader(f"Processing: {mission}")
        logger.info(f"Agency: {info['agency']}")
        logger.info(f"Data volume: {info['size_gb']} GB")
        
        request = generate_nasa_dsn_request(mission, info)
        filename = f"step003_{info['agency']}_request_{mission}.json"
        
        request_path = output_dir / filename
        with open(request_path, 'w') as f:
            json.dump(request, f, indent=2)
        
        generated[mission] = filename
        logger.success(f"Generated: {filename}")
        logger.add_output_file(request_path, f"DSN request for {mission}")
    
    duration = time.time() - start_time
    logger.section("SUMMARY")
    logger.info(f"Generated {len(generated)} request packages in: {output_dir}")
    logger.info(f"Total data volume: {sum(m['size_gb'] for m in missions.values()):.1f} GB")
    logger.info("Next step: Submit requests to NASA/ESA archives")
    logger.info("Timeline: 2-3 months for independent verification")
    logger.log_step_summary(duration, "SUCCESS")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
