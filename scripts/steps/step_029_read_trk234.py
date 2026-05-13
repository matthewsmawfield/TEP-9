#!/usr/bin/env python3
"""
Read real TRK-2-34 data using PyTrk234 library.
Extract Doppler measurements for TEP analysis.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json

from scripts.utils.dsn_pds_ingest import flyby_tracking_missions, load_perigee_datetime
from scripts.utils.dsn_tracking_discovery import discover_trk234_file, is_trk234_archive
from scripts.utils.step_logger import StepLogger
from scripts.utils.trk218_extract import extract_trk218_measurements
from scripts.utils.trk234_extract import extract_trk234_measurements


def extract_tracking_measurements(filepath: Path):
    if is_trk234_archive(filepath):
        return extract_trk234_measurements(filepath)
    return extract_trk218_measurements(filepath)


def main():
    """Read TRK-2-34 data and extract Doppler measurements."""
    logger = StepLogger("step_029_read_trk234", PROJECT_ROOT)
    logger.header("STEP 029: READ TRK-2-34 DATA FORMAT")

    try:
        mission_order = flyby_tracking_missions(PROJECT_ROOT)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        logger.log_step_summary(0, "FAILED")
        return 1

    trk_file = None
    selected_mission = None
    perigee = None

    for mission in mission_order:
        try:
            mission_perigee = load_perigee_datetime(PROJECT_ROOT, mission)
        except (FileNotFoundError, KeyError, ValueError) as exc:
            logger.warning(f"Skipping {mission}: perigee metadata unavailable ({exc})")
            continue

        candidate = discover_trk234_file(
            PROJECT_ROOT,
            mission=mission,
            perigee=mission_perigee,
            window_hours=48.0,
        )
        if candidate is None:
            continue
        try:
            measurements = extract_tracking_measurements(candidate)
        except Exception as exc:
            logger.warning(f"Skipping {candidate}: extraction failed ({exc})")
            continue

        if measurements:
            trk_file = candidate
            selected_mission = mission
            perigee = mission_perigee
            break

    if trk_file is None:
        output_dir = PROJECT_ROOT / 'results'
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / 'step029_trk234_extracted.json'

        no_data_result = {
            'status': 'NO_DATA_AVAILABLE',
            'message': 'No perigee-matched TRK-2-34 product found in local DSN cache',
            'note': 'Ingest perigee-matched NASA PDS tracking products before re-running Step 029',
            'missions_checked': mission_order,
            'data_points': 0,
        }

        with open(output_file, 'w', encoding='utf-8') as handle:
            json.dump(no_data_result, handle, indent=2)
        logger.add_output_file(output_file, "TRK-2-34 extraction status JSON")
        logger.warning(f"No perigee-matched TRK-2-34 product found - saved status to: {output_file}")
        logger.log_step_summary(0, "PARTIAL")
        return 0

    logger.info(f"Reading TRK-2-34 file: {trk_file}")
    logger.info(f"File exists: {trk_file.exists()}")
    
    if not trk_file.exists():
        # Save status indicating TRK-2-34 data not available (external data dependency)
        output_dir = PROJECT_ROOT / 'results'
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / 'step029_trk234_extracted.json'
        
        no_data_result = {
            'status': 'NO_DATA_AVAILABLE',
            'message': 'TRK-2-34 raw data must be downloaded from NASA PDS',
            'note': 'External DSN data not included in repository',
            'source_file': str(trk_file),
            'data_points': 0
        }
        
        with open(output_file, 'w') as f:
            json.dump(no_data_result, f, indent=2)
        logger.add_output_file(output_file, "TRK-2-34 extraction status JSON")
        
        logger.warning(f"TRK-2-34 file not found - saved status to: {output_file}")
        logger.info("This is external DSN data that must be downloaded separately")
        logger.log_step_summary(0, "PARTIAL")
        return 0  # Not a failure - external data dependency
    
    try:
        measurements = extract_tracking_measurements(trk_file)
        format_label = 'TRK-2-34' if is_trk234_archive(trk_file) else 'TRK-2-18'
        logger.info(f"Extracted {len(measurements)} {format_label} observables")

        output_dir = PROJECT_ROOT / 'results'
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / 'step029_trk234_extracted.json'

        if not measurements:
            no_data_result = {
                'status': 'NO_DATA_AVAILABLE',
                'message': 'TRK-2-34 file parsed but no radiometric records were found',
                'note': 'Requires valid TRK-2-34 tracking records in NASA PDS source file',
                'source_file': str(trk_file),
                'data_points': 0
            }
            with open(output_file, 'w') as f:
                json.dump(no_data_result, f, indent=2)
            logger.add_output_file(output_file, "TRK-2-34 extraction status JSON")
            logger.warning(f"No radiometric records found - saved status to: {output_file}")
            logger.log_step_summary(0, "PARTIAL")
            return 0

        result = {
            'status': 'SUCCESS',
            'data_source': f'NASA_PDS_{format_label}',
            'mission': selected_mission,
            'perigee_utc': perigee.isoformat() if perigee is not None else None,
            'analysis_mode': 'perigee_matched_flyby',
            'source_file': str(trk_file),
            'data_points': len(measurements),
            'derivation': (
                'sfdu_index values are sequential TRK-2-34 SFDU record positions from PyTrk234'
            ),
            'sample_measurements': measurements[:5],
        }
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)
        logger.add_output_file(output_file, "TRK-2-34 extracted Doppler JSON")

        logger.success(f"Saved extracted data to: {output_file}")
        for index, data_point in enumerate(measurements[:5], start=1):
            logger.info(f"  {index}. {data_point}")

        logger.log_step_summary(0, "SUCCESS")
        return 0

    except ImportError:
        logger.error("PyTrk234 library not installed")
        logger.info("Install with: pip install trk234")
        logger.log_step_summary(0, "FAILED")
        return 1
        
    except Exception as e:
        logger.error(f"Error reading TRK-2-34 file: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        logger.log_step_summary(0, "FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
