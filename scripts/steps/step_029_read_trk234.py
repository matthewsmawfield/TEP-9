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

from scripts.utils.step_logger import StepLogger


def main():
    """Read TRK-2-34 data and extract Doppler measurements."""
    logger = StepLogger("step_029_read_trk234", PROJECT_ROOT)
    logger.header("STEP 029: READ TRK-2-34 DATA FORMAT")
    
    # TRK-2-34 file path
    trk_file = PROJECT_ROOT / 'data' / 'raw' / 'dsn_tracking' / 'MESSENGER_2005' / 'NASA_PDS' / '151161600sc236dss15ddor_234.dat'
    
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
        import trk234
        
        # Read the TRK-2-34 file
        reader = trk234.Reader(str(trk_file))
        
        logger.info(f"Number of SFDUs: {len(reader.sfdu_list)}")
        
        # Extract Doppler data
        doppler_data = []
        
        for i, sfdu in enumerate(reader.sfdu_list):
            # Check if this SFDU has tracking data
            if hasattr(sfdu, 'trk_chdo') and sfdu.trk_chdo is not None:
                # Get tracking data
                trk = sfdu.trk_chdo
                
                # Extract relevant fields
                data_point = {
                    'sfdu_index': i,
                    'sfdutype': sfdu.label.sfdutype if hasattr(sfdu.label, 'sfdutype') else None,
                }
                
                # Add tracking data fields if available
                if hasattr(trk, 'recvtime'):
                    data_point['recvtime'] = str(trk.recvtime) if trk.recvtime else None
                if hasattr(trk, 'doppler'):
                    data_point['doppler'] = float(trk.doppler) if trk.doppler is not None else None
                if hasattr(trk, 'range'):
                    data_point['range'] = float(trk.range) if trk.range is not None else None
                if hasattr(trk, 'rxtonefreq'):
                    data_point['rxtonefreq'] = float(trk.rxtonefreq) if trk.rxtonefreq is not None else None
                
                # Only add if it has useful data
                if data_point.get('doppler') is not None or data_point.get('range') is not None:
                    doppler_data.append(data_point)
        
        logger.info(f"Extracted {len(doppler_data)} data points with Doppler/range data")
        
        # Save to JSON
        output_dir = PROJECT_ROOT / 'results'
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / 'step029_trk234_extracted.json'

        if not doppler_data:
            no_data_result = {
                'status': 'NO_DATA_AVAILABLE',
                'message': 'TRK-2-34 file parsed but no Doppler/range records were found',
                'note': 'Requires valid tracking records in NASA PDS source file',
                'source_file': str(trk_file),
                'data_points': 0
            }
            with open(output_file, 'w') as f:
                json.dump(no_data_result, f, indent=2)
            logger.add_output_file(output_file, "TRK-2-34 extraction status JSON")
            logger.warning(f"No Doppler/range records found - saved status to: {output_file}")
            logger.log_step_summary(0, "PARTIAL")
            return 0

        with open(output_file, 'w') as f:
            json.dump(doppler_data, f, indent=2)
        logger.add_output_file(output_file, "TRK-2-34 extracted Doppler JSON")
        
        logger.success(f"Saved extracted data to: {output_file}")
        
        # Log first few data points
        logger.info("\nFirst 5 data points:")
        for i, dp in enumerate(doppler_data[:5]):
            logger.info(f"  {i+1}. {dp}")
        
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
