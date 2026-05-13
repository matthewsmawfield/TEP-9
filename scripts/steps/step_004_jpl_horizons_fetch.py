#!/usr/bin/env python3
"""
Step 004: Fetch Real Trajectory Data from JPL Horizons

This step fetches real spacecraft trajectory data from JPL's official
Horizons system for Earth flyby analysis. This is REAL DATA, not simulated.

JPL Horizons provides:
- Range and range rate data
- Official spacecraft ephemerides
- High-precision trajectory information

Usage:
    python step_004_jpl_horizons_fetch.py --mission NEAR_1998
    python step_004_jpl_horizons_fetch.py --all
"""

import sys
import json
import argparse
import time
from pathlib import Path

# Add parent directory to path for imports
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.jpl_horizons.jpl_horizons_query import JPLHorizonsQuery
from scripts.utils.jpl_horizons.jpl_horizons_processor import JPLHorizonsProcessor
from scripts.utils.step_logger import StepLogger


class JPLHorizonsFetcher:
    """Fetch real trajectory data from JPL Horizons for flyby analysis"""
    
    def __init__(self):
        self.query = JPLHorizonsQuery()
        self.processor = JPLHorizonsProcessor()
        self.base_dir = Path(__file__).parent.parent.parent / 'data' / 'raw' / 'jpl_horizons'
        self.query_dir = Path(__file__).parent.parent / 'utils' / 'jpl_horizons' / 'queries'
        
        # Mission names and their JPL Horizons command names
        self.missions = {
            'NEAR_1998': 'NEAR_1998',
            'Galileo_1990': 'Galileo_1990',
            'Cassini_1999': 'Cassini_1999',
            'Rosetta_2005': 'Rosetta_2005',
            'MESSENGER_2005': 'MESSENGER_2005',
            'Juno_2013': 'Juno_2013',
            'Rosetta_2007': 'Rosetta_2007',
            'Rosetta_2009': 'Rosetta_2009',
            'Galileo_1992': 'Galileo_1992',
            'Stardust_2001': 'Stardust_2001',
            'OSIRIS-REx_2017': 'OSIRIS-REx_2017',
            'BepiColombo_2020': 'BepiColombo_2020',
            'BepiColombo_2021': 'BepiColombo_2021'
        }
    
    def fetch_mission(self, mission_name: str, logger: StepLogger) -> dict:
        """Fetch real trajectory data for a specific mission"""
        logger.progress(f"Fetching JPL Horizons data for {mission_name}")
        
        # Check if query file exists
        query_file = self.query_dir / f"{mission_name}.q"
        if not query_file.exists():
            logger.warning(f"Query file not found: {query_file}")
            return {'mission': mission_name, 'status': 'NO_QUERY_FILE', 'data': None}
        
        logger.debug(f"Using query file: {query_file}")
        
        # Create output directory
        output_dir = self.base_dir / mission_name
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Output directory: {output_dir}")
        
        try:
            # Fetch data from JPL Horizons
            logger.progress("Querying JPL Horizons...")
            response = self.query.fetch(str(query_file))
            logger.debug(f"Response length: {len(response)} characters")
            
            # Save raw response
            raw_file = output_dir / f"{mission_name}_raw_response.txt"
            with open(raw_file, 'w', encoding='utf-8') as f:
                f.write(response)
            logger.success(f"Saved raw response: {raw_file}")
            logger.add_output_file(raw_file, "JPL Horizons raw response")
            
            # Process to JSON
            logger.progress("Processing trajectory data...")
            json_file = output_dir / f"{mission_name}_trajectory.json"
            result = self.processor.process_to_json(response, str(json_file))
            logger.success(f"Saved processed data: {json_file}")
            logger.info(f"Data points: {len(result['timestamp'])}")
            logger.add_output_file(json_file, "Processed trajectory JSON")
            
            return {
                'mission': mission_name,
                'status': 'SUCCESS',
                'data_points': len(result['timestamp']),
                'raw_file': str(raw_file),
                'json_file': str(json_file)
            }
            
        except Exception as e:
            logger.error(f"Error: {e}")
            logger.debug(f"Exception type: {type(e).__name__}")
            return {'mission': mission_name, 'status': 'ERROR', 'error': str(e)}
    
    def fetch_all(self, logger: StepLogger):
        """Fetch real trajectory data for all missions"""
        logger.header("STEP 004: Fetch Real Trajectory Data from JPL Horizons")
        logger.header("STEP 004: FETCH ALL MISSIONS FROM JPL HORIZONS")
        logger.parameter("base_dir", self.base_dir)
        logger.parameter("query_dir", self.query_dir)
        logger.parameter("Number of missions", len(self.missions))
        
        results = []
        for mission_name in self.missions.keys():
            logger.subheader(f"Processing: {mission_name}")
            result = self.fetch_mission(mission_name, logger)
            results.append(result)
        
        # Summary
        logger.section("SUMMARY: JPL Horizons Data Fetch")
        success_count = 0
        for result in results:
            status = result['status']
            if status == 'SUCCESS':
                logger.success(f"{result['mission']:<20} {result['data_points']} data points")
                success_count += 1
            elif status == 'NO_QUERY_FILE':
                logger.warning(f"{result['mission']:<20} {status}")
            else:
                logger.error(f"{result['mission']:<20} {status}")
        
        # Save summary
        summary_file = self.base_dir / 'step004_fetch_summary.json'
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2)
        logger.success(f"Saved summary to: {summary_file}")
        logger.add_output_file(summary_file, "Fetch summary JSON")
        
        logger.info(f"Successfully fetched: {success_count}/{len(self.missions)} missions")
        return results


def main():
    parser = argparse.ArgumentParser(description='Fetch real trajectory data from JPL Horizons')
    parser.add_argument('--mission', type=str, help='Mission name (e.g., NEAR_1998)')
    parser.add_argument('--all', action='store_true', help='Fetch all missions')
    
    args = parser.parse_args()
    
    PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
    logger = StepLogger("step_004_jpl_horizons_fetch", PROJECT_ROOT)
    start_time = time.time()
    
    fetcher = JPLHorizonsFetcher()
    
    # Default to --all if no arguments provided (for pipeline execution)
    if args.mission:
        logger.header(f"STEP 010: Fetch {args.mission} from JPL Horizons")
        result = fetcher.fetch_mission(args.mission, logger)
        duration = time.time() - start_time
        logger.info(f"Result: {result}")
        logger.log_step_summary(duration, "SUCCESS" if result['status'] == 'SUCCESS' else "PARTIAL")
    else:
        # Default to fetching all missions for pipeline
        logger.header("STEP 010: FETCH ALL MISSIONS FROM JPL HORIZONS")
        results = fetcher.fetch_all(logger)
        duration = time.time() - start_time
        all_success = all(r.get('status') == 'SUCCESS' for r in results)
        logger.log_step_summary(duration, "SUCCESS" if all_success else "PARTIAL")


if __name__ == '__main__':
    main()
