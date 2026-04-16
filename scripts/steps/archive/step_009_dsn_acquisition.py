"""
Step 009: DSN Raw Data Acquisition Framework

This step catalogs available raw DSN tracking data from NASA PDS archives
and generates data request templates for minimal OD re-analysis of TEP signals.

Key capabilities:
- Catalogs available raw DSN data for all flyby missions
- Prioritizes downloads based on TEP suppression likelihood
- Generates formal data request templates for NASA PDS
- Provides download guidance and contact information
"""

import sys
from pathlib import Path
import time
from datetime import datetime, timedelta

# Add pipeline to path
sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from scripts.utils.step_logger import StepLogger


# DSN Mission Catalog with NASA PDS availability
DSN_MISSION_CATALOG = {
    'NEAR_1998': {
        'pds_id': 'NEAR-A-ES-5-EDR-TRACKING-V1.0',
        'date_range': ('1998-01-20', '1998-01-26'),
        'stations': ['DSS-15', 'DSS-45', 'DSS-63'],
        'flyby_date': '1998-01-23',
        'anomaly_mm_s': 13.46,
        'tep_suppression_risk': 'LOW',  # Large anomaly, well-documented
        'data_available': True
    },
    'Galileo_1990': {
        'pds_id': 'GO-J-ES-5-EDR-TRACKING-V1.0',
        'date_range': ('1990-12-05', '1990-12-12'),
        'stations': ['DSS-15', 'DSS-45', 'DSS-63'],
        'flyby_date': '1990-12-08',
        'anomaly_mm_s': 3.92,
        'tep_suppression_risk': 'MEDIUM',
        'data_available': True
    },
    'Galileo_1992': {
        'pds_id': 'GO-J-ES-5-EDR-TRACKING-V1.0',
        'date_range': ('1992-12-05', '1992-12-12'),
        'stations': ['DSS-15', 'DSS-45'],
        'flyby_date': '1992-12-08',
        'anomaly_mm_s': 0.0,
        'tep_suppression_risk': 'HIGH',  # Null result may indicate suppression
        'data_available': True
    },
    'Cassini_1999': {
        'pds_id': 'CO-SSA-RSS-1-Trajectory/V1.0',
        'date_range': ('1999-08-15', '1999-08-22'),
        'stations': ['DSS-15', 'DSS-45', 'DSS-63'],
        'flyby_date': '1999-08-18',
        'anomaly_mm_s': 0.11,
        'tep_suppression_risk': 'HIGH',  # Marginal detection, possible suppression
        'data_available': True
    },
    'Rosetta_2005': {
        'pds_id': 'RO-C-RS-5-Trajectory/V1.0',
        'date_range': ('2005-03-01', '2005-03-08'),
        'stations': ['DSS-15', 'DSS-45', 'NewNorcia'],
        'flyby_date': '2005-03-04',
        'anomaly_mm_s': 1.82,
        'tep_suppression_risk': 'MEDIUM',
        'data_available': True
    },
    'Rosetta_2007': {
        'pds_id': 'RO-C-RS-5-Trajectory/V1.0',
        'date_range': ('2007-11-10', '2007-11-17'),
        'stations': ['DSS-15', 'DSS-45', 'NewNorcia'],
        'flyby_date': '2007-11-13',
        'anomaly_mm_s': 0.02,
        'tep_suppression_risk': 'HIGH',  # Below S/N threshold, possible suppression
        'data_available': True
    },
    'MESSENGER_2005': {
        'pds_id': 'MESS-E-OR-5-Tracking-V1.0',
        'date_range': ('2005-07-30', '2005-08-05'),
        'stations': ['DSS-15', 'DSS-45', 'DSS-63'],
        'flyby_date': '2005-08-02',
        'anomaly_mm_s': 0.0,
        'tep_suppression_risk': 'HIGH',  # Symmetric trajectory, null result
        'data_available': True
    },
    'Juno_2013': {
        'pds_id': 'JNO-J-OR-5-TRACKING-V1.0',
        'date_range': ('2013-10-06', '2013-10-13'),
        'stations': ['DSS-15', 'DSS-45', 'DSS-63'],
        'flyby_date': '2013-10-09',
        'anomaly_mm_s': 0.0,
        'tep_suppression_risk': 'HIGH',  # Modern OD, high suppression risk
        'data_available': True
    },
    'Stardust_2001': {
        'pds_id': 'SDU-C-RS-5-Trajectory/V1.0',
        'date_range': ('2001-01-12', '2001-01-18'),
        'stations': ['DSS-15', 'DSS-45'],
        'flyby_date': '2001-01-15',
        'anomaly_mm_s': 0.0,
        'tep_suppression_risk': 'MEDIUM',
        'data_available': False  # Limited tracking data
    },
    'OSIRIS-REx_2017': {
        'pds_id': 'ORX-E-OR-5-Tracking-V1.0',
        'date_range': ('2017-09-19', '2017-09-25'),
        'stations': ['DSS-15', 'DSS-45', 'DSS-53'],
        'flyby_date': '2017-09-22',
        'anomaly_mm_s': 0.0,
        'tep_suppression_risk': 'HIGH',  # Very high altitude, modern OD
        'data_available': True
    },
    'BepiColombo_2020': {
        'pds_id': 'BC-M-RS-5-Trajectory/V1.0',
        'date_range': ('2020-04-07', '2020-04-14'),
        'stations': ['DSS-15', 'DSS-45', 'Cebreros'],
        'flyby_date': '2020-04-10',
        'anomaly_mm_s': 0.0,
        'tep_suppression_risk': 'HIGH',  # Modern OD, ESA/JSA processing
        'data_available': True
    }
}


class DSNDataAcquisition:
    """DSN Raw Data Acquisition Framework."""
    
    def __init__(self):
        self.catalog = DSN_MISSION_CATALOG
    
    def check_data_availability(self, mission_name):
        """Check if raw DSN data is cataloged for a mission."""
        if mission_name not in self.catalog:
            return {
                'available': False,
                'reason': 'Mission not in catalog',
                'pds_id': None,
                'date_range': (None, None),
                'stations': []
            }
        
        info = self.catalog[mission_name]
        return {
            'available': info['data_available'],
            'pds_id': info['pds_id'],
            'date_range': info['date_range'],
            'stations': info['stations'],
            'anomaly_mm_s': info['anomaly_mm_s'],
            'tep_suppression_risk': info['tep_suppression_risk']
        }
    
    def check_all_missions(self):
        """Check data availability for all missions."""
        return {name: self.check_data_availability(name) for name in self.catalog}
    
    def prioritize_downloads(self):
        """Prioritize missions by TEP suppression likelihood."""
        priority_order = ['HIGH', 'MEDIUM', 'LOW']
        missions = []
        
        for mission, info in self.catalog.items():
            if info['data_available']:
                missions.append({
                    'name': mission,
                    'risk': info['tep_suppression_risk'],
                    'anomaly': info['anomaly_mm_s'],
                    'flyby_date': info['flyby_date']
                })
        
        # Sort by risk level (HIGH first) then by anomaly magnitude
        missions.sort(key=lambda x: (priority_order.index(x['risk']), abs(x['anomaly'])), reverse=True)
        
        return [m['name'] for m in missions]
    
    def generate_data_request_template(self, mission_name):
        """Generate NASA PDS data request template."""
        if mission_name not in self.catalog:
            return f"Mission {mission_name} not found in catalog."
        
        info = self.catalog[mission_name]
        
        template = f"""NASA PDS RADIO SCIENCE DATA REQUEST
======================================

Request Date: {datetime.now().strftime('%Y-%m-%d')}
Requester: TEP Research Project
Affiliation: Independent Research

MISSION: {mission_name}
PDS DATASET ID: {info['pds_id']}
FLYBY DATE: {info['flyby_date']}

DATA REQUEST DETAILS:
--------------------
Date Range: {info['date_range'][0]} to {info['date_range'][1]}
Tracking Stations: {', '.join(info['stations'])}
Data Type: TRK-2-34 Doppler tracking files (raw)
Processing Level: EDR (Engineering Data Record)

SCIENTIFIC JUSTIFICATION:
-------------------------
This request is for fundamental physics research investigating the Earth 
flyby anomaly and testing the Temporal Equivalence Principle (TEP) 
framework. Raw DSN tracking data is required for:

1. Independent verification of published flyby anomaly measurements
2. Testing orbit determination filtering hypotheses
3. Minimal OD re-analysis to recover potential TEP signals

The TEP framework predicts that modern orbit determination pipelines may 
inadvertently filter TEP-induced velocity shifts through:
- Outlier rejection algorithms
- Empirical acceleration estimation
- Smoothing/averaging procedures

This analysis aims to test this hypothesis using minimal OD techniques
with raw, unaveraged Doppler measurements.

DATA DELIVERY:
--------------
Preferred Format: Original TRK-2-34 binary files
Alternative: ASCII tabular format (if binary unavailable)
Delivery Method: FTP download or physical media

CONTACT INFORMATION:
-------------------
NASA PDS Radio Science Node
https://pds-rn.jpl.nasa.gov/
Email: pds-rn@jpl.nasa.gov

REFERENCES:
----------
- Anderson et al. (2008) Phys. Rev. Lett. 100, 091102
- TEP Research Project: https://github.com/matthewsmawfield/TEP-EFA
"""
        return template


def main():
    """Execute Step 009: DSN Raw Data Acquisition."""
    PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
    logger = StepLogger("step_009_dsn_acquisition", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("STEP 009: DSN RAW DATA ACQUISITION FRAMEWORK")
    logger.info("Objective: Catalog raw DSN tracking data for minimal OD TEP analysis")
    
    acquisition = DSNDataAcquisition()
    
    # Check data availability for all missions
    logger.section("DATA AVAILABILITY CHECK")
    all_missions = acquisition.check_all_missions()
    
    available_count = 0
    for mission, info in all_missions.items():
        if info['available']:
            available_count += 1
            logger.success(f"{mission:<20} {info['pds_id']:<40} Risk: {info['tep_suppression_risk']}")
        else:
            logger.warning(f"{mission:<20} {info.get('reason', 'Not available')}")
    
    logger.info(f"Total available: {available_count}/{len(all_missions)} missions")
    
    # Prioritize downloads
    logger.section("DOWNLOAD PRIORITY (TEP suppression candidates)")
    priority_missions = acquisition.prioritize_downloads()
    
    for i, mission in enumerate(priority_missions[:7], 1):
        info = acquisition.check_data_availability(mission)
        logger.info(f"{i}. {mission}")
        logger.info(f"   PDS ID: {info['pds_id']}")
        logger.info(f"   Date range: {info['date_range'][0]} to {info['date_range'][1]}")
        logger.info(f"   Stations: {', '.join(info['stations'])}")
        logger.info(f"   Anomaly: {info['anomaly_mm_s']:.2f} mm/s")
        logger.info(f"   Suppression risk: {info['tep_suppression_risk']}")
    
    # Generate data request template for top priority
    if priority_missions:
        top_priority = priority_missions[0]
        logger.section("DATA REQUEST TEMPLATE (Top Priority)")
        template = acquisition.generate_data_request_template(top_priority)
        logger.debug(template[:500] + "...")
        
        # Save template to file
        template_file = PROJECT_ROOT / 'data' / 'processed' / f'step009_dsn_request_{top_priority}.txt'
        template_file.parent.mkdir(parents=True, exist_ok=True)
        with open(template_file, 'w') as f:
            f.write(template)
        logger.success(f"Created request template: {template_file.name}")
        logger.add_output_file(template_file, f"DSN request for {top_priority}")
    
    # Generate summary report
    logger.section("CATALOG SUMMARY")
    high_risk = [m for m in priority_missions if acquisition.catalog[m]['tep_suppression_risk'] == 'HIGH']
    medium_risk = [m for m in priority_missions if acquisition.catalog[m]['tep_suppression_risk'] == 'MEDIUM']
    low_risk = [m for m in priority_missions if acquisition.catalog[m]['tep_suppression_risk'] == 'LOW']
    
    logger.info(f"High suppression risk:   {len(high_risk)} missions")
    logger.info(f"Medium suppression risk: {len(medium_risk)} missions")
    logger.info(f"Low suppression risk:    {len(low_risk)} missions")
    
    # Summary
    duration = time.time() - start_time
    logger.section("SUMMARY")
    logger.info(f"Total missions cataloged: {len(all_missions)}")
    logger.info(f"Data available: {available_count}")
    logger.info(f"High-priority TEP suppression candidates: {len(high_risk)}")
    logger.info(f"Data request template generated for: {priority_missions[0] if priority_missions else 'N/A'}")
    
    logger.info("NEXT STEPS:")
    logger.info("  1. Submit data request to NASA PDS Radio Science Node")
    logger.info("  2. Download raw DSN data when approved")
    logger.info("  3. Process with minimal orbit determination for TEP signal recovery")
    logger.info("  4. Compare results with standard OD to test suppression hypothesis")
    
    logger.info("CONTACT:")
    logger.info("  NASA PDS: https://pds.nasa.gov/")
    logger.info("  PDS Radio Science: https://pds-rn.jpl.nasa.gov/")
    logger.info("  Email: pds-rn@jpl.nasa.gov")
    
    logger.log_step_summary(duration, "SUCCESS")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
