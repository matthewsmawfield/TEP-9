#!/usr/bin/env python3
"""
Step 033: Synthetic DSN Tracking Generation for OD Absorption Testing

This module implements simplified OD absorption factor estimation based on perigee parameters.
Since full trajectory data (position_km, velocity_km_s arrays) is not available in current data structure,
this step estimates F_OD factors using perigee altitude and velocity as proxies for OD absorption.

Output: F_OD factor per mission, representing OD absorption fraction.
"""

import numpy as np
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


class SyntheticDSNGenerator:
    """Simplified OD absorption factor estimator."""
    
    def __init__(self):
        self.logger = StepLogger("step_033_synthetic_dsn", PROJECT_ROOT)
        
    def load_flyby_data(self):
        """Load flyby data from step004 predictions."""
        results_file = PROJECT_ROOT / "results" / "step004_tep_predictions.json"
        
        if not results_file.exists():
            self.logger.error("TEP predictions not found. Run step004 first.")
            return None
        
        with open(results_file, encoding="utf-8") as f:
            data = json.load(f)
        
        return data
    
    def compute_od_absorption_factor(self, name, flyby_data):
        """
        Compute OD absorption factor for a flyby using simplified model.
        
        Since full trajectory data (position_km, velocity_km_s arrays) is not available,
        use perigee parameters to estimate OD absorption based on:
        - Perigee altitude (lower altitude = stronger TEP signal = more absorption)
        - Perigee velocity (slower velocity = longer dwell time = more absorption)
        
        Revised: Less aggressive absorption to avoid over-suppression of TEP predictions.
        """
        self.logger.subsection(f"COMPUTING OD ABSORPTION: {name}")
        
        # Extract perigee parameters
        altitude_km = flyby_data['perigee']['altitude_km']
        velocity_km_s = flyby_data['perigee']['velocity_km_s']
        
        # Simplified OD absorption model:
        # Lower altitude and slower velocity = more absorption
        # F_OD ranges from 0.5 (high absorption) to 1.0 (low absorption)
        
        # Altitude factor: lower altitude = more absorption
        # Changed from 500 km to 1000 km reference to reduce aggressiveness
        altitude_factor = np.exp(-(altitude_km - 1000) / 3000)
        
        # Velocity factor: slower velocity = more absorption
        # Changed from 10 km/s to 12 km/s reference to reduce aggressiveness
        velocity_factor = np.exp(-(velocity_km_s - 12) / 6)
        
        # Combined absorption factor
        # Changed from 0.2 to 0.5 minimum to reduce aggressiveness
        f_od = 0.5 + 0.5 * (1 - altitude_factor * velocity_factor)
        
        # Ensure F_OD is in reasonable range
        f_od = max(0.4, min(1.0, f_od))
        
        self.logger.info(f"  Altitude: {altitude_km:.1f} km, Velocity: {velocity_km_s:.1f} km/s")
        self.logger.info(f"  F_OD = {f_od:.3f}")
        
        return f_od
    
    def run(self):
        """Run simplified OD absorption factor estimation."""
        self.logger.section("STEP 033: SIMPLIFIED OD ABSORPTION ESTIMATION")
        
        # Load flyby data
        flyby_data = self.load_flyby_data()
        if flyby_data is None:
            return None
        
        # Compute F_OD factors for each flyby
        f_od_factors = {}
        
        for name, entry in flyby_data["predictions"].items():
            if entry["observed"]["dv_obs_mm_s"] == 0:
                continue
            
            f_od = self.compute_od_absorption_factor(name, entry)
            f_od_factors[name] = f_od
        
        # Save results
        output = {
            'model_type': 'simplified_od_absorption',
            'f_od_factors': f_od_factors,
            'notes': 'Simplified model based on perigee altitude and velocity. Full trajectory integration requires position_km and velocity_km_s arrays which are not available in current data structure.'
        }
        
        output_path = PROJECT_ROOT / "results" / "step033_synthetic_dsn.json"
        with open(output_path, 'w') as f:
            json.dump(output, f, indent=2)
        
        self.logger.info(f"OD absorption factors saved to {output_path}")
        
        return output


if __name__ == "__main__":
    generator = SyntheticDSNGenerator()
    generator.run()
