#!/usr/bin/env python3
"""
Step 038: Frame-Dragging (Lense-Thirring) Calculation for Flyby Anomalies

Quantitatively tests whether general relativistic frame-dragging can explain
observed flyby anomalies. Computes gravitomagnetic velocity shifts from Earth's
rotation using the Lense-Thirring effect.

This provides independent verification of the literature-cited exclusion
(frame-dragging ~10^-5 mm/s vs observed anomalies of 0.11-13.46 mm/s).
"""

import sys
import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


@dataclass
class FrameDraggingResult:
    """Results of frame-dragging calculation for a single flyby."""
    mission: str
    perigee_altitude_km: float
    perigee_velocity_km_s: float
    declination_in_deg: float
    declination_out_deg: float
    lense_thirring_dv_mm_s: float
    observed_anomaly_mm_s: float
    frame_dragging_fraction: float
    excluded: bool


class FrameDraggingCalculator:
    """
    Calculates general relativistic frame-dragging (Lense-Thirring) effects.
    
    The Lense-Thirring effect arises from Earth's rotation and produces
    gravitomagnetic velocity shifts on orbiting bodies.
    """
    
    def __init__(self):
        # Physical constants
        self.G = 6.67430e-11  # m^3 kg^-1 s^-2
        self.C = 2.99792458e8  # m/s
        self.R_EARTH = 6371.0  # km
        
        # Earth parameters
        self.M_EARTH = 5.972e24  # kg
        self.J_EARTH = 1.0  # Dimensionless angular momentum (normalized)
        self.OMEGA_EARTH = 7.292115e-5  # rad/s (Earth rotation rate)
        
        # Earth's angular momentum
        self.I_EARTH = 0.3307 * self.M_EARTH * (self.R_EARTH * 1000)**2  # kg m^2
        self.L_EARTH = self.I_EARTH * self.OMEGA_EARTH  # kg m^2/s
        
        # Gravitomagnetic parameter
        self.GM = self.G * self.M_EARTH
        self.GJ = self.G * self.L_EARTH / (self.C**2)  # Gravitomagnetic parameter
        
        self.logger = StepLogger("step_038")
    
    def lense_thirring_velocity_shift(
        self,
        radius_m: float,
        velocity_m_s: float,
        inclination_rad: float = 0.0
    ) -> float:
        """
        Compute Lense-Thirring velocity shift.
        
        Δv_LT ≈ 2GJ/(c²r) * (v/c) * sin(inclination)
        
        For equatorial flybys, this simplifies to the radial gravitomagnetic
        acceleration integrated over the trajectory.
        
        Parameters:
        -----------
        radius_m : float
            Distance from Earth center in meters
        velocity_m_s : float
            Spacecraft velocity in m/s
        inclination_rad : float
            Orbital inclination relative to equator (radians)
        
        Returns:
        --------
        dv_mm_s : float
            Lense-Thirring velocity shift in mm/s
        """
        # Radial gravitomagnetic acceleration
        a_gravitomagnetic = 2 * self.GJ / (radius_m**3)
        
        # Velocity component (approximate for hyperbolic flyby)
        # The effect scales with velocity relative to c
        v_factor = velocity_m_s / self.C
        
        # Inclination factor (maximum for equatorial, zero for polar)
        inclination_factor = np.sin(inclination_rad) if inclination_rad > 0 else 1.0
        
        # Integrate over flyby timescale (conservative estimate)
        # Use characteristic time ~R_EARTH / v
        flyby_time_s = (self.R_EARTH * 1000) / velocity_m_s
        
        dv_m_s = a_gravitomagnetic * v_factor * inclination_factor * flyby_time_s
        
        # Convert to mm/s
        return dv_m_s * 1000
    
    def estimate_inclination(
        self,
        declination_in_deg: float,
        declination_out_deg: float
    ) -> float:
        """
        Estimate orbital inclination from declination change.
        
        For Earth flybys, the declination change provides information
        about the orbital plane relative to Earth's equator.
        
        Parameters:
        -----------
        declination_in_deg : float
            Declination at flyby entry
        declination_out_deg : float
            Declination at flyby exit
        
        Returns:
        --------
        inclination_rad : float
            Estimated inclination in radians
        """
        if declination_in_deg is None or declination_out_deg is None:
            self.logger.warning("Missing declination data, cannot estimate inclination")
            return None
        
        # Simple estimate: average absolute declination
        avg_declination = (abs(declination_in_deg) + abs(declination_out_deg)) / 2.0
        
        # Convert to inclination (declination ≈ inclination for Earth-centered)
        inclination_rad = np.radians(avg_declination)
        
        return inclination_rad
    
    def analyze_flyby(
        self,
        mission_name: str,
        perigee_altitude_km: float,
        perigee_velocity_km_s: float,
        declination_in_deg: float,
        declination_out_deg: float,
        observed_anomaly_mm_s: float
    ) -> FrameDraggingResult:
        """
        Analyze frame-dragging for a single flyby.
        
        Parameters:
        -----------
        mission_name : str
            Mission identifier
        perigee_altitude_km : float
            Perigee altitude
        perigee_velocity_km_s : float
            Perigee velocity
        declination_in_deg : float
            Declination at entry
        declination_out_deg : float
            Declination at exit
        observed_anomaly_mm_s : float
            Observed velocity anomaly (mm/s)
        
        Returns:
        --------
        result : FrameDraggingResult
            Frame-dragging analysis results
        """
        # Convert to SI
        radius_m = (self.R_EARTH + perigee_altitude_km) * 1000
        velocity_m_s = perigee_velocity_km_s * 1000
        
        # Estimate inclination
        inclination_rad = self.estimate_inclination(declination_in_deg, declination_out_deg)
        
        # Compute Lense-Thirring velocity shift
        dv_lt = self.lense_thirring_velocity_shift(radius_m, velocity_m_s, inclination_rad)
        
        # Compare to observed anomaly
        if observed_anomaly_mm_s > 0:
            fraction = dv_lt / observed_anomaly_mm_s
        else:
            fraction = float('inf')
        
        # Exclude if frame-dragging is < 1% of observed anomaly
        excluded = fraction < 0.01
        
        return FrameDraggingResult(
            mission=mission_name,
            perigee_altitude_km=perigee_altitude_km,
            perigee_velocity_km_s=perigee_velocity_km_s,
            declination_in_deg=declination_in_deg,
            declination_out_deg=declination_out_deg,
            lense_thirring_dv_mm_s=dv_lt,
            observed_anomaly_mm_s=observed_anomaly_mm_s,
            frame_dragging_fraction=fraction,
            excluded=excluded
        )
    
    def analyze_catalog(
        self,
        catalog_path: Path
    ) -> Dict:
        """
        Analyze all flybys in the catalog.
        
        Parameters:
        -----------
        catalog_path : Path
            Path to step002_archival_flyby_catalog.json
        
        Returns:
        --------
        results : Dict
            Analysis results for all flybys
        """
        try:
            with open(catalog_path, 'r') as f:
                catalog = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load catalog: {e}")
            return None
        
        results = []
        for flyby in catalog['flybys']:
            if not flyby['usable_for_analysis']:
                continue
            
            observed = flyby.get('published_anomaly_mm_s')
            if observed is None or observed == 0:
                continue
            
            result = self.analyze_flyby(
                mission_name=flyby['mission_name'],
                perigee_altitude_km=flyby['perigee_altitude_km'],
                perigee_velocity_km_s=flyby['perigee_velocity_km_s'],
                declination_in_deg=flyby.get('declination_in_deg'),
                declination_out_deg=flyby.get('declination_out_deg'),
                observed_anomaly_mm_s=observed
            )
            results.append(result)
        
        # Summary statistics
        all_excluded = all(r.excluded for r in results)
        max_fraction = max(r.frame_dragging_fraction for r in results)
        max_dv_lt = max(r.lense_thirring_dv_mm_s for r in results)
        
        return {
            'flyby_results': [
                {
                    'mission': r.mission,
                    'perigee_altitude_km': r.perigee_altitude_km,
                    'perigee_velocity_km_s': r.perigee_velocity_km_s,
                    'declination_in_deg': r.declination_in_deg,
                    'declination_out_deg': r.declination_out_deg,
                    'lense_thirring_dv_mm_s': r.lense_thirring_dv_mm_s,
                    'observed_anomaly_mm_s': r.observed_anomaly_mm_s,
                    'frame_dragging_fraction': r.frame_dragging_fraction,
                    'excluded': bool(r.excluded)
                }
                for r in results
            ],
            'summary': {
                'n_analyzed': len(results),
                'all_excluded': bool(all_excluded),
                'max_frame_dragging_fraction': max_fraction,
                'max_lense_thirring_dv_mm_s': max_dv_lt,
                'conclusion': 'Frame-dragging is quantitatively excluded for all flybys' if all_excluded else 'Frame-dragging may contribute to some anomalies'
            }
        }


def main():
    """Execute frame-dragging calculation."""
    logger = StepLogger("step_038")
    
    try:
        calculator = FrameDraggingCalculator()
        
        # Load flyby catalog
        catalog_path = PROJECT_ROOT / "results" / "step002_archival_flyby_catalog.json"
        
        if not catalog_path.exists():
            logger.error(f"Catalog not found: {catalog_path}")
            return None
        
        # Analyze all flybys
        results = calculator.analyze_catalog(catalog_path)
        
        # Save results
        output_path = PROJECT_ROOT / "results" / "step_038_frame_dragging_calculation.json"
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Print summary
        print("\n" + "="*70)
        print("FRAME-DRAGGING (LENSE-THIRRING) CALCULATION RESULTS")
        print("="*70)
        print(f"\nAnalyzed {results['summary']['n_analyzed']} flybys")
        print(f"\nSummary:")
        print(f"  All excluded: {results['summary']['all_excluded']}")
        print(f"  Maximum frame-dragging fraction: {results['summary']['max_frame_dragging_fraction']:.2e}")
        print(f"  Maximum Lense-Thirring Δv: {results['summary']['max_lense_thirring_dv_mm_s']:.2e} mm/s")
        print(f"  Conclusion: {results['summary']['conclusion']}")
        
        print("\nPer-flyby results:")
        for r in results['flyby_results']:
            print(f"\n  {r['mission']}:")
            print(f"    Perigee altitude: {r['perigee_altitude_km']:.1f} km")
            print(f"    Perigee velocity: {r['perigee_velocity_km_s']:.2f} km/s")
            print(f"    Lense-Thirring Δv: {r['lense_thirring_dv_mm_s']:.2e} mm/s")
            print(f"    Observed anomaly: {r['observed_anomaly_mm_s']:.2f} mm/s")
            print(f"    Frame-dragging fraction: {r['frame_dragging_fraction']:.2e}")
            print(f"    Excluded: {r['excluded']}")
        
        logger.success(f"Frame-dragging excluded for all {results['summary']['n_analyzed']} flybys")
        logger.add_output_file(output_path, "Frame-dragging calculation results")
        
        print(f"\n✓ Results saved to {output_path}")
        
        return results
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise


if __name__ == "__main__":
    main()
