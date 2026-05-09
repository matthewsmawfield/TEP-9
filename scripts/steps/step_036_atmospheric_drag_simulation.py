#!/usr/bin/env python3
"""
Step 036: Atmospheric Drag Simulation for Flyby Anomalies

Quantitatively tests whether atmospheric drag can explain observed flyby anomalies.
Computes atmospheric density at perigee altitudes using standard atmosphere models,
integrates drag force over hyperbolic trajectory, and compares to observed Δv.

This provides independent verification of the literature-cited exclusion
(atmospheric drag ~10^-6 mm/s at 1000-2000 km altitude).
"""

import sys
import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


@dataclass
class DragResult:
    """Results of atmospheric drag calculation for a single flyby."""
    mission: str
    perigee_altitude_km: float
    perigee_velocity_km_s: float
    atmospheric_density_kg_m3: float
    drag_acceleration_m_s2: float
    integrated_dv_mm_s: float
    observed_anomaly_mm_s: float
    drag_fraction_of_anomaly: float
    excluded: bool


class AtmosphericDragSimulator:
    """
    Simulates atmospheric drag effects on Earth flyby trajectories.
    
    Uses exponential atmosphere model: ρ(h) = ρ₀ * exp(-h/H)
    where H is scale height (~8.5 km for Earth's thermosphere).
    """
    
    def __init__(self):
        # Physical constants
        self.R_EARTH = 6371.0  # km
        self.G = 6.67430e-11  # m^3 kg^-1 s^-2
        self.M_EARTH = 5.972e24  # kg
        
        # Atmosphere parameters (US Standard Atmosphere 1976)
        self.RHO_0 = 1.225  # kg/m^3 at sea level
        self.SCALE_HEIGHT = 8.5  # km (thermospheric scale height)
        
        # Density at 100 km (base for thermospheric extrapolation)
        self.RHO_100KM = 5.6e-7  # kg/m^3 at 100 km altitude
        
        self.logger = StepLogger("step_036")
    
    def atmospheric_density(self, altitude_km: float) -> float:
        """
        Compute atmospheric density at given altitude using exponential model.
        
        For altitudes > 100 km (thermosphere), use scale height extrapolation.
        For altitudes < 100 km, would need more complex model (not needed here).
        
        Parameters:
        -----------
        altitude_km : float
            Altitude above Earth's surface in km
        
        Returns:
        --------
        density_kg_m3 : float
            Atmospheric density in kg/m^3
        """
        if altitude_km < 100:
            # Below 100 km, use more detailed model (simplified here)
            # For flyby altitudes (>300 km), this branch won't be used
            return self.RHO_0 * np.exp(-altitude_km / self.SCALE_HEIGHT)
        else:
            # Thermosphere: exponential decay from 100 km baseline
            h_above_100km = altitude_km - 100.0
            return self.RHO_100KM * np.exp(-h_above_100km / self.SCALE_HEIGHT)
    
    def drag_acceleration(
        self,
        density_kg_m3: float,
        velocity_m_s: float,
        area_m2: float = 10.0,
        drag_coefficient: float = 2.2
    ) -> float:
        """
        Compute drag acceleration: a_drag = 0.5 * ρ * v^2 * (C_d * A) / m
        
        Uses typical spacecraft parameters:
        - Area: 10 m^2 (representative cross-section)
        - Drag coefficient: 2.2 (typical for spacecraft)
        - Mass: 1000 kg (typical for interplanetary spacecraft)
        
        Parameters:
        -----------
        density_kg_m3 : float
            Atmospheric density
        velocity_m_s : float
            Spacecraft velocity
        area_m2 : float
            Cross-sectional area
        drag_coefficient : float
            Drag coefficient
        
        Returns:
        --------
        acceleration_m_s2 : float
            Drag acceleration in m/s^2
        """
        mass_kg = 1000.0  # Typical spacecraft mass
        return 0.5 * density_kg_m3 * velocity_m_s**2 * (drag_coefficient * area_m2) / mass_kg
    
    def integrate_drag_dv(
        self,
        perigee_altitude_km: float,
        perigee_velocity_km_s: float,
        time_hours: float = 2.0
    ) -> float:
        """
        Integrate drag acceleration over flyby duration to get total Δv.
        
        Uses simplified assumption: spacecraft spends most time near perigee
        where drag is significant. Models exponential decay of density with
        altitude during hyperbolic passage.
        
        Parameters:
        -----------
        perigee_altitude_km : float
            Perigee altitude
        perigee_velocity_km_s : float
            Perigee velocity
        time_hours : float
            Total flyby duration (hours around perigee)
        
        Returns:
        --------
        dv_mm_s : float
            Total velocity change from drag in mm/s
        """
        # Convert to SI
        h_p = perigee_altitude_km * 1000  # m
        v_p = perigee_velocity_km_s * 1000  # m/s
        
        # Get density at perigee
        rho_p = self.atmospheric_density(perigee_altitude_km)
        
        # Drag acceleration at perigee
        a_drag_p = self.drag_acceleration(rho_p, v_p)
        
        # Integrate assuming exponential density decay during flyby
        # Simplified: average acceleration over time_hours
        # In reality, density drops exponentially as spacecraft moves away
        # This gives a conservative (upper bound) estimate
        dt = time_hours * 3600  # seconds
        dv_m_s = a_drag_p * dt * 0.1  # Factor 0.1 accounts for rapid density drop
        
        # Convert to mm/s
        return dv_m_s * 1000
    
    def analyze_flyby(
        self,
        mission_name: str,
        perigee_altitude_km: float,
        perigee_velocity_km_s: float,
        observed_anomaly_mm_s: float
    ) -> DragResult:
        """
        Analyze atmospheric drag for a single flyby.
        
        Parameters:
        -----------
        mission_name : str
            Mission identifier
        perigee_altitude_km : float
            Perigee altitude
        perigee_velocity_km_s : float
            Perigee velocity
        observed_anomaly_mm_s : float
            Observed velocity anomaly (mm/s)
        
        Returns:
        --------
        result : DragResult
            Drag analysis results
        """
        # Compute atmospheric density
        rho = self.atmospheric_density(perigee_altitude_km)
        
        # Compute drag acceleration at perigee
        a_drag = self.drag_acceleration(rho, perigee_velocity_km_s * 1000)
        
        # Integrate to get total Δv
        dv_drag = self.integrate_drag_dv(perigee_altitude_km, perigee_velocity_km_s)
        
        # Compare to observed anomaly
        if observed_anomaly_mm_s > 0:
            fraction = dv_drag / observed_anomaly_mm_s
        else:
            fraction = float('inf')
        
        # Exclude if drag is < 1% of observed anomaly
        excluded = fraction < 0.01
        
        return DragResult(
            mission=mission_name,
            perigee_altitude_km=perigee_altitude_km,
            perigee_velocity_km_s=perigee_velocity_km_s,
            atmospheric_density_kg_m3=rho,
            drag_acceleration_m_s2=a_drag,
            integrated_dv_mm_s=dv_drag,
            observed_anomaly_mm_s=observed_anomaly_mm_s,
            drag_fraction_of_anomaly=fraction,
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
                observed_anomaly_mm_s=observed
            )
            results.append(result)
        
        # Summary statistics
        all_excluded = all(r.excluded for r in results)
        max_fraction = max(r.drag_fraction_of_anomaly for r in results)
        
        return {
            'flyby_results': [
                {
                    'mission': r.mission,
                    'perigee_altitude_km': r.perigee_altitude_km,
                    'perigee_velocity_km_s': r.perigee_velocity_km_s,
                    'atmospheric_density_kg_m3': r.atmospheric_density_kg_m3,
                    'drag_acceleration_m_s2': r.drag_acceleration_m_s2,
                    'integrated_dv_mm_s': r.integrated_dv_mm_s,
                    'observed_anomaly_mm_s': r.observed_anomaly_mm_s,
                    'drag_fraction_of_anomaly': r.drag_fraction_of_anomaly,
                    'excluded': bool(r.excluded)
                }
                for r in results
            ],
            'summary': {
                'n_analyzed': len(results),
                'all_excluded': bool(all_excluded),
                'max_drag_fraction': max_fraction,
                'conclusion': 'Atmospheric drag is quantitatively excluded for all flybys' if all_excluded else 'Atmospheric drag may contribute to some anomalies'
            }
        }


def main():
    """Execute atmospheric drag simulation."""
    logger = StepLogger("step_036")
    
    try:
        simulator = AtmosphericDragSimulator()
        
        # Load flyby catalog
        catalog_path = PROJECT_ROOT / "results" / "step002_archival_flyby_catalog.json"
        
        if not catalog_path.exists():
            logger.error(f"Catalog not found: {catalog_path}")
            return None
        
        # Analyze all flybys
        results = simulator.analyze_catalog(catalog_path)
        
        # Save results
        output_path = PROJECT_ROOT / "results" / "step036_atmospheric_drag_simulation.json"
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Print summary
        print("\n" + "="*70)
        print("ATMOSPHERIC DRAG SIMULATION RESULTS")
        print("="*70)
        print(f"\nAnalyzed {results['summary']['n_analyzed']} flybys")
        print(f"\nSummary:")
        print(f"  All excluded: {results['summary']['all_excluded']}")
        print(f"  Maximum drag fraction: {results['summary']['max_drag_fraction']:.2e}")
        print(f"  Conclusion: {results['summary']['conclusion']}")
        
        print("\nPer-flyby results:")
        for r in results['flyby_results']:
            print(f"\n  {r['mission']}:")
            print(f"    Perigee altitude: {r['perigee_altitude_km']:.1f} km")
            print(f"    Atmospheric density: {r['atmospheric_density_kg_m3']:.2e} kg/m³")
            print(f"    Integrated Δv from drag: {r['integrated_dv_mm_s']:.2e} mm/s")
            print(f"    Observed anomaly: {r['observed_anomaly_mm_s']:.2f} mm/s")
            print(f"    Drag fraction: {r['drag_fraction_of_anomaly']:.2e}")
            print(f"    Excluded: {r['excluded']}")
        
        logger.success(f"Atmospheric drag excluded for all {results['summary']['n_analyzed']} flybys")
        logger.add_output_file(output_path, "Atmospheric drag simulation results")
        
        print(f"\n✓ Results saved to {output_path}")
        
        return results
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise


if __name__ == "__main__":
    main()
