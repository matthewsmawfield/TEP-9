"""
Saturation Model from Pulsar Data (Paper 14)

This module implements the saturation model from globular cluster pulsar analysis.
From Paper 14: TEP effects saturate at high densities rather than scaling linearly.

Key insight from Paper 14:
-- Suppressed density scaling: observed slope = 0.39 vs Newtonian = 0.72 (4.1σ rejection)
-- This indicates TEP effects saturate in dense environments
- Saturation threshold: ρ_T ≈ 20 g/cm³ (universal Temporal Topology density)

For flyby analysis, this means:
- Low-altitude flybys (higher density) may expe
- High-altitude flybys (lower density) are in unscreened regime
- The saturation behavior can explain heterogeneity in β fits
"""

import numpy as np
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger

# Saturation parameters from Paper 14
RHO_T = 20.0  # g/cm³ - universal Temporal Topology density
SATURATION_SLOPE_OBSERVED = 0.39  # Observed slope from pulsar data
SATURATION_SLOPE_NEWTONIAN = 0.72  # Expected Newtonian slope
SATURATION_SIGMA_REJECTION = 4.1  # Sigma level of rejection


class SaturationModel:
    """
    Implements saturation behavior from Paper 14 pulsar analysis.
    
    The model captures the observation that TEP effects saturate
    in dense environments rather than scaling indefinitely.
    """
    
    def __init__(self):
        self.logger = StepLogger("step_012_saturation_model", PROJECT_ROOT)
    
    def saturation_factor(self, density_g_cm3):
        """
        Calculate saturation factor based on density.
        
        From Paper 14: Effects saturate when ρ approaches ρ_c.
        The saturation factor S_sat = 1 / (1 + (ρ/ρ_c)^γ)
        where γ controls the saturation strength.
        
        For ρ << ρ_c: S_sat ≈ 1 (no saturation)
        For ρ >> ρ_c: S_sat ≈ (ρ_c/ρ)^γ (saturated)
        """
        if density_g_cm3 <= 0:
            return 1.0
        
        # Normalized density
        rho_norm = density_g_cm3 / RHO_T
        
        # Saturation strength (calibrated from pulsar data)
        gamma = 1.0  # Can be adjusted based on calibration
        
        # Saturation factor
        saturation_factor = 1.0 / (1.0 + rho_norm**gamma)
        
        return saturation_factor
    
    def local_density_at_altitude(self, altitude_km):
        """
        Calculate local density at given altitude.
        
        Uses piecewise exponential atmosphere model (NRLMSISE-00 approx).
        """
        if altitude_km < 0:
            return 5.515  # Mean Earth density (g/cm³)
            
        if altitude_km <= 100:
            return 1.225e-3 * np.exp(-altitude_km / 8.5)
        elif altitude_km <= 300:
            return 9.3e-12 * np.exp(-(altitude_km - 100) / 40.0)
        elif altitude_km <= 600:
            return 6.3e-14 * np.exp(-(altitude_km - 300) / 80.0)
        elif altitude_km <= 1000:
            return 1.5e-15 * np.exp(-(altitude_km - 600) / 150.0)
        else:
            return max(1.0e-17 * np.exp(-(altitude_km - 1000) / 300.0), 1e-20)
    
    def apply_saturation_to_beta(self, beta, altitude_km):
        """
        Apply saturation correction to effective coupling.
        
        β_sat = β × S_sat(ρ)
        
        This accounts for the suppressed density scaling observed
        in pulsar data (Paper 14).
        """
        density = self.local_density_at_altitude(altitude_km)
        saturation_factor = self.saturation_factor(density)
        
        return beta * saturation_factor
    
    def analyze_flyby_saturation(self):
        """
        Analyze saturation effects for all flybys.
        
        This determines which flybys are in the saturated regime
        and quantifies the expected reduction in β_eff.
        """
        self.logger.section("SATURATION MODEL ANALYSIS")
        
        # Load flyby data
        predictions_file = PROJECT_ROOT / 'results' / 'step004_tep_predictions.json'
        
        if not predictions_file.exists():
            self.logger.error("TEP predictions not found. Run step004 first.")
            return None
        
        try:
            with open(predictions_file) as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load predictions: {e}")
            return 1
        
        self.logger.subsection("FLYBY SATURATION ANALYSIS")
        
        results = []
        for name, pred in data['predictions'].items():
            altitude = pred['geometry']['altitude_km']
            beta_eff = pred['tep_predictions'].get('beta_eff', pred['tep_predictions'].get('beta_initial', 1e-4))
            
            # Calculate density and saturation
            density = self.local_density_at_altitude(altitude)
            saturation_factor = self.saturation_factor(density)
            beta_sat = beta_eff * saturation_factor
            
            self.logger.info(f"{name}:")
            self.logger.info(f"  Altitude: {altitude:.0f} km")
            self.logger.info(f"  Density: {density:.2e} g/cm³")
            self.logger.info(f"  Saturation factor: {saturation_factor:.2e}")
            self.logger.info(f"  β_eff (no saturation): {beta_eff:.2e}")
            self.logger.info(f"  β_eff (with saturation): {beta_sat:.2e}")
            
            results.append({
                'name': name,
                'altitude_km': altitude,
                'density_g_cm3': density,
                'saturation_factor': saturation_factor,
                'beta_eff_original': beta_eff,
                'beta_eff_saturated': beta_sat
            })
        
        # Summary statistics
        self.logger.subsection("SUMMARY STATISTICS")
        
        avg_saturation = np.mean([r['saturation_factor'] for r in results])
        self.logger.info(f"Average saturation factor: {avg_saturation:.2e}")
        
        saturated_count = sum(1 for r in results if r['saturation_factor'] < 0.5)
        self.logger.info(f"Flybys in saturated regime (S < 0.5): {saturated_count}/{len(results)}")
        
        return {
            'flyby_results': results,
            'average_saturation_factor': float(avg_saturation),
            'saturated_count': saturated_count,
            'total_count': len(results)
        }


def main():
    """Execute saturation model analysis."""
    logger = StepLogger("step_012_saturation_model", PROJECT_ROOT)
    logger.section("STEP 012: SATURATION MODEL (PAPER 14)")
    
    model = SaturationModel()
    
    # Analyze flyby saturation
    results = model.analyze_flyby_saturation()
    
    if results is None:
        logger.log_step_summary(0, "FAILED")
        return 1
    
    # Save results
    results_dir = PROJECT_ROOT / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = results_dir / 'step012_saturation_model.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.section("SAVING RESULTS")
    logger.info(f"Results saved to: {output_file}")
    
    logger.log_step_summary(0, "SUCCESS")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
