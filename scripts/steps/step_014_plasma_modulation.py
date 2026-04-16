"""
Plasma Modulation for TEP Scalar Force

This module implements plasma-dependent screening effects to address
the Cassini sign mismatch and other anomalies in flyby data.

Key features:
- Plasma density models along flyby trajectories
- Plasma-dependent screening factor
- Sign-reversal mechanisms in plasma
- Ionospheric and magnetospheric effects
"""

import numpy as np
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


class PlasmaModulationModel:
    """
    Plasma modulation model for TEP scalar force.
    
    The scalar field can be screened or enhanced by plasma density,
    potentially causing sign reversals in certain regimes.
    """
    
    def __init__(self):
        self.logger = StepLogger("step_014_plasma_modulation", PROJECT_ROOT)
        
        # Plasma density thresholds (cm^-3)
        self.n_crit = 1e4  # Critical density for screening reversal
        self.n_screen = 1e3  # Density where screening becomes significant
        
    def ionospheric_density(self, altitude_km, solar_activity=1.0):
        """
        Calculate ionospheric electron density at given altitude.
        
        Uses Chapman layer model simplified for Earth's ionosphere.
        """
        # Chapman layer parameters
        h_max = 300  # km - F2 layer peak
        scale_height = 50  # km
        n_max = 1e6 * solar_activity  # cm^-3 - peak density
        
        # Chapman function
        z = (altitude_km - h_max) / scale_height
        density = n_max * np.exp(0.5 * (1 - z - np.exp(-z)))
        
        return max(density, 1e-10)  # Minimum density
    
    def magnetospheric_density(self, altitude_km, l_shell=1.0):
        """
        Calculate magnetospheric plasma density at given altitude.
        
        Uses empirical model for plasma density in magnetosphere.
        """
        if altitude_km < 1000:
            # Inside plasmasphere
            n_eq = 1e3  # cm^-3 - equatorial density at L=1
            density = n_eq / (l_shell**4) * np.exp(-(altitude_km - 300) / 1000)
        else:
            # Outside plasmasphere
            n_eq = 0.1  # cm^-3 - much lower
            density = n_eq / (l_shell**4)
        
        return max(density, 1e-10)
    
    def total_plasma_density(self, altitude_km, solar_activity=1.0, l_shell=1.0):
        """
        Calculate total plasma density (ionospheric + magnetospheric).
        """
        n_iono = self.ionospheric_density(altitude_km, solar_activity)
        n_mag = self.magnetospheric_density(altitude_km, l_shell)
        
        # Combine densities
        n_total = n_iono + n_mag
        
        return n_total
    
    def plasma_screening_factor(self, plasma_density):
        """
        Calculate plasma screening factor.
        
        In high-density plasma, the scalar field can be screened
        or even experience sign reversal due to charge screening.
        
        S_plasma = 1 / (1 + (n/n_crit)^2)
        
        For n >> n_crit: strong screening (S << 1)
        For n << n_crit: no screening (S ≈ 1)
        """
        ratio = plasma_density / self.n_crit
        screening_factor = 1.0 / (1.0 + ratio**2)
        
        return screening_factor
    
    def plasma_sign_factor(self, plasma_density, solar_activity=1.0):
        """
        Calculate potential sign reversal factor.
        
        In certain plasma regimes, the scalar force can reverse sign
        due to charge screening effects or disformal coupling.
        
        Returns factor in [-1, 1] where negative indicates sign reversal.
        """
        ratio = plasma_density / self.n_crit
        
        # Sign reversal occurs at intermediate densities
        # Enhanced model includes solar activity dependence
        # Higher solar activity increases ionospheric density, promoting sign reversal
        
        # Critical density threshold depends on solar activity
        n_crit_effective = self.n_crit / solar_activity
        ratio_effective = plasma_density / n_crit_effective
        
        if ratio_effective < 0.05:
            sign_factor = 1.0
        elif ratio_effective < 0.5:
            # Transition region where sign can flip
            # Smooth sigmoid transition
            sign_factor = 1.0 - 2.0 * (ratio_effective - 0.05) / 0.45
        elif ratio_effective < 2.0:
            # Full sign reversal region
            sign_factor = -1.0
        else:
            # Very high density - may return to positive (screening dominates)
            sign_factor = -1.0 + 0.5 * (ratio_effective - 2.0) / 10.0
        
        return max(min(sign_factor, 1.0), -1.0)  # Clamp to [-1, 1]
    
    def calculate_plasma_effects(self, flyby_data):
        """
        Calculate plasma modulation effects for a flyby.
        
        Returns:
        - plasma_density_profile along trajectory
        - plasma_screening_factor
        - plasma_sign_factor
        """
        altitude_km = flyby_data['altitude_km']
        
        # Estimate solar activity based on date (simplified)
        # Cassini 1999: solar minimum
        # NEAR 1998: solar minimum
        # Galileo 1990: solar maximum
        solar_activity = self.estimate_solar_activity(flyby_data.get('name', ''))
        
        # Calculate plasma density
        n_plasma = self.total_plasma_density(altitude_km, solar_activity)
        
        # Calculate screening and sign factors
        S_plasma = self.plasma_screening_factor(n_plasma)
        sign_factor = self.plasma_sign_factor(n_plasma, solar_activity)
        
        return {
            'plasma_density_cm3': float(n_plasma),
            'plasma_screening_factor': float(S_plasma),
            'plasma_sign_factor': float(sign_factor),
            'solar_activity': solar_activity
        }
    
    def estimate_solar_activity(self, spacecraft_name):
        """
        Estimate solar activity level based on flyby date.
        
        Returns: solar_activity factor (1.0 = average, <1 = minimum, >1 = maximum)
        """
        # Simplified solar cycle estimation
        if 'Galileo_1990' in spacecraft_name:
            return 1.5  # Solar maximum
        elif 'Cassini' in spacecraft_name or 'NEAR' in spacecraft_name:
            return 0.5  # Solar minimum
        else:
            return 1.0  # Average


def main():
    """Execute plasma modulation analysis."""
    logger = StepLogger("step_014_plasma_modulation", PROJECT_ROOT)
    logger.section("STEP 014: PLASMA MODULATION ANALYSIS")
    
    model = PlasmaModulationModel()
    
    # Load flyby data from results folder
    predictions_file = PROJECT_ROOT / 'results' / 'step004_tep_predictions.json'
    
    if not predictions_file.exists():
        logger.error("TEP predictions not found. Run step004 first.")
        logger.log_step_summary(0, "FAILED")
        return 1
    
    with open(predictions_file) as f:
        data = json.load(f)
    
    logger.subsection("PLASMA DENSITY ANALYSIS")
    
    results = {}
    
    for name, pred in data['predictions'].items():
        if pred['observed']['dv_obs_mm_s'] != 0:
            flyby_data = {
                'name': name,
                'altitude_km': pred['perigee']['altitude_km'],
                'dv_obs_mm_s': pred['observed']['dv_obs_mm_s'],
                'dv_tep_mm_s': pred['tep_predictions']['dv_tep_mm_s']
            }
            
            plasma_effects = model.calculate_plasma_effects(flyby_data)
            results[name] = plasma_effects
            
            logger.info(f"{name}:")
            logger.info(f"  Altitude: {flyby_data['altitude_km']:.0f} km")
            logger.info(f"  Plasma density: {plasma_effects['plasma_density_cm3']:.2e} cm^-3")
            logger.info(f"  Plasma screening factor: {plasma_effects['plasma_screening_factor']:.3e}")
            logger.info(f"  Plasma sign factor: {plasma_effects['plasma_sign_factor']:.3f}")
            logger.info(f"  Solar activity: {plasma_effects['solar_activity']:.1f}")
    
    # Special analysis for Cassini sign mismatch
    if 'Cassini_1999' in results:
        logger.subsection("CASSINI SIGN MISMATCH ANALYSIS")
        cassini = results['Cassini_1999']
        
        logger.info("Cassini observed: +0.11 mm/s")
        logger.info("Cassini predicted (no plasma): -0.19 mm/s")
        logger.info(f"Plasma sign factor: {cassini['plasma_sign_factor']:.3f}")
        
        if cassini['plasma_sign_factor'] < 0:
            logger.info("CONCLUSION: Plasma modulation can explain sign reversal")
        else:
            logger.info("CONCLUSION: Plasma modulation insufficient for sign reversal")
    
    # Save results
    results_dir = PROJECT_ROOT / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = results_dir / 'step014_plasma_modulation_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.section("SAVING RESULTS")
    logger.info(f"Results saved to: {output_file}")
    
    logger.log_step_summary(0, "SUCCESS")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
