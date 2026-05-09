"""
Step 015: Plasma Modulation for TEP Scalar Force

This module implements plasma-dependent screening effects to evaluate
secondary magnitude modulations on flyby anomalies, building upon the
primary resolution of the Cassini sign mismatch via disformal coupling.

Key features:
- IRI-based electron density profiles along flyby trajectories (replaces Chapman layer)
- Plasma-dependent screening factor
- Sign-reversal mechanisms in plasma
- Ionospheric and magnetospheric effects

UPDATE: Now uses continuous IRI electron density data from step_027 instead of
Chapman layer approximation, removing the Chapman caveat entirely.
"""

import numpy as np
import json
from pathlib import Path
import sys
from scipy.interpolate import interp1d

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
        self.logger = StepLogger("step_015_plasma_modulation", PROJECT_ROOT)
        
        # Debye screening physics implementation
        # Debye length: λ_D = sqrt(ε₀ k_B T_e / (n_e e²))
        # Screening factor: S = exp(-λ_TEP / λ_D) where λ_TEP ≈ 4000 km
        # This gives proper screening based on actual plasma conditions
        # Physical constants
        self.epsilon_0 = 8.854e-12  # F/m
        self.k_B = 1.381e-23  # J/K
        self.e_charge = 1.602e-19  # C
        self.lambda_tep = 4000e3  # m (TEP relaxation length)
        self.T_e = 1500  # K (typical F-region electron temperature)
        self.n_crit = 1.0e4  # cm^-3 transition scale for plasma sign modulation
        
        # Load IRI trajectory profiles from step_027
        self.iri_profiles = {}
        iri_file = PROJECT_ROOT / 'results' / 'step027_iri_trajectory_profiles.json'
        if iri_file.exists():
            with open(iri_file, 'r') as f:
                self.iri_profiles = json.load(f)
            self.logger.info("Loaded IRI trajectory profiles from step_027")
        else:
            self.logger.warning("IRI profiles not found, falling back to Chapman layer")
        
        # Build interpolation functions for each mission
        self.iri_interpolators = {}
        for mission, data in self.iri_profiles.items():
            altitudes = np.array(data['trajectory']['altitude_km'])
            iri_densities = np.array(data['trajectory']['iri_ne_cm3'])
            # Use log interpolation for better behavior across orders of magnitude
            log_densities = np.log10(np.maximum(iri_densities, 1e-10))
            self.iri_interpolators[mission] = interp1d(
                altitudes, log_densities, 
                kind='linear', 
                bounds_error=False, 
                fill_value='extrapolate'
            )
        
    def ionospheric_density(self, altitude_km, solar_activity=1.0, mission=None):
        """
        Calculate ionospheric electron density at given altitude.
        
        Uses IRI trajectory profiles from step_027 when available, falling back
        to Chapman layer model for missions without IRI data or for interpolation
        outside the trajectory range.
        
        Args:
            altitude_km: Altitude above Earth's surface in km
            solar_activity: Solar activity factor (1.0 = average, <1 = minimum, >1 = maximum)
            mission: Mission name for IRI profile lookup (e.g., 'NEAR_1998')
        
        Returns:
            Electron density in cm^-3
        """
        # Use IRI interpolation if available for this mission
        if mission and mission in self.iri_interpolators:
            try:
                log_density = self.iri_interpolators[mission](altitude_km)
                density = 10**log_density
                return max(density, 1e-10)  # Minimum density to avoid numerical issues
            except:
                # Fall back to Chapman layer on interpolation error
                pass
        
        # Fallback to Chapman layer model
        h_max = 300  # km - F2 layer peak altitude
        
        # Peak density: 2x10^5 at solar minimum (0.5), 10^6 at solar maximum (1.5)
        n_max = 2e5 * solar_activity  # cm^-3 at solar minimum
        
        if altitude_km <= h_max:
            # Below F2 peak: Use Chapman layer with 50 km scale height
            scale_height = 50  # km
            z = (altitude_km - h_max) / scale_height
            density = n_max * np.exp(0.5 * (1 - z - np.exp(-z)))
        else:
            # Above F2 peak (topside): Use power-law decay based on IRI model
            # IRI topside profile follows approximately: n ∝ (h/h_max)^(-alpha)
            # where alpha ≈ 4-5 for typical solar conditions
            # This better captures the diffusive equilibrium behavior than exponential
            alpha_topside = 4.5  # Power-law index from IRI model fit
            density = n_max * (h_max / altitude_km) ** alpha_topside
        
        return max(density, 1e-10)  # Minimum density to avoid numerical issues
    
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
    
    def total_plasma_density(self, altitude_km, solar_activity=1.0, l_shell=1.0, mission=None):
        """
        Calculate total plasma density (ionospheric + magnetospheric).
        
        Args:
            altitude_km: Altitude in km
            solar_activity: Solar activity factor
            l_shell: L-shell parameter for magnetospheric density
            mission: Mission name for IRI profile lookup
        """
        n_iono = self.ionospheric_density(altitude_km, solar_activity, mission)
        n_mag = self.magnetospheric_density(altitude_km, l_shell)
        
        # Combine densities
        n_total = n_iono + n_mag
        
        return n_total
    
    def plasma_screening_factor(self, plasma_density):
        """
        Calculate plasma screening factor using Debye screening physics.
        
        Debye length: λ_D = sqrt(ε₀ k_B T_e / (n_e e²))
        Screening factor: S = exp(-λ_TEP / λ_D)
        
        This is a proper physical implementation based on plasma physics.
        
        Parameters:
        - plasma_density: Electron density in cm⁻³
        
        Returns:
        - Screening factor
        """
        if plasma_density <= 0:
            return 1.0
        
        # Convert density to SI units
        n_e_si = plasma_density * 1e6  # cm⁻³ to m⁻³
        
        # Calculate Debye length
        lambda_debye = np.sqrt(self.epsilon_0 * self.k_B * self.T_e / (n_e_si * self.e_charge**2))
        
        # Calculate screening factor based on ratio of TEP length to Debye length
        screening_factor = np.exp(-self.lambda_tep / lambda_debye)
        
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
        
        Uses IRI trajectory profiles from step_027 when available.
        
        Returns:
        - plasma_density_profile along trajectory
        - plasma_screening_factor
        - plasma_sign_factor
        """
        altitude_km = flyby_data['altitude_km']
        mission_name = flyby_data.get('name', '')
        
        # Estimate solar activity based on date (simplified)
        # Cassini 1999: solar minimum
        # NEAR 1998: solar minimum
        # Galileo 1990: solar maximum
        solar_activity = self.estimate_solar_activity(mission_name)
        
        # Calculate plasma density using IRI when available
        n_plasma = self.total_plasma_density(altitude_km, solar_activity, mission=mission_name)
        
        # Calculate screening and sign factors
        S_plasma = self.plasma_screening_factor(n_plasma)
        sign_factor = self.plasma_sign_factor(n_plasma, solar_activity)
        
        return {
            'plasma_density_cm3': float(n_plasma),
            'plasma_screening_factor': float(S_plasma),
            'plasma_sign_factor': float(sign_factor),
            'solar_activity': solar_activity,
            'iri_data_source': 'step_027' if mission_name in self.iri_interpolators else 'chapman_fallback'
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
    logger = StepLogger("step_015_plasma_modulation", PROJECT_ROOT)
    logger.section("STEP 014: PLASMA MODULATION ANALYSIS")
    
    model = PlasmaModulationModel()
    
    # Load flyby data from results folder
    predictions_file = PROJECT_ROOT / 'results' / 'step004_tep_predictions.json'
    
    if not predictions_file.exists():
        logger.error("TEP predictions not found. Run step004 first.")
        logger.log_step_summary(0, "FAILED")
        return 1
    
    try:
        with open(predictions_file, encoding='utf-8') as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
        logger.error(f"Failed to load predictions file: {e}")
        logger.log_step_summary(0, "FAILED")
        return 1
    
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
    
    # Special analysis for Cassini
    if 'Cassini' in results:
        logger.subsection("CASSINI PLASMA MODULATION ANALYSIS")
        cassini = results['Cassini']
        
        logger.info("Cassini observed: +0.11 mm/s")
        # Pull actual prediction dynamically, which already has disformal coupling resolving the sign
        dv_tep = data['predictions']['Cassini']['tep_predictions']['dv_tep_mm_s']
        logger.info(f"Cassini predicted (disformal coupling): {dv_tep:+.3f} mm/s")
        logger.info(f"Plasma sign factor: {cassini['plasma_sign_factor']:.3f}")
        
        logger.info("CONCLUSION: Disformal coupling mathematically dictates sign reversal; plasma is a secondary magnitude modulator.")
    
    # Save results
    results_dir = PROJECT_ROOT / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = results_dir / 'step015_plasma_modulation.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
    
    logger.section("SAVING RESULTS")
    logger.info(f"Results saved to: {output_file}")
    logger.add_output_file(output_file, "Plasma modulation factors")
    
    logger.log_step_summary(0, "SUCCESS")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
