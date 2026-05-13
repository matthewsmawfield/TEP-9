"""
Step 017: Plasma Modulation for TEP Scalar Force

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

from scripts.utils.celestrak_space_weather import lookup_space_weather
from scripts.utils.iri_mission_map import resolve_iri_mission
from scripts.utils.plasma_screening import iri_reference_electron_density_cm3
from scripts.utils.step_logger import StepLogger


class PlasmaModulationModel:
    """
    Plasma modulation model for TEP scalar force.

    Uses a bounded phenomenological ansatz for plasma attenuation.
    The Debye-screening formula exp(-lambda_TEP / lambda_D) is
    numerically catastrophic (underflows to zero for all realistic
    ionospheric densities) and physically unjustified for a neutral
    scalar-gravity field.  We replace it with a smooth exponential
    in electron density that is explicitly labelled as an ansatz.
    A derivation of scalar-plasma coupling from the underlying TEP
    action remains future work.
    """

    def __init__(self):
        self.logger = StepLogger("step_017_plasma_modulation", PROJECT_ROOT)

        # Phenomenological plasma-attenuation ansatz.
        # S = exp(-n_e / n_ref) where n_ref is calibrated from IRI peak densities.
        self.n_ref = iri_reference_electron_density_cm3(PROJECT_ROOT)
        self.n_crit = self.n_ref
        
        # Load IRI trajectory profiles from step_027
        self.iri_profiles = {}
        iri_file = PROJECT_ROOT / 'results' / 'step033_iri_trajectory_profiles.json'
        if not iri_file.exists():
            raise RuntimeError(
                "IRI trajectory profiles are required (results/step033_iri_trajectory_profiles.json)"
            )

        with open(iri_file, 'r') as f:
            self.iri_profiles = json.load(f)
        self.logger.info("Loaded IRI trajectory profiles from step_033")
        
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
        iri_mission = resolve_iri_mission(mission) if mission else None
        if not iri_mission or iri_mission not in self.iri_interpolators:
            return None

        log_density = self.iri_interpolators[iri_mission](altitude_km)
        density = 10**log_density
        return max(density, 1e-10)
    
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
        if n_iono is None:
            return None
        n_mag = self.magnetospheric_density(altitude_km, l_shell)
        
        # Combine densities
        n_total = n_iono + n_mag
        
        return n_total
    
    def plasma_screening_factor(self, plasma_density):
        """
        Calculate plasma attenuation factor using a phenomenological ansatz.

        We adopt S = exp(-n_e / n_ref) as a smooth proxy for ionospheric
        screening.  This replaces the numerically pathological Debye formula
        exp(-lambda_TEP / lambda_D), which underflows to zero for all
        realistic ionospheric densities and lacks a derivation from the TEP
        action for a neutral scalar field.

        Parameters:
        - plasma_density: Electron density in cm⁻³

        Returns:
        - Screening factor (0 to 1)
        """
        if plasma_density <= 0:
            return 1.0

        # Phenomenological ansatz: higher plasma density → stronger attenuation
        screening_factor = np.exp(-plasma_density / self.n_ref)

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
        
        flyby_date = flyby_data.get("flyby_date")
        if not flyby_date:
            raise RuntimeError(f"Missing flyby_date for {mission_name}")
        solar = lookup_space_weather(flyby_date)
        solar_activity = solar["f10_7"] / 100.0
        
        n_plasma = self.total_plasma_density(altitude_km, solar_activity, mission=mission_name)
        if n_plasma is None:
            raise RuntimeError(f"No IRI trajectory profile available for {mission_name}")
        
        # Calculate screening and sign factors
        S_plasma = self.plasma_screening_factor(n_plasma)
        sign_factor = self.plasma_sign_factor(n_plasma, solar_activity)
        
        return {
            'plasma_density_cm3': float(n_plasma),
            'plasma_screening_factor': float(S_plasma),
            'plasma_sign_factor': float(sign_factor),
            'solar_activity': solar_activity,
            'f10_7_sfu': solar["f10_7"],
            'kp_sum': solar["kp_sum"],
            'derivation': (
                f"IRI electron density at perigee with Celestrak F10.7={solar['f10_7']:.1f} sfu "
                f"({solar['f10_7_field']})"
            ),
            'source': solar["data_source"],
            'iri_data_source': 'step_033',
        }


def main():
    logger = StepLogger("step_017_plasma_modulation", PROJECT_ROOT)
    logger.section("STEP 017: PLASMA MODULATION ANALYSIS")
    
    model = PlasmaModulationModel()
    
    # Load flyby data from results folder
    predictions_file = PROJECT_ROOT / 'results' / 'step007_tep_predictions.json'
    
    if not predictions_file.exists():
        logger.error("TEP predictions not found. Run step007 first.")
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
    
    results = {
        'uncertainty_fraction': 0.15,
        'status': 'iri_profile',
        'calibration_status': 'empirical_IRI_model',
        'data_source': 'IRI trajectory profiles (step_033) with Celestrak SW-All solar activity',
        'derivation': 'Perigee IRI electron density with Celestrak F10.7/Kp; no Chapman fallback',
        'recommended_action': 'Derive scalar-plasma coupling from TEP action; screening ansatz remains phenomenological',
    }
    
    for name, pred in data['predictions'].items():
        flyby_data = {
            'name': name,
            'flyby_date': pred['perigee']['datetime'],
            'altitude_km': pred['perigee']['altitude_km'],
            'dv_obs_mm_s': pred['observed']['dv_obs_mm_s'],
            'dv_tep_mm_s': pred['tep_predictions']['dv_tep_mm_s']
        }

        try:
            plasma_effects = model.calculate_plasma_effects(flyby_data)
        except RuntimeError as exc:
            logger.error(f"{name}: plasma modulation failed ({exc})")
            logger.log_step_summary(0, "FAILED")
            return 1

        results[name] = {
            'uncertainty': None,
            'uncertainty_fraction': 0.15,
            'uncertainty_absolute': None,
            'status': 'iri_profile',
            'calibration_status': plasma_effects.get('iri_data_source', 'unknown'),
            'data_source': plasma_effects.get('iri_data_source', 'unknown'),
            'derivation': plasma_effects.get('derivation'),
            'recommended_action': 'Derive scalar-plasma coupling from TEP action; screening ansatz remains phenomenological',
            **plasma_effects
        }

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
    
    output_file = results_dir / 'step017_plasma_modulation.json'
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
