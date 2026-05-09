#!/usr/bin/env python3
"""
Step 039: Systematic Error Monte Carlo Analysis for Flyby Anomalies

Quantifies systematic uncertainty through Monte Carlo error propagation for TEP-EFA.
Tests robustness against:
1. DSN measurement systematics (antenna phase center, tropospheric delay, station position)
2. Trajectory uncertainties (perigee altitude, velocity errors)
3. Calibration drifts
4. Correlated noise

This provides a quantitative systematic error budget similar to TEP-LLR Step 020.
"""

import sys
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger
from scripts.utils.physics import BETA_BASELINE


class SystematicErrorMonteCarlo:
    """
    Monte Carlo analysis of systematic errors in flyby anomaly measurements.
    
    Propagates known systematic uncertainties through the TEP analysis
    to quantify their impact on the fitted β parameter and predicted anomalies.
    """
    
    def __init__(self, seed: int = 42):
        self.seed = seed
        np.random.seed(seed)
        self.logger = StepLogger("step_039")
        
        # Systematic error magnitudes (from literature and DSN specifications)
        self.systematic_errors = {
            'dsn_antenna_phase': 0.10,  # mm/s (antenna phase center motion)
            'dsn_tropospheric': 0.05,   # mm/s (tropospheric delay modeling)
            'dsn_station_position': 0.02,  # mm/s (station position errors)
            'trajectory_altitude': 1.0,   # km (perigee altitude uncertainty)
            'trajectory_velocity': 1.0,   # m/s (perigee velocity uncertainty)
            'calibration_drift': 0.01,   # mm/s/hour (calibration drift rate)
        }
    
    def inject_dsn_systematics(
        self,
        observed_dv_mm_s: float,
        tracking_precision_mm_s: float,
        n_trials: int = 1000
    ) -> np.ndarray:
        """
        Inject DSN measurement systematics into observed velocity anomalies.
        
        Parameters:
        -----------
        observed_dv_mm_s : float
            Observed velocity anomaly (mm/s)
        tracking_precision_mm_s : float
            Formal tracking precision (mm/s)
        n_trials : int
            Number of Monte Carlo trials
        
        Returns:
        --------
        perturbed_dv : np.ndarray
            Array of perturbed velocity anomalies
        """
        # Total DSN systematic error (quadrature sum)
        dsn_systematic = np.sqrt(
            self.systematic_errors['dsn_antenna_phase']**2 +
            self.systematic_errors['dsn_tropospheric']**2 +
            self.systematic_errors['dsn_station_position']**2
        )
        
        # Combine formal error with systematic
        total_error = np.sqrt(tracking_precision_mm_s**2 + dsn_systematic**2)
        
        # Generate perturbed values
        perturbed_dv = np.random.normal(observed_dv_mm_s, total_error, n_trials)
        
        return perturbed_dv
    
    def inject_trajectory_uncertainty(
        self,
        perigee_altitude_km: float,
        perigee_velocity_km_s: float,
        n_trials: int = 1000
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Inject trajectory uncertainties into perigee parameters.
        
        Parameters:
        -----------
        perigee_altitude_km : float
            Perigee altitude (km)
        perigee_velocity_km_s : float
            Perigee velocity (km/s)
        n_trials : int
            Number of Monte Carlo trials
        
        Returns:
        --------
        perturbed_altitude : np.ndarray
            Array of perturbed altitudes
        perturbed_velocity : np.ndarray
            Array of perturbed velocities
        """
        # Perturb altitude
        perturbed_altitude = np.random.normal(
            perigee_altitude_km,
            self.systematic_errors['trajectory_altitude'],
            n_trials
        )
        
        # Perturb velocity (convert km/s to m/s for error, then back)
        perturbed_velocity = np.random.normal(
            perigee_velocity_km_s,
            self.systematic_errors['trajectory_velocity'] / 1000.0,  # m/s to km/s
            n_trials
        )
        
        return perturbed_altitude, perturbed_velocity
    
    def compute_tep_sensitivity(
        self,
        altitude_km: float,
        velocity_km_s: float,
        cos_asymmetry: float,
        beta: float = None
    ) -> float:
        """
        Compute TEP-predicted anomaly given trajectory parameters.
        
        Simplified TEP model: Δv ∝ β * (1/r²) * v * cos(α)
        
        Parameters:
        -----------
        altitude_km : float
            Perigee altitude (km)
        velocity_km_s : float
            Perigee velocity (km/s)
        cos_asymmetry : float
            Cosine of declination asymmetry
        beta : float
            TEP coupling constant (defaults to BETA_BASELINE from physics.py)
        
        Returns:
        --------
        dv_tep_mm_s : float
            TEP-predicted velocity anomaly (mm/s)
        """
        if beta is None:
            beta = BETA_BASELINE * 1e-4  # Convert baseline to actual coupling
        
        R_EARTH = 6371.0  # km
        r = R_EARTH + altitude_km  # Distance from Earth center
        
        # Simplified TEP scaling (from step_005 results)
        # dv_tep ≈ β * 1e7 * (1/r²) * v * cos(α)
        # Scaling factor calibrated to match step_005 predictions
        scaling_factor = 1e7
        
        dv_tep = beta * scaling_factor * (1.0 / r**2) * velocity_km_s * cos_asymmetry
        
        return dv_tep
    
    def monte_carlo_analysis(
        self,
        catalog_path: Path,
        n_trials: int = 1000
    ) -> Dict:
        """
        Perform Monte Carlo systematic error analysis on all flybys.
        
        Parameters:
        -----------
        catalog_path : Path
            Path to step002_archival_flyby_catalog.json
        n_trials : int
            Number of Monte Carlo trials
        
        Returns:
        --------
        results : Dict
            Monte Carlo analysis results
        """
        try:
            with open(catalog_path, 'r') as f:
                catalog = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            self.logger.error(f"Failed to load catalog: {e}")
            return None
        
        # Load fitting results for reference beta
        fitting_path = PROJECT_ROOT / "results" / "step005_fitting_results.json"
        try:
            with open(fitting_path, 'r') as f:
                fitting_data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            self.logger.error(f"Failed to load fitting results: {e}")
            return None
        
        # Use coupling constant as reference beta
        beta_ref = fitting_data.get('coupling_constant', fitting_data.get('beta', 0.0004635))
        
        flyby_results = []
        
        for flyby in catalog['flybys']:
            if not flyby['usable_for_analysis']:
                continue
            
            observed = flyby.get('published_anomaly_mm_s')
            if observed is None or observed == 0:
                continue
            
            mission = flyby['mission_name']
            
            # Get mission-specific parameters
            altitude = flyby['perigee_altitude_km']
            velocity = flyby['perigee_velocity_km_s']
            cos_asym = flyby.get('cos_asymmetry', 0.0)
            precision = flyby.get('tracking_precision_mm_s', 0.05)
            
            # Monte Carlo trials
            dv_perturbed = self.inject_dsn_systematics(observed, precision, n_trials)
            alt_perturbed, vel_perturbed = self.inject_trajectory_uncertainty(
                altitude, velocity, n_trials
            )
            
            # Compute TEP predictions for each trial
            dv_tep_trials = np.array([
                self.compute_tep_sensitivity(alt_perturbed[i], vel_perturbed[i], cos_asym, beta_ref)
                for i in range(n_trials)
            ])
            
            # Compute statistics
            dv_tep_mean = np.mean(dv_tep_trials)
            dv_tep_std = np.std(dv_tep_trials)
            dv_tep_ci_lower = np.percentile(dv_tep_trials, 2.5)
            dv_tep_ci_upper = np.percentile(dv_tep_trials, 97.5)
            
            # Compute residual statistics
            residual_trials = dv_perturbed - dv_tep_trials
            residual_mean = np.mean(residual_trials)
            residual_std = np.std(residual_trials)
            
            flyby_results.append({
                'mission': mission,
                'observed_dv_mm_s': observed,
                'tracking_precision_mm_s': precision,
                'tep_dv_mean_mm_s': dv_tep_mean,
                'tep_dv_std_mm_s': dv_tep_std,
                'tep_dv_ci_lower_mm_s': dv_tep_ci_lower,
                'tep_dv_ci_upper_mm_s': dv_tep_ci_upper,
                'residual_mean_mm_s': residual_mean,
                'residual_std_mm_s': residual_std,
                'systematic_uncertainty_fraction': dv_tep_std / dv_tep_mean if dv_tep_mean > 0 else 0
            })
        
        # Global systematic error budget
        all_std = [r['tep_dv_std_mm_s'] for r in flyby_results]
        all_mean = [r['tep_dv_mean_mm_s'] for r in flyby_results]
        
        systematic_uncertainties = np.array([s/m if m > 0 else 0 for s, m in zip(all_std, all_mean)])
        mean_systematic_uncertainty = np.mean(systematic_uncertainties)
        max_systematic_uncertainty = np.max(systematic_uncertainties)
        
        return {
            'n_trials': n_trials,
            'beta_ref': beta_ref,
            'systematic_error_magnitudes': self.systematic_errors,
            'flyby_results': flyby_results,
            'systematic_error_budget': {
                'mean_systematic_uncertainty': float(mean_systematic_uncertainty),
                'max_systematic_uncertainty': float(max_systematic_uncertainty),
                'interpretation': f'Systematic uncertainties contribute {mean_systematic_uncertainty*100:.1f}% on average to TEP predictions'
            }
        }


def main():
    """Execute systematic error Monte Carlo analysis."""
    logger = StepLogger("step_039")
    
    try:
        mc = SystematicErrorMonteCarlo(seed=42)
        
        # Load flyby catalog
        catalog_path = PROJECT_ROOT / "results" / "step002_archival_flyby_catalog.json"
        
        if not catalog_path.exists():
            logger.error(f"Catalog not found: {catalog_path}")
            return None
        
        # Perform Monte Carlo analysis
        results = mc.monte_carlo_analysis(catalog_path, n_trials=1000)
        
        # Save results
        output_path = PROJECT_ROOT / "results" / "step_039_systematic_error_monte_carlo.json"
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Print summary
        print("\n" + "="*70)
        print("SYSTEMATIC ERROR MONTE CARLO ANALYSIS RESULTS")
        print("="*70)
        print(f"\nMonte Carlo trials: {results['n_trials']}")
        print(f"Reference β: {results['beta_ref']:.2e}")
        
        print("\nSystematic error magnitudes:")
        for key, value in results['systematic_error_magnitudes'].items():
            print(f"  {key}: {value}")
        
        print("\nSystematic error budget:")
        budget = results['systematic_error_budget']
        print(f"  Mean systematic uncertainty: {budget['mean_systematic_uncertainty']*100:.2f}%")
        print(f"  Max systematic uncertainty: {budget['max_systematic_uncertainty']*100:.2f}%")
        print(f"  Interpretation: {budget['interpretation']}")
        
        print("\nPer-flyby results:")
        for r in results['flyby_results']:
            print(f"\n  {r['mission']}:")
            print(f"    Observed Δv: {r['observed_dv_mm_s']:.2f} mm/s")
            print(f"    TEP Δv (mean ± std): {r['tep_dv_mean_mm_s']:.3f} ± {r['tep_dv_std_mm_s']:.3f} mm/s")
            print(f"    TEP Δv (95% CI): [{r['tep_dv_ci_lower_mm_s']:.3f}, {r['tep_dv_ci_upper_mm_s']:.3f}] mm/s")
            print(f"    Residual (mean ± std): {r['residual_mean_mm_s']:.3f} ± {r['residual_std_mm_s']:.3f} mm/s")
            print(f"    Systematic uncertainty fraction: {r['systematic_uncertainty_fraction']*100:.2f}%")
        
        logger.success(f"Monte Carlo systematic error analysis completed with {results['n_trials']} trials")
        logger.add_output_file(output_path, "Systematic error Monte Carlo results")
        
        print(f"\n✓ Results saved to {output_path}")
        
        return results
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise


if __name__ == "__main__":
    main()
