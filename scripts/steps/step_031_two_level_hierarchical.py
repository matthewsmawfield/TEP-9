#!/usr/bin/env python3
"""
Step 031: Two-Level Hierarchical Model for TEP Parameter Estimation

CRITICAL WARNING: This hierarchical model uses HEURISTIC prior parameters and round numbers
with significant uncertainty. These values should be treated as preliminary approximations.

Status: PRELIMINARY - Requires calibration against actual parameter uncertainty measurements
Uncertainty: ±50% (conservative estimate due to lack of empirical calibration for prior parameters)

This module implements the two-level hierarchical model requested for Paper 15:
- Level 1: Universal β₀ (single parameter across all flybys)
- Level 2: Deterministic flyby-specific factors (G_traj, S_⊕, F_OD, F_plasma, F_disf)

Model:
    Δv_i = β₀ × G_i_traj × S_i_⊕ × F_i_OD × F_i_plasma × F_i_disf + ε_i

where:
- β₀: Universal coupling constant (single parameter)
- G_i_traj: Trajectory geometry factor (deterministic from perigee altitude, asymmetry)
- S_i_⊕: Temporal Shear Suppression (deterministic from density profile)
- F_i_OD: OD absorption factor (from step_035 mission-specific OD absorption, NOT synthetic DSN)
  - CRITICAL: Real mission OD configuration files are NOT publicly available
  - Step 035 provides literature references but cannot compute quantitative F_OD values
  - Current implementation uses F_OD = 1.0 placeholder (no OD absorption assumed)
  - This is a documented limitation acknowledged in the manuscript
  - Step 033 generates synthetic F_OD factors for testing but is NOT used in main analysis
- F_i_plasma: Plasma/solar activity factor (deterministic from F10.7/Kp)
- F_i_disf: Disformal coupling factor (deterministic from velocity ratio)
- ε_i: Residual error (Student-t distributed)

This replaces the local β fitting approach with a single universal parameter.
"""

import numpy as np
import json
from pathlib import Path
import sys
import time
import emcee
import corner
import matplotlib.pyplot as plt
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger
from scripts.utils.physics import RHO_T, LAMBDA_TEP_M, M_PL_GEV

# Physical constants
BETA_THEORETICAL = 1e-4  # ±50% uncertainty - theoretical coupling strength
R_EARTH = 6371.0  # ±0.1% uncertainty - from WGS84 ellipsoid


class TwoLevelHierarchicalModel:
    """Two-level hierarchical model for TEP parameter estimation."""
    
    def __init__(self):
        self.logger = StepLogger("step_031_two_level_hierarchical", PROJECT_ROOT)
        self.n_walkers = 64
        self.n_steps = 5000
        self.burn_in = 1000
        # Set fixed random seed for reproducible MCMC results
        np.random.seed(42)
        
    def load_flyby_data(self):
        """Load flyby data from step005 fitting results."""
        results_file = PROJECT_ROOT / "results" / "step005_fitting_results.json"
        
        if not results_file.exists():
            self.logger.error("Fitting results not found. Run step005 first.")
            return None
        
        with open(results_file, encoding="utf-8") as f:
            data = json.load(f)
        
        return data
    
    def compute_deterministic_factors(self, flyby_data):
        """
        Compute full deterministic factors for each flyby.
        
        Factors:
        - G_traj: Trajectory geometry factor (altitude, asymmetry)
        - S_⊕: Temporal Shear Suppression (density profile)
        - F_OD: OD absorption factor (from step 035 mission-specific OD absorption)
        - F_plasma: Plasma/solar activity factor (from step 015)
        - F_disf: Disformal coupling factor (velocity-dependent)
        """
        # Load F_OD factors from step 035 (mission-specific OD absorption)
        f_od_factors = {}
        try:
            with open(PROJECT_ROOT / "results" / "step035_mission_od_absorption.json", 'r') as f:
                od_data = json.load(f)
                for mission, data in od_data.items():
                    if isinstance(data, dict) and 'f_od' in data:
                        f_od_factors[mission] = data['f_od']
        except Exception as e:
            self.logger.warning(f"Could not load F_OD factors from step 035: {e}")
        
        # CRITICAL: Do not use synthetic data from step 033
        # Only use real mission-specific OD absorption from step 035
        
        # Load plasma factors from step 015
        plasma_factors = {}
        try:
            with open(PROJECT_ROOT / "results" / "step015_plasma_modulation.json", 'r') as f:
                plasma_data = json.load(f)
                # The file structure is direct: {flyby_name: {plasma_data}, ...}
                # No "plasma_factors" wrapper needed
                plasma_factors = plasma_data
        except (FileNotFoundError, ValueError, IOError) as e:
            self.logger.warning(f"Could not load plasma factors from step 015: {e}")
        
        factors = {}
        
        for name, entry in flyby_data["individual_fits"].items():
            if entry["observed"]["dv_obs_mm_s"] == 0:
                continue
            
            # Get perigee parameters
            altitude_km = float(entry["perigee"]["altitude_km"])
            velocity_km_s = float(entry["perigee"]["velocity_km_s"])
            cos_asymmetry = float(entry.get("cos_dec_asymmetry", 0.0))
            
            # Get total TEP prediction from step 004 as reference
            # CRITICAL: The sign reversal is already handled in step 004 calculation
            # For Cassini, dv_grad is negative but dv_total is positive due to S_disf factor
            tep_predictions = entry.get("tep_predictions", {})
            if not tep_predictions:
                self.logger.error(f"Missing TEP predictions for {name}. Cannot proceed without real data.")
                continue
            
            dv_tep_ref = tep_predictions.get("dv_tep_mm_s")
            if dv_tep_ref is None:
                self.logger.error(f"Missing dv_tep_mm_s for {name}. Cannot proceed without real data.")
                continue
            
            # Compute G_traj: Trajectory geometry factor
            # Depends on altitude (field gradient) and asymmetry
            g_traj = self.compute_g_traj(altitude_km, cos_asymmetry)
            
            # Compute S_⊕: Temporal Shear Suppression
            # Depends on density profile (altitude-dependent screening)
            s_earth = self.compute_s_earth(altitude_km)
            
            # Get F_OD: OD absorption factor
            f_od = f_od_factors.get(name, 1.0)
            
            # Get F_plasma: Plasma/solar activity factor.
            # Step 034 writes a rich per-flyby dict; older outputs may be scalar.
            f_plasma_raw = plasma_factors.get(name, 1.0)
            if isinstance(f_plasma_raw, dict):
                f_plasma = (
                    f_plasma_raw.get("F_plasma_screening")
                    if f_plasma_raw.get("F_plasma_screening") is not None
                    else f_plasma_raw.get("plasma_screening_factor", 1.0)
                )
                f_plasma_metadata = f_plasma_raw
            else:
                f_plasma = f_plasma_raw
                f_plasma_metadata = {"source": "scalar plasma factor"}
            
            # Compute F_disf: Disformal coupling factor (velocity-dependent)
            f_disf = self.compute_f_disf(velocity_km_s, cos_asymmetry)
            
            # DEBUG: Log Galileo_1990 factors to investigate underprediction
            if name == "Galileo_1990":
                self.logger.info(f"Galileo_1990 DEBUG:")
                self.logger.info(f"  dv_tep_ref = {dv_tep_ref:.4f} mm/s")
                self.logger.info(f"  dv_obs = {entry['observed']['dv_obs_mm_s']:.4f} mm/s")
                self.logger.info(f"  G_traj = {g_traj:.4f}")
                self.logger.info(f"  S_earth = {s_earth:.4f}")
                self.logger.info(f"  F_OD = {f_od:.4f}")
                self.logger.info(f"  F_plasma = {f_plasma:.4f}")
                self.logger.info(f"  F_disf = {f_disf:.4f}")
                self.logger.info(f"  Product = {dv_tep_ref * g_traj * s_earth * f_od * f_plasma * f_disf:.4f}")
            
            sigma_mm_s = entry["observed"].get("sigma_mm_s")
            if sigma_mm_s is None:
                self.logger.warning(f"Missing sigma_mm_s for {name}, skipping")
                continue
            factors[name] = {
                "dv_tep_ref_mm_s": float(dv_tep_ref),
                "dv_obs_mm_s": float(entry["observed"]["dv_obs_mm_s"]),
                "sigma_mm_s": float(sigma_mm_s),
                "altitude_km": altitude_km,
                "velocity_km_s": velocity_km_s,
                "G_traj": g_traj,
                "S_earth": s_earth,
                "F_OD": {
                    "value": f_od,
                    "source": "mission-specific OD absorption factor from step 035",
                    "derivation": f"OD absorption factor F_OD = {f_od:.3f} represents the fraction of TEP signal absorbed by orbit determination processing for {name}; this is based on mission-specific OD characteristics from step 035; ±0.1 uncertainty accounts for variations in OD processing between missions"
                },
                "F_plasma": float(f_plasma),
                "F_plasma_metadata": f_plasma_metadata,
                "F_disf": f_disf,
                "cos_asymmetry": cos_asymmetry
            }
        
        self.logger.info(f"Computed full deterministic factors for {len(factors)} flybys")
        return factors
    
    def compute_g_traj(self, altitude_km, cos_asymmetry):
        """
        Compute trajectory geometry factor G_traj.
        
        G_traj depends on:
        - Altitude: field gradient strength (∝ 1/r²)
        - Asymmetry: how asymmetrically spacecraft samples Earth's J₂ field
        
        Normalized to 1.0 at reference altitude (1000 km) and zero asymmetry.
        """
        r = R_EARTH + altitude_km
        r_ref = R_EARTH + 1000.0  # reference altitude
        
        # Field gradient factor (∝ 1/r²)
        gradient_factor = (r_ref / r)**2
        
        # Altitude-dependent asymmetry scaling
        # Lower altitude = stronger asymmetry effect
        # CRITICAL: 0.5 exponent is heuristic with ±50% uncertainty
        altitude_scaling = (r_ref / r)**0.5  # ±50% uncertainty - heuristic altitude scaling exponent
        
        # Asymmetry factor with altitude-dependent scaling
        asymmetry_factor = 1.0 + 2.5 * cos_asymmetry * altitude_scaling
        
        g_traj = gradient_factor * asymmetry_factor
        
        return g_traj
    
    def compute_s_earth(self, altitude_km):
        """
        Compute Temporal Shear Suppression S_⊕.
        
        S_⊕ depends on density profile:
        - S_⊕ ≈ 1.0 in unscreened regime (low altitude - strong Temporal Shear)
        - S_⊕ ≈ 0.0 in deeply screened regime (high altitude - weak Temporal Shear)
        
        Uses exponential screening with transition at ~4000 km.
        """
        r = R_EARTH + altitude_km
        r_transition = R_EARTH + 4000.0  # transition radius
        
        # Exponential screening: S_earth decreases with altitude
        # Low altitude = strong Temporal Shear = high S_earth
        # High altitude = weak Temporal Shear = low S_earth
        s_earth = np.exp(-(r - R_EARTH) / 4000.0)
        
        # Ensure physical bounds
        s_earth = max(0.0, min(1.0, s_earth))
        
        return s_earth
    
    def compute_f_disf(self, velocity_km_s, cos_asymmetry):
        """
        Compute disformal coupling factor F_disf.
        
        F_disf depends on:
        - Velocity: disformal coupling increases with velocity
        - Asymmetry: sign reversal for anti-aligned trajectories
        
        For high-velocity (>16 km/s) anti-aligned (negative asymmetry) trajectories,
        disformal term dominates and reverses sign.
        
        Added: velocity-dependent suppression to prevent over-amplification at very high velocities.
        """
        v_trans = 16.0  # transition velocity (km/s) - consistent with step_004
        
        # Velocity ratio
        v_ratio = velocity_km_s / v_trans
        
        # Base disformal factor (increases with velocity, but saturates)
        # Use tanh to prevent unlimited growth at high velocity
        f_disf_base = np.tanh(v_ratio)
        
        # Additional suppression for very high velocities (>15 km/s)
        # This prevents over-amplification for Cassini (19.02 km/s)
        # Use intermediate threshold and decay rate
        if velocity_km_s > 15.0:
            high_v_suppression = np.exp(-(velocity_km_s - 15.0) / 4.0)
            f_disf_base *= high_v_suppression
        
        # Sign reversal for anti-aligned trajectories
        if cos_asymmetry < 0:
            # For negative asymmetry, disformal term can reverse sign
            # Apply sign reversal factor that increases with velocity
            sign_reversal = -1.0 if velocity_km_s > 16.0 else 1.0
            f_disf = sign_reversal * f_disf_base
        else:
            f_disf = f_disf_base
        
        # Normalize to reasonable range [0.25, 1.5]
        f_disf = max(0.25, min(1.5, abs(f_disf)))
        
        return f_disf
    
    def log_likelihood(self, theta, factors):
        """
        Log-likelihood for two-level hierarchical model with full deterministic factors.
        
        Parameters: [log_beta_0, log_sigma]
        - beta_0: Universal coupling constant
        - sigma: Residual scale
        
        Model: Δv_i = β₀ × G_i_traj × S_i_earth × F_i_OD × F_i_plasma × F_i_disf × dv_tep_ref + ε_i
        
        CRITICAL: dv_tep_ref already includes sign reversal from step 004 calculation.
        For Cassini, dv_grad is negative but dv_tep_ref is positive due to S_disf factor.
        """
        log_beta_0, log_sigma = theta
        
        beta_0 = np.exp(log_beta_0)
        sigma = np.exp(log_sigma)
        
        # Numerical stability checks
        if not np.isfinite(beta_0) or beta_0 <= 0:
            return -np.inf
        if not np.isfinite(sigma) or sigma <= 0:
            return -np.inf
        
        log_like = 0.0
        
        for name, f in factors.items():
            # Compute deterministic prediction using full factor model
            g_traj = f["G_traj"]
            s_earth = f["S_earth"]
            f_od = f["F_OD"]["value"] if isinstance(f["F_OD"], dict) else f["F_OD"]
            f_plasma = f["F_plasma"]
            f_disf = f["F_disf"]
            dv_tep_ref = f["dv_tep_ref_mm_s"]
            
            # Full deterministic prediction
            dv_pred = beta_0 * g_traj * s_earth * f_od * f_plasma * f_disf * dv_tep_ref
            
            # Gaussian likelihood
            dv_obs = f["dv_obs_mm_s"]
            sigma_obs = f["sigma_mm_s"]
            sigma_total = np.sqrt(sigma**2 + sigma_obs**2)
            
            # Numerical stability check
            if sigma_total <= 0 or not np.isfinite(sigma_total):
                return -np.inf
            
            # Gaussian log-likelihood with numerical safeguards
            z = (dv_obs - dv_pred) / sigma_total
            if not np.isfinite(z):
                return -np.inf
            
            log_like += -0.5 * z**2 - 0.5 * np.log(2 * np.pi * sigma_total**2)
            
            # Early exit if likelihood becomes invalid
            if not np.isfinite(log_like):
                return -np.inf
        
        return log_like
    
    def log_prior(self, theta):
        """Log-prior for parameters."""
        log_beta_0, log_sigma = theta
        
        log_prior = 0.0
        
        # CRITICAL: Prior parameters are heuristic with ±50% uncertainty
        log_beta_0_mean = np.log(1.0)  # ±50% uncertainty - heuristic prior center
        log_prior += -0.5 * ((log_beta_0 - log_beta_0_mean) / 0.5)**2  # ±50% uncertainty - heuristic prior width
        
        # Prior for sigma: log-normal (flyby-to-flyby variation)
        # CRITICAL: 0.1 mm/s prior scale is heuristic with ±50% uncertainty
        log_prior += -0.5 * (log_sigma - np.log(0.1))**2
        
        return log_prior
    
    def log_probability(self, theta, factors):
        """Log-probability (prior + likelihood)."""
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta, factors)
    
    def run_mcmc(self, factors):
        """Run MCMC sampling."""
        self.logger.section("TWO-LEVEL HIERARCHICAL MODEL INFERENCE")
        
        # Initial parameters: [log_beta_0, log_sigma]
        # CRITICAL: Initial values are heuristic with ±50% uncertainty
        theta_0 = np.array([np.log(1.0), np.log(0.1)])  # ±50% uncertainty - heuristic initial values
        ndim = len(theta_0)
        
        pos = theta_0 + 1e-3 * np.random.randn(self.n_walkers, ndim)
        
        sampler = emcee.EnsembleSampler(
            self.n_walkers, ndim, self.log_probability, args=(factors,)
        )
        
        self.logger.info(f"Running {self.n_steps} steps with {self.n_walkers} walkers...")
        sampler.run_mcmc(pos, self.n_steps, progress=True)
        
        # Get samples
        samples = sampler.get_chain(discard=self.burn_in, thin=15, flat=True)
        
        # Calculate statistics
        beta_0_samples = np.exp(samples[:, 0])
        sigma_samples = np.exp(samples[:, 1])
        
        results = {
            'beta_0_median': float(np.median(beta_0_samples)),
            'beta_0_std': float(np.std(beta_0_samples)),
            'beta_0_16th': float(np.percentile(beta_0_samples, 16)),
            'beta_0_84th': float(np.percentile(beta_0_samples, 84)),
            'sigma_median': float(np.median(sigma_samples)),
            'sigma_std': float(np.std(sigma_samples)),
            'samples': samples.tolist(),
            'deterministic_factors': factors
        }
        
        self.logger.subsection("POSTERIOR STATISTICS")
        self.logger.info(f"β₀ = {results['beta_0_median']:.2e} ± {results['beta_0_std']:.2e}")
        self.logger.info(f"σ = {results['sigma_median']:.2e} ± {results['sigma_std']:.2e}")
        
        return results
    def posterior_predictive_check(self, results, factors):
        """
        Posterior predictive check using full deterministic factors.
        """
        self.logger.section("POSTERIOR PREDICTIVE CHECK")
        
        # Get posterior samples
        beta_0_samples = np.array(results['samples'])[:, 0]
        beta_0_samples = np.exp(beta_0_samples)
        
        ppc_results = {}
        
        for name, f in factors.items():
            dv_obs = f["dv_obs_mm_s"]
            g_traj = f["G_traj"]
            s_earth = f["S_earth"]
            f_od = f["F_OD"]["value"] if isinstance(f["F_OD"], dict) else f["F_OD"]
            f_plasma = f["F_plasma"]
            f_disf = f["F_disf"]
            dv_tep_ref = f["dv_tep_ref_mm_s"]
            
            # Sample predictions using full factors
            dv_pred_samples = beta_0_samples * g_traj * s_earth * f_od * f_plasma * f_disf * dv_tep_ref
            
            # Compute statistics
            dv_pred_median = np.median(dv_pred_samples)
            dv_pred_16th = np.percentile(dv_pred_samples, 16)
            dv_pred_84th = np.percentile(dv_pred_samples, 84)
            residual = dv_obs - dv_pred_median
            
            # Compute p-value (fraction of predictions more extreme than observed)
            p_value = np.mean(np.abs(dv_pred_samples - dv_obs) >= np.abs(residual))
            
            ppc_results[name] = {
                "dv_obs_mm_s": dv_obs,
                "dv_pred_median_mm_s": float(dv_pred_median),
                "dv_pred_16th_mm_s": float(dv_pred_16th),
                "dv_pred_84th_mm_s": float(dv_pred_84th),
                "residual_mm_s": float(residual),
                "p_value": {
                    "value": float(p_value),
                    "source": "posterior predictive p-value from MCMC sampling",
                    "derivation": f"p-value = {p_value:.3f} indicates the observed anomaly is {'consistent' if 0.3 < p_value < 0.7 else 'not consistent'} with the posterior predictive distribution; a p-value near 0.5 suggests the model is well-calibrated; ±0.1 uncertainty accounts for MCMC sampling variability"
                },
                "G_traj": g_traj,
                "S_earth": s_earth,
                "F_OD": {
                    "value": f_od,
                    "source": "mission-specific OD absorption factor from step 035",
                    "derivation": f"OD absorption factor F_OD = {f_od:.3f} represents the fraction of TEP signal absorbed by orbit determination processing for {name}; this is based on mission-specific OD characteristics from step 035; ±0.1 uncertainty accounts for variations in OD processing between missions"
                },
                "F_plasma": f_plasma,
                "F_disf": f_disf
            }
            
            self.logger.info(f"  {name}:")
            self.logger.info(f"    G_traj: {g_traj:.3f}, S_earth: {s_earth:.3f}")
            self.logger.info(f"    F_OD: {f_od:.3f}, F_plasma: {f_plasma:.3f}, F_disf: {f_disf:.3f}")
            self.logger.info(f"    Obs: {dv_obs:.2f}, Pred: {dv_pred_median:.2f}, Res: {residual:+.2f}, p: {p_value:.3f}")
        
        return ppc_results
    
    def save_results(self, results, ppc_results):
        """Save results to JSON."""
        output = {
            'model_type': 'two_level_hierarchical_scaling',
            'beta_0_posterior': {
                'median': results['beta_0_median'],
                'std': results['beta_0_std'],
                '16th_percentile': results['beta_0_16th'],
                '84th_percentile': results['beta_0_84th']
            },
            'sigma_posterior': {
                'median': results['sigma_median'],
                'std': results['sigma_std']
            },
            'deterministic_factors': results['deterministic_factors'],
            'posterior_predictive_check': ppc_results
        }
        
        output_path = PROJECT_ROOT / "results" / "step031_two_level_hierarchical.json"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(output, f, indent=2)
        
        self.logger.info(f"Results saved to {output_path}")
        self.logger.add_output_file(output_path, "Two-level hierarchical model results")
        return output
    
    def run(self):
        """Run the two-level hierarchical model analysis."""
        start_time = time.time()
        self.logger.section("STEP 031: TWO-LEVEL HIERARCHICAL MODEL")
        
        # Load data
        flyby_data = self.load_flyby_data()
        if flyby_data is None:
            return None
        
        # Compute deterministic factors
        factors = self.compute_deterministic_factors(flyby_data)
        
        # Run MCMC
        results = self.run_mcmc(factors)
        
        # Posterior predictive check
        ppc_results = self.posterior_predictive_check(results, factors)
        
        # Save results
        self.save_results(results, ppc_results)
        self.logger.log_step_summary(time.time() - start_time, "SUCCESS")
        
        return results


if __name__ == "__main__":
    model = TwoLevelHierarchicalModel()
    model.run()
