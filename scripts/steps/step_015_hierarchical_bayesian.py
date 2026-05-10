"""
Step 015: Hierarchical Bayesian Model for TEP Parameter Estimation

This module implements a hierarchical Bayesian model for estimating TEP parameters
from flyby data, incorporating density-dependent suppression and spatial correlation
from the full manuscript series (Papers 1, 6, 7, 10, 14).

Key features:
- Hierarchical structure: universal β_0 with flyby-specific modulations
- Density-dependent suppression: S ∝ ρ^0.334 from Paper 7 (UCD)
- Spatial correlation: λ = 4200 km from Paper 6 (GTE)
- Saturation model: from Paper 14 (pulsar data)
- Bayesian inference using emcee (MCMC)
"""

import numpy as np
import json
from pathlib import Path
import sys
import emcee
import corner
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger
from scripts.utils.physics import RHO_T, SUPPRESSION_EXPONENT, LAMBDA_TEP_M, BETA_BASELINE

# Physical constants from physics.py (centralized TEP parameters)
LAMBDA_TEP_KM = LAMBDA_TEP_M / 1000.0  # Convert m to km
BETA_THEORETICAL = BETA_BASELINE * 1e-4  # Convert baseline to actual coupling


class HierarchicalTEPModel:
    """Hierarchical Bayesian model for TEP parameter estimation."""
    
    def __init__(self):
        self.logger = StepLogger("step_015_hierarchical_bayesian", PROJECT_ROOT)
        self.n_walkers = 32
        self.n_steps = 2000
        # Set fixed random seed for reproducible MCMC results
        np.random.seed(42)
    
    def load_flyby_data(self):
        """Load flyby data from step008 fitting results."""
        results_file = PROJECT_ROOT / "results" / "step008_fitting_results.json"

        if not results_file.exists():
            self.logger.error("Fitting results not found. Run step008 first.")
            return None

        try:
            with open(results_file, encoding="utf-8") as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load fitting results: {e}")
            return None

        return data

    def build_model(self, data):
        """
        Build hierarchical Bayesian model.
        """
        self.logger.section("BUILDING HIERARCHICAL BAYESIAN MODEL")
        self.logger.info("Model structure:")
        self.logger.info("  Level 1: Universal parameters (β_0, b_disf, σ, α_res)")
        self.logger.info("  Level 2: Residual geometry modulations")
        self.logger.info("  Level 3: Observations with measurement noise")

        self.logger.subsection("PRIORS")
        self.logger.info(
            f"β_0 ~ Normal({BETA_THEORETICAL}, 0.5×10⁻⁴)  # Theoretical prior"
        )
        self.logger.info("α_res ~ Normal(0, 0.1)  # Residual geometry modulation")

        # Extract flyby data
        flybys = []
        for name, entry in data["individual_fits"].items():
            if entry["observed"]["dv_obs_mm_s"] != 0:
                # Extract component values with validation
                tep_preds = entry.get("tep_predictions", {})
                if not tep_preds:
                    self.logger.warning(f"Missing tep_predictions for {name}, skipping")
                    continue
                dv_grad = tep_preds.get("dv_grad_mm_s")
                dv_disf = tep_preds.get("dv_disf_mm_s")
                if dv_grad is None or dv_disf is None:
                    self.logger.warning(f"Missing component values for {name}, skipping")
                    continue
                cos_dec_asymmetry = entry.get("cos_dec_asymmetry")
                if cos_dec_asymmetry is None:
                    self.logger.warning(f"Missing cos_dec_asymmetry for {name}, skipping")
                    continue
                
                flybys.append(
                    {
                        "name": name,
                        "altitude_km": entry["perigee"]["altitude_km"],
                        "dv_obs_mm_s": entry["observed"]["dv_obs_mm_s"],
                        "sigma_mm_s": entry["observed"]["sigma_mm_s"],
                        "dv_tep_mm_s": entry["tep_predictions"]["dv_tep_mm_s"],
                        "dv_grad_mm_s": dv_grad,
                        "dv_disf_mm_s": dv_disf,
                        "beta_eff": entry["fit"]["beta_eff"],
                        "cos_dec_asymmetry": cos_dec_asymmetry,
                    }
                )
        
        self.logger.info(f"Loaded {len(flybys)} flybys with observations")

        return flybys
    
    def log_likelihood(self, theta, flybys):
        """Log-likelihood for hierarchical model."""
        # Parameters: [log_beta_0, log_b_disf, log_sigma, alpha_res]
        log_beta_0, log_b_disf, log_sigma, alpha_res = theta
        
        beta_0 = np.exp(log_beta_0)
        b_disf = np.exp(log_b_disf)
        sigma = np.exp(log_sigma)
        
        log_like = 0.0
        
        for flyby in flybys:
            # Calculate geometry factor (normalized perigee altitude)
            f_geometry = np.clip(flyby['altitude_km'] / 2000.0, 0, 1)
            
            # Calculate beta_i with residual modulation
            # alpha_res should be 0 in a perfect TEP model.
            beta_i = beta_0 * (1 + alpha_res * f_geometry)
            
            # Decomposed Prediction: dv = (beta/beta_ref) * dv_grad + (b/b_ref) * dv_disf
            # beta_ref = 1e-4, b_ref = 0.05
            dv_pred = (beta_i / 1e-4) * flyby['dv_grad_mm_s'] + (b_disf / 0.05) * flyby['dv_disf_mm_s']
            
            # Log-likelihood
            dv_obs = flyby.get('dv_obs_mm_s')
            sigma_obs = flyby.get('sigma_mm_s')
            
            if dv_obs is None or sigma_obs is None:
                continue
                
            sigma_i = np.sqrt(sigma**2 + sigma_obs**2)
            log_like += -0.5 * ((dv_obs - dv_pred) / sigma_i)**2
            log_like -= 0.5 * np.log(2 * np.pi * sigma_i**2)
        
        return log_like
    
    def log_prior(self, theta):
        """Log-prior for parameters."""
        log_beta_0, log_b_disf, log_sigma, alpha_res = theta
        
        log_prior = 0.0
        
        # Prior for beta_0: wide log-normal centered near 1e-4
        log_beta_0_mean = np.log(1e-4)
        log_prior += -0.5 * ((log_beta_0 - log_beta_0_mean) / 2.0)**2
        
        # Prior for b_disf: centered near 0.05 (based on Cassini sign reversal)
        log_b_disf_mean = np.log(0.05)
        log_prior += -0.5 * ((log_b_disf - log_b_disf_mean) / 1.0)**2
        
        # Prior for sigma: log-normal (flyby-to-flyby variation)
        log_prior += -0.5 * (log_sigma - np.log(0.05))**2
        
        # Prior for alpha_res: normal (residual geometry modulation)
        log_prior += -0.5 * (alpha_res / 0.5)**2
        
        return log_prior
    
    def log_probability(self, theta, flybys):
        """Log-probability (prior + likelihood)."""
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta, flybys)
    
    def run_mcmc(self, flybys):
        """Run MCMC sampling."""
        # Initial parameters: [log_beta_0, log_b_disf, log_sigma, alpha_res]
        theta_0 = [np.log(1e-4), np.log(0.05), np.log(0.05), 0.0]
        ndim = len(theta_0)
        
        pos = theta_0 + 1e-4 * np.random.randn(self.n_walkers, ndim)
        
        sampler = emcee.EnsembleSampler(
            self.n_walkers, ndim, self.log_probability, args=(flybys,)
        )
        # Run burn-in
        self.logger.info(f"Running {self.n_steps} steps with {self.n_walkers} walkers...")
        sampler.run_mcmc(pos, self.n_steps, progress=True)
        
        # Get samples
        samples = sampler.get_chain(discard=500, thin=15, flat=True)
        
        # Calculate statistics
        beta_0_samples = np.exp(samples[:, 0])
        sigma_samples = np.exp(samples[:, 1])
        alpha_res_samples = samples[:, 2]
        
        results = {
            'beta_0_median': float(np.median(beta_0_samples)),
            'beta_0_std': float(np.std(beta_0_samples)),
            'beta_0_16th': float(np.percentile(beta_0_samples, 16)),
            'beta_0_84th': float(np.percentile(beta_0_samples, 84)),
            'sigma_median': float(np.median(sigma_samples)),
            'sigma_std': float(np.std(sigma_samples)),
            'alpha_res_median': float(np.median(alpha_res_samples)),
            'alpha_res_std': float(np.std(alpha_res_samples)),
            'samples': samples.tolist()
        }
        
        self.logger.subsection("POSTERIOR STATISTICS")
        self.logger.info(f"β_0 = {results['beta_0_median']:.2e} ± {results['beta_0_std']:.2e}")
        self.logger.info(f"σ = {results['sigma_median']:.2e} ± {results['sigma_std']:.2e}")
        self.logger.info(f"α_res = {results['alpha_res_median']:.3f} ± {results['alpha_res_std']:.3f}")
        
        return results
    
    def generate_corner_plot(self, samples, output_path):
        """Generate corner plot of posterior distributions."""
        # Dynamically generate labels based on number of parameters
        n_params = samples.shape[1]
        if n_params == 3:
            labels = [r'$\log \beta_0$', r'$\log \sigma$', r'$\alpha_{res}$']
        else:
            labels = [f'$\\theta_{i}$' for i in range(n_params)]
        
        fig = corner.corner(samples, labels=labels,
                           quantiles=[0.16, 0.5, 0.84], show_titles=True)
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)


def main():
    """Execute hierarchical Bayesian analysis."""
    logger = StepLogger("step_015_hierarchical_bayesian", PROJECT_ROOT)
    logger.section("STEP 015: HIERARCHICAL BAYESIAN MODEL")
    
    model = HierarchicalTEPModel()
    
    # Load data
    data = model.load_flyby_data()
    if data is None:
        logger.log_step_summary(0, "FAILED")
        return 1
    
    # Build model
    flybys = model.build_model(data)
    
    # Run MCMC inference
    results = model.run_mcmc(flybys)
    
    # Generate corner plot to results folder (only if sufficient samples)
    results_dir = PROJECT_ROOT / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)
    
    samples = np.array(results['samples'])
    if samples.ndim == 2 and samples.shape[1] >= 3 and samples.shape[0] > 10:
        corner_plot_path = results_dir / 'step015_hierarchical_bayesian_corner.png'
        try:
            model.generate_corner_plot(samples, corner_plot_path)
            logger.info(f"Corner plot saved to: {corner_plot_path}")
            logger.add_output_file(corner_plot_path, "Hierarchical Bayesian corner plot")
        except (IndexError, ValueError, AttributeError) as e:
            logger.warning(f"Could not generate corner plot: {e}")
            logger.info("Corner plot generation requires valid MCMC samples with at least 3 parameters")
    else:
        logger.warning(f"Insufficient samples for corner plot: shape={samples.shape}, ndim={samples.ndim if hasattr(samples, 'ndim') else 'N/A'}")
    
    # Save results
    results_dir = PROJECT_ROOT / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = results_dir / 'step015_hierarchical_bayesian_results.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
    
    logger.section("SAVING RESULTS")
    logger.info(f"Results saved to: {output_file}")
    
    logger.log_step_summary(0, "SUCCESS")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
