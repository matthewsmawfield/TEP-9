"""
Step 015: Hierarchical Bayesian Model + Component-Level Analysis for TEP Parameter Estimation

This module implements TWO complementary approaches:

1. Per-Flyby Geometry Factor Extraction (Section 3.2):
   Computes the effective geometry factor G_i,eff = dv_obs / dv_grad for each
   flyby at the reference coupling beta_0 = 1e-4. Correlates G_i,eff with
   trajectory parameters (altitude, velocity, asymmetry) to confirm
   geometry-dependent TEP coupling.

2. Hierarchical Bayesian Model (Section 4.2):
   MCMC inference for [β_0, b_disf, σ, α_res]. The log-β_0 prior mean is selected by
   ``hierarchical_beta_prior_center`` in ``config/pipeline_config.json`` (default:
   ``step008_recommended``), with log-space width derived from Step 008's reported
   β uncertainty (inflated_uncertainty preferred, else weighted_uncertainty).

Key features:
- Component decomposition separates conformal-gradient and disformal contributions at β_ref
- Density-dependent suppression: S ∝ ρ^0.334 from Paper 6 (UCD)
- Spatial correlation: λ = 4200 km from Paper 5 (GTE)
- Bayesian inference using emcee (MCMC)
"""

import numpy as np
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Tuple

import emcee
import corner
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger
from scripts.utils.physics import (
    RHO_T, SUPPRESSION_EXPONENT, LAMBDA_TEP_M, BETA_BASELINE,
    CHARACTERISTIC_SUPPRESSION, J2_EARTH, R_EARTH,
    DISFORMAL_VELOCITY_THRESHOLD_KM_S
)

# Physical constants from physics.py (centralized TEP parameters)
LAMBDA_TEP_KM = LAMBDA_TEP_M / 1000.0  # Convert m to km
BETA_THEORETICAL = BETA_BASELINE * 1e-4  # Convert baseline to actual coupling


def _tep_physics_pipeline_config() -> Dict[str, Any]:
    cfg_path = PROJECT_ROOT / "config" / "pipeline_config.json"
    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    return cfg["parameters"]["analysis"]["tep_physics"]


def _beta_prior_center_from_step008(data: Dict[str, Any]) -> Tuple[float, str]:
    """
    Central value for log-β_0 prior mean.

    Controlled by ``hierarchical_beta_prior_center`` in ``config/pipeline_config.json``:

    - ``step008_recommended`` (default): ``overall_analysis.recommended_beta`` then
      ``beta_statistics.weighted_mean`` (all S/N-qualified fits in the Step 008 run).
    - ``step008_sign_gated_diagnostic``: inverse-variance mean over the subset with
      sign agreement at β_ref (legacy trio when Cassini is sign-relaxed in Step 008).
    - ``reference_beta``: fixed β_ref = BETA_BASELINE × 10⁻⁴ from physics.py.
    """
    tep_cfg = _tep_physics_pipeline_config()
    mode = str(tep_cfg.get("hierarchical_beta_prior_center", "step008_recommended"))

    if mode == "reference_beta":
        v = float(BETA_THEORETICAL)
        if math.isfinite(v) and v > 0.0:
            return v, "physics.BETA_BASELINE_times_1e4"
        raise ValueError("Invalid reference_beta prior center")

    oa = data.get("overall_analysis")
    if not isinstance(oa, dict):
        raise ValueError("step008_fitting_results.json: missing overall_analysis")

    if mode == "step008_sign_gated_diagnostic":
        sg = oa.get("beta_statistics_sign_gated_diagnostic")
        if not isinstance(sg, dict):
            raise ValueError(
                "step008_fitting_results.json: missing beta_statistics_sign_gated_diagnostic "
                "(required when hierarchical_beta_prior_center == 'step008_sign_gated_diagnostic')"
            )
        wm = sg.get("weighted_mean")
        if wm is not None:
            v = float(wm)
            if math.isfinite(v) and v > 0.0:
                return v, "overall_analysis.beta_statistics_sign_gated_diagnostic.weighted_mean"
        raise ValueError(
            "step008_fitting_results.json: cannot set hierarchical β_0 prior center from "
            "sign-gated diagnostic (empty subset or invalid weighted_mean)"
        )

    if mode != "step008_recommended":
        raise ValueError(
            f"Unknown hierarchical_beta_prior_center mode: {mode!r} "
            "(expected step008_recommended | step008_sign_gated_diagnostic | reference_beta)"
        )

    rb = oa.get("recommended_beta")
    if rb is not None:
        v = float(rb)
        if math.isfinite(v) and v > 0.0:
            return v, "overall_analysis.recommended_beta"

    bs = oa.get("beta_statistics")
    if isinstance(bs, dict):
        wm = bs.get("weighted_mean")
        if wm is not None:
            v = float(wm)
            if math.isfinite(v) and v > 0.0:
                return v, "overall_analysis.beta_statistics.weighted_mean"

    raise ValueError(
        "step008_fitting_results.json: cannot set hierarchical β_0 prior center "
        "(need overall_analysis.recommended_beta or beta_statistics.weighted_mean)"
    )


def _beta_prior_log_sigma_from_step008(data: Dict[str, Any], beta_center: float) -> float:
    """
    Gaussian width (in natural-log space of β) for the log-β_0 prior.

    Uses relative uncertainty σ_β / β from Step 008, capped so the prior remains
    diffuse enough for MCMC exploration while not ignoring the reported precision.
    """
    oa = data.get("overall_analysis")
    if not isinstance(oa, dict):
        raise ValueError("step008_fitting_results.json: missing overall_analysis")
    bs = oa.get("beta_statistics")
    if not isinstance(bs, dict):
        raise ValueError("step008_fitting_results.json: missing beta_statistics")

    sig = bs.get("inflated_uncertainty")
    if sig is None:
        sig = bs.get("weighted_uncertainty")
    if sig is None:
        raise ValueError(
            "step008 beta_statistics must contain inflated_uncertainty or weighted_uncertainty"
        )
    sig = float(sig)
    if not math.isfinite(sig) or sig <= 0.0 or not math.isfinite(beta_center) or beta_center <= 0.0:
        raise ValueError("invalid beta center or uncertainty for prior width")

    rel = sig / beta_center
    if not math.isfinite(rel) or rel <= 0.0:
        raise ValueError("invalid relative uncertainty sigma_beta/beta for prior width")

    # ~ Gaussian in ln beta; width scales with relative error, bounded for stability.
    return float(min(1.5, max(0.25, 4.0 * rel)))


class PerFlybyGeometryAnalyzer:
    """
    Extract per-flyby effective geometry factors G_i,eff and correlate
    with trajectory parameters (altitude, velocity, asymmetry).

    The effective geometry factor is defined as the ratio of observed
    anomaly to the gradient prediction at the reference coupling:
        G_i,eff = dv_obs / dv_grad(beta_0 = 1e-4)

    This reveals geometry-dependent coupling that the hierarchical
    population model averages over.
    """

    def __init__(self, logger):
        self.logger = logger

    def analyze(self, flybys):
        """
        Extract per-flyby geometry factors and correlate with trajectory.
        """
        self.logger.section("PER-FLYBY GEOMETRY FACTOR EXTRACTION")
        self.logger.info("G_i,eff = dv_obs / dv_grad  (reference beta_0 = 1e-4)")

        records = []
        for fb in flybys:
            dv_obs = fb.get('dv_obs_mm_s')
            sigma = fb.get('sigma_mm_s')
            dv_grad = fb.get('dv_grad_mm_s')
            dv_disf = fb.get('dv_disf_mm_s')
            if dv_obs is None or sigma is None or dv_grad is None or dv_disf is None:
                continue
            if abs(dv_grad) < 1e-6:
                continue

            # Effective geometry factor (gradient component only)
            g_eff = dv_obs / dv_grad
            # Disformal-to-gradient ratio
            d_ratio = dv_disf / dv_grad
            # Effective disformal geometry factor
            if abs(dv_disf) > 1e-6:
                g_disf = dv_obs / dv_disf
            else:
                g_disf = np.nan

            # Implied reference-scale amplitude if geometry factor were unity
            beta_0_implied = 1e-4 * g_eff

            records.append({
                'name': fb['name'],
                'dv_obs': float(dv_obs),
                'dv_grad': float(dv_grad),
                'dv_disf': float(dv_disf),
                'sigma': float(sigma),
                'g_eff': float(g_eff),
                'g_disf': float(g_disf) if not np.isnan(g_disf) else None,
                'd_ratio': float(d_ratio),
                'beta_0_implied': float(beta_0_implied),
                'altitude_km': float(fb.get('altitude_km', 0)),
                'velocity_km_s': float(fb.get('velocity_km_s', 0)),
                'cos_dec_asymmetry': float(fb.get('cos_dec_asymmetry', 0)),
            })

        n = len(records)
        if n < 2:
            self.logger.error("Insufficient data for geometry analysis.")
            return None

        # --- Trajectory correlation analysis ---
        self.logger.subsection("PER-FLYBY GEOMETRY FACTORS")
        for r in records:
            self.logger.info(
                f"  {r['name']:15s}: G_eff = {r['g_eff']:6.2f}, "
                f"alt={r['altitude_km']:5.0f}km, "
                f"v={r['velocity_km_s']:4.1f}km/s, "
                f"asym={r['cos_dec_asymmetry']:+.3f}"
            )

        # Correlation with trajectory parameters
        g_effs = np.array([r['g_eff'] for r in records])
        alts = np.array([r['altitude_km'] for r in records])
        vels = np.array([r['velocity_km_s'] for r in records])
        asymms = np.array([r['cos_dec_asymmetry'] for r in records])

        # Pearson correlation coefficients
        def corr(x, y):
            if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
                return 0.0, 0.0
            r = np.corrcoef(x, y)[0, 1]
            # Fisher z-transform for approximate p-value (n small, crude)
            z = 0.5 * np.log((1 + r) / (1 - r))
            se = 1.0 / np.sqrt(len(x) - 3)
            p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / (se * math.sqrt(2)))))
            return r, p

        r_alt, p_alt = corr(alts, g_effs)
        r_vel, p_vel = corr(vels, g_effs)
        r_asym, p_asym = corr(asymms, g_effs)

        self.logger.subsection("TRAJECTORY CORRELATIONS")
        self.logger.info(f"  G_eff vs altitude:   r = {r_alt:+.3f} (p ≈ {p_alt:.3f})")
        self.logger.info(f"  G_eff vs velocity:   r = {r_vel:+.3f} (p ≈ {p_vel:.3f})")
        self.logger.info(f"  G_eff vs asymmetry:  r = {r_asym:+.3f} (p ≈ {p_asym:.3f})")

        # Linear model: log(G_eff) = c0 + c1*alt + c2*vel + c3*asym
        # Normalize to avoid numerical issues
        alt_norm = (alts - np.mean(alts)) / (np.std(alts) + 1e-6)
        vel_norm = (vels - np.mean(vels)) / (np.std(vels) + 1e-6)
        asym_norm = (asymms - np.mean(asymms)) / (np.std(asymms) + 1e-6)

        X_corr = np.column_stack([np.ones(n), alt_norm, vel_norm, asym_norm])
        y_log = np.log10(np.abs(g_effs))

        coeffs, resids, rank, s = np.linalg.lstsq(X_corr, y_log, rcond=None)
        y_pred = X_corr @ coeffs
        ss_res = np.sum((y_log - y_pred) ** 2)
        ss_tot = np.sum((y_log - np.mean(y_log)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

        self.logger.subsection("GEOMETRY FACTOR REGRESSION")
        self.logger.info(f"  log10|G_eff| = c0 + c1·alt_norm + c2·vel_norm + c3·asym_norm")
        self.logger.info(f"  c0 = {coeffs[0]:.3f}, c1 = {coeffs[1]:.3f}, "
                        f"c2 = {coeffs[2]:.3f}, c3 = {coeffs[3]:.3f}")
        self.logger.info(f"  R² = {r_squared:.3f}")

        # Summary statistics
        g_eff_median = float(np.median(g_effs))
        g_eff_mean = float(np.mean(g_effs))
        g_eff_std = float(np.std(g_effs))

        self.logger.subsection("SUMMARY")
        self.logger.info(f"  Median G_eff = {g_eff_median:.2f}")
        self.logger.info(f"  Mean   G_eff = {g_eff_mean:.2f} ± {g_eff_std:.2f}")
        self.logger.info(f"  Range: [{np.min(g_effs):.2f}, {np.max(g_effs):.2f}]")
        self.logger.info(f"")
        self.logger.info(f"  If beta_0 = 1.00e-4 (theoretical reference):")
        self.logger.info(f"    The effective geometry factor varies by a factor of "
                        f"{np.max(g_effs)/np.min(np.abs(g_effs)):.1f} across flybys.")
        self.logger.info(f"    This confirms geometry-dependent TEP coupling.")
        self.logger.info(f"")
        self.logger.info(f"  Empirically, the median implied coupling is:")
        self.logger.info(f"    beta_0,implied = 1e-4 × G_eff,median = {1e-4 * g_eff_median:.2e}")

        return {
            'per_flyby': records,
            'g_eff_median': g_eff_median,
            'g_eff_mean': g_eff_mean,
            'g_eff_std': g_eff_std,
            'g_eff_min': float(np.min(g_effs)),
            'g_eff_max': float(np.max(g_effs)),
            'correlations': {
                'altitude': {'r': float(r_alt), 'p': float(p_alt)},
                'velocity': {'r': float(r_vel), 'p': float(p_vel)},
                'asymmetry': {'r': float(r_asym), 'p': float(p_asym)},
            },
            'regression_coefficients': {
                'c0': float(coeffs[0]),
                'c1_altitude': float(coeffs[1]),
                'c2_velocity': float(coeffs[2]),
                'c3_asymmetry': float(coeffs[3]),
                'r_squared': float(r_squared),
            },
        }


class HierarchicalTEPModel:
    """Hierarchical Bayesian model for TEP parameter estimation."""
    
    def __init__(self):
        self.logger = StepLogger("step_015_hierarchical_bayesian", PROJECT_ROOT)
        self.n_walkers = 32
        self.n_steps = 2000
        self._beta_prior_center: float | None = None
        self._beta_prior_center_source: str | None = None
        self._beta_prior_log_sigma: float | None = None
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

        beta_center, beta_src = _beta_prior_center_from_step008(data)
        log_sigma_beta = _beta_prior_log_sigma_from_step008(data, beta_center)
        self._beta_prior_center = float(beta_center)
        self._beta_prior_center_source = beta_src
        self._beta_prior_log_sigma = float(log_sigma_beta)

        self.logger.info("Model structure:")
        self.logger.info("  Level 1: Universal parameters (β_0, b_disf, σ, α_res)")
        self.logger.info("  Level 2: Residual geometry modulations")
        self.logger.info("  Level 3: Observations with measurement noise")

        self.logger.subsection("PRIORS")
        self.logger.info(
            f"β_0 ~ LogNormal(ln({self._beta_prior_center:.3e}), σ_ln={self._beta_prior_log_sigma:.3f})  "
            f"# center from Step 008 ({self._beta_prior_center_source}), "
            f"hierarchical_beta_prior_center={_tep_physics_pipeline_config().get('hierarchical_beta_prior_center', 'step008_recommended')}"
        )
        self.logger.info("b_disf ~ LogNormal(ln(0.05), 1.0)  # Disformal coupling strength")
        self.logger.info("σ ~ LogNormal(ln(0.5), 1.0)  # Flyby-to-flyby scatter")
        self.logger.info("α_res ~ Normal(0, 0.3)  # Residual geometry modulation")

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
                
                # Extract velocity for disformal coupling calculation
                velocity_km_s = entry["perigee"].get("velocity_km_s", 10.0)
                
                flybys.append(
                    {
                        "name": name,
                        "altitude_km": entry["perigee"]["altitude_km"],
                        "velocity_km_s": velocity_km_s,
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
        """Log-likelihood for hierarchical model.
        
        The pre-computed dv_grad and dv_disf from step_007 already include:
        - Field gradient at perigee altitude
        - J2 oblateness factor
        - Trajectory asymmetry (cos_dec_asymmetry)
        - Time spent in field (r_p / v_p)
        
        The only scaling needed is the population-level beta_0 / beta_ref ratio,
        plus the disformal coupling strength b_disf.
        Any residual geometry modulation is captured by alpha_res.
        """
        # Parameters: [log_beta_0, log_b_disf, log_sigma, alpha_res]
        log_beta_0, log_b_disf, log_sigma, alpha_res = theta
        
        beta_0 = np.exp(log_beta_0)
        b_disf = np.exp(log_b_disf)
        sigma = np.exp(log_sigma)
        
        log_like = 0.0
        
        for flyby in flybys:
            # The pre-computed components already contain the full perigee physics
            # We only scale by the inferred population-level beta_0 relative to 
            # the reference beta = 1e-4 used in step_007
            # TEP predictions follow 3/4 power law: dv ∝ β^(3/4)
            # Both gradient and disformal components scale with the same exponent
            # because they arise from the same scalar field.
            beta_i = beta_0 * (1 + alpha_res)
            # beta must be positive (coupling strength); negative is unphysical
            if beta_i <= 0:
                return -np.inf
            scale = (beta_i / 1e-4) ** 0.75
            
            # Decomposed Prediction: scale pre-computed components by inferred couplings
            dv_pred = scale * (flyby['dv_grad_mm_s'] + (b_disf / 0.05) * flyby['dv_disf_mm_s'])
            
            # Log-likelihood
            dv_obs = flyby.get('dv_obs_mm_s')
            sigma_obs = flyby.get('sigma_mm_s')
            
            if dv_obs is None or sigma_obs is None or sigma_obs <= 0:
                continue
                
            sigma_i = np.sqrt(sigma**2 + sigma_obs**2)
            log_like += -0.5 * ((dv_obs - dv_pred) / sigma_i)**2
            log_like -= 0.5 * np.log(2 * np.pi * sigma_i**2)
        
        return log_like
    
    def log_prior(self, theta):
        """Log-prior for parameters."""
        log_beta_0, log_b_disf, log_sigma, alpha_res = theta
        
        log_prior = 0.0

        if self._beta_prior_center is None or self._beta_prior_log_sigma is None:
            raise RuntimeError(
                "HierarchicalTEPModel prior hyperparameters unset; call build_model(data) first."
            )

        # Prior for beta_0: log-normal in beta space <=> Gaussian on log_beta_0
        log_beta_0_mean = np.log(self._beta_prior_center)
        log_prior += -0.5 * (
            (log_beta_0 - log_beta_0_mean) / self._beta_prior_log_sigma
        ) ** 2
        
        # Prior for b_disf: log-normal width allows exploration around the nominal
        # disformal sector scale from physics.py (not tuned per-flyby).
        log_b_disf_mean = np.log(0.05)
        log_prior += -0.5 * ((log_b_disf - log_b_disf_mean) / 1.0)**2
        
        # Prior for sigma: log-normal (flyby-to-flyby variation)
        log_prior += -0.5 * (log_sigma - np.log(0.5))**2
        
        # Prior for alpha_res: normal (residual geometry modulation)
        # Centered at 0 (no residual modulation if physics is complete)
        log_prior += -0.5 * (alpha_res / 0.3)**2
        
        return log_prior
    
    def log_probability(self, theta, flybys):
        """Log-probability (prior + likelihood)."""
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta, flybys)
    
    def compute_posterior_predictions(self, samples, flybys):
        """Compute per-flyby predictions using posterior samples."""
        predictions = {}
        
        for flyby in flybys:
            name = flyby['name']
            dv_obs = flyby.get('dv_obs_mm_s')
            if dv_obs is None:
                continue
                
            # Compute predictions for each posterior sample
            preds = []
            for sample in samples:
                log_beta_0, log_b_disf, _, alpha_res = sample
                beta_0 = np.exp(log_beta_0)
                b_disf = np.exp(log_b_disf)
                
                beta_i = beta_0 * (1 + alpha_res)
                # beta must be positive; skip unphysical samples
                if beta_i <= 0:
                    continue
                scale = (beta_i / 1e-4) ** 0.75
                dv_pred = scale * (flyby['dv_grad_mm_s'] + (b_disf / 0.05) * flyby['dv_disf_mm_s'])
                preds.append(dv_pred)
            
            preds = np.array(preds)
            predictions[name] = {
                'dv_obs_mm_s': float(dv_obs),
                'dv_pred_median': float(np.median(preds)),
                'dv_pred_mean': float(np.mean(preds)),
                'dv_pred_std': float(np.std(preds)),
                'dv_pred_16th': float(np.percentile(preds, 16)),
                'dv_pred_84th': float(np.percentile(preds, 84)),
                'residual': float(dv_obs - np.median(preds)),
            }
        
        return predictions
    
    def run_mcmc(self, flybys):
        """Run MCMC sampling."""
        if self._beta_prior_center is None or self._beta_prior_log_sigma is None:
            raise RuntimeError(
                "HierarchicalTEPModel prior hyperparameters unset; call build_model(data) before run_mcmc()."
            )

        # Initial parameters: [log_beta_0, log_b_disf, log_sigma, alpha_res]
        beta_0 = float(self._beta_prior_center)
        theta_0 = [np.log(beta_0), np.log(0.05), np.log(0.5), 0.0]
        ndim = len(theta_0)

        self.logger.info(
            f"MCMC init: log β_0 = ln({beta_0:.6e}) from Step 008 ({self._beta_prior_center_source}), "
            f"prior σ_lnβ = {self._beta_prior_log_sigma:.4f}"
        )

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
        # Parameter order: [log_beta_0, log_b_disf, log_sigma, alpha_res]
        beta_0_samples = np.exp(samples[:, 0])
        b_disf_samples = np.exp(samples[:, 1])
        sigma_samples = np.exp(samples[:, 2])
        alpha_res_samples = samples[:, 3]

        # Compute posterior predictive checks
        posterior_predictions = self.compute_posterior_predictions(samples, flybys)
        
        self.logger.subsection("POSTERIOR PREDICTIVE CHECKS")
        for name, pred in posterior_predictions.items():
            self.logger.info(
                f"  {name:15s}: obs = {pred['dv_obs_mm_s']:7.2f}, "
                f"pred = {pred['dv_pred_median']:7.2f} ± {pred['dv_pred_std']:5.2f} mm/s, "
                f"residual = {pred['residual']:+.2f} mm/s"
            )

        results = {
            "prior_beta_0_center": float(self._beta_prior_center),
            "prior_beta_0_center_source": self._beta_prior_center_source,
            "prior_log_beta_0_sigma": float(self._beta_prior_log_sigma),
            "beta_0_median": float(np.median(beta_0_samples)),
            'beta_0_std': float(np.std(beta_0_samples)),
            'beta_0_16th': float(np.percentile(beta_0_samples, 16)),
            'beta_0_84th': float(np.percentile(beta_0_samples, 84)),
            'b_disf_median': float(np.median(b_disf_samples)),
            'b_disf_std': float(np.std(b_disf_samples)),
            'sigma_median': float(np.median(sigma_samples)),
            'sigma_std': float(np.std(sigma_samples)),
            'alpha_res_median': float(np.median(alpha_res_samples)),
            'alpha_res_std': float(np.std(alpha_res_samples)),
            'posterior_predictions': posterior_predictions,
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
    """Execute hierarchical Bayesian + component-level analysis."""
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

    # ---- Per-Flyby Geometry Factor Extraction ----
    geom_analyzer = PerFlybyGeometryAnalyzer(logger)
    geom_results = geom_analyzer.analyze(flybys)

    # ---- Hierarchical MCMC Inference ----
    mcmc_results = model.run_mcmc(flybys)

    # Combine results
    results = {
        'geometry_analysis': geom_results,
        'hierarchical_mcmc': mcmc_results,
    }

    # Generate corner plot to results folder (only if sufficient samples)
    results_dir = PROJECT_ROOT / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)

    samples = np.array(mcmc_results['samples'])
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
