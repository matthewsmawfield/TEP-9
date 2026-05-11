#!/usr/bin/env python3
"""
Step 026: Stable Model Comparison with Regularized Likelihood

Four-tier model comparison framework addressing reviewer feedback on
Bayes factor overstatement and parameter-count transparency.

Model Tiers:
------------
1. Null (M_0): 0 parameters, predicts Δv = 0 for all flybys.

2. Anderson Empirical (M_A): 2 parameters. Fitted amplitude on trajectory
   asymmetry (cos δ_in - cos δ_out) plus offset. This is the empirical
   latitude/asymmetry formula from Anderson et al. (2008) simplified to
   asymmetry-only because perigee latitude data are not catalogued.

3. TEP Restricted (M_T^res): 1 parameter. Universal coupling β is the sole
   fitted parameter. All other TEP quantities are pre-specified:
   - λ_TEP ≈ 4000 km (from GNSS atomic clock correlations, Step 016)
   - S_⊕ ≈ 0.35 (from UCD soliton first-principles, Step 010)
   - v_trans ≈ 16.8 km/s (from TEP field equations)
   - Geometry (altitude, velocity, declinations) from JPL Horizons ephemerides

4. TEP Flexible (M_T^flex): 3 parameters. β (conformal coupling),
   b_disf (disformal amplitude), and offset (residual modulation capturing
   plasma screening, OD absorption, and unmodeled systematics).

The restricted model is the scientifically important one because all
parameters except β are pre-specified from independent measurements or
first-principles derivations. The Bayes factor reported for TEP vs Null
refers strictly to the restricted model.

Fixes from prior version:
-------------------------
- Replace the overfitted 4-parameter independent-β_i "Empirical" model with
  the Anderson asymmetry formula.
- Rename "TEP" to "TEP restricted" and add "TEP flexible" as a separate tier.
- Explicitly document pre-specified vs fitted parameters.
- Use log-likelihood and AICc/BIC corrections for stability.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


def convert_to_native_types(obj):
    """Recursively convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: convert_to_native_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_native_types(item) for item in obj]
    elif isinstance(obj, (np.integer, np.floating)):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        return convert_to_native_types(obj.tolist())
    else:
        return obj


class StableModelComparison:
    """
    Stable model comparison with regularized likelihood.
    """
    
    def __init__(self, logger):
        self.logger = logger
    
    def load_data(self):
        """Load TEP predictions, observations, and geometry components."""
        pred_file = PROJECT_ROOT / 'results' / 'step007_tep_predictions.json'
        
        if not pred_file.exists():
            self.logger.error(f"TEP predictions not found: {pred_file}")
            return None
        
        with open(pred_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Extract flybys with S/N >= 2
        flybys = []
        for name, pred in data['predictions'].items():
            dv_obs = pred['observed']['dv_obs_mm_s']
            dv_unc = pred['observed'].get('sigma_mm_s')
            
            if dv_obs is None or dv_unc is None:
                continue
            
            snr = abs(dv_obs) / dv_unc if dv_unc > 0 else 0
            if snr < 2:
                continue
            
            geometry = pred.get('geometry', {})
            tep_pred = pred.get('tep_predictions', {})
            
            flybys.append({
                'name': name,
                'dv_obs': dv_obs,
                'dv_unc': dv_unc,
                'dv_pred_base': tep_pred.get('dv_tep_mm_s', 0.0),
                'dv_grad_base': tep_pred.get('dv_grad_mm_s', 0.0),
                'dv_disf_base': tep_pred.get('dv_disf_mm_s', 0.0),
                'cos_asymmetry': geometry.get('cos_dec_asymmetry', 0.0),
                'snr': snr
            })
        
        return flybys
    
    def _gauss_loglike(self, residual, sigma_total):
        """Single Gaussian log-likelihood term."""
        return -0.5 * ((residual / sigma_total)**2) - 0.5 * np.log(2 * np.pi * sigma_total**2)

    def log_likelihood_anderson(self, flybys, A, B, sigma_sys=0.0):
        """
        Anderson empirical model: Δv = A * cos_asymmetry + B.
        
        This captures the core empirical finding of Anderson et al. (2008)
        that the anomaly correlates with trajectory asymmetry. Perigee
        latitude is not included because it is not catalogued.
        """
        log_like = 0.0
        for fb in flybys:
            dv_pred = A * fb['cos_asymmetry'] + B
            sigma_total = np.sqrt(fb['dv_unc']**2 + sigma_sys**2)
            residual = fb['dv_obs'] - dv_pred
            log_like += self._gauss_loglike(residual, sigma_total)
        return log_like

    def fit_anderson(self, flybys, sigma_sys=0.0):
        """Least-squares fit of Anderson empirical model to data."""
        X = np.array([[fb['cos_asymmetry'], 1.0] for fb in flybys])
        y = np.array([fb['dv_obs'] for fb in flybys])
        # Weighted least squares with measurement uncertainty
        w = np.array([1.0 / (fb['dv_unc']**2 + sigma_sys**2) for fb in flybys])
        W = np.diag(w)
        beta_hat = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ y, rcond=None)[0]
        A_fit, B_fit = beta_hat[0], beta_hat[1]
        log_like = self.log_likelihood_anderson(flybys, A_fit, B_fit, sigma_sys)
        return A_fit, B_fit, log_like

    def log_likelihood_tep_restricted(self, flybys, beta, sigma_sys=0.0):
        """
        TEP Restricted: 1 fitted parameter (β).
        
        All other quantities are pre-specified:
        - λ_TEP ≈ 4000 km (GNSS Step 016)
        - S_⊕ ≈ 0.35 (UCD Step 010)
        - v_trans ≈ 16.8 km/s (TEP field equations)
        - Geometry from JPL Horizons
        """
        log_like = 0.0
        for fb in flybys:
            dv_pred = fb['dv_pred_base'] * (beta / 1e-4)
            sigma_total = np.sqrt(fb['dv_unc']**2 + sigma_sys**2)
            residual = fb['dv_obs'] - dv_pred
            log_like += self._gauss_loglike(residual, sigma_total)
        return log_like

    def fit_tep_restricted(self, flybys, sigma_sys=0.0):
        """Weighted least-squares fit of β for TEP restricted."""
        num = sum(fb['dv_obs'] * fb['dv_pred_base'] / (fb['dv_unc']**2 + sigma_sys**2) for fb in flybys)
        den = sum(fb['dv_pred_base']**2 / (fb['dv_unc']**2 + sigma_sys**2) for fb in flybys)
        beta_fit = 1e-4 * (num / den) if den > 0 else 1e-4
        log_like = self.log_likelihood_tep_restricted(flybys, beta_fit, sigma_sys)
        return beta_fit, log_like

    def log_likelihood_tep_flexible(self, flybys, beta, b_disf, offset, sigma_sys=0.0):
        """
        TEP Flexible: 3 fitted parameters (β, b_disf, offset).
        
        The prediction is:
            Δv = β * dv_grad_base + b_disf * dv_disf_base + offset
        
        b_disf allows the disformal amplitude to vary independently,
        and offset captures residual modulation (plasma, OD, etc.).
        """
        log_like = 0.0
        for fb in flybys:
            dv_pred = beta * fb['dv_grad_base'] + b_disf * fb['dv_disf_base'] + offset
            sigma_total = np.sqrt(fb['dv_unc']**2 + sigma_sys**2)
            residual = fb['dv_obs'] - dv_pred
            log_like += self._gauss_loglike(residual, sigma_total)
        return log_like

    def fit_tep_flexible(self, flybys, sigma_sys=0.0):
        """Weighted least-squares fit of β, b_disf, offset for TEP flexible."""
        X = np.array([[fb['dv_grad_base'], fb['dv_disf_base'], 1.0] for fb in flybys])
        y = np.array([fb['dv_obs'] for fb in flybys])
        w = np.array([1.0 / (fb['dv_unc']**2 + sigma_sys**2) for fb in flybys])
        W = np.diag(w)
        beta_hat = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ y, rcond=None)[0]
        beta_fit, b_disf_fit, offset_fit = beta_hat[0], beta_hat[1], beta_hat[2]
        log_like = self.log_likelihood_tep_flexible(flybys, beta_fit, b_disf_fit, offset_fit, sigma_sys)
        return beta_fit, b_disf_fit, offset_fit, log_like

    def log_likelihood_null(self, flybys, sigma_sys=0.0):
        """
        Null model: 0 parameters, predicts Δv = 0.
        """
        log_like = 0.0
        for fb in flybys:
            sigma_total = np.sqrt(fb['dv_unc']**2 + sigma_sys**2)
            residual = fb['dv_obs']
            log_like += self._gauss_loglike(residual, sigma_total)
        return log_like
    
    def compute_aic_bic(self, log_like, n_params, n_data):
        """
        Compute AIC, AICc, and BIC with proper corrections.
        
        AIC = -2*log(L) + 2*k
        AICc = AIC + 2*k*(k+1)/(n-k-1)  (small sample correction)
        BIC = -2*log(L) + k*ln(n)
        """
        aic = -2 * log_like + 2 * n_params
        
        # AICc correction for small samples
        if n_data - n_params - 1 > 0:
            aic_correction = 2 * n_params * (n_params + 1) / (n_data - n_params - 1)
            aicc = aic + aic_correction
        else:
            aicc = float('inf')  # Not enough data
        
        bic = -2 * log_like + n_params * np.log(n_data)
        
        return aic, aicc, bic
    
    def bayes_factor_approx(self, log_like_1, log_like_2, k1, k2, n):
        """
        Approximate Bayes factor using BIC approximation.
        
        BF ≈ exp((BIC2 - BIC1)/2)
        
        This is a large-sample approximation that's numerically stable.
        """
        bic1 = -2 * log_like_1 + k1 * np.log(n)
        bic2 = -2 * log_like_2 + k2 * np.log(n)
        
        delta_bic = bic2 - bic1
        log_bf = delta_bic / 2
        
        # Clip to avoid overflow
        log_bf = np.clip(log_bf, -50, 50)
        
        return np.exp(log_bf), delta_bic
    
    def stable_model_comparison(self, flybys):
        """
        Four-tier model comparison with regularized likelihood.
        """
        self.logger.section("FOUR-TIER MODEL COMPARISON")
        
        n_data = len(flybys)
        self.logger.info(f"Using {n_data} flybys with S/N >= 2")
        
        # Load corrected uncertainty to get systematic component
        uncertainty_file = PROJECT_ROOT / 'results' / 'step025_corrected_uncertainty.json'
        with open(uncertainty_file, 'r', encoding='utf-8') as f:
            uncertainty_data = json.load(f)
        
        heterogeneity_rel_unc = uncertainty_data['corrected_uncertainty_budget']['heterogeneity_relative_uncertainty']
        obs_velocities = [fb['dv_obs'] for fb in flybys]
        mean_obs_velocity = np.mean(np.abs(obs_velocities))
        sigma_sys = mean_obs_velocity * heterogeneity_rel_unc
        self.logger.info(f"Systematic uncertainty (heterogeneity): {sigma_sys:.4f} mm/s")
        
        # --- Tier 1: Null ---
        log_like_null = self.log_likelihood_null(flybys, sigma_sys)
        aic_null, aicc_null, bic_null = self.compute_aic_bic(log_like_null, 0, n_data)
        
        # --- Tier 2: Anderson Empirical ---
        A_fit, B_fit, log_like_anderson = self.fit_anderson(flybys, sigma_sys)
        aic_anderson, aicc_anderson, bic_anderson = self.compute_aic_bic(log_like_anderson, 2, n_data)
        
        # --- Tier 3: TEP Restricted ---
        beta_fit, log_like_tep_res = self.fit_tep_restricted(flybys, sigma_sys)
        aic_tep_res, aicc_tep_res, bic_tep_res = self.compute_aic_bic(log_like_tep_res, 1, n_data)
        
        # --- Tier 4: TEP Flexible ---
        beta_flex, b_disf_flex, offset_flex, log_like_tep_flex = self.fit_tep_flexible(flybys, sigma_sys)
        aic_tep_flex, aicc_tep_flex, bic_tep_flex = self.compute_aic_bic(log_like_tep_flex, 3, n_data)
        
        self.logger.info("")
        self.logger.info("Fitted parameters:")
        self.logger.info(f"  Anderson: A={A_fit:.3f}, B={B_fit:.3f}")
        self.logger.info(f"  TEP restricted: beta={beta_fit:.2e}")
        self.logger.info(f"  TEP flexible: beta={beta_flex:.2e}, b_disf={b_disf_flex:.2e}, offset={offset_flex:.3f}")
        
        self.logger.info("")
        self.logger.info("Log-likelihoods:")
        self.logger.info(f"  Null:        {log_like_null:.2f}")
        self.logger.info(f"  Anderson:    {log_like_anderson:.2f}")
        self.logger.info(f"  TEP res:     {log_like_tep_res:.2f}")
        self.logger.info(f"  TEP flex:    {log_like_tep_flex:.2f}")
        
        self.logger.info("")
        self.logger.info("Information Criteria:")
        self.logger.info(f"  Null:        AIC={aic_null:.1f}, BIC={bic_null:.1f}")
        self.logger.info(f"  Anderson:    AIC={aic_anderson:.1f}, BIC={bic_anderson:.1f}")
        self.logger.info(f"  TEP res:     AIC={aic_tep_res:.1f}, BIC={bic_tep_res:.1f}")
        self.logger.info(f"  TEP flex:    AIC={aic_tep_flex:.1f}, BIC={bic_tep_flex:.1f}")
        
        # Bayes factors (all relative to Null)
        bf_anderson_null, dbic_anderson_null = self.bayes_factor_approx(
            log_like_anderson, log_like_null, 2, 0, n_data)
        bf_tep_res_null, dbic_tep_res_null = self.bayes_factor_approx(
            log_like_tep_res, log_like_null, 1, 0, n_data)
        bf_tep_flex_null, dbic_tep_flex_null = self.bayes_factor_approx(
            log_like_tep_flex, log_like_null, 3, 0, n_data)
        
        # Anderson vs TEP restricted
        bf_tep_res_anderson, dbic_tep_res_anderson = self.bayes_factor_approx(
            log_like_tep_res, log_like_anderson, 1, 2, n_data)
        
        self.logger.info("")
        self.logger.info("Bayes Factors (vs Null):")
        self.logger.info(f"  Anderson vs Null:       {bf_anderson_null:.2f} (ΔBIC={dbic_anderson_null:.1f})")
        self.logger.info(f"  TEP restricted vs Null: {bf_tep_res_null:.2f} (ΔBIC={dbic_tep_res_null:.1f})")
        self.logger.info(f"  TEP flexible vs Null:   {bf_tep_flex_null:.2f} (ΔBIC={dbic_tep_flex_null:.1f})")
        self.logger.info(f"  TEP res vs Anderson:    {bf_tep_res_anderson:.2f} (ΔBIC={dbic_tep_res_anderson:.1f})")
        
        # Model selection based on BIC
        bic_values = {
            'Null': bic_null,
            'Anderson': bic_anderson,
            'TEP_restricted': bic_tep_res,
            'TEP_flexible': bic_tep_flex
        }
        best_model_bic = min(bic_values, key=bic_values.get)
        
        # Akaike weights (TEP restricted vs Null vs Anderson)
        aicc_comp = {'TEP_restricted': aicc_tep_res, 'Null': aicc_null, 'Anderson': aicc_anderson}
        finite_aicc = {k: v for k, v in aicc_comp.items() if np.isfinite(v)}
        if finite_aicc:
            min_aicc = min(finite_aicc.values())
            delta_aicc = {k: v - min_aicc for k, v in finite_aicc.items()}
            sum_exp = sum(np.exp(-0.5 * da) for da in delta_aicc.values())
            akaike_weights = {k: np.exp(-0.5 * da) / sum_exp for k, da in delta_aicc.items()}
        else:
            akaike_weights = {'TEP_restricted': 1.0/3, 'Null': 1.0/3, 'Anderson': 1.0/3}
        
        results = {
            'n_data': n_data,
            'systematic_uncertainty_mm_s': float(sigma_sys),
            'log_likelihoods': {
                'Null': float(log_like_null),
                'Anderson': float(log_like_anderson),
                'TEP_restricted': float(log_like_tep_res),
                'TEP_flexible': float(log_like_tep_flex)
            },
            'fitted_parameters': {
                'Anderson': {'A': float(A_fit), 'B': float(B_fit)},
                'TEP_restricted': {'beta': float(beta_fit)},
                'TEP_flexible': {
                    'beta': float(beta_flex),
                    'b_disf': float(b_disf_flex),
                    'offset': float(offset_flex)
                }
            },
            'information_criteria': {
                'Null': {'AIC': float(aic_null), 'AICc': float(aicc_null), 'BIC': float(bic_null)},
                'Anderson': {'AIC': float(aic_anderson), 'AICc': float(aicc_anderson), 'BIC': float(bic_anderson)},
                'TEP_restricted': {'AIC': float(aic_tep_res), 'AICc': float(aicc_tep_res), 'BIC': float(bic_tep_res)},
                'TEP_flexible': {'AIC': float(aic_tep_flex), 'AICc': float(aicc_tep_flex), 'BIC': float(bic_tep_flex)}
            },
            'bayes_factors': {
                'Anderson_vs_Null': float(bf_anderson_null),
                'TEP_restricted_vs_Null': float(bf_tep_res_null),
                'TEP_flexible_vs_Null': float(bf_tep_flex_null),
                'TEP_restricted_vs_Anderson': float(bf_tep_res_anderson),
                'delta_BIC_Anderson_vs_Null': float(dbic_anderson_null),
                'delta_BIC_TEP_restricted_vs_Null': float(dbic_tep_res_null),
                'delta_BIC_TEP_flexible_vs_Null': float(dbic_tep_flex_null),
                'delta_BIC_TEP_restricted_vs_Anderson': float(dbic_tep_res_anderson)
            },
            'model_selection': {
                'best_model_BIC': best_model_bic,
                'akaike_weights': akaike_weights,
                'interpretation': f"BIC selects {best_model_bic} model"
            },
            'pre_specified_parameters': {
                'lambda_TEP_km': 4000,
                'lambda_TEP_source': 'GNSS atomic clock correlations (Step 016)',
                'S_earth': 0.35,
                'S_earth_source': 'UCD soliton first-principles (Step 010)',
                'v_trans_km_s': 16.8,
                'v_trans_source': 'TEP field equations',
                'geometry_source': 'JPL Horizons ephemerides'
            }
        }
        
        return results


def main():
    """Execute stable model comparison."""
    logger = StepLogger("step_026_stable_model_comparison", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("STEP 026: FOUR-TIER MODEL COMPARISON")
    
    logger.section("FRAMEWORK")
    logger.info("Tiers: Null (0 param) | Anderson (2 param) | TEP restricted (1 param) | TEP flexible (3 param)")
    logger.info("TEP restricted is the primary evidence claim: all parameters except beta are pre-specified.")
    
    comparison = StableModelComparison(logger)
    
    # Load data
    flybys = comparison.load_data()
    if flybys is None:
        logger.error("Failed to load data")
        return 1
    
    logger.info(f"Loaded {len(flybys)} flybys with S/N >= 2")
    
    # Perform stable model comparison
    results = comparison.stable_model_comparison(flybys)
    
    # Save results
    output_file = PROJECT_ROOT / 'results' / 'step026_stable_model_comparison.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(convert_to_native_types(results), f, indent=2)
    
    logger.success(f"Stable model comparison complete. Saved to {output_file}")
    logger.add_output_file(output_file, "Stable model comparison results")
    
    duration = time.time() - start_time
    logger.log_step_summary(duration, "SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
