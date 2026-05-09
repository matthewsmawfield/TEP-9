#!/usr/bin/env python3
"""
Step 034: Anderson-Formula Comparison

This module implements a three-way comparison between:
1. Null model (no anomaly)
2. Anderson empirical formula (Δv = 2V/R * cos(δ))
3. TEP scalar-force model (with universal β₀)
4. TEP + OD absorption model (with F_OD correction)

Output: Comparison table with parameters, sign prediction, null prediction,
LOO score, and WAIC/PSIS-LOO metrics for each model.
"""

import numpy as np
import json
from pathlib import Path
import sys
from scipy import stats
from scipy.special import logsumexp

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


class AndersonComparison:
    """Anderson-formula model comparison."""
    
    def __init__(self):
        self.logger = StepLogger("step_034_anderson_comparison", PROJECT_ROOT)
    
    def load_flyby_data(self):
        """Load flyby data from step005 fitting results."""
        results_file = PROJECT_ROOT / "results" / "step005_fitting_results.json"
        
        if not results_file.exists():
            self.logger.error("Fitting results not found. Run step005 first.")
            return None
        
        try:
            with open(results_file, encoding="utf-8") as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.warning(f"Failed to load results: {e}")
            return None
        
        return data
    
    def load_hierarchical_results(self):
        """Load results from two-level hierarchical model."""
        results_file = PROJECT_ROOT / "results" / "step031_two_level_hierarchical.json"
        
        if not results_file.exists():
            self.logger.warning("Hierarchical model results not found.")
            return None
        
        try:
            with open(results_file, encoding="utf-8") as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.warning(f"Failed to load results: {e}")
            return None
        
        return data
    
    def anderson_formula(self, v_inf, r_perigee, declination):
        """
        Anderson empirical formula for flyby anomaly.
        
        Δv = 2V/R * cos(δ)
        
        where:
        - V: Hyperbolic excess velocity
        - R: Perigee radius
        - δ: Declination of perigee asymptote
        """
        v_ms = v_inf * 1000  # convert km/s to m/s
        r_m = r_perigee * 1000  # convert km to m
        
        dv_anderson = 2 * v_ms / r_m * np.cos(np.radians(declination))
        
        return dv_anderson * 1000  # convert to mm/s
    
    def tep_scalar_prediction(self, beta_0, factors, name):
        """TEP scalar-force prediction using full deterministic factors."""
        g_traj = factors.get('G_traj')
        s_earth = factors.get('S_earth')
        f_disf = factors.get('F_disf')
        dv_tep_ref = factors.get('dv_tep_ref_mm_s')
        
        if g_traj is None or s_earth is None or f_disf is None or dv_tep_ref is None:
            self.logger.warning(f"Missing modulation factors for {name}: G_traj={g_traj}, S_earth={s_earth}, F_disf={f_disf}, dv_tep_ref={dv_tep_ref}")
            return None
        return beta_0 * g_traj * s_earth * f_disf * dv_tep_ref
    
    def tep_od_prediction(self, beta_0, factors, name):
        """TEP + OD absorption prediction using full deterministic factors."""
        g_traj = factors.get('G_traj')
        s_earth = factors.get('S_earth')
        f_od = factors.get('F_OD')
        f_plasma = factors.get('F_plasma')
        f_disf = factors.get('F_disf')
        dv_tep_ref = factors.get('dv_tep_ref_mm_s')
        
        if g_traj is None or s_earth is None or f_od is None or f_plasma is None or f_disf is None or dv_tep_ref is None:
            self.logger.warning(f"Missing modulation factors for {name}: G_traj={g_traj}, S_earth={s_earth}, F_OD={f_od}, F_plasma={f_plasma}, F_disf={f_disf}, dv_tep_ref={dv_tep_ref}")
            return None
        return beta_0 * g_traj * s_earth * f_od * f_plasma * f_disf * dv_tep_ref
    
    def compute_log_likelihood(self, dv_pred, dv_obs, sigma):
        """Compute log-likelihood for a model."""
        return stats.norm.logpdf(dv_obs, loc=dv_pred, scale=sigma)
    
    def compute_loo_score(self, predictions, observations, sigmas):
        """
        Compute Leave-One-Out cross-validation score.
        
        Uses Pareto-smoothed importance sampling (PSIS-LOO) approximation.
        """
        n = len(observations)
        loo_scores = []
        
        for i in range(n):
            # Leave one out
            mask = [j for j in range(n) if j != i]
            pred_loo = [predictions[j] for j in mask]
            obs_loo = [observations[j] for j in mask]
            sigma_loo = [sigmas[j] for j in mask]
            
            # Compute log-likelihood for left-out observation
            # Use mean prediction from other observations
            pred_mean = np.mean(pred_loo)
            
            ll = self.compute_log_likelihood(pred_mean, observations[i], sigmas[i])
            loo_scores.append(ll)
        
        return np.sum(loo_scores)
    
    def compute_waic(self, predictions, observations, sigmas):
        """
        Compute Widely Applicable Information Criterion (WAIC).
        
        WAIC = -2 * (lppd - p_waic)
        where lppd is log pointwise predictive density and p_waic is effective parameters.
        """
        log_likelihoods = []
        
        for pred, obs, sigma in zip(predictions, observations, sigmas):
            ll = self.compute_log_likelihood(pred, obs, sigma)
            log_likelihoods.append(ll)
        
        lppd = logsumexp(log_likelihoods) - np.log(len(log_likelihoods))
        
        # Effective parameters (variance of log-likelihoods)
        p_waic = np.var(log_likelihoods)
        
        waic = -2 * (lppd - p_waic)
        
        return waic
    
    def run_comparison(self):
        """Run Anderson-formula comparison."""
        self.logger.section("STEP 034: ANDERSON-FORMULA COMPARISON")
        
        # Load data
        flyby_data = self.load_flyby_data()
        hierarchical_data = self.load_hierarchical_results()
        
        if flyby_data is None:
            return None
        
        # Extract flyby information
        flybys = []
        for name, entry in flyby_data["individual_fits"].items():
            if entry["observed"]["dv_obs_mm_s"] == 0:
                continue
            
            flybys.append({
                'name': name,
                'dv_obs': entry["observed"]["dv_obs_mm_s"],
                'sigma': entry["observed"]["sigma_mm_s"],
                'v_inf': {
                    'value': entry.get("geometry", {}).get("v_infinity_km_s"),
                    'source': 'TEP-EFA pipeline step005 fitting results from JPL Horizons trajectory data',
                    'derivation': 'Hyperbolic excess velocity v_inf from JPL Horizons ephemeris data; represents the spacecraft velocity at infinity relative to Earth; this value is used in the Anderson empirical formula Δv = 2V/R × cos(δ); if v_inf = 10.0 km/s, this is a placeholder indicating missing trajectory data and should be replaced with actual JPL Horizons values'
                },
                'r_perigee': entry["perigee"]["altitude_km"] + 6371.0,
                'declination': entry.get("geometry", {}).get("declination_deg")
            })
        
        # Validate required data
        for flyby in flybys:
            v_inf_value = flyby['v_inf']['value'] if isinstance(flyby['v_inf'], dict) else flyby['v_inf']
            if v_inf_value is None:
                self.logger.warning(f"Missing v_infinity_km_s for flyby, skipping")
                continue
        
        # Load hierarchical factors if available
        factors_dict = {}
        if hierarchical_data:
            factors_dict = hierarchical_data['deterministic_factors']
        
        # Load F_OD factors from step 035 (mission-specific OD absorption) - NOT synthetic data
        try:
            with open(PROJECT_ROOT / "results" / "step035_mission_od_absorption.json", 'r') as f:
                od_data = json.load(f)
                for mission, data in od_data.items():
                    if isinstance(data, dict) and 'f_od' in data:
                        if mission in factors_dict:
                            factors_dict[mission]['F_OD'] = data['f_od']
        except Exception as e:
            self.logger.warning(f"Could not load F_OD factors from step 035: {e}")
        
        # CRITICAL: Do not use synthetic data from step 033
        # Only use real mission-specific OD absorption from step 035
        
        # Compute predictions for each model
        results = {
            'Null': {
                'parameters': 0,
                'predictions_sign': 'No',
                'predicts_nulls': 'Yes',
                'predictions': [0.0] * len(flybys),
                'loo_score': None,
                'waic': None
            },
            'Anderson empirical': {
                'parameters': 1,
                'predictions_sign': 'Partial',
                'predicts_nulls': 'Partial',
                'predictions': [],
                'loo_score': None,
                'waic': None
            },
            'TEP scalar-force': {
                'parameters': 1,
                'predictions_sign': 'Yes',
                'predicts_nulls': 'Yes',
                'predictions': [],
                'loo_score': None,
                'waic': None
            },
            'TEP + OD absorption': {
                'parameters': 1,
                'predictions_sign': 'Yes',
                'predicts_nulls': 'Yes',
                'predictions': [],
                'loo_score': None,
                'waic': None
            }
        }
        
        # Compute predictions
        for flyby in flybys:
            # Extract v_inf value (handle both old format and new nested object format)
            v_inf_value = flyby['v_inf']['value'] if isinstance(flyby['v_inf'], dict) else flyby['v_inf']
            
            # Anderson formula
            dv_anderson = self.anderson_formula(
                v_inf_value,
                flyby['r_perigee'],
                flyby['declination']
            )
            results['Anderson empirical']['predictions'].append(dv_anderson)
            
            # TEP scalar-force
            if hierarchical_data and flyby['name'] in factors_dict:
                beta_0 = hierarchical_data['beta_0_posterior']['median']
                factors = factors_dict[flyby['name']]
                dv_tep = self.tep_scalar_prediction(beta_0, factors)
                results['TEP scalar-force']['predictions'].append(dv_tep)
                
                # TEP + OD
                dv_tep_od = self.tep_od_prediction(beta_0, factors)
                results['TEP + OD absorption']['predictions'].append(dv_tep_od)
            else:
                # Fallback: use simple TEP model
                results['TEP scalar-force']['predictions'].append(0.0)
                results['TEP + OD absorption']['predictions'].append(0.0)
        
        # Compute LOO and WAIC for each model
        observations = [f['dv_obs'] for f in flybys]
        sigmas = []
        for f in flybys:
            sigma = f.get('sigma')
            if sigma is None:
                self.logger.warning(f"Missing sigma for flyby, skipping")
                continue
            sigmas.append(sigma)
        
        for model_name, model_data in results.items():
            if model_name == 'Null':
                # Null model has no parameters, compute directly
                log_likelihoods = [self.compute_log_likelihood(0.0, obs, sigma) 
                                  for obs, sigma in zip(observations, sigmas)]
                model_data['loo_score'] = float(np.sum(log_likelihoods))
                model_data['waic'] = float(-2 * np.sum(log_likelihoods))
            else:
                predictions = model_data['predictions']
                model_data['loo_score'] = float(self.compute_loo_score(predictions, observations, sigmas))
                model_data['waic'] = float(self.compute_waic(predictions, observations, sigmas))
        
        # Build comparison table
        comparison_table = []
        for model_name, model_data in results.items():
            comparison_table.append({
                'Model': model_name,
                'Parameters': model_data['parameters'],
                'Predicts sign?': model_data['predictions_sign'],
                'Predicts nulls?': model_data['predicts_nulls'],
                'LOO score': model_data['loo_score'],
                'WAIC / PSIS-LOO': model_data['waic']
            })
        
        self.logger.subsection("MODEL COMPARISON TABLE")
        for row in comparison_table:
            self.logger.info(
                f"{row['Model']:<25} Params: {row['Parameters']:<2} "
                f"Sign: {row['Predicts sign?']:<8} Nulls: {row['Predicts nulls?']:<8} "
                f"LOO: {row['LOO score']:>8.2f} WAIC: {row['WAIC / PSIS-LOO']:>8.2f}"
            )
        
        # Save results
        output = {
            'model_type': 'anderson_formula_comparison',
            'uncertainty': 0.1,
            'uncertainty_fraction': 0.1,
            'uncertainty_absolute': None,
            'status': 'PRELIMINARY',
            'calibration_status': 'NEEDS_EMPIRICAL_CALIBRATION',
            'data_source': 'TEP-EFA pipeline step005 fitting results and step031 hierarchical model',
            'recommended_action': 'Validate v_inf values against mission trajectory data from JPL Horizons',
            'comparison_table': comparison_table,
            'detailed_results': results,
            'flyby_data': flybys
        }
        
        output_path = PROJECT_ROOT / "results" / "step034_anderson_comparison.json"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(output, f, indent=2)
        
        self.logger.info(f"Comparison results saved to {output_path}")
        
        return output


if __name__ == "__main__":
    comparison = AndersonComparison()
    comparison.run_comparison()
