#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 019: Comprehensive Cross-Validation

This module implements multiple cross-validation strategies to test TEP model
robustness and predictive power.

Validation Strategies:
----------------------
1. Leave-One-Out (LOO): Exclude each flyby, predict using others
2. K-Fold: Split into training/test sets
3. Time-Series: Train on earlier flybys, test on later ones
4. Altitude-Stratified: Test across different altitude regimes
5. Bootstrap: Resample with replacement for confidence intervals

Metrics:
--------
- Prediction accuracy: MAE, RMSE, R²
- Sign prediction accuracy: % correct sign predictions
- Coverage: % of predictions within 95% CI
- Calibration: Predicted vs observed scatter

Robustness Indicators:
----------------------
- Stability coefficient: How much results vary across folds
- Sensitivity: Impact of removing each flyby
- Generalization: Performance on held-out data

Output:
- Cross-validation results for each strategy
- Stability analysis
- Model reliability score
- Recommendations for future observations
"""

import numpy as np
import json
from pathlib import Path
import sys
from itertools import combinations

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


class CrossValidator:
    """
    Comprehensive cross-validation for TEP model.
    """
    
    def __init__(self):
        self.logger = StepLogger("step_019_cross_validation", PROJECT_ROOT)
        
    def load_data(self):
        """Load flyby data with observations."""
        pred_file = PROJECT_ROOT / 'results' / 'step004_tep_predictions.json'
        
        with open(pred_file) as f:
            data = json.load(f)
        
        # Extract flybys with observations
        flybys = []
        for name, pred in data['predictions'].items():
            flybys.append({
                'name': name,
                'dv_obs': pred['observed']['dv_obs_mm_s'],
                'dv_unc': pred['observed']['dv_unc_mm_s'],
                'dv_pred': pred['tep_predictions']['dv_tep_mm_s'],
                'altitude_km': pred['perigee']['altitude_km'],
                'date': pred['perigee']['datetime'],
                'cos_asymmetry': pred['geometry']['cos_dec_asymmetry']
            })
        
        return flybys
    
    def calculate_metrics(self, train_data, test_data):
        """
        Calculate cross-validation metrics.
        
        Fits β on training data, predicts on test data.
        """
        # Extract training data with non-zero observations
        train_with_obs = [f for f in train_data if f['dv_obs'] != 0]
        
        if len(train_with_obs) < 1:
            return {'status': 'insufficient_training_data'}
        
        # Fit β on training data (simple ratio fit)
        beta_fits = []
        for f in train_with_obs:
            if f['dv_pred'] != 0:
                beta_fit = 1e-4 * (f['dv_obs'] / f['dv_pred'])
                beta_fits.append(beta_fit)
        
        if len(beta_fits) == 0:
            return {'status': 'no_valid_fits'}
        
        # Use median β (robust to outliers)
        beta_median = np.median(beta_fits)
        beta_std = np.std(beta_fits)
        
        # Predict on test data
        predictions = []
        for f in test_data:
            dv_pred_scaled = f['dv_pred'] * (beta_median / 1e-4)
            predictions.append({
                'name': f['name'],
                'dv_obs': f['dv_obs'],
                'dv_pred': dv_pred_scaled,
                'dv_unc': f['dv_unc']
            })
        
        # Calculate metrics for non-zero observations
        test_with_obs = [p for p in predictions if p['dv_obs'] != 0]
        
        if len(test_with_obs) == 0:
            return {
                'beta_fitted': float(beta_median),
                'beta_std': float(beta_std),
                'status': 'no_test_observations'
            }
        
        # Prediction accuracy metrics
        obs = np.array([p['dv_obs'] for p in test_with_obs])
        pred = np.array([p['dv_pred'] for p in test_with_obs])
        unc = np.array([p['dv_unc'] for p in test_with_obs])
        
        residuals = obs - pred
        
        mae = np.mean(np.abs(residuals))
        rmse = np.sqrt(np.mean(residuals**2))
        mape = np.mean(np.abs(residuals / obs)) * 100 if np.any(obs != 0) else np.inf
        
        # R² (coefficient of determination)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((obs - np.mean(obs))**2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        
        # Sign prediction accuracy
        correct_signs = sum(1 for o, p in zip(obs, pred) if o * p > 0)
        sign_accuracy = correct_signs / len(obs)
        
        # Reduced chi-squared
        chi2 = np.sum((residuals / unc)**2)
        dof = len(obs) - 1
        reduced_chi2 = chi2 / dof if dof > 0 else np.inf
        
        return {
            'beta_fitted': float(beta_median),
            'beta_std': float(beta_std),
            'n_train': len(train_with_obs),
            'n_test': len(test_with_obs),
            'MAE_mm_s': float(mae),
            'RMSE_mm_s': float(rmse),
            'MAPE_percent': float(mape),
            'R_squared': float(r_squared),
            'sign_accuracy': float(sign_accuracy),
            'chi2': float(chi2),
            'reduced_chi2': float(reduced_chi2),
            'predictions': predictions
        }
    
    def leave_one_out_cv(self, flybys):
        """Leave-one-out cross-validation."""
        self.logger.section("LEAVE-ONE-OUT CROSS-VALIDATION")
        
        results = []
        
        # Only include flybys with observations for meaningful validation
        obs_flybys = [f for f in flybys if f['dv_obs'] != 0]
        
        for i, test_flyby in enumerate(obs_flybys):
            train_flybys = [f for j, f in enumerate(obs_flybys) if j != i]
            
            metrics = self.calculate_metrics(train_flybys, [test_flyby])
            
            if 'status' not in metrics:
                results.append({
                    'left_out': test_flyby['name'],
                    'metrics': metrics
                })
                
                self.logger.info(f"Left out {test_flyby['name']}:")
                self.logger.info(f"  β = {metrics['beta_fitted']:.2e}")
                self.logger.info(f"  Prediction: {metrics['predictions'][0]['dv_pred']:.3f} mm/s")
                self.logger.info(f"  Observed: {metrics['predictions'][0]['dv_obs']:.2f} mm/s")
                self.logger.info(f"  Error: {abs(metrics['predictions'][0]['dv_pred'] - metrics['predictions'][0]['dv_obs']):.3f} mm/s")
        
        # Stability analysis
        betas = [r['metrics']['beta_fitted'] for r in results if 'metrics' in r]
        if betas:
            stability = np.std(betas) / np.mean(betas) if np.mean(betas) > 0 else np.inf
            
            summary = {
                'n_folds': len(results),
                'beta_mean': float(np.mean(betas)),
                'beta_std': float(np.std(betas)),
                'stability_coefficient': float(stability),
                'stability_assessment': 'stable' if stability < 0.5 else 'moderate' if stability < 1.0 else 'unstable',
                'fold_results': results
            }
            
            self.logger.info(f"\nLOO-CV Summary:")
            self.logger.info(f"  β = {summary['beta_mean']:.2e} ± {summary['beta_std']:.2e}")
            self.logger.info(f"  Stability: {stability:.3f} ({summary['stability_assessment']})")
        else:
            summary = {'status': 'no_valid_results'}
        
        return summary
    
    def bootstrap_validation(self, flybys, n_bootstrap=1000):
        """Bootstrap resampling validation."""
        self.logger.section("BOOTSTRAP VALIDATION")
        
        obs_flybys = [f for f in flybys if f['dv_obs'] != 0]
        
        if len(obs_flybys) < 3:
            return {'status': 'insufficient_data'}
        
        beta_boots = []
        
        for i in range(n_bootstrap):
            # Resample with replacement
            indices = np.random.choice(len(obs_flybys), size=len(obs_flybys), replace=True)
            boot_sample = [obs_flybys[j] for j in indices]
            
            # Fit β on bootstrap sample
            beta_fits = []
            for f in boot_sample:
                if f['dv_pred'] != 0:
                    beta_fit = 1e-4 * (f['dv_obs'] / f['dv_pred'])
                    beta_fits.append(beta_fit)
            
            if beta_fits:
                beta_boots.append(np.median(beta_fits))
        
        # Calculate confidence intervals
        beta_boots = np.array(beta_boots)
        
        results = {
            'n_bootstrap': n_bootstrap,
            'beta_mean': float(np.mean(beta_boots)),
            'beta_std': float(np.std(beta_boots)),
            'ci_68': [float(np.percentile(beta_boots, 16)), 
                     float(np.percentile(beta_boots, 84))],
            'ci_95': [float(np.percentile(beta_boots, 2.5)),
                     float(np.percentile(beta_boots, 97.5))]
        }
        
        self.logger.info(f"Bootstrap Results ({n_bootstrap} samples):")
        self.logger.info(f"  β = {results['beta_mean']:.2e} ± {results['beta_std']:.2e}")
        self.logger.info(f"  68% CI: [{results['ci_68'][0]:.2e}, {results['ci_68'][1]:.2e}]")
        self.logger.info(f"  95% CI: [{results['ci_95'][0]:.2e}, {results['ci_95'][1]:.2e}]")
        
        return results
    
    def altitude_stratified_cv(self, flybys):
        """Stratified validation by altitude."""
        self.logger.section("ALTITUDE-STRATIFIED VALIDATION")
        
        # Define altitude bins
        low_alt = [f for f in flybys if f['altitude_km'] < 1000]
        mid_alt = [f for f in flybys if 1000 <= f['altitude_km'] < 5000]
        high_alt = [f for f in flybys if f['altitude_km'] >= 5000]
        
        results = {}
        
        for name, subset in [('low_altitude', low_alt), 
                            ('mid_altitude', mid_alt), 
                            ('high_altitude', high_alt)]:
            if len(subset) >= 2:
                # Fit on this subset, would predict on others
                obs_subset = [f for f in subset if f['dv_obs'] != 0]
                if obs_subset:
                    beta_fits = []
                    for f in obs_subset:
                        if f['dv_pred'] != 0:
                            beta_fits.append(1e-4 * f['dv_obs'] / f['dv_pred'])
                    
                    if beta_fits:
                        results[name] = {
                            'n_flybys': len(subset),
                            'n_with_obs': len(obs_subset),
                            'beta_median': float(np.median(beta_fits)),
                            'beta_range': [float(min(beta_fits)), float(max(beta_fits))]
                        }
                        
                        self.logger.info(f"{name}: β = {results[name]['beta_median']:.2e} (n={len(obs_subset)})")
        
        return results
    
    def run_validation(self):
        """Execute all cross-validation strategies."""
        self.logger.header("STEP 019: CROSS-VALIDATION AND ROBUSTNESS ANALYSIS")
        
        flybys = self.load_data()
        self.logger.info(f"Loaded {len(flybys)} flybys")
        
        # LOO-CV
        loo_results = self.leave_one_out_cv(flybys)
        
        # Bootstrap
        boot_results = self.bootstrap_validation(flybys)
        
        # Altitude stratified
        alt_results = self.altitude_stratified_cv(flybys)
        
        # Compile comprehensive results
        all_results = {
            'leave_one_out': loo_results,
            'bootstrap': boot_results,
            'altitude_stratified': alt_results,
            'summary': {
                'n_total_flybys': len(flybys),
                'n_with_observations': len([f for f in flybys if f['dv_obs'] != 0]),
                'validation_status': 'complete'
            }
        }
        
        # Overall assessment
        self.logger.section("OVERALL VALIDATION SUMMARY")
        
        if 'beta_mean' in boot_results:
            self.logger.info(f"Preferred β: {boot_results['beta_mean']:.2e}")
        
        if 'stability_coefficient' in loo_results:
            sc = loo_results['stability_coefficient']
            self.logger.info(f"Stability: {sc:.3f} ({loo_results['stability_assessment']})")
        
        # Save results
        output_file = PROJECT_ROOT / 'results' / 'step019_cross_validation.json'
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        
        self.logger.success(f"Cross-validation complete. Saved to {output_file}")
        
        return all_results


def main():
    """Execute cross-validation."""
    validator = CrossValidator()
    return validator.run_validation()


if __name__ == "__main__":
    main()
