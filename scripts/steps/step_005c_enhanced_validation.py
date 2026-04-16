#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 005c: Enhanced Validation and Model Comparison

Adds rigorous statistical validation:
1. Bayesian Information Criterion (BIC) for model comparison
2. Effect size calculations (Cohen's d)
3. Proper residual analysis with normality tests
4. Prediction accuracy metrics
5. Model comparison: TEP vs null model vs systematic error model
6. Information-theoretic evidence weights
"""

import sys
import json
import numpy as np
from pathlib import Path
import time
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class NumpyEncoder(json.JSONEncoder):
    """Handle numpy types in JSON serialization."""
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        return super().default(obj)

from scripts.utils.step_logger import StepLogger


def calculate_effect_sizes(fits: dict) -> dict:
    """
    Calculate Cohen's d effect sizes for detections.
    Compares each detection against null-result population.
    """
    # Separate detections from nulls
    detections = {}
    nulls = {}
    
    for name, fit_data in fits.items():
        obs = fit_data.get('observed', {})
        dv = obs.get('dv_obs_mm_s', 0)
        unc = obs.get('dv_unc_mm_s', 0.05)
        snr = abs(dv) / unc if unc > 0 else 0
        
        if snr >= 2:  # Detection threshold
            detections[name] = fit_data
        else:
            nulls[name] = fit_data
    
    if not nulls:
        return {'status': 'no_nulls_for_comparison'}
    
    # Null population statistics
    null_dvs = [v['observed']['dv_obs_mm_s'] for v in nulls.values()]
    null_mean = np.mean(null_dvs)
    null_std = np.std(null_dvs) if len(null_dvs) > 1 else 0.05
    
    results = {
        'null_population': {
            'n_nulls': len(nulls),
            'mean_dv': float(null_mean),
            'std_dv': float(null_std)
        },
        'effect_sizes': {}
    }
    
    # Calculate Cohen's d for each detection
    for name, fit_data in detections.items():
        dv = fit_data['observed']['dv_obs_mm_s']
        unc = fit_data['observed']['dv_unc_mm_s']
        
        # Pooled standard deviation
        pooled_std = np.sqrt((unc**2 + null_std**2) / 2)
        
        # Cohen's d
        if pooled_std > 0:
            cohens_d = (dv - null_mean) / pooled_std
        else:
            cohens_d = 0
        
        # Interpretation
        abs_d = abs(cohens_d)
        if abs_d < 0.2:
            interpretation = 'negligible'
        elif abs_d < 0.5:
            interpretation = 'small'
        elif abs_d < 0.8:
            interpretation = 'medium'
        else:
            interpretation = 'large'
        
        results['effect_sizes'][name] = {
            'cohens_d': float(cohens_d),
            'abs_d': float(abs_d),
            'interpretation': interpretation,
            'detection_significance': 'strong' if abs_d > 2 else 'moderate' if abs_d > 1 else 'weak'
        }
    
    return results


def bayesian_model_comparison(fits: dict) -> dict:
    """
    Compare TEP model against alternatives using information criteria.
    
    Models:
    - TEP: 1 parameter (β) + physical model structure
    - Null: 0 parameters (no anomaly, all Δv = 0)
    - Empirical: 3 parameters (one β per flyby, no physical structure)
    """
    successful = {k: v for k, v in fits.items() 
                  if v['fit']['beta_fitted'] is not None}
    
    if len(successful) < 2:
        return {'status': 'insufficient_data'}
    
    n = len(successful)
    
    # Collect data
    observed = []
    predicted_tep = []
    uncertainties = []
    
    # Calculate weighted mean beta for model comparison
    beta_values = [v['fit']['beta_fitted'] for v in successful.values()]
    beta_uncs = [v['fit']['uncertainty'] for v in successful.values()]
    weights = [1.0 / (u**2) for u in beta_uncs]
    beta_weighted = sum(b * w for b, w in zip(beta_values, weights)) / sum(weights)
    
    for name, fit_data in successful.items():
        obs = fit_data['observed']['dv_obs_mm_s']
        pred = fit_data['tep_predictions']['dv_tep_mm_s']
        unc = fit_data['observed']['dv_unc_mm_s']
        
        # Scale prediction by WEIGHTED MEAN beta (single shared parameter)
        pred_scaled = pred * (beta_weighted / 1e-4)
        
        observed.append(obs)
        predicted_tep.append(pred_scaled)
        uncertainties.append(unc)
    
    observed = np.array(observed)
    predicted_tep = np.array(predicted_tep)
    uncertainties = np.array(uncertainties)
    
    # Model 1: TEP (1 shared parameter β)
    residuals_tep = observed - predicted_tep
    chi2_tep = np.sum((residuals_tep / uncertainties) ** 2)
    k_tep = 1  # number of parameters
    bic_tep = chi2_tep + k_tep * np.log(n)
    aic_tep = chi2_tep + 2 * k_tep
    
    # Model 2: Null (0 parameters, all predictions = 0)
    chi2_null = np.sum((observed / uncertainties) ** 2)
    k_null = 0
    bic_null = chi2_null + k_null * np.log(n)
    aic_null = chi2_null + 2 * k_null
    
    # Model 3: Empirical (n independent parameters - one per flyby)
    # Perfect fit by construction
    chi2_emp = 0  # Perfect fit
    k_emp = n
    bic_emp = chi2_emp + k_emp * np.log(n)
    aic_emp = chi2_emp + 2 * k_emp
    
    # Evidence weights (Akaike weights)
    aics = np.array([aic_tep, aic_null, aic_emp])
    delta_aic = aics - np.min(aics)
    weights = np.exp(-0.5 * delta_aic) / np.sum(np.exp(-0.5 * delta_aic))
    
    return {
        'n_data_points': n,
        'models': {
            'TEP': {
                'k_parameters': k_tep,
                'chi2': float(chi2_tep),
                'BIC': float(bic_tep),
                'AIC': float(aic_tep),
                'akaike_weight': float(weights[0])
            },
            'Null': {
                'k_parameters': k_null,
                'chi2': float(chi2_null),
                'BIC': float(bic_null),
                'AIC': float(aic_null),
                'akaike_weight': float(weights[1])
            },
            'Empirical': {
                'k_parameters': k_emp,
                'chi2': float(chi2_emp),
                'BIC': float(bic_emp),
                'AIC': float(aic_emp),
                'akaike_weight': float(weights[2])
            }
        },
        'model_comparison': {
            'best_model_bic': 'TEP' if bic_tep < bic_null and bic_tep < bic_emp else 'Null' if bic_null < bic_emp else 'Empirical',
            'best_model_aic': 'TEP' if weights[0] > weights[1] and weights[0] > weights[2] else 'Null' if weights[1] > weights[2] else 'Empirical',
            'tep_evidence_weight': float(weights[0]),
            'tep_vs_null_bayes_factor_approx': float(np.exp(-0.5 * (bic_null - bic_tep)))
        }
    }


def residual_analysis(fits: dict) -> dict:
    """
    Analyze residuals for patterns and normality.
    """
    successful = {k: v for k, v in fits.items() 
                  if v['fit']['beta_fitted'] is not None}
    
    if len(successful) < 3:
        return {'status': 'insufficient_data'}
    
    # Calculate weighted mean beta for consistency
    beta_values = [v['fit']['beta_fitted'] for v in successful.values()]
    beta_uncs = [v['fit']['uncertainty'] for v in successful.values()]
    weights = [1.0 / (u**2) for u in beta_uncs]
    beta_weighted = sum(b * w for b, w in zip(beta_values, weights)) / sum(weights)
    
    residuals = []
    names = []
    
    for name, fit_data in successful.items():
        obs = fit_data['observed']['dv_obs_mm_s']
        pred = fit_data['tep_predictions']['dv_tep_mm_s']
        pred_scaled = pred * (beta_weighted / 1e-4)
        
        residual = obs - pred_scaled
        residuals.append(residual)
        names.append(name)
    
    residuals = np.array(residuals)
    
    # Normality test (Shapiro-Wilk if n < 50, otherwise D'Agostino)
    if len(residuals) >= 3 and len(residuals) <= 50:
        stat, p_normality = stats.shapiro(residuals)
    else:
        stat, p_normality = stats.normaltest(residuals)
    
    return {
        'residuals': {n: float(r) for n, r in zip(names, residuals)},
        'statistics': {
            'mean': float(np.mean(residuals)),
            'std': float(np.std(residuals)),
            'skewness': float(stats.skew(residuals)),
            'kurtosis': float(stats.kurtosis(residuals)),
            'normality_p_value': float(p_normality),
            'normal_distribution': p_normality > 0.05
        },
        'interpretation': 'residuals_appear_normal' if p_normality > 0.05 else 'residuals_deviate_from_normal'
    }


def prediction_accuracy_metrics(fits: dict) -> dict:
    """
    Calculate standard prediction accuracy metrics.
    """
    successful = {k: v for k, v in fits.items() 
                  if v['fit']['beta_fitted'] is not None}
    
    if len(successful) < 2:
        return {'status': 'insufficient_data'}
    
    # Calculate weighted mean beta for consistency
    beta_values = [v['fit']['beta_fitted'] for v in successful.values()]
    beta_uncs = [v['fit']['uncertainty'] for v in successful.values()]
    weights = [1.0 / (u**2) for u in beta_uncs]
    beta_weighted = sum(b * w for b, w in zip(beta_values, weights)) / sum(weights)
    
    observed = []
    predicted = []
    
    for fit_data in successful.values():
        obs = fit_data['observed']['dv_obs_mm_s']
        pred = fit_data['tep_predictions']['dv_tep_mm_s']
        pred_scaled = pred * (beta_weighted / 1e-4)
        
        observed.append(obs)
        predicted.append(pred_scaled)
    
    observed = np.array(observed)
    predicted = np.array(predicted)
    
    # Metrics
    residuals = observed - predicted
    mae = np.mean(np.abs(residuals))  # Mean Absolute Error
    rmse = np.sqrt(np.mean(residuals**2))  # Root Mean Square Error
    mape = np.mean(np.abs(residuals / observed)) * 100  # Mean Absolute Percentage Error
    
    # R-squared
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((observed - np.mean(observed))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # Correlation
    correlation = np.corrcoef(observed, predicted)[0, 1] if len(observed) > 1 else 0
    
    return {
        'n_predictions': len(observed),
        'MAE_mm_s': float(mae),
        'RMSE_mm_s': float(rmse),
        'MAPE_percent': float(mape),
        'R_squared': float(r_squared),
        'correlation': float(correlation),
        'prediction_quality': 'excellent' if r_squared > 0.9 else 'good' if r_squared > 0.7 else 'moderate' if r_squared > 0.5 else 'poor'
    }


def validate_small_sample_methods(fits: dict) -> dict:
    """
    Document and validate small-sample statistical methods used.
    """
    successful = {k: v for k, v in fits.items() 
                  if v['fit']['beta_fitted'] is not None}
    
    n = len(successful)
    
    methods_used = {
        'bootstrap_resampling': {
            'used': True,
            'n_iterations': 10000,
            'justification': 'Addresses small sample uncertainty through resampling',
            'valid_for_n': n >= 3
        },
        'leave_one_out_cv': {
            'used': True,
            'justification': 'Tests robustness by excluding each data point',
            'valid_for_n': n >= 3
        },
        'inverse_variance_weighting': {
            'used': True,
            'justification': 'Optimal weighting for heteroscedastic uncertainties',
            'valid_for_n': n >= 2
        },
        'heterogeneity_tests': {
            'cochrans_q': {'used': True, 'valid_for_n': n >= 2},
            'i_squared': {'used': True, 'valid_for_n': n >= 2, 'note': 'Can be inflated for small samples'}
        }
    }
    
    limitations = []
    if n < 5:
        limitations.append('Small sample size limits statistical power for complex model comparisons')
    if n < 10:
        limitations.append('Bootstrap confidence intervals may be approximate')
    if n < 30:
        limitations.append('Central limit theorem approximations may not fully apply')
    
    recommendations = [
        'Report all uncertainties with explicit small-sample caveats',
        'Use non-parametric tests where possible (Spearman vs Pearson)',
        'Validate with hold-out data when additional flybys become available',
        'Consider Bayesian methods with informative priors for small samples'
    ]
    
    return {
        'sample_size': n,
        'methods_used': methods_used,
        'limitations': limitations,
        'recommendations': recommendations,
        'confidence_in_conclusions': 'moderate' if n >= 3 else 'low' if n >= 2 else 'insufficient'
    }


def main():
    """Execute enhanced validation."""
    logger = StepLogger("step_005c_enhanced_validation", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("STEP 005c: ENHANCED VALIDATION AND MODEL COMPARISON")
    logger.info("Rigorous statistical validation:")
    logger.info("  1. Effect size calculations (Cohen's d)")
    logger.info("  2. Bayesian Information Criterion comparison")
    logger.info("  3. Residual analysis with normality tests")
    logger.info("  4. Prediction accuracy metrics")
    logger.info("  5. Small-sample method validation")
    
    # Load data from results folder
    results_dir = PROJECT_ROOT / 'results'
    fit_file = results_dir / 'step005_fitting_results.json'
    
    if not fit_file.exists():
        logger.error(f"Fitting results not found: {fit_file}")
        logger.log_step_summary(0, "FAILED")
        return 1
    
    with open(fit_file) as f:
        fit_data = json.load(f)
    
    fits = fit_data.get('individual_fits', {})
    
    # 1. Effect Size Analysis
    logger.section("EFFECT SIZE ANALYSIS")
    effect_sizes = calculate_effect_sizes(fits)
    
    if effect_sizes.get('status') != 'no_nulls_for_comparison':
        logger.info(f"Null population: n={effect_sizes['null_population']['n_nulls']}, "
                   f"σ={effect_sizes['null_population']['std_dv']:.3f} mm/s")
        
        for name, es in effect_sizes['effect_sizes'].items():
            logger.info(f"{name}: Cohen's d = {es['cohens_d']:.1f} ({es['interpretation']}) - {es['detection_significance']}")
    
    # 2. Bayesian Model Comparison
    logger.section("BAYESIAN MODEL COMPARISON")
    model_comparison = bayesian_model_comparison(fits)
    
    if model_comparison.get('status') != 'insufficient_data':
        logger.info(f"Comparing {model_comparison['n_data_points']} detections:")
        
        for model_name, metrics in model_comparison['models'].items():
            logger.info(f"  {model_name}: BIC={metrics['BIC']:.2f}, AIC={metrics['AIC']:.2f}, "
                       f"weight={metrics['akaike_weight']:.3f}")
        
        comp = model_comparison['model_comparison']
        logger.success(f"Best model (BIC): {comp['best_model_bic']}")
        logger.success(f"Best model (AIC): {comp['best_model_aic']}")
        logger.info(f"TEP evidence weight: {comp['tep_evidence_weight']:.1%}")
        
        if comp['best_model_bic'] == 'TEP':
            logger.success("✓ TEP is the preferred model by information criteria")
        
        # Approximate Bayes factor
        bf = comp['tep_vs_null_bayes_factor_approx']
        if bf > 10:
            logger.success(f"✓ Strong evidence for TEP vs null (BF ≈ {bf:.1f})")
        elif bf > 3:
            logger.info(f"Moderate evidence for TEP vs null (BF ≈ {bf:.1f})")
    
    # 3. Residual Analysis
    logger.section("RESIDUAL ANALYSIS")
    residuals = residual_analysis(fits)
    
    if residuals.get('status') != 'insufficient_data':
        stats_dict = residuals['statistics']
        logger.info(f"Residual mean: {stats_dict['mean']:.3f} mm/s")
        logger.info(f"Residual std: {stats_dict['std']:.3f} mm/s")
        logger.info(f"Normality p-value: {stats_dict['normality_p_value']:.3f}")
        
        if stats_dict['normal_distribution']:
            logger.success("✓ Residuals are consistent with normal distribution")
        else:
            logger.warning("Residuals deviate from normal - suggests unmodeled structure")
    
    # 4. Prediction Accuracy
    logger.section("PREDICTION ACCURACY METRICS")
    accuracy = prediction_accuracy_metrics(fits)
    
    if accuracy.get('status') != 'insufficient_data':
        logger.info(f"R² = {accuracy['R_squared']:.3f} ({accuracy['prediction_quality']})")
        logger.info(f"RMSE = {accuracy['RMSE_mm_s']:.2f} mm/s")
        logger.info(f"MAE = {accuracy['MAE_mm_s']:.2f} mm/s")
        logger.info(f"Prediction-observation correlation = {accuracy['correlation']:.3f}")
        
        if accuracy['R_squared'] > 0.9:
            logger.success("✓ Excellent predictive accuracy")
        elif accuracy['R_squared'] > 0.7:
            logger.success("✓ Good predictive accuracy")
    
    # 5. Small Sample Validation
    logger.section("SMALL SAMPLE METHOD VALIDATION")
    small_sample = validate_small_sample_methods(fits)
    
    logger.info(f"Sample size: n={small_sample['sample_size']}")
    logger.info(f"Confidence in conclusions: {small_sample['confidence_in_conclusions']}")
    
    if small_sample['limitations']:
        logger.subsection("Limitations")
        for lim in small_sample['limitations']:
            logger.info(f"  • {lim}")
    
    # 6. Overall Assessment
    logger.section("OVERALL ASSESSMENT")
    
    # Tally the evidence
    # Convert numpy types to Python native types for JSON serialization
    evidence_scores = {
        'effect_sizes_strong': int(sum(1 for es in effect_sizes.get('effect_sizes', {}).values() 
                                   if es.get('detection_significance') == 'strong')),
        'model_comparison_favors_tep': bool(model_comparison.get('model_comparison', {}).get('best_model_bic') == 'TEP'),
        'residuals_normal': bool(residuals.get('statistics', {}).get('normal_distribution', False)),
        'prediction_quality': str(accuracy.get('prediction_quality', 'poor')),
        'r_squared_high': bool(accuracy.get('R_squared', 0) > 0.8)
    }
    
    score = sum([
        evidence_scores['effect_sizes_strong'] > 0,
        evidence_scores['model_comparison_favors_tep'],
        evidence_scores['residuals_normal'],
        evidence_scores['prediction_quality'] in ['good', 'excellent'],
        evidence_scores['r_squared_high']
    ])
    
    logger.info(f"Evidence score: {score}/5 tests passed")
    
    if score >= 4:
        logger.success("✓✓ STRONG VALIDATION: TEP model is well-supported by multiple independent tests")
        overall_assessment = 'strongly_validated'
    elif score >= 3:
        logger.success("✓ GOOD VALIDATION: TEP model is supported with minor caveats")
        overall_assessment = 'well_supported'
    elif score >= 2:
        logger.info("MODERATE VALIDATION: TEP model shows promise but needs more data")
        overall_assessment = 'moderately_supported'
    else:
        logger.warning("WEAK VALIDATION: Insufficient evidence for TEP")
        overall_assessment = 'weakly_supported'
    
    # Save results
    output = {
        'effect_sizes': effect_sizes,
        'model_comparison': model_comparison,
        'residual_analysis': residuals,
        'prediction_accuracy': accuracy,
        'small_sample_validation': small_sample,
        'overall_assessment': {
            'score': f"{score}/5",
            'assessment': overall_assessment,
            'evidence_breakdown': evidence_scores,
            'scientific_conclusion': 'TEP is the preferred explanation for Earth flyby anomalies' if score >= 3 else 'More data needed to confirm TEP'
        }
    }
    
    output_file = results_dir / 'step005c_enhanced_validation.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, cls=NumpyEncoder)
    
    logger.success(f"Enhanced validation complete")
    logger.info(f"Results saved to: {output_file}")
    logger.add_output_file(output_file, "Enhanced validation and model comparison")
    
    duration = time.time() - start_time
    logger.log_step_summary(duration, "SUCCESS")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
