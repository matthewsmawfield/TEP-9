#!/usr/bin/env python3
"""
Statistical utilities for TEP-EFA analysis.
Ported from TEP-LLR Research-Grade Hardening Suite.
"""

import numpy as np
from scipy import linalg, stats
from typing import Dict

def robust_regression(y: np.ndarray, X: np.ndarray, weights: np.ndarray = None, 
                      scale_errors_by_birge: bool = True) -> Dict:
    """
    Perform numerically stable weighted linear regression using QR decomposition.
    
    Formula:
    --------
    Solve for β: y = Xβ + ε
    Weighted version: W^{1/2}y = W^{1/2}Xβ + ε'
    Using QR: Rβ = Qᵀ(W^{1/2}y)
    
    Covariance:
    -----------
    Σ_β = (XᵀWX)⁻¹ σ² = R⁻¹(R⁻ᵀ) σ²
    where σ² is the unbiased variance estimate (Birge-scaled if enabled).
    
    Parameters:
    -----------
    y : np.ndarray, shape (n,)
        Dependent variable (residuals)
    X : np.ndarray, shape (n, k)
        Design matrix (including intercept column if needed)
    weights : np.ndarray, shape (n,), optional
        Observation weights (1/σ²). If None, uniform weighting is used.
    scale_errors_by_birge : bool, default True
        If True, scales formal errors by the Birge Ratio max(1.0, sqrt(chi2_red)).
        
    Returns:
    --------
    Dict containing:
        - coefficients: np.ndarray (β)
        - errors: np.ndarray (scaled standard errors)
        - chi2_red: Reduced chi-squared
        - birge_ratio: Birge Ratio (sqrt of chi2_red)
        - condition_number: Matrix condition number κ(R)
        - n_obs: Number of observations
        - dof: Degrees of freedom
        - rss: Residual sum of squares
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    n, k = X.shape
    dof = n - k
    
    if n <= k:
        return {'coefficients': np.full(k, np.nan), 'errors': np.full(k, np.nan), 
                'chi2_red': np.nan, 'birge_ratio': np.nan, 'condition_number': np.nan}
        
    if weights is None:
        weights = np.ones(n)
    
    # Weighted transformation
    sqrt_w = np.sqrt(weights)
    yw = y * sqrt_w
    Xw = X * sqrt_w[:, np.newaxis]
    
    # QR Decomposition for stability
    # Xw = Q @ R
    Q, R = linalg.qr(Xw, mode='economic')
    
    # Check condition number
    # Stability threshold: kappa < 1e12 for float64
    s = linalg.svdvals(R)
    cond = s[0] / s[-1] if s[-1] > 0 else np.inf
    
    if cond > 1e12:
        return {
            'coefficients': np.full(k, np.nan), 
            'errors': np.full(k, np.nan), 
            'chi2_red': np.nan, 
            'birge_ratio': np.nan, 
            'condition_number': cond,
            'status': 'SINGULAR'
        }
    
    # Solve R β = Qᵀ yw
    qty = Q.T @ yw
    try:
        beta = linalg.solve_triangular(R, qty)
        
        # Residuals and Statistics
        y_pred = X @ beta
        residuals = y - y_pred
        rss = np.sum(weights * residuals**2)
        mse = rss / dof
        
        # Formal Covariance Matrix: (XᵀWX)⁻¹ = (RᵀR)⁻¹
        # (RᵀR)⁻¹ = R⁻¹ (R⁻ᵀ)
        R_inv = linalg.inv(R)
        cov_raw = (R_inv @ R_inv.T)
        errors_raw = np.sqrt(np.diag(cov_raw))
        
        # Birge Scaling (only scale up if chi2_red > 1)
        chi2_red = rss / dof
        birge_ratio = np.sqrt(chi2_red)
        scaling_factor = max(1.0, birge_ratio) if scale_errors_by_birge else 1.0
        
        errors = errors_raw * scaling_factor
        cov_scaled = cov_raw * (scaling_factor**2)
        
        return {
            'coefficients': beta,
            'errors': errors,
            'chi2_red': chi2_red,
            'birge_ratio': birge_ratio,
            'condition_number': cond,
            'n_obs': n,
            'dof': dof,
            'rss': rss,
            'mse': chi2_red,
            'cov': cov_scaled
        }
    except (linalg.LinAlgError, ValueError):
        return {'coefficients': np.full(k, np.nan), 'errors': np.full(k, np.nan), 
                'chi2_red': np.nan, 'birge_ratio': np.nan, 'condition_number': cond}

def detect_outliers_iqr(residuals: np.ndarray, multiplier: float = 3.0) -> np.ndarray:
    """Detect outliers using IQR method."""
    q25 = np.percentile(residuals, 25)
    q75 = np.percentile(residuals, 75)
    iqr = q75 - q25
    lower_bound = q25 - multiplier * iqr
    upper_bound = q75 + multiplier * iqr
    return (residuals < lower_bound) | (residuals > upper_bound)

def detect_outliers_sigma(residuals: np.ndarray, sigma_threshold: float = 5.0) -> np.ndarray:
    """Detect outliers using sigma (standard deviation) method."""
    median = np.median(residuals)
    mad = np.median(np.abs(residuals - median))
    sigma = 1.4826 * mad
    threshold = sigma_threshold * sigma
    return np.abs(residuals - median) > threshold
def weighted_mean(values: np.ndarray, uncertainties: np.ndarray, scale_errors_by_birge: bool = True) -> Dict:
    """
    Compute inverse-variance weighted mean with optional Birge scaling.
    
    Parameters:
    -----------
    values : np.ndarray
        Array of values to average
    uncertainties : np.ndarray
        Array of 1-sigma uncertainties
    scale_errors_by_birge : bool, default True
        If True, scales the final error by sqrt(chi2_red) if chi2_red > 1.
        
    Returns:
    --------
    Dict containing:
        - mean: Weighted mean
        - error: Standard error (possibly scaled)
        - chi2_red: Reduced chi-squared
        - birge_ratio: Birge Ratio
        - n_obs: Number of observations
    """
    values = np.asarray(values)
    uncertainties = np.asarray(uncertainties)
    weights = 1.0 / (uncertainties**2)
    
    n = len(values)
    if n == 0:
        return {'mean': np.nan, 'error': np.nan, 'chi2_red': np.nan}
    
    sum_w = np.sum(weights)
    mean = np.sum(values * weights) / sum_w
    
    # Internal error (formal uncertainty)
    error_int = 1.0 / np.sqrt(sum_w)
    
    # External error (from scatter)
    chi2 = np.sum(((values - mean) / uncertainties)**2)
    dof = n - 1
    chi2_red = chi2 / dof if dof > 0 else 0.0
    
    birge_ratio = np.sqrt(chi2_red)
    scaling_factor = max(1.0, birge_ratio) if scale_errors_by_birge else 1.0
    
    error = error_int * scaling_factor
    
    return {
        'mean': mean,
        'error': error,
        'chi2_red': chi2_red,
        'birge_ratio': birge_ratio,
        'n_obs': n,
        'dof': dof,
        'chi2': chi2
    }


def permutation_spearman_test(
    x: np.ndarray,
    y: np.ndarray,
    n_permutations: int = 100000,
    alternative: str = 'greater',
    random_seed: int = 42,
) -> Dict:
    """
    Permutation test for Spearman rank correlation.

    Tests the null hypothesis that the observed rank correlation between
    x and y could arise from random pairing. This is the rigorous
    non-parametric test for the geometry-correlation smoking gun:
    systematic errors have no mechanism to produce correlations between
    anomaly magnitude and trajectory asymmetry, so a significant
    permutation p-value directly constrains the systematic-error
    alternative.

    Parameters
    ----------
    x : np.ndarray
        Predictor values (e.g., trajectory asymmetry).
    y : np.ndarray
        Response values (e.g., observed anomaly magnitude).
    n_permutations : int
        Number of random permutations for the null distribution.
    alternative : str
        'greater' (one-sided: rho > 0), 'two-sided', or 'less'.
    random_seed : int
        Seed for reproducibility.

    Returns
    -------
    Dict with keys:
        - observed_rho: observed Spearman correlation
        - p_value: permutation p-value
        - null_median: median of null distribution
        - null_95ci: [2.5th, 97.5th] percentiles of null
        - n_permutations: number of permutations used
        - n_samples: sample size
        - alternative: test direction
    """
    rng = np.random.default_rng(random_seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)

    if n < 3:
        return {
            'observed_rho': np.nan,
            'p_value': np.nan,
            'null_median': np.nan,
            'null_95ci': [np.nan, np.nan],
            'n_permutations': 0,
            'n_samples': n,
            'alternative': alternative,
            'warning': 'Need at least 3 points for rank correlation',
        }

    observed_rho, _ = stats.spearmanr(x, y)

    null_dist = np.empty(n_permutations)
    for i in range(n_permutations):
        y_perm = rng.permutation(y)
        rho_perm, _ = stats.spearmanr(x, y_perm)
        null_dist[i] = rho_perm

    if alternative == 'greater':
        p_value = np.mean(null_dist >= observed_rho)
    elif alternative == 'less':
        p_value = np.mean(null_dist <= observed_rho)
    else:
        p_value = np.mean(np.abs(null_dist) >= np.abs(observed_rho))

    return {
        'observed_rho': float(observed_rho),
        'p_value': float(p_value),
        'null_median': float(np.median(null_dist)),
        'null_95ci': [float(np.percentile(null_dist, 2.5)), float(np.percentile(null_dist, 97.5))],
        'n_permutations': n_permutations,
        'n_samples': n,
        'alternative': alternative,
    }


def exact_spearman_pvalue(x: np.ndarray, y: np.ndarray) -> Dict:
    """
    Exact Spearman p-value via enumeration of all permutations.

    Only feasible for n <= 10. For n <= 8 this is fast; for n=9 or 10
    it may take a few seconds. For larger n, use permutation_spearman_test.

    Returns
    -------
    Dict with observed_rho, exact_p_value (one-sided, rho >= observed),
    and n_permutations (= n!).
    """
    from itertools import permutations

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)

    if n > 10:
        return {
            'observed_rho': np.nan,
            'exact_p_value': np.nan,
            'n_permutations': 0,
            'n_samples': n,
            'warning': f'n={n} too large for exact enumeration; use permutation_spearman_test',
        }

    observed_rho, _ = stats.spearmanr(x, y)

    count_ge = 0
    total = 0
    for perm in permutations(y):
        rho_perm, _ = stats.spearmanr(x, np.array(perm))
        if rho_perm >= observed_rho:
            count_ge += 1
        total += 1

    return {
        'observed_rho': float(observed_rho),
        'exact_p_value': float(count_ge / total),
        'n_permutations': total,
        'n_samples': n,
    }


def juno_falsification_power_analysis(
    predicted_dv_mm_s: float,
    prediction_uncertainty_mm_s: float,
    measurement_precision_mm_s: float,
    alpha: float = 0.05,
    power_targets: tuple[float, ...] = (0.5, 0.8, 0.9, 0.95),
) -> Dict:
    """
    Statistical power analysis for the Juno falsification test.

    Computes the power to detect the TEP-predicted Juno velocity shift
    given the measurement precision and the between-flyby prediction
    uncertainty (random-effects scatter).

    The total uncertainty is sigma_total = sqrt(sigma_pred^2 + sigma_meas^2).
    Power is computed for a one-sided z-test of H0: mu = 0 vs H1: mu > 0
    at significance level alpha.

    Also reports the Minimum Detectable Effect (MDE) -- the smallest
    true effect detectable at the specified power targets.

    Parameters
    ----------
    predicted_dv_mm_s : float
        TEP predicted velocity shift for Juno (mm/s).
    prediction_uncertainty_mm_s : float
        Random-effects prediction uncertainty (between-flyby scatter, mm/s).
    measurement_precision_mm_s : float
        Single-measurement DSN precision (mm/s).
    alpha : float
        Significance level for the test.
    power_targets : tuple
        Power targets at which to report MDE.

    Returns
    -------
    Dict with:
        - power: power to detect predicted_dv_mm_s at given noise levels
        - sigma_total: combined uncertainty (mm/s)
        - effect_size_d: Cohen's d = predicted_dv / sigma_total
        - mde_mm_s: dict mapping power targets to minimum detectable effects
        - falsification_precision: measurement precision needed for
          the observed null to falsify TEP at 3-sigma (i.e. the
          prediction would need to be > 3*sigma away from zero)
    """
    sigma_total = float(np.sqrt(
        prediction_uncertainty_mm_s**2 + measurement_precision_mm_s**2
    ))

    if sigma_total <= 0 or predicted_dv_mm_s <= 0:
        return {
            'power': float('nan'),
            'sigma_total_mm_s': sigma_total,
            'effect_size_d': float('nan'),
            'mde_mm_s': {},
            'falsification_precision_mm_s': float('nan'),
            'warning': 'Invalid input: predicted_dv must be > 0 and sigma_total must be > 0',
        }

    effect_size_d = predicted_dv_mm_s / sigma_total

    # One-sided z-test power
    z_alpha = stats.norm.ppf(1 - alpha)
    z_beta = effect_size_d - z_alpha
    power = float(stats.norm.cdf(z_beta))

    # Minimum Detectable Effect at each power target
    mde = {}
    for p in power_targets:
        if p <= alpha:
            continue
        z_beta_mde = stats.norm.ppf(p)
        mde[p] = float((z_alpha + z_beta_mde) * sigma_total)

    # Falsification precision: what measurement precision would make
    # the predicted effect falsifiable at 3-sigma given the prediction uncertainty?
    # We need predicted_dv > 3 * sqrt(sigma_pred^2 + sigma_meas^2)
    # => sigma_meas < sqrt((predicted_dv/3)^2 - sigma_pred^2)
    threshold = predicted_dv_mm_s / 3.0
    if threshold > prediction_uncertainty_mm_s:
        falsification_precision = float(np.sqrt(threshold**2 - prediction_uncertainty_mm_s**2))
    else:
        falsification_precision = 0.0  # Already falsified even with zero measurement noise

    return {
        'predicted_dv_mm_s': float(predicted_dv_mm_s),
        'prediction_uncertainty_mm_s': float(prediction_uncertainty_mm_s),
        'measurement_precision_mm_s': float(measurement_precision_mm_s),
        'sigma_total_mm_s': sigma_total,
        'effect_size_d': float(effect_size_d),
        'alpha': alpha,
        'power': power,
        'mde_mm_s': mde,
        'falsification_precision_mm_s': falsification_precision,
        'interpretation': (
            f"At sigma_total={sigma_total:.3f} mm/s, the predicted "
            f"{predicted_dv_mm_s:.2f} mm/s effect has Cohen's d={effect_size_d:.2f} "
            f"and one-sided power={power:.3f}."
        ),
    }
