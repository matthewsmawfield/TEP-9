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
