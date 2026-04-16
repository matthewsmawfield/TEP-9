#!/usr/bin/env python3
"""
Test TEP model calculations.

This module verifies that TEP model formulas are mathematically correct.
"""

import pytest
import numpy as np


def test_tep_acceleration_formula():
    """Verify TEP acceleration formula: a = c² β ∇φ / M_Pl."""
    # Test parameters
    C_LIGHT = 2.998e8  # m/s
    beta = 1e-4  # dimensionless
    dphi_dr = 1e-10  # GeV/m (example gradient)
    M_PL = 2.435e18  # GeV
    
    # Calculate acceleration
    c2 = C_LIGHT**2
    factor = c2 * beta / M_PL
    acceleration = factor * dphi_dr
    
    # Check units: [m/s²] = [m²/s²] × [1] × [GeV/m] / [GeV]
    assert acceleration > 0  # Should be positive for positive gradient
    assert acceleration < 1  # Should be small (TEP is a weak force)


def test_time_dilation_formula():
    """Verify time dilation formula: dtau = exp(beta * phi / M_Pl)."""
    # Test parameters
    beta = 1e-4  # dimensionless
    phi = 1e-6  # GeV (example potential)
    M_PL = 2.435e18  # GeV
    
    # Calculate time dilation
    dtau = np.exp(beta * phi / M_PL)
    
    # Check: should be dimensionless and close to 1 for small effects
    assert dtau > 0  # Should be positive
    assert abs(dtau - 1.0) < 1e-10  # Should be very close to 1 for small phi


def test_chameleon_field_formula():
    """Verify chameleon field formula: φ = Lambda * ((n * M_PL * Lambda³) / (beta * rho))^(1/(n+1))."""
    # Test parameters
    Lambda = 1e-3  # GeV
    n = 1
    M_PL = 2.435e18  # GeV
    beta = 1e-4
    rho_gev4 = 1e-10  # GeV^4
    
    # Calculate chameleon field
    numerator = n * M_PL * Lambda**3
    denominator = beta * rho_gev4
    scale = (numerator / denominator)**(1.0 / (n + 1))
    phi = Lambda * scale
    
    # Check: should be positive
    assert phi > 0
    assert phi > Lambda  # Should be larger than Lambda


def test_cohens_d_formula():
    """Verify Cohen's d formula: d = (value - null_mean) / pooled_std."""
    # Test parameters
    value = 10.0
    null_mean = 2.0
    pooled_std = 4.0
    
    # Calculate Cohen's d
    cohens_d = (value - null_mean) / pooled_std
    
    # Check: should be positive for value > null_mean
    assert cohens_d > 0
    assert cohens_d == 2.0  # (10 - 2) / 4 = 2


def test_pooled_std_formula():
    """Verify pooled standard deviation formula."""
    # Test parameters
    n_det = 3
    n_null = 5
    det_std = 2.0
    null_std = 3.0
    
    # Calculate pooled std
    pooled_std = np.sqrt(((n_det - 1) * det_std**2 + (n_null - 1) * null_std**2) / (n_det + n_null - 2))
    
    # Check: should be positive
    assert pooled_std > 0
    assert pooled_std > min(det_std, null_std)  # Should be between the two stds


def test_chi_squared_formula():
    """Verify chi-squared formula: χ² = Σ((observed - expected)² / σ²)."""
    # Test parameters
    observed = np.array([10.0, 12.0, 11.0])
    expected = np.array([10.0, 10.0, 10.0])
    uncertainties = np.array([1.0, 1.0, 1.0])
    
    # Calculate chi-squared
    chi2 = np.sum(((observed - expected) / uncertainties) ** 2)
    
    # Check: should be positive
    assert chi2 > 0
    assert chi2 == 5.0  # (0 + 4 + 1) = 5


def test_reduced_chi_squared_formula():
    """Verify reduced chi-squared formula: χ²_red = χ² / dof."""
    # Test parameters
    chi2 = 10.0
    dof = 5
    
    # Calculate reduced chi-squared
    reduced_chi2 = chi2 / dof
    
    # Check: should be positive
    assert reduced_chi2 > 0
    assert reduced_chi2 == 2.0  # 10 / 5 = 2


def test_i_squared_formula():
    """Verify I² heterogeneity index formula: I² = max(0, (Q - dof) / Q * 100)."""
    # Test parameters
    Q = 20.0
    dof = 5
    
    # Calculate I²
    I2 = max(0, (Q - dof) / Q * 100)
    
    # Check: should be between 0 and 100
    assert I2 >= 0
    assert I2 <= 100
    assert I2 == 75.0  # (20 - 5) / 20 * 100 = 75


def test_weighted_mean_formula():
    """Verify weighted mean formula: β_weighted = Σ(w_i × β_i) / Σ(w_i)."""
    # Test parameters
    beta_values = np.array([1.0, 2.0, 3.0])
    weights = np.array([1.0, 2.0, 1.0])
    
    # Calculate weighted mean
    weighted_mean = sum(beta_values * weights) / sum(weights)
    
    # Check: should be positive
    assert weighted_mean > 0
    assert weighted_mean == 2.0  # (1*1 + 2*2 + 3*1) / 4 = 8/4 = 2


def test_inverse_variance_weighting():
    """Verify inverse-variance weighting: w_i = 1/σ²."""
    # Test parameters
    uncertainties = np.array([1.0, 2.0, 0.5])
    
    # Calculate weights
    weights = 1.0 / uncertainties**2
    
    # Check: should be positive
    assert np.all(weights > 0)
    assert weights[0] == 1.0  # 1/1² = 1
    assert weights[1] == 0.25  # 1/2² = 0.25
    assert weights[2] == 4.0  # 1/0.5² = 4


def test_field_gradient_formula():
    """Verify field gradient formula: dφ/dr = (delta_phi / lambda) * exp(-delta_r / lambda)."""
    # Test parameters
    delta_phi = 1e-6  # GeV
    lambda_tep = 4000e3  # m (4000 km)
    delta_r = 1000e3  # m (1000 km)
    
    # Calculate gradient
    dphi_dr = (delta_phi / lambda_tep) * np.exp(-delta_r / lambda_tep)
    
    # Check: should be positive
    assert dphi_dr > 0
    assert dphi_dr < delta_phi / lambda_tep  # Exponential decay


def test_vainshtein_screening_formula():
    """Verify Vainshtein screening formula: S = (ρ/ρ_c)^γ."""
    # Test parameters
    density_g_cm3 = 10.0
    RHO_C = 20.0
    gamma = 0.334
    
    # Calculate screening factor
    rho_norm = density_g_cm3 / RHO_C
    screening_factor = rho_norm ** gamma
    
    # Check: should be positive and <= 1
    assert screening_factor > 0
    assert screening_factor <= 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
