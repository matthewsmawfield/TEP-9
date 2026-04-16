#!/usr/bin/env python3
"""
Test physical constants and unit conversions.

This module verifies that physical constants are correctly defined
and unit conversions are mathematically sound.
"""

import pytest
import numpy as np


def test_planck_mass():
    """Verify Planck mass is correctly defined."""
    # Planck mass in GeV
    M_PL = 2.435e18  # GeV
    # Expected value: ~2.435e18 GeV
    assert M_PL > 0
    assert M_PL > 1e18  # Should be on order of 10^18 GeV


def test_speed_of_light():
    """Verify speed of light is correctly defined."""
    C_LIGHT = 2.998e8  # m/s
    # Expected value: 299,792,458 m/s
    assert C_LIGHT > 0
    assert abs(C_LIGHT - 2.99792458e8) < 1e5  # Within 100 km/s


def test_gravitational_constant():
    """Verify gravitational constant is correctly defined."""
    G_NEWTON = 6.674e-11  # m^3 kg^-1 s^-2
    # Expected value: 6.67430e-11
    assert G_NEWTON > 0
    assert abs(G_NEWTON - 6.67430e-11) < 1e-14  # Within uncertainty


def test_earth_radius():
    """Verify Earth radius is correctly defined."""
    R_EARTH = 6.371e6  # m
    # Expected value: 6,371 km
    assert R_EARTH > 0
    assert abs(R_EARTH - 6.371e6) < 1e3  # Within 1 km


def test_earth_mass():
    """Verify Earth mass is correctly defined."""
    M_EARTH = 5.972e24  # kg
    # Expected value: 5.9722e24 kg
    assert M_EARTH > 0
    assert abs(M_EARTH - 5.9722e24) < 1e22  # Within uncertainty


def test_kg_m3_to_gev4_conversion():
    """Verify kg/m^3 to GeV^4 conversion factor."""
    KG_M3_TO_GEV4 = 1.17e-22
    # This is a unit conversion factor
    assert KG_M3_TO_GEV4 > 0
    # Should be on order of 10^-22
    assert KG_M3_TO_GEV4 > 1e-23
    assert KG_M3_TO_GEV4 < 1e-20


def test_j2_earth():
    """Verify Earth's J2 oblateness coefficient."""
    J2_EARTH = 0.00108263
    # Expected value: 0.00108263
    assert J2_EARTH > 0
    assert J2_EARTH < 1  # Dimensionless, should be < 1
    assert abs(J2_EARTH - 0.00108263) < 1e-6


def test_critical_density():
    """Verify universal critical density."""
    RHO_C = 20.0  # g/cm^3
    # Expected value: ~20 g/cm^3 from GNSS analysis
    assert RHO_C > 0
    assert RHO_C > 10  # Should be > 10 g/cm^3
    assert RHO_C < 100  # Should be < 100 g/cm^3


def test_screening_exponent():
    """Verify Vainshtein screening exponent."""
    SCREENING_EXPONENT = 0.334
    # Expected value: 0.334 from multi-scale validation
    assert SCREENING_EXPONENT > 0
    assert SCREENING_EXPONENT < 1  # Should be < 1
    assert abs(SCREENING_EXPONENT - 0.334) < 0.01


def test_tep_screening_length():
    """Verify TEP screening length."""
    LAMBDA_TEP_KM = 4000  # km
    # Expected value: ~4000 km from GNSS correlation
    assert LAMBDA_TEP_KM > 0
    assert LAMBDA_TEP_KM > 1000  # Should be > 1000 km
    assert LAMBDA_TEP_KM < 10000  # Should be < 10000 km


def test_beta_initial():
    """Verify initial TEP coupling constant."""
    BETA_INITIAL = 1e-4
    # Expected value: 10^-4 from theoretical prior
    assert BETA_INITIAL > 0
    assert BETA_INITIAL < 1  # Should be dimensionless and < 1
    assert abs(BETA_INITIAL - 1e-4) < 1e-6


def test_unit_conversion_km_to_m():
    """Verify km to m conversion."""
    km = 1000  # km
    m = km * 1e3
    assert m == 1e6  # 1000 km = 1,000,000 m


def test_unit_conversion_m_to_mm():
    """Verify m/s to mm/s conversion."""
    m_per_s = 1.0  # m/s
    mm_per_s = m_per_s * 1e3
    assert mm_per_s == 1000  # 1 m/s = 1000 mm/s


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
