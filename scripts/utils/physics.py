"""
Physics constants and fundamental units for the Temporal Equivalence Principle (TEP) framework.
Standardized for Jakarta v0.8.

References:
- Jakarta v0.8 Manuscript (2026)
- NIST CODATA 2018
- WGS84 Geodetic System
- GNSS atomic clock correlations (Paper 6, UCD)

CRITICAL NOTE: Some TEP theory parameters are empirically calibrated from GNSS data
and have significant systematic uncertainty. These are clearly marked below.
"""

from typing import Final

# Fundamental Constants
C_LIGHT: Final[float] = 299792458.0  # Speed of light in vacuum [m/s]
G_NEWTON: Final[float] = 6.67430e-11  # Newton's gravitational constant [m³ kg⁻¹ s⁻²]
H_BAR: Final[float] = 1.054571817e-34  # Reduced Planck constant [J s]

# Planck Units (Natural Units for TEP Scalar Field)
# Standardized to Reduced Planck Mass for Jakarta v0.8 consistency
M_PL_GEV: Final[float] = 2.435e18  # Reduced Planck Mass [GeV]
M_PL_KG: Final[float] = 4.341e-9    # Reduced Planck Mass [kg]

# Earth Parameters (WGS84)
R_EARTH: Final[float] = 6371000.0   # Mean Earth radius [m]
M_EARTH: Final[float] = 5.9722e24   # Earth mass [kg]
GM_EARTH: Final[float] = 3.986004418e14  # Earth gravitational parameter [m³ s⁻²]
J2_EARTH: Final[float] = 1.08263e-3  # Earth's second dynamic form factor (J2)
J3_EARTH: Final[float] = -2.53266e-6  # Earth's third dynamic form factor (J3, EGM96)
J4_EARTH: Final[float] = -1.61099e-6  # Earth's fourth dynamic form factor (J4, EGM96)

# Astronomy Units
AU_METERS: Final[float] = 1.495978707e11  # Astronomical Unit [m]
PC_METERS: Final[float] = 3.08567758e16   # Parsec [m]
LY_METERS: Final[float] = 9.46073047e15   # Light year [m]

# TEP Baseline Theory Parameters (Jakarta v0.8)
LAMBDA_BASELINE_GEV: Final[float] = 1.0e-2  # 10 MeV scale [GeV]
BETA_BASELINE: Final[float] = 1.0           # Unit scalar coupling
N_TOPOLOGY: Final[int] = 3                   # Continuous gradient suppression index

# TEP Field Relaxation and Screening (Empirical Anchors)
# CRITICAL: These parameters are empirically calibrated from GNSS atomic clock correlations
# and have significant systematic uncertainty. They are NOT fundamental constants.
LAMBDA_TEP_M: Final[float] = 4200000.0       # Unified relaxation length [m] (±15% from GNSS calibration)
R_TRANSITION_M: Final[float] = 4146000.0     # PREM-derived transition radius [m] (±10% from PREM model)
RHO_T: Final[float] = 20.0                   # Temporal Topology saturation density [g/cm³] (±40% from GNSS calibration)
RHO_T_ERROR: Final[float] = 8.0             # 40% systematic uncertainty from GNSS calibration
RHO_T_SOURCE: Final[str] = "GNSS atomic clock correlations (Paper 6, UCD) - ±40% uncertainty"
LAMBDA_TEP_UNCERTAINTY: Final[float] = 0.15  # ±15% systematic uncertainty
R_TRANSITION_UNCERTAINTY: Final[float] = 0.10  # ±10% systematic uncertainty
SUPPRESSION_EXPONENT: Final[float] = 0.334    # Theoretical density scaling (≈ 1/3, ±0.05 from multi-scale validation)

# Disformal Coupling Parameters (Yogyakarta v0.1)
# CRITICAL: These are theoretical parameters with limited empirical validation
DISFORMAL_COUPLING_STRENGTH: Final[float] = 0.05  # α_B - Disformal coupling strength (±50% uncertainty - theoretical)
DISFORMAL_VELOCITY_THRESHOLD_KM_S: Final[float] = 16.8  # v_trans - Velocity threshold [km/s] (±20% uncertainty - theoretical)
DISFORMAL_COUPLING_UNCERTAINTY: Final[float] = 0.5  # ±50% uncertainty (theoretical)
DISFORMAL_VELOCITY_UNCERTAINTY: Final[float] = 0.2  # ±20% uncertainty (theoretical)

# Legacy alias for backward compatibility
RHO_TEMPORAL_TOPOLOGY_G_CM3: Final[float] = RHO_T  # Universal Temporal Topology density [g/cm³] (deprecated, use RHO_T)

# Screening factors (Paper 6 UCD + Paper 15 EFA)
# S_⊕: surface gradient-suppression ratio used in flyby β_eff = β × S_⊕
# S_ucd: embedding depth R_sol / R_⊕ used in Paper 6 (distinct from S_⊕)
CHARACTERISTIC_SUPPRESSION: Final[float] = (R_EARTH - R_TRANSITION_M) / R_EARTH
UCD_EMBEDDING_FACTOR: Final[float] = R_TRANSITION_M / R_EARTH

# Conversion Factors
KG_M3_TO_GEV4: Final[float] = 4.318e-21  # kg/m³ to GeV⁴ conversion (Natural Units)


def screened_beta(beta: float, surface_suppression: float = CHARACTERISTIC_SUPPRESSION) -> float:
    """Map universal conformal β to Earth-surface screened β_eff."""
    return beta * surface_suppression


def ppn_gamma_deviation(beta_eff: float) -> float:
    """Magnitude |γ − 1| used for Cassini checks (Jakarta v0.8 Sec. 7: γ − 1 = −2 α_eff² in DEF; here |γ−1| ≈ 2 β_eff² with β_eff the screened dimensionless coupling for A = exp(β φ/M_Pl))."""
    return 2.0 * beta_eff**2


def validate_screened_coupling(
    beta: float,
    beta_eff: float,
    *,
    surface_suppression: float = CHARACTERISTIC_SUPPRESSION,
    rtol: float = 0.02,
) -> None:
    expected = screened_beta(beta, surface_suppression)
    scale = max(abs(beta_eff), 1e-30)
    if abs(expected - beta_eff) / scale > rtol:
        raise ValueError(
            "β_eff is inconsistent with β × S_⊕: "
            f"beta={beta:.6e}, beta_eff={beta_eff:.6e}, "
            f"expected={expected:.6e}, S_earth={surface_suppression:.6f}"
        )

def get_tep_metadata() -> dict:
    """Return theory version and metadata for data provenance."""
    return {
        "theory_version": "Jakarta v0.8",
        "paradigm": "Temporal Shear Suppression (Continuous Gradient)",
        "coupling_convention": "A(phi) = exp(beta phi / M_Pl)",
        "mpl_definition": "Reduced Planck Mass (2.435e18 GeV)",
        "standard_constants": {
            "Lambda_TEP": f"{LAMBDA_TEP_M/1e3} km",
            "R_transition": f"{R_TRANSITION_M/1e3} km",
            "rho_T": f"{RHO_T} g/cm³",
            "rho_T_error": f"{RHO_T_ERROR} g/cm³",
            "rho_T_source": RHO_T_SOURCE,
            "S_earth": f"{CHARACTERISTIC_SUPPRESSION:.3f}",
            "S_ucd_embedding": f"{UCD_EMBEDDING_FACTOR:.3f}",
        }
    }
