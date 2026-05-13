#!/usr/bin/env python3
"""
Step 007: TEP Scalar Force Model (Analytical Reference)

This module provides the analytical framework for the TEP scalar force.
NOTE: This script uses the perigee approximation (impulse approximation) 
for preliminary analysis. Rigorous 3D path integration is performed in 
Step 009 (trajectory_integration.py) to account for directional 
heterogeneity and exact ∇φ sampling.
"""

import datetime
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.physics import (
    C_LIGHT,
    G_NEWTON,
    R_EARTH,
    M_EARTH,
    J2_EARTH,
    M_PL_GEV as M_PL,
    CHARACTERISTIC_SUPPRESSION,
    BETA_BASELINE,
    KG_M3_TO_GEV4,
    LAMBDA_BASELINE_GEV,
    LAMBDA_TEP_M,
    R_TRANSITION_M,
    RHO_T,
    SUPPRESSION_EXPONENT,
    DISFORMAL_COUPLING_STRENGTH,
    DISFORMAL_VELOCITY_THRESHOLD_KM_S,
    get_tep_metadata
)

# NOTE: Local specialized density/geoid classes remain in enhanced_physics
from scripts.utils.enhanced_physics import (
    EarthDensityModel,
    EarthGeoidModel,
)

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep_geometry_envelope import (
    compose_geometry_envelope,
    disformal_envelope_factor,
    j2_only_bracket,
    zonal_harmonic_bracket,
)

# TEP Theoretical Framework (Standardized for Jakarta v0.8)
# These constants are centrally managed in scripts/utils/physics.py
LAMBDA_TEP_KM = LAMBDA_TEP_M / 1e3
R_TRANSITION_KM = R_TRANSITION_M / 1e3
# Note: CHARACTERISTIC_SUPPRESSION (S_⊕) is the UCD surface gradient-suppression ratio
RELAXATION_EXPONENT = SUPPRESSION_EXPONENT

# Load config
try:
    with open(
        PROJECT_ROOT / "config" / "pipeline_config.json", "r", encoding="utf-8"
    ) as f:
        config = json.load(f)
    tep_config = config["parameters"]["analysis"]["tep_physics"]
    N_TOPOLOGY = tep_config["n_temporal_shear_suppression"]
    LAMBDA_GEV = tep_config["lambda_gev"]
    BETA_INITIAL = BETA_BASELINE * 1e-4  # Unified Yogyakarta anchor from physics.py
    GEOMETRIC_SCREENING_FACTOR = tep_config.get(
        "geometric_screening_factor", CHARACTERISTIC_SUPPRESSION
    )
except (FileNotFoundError, json.JSONDecodeError, KeyError, PermissionError) as e:
    # Use standardized Jakarta v0.8.0 defaults from physics.py
    N_TOPOLOGY = 3
    LAMBDA_GEV = 0.01
    BETA_INITIAL = BETA_BASELINE * 1e-4
    GEOMETRIC_SCREENING_FACTOR = CHARACTERISTIC_SUPPRESSION

GM_EARTH = 3.986004418e14  # m³/s²

# Disformal coupling parameters (from physics.py)
# DISFORMAL_COUPLING_STRENGTH and DISFORMAL_VELOCITY_THRESHOLD_KM_S are imported from physics.py

# CMB dipole parameters (Planck 2018) — scalar field rest-frame anchor
CMB_DIPOLE_RA_DEG = 167.94
CMB_DIPOLE_DEC_DEG = -6.93
CMB_DIPOLE_VELOCITY_KM_S = 369.82

# Earth orbital parameters for CMB-frame velocity computation
EARTH_ORBITAL_ECCENTRICITY = 0.0167086
EARTH_ORBITAL_SEMI_MAJOR_AU = 1.000001018
EARTH_ORBITAL_LONGITUDE_OF_PERIHELION_DEG = 102.937348
AU_M = 1.495978707e11
M_SUN_KG = 1.98847e30
G_MKS = 6.67430e-11


def _julian_day(dt):
    y, m, d = dt.year, dt.month, dt.day
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    JD = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5
    frac = (dt.hour + dt.minute / 60.0 + dt.second / 3600.0) / 24.0
    return JD + frac


def _earth_orbital_velocity_equatorial(dt):
    """
    Compute Earth's heliocentric orbital velocity vector in J2000 equatorial frame.
    Returns (vx, vy, vz) in km/s.
    """
    JD = _julian_day(dt)
    T = (JD - 2451545.0) / 36525.0

    L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T * T
    L0 = L0 % 360.0

    M = math.radians((357.52911 + 35999.05029 * T - 0.0001537 * T * T) % 360.0)
    e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T * T

    C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * math.sin(M)
    C += (0.019993 - 0.000101 * T) * math.sin(2 * M)
    C += 0.000289 * math.sin(3 * M)

    true_lon_sun = (L0 + C) % 360.0
    true_lon_earth = (true_lon_sun + 180.0) % 360.0
    true_lon_earth_rad = math.radians(true_lon_earth)

    nu_rad = math.radians((true_lon_earth - EARTH_ORBITAL_LONGITUDE_OF_PERIHELION_DEG) % 360.0)
    a = EARTH_ORBITAL_SEMI_MAJOR_AU
    r_au = a * (1 - e * e) / (1 + e * math.cos(nu_rad))

    mu_km3_s2 = 1.32712440018e11
    r_km = r_au * AU_M / 1000.0
    p_km = a * AU_M / 1000.0 * (1 - e * e)
    h = math.sqrt(mu_km3_s2 * p_km)

    v_r = (mu_km3_s2 / h) * e * math.sin(nu_rad)
    v_t = (mu_km3_s2 / h) * (1 + e * math.cos(nu_rad))

    vx_ecl = v_r * math.cos(true_lon_earth_rad) - v_t * math.sin(true_lon_earth_rad)
    vy_ecl = v_r * math.sin(true_lon_earth_rad) + v_t * math.cos(true_lon_earth_rad)
    vz_ecl = 0.0

    eps = math.radians(23.439281)
    vx_eq = vx_ecl
    vy_eq = vy_ecl * math.cos(eps) - vz_ecl * math.sin(eps)
    vz_eq = vy_ecl * math.sin(eps) + vz_ecl * math.cos(eps)

    return vx_eq, vy_eq, vz_eq


def _cmb_dipole_unit_vector():
    ra = math.radians(CMB_DIPOLE_RA_DEG)
    dec = math.radians(CMB_DIPOLE_DEC_DEG)
    return np.array([
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec)
    ])


def _sc_velocity_unit_vector_equatorial(dec_in_deg):
    dec = math.radians(dec_in_deg)
    return np.array([math.cos(dec), 0.0, math.sin(dec)])


def compute_sc_cmb_cos_theta(dec_in_deg, state_vectors=None, mission_name=None):
    """
    Compute cos(theta) between SC velocity direction and CMB dipole apex.

    Uses 3D state vectors from JPL Horizons (step_040a) if available.
    Falls back to declination-only approximation otherwise.
    """
    n_cmb = _cmb_dipole_unit_vector()

    # Try 3D state vectors first
    if state_vectors and mission_name:
        vec = state_vectors.get(mission_name)
        if vec is not None and "ux_ecl" in vec:
            ux_ecl = vec["ux_ecl"]
            uy_ecl = vec["uy_ecl"]
            uz_ecl = vec["uz_ecl"]
            eps = math.radians(23.439281)
            ux_eq = ux_ecl
            uy_eq = uy_ecl * math.cos(eps) + uz_ecl * math.sin(eps)
            uz_eq = -uy_ecl * math.sin(eps) + uz_ecl * math.cos(eps)
            v_sc_hat = np.array([ux_eq, uy_eq, uz_eq])
            return float(np.dot(v_sc_hat, n_cmb))

    # Fallback to declination-only
    if dec_in_deg is not None:
        v_sc_hat = _sc_velocity_unit_vector_equatorial(dec_in_deg)
        return float(np.dot(v_sc_hat, n_cmb))

    return None


def compute_cmb_frame_speed_kms(flyby_date_str, v_sc_kms, dec_in_deg):
    """
    Compute the total CMB-frame speed for a flyby.

    In TEP, the disformal coupling depends on velocity relative to the scalar
    field rest frame. If the CMB frame approximates the scalar rest frame, the
    effective velocity is the vector sum:

        v_cmb = v_CMB_dipole + v_Earth_orbit + v_SC

    Returns the magnitude |v_cmb| in km/s, or None if flyby_date_str is empty
    (caller should fall back to geocentric velocity).
    """
    if not flyby_date_str:
        return None
    try:
        dt = datetime.datetime.strptime(flyby_date_str, "%Y-%m-%d")
    except ValueError:
        try:
            dt = datetime.datetime.strptime(flyby_date_str[:10], "%Y-%m-%d")
        except ValueError:
            return None

    n_cmb = _cmb_dipole_unit_vector()
    v_cmb_bulk = CMB_DIPOLE_VELOCITY_KM_S * n_cmb

    v_earth = np.array(_earth_orbital_velocity_equatorial(dt))

    v_sc_hat = _sc_velocity_unit_vector_equatorial(dec_in_deg)
    v_sc = v_sc_hat * v_sc_kms

    v_total = v_cmb_bulk + v_earth + v_sc
    return float(np.linalg.norm(v_total))


def _vinf_declinations_from_ephemeris(trajectory_data):
    ephemeris = trajectory_data.get("ephemeris", [])
    if not ephemeris or len(ephemeris) < 10:
        raise ValueError("Insufficient ephemeris data for trajectory reconstruction.")
    ranges = [p["range_km"] for p in ephemeris]
    perigee_idx = int(np.argmin(ranges))
    p = ephemeris[perigee_idx]
    r_vec = np.array([p["x_km"], p["y_km"], p["z_km"]]) * 1e3
    v_vec = np.array([p["vx_km_s"], p["vy_km_s"], p["vz_km_s"]]) * 1e3
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    e_vec = np.cross(v_vec, h_vec) / GM_EARTH - r_vec / r
    e = np.linalg.norm(e_vec)
    if e < 1.001:
        raise ValueError(f"Orbit is not hyperbolic (e={e:.4f}). declination extraction requires hyperbolic flyby.")
    x_hat = e_vec / e
    h_hat = h_vec / h
    y_hat = np.cross(h_hat, x_hat)
    sqrt_e2_1 = math.sqrt(e**2 - 1)
    v_in_hat = (x_hat + sqrt_e2_1 * y_hat) / e
    v_out_hat = (-x_hat + sqrt_e2_1 * y_hat) / e
    dec_in = math.degrees(math.asin(np.clip(v_in_hat[2], -1, 1)))
    dec_out = math.degrees(math.asin(np.clip(v_out_hat[2], -1, 1)))
    return dec_in, dec_out


def extract_trajectory_geometry(flyby_data, spacecraft_name=""):
    """
    Extract trajectory geometry from flyby data.
    
    Priority:
    1. If ``data/trajectories/{name}_ephemeris.csv`` exists, load it (required to succeed).
    2. Otherwise use published declinations from the flyby catalog.
    """
    name = spacecraft_name or flyby_data.get("mission_name", "Unknown")
    
    # Try ephemeris reconstruction first
    ephemeris_path = PROJECT_ROOT / "data" / "trajectories" / f"{name.lower()}_ephemeris.csv"
    
    dec_in_deg = flyby_data.get("declination_in_deg")
    dec_out_deg = flyby_data.get("declination_out_deg")
    declination_source = "catalog_published"
    
    ephemeris = []
    perigee_lat = 0.0
    
    if ephemeris_path.exists():
        try:
            trajectory_data = _load_ephemeris(ephemeris_path)
            dec_in_deg, dec_out_deg = _vinf_declinations_from_ephemeris(trajectory_data)
            declination_source = "reconstructed_ephemeris"
            ephemeris = trajectory_data.get("ephemeris", [])
            
            # Extract perigee latitude if ephemeris is available
            if ephemeris:
                ranges = [p["range_km"] for p in ephemeris]
                p_idx = int(np.argmin(ranges))
                p = ephemeris[p_idx]
                r_p = math.sqrt(p["x_km"]**2 + p["y_km"]**2 + p["z_km"]**2)
                if r_p > 0:
                    perigee_lat = math.asin(p["z_km"] / r_p)
        except Exception as e:
            raise RuntimeError(
                f"Ephemeris file is present but could not be used for {name}: {ephemeris_path}"
            ) from e

    if dec_in_deg is None or dec_out_deg is None:
        raise ValueError(f"Insufficient geometry data for {name} (no ephemeris or catalog declinations).")

    v_p = flyby_data.get("perigee_velocity_km_s")
    alt_p = flyby_data.get("perigee_altitude_km")
    
    if v_p is None or alt_p is None:
        raise ValueError("Missing critical perigee parameters (velocity or altitude).")

    dec_in_rad = math.radians(dec_in_deg)
    dec_out_rad = math.radians(dec_out_deg)
    # Use catalog cos_asymmetry if present (from manuscript), otherwise calculate
    cos_dec_asymmetry = flyby_data.get("cos_asymmetry")
    if cos_dec_asymmetry is None:
        cos_dec_asymmetry = math.cos(dec_in_rad) - math.cos(dec_out_rad)

    flyby_date = flyby_data.get("flyby_date", "")

    return {
        "dec_in_deg": dec_in_deg,
        "dec_out_deg": dec_out_deg,
        "dec_in_rad": dec_in_rad,
        "dec_out_rad": dec_out_rad,
        "cos_dec_asymmetry": cos_dec_asymmetry,
        "declination_source": declination_source,
        "v_perigee_m_s": v_p * 1e3,
        "altitude_km": alt_p,
        "perigee_latitude_rad": perigee_lat,
        "ephemeris": ephemeris,
        "flyby_date": flyby_date,
    }


class TEPTemporalTopologyModel:
    def __init__(
        self,
        beta=BETA_INITIAL,
        Lambda_GeV=LAMBDA_GEV,
        n=N_TOPOLOGY,
        characteristic_suppression=CHARACTERISTIC_SUPPRESSION,
    ):
        self.beta = beta
        self.Lambda = Lambda_GeV
        self.n = n
        self.characteristic_suppression = characteristic_suppression
        self.lambda_tep = LAMBDA_TEP_M
        
        # =====================================================================
        # SELF-CONSISTENT FIELD PARAMETERS (Jakarta v0.8 field equations)
        # =====================================================================
        # phi is computed from the Jakarta v0.8 / EFA field minimum (Temporal Shear Suppression):
        #   phi(rho) = Lambda * [n * Lambda^(n+4) * M_Pl / (2 * beta * rho_GeV4)]^(1/(n+1))
        # Matter coupling supplies the factor 2 in the denominator (Einstein-frame
        # TSS implementation shared with Step 011 / Step 019; EFA v0.1 Yogyakarta text).
        # with n=3 (Temporal Topology index), Lambda=10 MeV.
        #
        # Earth interior density ~5515 kg/m³ gives phi_earth ~ 2.4×10⁴ GeV.
        # Interplanetary vacuum ~1e-20 kg/m³ gives phi_space ~ 2.0×10¹⁰ GeV.
        # These are computed self-consistently for the reference beta=1e-4.
        # =====================================================================
        
        self.phi_earth = self._phi_of_rho(5515.0)
        self.phi_surface = self._phi_of_rho(2700.0)
        self.phi_space = self._phi_of_rho(1e-20)
        self.delta_phi = self.phi_space - self.phi_earth

    def _phi_of_rho(self, rho_kg_m3):
        rho_gev4 = rho_kg_m3 * KG_M3_TO_GEV4
        if rho_gev4 <= 0 or self.beta <= 0:
            return self.Lambda * 1e9
        # Jakarta v0.8 consistency: use factor of 2 in denominator for field minimum
        numerator = self.n * (self.Lambda ** (4 + self.n)) * M_PL
        denominator = 2.0 * self.beta * rho_gev4
        if denominator <= 0:
            return self.Lambda * 1e9
        scale = (numerator / denominator) ** (1.0 / (self.n + 1))
        return self.Lambda * scale

    def phi(self, r):
        if r <= R_EARTH:
            return self.phi_earth
        delta_r = r - R_EARTH
        frac = 1.0 - np.exp(-delta_r / self.lambda_tep)
        return self.phi_earth + self.delta_phi * frac

    def field_gradient(self, r):
        if r <= R_EARTH:
            return 0.0
        delta_r = r - R_EARTH
        return (self.delta_phi / self.lambda_tep) * np.exp(-delta_r / self.lambda_tep)

    def geometry_modulation_factors(
        self, altitude_km, latitude_deg, velocity_km_s, plasma_density_cm3=1000
    ):
        f_inclination = 1.0 + 0.15 * abs(np.sin(np.radians(latitude_deg)))
        f_j2 = (1.0 - 0.00054 * np.cos(np.radians(latitude_deg)) ** 2) * np.exp(
            -altitude_km / 2000.0
        )
        f_plasma = (1.0 + plasma_density_cm3 / 5000.0) ** (-0.3)
        v_threshold = DISFORMAL_VELOCITY_THRESHOLD_KM_S
        f_velocity = (
            (v_threshold / velocity_km_s) ** 4.0 if velocity_km_s > v_threshold else 1.0
        )
        f_total_core = f_inclination * f_j2 * f_plasma * f_velocity
        return {
            "f_inclination": f_inclination,
            "f_j2": f_j2,
            "f_plasma": f_plasma,
            "f_velocity": f_velocity,
            "f_total_core": f_total_core,
        }

    def disformal_modulation_factor(
        self, v_sc_m_s, cos_asymmetry, v_cmb_frame_kms=None, sc_cmb_cos_theta=None
    ):
        return disformal_envelope_factor(
            v_sc_m_s,
            cos_asymmetry,
            v_cmb_frame_kms=v_cmb_frame_kms,
            sc_cmb_cos_theta=sc_cmb_cos_theta,
        )

    def _plasma_density_cm3(self, altitude_km):
        if altitude_km < 1000:
            return 10000.0 * np.exp(-altitude_km / 200.0)
        return 100.0

    def _geometry_context(self, geometry, plasma_effects=None):
        altitude_km = geometry["altitude_km"]
        v_p = geometry["v_perigee_m_s"]
        cos_asymmetry = geometry["cos_dec_asymmetry"]
        perigee_lat_deg = np.degrees(geometry["perigee_latitude_rad"])
        velocity_km_s = v_p / 1000.0
        r_p = altitude_km * 1e3 + R_EARTH

        if plasma_effects is None:
            raise ValueError("Missing perigee plasma effects for deterministic geometry envelope.")
        plasma_density = plasma_effects["plasma_density_cm3"]
        plasma_screening = plasma_effects["plasma_screening_factor"]
        plasma_sign = plasma_effects["plasma_sign_factor"]

        modulation = self.geometry_modulation_factors(
            altitude_km, perigee_lat_deg, velocity_km_s, plasma_density
        )
        envelope = compose_geometry_envelope(
            altitude_km=altitude_km,
            latitude_deg=perigee_lat_deg,
            velocity_km_s=velocity_km_s,
            cos_asymmetry=cos_asymmetry,
            plasma_density_cm3=plasma_density,
            modulation=modulation,
            plasma_screening_factor=plasma_screening,
            plasma_sign_factor=plasma_sign,
        )

        v_cmb = compute_cmb_frame_speed_kms(
            geometry.get("flyby_date", ""),
            velocity_km_s,
            geometry.get("dec_in_deg", 0.0),
        )
        sc_cmb_cos = compute_sc_cmb_cos_theta(
            geometry.get("dec_in_deg"),
            geometry.get("state_vectors"),
            geometry.get("mission_name"),
        )
        s_disf = self.disformal_modulation_factor(
            v_p, cos_asymmetry, v_cmb, sc_cmb_cos
        )

        return {
            "r_p": r_p,
            "v_p": v_p,
            "cos_asymmetry": cos_asymmetry,
            "perigee_lat_deg": perigee_lat_deg,
            "harmonic_bracket": zonal_harmonic_bracket(perigee_lat_deg, r_p),
            "j2_bracket": j2_only_bracket(perigee_lat_deg, r_p),
            "envelope": envelope,
            "s_disf": s_disf,
        }

    def tep_velocity_shift(
        self, geometry, plasma_effects=None
    ):
        ephemeris = geometry.get("ephemeris", [])
        ctx = self._geometry_context(geometry, plasma_effects)
        r_p = ctx["r_p"]
        v_p = ctx["v_p"]
        cos_asymmetry = ctx["cos_asymmetry"]
        envelope_factor = ctx["envelope"]["geometry_envelope"]
        s_disf = ctx["s_disf"]

        if not ephemeris:
            beta_eff = self.effective_coupling(r_p)
            dphi_dr = self.field_gradient(r_p)
            dv_base = (
                beta_eff
                * (C_LIGHT**2)
                * dphi_dr
                / M_PL
                * (r_p / v_p)
                * ctx["harmonic_bracket"]
                * cos_asymmetry
                * envelope_factor
                * 1e3
            )
            return dv_base * s_disf

        integral_F_dt = 0.0
        for i in range(len(ephemeris) - 1):
            p1, p2 = ephemeris[i], ephemeris[i + 1]
            try:
                t1 = datetime.datetime.strptime(
                    p1["datetime"].replace("A.D. ", "").strip(), "%Y-%b-%d %H:%M:%S.%f"
                )
                t2 = datetime.datetime.strptime(
                    p2["datetime"].replace("A.D. ", "").strip(), "%Y-%b-%d %H:%M:%S.%f"
                )
            except ValueError:
                t1 = datetime.datetime.strptime(
                    p1["datetime"].replace("A.D. ", "").strip(), "%Y-%b-%d %H:%M:%S"
                )
                t2 = datetime.datetime.strptime(
                    p2["datetime"].replace("A.D. ", "").strip(), "%Y-%b-%d %H:%M:%S"
                )
            dt = (t2 - t1).total_seconds()
            if dt <= 0:
                continue
            r = np.linalg.norm(
                [
                    (p1["x_km"] + p2["x_km"]) / 2 * 1e3,
                    (p1["y_km"] + p2["y_km"]) / 2 * 1e3,
                    (p1["z_km"] + p2["z_km"]) / 2 * 1e3,
                ]
            )
            if r <= R_EARTH:
                continue
            lat_deg = math.degrees(math.asin(((p1["z_km"] + p2["z_km"]) / 2) / (r / 1e3)))
            bracket = zonal_harmonic_bracket(lat_deg, r)
            integral_F_dt += (
                self.effective_coupling(r)
                * (C_LIGHT**2)
                * self.field_gradient(r)
                / M_PL
                * bracket
                * cos_asymmetry
                * envelope_factor
                * dt
            )

        return integral_F_dt * s_disf * 1e3

    def ambient_relaxation_factor(self, density_g_cm3):
        if density_g_cm3 <= 0:
            return 1.0
        rho_norm = density_g_cm3 / RHO_T
        return min(rho_norm ** (-RELAXATION_EXPONENT), 1.0)

    def local_density_at_altitude(self, altitude_km):
        if altitude_km < 0:
            return 5.515
        if altitude_km > 1000:
            return 1e-20
        return 1.225e-3 * np.exp(-altitude_km / 8.5)

    def effective_coupling(self, r):
        """
        Continuous density-driven screening via Temporal Shear suppression (v0.8).
        Aligned with Jakarta v0.8: A(phi) = exp(beta*phi/M_Pl).
        """
        altitude_km = (r - R_EARTH) / 1e3
        return self.beta * self.ambient_relaxation_factor(
            self.local_density_at_altitude(altitude_km)
        )


def main():
    """Execute TEP model predictions for all archival flybys."""
    logger = StepLogger("step_007_tep_model", PROJECT_ROOT)
    start_time = time.time()

    logger.header("STEP 007: TEP TEMPORAL TOPOLOGY MODEL")

    # Load flyby catalog
    catalog_path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
    if not catalog_path.exists():
        logger.error(f"Archival catalog not found: {catalog_path}")
        return 1

    with open(catalog_path, "r", encoding="utf-8") as f:
        catalog = json.load(f)

    flybys = catalog.get("flybys", [])
    logger.info(f"Loaded {len(flybys)} flybys from archival catalog")

    # Use theoretical reference coupling BETA_INITIAL = 1e-04.
    # Step_008 fits empirical scaling factors relative to this reference.
    # Loading fitted beta here would break consistency with step_008's formula.
    beta_ref = BETA_INITIAL
    logger.info(f"Using reference beta = {beta_ref:.4e}")

    model = TEPTemporalTopologyModel(beta=beta_ref)
    predictions = {}

    from scripts.steps.step_017_plasma_modulation import PlasmaModulationModel

    plasma_model = PlasmaModulationModel()

    # Load 3D state vectors for CMB alignment computation
    state_vectors_path = PROJECT_ROOT / "results" / "step040a_3d_state_vectors.json"
    state_vectors = {}
    if state_vectors_path.exists():
        with open(state_vectors_path, "r", encoding="utf-8") as f:
            state_vectors = json.load(f)
        logger.info(f"Loaded {len(state_vectors)} 3D state vectors for CMB alignment")

    for flyby in flybys:
        name = flyby["mission_name"]
        logger.subheader(f"Predicting: {name}")

        try:
            # Extract geometry
            geometry = extract_trajectory_geometry(flyby, name)
            geometry["state_vectors"] = state_vectors
            geometry["mission_name"] = name

            perigee_vector = state_vectors.get(name)
            if perigee_vector and geometry["perigee_latitude_rad"] == 0.0:
                geometry["perigee_latitude_rad"] = math.radians(perigee_vector["dec_deg"])

            plasma_effects = plasma_model.calculate_plasma_effects(
                {
                    "name": name,
                    "flyby_date": flyby.get("flyby_date", geometry.get("flyby_date", "")),
                    "altitude_km": geometry["altitude_km"],
                }
            )

            ctx = model._geometry_context(geometry, plasma_effects)
            dv_tep = model.tep_velocity_shift(geometry, plasma_effects)
        except (ValueError, FileNotFoundError, RuntimeError) as e:
            logger.warning(f"  Skipping {name}: {str(e)}")
            continue
        except Exception as e:
            logger.error(f"  Unexpected error modeling {name}: {str(e)}")
            raise

        r_p = ctx["r_p"]
        beta_eff = model.effective_coupling(r_p)
        dphi_dr = model.field_gradient(r_p)
        v_p = ctx["v_p"]
        cos_asym_factor = ctx["cos_asymmetry"]
        envelope_factor = ctx["envelope"]["geometry_envelope"]

        dv_grad = dv_tep / ctx["s_disf"] if ctx["s_disf"] else dv_tep
        dv_total = dv_tep
        dv_disf = dv_total - dv_grad

        predictions[name] = {
            "spacecraft": name,
            "perigee": {
                "altitude_km": flyby["perigee_altitude_km"],
                "velocity_km_s": flyby["perigee_velocity_km_s"],
                "datetime": flyby.get("perigee_time", flyby["flyby_date"]),
            },
            "observed": {
                "dv_obs_mm_s": flyby.get("published_anomaly_mm_s"),
                "sigma_mm_s": flyby.get("published_anomaly_uncertainty_mm_s"),
            },
            "tep_predictions": {
                "dv_tep_mm_s": dv_tep,
                "dv_grad_mm_s": dv_grad,
                "dv_disf_mm_s": dv_disf,
                "model_version": "v5.4-Yogyakarta-geometry-envelope",
                "topology_n": N_TOPOLOGY,
                "beta_reference": beta_ref,
                "geometry_envelope": ctx["envelope"],
            },
            "geometry": geometry,
        }

        obs_dv = flyby.get('published_anomaly_mm_s')
        if obs_dv is None:
            logger.warning(f"  No published anomaly for {name}. Using N/A.")
            
        logger.info(f"  Predicted Δv: {dv_tep:.4f} mm/s")
        logger.info(
            f"  Observed Δv:  {('N/A' if obs_dv is None else f'{obs_dv:.2f}')} mm/s"
        )

    # Save predictions
    results_dir = PROJECT_ROOT / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    output_path = results_dir / "step007_tep_predictions.json"

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump({"predictions": predictions}, f, indent=2)

    logger.success(f"Saved {len(predictions)} predictions to: {output_path}")
    logger.add_output_file(output_path, "TEP predictions")

    duration = time.time() - start_time
    logger.log_step_summary(duration, "SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
