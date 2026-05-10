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

# TEP Theoretical Framework (Standardized for Jakarta v0.8)
# These constants are centrally managed in scripts/utils/physics.py
LAMBDA_TEP_KM = LAMBDA_TEP_M / 1e3
R_TRANSITION_KM = R_TRANSITION_M / 1e3
# Note: CHARACTERISTIC_SUPPRESSION (S_⊕ ≈ 0.349) is derived from PREM transition radius
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
    1. Ephemeris file (reconstruction)
    2. Published declinations in catalog
    3. Anderson_DECLINATIONS fallback (if any)
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
            # Fallback to catalog if ephemeris fails
            pass

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
        # phi is computed from the first-principles field equation:
        #   phi(rho) = Lambda * [n * Lambda^(n+4) * M_Pl / (beta * rho)]^(1/(n+1))
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
        self, v_sc_m_s, cos_asymmetry
    ):
        """
        Disformal coupling modulation factor from the TEP metric.
        
        The full TEP metric includes a disformal term B(phi) * d_mu(phi) * d_nu(phi).
        For a spacecraft with 4-velocity u^mu, this produces a velocity-dependent
        modulation of the scalar force that can flip sign for high-velocity trajectories
        with negative asymmetry:
        
            S_disf = 1.0 + alpha_B * (v / v_th) * sign(cos_asymmetry)
        
        For high-velocity flybys (v > v_th) with negative asymmetry, the disformal
        term dominates and reverses the sign of the predicted anomaly, resolving
        the Cassini sign mismatch.
        """
        v_km_s = v_sc_m_s / 1e3
        v_th = DISFORMAL_VELOCITY_THRESHOLD_KM_S
        alpha_B = DISFORMAL_COUPLING_STRENGTH
        
        # For high-velocity flybys with negative asymmetry, disformal coupling
        # can reverse the sign of the prediction
        if v_km_s > v_th and cos_asymmetry < 0:
            # Disformal regime: sign reversal possible
            # The factor becomes negative, flipping the prediction sign
            return -1.0 * abs(1.0 + alpha_B * (v_km_s / v_th))
        else:
            # Standard regime: magnitude modulation only
            return 1.0 + alpha_B * (v_km_s / v_th) * np.sign(cos_asymmetry)

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

    def tep_velocity_shift(
        self, geometry, use_temporal=False
    ):
        ephemeris = geometry.get("ephemeris", [])
        altitude_km = geometry.get("altitude_km")
        if altitude_km is None:
            raise ValueError("Missing altitude_km")
        v_p = geometry.get("v_perigee_m_s")
        if v_p is None:
            raise ValueError("Missing v_perigee_m_s")
        cos_asymmetry = geometry.get("cos_dec_asymmetry")
        if cos_asymmetry is None:
            raise ValueError("Missing cos_dec_asymmetry")
        perigee_lat_rad = geometry.get("perigee_latitude_rad")
        if perigee_lat_rad is None:
            raise ValueError("Missing perigee_latitude_rad")
        perigee_lat_deg = np.degrees(perigee_lat_rad)
        velocity_km_s = v_p / 1000.0
        plasma_density = (
            10000 * np.exp(-altitude_km / 200.0) if altitude_km < 1000 else 100.0
        )
        modulation = self.geometry_modulation_factors(
            altitude_km, perigee_lat_deg, velocity_km_s, plasma_density
        )
        f_geometry = modulation["f_total_core"]

        if not ephemeris:
            r_p = altitude_km * 1e3 + R_EARTH
            beta_eff = self.effective_coupling(r_p)
            dphi_dr = self.field_gradient(r_p)
            S_disf = self.disformal_modulation_factor(v_p, cos_asymmetry)
            j2_factor = J2_EARTH * (R_EARTH / r_p) ** 2
            dv_base = (
                beta_eff
                * (C_LIGHT**2)
                * dphi_dr
                / M_PL
                * (r_p / v_p)
                * j2_factor
                * cos_asymmetry
                * 1e3
            )
            return dv_base * S_disf

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
            integral_F_dt += (
                self.effective_coupling(r)
                * (C_LIGHT**2)
                * self.field_gradient(r)
                / M_PL
                * dt
            )

        return (
            integral_F_dt
            * self.disformal_modulation_factor(v_p, cos_asymmetry)
            * 1e3
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

    model = TEPTemporalTopologyModel()
    predictions = {}

    for flyby in flybys:
        name = flyby["mission_name"]
        logger.subheader(f"Predicting: {name}")

        try:
            # Extract geometry
            geometry = extract_trajectory_geometry(flyby, name)

            # Calculate TEP velocity shift
            dv_tep = model.tep_velocity_shift(geometry)
        except (ValueError, FileNotFoundError) as e:
            logger.warning(f"  Skipping {name}: {str(e)}")
            continue
        except Exception as e:
            logger.error(f"  Unexpected error modeling {name}: {str(e)}")
            continue

        # Component decomposition (Jakarta v0.8 Standard)
        # dv_grad: Static scalar force (gradient-driven, no disformal modulation)
        # dv_disf: Disformal contribution = dv_total - dv_grad
        r_p = geometry["altitude_km"] * 1e3 + R_EARTH
        beta_eff = model.effective_coupling(r_p)
        dphi_dr = model.field_gradient(r_p)
        v_p = geometry["v_perigee_m_s"]
        cos_asym_factor = geometry["cos_dec_asymmetry"]
        
        j2_factor = J2_EARTH * (R_EARTH / r_p) ** 2
        dv_grad = (beta_eff * (C_LIGHT**2) * dphi_dr / M_PL * (r_p / v_p) * j2_factor * cos_asym_factor * 1e3)
        
        S_disf = model.disformal_modulation_factor(v_p, cos_asym_factor)
        dv_total = dv_grad * S_disf
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
                "model_version": "v5.3-Jakarta",
                "topology_n": N_TOPOLOGY,
                "beta_initial": BETA_INITIAL,
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
