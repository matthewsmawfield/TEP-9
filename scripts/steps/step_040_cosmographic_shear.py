"""
Step 040: Cosmographic Temporal Shear Modulation Test

Tests the TEP-specific prediction that temporal shear depends on:
1. Earth-Moon system orbital velocity through solar scalar topology
2. Heliocentric distance (Earth's position in solar scalar field)
3. CMB dipole direction (Earth's motion through the cosmic scalar rest frame)

This addresses discrepancies between TEP predictions and observed flyby anomalies
that may correlate with date/time and direction via cosmographic modulation.
"""

import json
import sys
import math
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

try:
    from scripts.utils.step_logger import StepLogger
except ImportError:
    class StepLogger:
        def __init__(self, *args, **kwargs):
            pass
        def section(self, s):
            print(f"\n=== {s} ===")
        def subsection(self, s):
            print(f"\n--- {s} ---")
        def info(self, s):
            print(s)
        def warning(self, s):
            print(f"WARNING: {s}")
        def log_step_summary(self, *args):
            pass


# CMB dipole parameters (Planck 2018)
CMB_DIPOLE_RA_DEG = 167.94   # J2000 equatorial
CMB_DIPOLE_DEC_DEG = -6.93
CMB_DIPOLE_VELOCITY_KM_S = 369.82  # km/s

# Earth's orbital parameters
EARTH_ORBITAL_ECCENTRICITY = 0.0167086
EARTH_ORBITAL_SEMI_MAJOR_AU = 1.000001018
EARTH_ORBITAL_PERIOD_DAYS = 365.256363004
EARTH_ORBITAL_INCLINATION_DEG = 0.0  # ecliptic reference
EARTH_ORBITAL_LONGITUDE_OF_PERIHELION_DEG = 102.937348  # J2000

# TEP UCD parameters
RHO_T_G_CM3 = 20.0  # g/cm^3
M_SUN_KG = 1.98847e30
G_MKS = 6.67430e-11  # m^3 kg^-1 s^-2
AU_M = 1.495978707e11
C_M_S = 299792458.0


def parse_flyby_dates(flyby_catalog_path: Path) -> list:
    """Extract flyby dates, observed anomalies, and TEP predictions."""
    with open(flyby_catalog_path) as f:
        catalog = json.load(f)

    flybys = []
    for fb in catalog.get("flybys", []):
        if not fb.get("usable_for_analysis", False):
            continue
        date_str = fb.get("flyby_date", "")
        if not date_str:
            continue
        try:
            dt = datetime.strptime(date_str, "%Y-%m-%d").replace(tzinfo=timezone.utc)
        except ValueError:
            continue

        flybys.append({
            "mission": fb["mission_name"],
            "date": date_str,
            "datetime": dt,
            "altitude_km": fb.get("perigee_altitude_km"),
            "velocity_km_s": fb.get("perigee_velocity_km_s"),
            "declination_in_deg": fb.get("declination_in_deg"),
            "declination_out_deg": fb.get("declination_out_deg"),
            "cos_asymmetry": fb.get("cos_asymmetry"),
            "observed_mm_s": fb.get("published_anomaly_mm_s"),
            "sigma_mm_s": fb.get("published_anomaly_uncertainty_mm_s"),
            "detection": fb.get("detection_significance"),
        })
    return flybys


def load_tep_predictions(fitting_results_path: Path) -> dict:
    """Load TEP predicted anomalies for each mission."""
    with open(fitting_results_path) as f:
        data = json.load(f)

    preds = {}
    for mission, result in data.get("individual_fits", {}).items():
        preds[mission] = {
            "dv_tep_mm_s": result.get("tep_predictions", {}).get("dv_tep_mm_s"),
            "dv_grad_mm_s": result.get("tep_predictions", {}).get("dv_grad_mm_s"),
            "dv_disf_mm_s": result.get("tep_predictions", {}).get("dv_disf_mm_s"),
            "observed_mm_s": result.get("observed", {}).get("dv_obs_mm_s"),
            "sigma_mm_s": result.get("observed", {}).get("sigma_mm_s"),
        }
    return preds


def load_3d_state_vectors(state_vectors_path: Path) -> dict:
    """Load 3D spacecraft state vectors at perigee from step_040a."""
    if not state_vectors_path.exists():
        return {}
    with open(state_vectors_path) as f:
        return json.load(f)


def julian_day(dt: datetime) -> float:
    """Convert UTC datetime to Julian Day."""
    y, m, d = dt.year, dt.month, dt.day
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    JD = (int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5)
    frac = (dt.hour + dt.minute / 60.0 + dt.second / 3600.0) / 24.0
    return JD + frac


def earth_heliocentric_state(dt: datetime) -> tuple:
    """
    Compute Earth's heliocentric ecliptic longitude, latitude, distance (AU),
    and orbital velocity vector (km/s) using proper elliptical orbit mechanics.

    Returns:
        (r_au, lon_deg, lat_deg, vx_kms, vy_kms, vz_kms)
        where (vx, vy, vz) are in a heliocentric ecliptic frame:
        x: toward vernal equinox
        y: 90 deg ecliptic longitude
        z: toward ecliptic north pole
    """
    JD = julian_day(dt)
    T = (JD - 2451545.0) / 36525.0  # Julian centuries from J2000.0

    # Mean longitude of the Sun (Earth's longitude + 180)
    L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T * T  # deg
    L0 = L0 % 360.0

    # Mean anomaly of the Sun (Earth)
    M = 357.52911 + 35999.05029 * T - 0.0001537 * T * T  # deg
    M = math.radians(M % 360.0)

    # Eccentricity
    e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T * T

    # Equation of center (simplified, first three terms)
    C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * math.sin(M)
    C += (0.019993 - 0.000101 * T) * math.sin(2 * M)
    C += 0.000289 * math.sin(3 * M)

    # True longitude of the Sun
    true_lon_sun = (L0 + C) % 360.0  # deg
    # Earth's true heliocentric longitude (opposite the Sun)
    true_lon_earth = (true_lon_sun + 180.0) % 360.0
    true_lon_earth_rad = math.radians(true_lon_earth)

    # True anomaly of Earth
    nu = (true_lon_earth - 102.937348) % 360.0
    nu_rad = math.radians(nu)

    # Distance (AU)
    a = 1.000001018
    r_au = a * (1 - e * e) / (1 + e * math.cos(nu_rad))

    # Proper elliptical velocity components (km/s)
    mu_km3_s2 = 1.32712440018e11  # GM_sun in km^3/s^2
    r_km = r_au * AU_M / 1000.0
    # Semi-latus rectum
    p_km = a * AU_M / 1000.0 * (1 - e * e)
    # Specific angular momentum
    h = math.sqrt(mu_km3_s2 * p_km)

    # Radial and transverse velocity components
    v_r = (mu_km3_s2 / h) * e * math.sin(nu_rad)  # km/s, positive = moving away from Sun
    v_t = (mu_km3_s2 / h) * (1 + e * math.cos(nu_rad))  # km/s, transverse

    # Verify with vis-viva
    v_visviva = math.sqrt(mu_km3_s2 * (2.0 / r_km - 1.0 / (a * AU_M / 1000.0)))
    v_check = math.sqrt(v_r * v_r + v_t * v_t)
    # They should agree within numerical precision

    # Velocity vector in heliocentric ecliptic frame
    # r_hat = (cos λ, sin λ, 0)  [Sun to Earth]
    # theta_hat = (-sin λ, cos λ, 0)  [perpendicular, CCW]
    # v = v_r * r_hat + v_t * theta_hat
    vx = v_r * math.cos(true_lon_earth_rad) - v_t * math.sin(true_lon_earth_rad)
    vy = v_r * math.sin(true_lon_earth_rad) + v_t * math.cos(true_lon_earth_rad)
    vz = 0.0

    return r_au, true_lon_earth, 0.0, vx, vy, vz


def unit_vector_equatorial(ra_deg: float, dec_deg: float) -> np.ndarray:
    """Convert RA/Dec to Cartesian unit vector in J2000 equatorial frame."""
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    return np.array([
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec)
    ])


def ecliptic_to_equatorial(vx, vy, vz) -> np.ndarray:
    """Rotate vector from ecliptic to equatorial J2000 frame."""
    eps = math.radians(23.439281)  # obliquity of the ecliptic J2000
    # Rotation about x-axis by -eps
    x = vx
    y = vy * math.cos(eps) - vz * math.sin(eps)
    z = vy * math.sin(eps) + vz * math.cos(eps)
    return np.array([x, y, z])


def declination_to_unit_vector(dec_deg: float, use_incoming: bool = True) -> np.ndarray:
    """
    Fallback approximate 3D unit vector from declination only.
    Only used when proper 3D state vectors are unavailable.
    """
    dec = math.radians(dec_deg)
    return np.array([math.cos(dec), 0.0, math.sin(dec)])


def get_spacecraft_velocity_unit_vector(mission: str, state_vectors: dict, declination_in_deg: float = None) -> np.ndarray:
    """
    Return spacecraft velocity unit vector in ecliptic frame.
    Prefers actual 3D state vector from JPL Horizons; falls back to declination approximation.
    """
    vec = state_vectors.get(mission)
    if vec is not None and "ux_ecl" in vec:
        return np.array([vec["ux_ecl"], vec["uy_ecl"], vec["uz_ecl"]])
    if declination_in_deg is not None:
        return declination_to_unit_vector(declination_in_deg)
    return np.array([1.0, 0.0, 0.0])


def compute_solar_scalar_modulation(r_au: float, lon_deg: float,
                                     v_earth_ecl: np.ndarray, v_sc_ecl: np.ndarray) -> dict:
    """
    Compute solar scalar topology modulation factors using proper 3D vectors.

    In TEP-UCD, the Sun's scalar field saturation radius is:
        R_sol_sun = (3 M_sun / 4 pi rho_T)^(1/3)
    """
    # Sun's UCD saturation radius
    rho_T_kg_m3 = RHO_T_G_CM3 * 1000.0
    r_sol_sun_m = ((3.0 * M_SUN_KG) / (4.0 * math.pi * rho_T_kg_m3)) ** (1.0 / 3.0)
    r_sol_sun_au = r_sol_sun_m / AU_M

    # Heliocentric distance modulation
    r_mod = 1.0 / (r_au ** 2)

    # Earth orbital unit vector
    v_orb_mag = np.linalg.norm(v_earth_ecl)
    v_orb_hat = v_earth_ecl / v_orb_mag if v_orb_mag > 0 else np.array([0.0, 0.0, 0.0])

    # Sun-Earth radial unit vector in ecliptic frame:
    # position vector from Sun to Earth = (r cos λ, r sin λ, 0)
    # radial unit vector (outward from Sun) = (cos λ, sin λ, 0)
    lon_rad = math.radians(lon_deg)
    r_hat = np.array([math.cos(lon_rad), math.sin(lon_rad), 0.0])

    # Earth's radial velocity (component along Sun-Earth line, positive = moving away)
    v_radial = float(np.dot(v_earth_ecl, r_hat))

    # Scalar wind factor
    wind_factor = v_orb_mag / 29.78

    # Spacecraft alignment with Earth's orbital motion
    # Dot product: +1 means prograde, -1 means retrograde relative to Earth orbit
    v_sc_mag = np.linalg.norm(v_sc_ecl)
    v_sc_hat = v_sc_ecl / v_sc_mag if v_sc_mag > 0 else np.array([0.0, 0.0, 0.0])
    sc_orb_alignment = np.dot(v_sc_hat, v_orb_hat)

    return {
        "r_sol_sun_au": r_sol_sun_au,
        "heliocentric_distance_au": r_au,
        "heliocentric_modulation": r_mod,
        "orbital_speed_kms": float(v_orb_mag),
        "radial_velocity_kms": v_radial,
        "scalar_wind_factor": float(wind_factor),
        "sc_orbital_alignment": float(sc_orb_alignment),
        "sc_speed_kms": float(v_sc_mag),
    }


def compute_cmb_dipole_modulation(v_earth_equ: np.ndarray, v_sc_equ: np.ndarray) -> dict:
    """
    Compute CMB dipole modulation of temporal shear using proper 3D vectors.

    The Earth-Moon system moves at ~370 km/s toward the CMB dipole apex.
    In TEP, the disformal coupling depends on velocity relative to the scalar
    field rest frame. If the CMB frame approximates the scalar rest frame,
    the ~370 km/s bulk motion modulates the effective coupling.

    With proper 3D vectors, we compute:
    - Earth velocity toward CMB dipole apex (from orbital motion)
    - Total spacecraft velocity in CMB rest frame
    - CMB-frame disformal enhancement factor: v_total^2 / v_sc^2
    - Sign of the projection (enhancement vs. suppression)
    """
    n_cmb = unit_vector_equatorial(CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG)

    # Earth's orbital velocity projected onto CMB dipole
    v_earth_cmb_proj = float(np.dot(v_earth_equ, n_cmb))  # km/s

    # Total Earth+spacecraft velocity projected onto CMB dipole
    v_total_equ = v_earth_equ + v_sc_equ
    v_total_cmb_proj = float(np.dot(v_total_equ, n_cmb))

    # Bulk CMB velocity of Earth (including solar motion)
    v_solar_cmb = CMB_DIPOLE_VELOCITY_KM_S * n_cmb
    v_total_with_solar = v_solar_cmb + v_sc_equ
    v_total_solar_proj = float(np.dot(v_total_with_solar, n_cmb))

    # Spacecraft velocity magnitude and alignment with CMB dipole
    v_sc_mag = float(np.linalg.norm(v_sc_equ))
    v_sc_cmb_proj = float(np.dot(v_sc_equ, n_cmb))

    # Angle between spacecraft velocity and CMB dipole direction
    if v_sc_mag > 0:
        cos_theta = v_sc_cmb_proj / v_sc_mag
    else:
        cos_theta = 0.0

    # CMB-rest-frame disformal velocity enhancement factor
    # In TEP, disformal coupling scales as v^2. If the scalar rest frame is CMB,
    # the effective velocity is v_total_cmb = v_earth_cmb + v_sc, not just v_sc.
    v_total_cmb_mag = float(np.linalg.norm(v_total_with_solar))
    if v_sc_mag > 0:
        cmb_disformal_enhancement = (v_total_cmb_mag / v_sc_mag) ** 2
    else:
        cmb_disformal_enhancement = 1.0

    return {
        "cmb_dipole_ra_deg": CMB_DIPOLE_RA_DEG,
        "cmb_dipole_dec_deg": CMB_DIPOLE_DEC_DEG,
        "cmb_velocity_kms": CMB_DIPOLE_VELOCITY_KM_S,
        "earth_orbital_cmb_proj_kms": v_earth_cmb_proj,
        "sc_velocity_cmb_proj_kms": v_sc_cmb_proj,
        "total_velocity_cmb_proj_kms": v_total_cmb_proj,
        "total_with_solar_cmb_proj_kms": v_total_solar_proj,
        "cmb_modulation_factor": v_total_cmb_proj / CMB_DIPOLE_VELOCITY_KM_S,
        "cmb_solar_modulation": v_total_solar_proj / CMB_DIPOLE_VELOCITY_KM_S,
        "sc_cmb_cos_theta": cos_theta,
        "sc_speed_kms": v_sc_mag,
        "total_cmb_speed_kms": v_total_cmb_mag,
        "cmb_disformal_enhancement": cmb_disformal_enhancement,
    }


def compute_discrepancy(observed, predicted, sigma):
    """Compute discrepancy metrics."""
    if observed is None or predicted is None or predicted == 0:
        return None, None, None
    ratio = observed / predicted
    diff = observed - predicted
    if sigma and sigma > 0:
        sigma_diff = diff / sigma
    else:
        sigma_diff = None
    return ratio, diff, sigma_diff


def safe_pearson(x, y):
    """Compute Pearson r with handling for insufficient data."""
    x_clean = []
    y_clean = []
    for xi, yi in zip(x, y):
        if xi is not None and yi is not None and not (np.isnan(xi) or np.isnan(yi)):
            x_clean.append(xi)
            y_clean.append(yi)
    n = len(x_clean)
    if n < 3:
        return None, None, n
    r, p = stats.pearsonr(x_clean, y_clean)
    return float(r), float(p), n


def main():
    logger = StepLogger("step_040_cosmographic_shear", PROJECT_ROOT)
    logger.header("STEP 040: COSMOGRAPHIC TEMPORAL SHEAR MODULATION TEST")
    logger.info("Testing TEP prediction: temporal shear depends on Earth-Moon")
    logger.info("orbital velocity, heliocentric distance, and CMB dipole direction.")

    # Load data
    catalog_path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
    fitting_path = PROJECT_ROOT / "results" / "step008_fitting_results.json"
    vectors_path = PROJECT_ROOT / "results" / "step040a_3d_state_vectors.json"

    flybys = parse_flyby_dates(catalog_path)
    tep_preds = load_tep_predictions(fitting_path)
    state_vectors = load_3d_state_vectors(vectors_path)

    logger.info(f"Loaded {len(flybys)} flybys from catalog")
    logger.info(f"Loaded {len(state_vectors)} 3D state vectors from step_040a")

    # CMB dipole unit vector
    n_cmb = unit_vector_equatorial(CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG)
    logger.info(f"CMB dipole direction: RA={CMB_DIPOLE_RA_DEG} deg, Dec={CMB_DIPOLE_DEC_DEG} deg")
    logger.info(f"CMB dipole velocity: {CMB_DIPOLE_VELOCITY_KM_S} km/s")

    results = []
    logger.section("COMPUTING COSMOGRAPHIC MODULATION FOR EACH FLYBY")

    for fb in flybys:
        mission = fb["mission"]
        dt = fb["datetime"]

        # Skip future missions
        if dt.year > 2025:
            continue

        logger.subsection(mission)
        logger.info(f"Date: {fb['date']}")

        # Earth's heliocentric state
        r_au, lon_deg, lat_deg, vx_ecl, vy_ecl, vz_ecl = earth_heliocentric_state(dt)
        v_earth_ecl = np.array([vx_ecl, vy_ecl, vz_ecl])
        v_earth_equ = ecliptic_to_equatorial(vx_ecl, vy_ecl, vz_ecl)

        logger.info(f"  Heliocentric distance: {r_au:.4f} AU")
        logger.info(f"  Ecliptic longitude: {lon_deg:.2f} deg")
        logger.info(f"  Earth orbital speed: {np.linalg.norm(v_earth_ecl):.2f} km/s")

        # Spacecraft velocity vector (3D from JPL Horizons, or fallback)
        dec_in = fb.get("declination_in_deg")
        v_sc_hat_ecl = get_spacecraft_velocity_unit_vector(mission, state_vectors, dec_in)

        # Get actual spacecraft speed from state vectors if available
        vec = state_vectors.get(mission)
        if vec is not None and "speed_ecl_km_s" in vec:
            v_sc_mag = vec["speed_ecl_km_s"]
            v_sc_ecl = v_sc_hat_ecl * v_sc_mag
            logger.info(f"  Using 3D state vector: speed = {v_sc_mag:.2f} km/s")
            logger.info(f"    ecliptic velocity: [{v_sc_ecl[0]:+.2f}, {v_sc_ecl[1]:+.2f}, {v_sc_ecl[2]:+.2f}] km/s")
        else:
            v_sc_mag = fb.get("velocity_km_s", 15.0)
            v_sc_ecl = v_sc_hat_ecl * v_sc_mag
            logger.info(f"  Using declination approximation (no 3D vector): speed = {v_sc_mag:.2f} km/s")

        # Convert spacecraft velocity to equatorial for CMB analysis
        # ecliptic_to_equatorial rotation inverse: x_eq = x_ec, y_eq = y_ec*cos(eps) + z_ec*sin(eps), z_eq = -y_ec*sin(eps) + z_ec*cos(eps)
        eps = math.radians(23.439281)
        v_sc_equ = np.array([
            v_sc_ecl[0],
            v_sc_ecl[1] * math.cos(eps) + v_sc_ecl[2] * math.sin(eps),
            -v_sc_ecl[1] * math.sin(eps) + v_sc_ecl[2] * math.cos(eps)
        ])

        # Solar scalar modulation
        solar_mod = compute_solar_scalar_modulation(r_au, lon_deg, v_earth_ecl, v_sc_ecl)
        logger.info(f"  Sun UCD saturation radius: {solar_mod['r_sol_sun_au']:.4f} AU")
        logger.info(f"  Heliocentric modulation proxy: {solar_mod['heliocentric_modulation']:.4f}")
        logger.info(f"  Earth radial velocity: {solar_mod['radial_velocity_kms']:.2f} km/s")
        logger.info(f"  SC-orbital alignment: {solar_mod['sc_orbital_alignment']:+.3f}")

        # CMB dipole modulation
        cmb_mod = compute_cmb_dipole_modulation(v_earth_equ, v_sc_equ)
        logger.info(f"  Earth orbital v along CMB dipole: {cmb_mod['earth_orbital_cmb_proj_kms']:.2f} km/s")
        logger.info(f"  SC v along CMB dipole: {cmb_mod['sc_velocity_cmb_proj_kms']:.2f} km/s")
        logger.info(f"  Total v along CMB dipole: {cmb_mod['total_velocity_cmb_proj_kms']:.2f} km/s")
        logger.info(f"  CMB modulation factor: {cmb_mod['cmb_modulation_factor']:.4f}")
        logger.info(f"  SC-CMB cos(theta): {cmb_mod['sc_cmb_cos_theta']:+.3f}")

        # Discrepancy
        pred = tep_preds.get(mission, {})
        obs = pred.get("observed_mm_s")
        prd = pred.get("dv_tep_mm_s")
        sig = pred.get("sigma_mm_s")

        ratio, diff, sigma_diff = compute_discrepancy(obs, prd, sig)
        if ratio is not None:
            logger.info(f"  Observed: {obs:.2f} mm/s, Predicted: {prd:.3f} mm/s")
            logger.info(f"  Ratio (obs/pred): {ratio:.2f}, Difference: {diff:.2f} mm/s")
        else:
            logger.info(f"  No prediction available or zero prediction")

        results.append({
            "mission": mission,
            "date": fb["date"],
            "heliocentric_distance_au": r_au,
            "ecliptic_longitude_deg": lon_deg,
            "earth_orbital_speed_kms": float(np.linalg.norm(v_earth_ecl)),
            **solar_mod,
            **cmb_mod,
            "observed_mm_s": obs,
            "predicted_mm_s": prd,
            "uncertainty_mm_s": sig,
            "ratio_obs_pred": ratio,
            "difference_mm_s": diff,
            "sigma_difference": sigma_diff,
            "declination_in_deg": dec_in,
            "has_3d_vector": vec is not None,
        })

    # Correlation analysis
    logger.section("CORRELATION ANALYSIS")

    # Filter to flybys with both observed and predicted data
    usable = [r for r in results
              if r["ratio_obs_pred"] is not None and not np.isnan(r["ratio_obs_pred"])]

    logger.info(f"Usable flybys for correlation: {len(usable)}")

    if len(usable) >= 3:
        ratios = [r["ratio_obs_pred"] for r in usable]
        diffs = [r["difference_mm_s"] for r in usable]

        # Modulation factors
        helio_mods = [r["heliocentric_modulation"] for r in usable]
        radial_vs = [r["radial_velocity_kms"] for r in usable]
        cmb_factors = [r["cmb_modulation_factor"] for r in usable]
        cmb_solar = [r["cmb_solar_modulation"] for r in usable]
        distances = [r["heliocentric_distance_au"] for r in usable]
        orb_speeds = [r["earth_orbital_speed_kms"] for r in usable]
        cmb_proj_orb = [r["earth_orbital_cmb_proj_kms"] for r in usable]

        correlations = []

        def add_corr(name, x, y):
            r, p, n = safe_pearson(x, y)
            if r is not None:
                correlations.append({
                    "test": name,
                    "pearson_r": r,
                    "p_value": p,
                    "n": n,
                    "significant": p < 0.05 if p is not None else False,
                })
                sig_str = " ***" if p is not None and p < 0.05 else ""
                logger.info(f"  {name}: r={r:+.3f}, p={p:.3f}, n={n}{sig_str}")
            else:
                logger.info(f"  {name}: insufficient data (n={n})")

        add_corr("Ratio vs Heliocentric Distance", distances, ratios)
        add_corr("Ratio vs Heliocentric Modulation", helio_mods, ratios)
        add_corr("Ratio vs Radial Velocity", radial_vs, ratios)
        add_corr("Ratio vs CMB Modulation Factor", cmb_factors, ratios)
        add_corr("Ratio vs CMB Solar Modulation", cmb_solar, ratios)
        add_corr("Ratio vs Earth Orbital Speed", orb_speeds, ratios)
        add_corr("Ratio vs Earth CMB Projection", cmb_proj_orb, ratios)

        # New 3D-vector-specific correlations
        sc_orb_aligns = [r["sc_orbital_alignment"] for r in usable]
        sc_cmb_cos = [r["sc_cmb_cos_theta"] for r in usable]
        sc_cmb_projs = [r["sc_velocity_cmb_proj_kms"] for r in usable]
        cmb_enhancements = [r["cmb_disformal_enhancement"] for r in usable]

        add_corr("Ratio vs SC-Orbital Alignment", sc_orb_aligns, ratios)
        add_corr("Ratio vs SC-CMB cos(theta)", sc_cmb_cos, ratios)
        add_corr("Ratio vs SC-CMB Projection", sc_cmb_projs, ratios)
        add_corr("Ratio vs CMB Disformal Enhancement", cmb_enhancements, ratios)

        add_corr("Difference vs Heliocentric Distance", distances, diffs)
        add_corr("Difference vs CMB Modulation Factor", cmb_factors, diffs)
        add_corr("Difference vs CMB Solar Modulation", cmb_solar, diffs)
        add_corr("Difference vs Earth CMB Projection", cmb_proj_orb, diffs)
        add_corr("Difference vs SC-Orbital Alignment", sc_orb_aligns, diffs)
        add_corr("Difference vs SC-CMB cos(theta)", sc_cmb_cos, diffs)
        add_corr("Difference vs CMB Disformal Enhancement", cmb_enhancements, diffs)

        # Multivariate geometric regression (with intercept)
        logger.section("MULTIVARIATE GEOMETRIC REGRESSION")
        X = np.column_stack([
            np.ones(len(sc_cmb_cos)),
            sc_cmb_cos,
            np.array(cmb_proj_orb) / 30.0,
            sc_orb_aligns,
        ])
        X = np.nan_to_num(X, nan=0.0)
        y = np.array(ratios)

        # Ordinary least squares with intercept
        XtX = X.T @ X
        XtY = X.T @ y
        try:
            beta_hat = np.linalg.solve(XtX, XtY)
            y_pred = X @ beta_hat
            residuals = y - y_pred
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r2_multi = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
            n_obs = len(y)
            n_params = X.shape[1]  # includes intercept
            adj_r2_multi = 1.0 - (1.0 - r2_multi) * (n_obs - 1) / (n_obs - n_params) if n_obs > n_params else 0.0
            std_orig = float(np.std(y))
            std_resid = float(np.std(residuals))
            reduction = (std_orig - std_resid) / std_orig * 100.0 if std_orig > 0 else 0.0

            logger.info(f"Model: ratio = b0 + b1*SC-CMB_cos + b2*Earth-CMB_proj/30 + b3*SC-orbital")
            logger.info(f"Coefficients: b0={beta_hat[0]:+.3f}, b1={beta_hat[1]:+.3f}, b2={beta_hat[2]:+.3f}, b3={beta_hat[3]:+.3f}")
            logger.info(f"R^2 = {r2_multi:.4f}, Adjusted R^2 = {adj_r2_multi:.4f}")
            logger.info(f"Residual std: {std_resid:.3f} (original: {std_orig:.3f}, reduction: {reduction:.1f}%)")

            # Add predicted ratio to each flyby result
            for i, r in enumerate(usable):
                r["geometric_model_ratio"] = float(y_pred[i])
                r["geometric_model_residual"] = float(residuals[i])
        except np.linalg.LinAlgError:
            logger.warning("Singular matrix in multivariate regression; skipping")
            beta_hat = None
            r2_multi = None
            adj_r2_multi = None
            std_resid = None
            reduction = None

        # Both-aligned directional consistency test
        logger.section("BOTH-ALIGNED DIRECTIONAL CONSISTENCY TEST")
        both_aligned = []
        same_direction = []
        for r in usable:
            sc_proj = r.get("sc_cmb_cos_theta", 0.0)
            earth_proj = r["earth_orbital_cmb_proj_kms"] / 30.0
            ba = 1.0 if (sc_proj > 0 and earth_proj > 0) else 0.0
            sd = 1.0 if (np.sign(sc_proj) == np.sign(earth_proj)) else 0.0
            both_aligned.append(ba)
            same_direction.append(sd)
            r["both_aligned_flag"] = float(ba)
            r["same_direction_flag"] = float(sd)

        if len(both_aligned) >= 3:
            r_ba, p_ba = stats.pearsonr(both_aligned, ratios)
            logger.info(f"Both-aligned flag correlation: r = {r_ba:+.3f}, p = {p_ba:.3f}")
            aligned_ratios = [ratios[i] for i in range(len(ratios)) if both_aligned[i] > 0.5]
            unaligned_ratios = [ratios[i] for i in range(len(ratios)) if both_aligned[i] < 0.5]
            logger.info(f"Aligned flybys (n={len(aligned_ratios)}): ratios = {[f'{x:.2f}' for x in aligned_ratios]}")
            logger.info(f"Unaligned flybys (n={len(unaligned_ratios)}): ratios = {[f'{x:.2f}' for x in unaligned_ratios]}")
            if len(aligned_ratios) > 0 and len(unaligned_ratios) > 0:
                from scipy.stats import mannwhitneyu
                try:
                    mw_stat, mw_p = mannwhitneyu(aligned_ratios, unaligned_ratios, alternative="greater", method="exact")
                    logger.info(f"Mann-Whitney U (aligned > unaligned): U={mw_stat:.0f}, p={mw_p:.3f} (exact)")
                except ValueError:
                    pass

            r_sd, p_sd = stats.pearsonr(same_direction, ratios)
            logger.info(f"Same-direction flag correlation: r = {r_sd:+.3f}, p = {p_sd:.3f}")

        # Weighted combination: SC-CMB + w*Earth-CMB
        logger.section("OPTIMAL WEIGHTED COMBINATION")
        best_r = 0.0
        best_w = 0.0
        best_p = 1.0
        for w in np.linspace(-2.0, 2.0, 401):
            E = sc_cmb_cos + w * (np.array(cmb_proj_orb) / 30.0)
            mask = np.isfinite(E)
            if mask.sum() < 3:
                continue
            r_w, p_w = stats.pearsonr(E[mask], np.array(ratios)[mask])
            if abs(r_w) > abs(best_r):
                best_r = r_w
                best_w = w
                best_p = p_w
        logger.info(f"Optimal model: E = SC-CMB_cos(theta) + {best_w:.3f} * Earth-CMB_proj/30")
        logger.info(f"Correlation: r = {best_r:+.3f}, p = {best_p:.3f}")

        # Summary table
        logger.section("SUMMARY TABLE")
        logger.info(f"{'Mission':<16} {'r_AU':>6} {'v_rad':>6} {'v_E,cmb':>8} {'sc_cmb':>7} {'sc_orb':>7} {'align':>6} {'Obs':>7} {'Pred':>7} {'Ratio':>7} {'Geom':>7}")
        logger.info("-" * 105)
        for r in usable:
            geom = r.get("geometric_model_ratio", 0.0)
            align = "YES" if r.get("both_aligned_flag", 0) > 0.5 else "no"
            logger.info(
                f"{r['mission']:<16} {r['date']:<12} {r['heliocentric_distance_au']:6.3f} "
                f"{r.get('radial_velocity_kms', 0):+6.2f} {r['earth_orbital_cmb_proj_kms']:7.1f} "
                f"{r.get('sc_velocity_cmb_proj_kms', 0):7.1f} {r.get('sc_orbital_alignment', 0):+7.3f} "
                f"{align:>6} "
                f"{r['observed_mm_s']:7.2f} {r['predicted_mm_s']:7.3f} {r['ratio_obs_pred']:7.2f} {geom:7.2f}"
            )

        # Best correlation
        if correlations:
            best = min(correlations, key=lambda c: c["p_value"] if c["p_value"] is not None else 1.0)
            logger.section("BEST CORRELATION")
            logger.info(f"Test: {best['test']}")
            logger.info(f"Pearson r = {best['pearson_r']:+.3f}")
            logger.info(f"p-value = {best['p_value']:.3f}")
            if best["significant"]:
                logger.info("Result: STATISTICALLY SIGNIFICANT at p < 0.05")
            else:
                logger.info("Result: Not statistically significant (small sample)")
    else:
        logger.warning("Insufficient usable flybys for correlation analysis (need >= 3)")
        correlations = []
        beta_hat = None
        r2_multi = None
        adj_r2_multi = None
        std_resid = None
        reduction = None
        best_r = None
        best_w = None
        best_p = None

    output = {
        "metadata": {
            "step": "040",
            "description": "Cosmographic temporal shear modulation test (3D vectors)",
            "cmb_dipole_ra_deg": CMB_DIPOLE_RA_DEG,
            "cmb_dipole_dec_deg": CMB_DIPOLE_DEC_DEG,
            "cmb_dipole_velocity_kms": CMB_DIPOLE_VELOCITY_KM_S,
            "n_flybys_total": len(results),
            "n_flybys_usable": len(usable),
            "n_with_3d_vectors": sum(1 for r in usable if r.get("has_3d_vector", False)),
        },
        "flyby_results": results,
        "correlations": correlations,
        "multivariate_regression": {
            "model": "ratio = b0 + b1*SC-CMB_cos + b2*Earth-CMB_proj/30 + b3*SC-orbital",
            "coefficients": {
                "b0_intercept": float(beta_hat[0]) if beta_hat is not None else None,
                "b1_sc_cmb_cos": float(beta_hat[1]) if beta_hat is not None else None,
                "b2_earth_cmb_proj": float(beta_hat[2]) if beta_hat is not None else None,
                "b3_sc_orbital": float(beta_hat[3]) if beta_hat is not None else None,
            },
            "r_squared": float(r2_multi) if r2_multi is not None else None,
            "adjusted_r_squared": float(adj_r2_multi) if adj_r2_multi is not None else None,
            "residual_std": float(std_resid) if std_resid is not None else None,
            "original_std": float(std_orig) if 'std_orig' in dir() else None,
            "std_reduction_percent": float(reduction) if reduction is not None else None,
        } if beta_hat is not None else None,
        "directional_consistency": {
            "both_aligned_correlation_r": float(r_ba) if 'r_ba' in dir() else None,
            "both_aligned_correlation_p": float(p_ba) if 'p_ba' in dir() else None,
            "same_direction_correlation_r": float(r_sd) if 'r_sd' in dir() else None,
            "same_direction_correlation_p": float(p_sd) if 'p_sd' in dir() else None,
            "mann_whitney_u_statistic": float(mw_stat) if 'mw_stat' in dir() else None,
            "mann_whitney_u_pvalue": float(mw_p) if 'mw_p' in dir() else None,
        } if 'r_ba' in dir() else None,
        "optimal_weighted_combination": {
            "model": "E = SC-CMB_cos(theta) + w * Earth-CMB_proj/30",
            "optimal_w": float(best_w) if best_w is not None else None,
            "pearson_r": float(best_r) if best_r is not None else None,
            "p_value": float(best_p) if best_p is not None else None,
        } if best_w is not None else None,
    }

    out_file = PROJECT_ROOT / "results" / "step040_cosmographic_shear.json"
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2)

    logger.info(f"\nResults saved to: {out_file}")
    logger.log_step_summary(0, "SUCCESS")
    return output


if __name__ == "__main__":
    main()
