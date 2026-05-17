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
import os
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
    """Load 3D spacecraft state vectors at perigee from step_038."""
    if not state_vectors_path.exists():
        return {}
    with open(state_vectors_path) as f:
        return json.load(f)


def icrs_cartesian_to_ecliptic_cartesian(x: float, y: float, z: float) -> np.ndarray:
    """
    Rotate ICRS / J2000 equatorial Cartesian components into ecliptic coordinates
    using the same obliquity convention as `step_038_extract_3d_vectors.equatorial_to_ecliptic`.
    """
    eps = math.radians(23.439281)
    xe = x
    ye = y * math.cos(eps) + z * math.sin(eps)
    ze = -y * math.sin(eps) + z * math.cos(eps)
    return np.array([xe, ye, ze], dtype=float)


def earth_sun_ephemeris_state(dt: datetime, logger, quiet: bool = False) -> dict:
    """
    Earth and Sun barycentric state in ICRS via Astropy, aligned with Horizons-grade
    inner-planet ephemeris (default: `builtin`; override with TEP_SOLAR_SYSTEM_EPHEMERIS).

    Returns heliocentric distance and ecliptic longitude from (r_Earth − r_Sun),
    heliocentric Earth velocity (v_Earth − v_Sun) in ecliptic km/s for solar modulation,
    Earth barycentric velocity in ICRS km/s for adding to geocentric spacecraft
    velocity from Horizons (consistent inertial chain), and Earth/Sun barycentric
    positions in ICRS (AU) for trajectory mapping in Step 042.
    """
    from astropy.time import Time
    from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
    import astropy.units as u

    t = Time(dt, scale="utc")
    ephem = os.environ.get("TEP_SOLAR_SYSTEM_EPHEMERIS", "de440").strip()
    if not ephem:
        ephem = "de440"

    try:
        with solar_system_ephemeris.set(ephem):
            p_e, v_e = get_body_barycentric_posvel("earth", t)
            p_s, v_s = get_body_barycentric_posvel("sun", t)
    except Exception as exc:
        raise RuntimeError(
            f"Step 040 Earth/Sun ephemeris failed (TEP_SOLAR_SYSTEM_EPHEMERIS={ephem!r}). "
            "Use 'builtin' (ships with Astropy) or another Astropy-supported ephemeris name."
        ) from exc

    pos_e = np.asarray(p_e.get_xyz().to(u.au).value, dtype=float).reshape(3)
    pos_s = np.asarray(p_s.get_xyz().to(u.au).value, dtype=float).reshape(3)
    vel_e = np.asarray(v_e.get_xyz().to(u.km / u.s).value, dtype=float).reshape(3)
    vel_s = np.asarray(v_s.get_xyz().to(u.km / u.s).value, dtype=float).reshape(3)

    r_rel_icrs = pos_e - pos_s
    r_rel_ecl = icrs_cartesian_to_ecliptic_cartesian(*r_rel_icrs)
    r_au = float(np.linalg.norm(r_rel_ecl))
    lon_deg = math.degrees(math.atan2(r_rel_ecl[1], r_rel_ecl[0])) % 360.0

    v_rel_icrs = vel_e - vel_s
    v_rel_ecl = icrs_cartesian_to_ecliptic_cartesian(*v_rel_icrs)

    if not quiet:
        logger.info(f"  Earth–Sun ephemeris ({ephem}): r = {r_au:.4f} AU, λ_ecl = {lon_deg:.2f} deg")
        logger.info(f"  Earth barycentric speed (ICRS): {float(np.linalg.norm(vel_e)):.2f} km/s")
        logger.info(f"  Earth helio-relative speed (ecliptic): {float(np.linalg.norm(v_rel_ecl)):.2f} km/s")

    return {
        "ephemeris": ephem,
        "r_au": r_au,
        "ecliptic_longitude_deg": lon_deg,
        "earth_barycentric_vel_kms_icrs": vel_e,
        "earth_helio_vel_kms_ecliptic": v_rel_ecl,
        "sun_barycentric_vel_kms_icrs": vel_s,
        "earth_barycentric_pos_au_icrs": pos_e.tolist(),
        "sun_barycentric_pos_au_icrs": pos_s.tolist(),
    }


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


def get_spacecraft_velocity_unit_vector(mission: str, state_vectors: dict) -> np.ndarray:
    """Return spacecraft velocity unit vector in ecliptic frame from step_038."""
    vec = state_vectors.get(mission)
    if vec is None or "ux_ecl" not in vec:
        raise RuntimeError(
            f"Missing JPL Horizons 3D state vector for {mission}; "
            "declination-only geometry is not used"
        )
    return np.array([vec["ux_ecl"], vec["uy_ecl"], vec["uz_ecl"]])


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
    - Earth barycentric velocity projected onto the CMB dipole apex (ICRS, km/s;
      added to Horizons geocentric spacecraft velocity for an inertial sum)
    - Total Earth+spacecraft velocity projected onto the CMB dipole (without the
      fixed Planck bulk template term; see bulk row in returned dict)
    - CMB-frame disformal enhancement factor: v_total^2 / v_sc^2
    - Sign of the projection (enhancement vs. suppression)
    """
    n_cmb = unit_vector_equatorial(CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG)

    # Earth's velocity projected onto CMB dipole (field name retained: historically
    # analytic orbital velocity; now barycentric ICRS from Astropy when Step 040 runs).
    v_earth_cmb_proj = float(np.dot(v_earth_equ, n_cmb))  # km/s

    # Total Earth+spacecraft velocity projected onto CMB dipole
    v_total_equ = v_earth_equ + v_sc_equ
    v_total_cmb_proj = float(np.dot(v_total_equ, n_cmb))

    # Bulk CMB velocity of Earth (including solar motion)
    v_solar_cmb = CMB_DIPOLE_VELOCITY_KM_S * n_cmb
    v_total_with_solar = v_solar_cmb + v_earth_equ + v_sc_equ
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


def _extract_geometry_arrays(rows: list) -> dict:
    """Modulation factors aligned with `rows` order (for correlation y-vector)."""
    return {
        "helio_mods": [r["heliocentric_modulation"] for r in rows],
        "radial_vs": [r["radial_velocity_kms"] for r in rows],
        "cmb_factors": [r["cmb_modulation_factor"] for r in rows],
        "cmb_solar": [r["cmb_solar_modulation"] for r in rows],
        "distances": [r["heliocentric_distance_au"] for r in rows],
        "orb_speeds": [r["earth_orbital_speed_kms"] for r in rows],
        "cmb_proj_orb": [r["earth_orbital_cmb_proj_kms"] for r in rows],
        "sc_orb_aligns": [r["sc_orbital_alignment"] for r in rows],
        "sc_cmb_cos": [r["sc_cmb_cos_theta"] for r in rows],
        "sc_cmb_projs": [r["sc_velocity_cmb_proj_kms"] for r in rows],
        "cmb_enhancements": [r["cmb_disformal_enhancement"] for r in rows],
    }


def ratio_geometry_correlations(rows: list, y_values: list, y_label: str) -> list:
    """
    Pearson correlations between each cosmographic proxy and a per-flyby scalar y
    (same length as rows). y_label prefixes each test name in outputs.
    """
    if len(rows) < 3 or len(rows) != len(y_values):
        return []
    g = _extract_geometry_arrays(rows)
    correlations = []

    def add_corr(name, x, y):
        r, p, n = safe_pearson(x, y)
        if r is not None:
            correlations.append({
                "test": f"{y_label} vs {name}",
                "pearson_r": r,
                "p_value": p,
                "n": n,
                "significant": p < 0.05 if p is not None else False,
            })

    y = y_values
    add_corr("Heliocentric Distance", g["distances"], y)
    add_corr("Heliocentric Modulation", g["helio_mods"], y)
    add_corr("Radial Velocity", g["radial_vs"], y)
    add_corr("CMB Modulation Factor", g["cmb_factors"], y)
    add_corr("CMB Solar Modulation", g["cmb_solar"], y)
    add_corr("Earth Orbital Speed", g["orb_speeds"], y)
    add_corr("Earth CMB Projection", g["cmb_proj_orb"], y)
    add_corr("SC-Orbital Alignment", g["sc_orb_aligns"], y)
    add_corr("SC-CMB cos(theta)", g["sc_cmb_cos"], y)
    add_corr("SC-CMB Projection", g["sc_cmb_projs"], y)
    add_corr("CMB Disformal Enhancement", g["cmb_enhancements"], y)
    return correlations


def max_abs_pearson_weighted_E_vs_y(
    sc_cmb_cos: list,
    cmb_proj_orb: list,
    y: np.ndarray,
    w_min: float = -2.0,
    w_max: float = 2.0,
    n_w: int = 201,
) -> tuple:
    """Scan w in E = cos(theta) + w * (Earth_CMB_proj/30); return (max_abs_r, best_w, r_at_best_w)."""
    sc = np.asarray(sc_cmb_cos, dtype=float)
    ep = np.asarray(cmb_proj_orb, dtype=float)
    yv = np.asarray(y, dtype=float)
    mask = np.isfinite(sc) & np.isfinite(ep) & np.isfinite(yv)
    sc, ep, yv = sc[mask], ep[mask], yv[mask]
    if len(yv) < 3:
        return 0.0, 0.0, 0.0
    base = ep / 30.0
    best_abs = 0.0
    best_w = 0.0
    best_r = 0.0
    for w in np.linspace(w_min, w_max, n_w):
        e_vec = sc + w * base
        if not np.all(np.isfinite(e_vec)):
            continue
        r_w = float(np.corrcoef(e_vec, yv)[0, 1])
        if abs(r_w) > best_abs:
            best_abs = abs(r_w)
            best_w = float(w)
            best_r = r_w
    return best_abs, best_w, best_r


def permutation_pvalue_max_abs_r_weighted(
    sc_cmb_cos: list,
    cmb_proj_orb: list,
    y: np.ndarray,
    n_perm: int = 4000,
    seed: int = 42,
) -> float:
    """
    Permutation null: randomly reassign y across flybys (destroying geometry pairing).
    p = (1 + #{|r_max_perm| >= |r_max_obs|}) / (n_perm + 1).
    """
    rng = np.random.default_rng(seed)
    y_obs = np.asarray(y, dtype=float)
    obs_max, _w0, _r0 = max_abs_pearson_weighted_E_vs_y(sc_cmb_cos, cmb_proj_orb, y_obs)
    exceed = 0
    for _ in range(n_perm):
        yp = rng.permutation(y_obs)
        m, _w, _r = max_abs_pearson_weighted_E_vs_y(sc_cmb_cos, cmb_proj_orb, yp)
        if m >= obs_max - 1e-15:
            exceed += 1
    return (exceed + 1.0) / (n_perm + 1.0)


def main():
    logger = StepLogger("step_040_cosmographic_shear", PROJECT_ROOT)
    logger.header("STEP 040: COSMOGRAPHIC TEMPORAL SHEAR MODULATION TEST")
    logger.info("Testing TEP prediction: temporal shear depends on Earth-Moon")
    logger.info("orbital velocity, heliocentric distance, and CMB dipole direction.")

    # Load data
    catalog_path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
    fitting_path = PROJECT_ROOT / "results" / "step008_fitting_results.json"
    vectors_path = PROJECT_ROOT / "results" / "step038_3d_state_vectors.json"

    flybys = parse_flyby_dates(catalog_path)
    tep_preds = load_tep_predictions(fitting_path)
    state_vectors = load_3d_state_vectors(vectors_path)

    logger.info(f"Loaded {len(flybys)} flybys from catalog")
    logger.info(f"Loaded {len(state_vectors)} 3D state vectors from step_038")

    results = []
    ephemeris_backend = os.environ.get("TEP_SOLAR_SYSTEM_EPHEMERIS", "builtin").strip() or "builtin"
    logger.section("COMPUTING COSMOGRAPHIC MODULATION FOR EACH FLYBY")

    for fb in flybys:
        mission = fb["mission"]
        dt = fb["datetime"]

        # Skip future missions
        if dt.year > 2025:
            continue

        if mission not in state_vectors:
            logger.warning(
                f"Skipping {mission}: no JPL Horizons 3D state vector in step_038 output"
            )
            continue

        logger.subsection(mission)
        logger.info(f"Date: {fb['date']}")

        # Earth–Sun geometry and Earth velocity: Astropy JPL-grade ephemeris (ICRS),
        # consistent with adding Horizons geocentric spacecraft velocity for CMB-frame sums.
        ephem = earth_sun_ephemeris_state(dt, logger)
        ephemeris_backend = str(ephem["ephemeris"])
        r_au = ephem["r_au"]
        lon_deg = ephem["ecliptic_longitude_deg"]
        v_earth_ecl = np.asarray(ephem["earth_helio_vel_kms_ecliptic"], dtype=float)
        v_earth_equ = np.asarray(ephem["earth_barycentric_vel_kms_icrs"], dtype=float)

        vec = state_vectors[mission]
        v_sc_hat_ecl = get_spacecraft_velocity_unit_vector(mission, state_vectors)
        v_sc_mag = vec["speed_ecl_km_s"]
        v_sc_ecl = v_sc_hat_ecl * v_sc_mag
        logger.info(f"  Using 3D state vector: speed = {v_sc_mag:.2f} km/s")
        logger.info(f"    ecliptic velocity: [{v_sc_ecl[0]:+.2f}, {v_sc_ecl[1]:+.2f}, {v_sc_ecl[2]:+.2f}] km/s")

        # Equatorial velocity for CMB dot products must match the same J2000
        # obliquity rotation as Earth's velocity. step_038 already stores
        # Horizons-derived (vx, vy, vz)_equatorial; use those to avoid applying
        # the inverse rotation by mistake (which scrambled SC-CMB geometry).
        if all(k in vec for k in ("vx_km_s", "vy_km_s", "vz_km_s")):
            v_sc_equ = np.array([vec["vx_km_s"], vec["vy_km_s"], vec["vz_km_s"]])
        else:
            v_sc_equ = ecliptic_to_equatorial(float(v_sc_ecl[0]), float(v_sc_ecl[1]), float(v_sc_ecl[2]))

        # Solar scalar modulation
        solar_mod = compute_solar_scalar_modulation(r_au, lon_deg, v_earth_ecl, v_sc_ecl)
        logger.info(f"  Sun UCD saturation radius: {solar_mod['r_sol_sun_au']:.4f} AU")
        logger.info(f"  Heliocentric modulation proxy: {solar_mod['heliocentric_modulation']:.4f}")
        logger.info(f"  Earth radial velocity: {solar_mod['radial_velocity_kms']:.2f} km/s")
        logger.info(f"  SC-orbital alignment: {solar_mod['sc_orbital_alignment']:+.3f}")

        # CMB dipole modulation
        cmb_mod = compute_cmb_dipole_modulation(v_earth_equ, v_sc_equ)
        logger.info(f"  Earth barycentric v along CMB dipole: {cmb_mod['earth_orbital_cmb_proj_kms']:.2f} km/s")
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
            "has_3d_vector": True,
        })

    # Correlation analysis
    logger.section("CORRELATION ANALYSIS")

    # Filter to flybys with both observed and predicted data
    usable = [r for r in results
              if r["ratio_obs_pred"] is not None and not np.isnan(r["ratio_obs_pred"])]

    logger.info(f"Usable flybys for correlation: {len(usable)}")

    correlations = []
    sensitivity = None
    perm_p_optimal = None

    if len(usable) >= 3:
        ratios = [r["ratio_obs_pred"] for r in usable]
        diffs = [r["difference_mm_s"] for r in usable]

        cmb_proj_orb = [r["earth_orbital_cmb_proj_kms"] for r in usable]
        sc_orb_aligns = [r["sc_orbital_alignment"] for r in usable]
        sc_cmb_cos = [r["sc_cmb_cos_theta"] for r in usable]

        correlations = ratio_geometry_correlations(usable, ratios, "Ratio")
        correlations.extend(ratio_geometry_correlations(usable, diffs, "Difference"))

        for c in correlations:
            sig_str = " ***" if c.get("p_value") is not None and c["p_value"] < 0.05 else ""
            logger.info(f"  {c['test']}: r={c['pearson_r']:+.3f}, p={c['p_value']:.3f}, n={c['n']}{sig_str}")

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

        r_ba = p_ba = r_sd = p_sd = mw_stat = mw_p = None
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
        _abs_max, best_w, best_r = max_abs_pearson_weighted_E_vs_y(
            sc_cmb_cos, cmb_proj_orb, np.array(ratios, dtype=float)
        )
        e_best = np.asarray(sc_cmb_cos, dtype=float) + float(best_w) * (np.asarray(cmb_proj_orb, dtype=float) / 30.0)
        y_arr = np.asarray(ratios, dtype=float)
        mask_b = np.isfinite(e_best) & np.isfinite(y_arr)
        if mask_b.sum() >= 3:
            _rb, best_p = stats.pearsonr(e_best[mask_b], y_arr[mask_b])
            best_p = float(best_p)
        else:
            best_p = 1.0
        logger.info(f"Optimal model: E = SC-CMB_cos(theta) + {best_w:.3f} * Earth-CMB_proj/30")
        logger.info(f"Correlation: r = {best_r:+.3f}, p = {best_p:.3f}")
        perm_p_optimal = permutation_pvalue_max_abs_r_weighted(
            sc_cmb_cos, cmb_proj_orb, y_arr, n_perm=4000, seed=42
        )
        logger.info(f"Permutation p (max |r| over w grid, exchangeable ratios, n_perm=4000): {perm_p_optimal:.4f}")

        logger.section("SENSITIVITY CHECKS")
        sensitivity = {
            "ratio_excluding_missions": [],
            "alternate_endpoints": [],
            "notes": (
                "Rosetta 2007 exclusion: near-zero TEP prediction inflates obs/pred ratio. "
                "signed_log_ratio: |ratio|>0 only. sigma_z: published σ only."
            ),
        }
        ex_rows = [r for r in usable if r["mission"] != "Rosetta_2007"]
        if len(ex_rows) >= 3:
            ry_ex = [float(r["ratio_obs_pred"]) for r in ex_rows]
            cor_ex = ratio_geometry_correlations(ex_rows, ry_ex, "Ratio")
            sensitivity["ratio_excluding_missions"].append({
                "excluded": ["Rosetta_2007"],
                "rationale": "Near-zero TEP prediction inflates obs/pred ratio",
                "n": len(ex_rows),
                "correlations": cor_ex,
            })
            enh = [c for c in cor_ex if c["test"] == "Ratio vs CMB Disformal Enhancement"]
            if enh:
                logger.info(
                    f"  Excl. Rosetta 2007 (n={len(ex_rows)}): "
                    f"CMB disformal enhancement vs ratio r={enh[0]['pearson_r']:+.3f}, p={enh[0]['p_value']:.3f}"
                )

        sl_rows = [
            r for r in usable
            if r.get("ratio_obs_pred") is not None
            and float(r["ratio_obs_pred"]) != 0.0
            and np.isfinite(float(r["ratio_obs_pred"]))
        ]
        sl_y = [
            float(np.sign(float(r["ratio_obs_pred"])) * math.log10(max(abs(float(r["ratio_obs_pred"])), 1e-300)))
            for r in sl_rows
        ]
        if len(sl_rows) >= 3:
            cor_sl = ratio_geometry_correlations(sl_rows, sl_y, "signed_log_ratio")
            sensitivity["alternate_endpoints"].append({
                "name": "signed_log10_abs_ratio",
                "definition": "sign(ratio)*log10(max(|ratio|,1e-300)); |ratio|>0 only",
                "n": len(sl_rows),
                "missions": [r["mission"] for r in sl_rows],
                "correlations": cor_sl,
            })
            logger.info(f"  Signed log|ratio| endpoint: n={len(sl_rows)}")

        sz_rows = [
            r for r in usable
            if r.get("sigma_difference") is not None and np.isfinite(float(r["sigma_difference"]))
        ]
        sz_y = [float(r["sigma_difference"]) for r in sz_rows]
        if len(sz_rows) >= 3:
            cor_sz = ratio_geometry_correlations(sz_rows, sz_y, "sigma_z")
            sensitivity["alternate_endpoints"].append({
                "name": "sigma_normalized_residual",
                "definition": "(observed - predicted) / published σ where σ > 0",
                "n": len(sz_rows),
                "missions": [r["mission"] for r in sz_rows],
                "correlations": cor_sz,
            })
            logger.info(f"  σ-normalized residual endpoint: n={len(sz_rows)}")

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
        std_orig = None
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
            "uncertainty": None,
            "uncertainty_fraction": None,
            "uncertainty_absolute": None,
            "status": "cosmographic_test",
            "calibration_status": "horizons_3d_vectors_and_catalog_dates",
            "earth_sun_ephemeris": ephemeris_backend,
            "data_source": "step038 JPL Horizons raw responses and step003 flyby catalog",
            "recommended_action": "Extend Horizons coverage before drawing cosmographic conclusions for skipped catalog flybys",
            "derivation": "Counts of catalog flybys with geocentric Horizons 3D vectors at perigee",
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
            "original_std": float(std_orig) if std_orig is not None else None,
            "std_reduction_percent": float(reduction) if reduction is not None else None,
        } if beta_hat is not None else None,
        "directional_consistency": {
            "both_aligned_correlation_r": float(r_ba) if r_ba is not None else None,
            "both_aligned_correlation_p": float(p_ba) if p_ba is not None else None,
            "same_direction_correlation_r": float(r_sd) if r_sd is not None else None,
            "same_direction_correlation_p": float(p_sd) if p_sd is not None else None,
            "mann_whitney_u_statistic": float(mw_stat) if mw_stat is not None else None,
            "mann_whitney_u_pvalue": float(mw_p) if mw_p is not None else None,
        } if r_ba is not None else None,
        "optimal_weighted_combination": {
            "model": "E = SC-CMB_cos(theta) + w * Earth-CMB_proj/30",
            "w_grid_min": -2.0,
            "w_grid_max": 2.0,
            "w_grid_points": 201,
            "optimal_w": float(best_w) if best_w is not None else None,
            "pearson_r": float(best_r) if best_r is not None else None,
            "p_value": float(best_p) if best_p is not None else None,
            "n_permutations_max_abs_r": 4000,
            "permutation_random_seed": 42,
            "permutation_p_value_max_abs_r": float(perm_p_optimal) if perm_p_optimal is not None else None,
        } if best_w is not None else None,
        "sensitivity": sensitivity,
    }

    out_file = PROJECT_ROOT / "results" / "step040_cosmographic_shear.json"
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2)

    logger.info(f"\nResults saved to: {out_file}")
    logger.log_step_summary(0, "SUCCESS")
    return output


if __name__ == "__main__":
    main()
