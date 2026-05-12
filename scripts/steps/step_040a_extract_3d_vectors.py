"""
Step 040a: Extract 3D State Vectors from JPL Horizons Raw Responses

Parses existing JPL Horizons raw response files (which contain RA/Dec in HMS format,
range, and range-rate) to reconstruct full 3D position and velocity vectors
at perigee in geocentric equatorial and ecliptic coordinates.

No new web requests needed — works from cached raw data.
"""

import json
import sys
import re
import math
from pathlib import Path
from datetime import datetime

import numpy as np

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

# Obliquity of the ecliptic J2000
OBLIQUITY_RAD = math.radians(23.439281)


def hms_to_degrees(h, m, s):
    """Convert hours/minutes/seconds to decimal degrees."""
    sign = -1 if h < 0 else 1
    return sign * (abs(h) + m / 60.0 + s / 3600.0) * 15.0


def dms_to_degrees(d, m, s):
    """Convert degrees/minutes/seconds to decimal degrees."""
    sign = -1 if d < 0 else 1
    return sign * (abs(d) + m / 60.0 + s / 60.0)


def parse_hms(hms_str):
    """Parse 'HH MM SS.SS' to decimal hours."""
    parts = hms_str.strip().split()
    if len(parts) == 3:
        h, m, s = float(parts[0]), float(parts[1]), float(parts[2])
        return hms_to_degrees(h, m, s)
    return None


def parse_dms(dms_str):
    """Parse '+DD MM SS.S' to decimal degrees."""
    parts = dms_str.strip().split()
    if len(parts) == 3:
        d, m, s = float(parts[0]), float(parts[1]), float(parts[2])
        return dms_to_degrees(d, m, s)
    return None


def parse_raw_response(raw_path: Path) -> dict:
    """
    Parse JPL Horizons raw response CSV with RA/Dec in HMS format.

    Returns dict with arrays: timestamp, ra_deg, dec_deg, range_km, range_rate_kms
    """
    with open(raw_path, 'r', encoding='utf-8') as f:
        text = f.read()

    # Find data between $$SOE and $$EOE
    soe = text.find('$$SOE')
    eoe = text.find('$$EOE')
    if soe == -1 or eoe == -1:
        return None

    data_text = text[soe + 5:eoe].strip()
    lines = data_text.split('\n')

    timestamps = []
    ra_deg_list = []
    dec_deg_list = []
    range_km_list = []
    range_rate_kms_list = []

    for line in lines:
        line = line.strip()
        if not line or line.startswith('$'):
            continue
        # Split by comma
        parts = [p.strip() for p in line.split(',')]
        if len(parts) < 7:
            continue

        try:
            ts_str = parts[0]
            # Parse date like "1998-Jan-23 05:00"
            dt = datetime.strptime(ts_str, "%Y-%b-%d %H:%M")

            ra_str = parts[3].strip()
            dec_str = parts[4].strip()
            lt_str = parts[5].strip()
            delta_str = parts[6].strip()
            deldot_str = parts[7].strip() if len(parts) > 7 else '0'

            ra_deg = parse_hms(ra_str)
            dec_deg = parse_dms(dec_str)
            range_km = float(delta_str)
            range_rate_kms = float(deldot_str)

            if ra_deg is not None and dec_deg is not None:
                timestamps.append(dt)
                ra_deg_list.append(ra_deg)
                dec_deg_list.append(dec_deg)
                range_km_list.append(range_km)
                range_rate_kms_list.append(range_rate_kms)
        except Exception:
            continue

    return {
        "timestamps": timestamps,
        "ra_deg": np.array(ra_deg_list),
        "dec_deg": np.array(dec_deg_list),
        "range_km": np.array(range_km_list),
        "range_rate_kms": np.array(range_rate_kms_list),
    }


def compute_3d_state_vectors(data: dict) -> dict:
    """
    Compute full 3D position and velocity vectors in geocentric equatorial frame.

    Position: r_vec = r * [cos(dec)cos(ra), cos(dec)sin(ra), sin(dec)]
    Velocity: decomposed into radial + transverse from dRA/dt, dDec/dt, and range-rate.
    """
    ts = data["timestamps"]
    ra = np.radians(data["ra_deg"])
    dec = np.radians(data["dec_deg"])
    r = data["range_km"] * 1000.0  # meters
    rdot = data["range_rate_kms"] * 1000.0  # m/s

    n = len(ts)
    if n < 3:
        return None

    # Time differences in seconds
    dt = np.zeros(n)
    for i in range(1, n):
        dt[i] = (ts[i] - ts[i - 1]).total_seconds()
    dt[0] = dt[1] if n > 1 else 60.0

    # Angular rates (central differences)
    dra = np.zeros(n)
    ddec = np.zeros(n)
    for i in range(1, n - 1):
        dra[i] = (ra[i + 1] - ra[i - 1]) / (dt[i] + dt[i + 1])
        ddec[i] = (dec[i + 1] - dec[i - 1]) / (dt[i] + dt[i + 1])
    dra[0] = (ra[1] - ra[0]) / dt[1]
    dra[-1] = (ra[-1] - ra[-2]) / dt[-1]
    ddec[0] = (dec[1] - dec[0]) / dt[1]
    ddec[-1] = (dec[-1] - dec[-2]) / dt[-1]

    # Unit vectors in equatorial frame
    cos_dec = np.cos(dec)
    sin_dec = np.sin(dec)
    cos_ra = np.cos(ra)
    sin_ra = np.sin(ra)

    # Position vector (m)
    rx = r * cos_dec * cos_ra
    ry = r * cos_dec * sin_ra
    rz = r * sin_dec

    # Velocity components
    # Radial: along r_vec
    # Transverse: perpendicular to r_vec in RA and Dec directions
    # v = rdot * r_hat + r * cos(dec) * dRA/dt * phi_hat + r * dDec/dt * theta_hat
    # In Cartesian:
    # v_x = rdot*cos(dec)*cos(ra) - r*cos(dec)*dra*sin(ra) - r*ddec*cos(ra)*sin(dec)
    # v_y = rdot*cos(dec)*sin(ra) + r*cos(dec)*dra*cos(ra) - r*ddec*sin(ra)*sin(dec)
    # v_z = rdot*sin(dec) + r*ddec*cos(dec)

    vx = rdot * cos_dec * cos_ra - r * cos_dec * dra * sin_ra - r * ddec * cos_ra * sin_dec
    vy = rdot * cos_dec * sin_ra + r * cos_dec * dra * cos_ra - r * ddec * sin_ra * sin_dec
    vz = rdot * sin_dec + r * ddec * cos_dec

    # Speed check
    speed = np.sqrt(vx**2 + vy**2 + vz**2)

    return {
        "timestamps": ts,
        "ra_rad": ra,
        "dec_rad": dec,
        "r_m": r,
        "rx_m": rx,
        "ry_m": ry,
        "rz_m": rz,
        "vx_m_s": vx,
        "vy_m_s": vy,
        "vz_m_s": vz,
        "speed_m_s": speed,
        "rdot_m_s": rdot,
        "dra_rad_s": dra,
        "ddec_rad_s": ddec,
    }


def extract_perigee_state(state: dict) -> dict:
    """Find index of minimum range and return state at perigee."""
    idx = np.argmin(state["r_m"])
    ts = state["timestamps"][idx]

    return {
        "datetime_utc": ts.strftime("%Y-%m-%d %H:%M:%S"),
        "perigee_index": int(idx),
        "ra_deg": float(np.degrees(state["ra_rad"][idx])),
        "dec_deg": float(np.degrees(state["dec_rad"][idx])),
        "range_km": float(state["r_m"][idx] / 1000.0),
        "altitude_km": float(state["r_m"][idx] / 1000.0 - 6378.137),
        "rx_km": float(state["rx_m"][idx] / 1000.0),
        "ry_km": float(state["ry_m"][idx] / 1000.0),
        "rz_km": float(state["rz_m"][idx] / 1000.0),
        "vx_km_s": float(state["vx_m_s"][idx] / 1000.0),
        "vy_km_s": float(state["vy_m_s"][idx] / 1000.0),
        "vz_km_s": float(state["vz_m_s"][idx] / 1000.0),
        "speed_km_s": float(state["speed_m_s"][idx] / 1000.0),
        "rdot_km_s": float(state["rdot_m_s"][idx] / 1000.0),
    }


def equatorial_to_ecliptic(vx, vy, vz):
    """Rotate vector from equatorial to ecliptic frame (J2000)."""
    eps = OBLIQUITY_RAD
    x = vx
    y = vy * math.cos(eps) + vz * math.sin(eps)
    z = -vy * math.sin(eps) + vz * math.cos(eps)
    return x, y, z


def unit_vector(vx, vy, vz):
    """Return normalized vector."""
    mag = math.sqrt(vx**2 + vy**2 + vz**2)
    if mag == 0:
        return 0, 0, 0
    return vx / mag, vy / mag, vz / mag


def main():
    logger = StepLogger("step_040a_extract_3d_vectors", PROJECT_ROOT)
    logger.header("STEP 040a: EXTRACT 3D STATE VECTORS FROM JPL HORIZONS RAW DATA")

    raw_dir = PROJECT_ROOT / "data" / "raw" / "jpl_horizons"
    missions = {
        "NEAR_1998": "NEAR",
        "Galileo_1990": "Galileo_1990",
        "Galileo_1992": "Galileo_1992",
        "Cassini_1999": "Cassini",
        "Rosetta_2005": "Rosetta_2005",
        "Rosetta_2007": "Rosetta_2007",
        "Rosetta_2009": "Rosetta_2009",
        "MESSENGER_2005": "MESSENGER_2005",
        "Juno_2013": "Juno",
        "Stardust_2001": "Stardust_2001",
    }

    results = {}
    success_count = 0
    fail_count = 0

    for folder, mission_name in missions.items():
        raw_path = raw_dir / folder / f"{folder}_raw_response.txt"
        if not raw_path.exists():
            logger.warning(f"No raw response for {mission_name}: {raw_path}")
            fail_count += 1
            continue

        logger.subsection(f"Processing {mission_name}")

        data = parse_raw_response(raw_path)
        if data is None or len(data["timestamps"]) == 0:
            logger.warning(f"  Could not parse raw response for {mission_name}")
            fail_count += 1
            continue

        logger.info(f"  Parsed {len(data['timestamps'])} data points")

        state = compute_3d_state_vectors(data)
        if state is None:
            logger.warning(f"  Insufficient data for 3D state computation")
            fail_count += 1
            continue

        perigee = extract_perigee_state(state)

        # Convert velocity to ecliptic for solar/CMB analysis
        vx_eq = perigee["vx_km_s"]
        vy_eq = perigee["vy_km_s"]
        vz_eq = perigee["vz_km_s"]
        vx_ec, vy_ec, vz_ec = equatorial_to_ecliptic(vx_eq, vy_eq, vz_eq)

        perigee["vx_ecl_km_s"] = vx_ec
        perigee["vy_ecl_km_s"] = vy_ec
        perigee["vz_ecl_km_s"] = vz_ec
        perigee["speed_ecl_km_s"] = math.sqrt(vx_ec**2 + vy_ec**2 + vz_ec**2)

        # Unit vectors
        ux_ec, uy_ec, uz_ec = unit_vector(vx_ec, vy_ec, vz_ec)
        perigee["ux_ecl"] = ux_ec
        perigee["uy_ecl"] = uy_ec
        perigee["uz_ecl"] = uz_ec

        logger.info(f"  Perigee: RA={perigee['ra_deg']:.2f} deg, Dec={perigee['dec_deg']:.2f} deg")
        logger.info(f"  Range: {perigee['range_km']:.1f} km, Altitude: {perigee['altitude_km']:.1f} km")
        logger.info(f"  Speed (equatorial): {perigee['speed_km_s']:.2f} km/s")
        logger.info(f"  Velocity (ecliptic): [{vx_ec:+.2f}, {vy_ec:+.2f}, {vz_ec:+.2f}] km/s")

        results[mission_name] = perigee
        success_count += 1

    logger.section("SUMMARY")
    logger.info(f"Successfully processed: {success_count}")
    logger.info(f"Failed: {fail_count}")

    out_file = PROJECT_ROOT / "results" / "step040a_3d_state_vectors.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Saved 3D state vectors to: {out_file}")
    logger.log_step_summary(0, "SUCCESS")
    return results


if __name__ == "__main__":
    main()
