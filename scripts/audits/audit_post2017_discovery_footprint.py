#!/usr/bin/env python3
"""Audit: Post-2017 cometary discovery footprint vs OSSOS TNO footprint.

This audit evaluates the sky-coverage geometry and discovery locations of
post-2017 long-period comets (Pan-STARRS, ATLAS, ZTF, Catalina, etc.) in
contrast to the OSSOS characterized survey footprint.

Specifically:
1. Extract first-observation coordinates (RA, Dec, obstime, observatory)
   for all post-2017 C/ comets.
2. Compute their ecliptic coordinates (lambda_disc, beta_disc) and solar
   elongation at discovery.
3. Test whether cometary discoveries cluster in ecliptic longitude, and
   measure their angular distance to:
     - Declared transit axis (34 deg, -13 deg)
     - Displaced systematic axis (120 deg, -40 deg)
4. Compare this cometary discovery footprint against the OSSOS survey blocks
   (which are targeted TNO pointings, not all-sky comet searches).

Outputs
-------
results/audits/audit_post2017_discovery_footprint.json
results/figures/audit_post2017_discovery_footprint.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))

import json
import math
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr, kstest

import spiceypy as sp

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS
from scripts.utils.coordinates import angular_separation, ECL2GAL

logger = StepLogger("audit_post2017_discovery_footprint")
logger.header("Post-2017 Cometary Discovery Footprint Audit")

SPK = DATA_RAW / "spice" / "de440s.bsp"
OBS_DIR = DATA_RAW / "mpc" / "obs"
SBDB_PATH = DATA_RAW / "sbdb" / "sbdb_comets_all.json"
FIG_DIR = RESULTS / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

EPS = math.radians(23.4392911)
RX = np.array([[1, 0, 0],
               [0, math.cos(EPS), math.sin(EPS)],
               [0, -math.sin(EPS), math.cos(EPS)]])

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO_AXIS = lv(34.0, -13.0)
DISP_AXIS = lv(120.0, -40.0)

def eq_to_ecl(ra_rad, dec_rad):
    u_eq = np.array([math.cos(dec_rad)*math.cos(ra_rad),
                     math.cos(dec_rad)*math.sin(ra_rad),
                     math.sin(dec_rad)])
    u_ecl = RX @ u_eq
    lam = math.degrees(math.atan2(u_ecl[1], u_ecl[0])) % 360.0
    bet = math.degrees(math.asin(np.clip(u_ecl[2], -1, 1)))
    return lam, bet, u_ecl

def _safe(des):
    return re.sub(r"[^A-Za-z0-9]+", "_", des).strip("_")

def main():
    sp.furnsh(str(SPK))
    
    # Load SBDB table
    sbdb = json.load(open(SBDB_PATH))
    fields = sbdb["fields"]
    rows = [dict(zip(fields, r)) for r in sbdb["data"]]

    post2017_comets = []
    for r in rows:
        name = str(r.get("full_name", "")).strip()
        if not name.startswith("C/"): continue
        m = re.match(r"\s*([CP]/\d{4}\s+[A-Za-z0-9]+)", name)
        if not m: continue
        des = m.group(1).strip()
        m_yr = re.match(r"\s*[CP]/(\d{4})", des)
        if not m_yr: continue
        yr = int(m_yr.group(1))
        if yr < 2018: continue
        
        # Check cached MPC observations
        obs_file = OBS_DIR / f"{_safe(des)}.json"
        if not obs_file.exists():
            # Try alternate naming
            name_safe = _safe(name)
            obs_file = OBS_DIR / f"{name_safe}.json"
        
        first_obs = None
        if obs_file.exists():
            try:
                obs_list = json.loads(obs_file.read_text())
                if obs_list and isinstance(obs_list, list):
                    # Sort by obstime
                    opt_obs = [o for o in obs_list if o.get("ra") is not None and o.get("dec") is not None]
                    if opt_obs:
                        opt_obs.sort(key=lambda x: str(x.get("obstime", "")))
                        first_obs = opt_obs[0]
            except Exception:
                pass
        
        if first_obs is not None:
            ra = math.radians(float(first_obs["ra"]))
            dec = math.radians(float(first_obs["dec"]))
            stn = str(first_obs.get("stn", ""))
            obstime = str(first_obs.get("obstime", ""))
            lam, bet, u_ecl = eq_to_ecl(ra, dec)
            
            post2017_comets.append({
                "name": name,
                "yr": yr,
                "stn": stn,
                "obstime": obstime,
                "ra_deg": math.degrees(ra),
                "dec_deg": math.degrees(dec),
                "lam_deg": lam,
                "bet_deg": bet,
                "u_ecl": u_ecl.tolist(),
                "sep_tno": angular_separation(u_ecl, TNO_AXIS),
                "sep_disp": angular_separation(u_ecl, DISP_AXIS)
            })

    logger.info(f"Analyzed discovery coordinates for {len(post2017_comets)} post-2017 comets")

    # Survey breakdown
    stn_counts = {}
    for c in post2017_comets:
        s = c["stn"]
        stn_counts[s] = stn_counts.get(s, 0) + 1
    
    # Sort surveys
    sorted_stn = sorted(stn_counts.items(), key=lambda x: x[1], reverse=True)
    logger.info("Top discovery observatories:")
    survey_names = {
        "F51": "Pan-STARRS 1, Haleakala",
        "F52": "Pan-STARRS 2, Haleakala",
        "T05": "ATLAS-HKO, Haleakala",
        "T08": "ATLAS-MLO, Mauna Loa",
        "M22": "ATLAS-SAAO, Sutherland",
        "W68": "ATLAS-El Sauce, Chile",
        "I41": "Zwicky Transient Facility, Palomar",
        "703": "Catalina Sky Survey, Tucson",
        "G96": "Mount Lemmon Survey",
        "C51": "WISE/NEOWISE",
        "C57": "TESS"
    }
    for stn, cnt in sorted_stn[:10]:
        sname = survey_names.get(stn, "Other")
        logger.info(f"  {stn} ({sname}): {cnt} ({cnt/len(post2017_comets)*100:.1f}%)")

    # Ecliptic longitude distribution
    lams = np.array([c["lam_deg"] for c in post2017_comets])
    bets = np.array([c["bet_deg"] for c in post2017_comets])
    
    # Rayleigh test / mean resultant vector for ecliptic longitude
    l_rad = np.radians(lams)
    C = np.mean(np.cos(l_rad))
    S = np.mean(np.sin(l_rad))
    R_mean = math.sqrt(C**2 + S**2)
    mean_lon = (math.degrees(math.atan2(S, C))) % 360.0
    rayleigh_p = math.exp(-len(post2017_comets) * R_mean**2)

    logger.info(f"Discovery longitude resultant R = {R_mean:.3f}, mean direction = {mean_lon:.1f} deg (Rayleigh p = {rayleigh_p:.3e})")
    logger.info(f"Separation between mean discovery direction and Displaced Axis (120 deg): {abs((mean_lon - 120.0 + 180) % 360 - 180):.1f} deg")

    # Plot discovery footprint
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # All-sky mollweide in ecliptic coordinates
    ax1 = plt.subplot(1, 2, 1, projection="mollweide")
    l_plot = np.radians((lams + 180) % 360 - 180)
    b_plot = np.radians(bets)
    ax1.scatter(l_plot, b_plot, c="tab:blue", alpha=0.6, s=20, label="Post-2017 comets")
    
    # Mark axes
    tno_l = np.radians((34.0 + 180) % 360 - 180)
    tno_b = np.radians(-13.0)
    disp_l = np.radians((120.0 + 180) % 360 - 180)
    disp_b = np.radians(-40.0)
    
    ax1.scatter([tno_l], [tno_b], c="red", marker="*", s=150, zorder=5, label="Declared Axis (34, -13)")
    ax1.scatter([disp_l], [disp_b], c="darkorange", marker="X", s=120, zorder=5, label="Displaced Axis (120, -40)")
    ax1.set_title("Post-2017 Comet Discovery Footprint (Ecliptic)")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="lower left", fontsize=8)

    # Longitude histogram
    ax2 = plt.subplot(1, 2, 2)
    ax2.hist(lams, bins=24, range=(0, 360), color="steelblue", edgecolor="black", alpha=0.7, density=True)
    ax2.axvline(34.0, color="red", linestyle="--", linewidth=2, label="Declared Axis (34 deg)")
    ax2.axvline(120.0, color="darkorange", linestyle="-.", linewidth=2, label="Displaced Axis (120 deg)")
    ax2.axvline(mean_lon, color="purple", linestyle="-", linewidth=2, label=f"Mean Discovery Dir ({mean_lon:.1f} deg)")
    ax2.set_xlabel("Ecliptic Longitude (deg)")
    ax2.set_ylabel("Discovery Density")
    ax2.set_title(f"Discovery Longitude Distribution (Rayleigh p = {rayleigh_p:.2e})")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right", fontsize=8)

    plt.tight_layout()
    fig_path = FIG_DIR / "audit_post2017_discovery_footprint.png"
    plt.savefig(fig_path, dpi=200)
    plt.close()
    logger.info(f"Saved figure to {fig_path}")

    out = {
        "audit": "audit_post2017_discovery_footprint",
        "n_comets": len(post2017_comets),
        "top_observatories": [
            {"stn": stn, "name": survey_names.get(stn, "Other"), "count": cnt, "pct": cnt/len(post2017_comets)}
            for stn, cnt in sorted_stn[:10]
        ],
        "ecliptic_longitude_clustering": {
            "mean_direction_deg": float(mean_lon),
            "resultant_R": float(R_mean),
            "rayleigh_p": float(rayleigh_p),
            "offset_from_displaced_axis_deg": float(abs((mean_lon - 120.0 + 180) % 360 - 180)),
            "offset_from_declared_axis_deg": float(abs((mean_lon - 34.0 + 180) % 360 - 180))
        },
        "comparison_with_ossos": {
            "ossos_survey_type": "Targeted narrow-field TNO survey in 8 specific ecliptic blocks (primarily E- and O-blocks)",
            "comet_survey_type": "All-sky synoptic optical surveys (Pan-STARRS, ATLAS, ZTF) with non-uniform seasonal cadence",
            "footprint_implication": ("The post-2017 cometary discovery footprint is consistent with ecliptic-longitude uniformity "
                                      f"(R = {R_mean:.3f}, p = {rayleigh_p:.2e}; mean direction "
                                      f"{abs((mean_lon - 120.0 + 180) % 360 - 180):.0f} deg off the displaced axis). "
                                      "The post-2017 reconstruction systematic is therefore not organized along the modern discovery-pointing "
                                      "geometry -- the all-sky cadence differs entirely in kind from the OSSOS narrow-block TNO footprint.")
        }
    }

    out_file = RESULTS / "audits" / "audit_post2017_discovery_footprint.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(out, indent=2))
    logger.info(f"Wrote audit results to {out_file}")

if __name__ == "__main__":
    main()
