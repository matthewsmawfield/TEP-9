#!/usr/bin/env python3
"""Audit: Ephemeris-version sensitivity (DE440s vs DE430).

This audit tests whether the cometary boundary reconstruction and the
in-cap vs out-of-cap rotation gap are sensitive to the planetary ephemeris
version.

CODE and Warsaw catalogue solutions were computed in the DE405/DE430 era,
whereas the primary TEP-9 baseline re-integrates orbits using DE440s.
This test integrates the class-1 CODE comets under both DE440s and DE430
(REBOUND/IAS15) to determine:
  1. The direct per-leg state and direction difference (|p_DE440 - p_DE430|).
  2. The stability of the simulated rotation drot_sim.
  3. The stability of the in-cap vs out-of-cap rotation gap and p-value.

Outputs
-------
results/audits/audit_ephemeris_de430_vs_de440.json
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))

import csv
import hashlib
import json
import math
import re
import urllib.request
import numpy as np
from html.parser import HTMLParser
from scipy.stats import mannwhitneyu, spearmanr

import rebound
import spiceypy as sp

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS
from scripts.utils.coordinates import angular_separation

logger = StepLogger("audit_ephemeris_de430_vs_de440")
logger.header("Ephemeris-Version Sensitivity Audit: DE440s vs DE430")

DE430_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp"
SPK_440 = DATA_RAW / "spice" / "de440s.bsp"
SPK_430 = DATA_RAW / "spice" / "de430.bsp"
PROV_FILE = DATA_RAW / "spice" / "provenance.json"

AU_KM   = 149597870.7
DAY_YR  = 365.25
EPS     = math.radians(23.4392911)
RX      = np.array([[1, 0, 0],
                    [0, math.cos(EPS), math.sin(EPS)],
                    [0, -math.sin(EPS), math.cos(EPS)]])

GM_SUN = 1.32712440018e11
GM = {"1": 2.2031868551e4, "2": 3.2485859200e5, "3": 4.0350323562e5,
      "4": 4.2828375814e4, "5": 1.2671276480e8, "6": 3.7940626000e7,
      "7": 5.7945490100e6, "8": 6.8365271006e6, "9": 1.0868657e3}
PLANET_IDS = list(GM.keys())

MU = 4 * math.pi ** 2
R_STOP = 250.0
T_MAX  = 20000.0
DT_OUT = 1.0

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO = lv(34.0, -13.0)

def ensure_de430():
    if not SPK_430.exists():
        logger.info(f"Downloading DE430 from {DE430_URL} (~115 MB)...")
        urllib.request.urlretrieve(DE430_URL, SPK_430)
        logger.info("Download complete.")
    
    sha256 = hashlib.sha256(SPK_430.read_bytes()).hexdigest()
    logger.info(f"DE430 SHA-256: {sha256}")
    
    # Update provenance
    if PROV_FILE.exists():
        prov = json.loads(PROV_FILE.read_text())
    else:
        prov = {"files": {}}
    if "files" not in prov:
        prov["files"] = {}
    prov["files"]["de430.bsp"] = {
        "url": DE430_URL,
        "bytes": SPK_430.stat().st_size,
        "sha256": sha256,
        "note": "JPL DE430 ephemeris, coverage 1550-2650"
    }
    PROV_FILE.write_text(json.dumps(prov, indent=2))
    return sha256

class TP(HTMLParser):
    def __init__(self):
        super().__init__(); self.rows = []; self.cur = []; self.buf = ""; self.in_td = False
    def handle_starttag(self, t, a):
        if t == "tr": self.cur = []
        elif t == "td": self.in_td = True; self.buf = ""
    def handle_endtag(self, t):
        if t == "td": self.in_td = False; self.cur.append(self.buf.strip())
        elif t == "tr" and self.cur: self.rows.append(self.cur)
    def handle_data(self, d):
        if self.in_td: self.buf += d

def parse_code(path):
    p = TP(); p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14: continue
        try:
            out[r[0].strip()] = dict(desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                T=r[7].strip(),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out

def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return angular_separation(a, b)

def jd_tt(y, m, d):
    if m <= 2: y -= 1; m += 12
    A = y // 100; B = 2 - A + A // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5

def perihelion_et(T):
    m = re.match(r"(\d{4})\s+(\d{1,2})\s+(\d+\.?\d*)", T)
    if not m:
        return None
    return (jd_tt(int(m.group(1)), int(m.group(2)), float(m.group(3))) - 2451545.0) * 86400.0

def body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return RX @ np.array(st[:3]) / AU_KM, RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR

def init_sim(et):
    sim = rebound.Simulation()
    sim.G = MU
    ps, vs = body_state("10", et)
    sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2], m=1.0)
    for b in PLANET_IDS:
        pp, vv = body_state(b, et)
        sim.add(x=pp[0], y=pp[1], z=pp[2], vx=vv[0], vy=vv[1], vz=vv[2],
                m=GM[b] / GM_SUN)
    return sim, ps, vs

def state_at_periapsis(ro):
    w_, O_, i_ = map(math.radians, (ro["w"], ro["Om"], ro["i"]))
    ph = perih_dir(w_, O_, i_)
    hh = np.array([math.sin(i_) * math.sin(O_), -math.sin(i_) * math.cos(O_), math.cos(i_)])
    th = np.cross(hh, ph); th /= np.linalg.norm(th)
    rvec = ro["q"] * ph
    vvec = math.sqrt(MU * (1 + ro["e"]) / ro["q"]) * th
    return rvec, vvec

def boundary_orbit(r_rel, v_rel, mtot):
    mu = MU * mtot
    r = np.linalg.norm(r_rel)
    h = np.cross(r_rel, v_rel)
    evec = (np.cross(v_rel, h) / mu) - r_rel / r
    en = np.linalg.norm(evec)
    phat = evec / en if en > 1e-12 else r_rel / r
    E = v_rel.dot(v_rel) / 2 - mu / r
    return phat, -2 * E / mu * 1e6

def integrate_leg(ro, et, direction):
    sim, ps, vs = init_sim(et)
    rvec, vvec = state_at_periapsis(ro)
    sim.add(x=rvec[0] + ps[0], y=rvec[1] + ps[1], z=rvec[2] + ps[2],
            vx=vvec[0] + vs[0], vy=vvec[1] + vs[1], vz=vvec[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    denc = np.full(sim.N - 1, np.inf)
    t = direction * DT_OUT
    while abs(t) < T_MAX:
        sim.integrate(t, exact_finish_time=0)
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y, p[nc].z - p[0].z])
        for j in range(1, sim.N - 1):
            d = math.sqrt((p[nc].x - p[j].x) ** 2 + (p[nc].y - p[j].y) ** 2 +
                          (p[nc].z - p[j].z) ** 2)
            if d < denc[j - 1]: denc[j - 1] = d
        if np.linalg.norm(r_rel) >= R_STOP:
            break
        t += direction * DT_OUT
    else:
        return None
    mtot = sum(pp.m for pp in sim.particles)
    rb = np.zeros(3); vb = np.zeros(3)
    for pp in sim.particles:
        rb += pp.m * np.array([pp.x, pp.y, pp.z])
        vb += pp.m * np.array([pp.vx, pp.vy, pp.vz])
    rb /= mtot; vb /= mtot
    p = sim.particles
    r_rel = np.array([p[nc].x, p[nc].y, p[nc].z]) - rb
    v_rel = np.array([p[nc].vx, p[nc].vy, p[nc].vz]) - vb
    phat, aa = boundary_orbit(r_rel, v_rel, mtot)
    return dict(phat=phat, aa=aa, denc=float(denc.min()), t_years=float(t))

def resid_cap(v, th, X_extra):
    inc = th < 60
    y = np.log10(v)
    X = np.column_stack([np.ones(len(y))] + X_extra)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    u = mannwhitneyu(resid[inc], resid[~inc], alternative="greater")
    rho, p = spearmanr(th, resid)
    return {"p": float(u.pvalue), "rho": float(rho), "p_2sided": float(p),
            "med_resid_in": float(np.median(resid[inc])),
            "med_resid_out": float(np.median(resid[~inc])),
            "gap": float(np.median(resid[inc]) - np.median(resid[~inc]))}

def main():
    sha430 = ensure_de430()
    
    osc  = parse_code(str(DATA_RAW / "code" / "code_osculating.html"))
    orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))
    fut  = parse_code(str(DATA_RAW / "code" / "code_future.html"))
    warsaw = {l[5:17].strip() for l in open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")) if len(l) > 115}

    sample = []
    for k, ro in orig.items():
        if k not in fut: continue
        if not (0 < ro["aa"] < 100 and ro["cls"] in ("1a", "1a+", "1b")) or k in warsaw:
            continue
        if k not in osc: continue
        et = perihelion_et(osc[k]["T"])
        if et is None: continue
        sample.append((k, ro, osc[k], et))

    logger.info(f"sample: {len(sample)} class-1 CODE comets")

    # Run under DE440s
    logger.info("Integrating under DE440s...")
    sp.kclear()
    sp.furnsh(str(SPK_440))
    res_440 = {}
    for k, ro, oo, et in sample:
        rb = integrate_leg(oo, et, -1)
        rf = integrate_leg(oo, et, +1)
        if rb and rf:
            res_440[k] = {"rb": rb, "rf": rf,
                          "drot": sep(rb["phat"], rf["phat"]),
                          "denc": min(rb["denc"], rf["denc"]),
                          "daa": rf["aa"] - rb["aa"]}

    # Run under DE430
    logger.info("Integrating under DE430...")
    sp.kclear()
    sp.furnsh(str(SPK_430))
    res_430 = {}
    for k, ro, oo, et in sample:
        rb = integrate_leg(oo, et, -1)
        rf = integrate_leg(oo, et, +1)
        if rb and rf:
            res_430[k] = {"rb": rb, "rf": rf,
                          "drot": sep(rb["phat"], rf["phat"]),
                          "denc": min(rb["denc"], rf["denc"]),
                          "daa": rf["aa"] - rb["aa"]}

    # Direct comparison on shared successful integrations
    common = sorted(set(res_440.keys()) & set(res_430.keys()))
    logger.info(f"Successfully integrated {len(common)} comets under both ephemerides")

    diff_in_deg = []
    diff_out_deg = []
    diff_drot_deg = []
    
    rows_440 = []
    rows_430 = []

    for k in common:
        r440 = res_440[k]
        r430 = res_430[k]
        ro = orig[k]
        oo = osc[k]
        po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
        th = sep(-po, TNO)

        din = sep(r440["rb"]["phat"], r430["rb"]["phat"])
        dout = sep(r440["rf"]["phat"], r430["rf"]["phat"])
        ddrot = abs(r440["drot"] - r430["drot"])

        diff_in_deg.append(din)
        diff_out_deg.append(dout)
        diff_drot_deg.append(ddrot)

        rows_440.append(dict(desig=k, theta=th, drot=r440["drot"],
                             daa=r440["daa"], denc=r440["denc"], q=ro["q"], i=ro["i"]))
        rows_430.append(dict(desig=k, theta=th, drot=r430["drot"],
                             daa=r430["daa"], denc=r430["denc"], q=ro["q"], i=ro["i"]))

    def eval_stats(rows):
        th = np.array([r["theta"] for r in rows])
        v = np.array([r["drot"] for r in rows])
        K = np.abs(np.array([r["daa"] for r in rows]))
        Q = np.array([r["q"] for r in rows])
        I = np.array([r["i"] for r in rows])
        D = np.array([r["denc"] for r in rows])
        
        # 1-covariate baseline: log10(|daa|+1)
        b1 = resid_cap(v, th, [np.log10(K + 1.0)])
        # 4-covariate baseline: log10(|daa|+1), log10(denc), q, i
        b4 = resid_cap(v, th, [np.log10(K + 1.0), np.log10(np.clip(D, 1e-3, None)), Q, I])
        
        # Raw cap contrast
        inc = th < 60
        u_raw = mannwhitneyu(v[inc], v[~inc], alternative="greater")
        
        return {
            "n": len(rows),
            "n_in": int(inc.sum()),
            "raw_med_in": float(np.median(v[inc])),
            "raw_med_out": float(np.median(v[~inc])),
            "raw_gap": float(np.median(v[inc]) - np.median(v[~inc])),
            "raw_p": float(u_raw.pvalue),
            "b1_kick": b1,
            "b4_covariates": b4
        }

    stats_440 = eval_stats(rows_440)
    stats_430 = eval_stats(rows_430)

    out = {
        "audit": "audit_ephemeris_de430_vs_de440",
        "n_comets": len(common),
        "de430_sha256": sha430,
        "direct_comparison": {
            "inbound_asymptote_diff_deg": {
                "median": float(np.median(diff_in_deg)),
                "max": float(np.max(diff_in_deg)),
                "p95": float(np.percentile(diff_in_deg, 95))
            },
            "outbound_asymptote_diff_deg": {
                "median": float(np.median(diff_out_deg)),
                "max": float(np.max(diff_out_deg)),
                "p95": float(np.percentile(diff_out_deg, 95))
            },
            "rotation_diff_deg": {
                "median": float(np.median(diff_drot_deg)),
                "max": float(np.max(diff_drot_deg)),
                "p95": float(np.percentile(diff_drot_deg, 95))
            }
        },
        "stats_de440s": stats_440,
        "stats_de430": stats_430,
        "stability_summary": {
            "raw_gap_de440": stats_440["raw_gap"],
            "raw_gap_de430": stats_430["raw_gap"],
            "raw_gap_diff": stats_440["raw_gap"] - stats_430["raw_gap"],
            "b4_p_de440": stats_440["b4_covariates"]["p"],
            "b4_p_de430": stats_430["b4_covariates"]["p"],
            "b4_gap_de440": stats_440["b4_covariates"]["gap"],
            "b4_gap_de430": stats_430["b4_covariates"]["gap"],
            "is_stable": bool(abs(stats_440["raw_gap"] - stats_430["raw_gap"]) < 0.005 and
                              abs(stats_440["b4_covariates"]["gap"] - stats_430["b4_covariates"]["gap"]) < 0.005)
        }
    }

    out_file = RESULTS / "audits" / "audit_ephemeris_de430_vs_de440.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(out, indent=2))
    logger.info(f"Wrote audit results to {out_file}")
    
    logger.info(f"Stability Verdict: {'STABLE' if out['stability_summary']['is_stable'] else 'SENSITIVE'}")
    logger.info(f"  DE440s b4 gap: {stats_440['b4_covariates']['gap']:.4f} dex (p={stats_440['b4_covariates']['p']:.4g})")
    logger.info(f"  DE430  b4 gap: {stats_430['b4_covariates']['gap']:.4f} dex (p={stats_430['b4_covariates']['p']:.4g})")
    logger.info(f"  Median per-comet rotation shift: {np.median(diff_drot_deg):.6f} deg")

if __name__ == "__main__":
    main()
