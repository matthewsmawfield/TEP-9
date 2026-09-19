"""step_062: Planet Nine insertion test -- does the published perturber
produce the observed comet anomaly on real trajectories?

step_042 established that subtracting the measured planetary baseline
leaves the in-cap anomaly undiminished.  This step asks the inverse
question for the specific conventional hypothesis: insert the Brown &
Batygin (2021, AJ 162, 219) Planet Nine into the same backward
integrations and measure what it does to each comet's reconstructed
orbit.

Two published realizations are tested:

  bb21_ml   maximum-likelihood model: m9 = 5.0 M_earth, a9 = 300 AU,
            e9 = 0.15, i9 = 17 deg, varpi9 = 254 deg, Om9 = 108 deg
            (BB21 Section on maximum likelihood, "nearly a local peak
            in every dimension")
  bb21_med  marginalized-median model: m9 = 6.9 M_earth, a9 = 461 AU,
            e9 = 0.30, i9 = 15.6 deg, varpi9 = 246.7 deg, Om9 = 96.9 deg
            (BB21 corner-plot medians)

P9's mean anomaly is unconstrained by the secular analysis (BB21 set
M9 = 0 and marginalized the position), so each model is run on a
four-point mean-anomaly grid {0, 90, 180, 270 deg}; M9 = 180 deg
(near apoapsis) is the residence-time-weighted most probable position
for an eccentric orbit.  This brackets every plausible placement.

For each class-1 CODE comet (identical selection to step_040/042) the
osculating elements at the catalogued perihelion epoch are converted
to a heliocentric state and integrated backward through Sun +
planetary-system barycentres (DE440s) to 255 AU heliocentric, twice:
with and without the P9 particle (REBOUND/IAS15).  At the boundary the
osculating periapsis direction and barycentric 1/a are recovered for
each run.  The paired differences

  drot_p9  = angle between boundary periapsis directions, P9 vs no-P9
  dk_p9    = |1/a_out(P9) - 1/a_out(noP9)|   [10^-6 AU^-1]
  denc_p9  = minimum comet-P9 distance [AU]

isolate the perturbation P9 would actually inject into the catalogued
original-orbit reconstruction.  The observed discrepancy d_of is the
orig->fut perihelion-direction residual used by every Phase-3 step.

Tests per (model, M9) cell and pooled:
  1. magnitude: median/max drot_p9 against the observed in-cap median
     d_of -- the trajectory-level impulse budget;
  2. targeting: Spearman(drot_p9, d_of) and in/out-cap Mann-Whitney on
     drot_p9 -- a causal perturber must hit the discrepant comets;
  3. absorption: log(d_of) ~ drot_p9 + |kick_e| + q residual cap
     contrast -- if P9 explained the signal the residual collapses.

Outputs
-------
results/step_b27_planet_nine_insertion.json
results/step_b27_planet_nine_insertion.csv   (per-comet P9 kicks)
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_062_planet_nine_insertion")
tee_stdout(logger)
logger.header("Planet Nine insertion test (REBOUND/IAS15 + DE440s)")

import csv
import json
import math
import re
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
from html.parser import HTMLParser

import rebound
import spiceypy as sp

# ------------------------------------------------------------------
# CODE table parser (identical to step_042)
# ------------------------------------------------------------------

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
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO = lv(34, -13)   # same axis convention as every Phase-3 transit step

# ------------------------------------------------------------------
# Ephemeris setup (identical to step_042)
# ------------------------------------------------------------------

SPK = DATA_RAW / "spice" / "de440s.bsp"
if not SPK.exists():
    raise FileNotFoundError(f"JPL ephemeris missing: {SPK}")
sp.furnsh(str(SPK))

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
BODY_NAME = {"1": "MeB", "2": "VB", "3": "EMB", "4": "MaB", "5": "JuB",
             "6": "SaB", "7": "UB", "8": "NB", "9": "PlB"}
PLANET_IDS = list(GM.keys())

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

# ------------------------------------------------------------------
# Planet Nine realizations (Brown & Batygin 2021, AJ 162, 219)
# ------------------------------------------------------------------

M_EARTH_SUN = 3.0035e-6          # M_earth / M_sun

P9_MODELS = {
    "bb21_ml":  dict(m9=5.0, a9=300.0, e9=0.15, i9=17.0,
                     varpi9=254.0, Om9=108.0,
                     note="BB21 maximum-likelihood point"),
    "bb21_med": dict(m9=6.9, a9=461.0, e9=0.30, i9=15.6,
                     varpi9=246.7, Om9=96.9,
                     note="BB21 marginalized-posterior medians"),
}
M9_GRID = [0.0, 90.0, 180.0, 270.0]   # deg; M9 unconstrained by BB21

def p9_state(model, m9_deg, et):
    """Heliocentric ecliptic state of P9 at epoch et (rebound elements
    conversion).  Mean anomaly propagated from J2000 by Kepler's law."""
    w9 = model["varpi9"] - model["Om9"]
    P9_yr = model["a9"] ** 1.5                     # Kepler, solar mass
    et_yr = et / (86400.0 * DAY_YR)
    M_at_et = math.radians(m9_deg) + (2 * math.pi / P9_yr) * et_yr
    sim = rebound.Simulation(); sim.G = MU
    sim.add(m=1.0)
    sim.add(a=model["a9"], e=model["e9"], inc=math.radians(model["i9"]),
            Omega=math.radians(model["Om9"]), omega=math.radians(w9),
            M=M_at_et, m=model["m9"] * M_EARTH_SUN)
    p = sim.particles[1]
    return (np.array([p.x, p.y, p.z]), np.array([p.vx, p.vy, p.vz]))

# ------------------------------------------------------------------
# REBOUND integration (extends step_042 with the P9 particle)
# ------------------------------------------------------------------

MU = 4 * math.pi ** 2
R_STOP = 255.0
T_MAX  = -20000.0
DT_OUT = 1.0

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
    """Osculating periapsis direction and barycentric 1/a [1e-6 AU^-1]
    from a heliocentric state relative to total mass mtot."""
    mu = MU * mtot
    r = np.linalg.norm(r_rel)
    h = np.cross(r_rel, v_rel)
    evec = (np.cross(v_rel, h) / mu) - r_rel / r
    en = np.linalg.norm(evec)
    phat = evec / en if en > 1e-12 else r_rel / r
    E = v_rel.dot(v_rel) / 2 - mu / r
    aa = -2 * E / mu * 1e6
    return phat, aa

def integrate(ro, et, p9=None):
    """Backward-integrate one comet to 255 AU.  p9 = (pos, vel, mass)
    heliocentric ecliptic for the inserted perturber, or None.

    Returns dict(phat, aa_out, denc_planet, denc_p9) or None."""
    sim, ps, vs = init_sim(et)
    np9 = -1
    if p9 is not None:
        sim.add(x=p9[0][0], y=p9[0][1], z=p9[0][2],
                vx=p9[1][0], vy=p9[1][1], vz=p9[1][2], m=p9[2])
        np9 = sim.N - 1
    rvec, vvec = state_at_periapsis(ro)
    sim.add(x=rvec[0] + ps[0], y=rvec[1] + ps[1], z=rvec[2] + ps[2],
            vx=vvec[0] + vs[0], vy=vvec[1] + vs[1], vz=vvec[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    denc = np.full(sim.N - 1, np.inf)
    t = -DT_OUT
    # t is negative for the backward leg; compare elapsed time, not signed t.
    while abs(t) < abs(T_MAX):
        sim.integrate(t, exact_finish_time=0)
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y, p[nc].z - p[0].z])
        for j in range(1, sim.N - 1):
            d = math.sqrt((p[nc].x - p[j].x) ** 2 + (p[nc].y - p[j].y) ** 2 +
                          (p[nc].z - p[j].z) ** 2)
            if d < denc[j - 1]: denc[j - 1] = d
        if np.linalg.norm(r_rel) >= R_STOP:
            break
        t -= DT_OUT
    else:
        return None
    mtot = sum(pp.m for pp in sim.particles)
    p = sim.particles
    r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y, p[nc].z - p[0].z])
    v_rel = np.array([p[nc].vx - p[0].vx, p[nc].vy - p[0].vy, p[nc].vz - p[0].vz])
    phat, aa_out = boundary_orbit(r_rel, v_rel, mtot)
    denc_p9 = float(denc[np9 - 1]) if np9 > 0 else float("nan")
    return dict(phat=phat, aa_out=aa_out,
                denc_planet=float(np.min(denc[:9])), denc_p9=denc_p9)

# ------------------------------------------------------------------
# Sample: identical selection to step_040 / step_042
# ------------------------------------------------------------------

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

# ------------------------------------------------------------------
# Integrate: one no-P9 reference + each (model, M9) cell
# ------------------------------------------------------------------

cells = [(mname, m9deg) for mname in P9_MODELS for m9deg in M9_GRID]
rows = []
failed = []
for k, ro, oo, et in sample:
    try:
        base = integrate(oo, et, None)
        if base is None:
            failed.append(k); continue
        po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
        pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
        rec = dict(desig=k, q=ro["q"], cls=ro["cls"],
                   theta=sep(-po, TNO), d_of=sep(-po, -pf),
                   aa_base=base["aa_out"], kick_e=oo["aa"] - base["aa_out"])
        for mname, m9deg in cells:
            model = P9_MODELS[mname]
            ps9, vs9 = p9_state(model, m9deg, et)
            rp = integrate(oo, et, (ps9, vs9, model["m9"] * M_EARTH_SUN))
            if rp is None:
                rec[f"{mname}_M{int(m9deg)}_drot"] = float("nan")
                rec[f"{mname}_M{int(m9deg)}_dk"] = float("nan")
                rec[f"{mname}_M{int(m9deg)}_denc"] = float("nan")
                continue
            rec[f"{mname}_M{int(m9deg)}_drot"] = sep(base["phat"], rp["phat"])
            rec[f"{mname}_M{int(m9deg)}_dk"] = abs(rp["aa_out"] - base["aa_out"])
            rec[f"{mname}_M{int(m9deg)}_denc"] = rp["denc_p9"]
        rows.append(rec)
    except Exception as e:
        failed.append(k); logger.warning(f"{k}: {e}")

logger.info(f"integrated {len(rows)} comets x {1 + len(cells)} runs; {len(failed)} failed")

# ------------------------------------------------------------------
# Statistics per (model, M9) cell and pooled
# ------------------------------------------------------------------

def cell_stats(key_rot, key_dk, sub):
    th = np.array([r["theta"] for r in sub]); v = np.array([r["d_of"] for r in sub])
    drot = np.array([r[key_rot] for r in sub]); dk = np.array([r[key_dk] for r in sub])
    K = np.abs(np.array([r["kick_e"] for r in sub])); Q = np.array([r["q"] for r in sub])
    ok = np.isfinite(drot)
    th, v, drot, dk, K, Q = th[ok], v[ok], drot[ok], dk[ok], K[ok], Q[ok]
    inc = th < 60
    out = {"n": int(len(v)), "n_in": int(inc.sum()),
           "med_drot": float(np.median(drot)),
           "med_drot_in": float(np.median(drot[inc])),
           "med_drot_out": float(np.median(drot[~inc])),
           "max_drot": float(np.max(drot)),
           "med_dk": float(np.median(dk)),
           "med_dk_in": float(np.median(dk[inc])),
           "obs_med_dof_in": float(np.median(v[inc])),
           "rotation_shortfall": float(np.median(v[inc]) / np.median(drot)) if np.median(drot) > 0 else float("inf"),
           "frac_drot_gt_dof": float(np.mean(drot > v))}
    u = mannwhitneyu(drot[inc], drot[~inc], alternative="greater")
    out["drot_in_gt_out_p"] = float(u.pvalue)
    rho, p = spearmanr(drot, v)
    out["spearman_drot_dof"] = {"rho": float(rho), "p_2sided": float(p)}
    y = np.log(v)
    X = np.column_stack([np.ones(len(y)), drot, K, Q])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    u2 = mannwhitneyu(resid[inc], resid[~inc], alternative="greater")
    out["residual_cap_after_p9"] = {"p": float(u2.pvalue)}
    return out

results = {"cells": {}, "pooled": {}}
for mname, m9deg in cells:
    kr, kd = f"{mname}_M{int(m9deg)}_drot", f"{mname}_M{int(m9deg)}_dk"
    results["cells"][f"{mname}_M{int(m9deg)}"] = {
        "model": {k2: P9_MODELS[mname][k2] for k2 in ("m9", "a9", "e9", "i9", "varpi9", "Om9")},
        "M9_deg": m9deg,
        "matched": cell_stats(kr, kd, [r for r in rows if r["q"] < 3.1]),
        "all_c1":  cell_stats(kr, kd, rows)}

# pooled: worst-case (max across cells) P9 rotation per comet
for rec in rows:
    rec["drot_max"] = float(np.nanmax([rec[f"{m}_M{int(g)}_drot"]
                                     for m, g in cells]))
    rec["dk_max"] = float(np.nanmax([rec[f"{m}_M{int(g)}_dk"]
                                   for m, g in cells]))
results["pooled"] = {
    "note": "per-comet maximum P9-attributable rotation across all "
            "(model, M9) cells -- the most generous bound",
    "matched": cell_stats("drot_max", "dk_max", [r for r in rows if r["q"] < 3.1]),
    "all_c1":  cell_stats("drot_max", "dk_max", rows)}

results["meta"] = {
    "n_all_c1": len(rows), "n_failed": len(failed), "failed": failed,
    "models": {m: {k2: P9_MODELS[m][k2] for k2 in
                   ("m9", "a9", "e9", "i9", "varpi9", "Om9", "note")}
               for m in P9_MODELS},
    "M9_grid_deg": M9_GRID,
    "integrator": "ias15", "boundary_AU": R_STOP,
    "reference": "Brown & Batygin 2021, AJ, 162, 219",
    "ephemeris": "data/raw/spice/de440s.bsp (provenance pinned)"}

out = str(RESULTS / "step_b27_planet_nine_insertion.json")
json.dump(results, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b27_planet_nine_insertion.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

for cell, d in results["cells"].items():
    a = d["all_c1"]
    logger.info(f"{cell}: med drot={a['med_drot']:.5f} deg "
                f"(in {a['med_drot_in']:.5f} / out {a['med_drot_out']:.5f}) "
                f"max {a['max_drot']:.4f}; shortfall x{a['rotation_shortfall']:.0f}; "
                f"resid cap p={a['residual_cap_after_p9']['p']:.4f}")
p = results["pooled"]["all_c1"]
logger.info(f"pooled max: med drot={p['med_drot']:.5f} deg, max {p['max_drot']:.4f}, "
            f"shortfall x{p['rotation_shortfall']:.0f}, resid p={p['residual_cap_after_p9']['p']:.4f}")
print("wrote", out)
print("wrote", csv_out)
