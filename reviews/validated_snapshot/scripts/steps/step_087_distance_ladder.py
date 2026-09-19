#!/usr/bin/env python3
"""Step 087 -- perturber distance ladder: what mass is required at
any distance?

Step 083 fixed the perturber's geometry to the published
Brown & Batygin (2021) median location and showed the required mass
is brown-dwarf scale (~4300 M_earth in-cap median).  That leaves one
escape route: a perturber elsewhere -- closer, where a smaller mass
delivers the same impulse, or farther, where surveys are shallower.
This step closes it by measuring the required mass as a function of
perturber semimajor axis, holding the anti-aligned shepherding
geometry (e9 = 0.30, i9 = 15.6 deg, varpi9 = 246.7 deg, Om9 = 96.9
deg, M9 = 180 deg residence-weighted) and the fiducial 6.9 M_earth
mass.  Because the injected rotation scales linearly in mass
(step 083, alpha = 1.00), a single rung at each distance yields the
required mass directly:

    m_req(a9) = 6.9 * (d_of / drot_6.9(a9))

Rungs: a9 = 250, 300, 461, 700, 1000 AU.  Deliverables: the
required-mass distribution per rung in and out of the cap, the
distance-scaling exponent gamma of m_req ~ a9^gamma, and the
comparison against survey exclusion lines (IRAS/WISE all-sky limits
of order Jupiter masses at these distances; BB21 nominal 5-10
M_earth).  A point mass that must sit at any of these distances at
~10^3 M_earth is excluded by every relevant survey; a perturber
closer still enters the planetary-ephemeris bounded regime.

Outputs
-------
results/step_b52_distance_ladder.json
results/step_b52_distance_ladder.csv   (per-comet ladder)
results/figures/step_b52_distance_ladder.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_087_distance_ladder")
tee_stdout(logger)
logger.header("Perturber distance ladder (REBOUND/IAS15 + DE440s)")

import csv
import json
import math
import re
import numpy as np
from scipy.stats import spearmanr
from html.parser import HTMLParser

import rebound
import spiceypy as sp

# ------------------------------------------------------------------
# CODE table parser (identical to step_062/083)
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

TNO = lv(34, -13)   # transit-axis convention (Phase-3 steps)

# ------------------------------------------------------------------
# Ephemeris setup (identical to step_062/083)
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
# Perturber distance rungs (BB21 median geometry, m9 = 6.9)
# ------------------------------------------------------------------

M_EARTH_SUN = 3.0035e-6
M9_EARTH = 6.9
A9_LADDER = [250.0, 300.0, 461.0, 700.0, 1000.0]
GEO = dict(e9=0.30, i9=15.6, varpi9=246.7, Om9=96.9)
M9_USE = 180.0

MU = 4 * math.pi ** 2
R_STOP = 255.0
T_MAX  = -20000.0
DT_OUT = 1.0

def p9_state(a9, m9_deg, et):
    w9 = GEO["varpi9"] - GEO["Om9"]
    P9_yr = a9 ** 1.5
    et_yr = et / (86400.0 * DAY_YR)
    M_at_et = math.radians(m9_deg) + (2 * math.pi / P9_yr) * et_yr
    sim = rebound.Simulation(); sim.G = MU
    sim.add(m=1.0)
    sim.add(a=a9, e=GEO["e9"], inc=math.radians(GEO["i9"]),
            Omega=math.radians(GEO["Om9"]), omega=math.radians(w9),
            M=M_at_et, m=M9_EARTH * M_EARTH_SUN)
    p = sim.particles[1]
    return (np.array([p.x, p.y, p.z]), np.array([p.vx, p.vy, p.vz]))

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
    aa = -2 * E / mu * 1e6
    return phat, aa

def integrate(ro, et, p9=None):
    sim, ps, vs = init_sim(et)
    if p9 is not None:
        sim.add(x=p9[0][0], y=p9[0][1], z=p9[0][2],
                vx=p9[1][0], vy=p9[1][1], vz=p9[1][2], m=p9[2])
    rvec, vvec = state_at_periapsis(ro)
    sim.add(x=rvec[0] + ps[0], y=rvec[1] + ps[1], z=rvec[2] + ps[2],
            vx=vvec[0] + vs[0], vy=vvec[1] + vs[1], vz=vvec[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    t = -DT_OUT
    # t is negative for the backward leg; compare elapsed time, not signed t.
    while abs(t) < abs(T_MAX):
        sim.integrate(t, exact_finish_time=0)
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y, p[nc].z - p[0].z])
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
    return dict(phat=phat, aa_out=aa_out)

# ------------------------------------------------------------------
# Sample: identical selection to step_062/083
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
# Integrate: baseline + one rung per distance per comet
# ------------------------------------------------------------------

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
                   theta=sep(-po, TNO), d_of=sep(-po, -pf))
        for a9 in A9_LADDER:
            key = f"a{int(a9)}"
            ps9, vs9 = p9_state(a9, M9_USE, et)
            rp = integrate(oo, et, (ps9, vs9, M9_EARTH * M_EARTH_SUN))
            rec[f"{key}_drot"] = sep(base["phat"], rp["phat"]) if rp else float("nan")
            rec[f"{key}_mreq"] = (M9_EARTH * rec["d_of"] / rec[f"{key}_drot"]
                                  if np.isfinite(rec[f"{key}_drot"]) and rec[f"{key}_drot"] > 0
                                  else float("nan"))
        rows.append(rec)
    except Exception as e:
        failed.append(k); logger.warning(f"{k}: {e}")

logger.info(f"integrated {len(rows)} comets x {1 + len(A9_LADDER)} runs; {len(failed)} failed")

# ------------------------------------------------------------------
# Required mass vs distance
# ------------------------------------------------------------------

inc = np.array([r["theta"] < 60 for r in rows])
dof = np.array([r["d_of"] for r in rows])

res = {"meta": {
    "n": len(rows), "n_failed": len(failed), "failed": failed,
    "geometry": {"e9": GEO["e9"], "i9": GEO["i9"],
                 "varpi9": GEO["varpi9"], "Om9": GEO["Om9"],
                 "M9_deg": M9_USE, "m9_earth": M9_EARTH,
                 "note": "BB21 median apsidal geometry; residence-"
                         "weighted most probable phase"},
    "a9_ladder_AU": A9_LADDER,
    "anchor_scaling": "alpha = 1.00 measured on the mass ladder "
                      "(step 083) -- m_req(a9) = 6.9*d_of/drot_6.9(a9)",
    "integrator": "ias15", "boundary_AU": R_STOP,
    "ephemeris": "data/raw/spice/de440s.bsp (provenance pinned)"},
    "rungs": {}}

for a9 in A9_LADDER:
    key = f"a{int(a9)}"
    dr = np.array([r[f"{key}_drot"] for r in rows])
    mr = np.array([r[f"{key}_mreq"] for r in rows])
    ok = np.isfinite(mr)
    res["rungs"][key] = {
        "a9_AU": a9,
        "med_drot": float(np.nanmedian(dr)),
        "m_req_median_in": float(np.nanmedian(mr[ok & inc])),
        "m_req_p05_in": float(np.nanpercentile(mr[ok & inc], 5)),
        "m_req_median_out": float(np.nanmedian(mr[ok & ~inc])),
        "frac_in_req_gt_1000": float(np.mean(mr[ok & inc] > 1000)),
        "frac_in_req_gt_100": float(np.mean(mr[ok & inc] > 100)),
        "frac_in_req_gt_318": float(np.mean(mr[ok & inc] > 318)),
        "note": "frac>318 = fraction of in-cap comets requiring more "
                "than a Jupiter mass at this distance -- the IRAS/WISE "
                "all-sky exclusion scale"}
    a = res["rungs"][key]
    logger.info(f"a9={a9:5.0f} AU: med drot={a['med_drot']:.5f} deg | "
                f"required mass in-cap med {a['m_req_median_in']:.0f} "
                f"M_earth, p05 {a['m_req_p05_in']:.0f}, "
                f">M_Jup {a['frac_in_req_gt_318']:.2f}")

# distance-scaling exponents on the in-cap medians
a9v = np.array(A9_LADDER)
mr_in = np.array([res["rungs"][f"a{int(a9)}"]["m_req_median_in"] for a9 in A9_LADDER])
dr_med = np.array([res["rungs"][f"a{int(a9)}"]["med_drot"] for a9 in A9_LADDER])
ok = np.isfinite(mr_in) & (mr_in > 0)
if ok.sum() >= 3:
    g = np.polyfit(np.log(a9v[ok]), np.log(mr_in[ok]), 1)
    res["distance_scaling"] = {
        "gamma_mreq_vs_a9": float(g[0]),
        "beta_drot_vs_a9": float(np.polyfit(np.log(a9v), np.log(dr_med), 1)[0]),
        "note": "gamma: required mass ~ a9^gamma on in-cap medians; "
                "beta: injected rotation ~ a9^beta on the 6.9 "
                "M_earth rungs; gamma ~ -beta by construction under "
                "linear mass scaling"}
    logger.info(f"distance scaling: m_req ~ a9^{g[0]:+.2f}; "
                f"drot ~ a9^{res['distance_scaling']['beta_drot_vs_a9']:+.2f}")

# ------------------------------------------------------------------
# Write
# ------------------------------------------------------------------

out = str(RESULTS / "step_b52_distance_ladder.json")
json.dump(res, open(out, "w"), indent=1, default=float)
csv_out = str(RESULTS / "step_b52_distance_ladder.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("wrote", out)
print("wrote", csv_out)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

ax = axes[0]
med_in = [res["rungs"][f"a{int(a)}"]["m_req_median_in"] for a in A9_LADDER]
p05_in = [res["rungs"][f"a{int(a)}"]["m_req_p05_in"] for a in A9_LADDER]
med_out = [res["rungs"][f"a{int(a)}"]["m_req_median_out"] for a in A9_LADDER]
ax.plot(A9_LADDER, med_in, "o-", c="crimson", label="in-cap median")
ax.plot(A9_LADDER, p05_in, "s--", c="crimson", alpha=0.5, label="in-cap 5th pct")
ax.plot(A9_LADDER, med_out, "o-", c="0.6", label="out-of-cap median")
ax.axhline(318, color="steelblue", ls=":", lw=1.2, label="1 $M_{Jup}$")
ax.axhline(6.9, color="k", ls="--", lw=1, label="BB21 nominal mass")
ax.axhspan(318, 1e5, color="steelblue", alpha=0.06)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("perturber semimajor axis $a_9$ (AU)")
ax.set_ylabel("required mass ($M_\\oplus$)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("required point mass vs perturber distance", fontsize=10)

ax = axes[1]
dr_in = []
for a9 in A9_LADDER:
    key = f"a{int(a9)}"
    dr_in.append(np.nanmedian([r[key + "_drot"] for r in rows if r["theta"] < 60]))
ax.plot(A9_LADDER, dr_in, "o-", c="teal", label="injected rotation (6.9 $M_\\oplus$)")
ax.axhline(float(np.median(dof[inc])), color="crimson", ls="--", lw=1.2,
           label="observed in-cap median $\\delta\\theta$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("perturber semimajor axis $a_9$ (AU)")
ax.set_ylabel("boundary rotation (deg)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("what 6.9 $M_\\oplus$ can do vs what is observed", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b52_distance_ladder.png", dpi=150)
print(f"wrote {FIG / 'step_b52_distance_ladder.png'}")
