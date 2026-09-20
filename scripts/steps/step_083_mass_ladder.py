#!/usr/bin/env python3
"""Step 083 -- perturber mass ladder: what point mass would the
transit anomaly require?

step_062 showed that the two published Brown & Batygin (2021)
realizations inject <= 0.003 deg of boundary rotation against the
>= 0.14 deg observed -- a shortfall of roughly two orders of
magnitude in angle.  This step converts that shortfall into a mass
statement.  The gravitational impulse a distant perturber delivers to
a comet's reconstructed orbit scales with the perturber's mass, so a
mass ladder at fixed geometry measures both the scaling exponent and
the mass a point-mass explanation would actually require.

Cells
-----
bb21_med_m7    m9 =   6.9 M_earth (the published realization, M9 = 180)
bb21_med_m20   m9 =  20   "
bb21_med_m50   m9 =  50   "
bb21_med_m150  m9 = 150   "
sct25_px       m9 =   4.4 M_earth, a = 290 AU, e = 0.29, i = 6.8 deg --
               the Siraj, Chyba & Tremaine (2025, ApJ 978, 139)
               quadrivariate best fit, placed in the anti-aligned
               shepherding configuration (varpi9 = 246.7, Om9 = 96.9,
               the same apsidal geometry as bb21_med; the fit
               constrains m, a, e, i only) at M9 = 180.

All cells use the residence-time-weighted most probable position
M9 = 180 deg.  The m9 = 6.9 anchor is taken from the step_b27 CSV
(deterministic integrator, identical inputs) so each comet needs one
no-P9 baseline plus four new integrations.

Per comet the periapsis-direction rotation drot is measured at each
mass; a log-log fit on the bb21_med ladder gives the scaling exponent
alpha (drot ~ m^alpha), and the required mass

    m_req = 6.9 * (d_of / drot_6.9)^(1/alpha)

is the point mass at the published location that would reproduce that
comet's observed orig->future rotation d_of.  Deliverables: the
scaling exponent, the required-mass distribution in and out of the
cap, and the same targeting/absorption statistics as step_062 for the
Siraj realization.

Context for interpretation (recorded in the JSON): Pan-STARRS1 +
ZTF + DES already exclude ~78 per cent of the BB21 parameter space
(Brown, Holman & Batygin 2024, AJ 167, 146); IRAS/WISE all-sky limits bound any
roughly Jupiter-mass-or-larger body well inside ~10^4 AU; and
planetary-ranging ephemerides bound a ~300-500 AU perturber at a few
tens of Earth masses.  A required mass of order 10^3-10^4 M_earth is
not a planet -- it is a brown dwarf that every relevant survey
excludes.

Outputs
-------
results/step_b48_mass_ladder.json
results/step_b48_mass_ladder.csv   (per-comet ladder)
results/figures/supplementary/step_b48_mass_ladder.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_083_mass_ladder")
tee_stdout(logger)
logger.header("Perturber mass ladder (REBOUND/IAS15 + DE440s)")

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
# CODE table parser (identical to step_062)
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
# Ephemeris setup (identical to step_062)
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
# Perturber realizations
# ------------------------------------------------------------------

M_EARTH_SUN = 3.0035e-6

P9_MODELS = {
    "bb21_med_m20":  dict(m9=20.0,  a9=461.0, e9=0.30, i9=15.6,
                          varpi9=246.7, Om9=96.9,
                          note="BB21 median geometry, mass ladder rung"),
    "bb21_med_m50":  dict(m9=50.0,  a9=461.0, e9=0.30, i9=15.6,
                          varpi9=246.7, Om9=96.9,
                          note="BB21 median geometry, mass ladder rung"),
    "bb21_med_m150": dict(m9=150.0, a9=461.0, e9=0.30, i9=15.6,
                          varpi9=246.7, Om9=96.9,
                          note="BB21 median geometry, mass ladder rung"),
    "sct25_px":      dict(m9=4.4,   a9=290.0, e9=0.29, i9=6.8,
                          varpi9=246.7, Om9=96.9,
                          note="Siraj, Chyba & Tremaine 2025 ApJ 978 139 "
                               "quadrivariate best fit (m,a,e,i); apsidal "
                               "geometry set to the anti-aligned "
                               "shepherding configuration"),
}
M9_USE = 180.0     # residence-time-weighted most probable position

MU = 4 * math.pi ** 2
R_STOP = 255.0
T_MAX  = -20000.0
DT_OUT = 1.0

def p9_state(model, m9_deg, et):
    w9 = model["varpi9"] - model["Om9"]
    P9_yr = model["a9"] ** 1.5
    et_yr = et / (86400.0 * DAY_YR)
    M_at_et = math.radians(m9_deg) + (2 * math.pi / P9_yr) * et_yr
    sim = rebound.Simulation(); sim.G = MU
    sim.add(m=1.0)
    sim.add(a=model["a9"], e=model["e9"], inc=math.radians(model["i9"]),
            Omega=math.radians(model["Om9"]), omega=math.radians(w9),
            M=M_at_et, m=model["m9"] * M_EARTH_SUN)
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
# Sample: identical selection to step_062
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

# anchor rotations at m9 = 6.9 (bb21_med) from the step_b27 CSV --
# deterministic integrator, identical inputs.  Two anchors: the
# residence-weighted M9 = 180 column used throughout the ladder, and
# the per-comet maximum across all four M9 placements -- the most
# favourable perturber position, used for the conservative bound.
anchor = {}
with open(RESULTS / "step_b27_planet_nine_insertion.csv") as f:
    for r in csv.DictReader(f):
        anchor[r["desig"]] = dict(
            drot_6p9=float(r["bb21_med_M180_drot"]),
            drot_gen=max(float(r[f"bb21_med_M{g}_drot"])
                         for g in (0, 90, 180, 270)),
            dk_6p9=float(r["bb21_med_M180_dk"]))

# ------------------------------------------------------------------
# Integrate: baseline + 4 new cells per comet
# ------------------------------------------------------------------

cells = list(P9_MODELS.keys())


def _run_comet(k, ro, oo, et):
    base = integrate(oo, et, None)
    if base is None:
        return None
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    rec = dict(desig=k, q=ro["q"], cls=ro["cls"],
               theta=sep(-po, TNO), d_of=sep(-po, -pf),
               drot_6p9=anchor.get(k, {}).get("drot_6p9", float("nan")),
               drot_gen=anchor.get(k, {}).get("drot_gen", float("nan")))
    for mname in cells:
        model = P9_MODELS[mname]
        ps9, vs9 = p9_state(model, M9_USE, et)
        rp = integrate(oo, et, (ps9, vs9, model["m9"] * M_EARTH_SUN))
        if rp is None:
            rec[f"{mname}_drot"] = float("nan")
            rec[f"{mname}_denc"] = float("nan")
            continue
        rec[f"{mname}_drot"] = sep(base["phat"], rp["phat"])
        rec[f"{mname}_denc"] = rp["denc_p9"]
    return rec


def _work_comet(job):
    k, ro, oo, et = job
    try:
        return _run_comet(k, ro, oo, et)
    except Exception as e:
        return ("__error__", k, str(e))


from scripts.utils.parallel import default_workers as _default_workers

if "--workers" in _sys.argv:
    N_WORK = max(1, int(_sys.argv[_sys.argv.index("--workers") + 1]))
else:
    N_WORK = _default_workers()

if N_WORK > 1 and len(sample) > 1:
    import multiprocessing as mp
    ctx = mp.get_context("fork") if _sys.platform != "win32" \
        else mp.get_context("spawn")

    def _init():
        # forked children inherit the parent's BSP fd; concurrent spkezr
        # reads through a shared descriptor corrupt each other -- reopen
        # the kernel so each worker holds its own file handle.
        sp.kclear()
        sp.furnsh(str(SPK))

    with ctx.Pool(min(N_WORK, len(sample)), initializer=_init) as pool:
        recs = pool.map(_work_comet, sample)
else:
    recs = [_work_comet(s) for s in sample]

rows = []
failed = []
for s, rec in zip(sample, recs):
    if rec is None:
        failed.append(s[0])
    elif isinstance(rec, tuple) and rec[0] == "__error__":
        failed.append(rec[1]); logger.warning(f"{rec[1]}: {rec[2]}")
    else:
        rows.append(rec)

logger.info(f"integrated {len(rows)} comets x {1 + len(cells)} runs; {len(failed)} failed")

# ------------------------------------------------------------------
# Scaling exponent and required mass
# ------------------------------------------------------------------

LADDER_M = np.array([6.9, 20.0, 50.0, 150.0])
def ladder_vals(rec):
    return np.array([rec["drot_6p9"], rec["bb21_med_m20_drot"],
                     rec["bb21_med_m50_drot"], rec["bb21_med_m150_drot"]])

for rec in rows:
    v = ladder_vals(rec)
    ok = np.isfinite(v) & (v > 0)
    if ok.sum() >= 3:
        a_fit = np.polyfit(np.log(LADDER_M[ok]), np.log(v[ok]), 1)
        rec["alpha"] = float(a_fit[0])
        # required mass to reach the observed rotation
        if np.isfinite(rec["drot_6p9"]) and rec["drot_6p9"] > 0:
            rec["m_req"] = float(6.9 * (rec["d_of"] / rec["drot_6p9"]) ** (1.0 / a_fit[0]))
        else:
            rec["m_req"] = float("nan")
        # conservative bound: most favourable M9 placement
        if np.isfinite(rec["drot_gen"]) and rec["drot_gen"] > 0:
            rec["m_req_gen"] = float(6.9 * (rec["d_of"] / rec["drot_gen"]) ** (1.0 / a_fit[0]))
        else:
            rec["m_req_gen"] = float("nan")
    else:
        rec["alpha"] = float("nan"); rec["m_req"] = float("nan")
        rec["m_req_gen"] = float("nan")

inc = np.array([r["theta"] < 60 for r in rows])
mr  = np.array([r["m_req"] for r in rows])
mrg = np.array([r["m_req_gen"] for r in rows])
al  = np.array([r["alpha"] for r in rows])
dof = np.array([r["d_of"] for r in rows])

okm = np.isfinite(mr)
okg = np.isfinite(mrg)
res = {"meta": {
    "n": len(rows), "n_failed": len(failed), "failed": failed,
    "models": {m: {k2: P9_MODELS[m][k2] for k2 in
                   ("m9","a9","e9","i9","varpi9","Om9","note")}
               for m in P9_MODELS},
    "anchor": "bb21_med m9=6.9 M9=180 drot taken from "
              "step_b27_planet_nine_insertion.csv (deterministic run)",
    "M9_deg": M9_USE, "integrator": "ias15", "boundary_AU": R_STOP,
    "references": ["Brown & Batygin 2021, AJ, 162, 219",
                   "Siraj, Chyba & Tremaine 2025, ApJ, 978, 139",
                   "Brown, Holman & Batygin 2024, AJ, 167, 146 (PS1: ~78% of BB21 "
                   "space excluded by PS1+ZTF+DES)",
                   "Phan et al. 2025, PASA (IRAS+AKARI far-IR search)"],
    "ephemeris": "data/raw/spice/de440s.bsp (provenance pinned)"}}

res["scaling"] = {
    "alpha_median": float(np.nanmedian(al)),
    "alpha_i16": float(np.nanpercentile(al, 16)),
    "alpha_i84": float(np.nanpercentile(al, 84)),
    "note": "per-comet log-log slope of boundary rotation vs perturber "
            "mass across the bb21_med ladder (4 rungs); ~1 = linear "
            "impulse scaling"}
logger.info(f"scaling exponent alpha = {res['scaling']['alpha_median']:.2f} "
            f"[{res['scaling']['alpha_i16']:.2f}, {res['scaling']['alpha_i84']:.2f}]")

res["required_mass"] = {
    "n": int(okm.sum()),
    "m_req_median_all": float(np.nanmedian(mr)),
    "m_req_median_in": float(np.nanmedian(mr[okm & inc])),
    "m_req_median_out": float(np.nanmedian(mr[okm & ~inc])),
    "m_req_p05_in": float(np.nanpercentile(mr[okm & inc], 5)),
    "m_req_min": float(np.nanmin(mr)),
    "frac_in_req_gt_1000": float(np.mean(mr[okm & inc] > 1000)),
    "frac_in_req_gt_100": float(np.mean(mr[okm & inc] > 100)),
    "m_req_gen_median_in": float(np.nanmedian(mrg[okg & inc])),
    "m_req_gen_p05_in": float(np.nanpercentile(mrg[okg & inc], 5)),
    "frac_in_gen_gt_1000": float(np.mean(mrg[okg & inc] > 1000)),
    "frac_in_gen_gt_100": float(np.mean(mrg[okg & inc] > 100)),
    "note": "m_req [M_earth] = 6.9 * (d_of/drot)^(1/alpha): the point "
            "mass at the BB21 median location that would reproduce each "
            "comet's observed rotation; m_req uses the residence-"
            "weighted M9=180 anchor, m_req_gen the most favourable of "
            "the four M9 placements -- the conservative bound"}
logger.info(f"required mass: median {res['required_mass']['m_req_median_all']:.0f} M_earth; "
            f"in-cap {res['required_mass']['m_req_median_in']:.0f}; "
            f"5th pct {res['required_mass']['m_req_p05_in']:.0f}; "
            f"frac>1000 {res['required_mass']['frac_in_req_gt_1000']:.2f}")

# ------------------------------------------------------------------
# Siraj realization: targeting + absorption statistics (step_062 style)
# ------------------------------------------------------------------

def cell_stats(key_rot, sub):
    th = np.array([r["theta"] for r in sub]); v = np.array([r["d_of"] for r in sub])
    drot = np.array([r[key_rot] for r in sub])
    ok = np.isfinite(drot)
    th, v, drot = th[ok], v[ok], drot[ok]
    inc_ = th < 60
    out = {"n": int(len(v)), "n_in": int(inc_.sum()),
           "med_drot": float(np.median(drot)),
           "med_drot_in": float(np.median(drot[inc_])) if inc_.sum() else float("nan"),
           "med_drot_out": float(np.median(drot[~inc_])) if (~inc_).sum() else float("nan"),
           "max_drot": float(np.max(drot)),
           "obs_med_dof_in": float(np.median(v[inc_])),
           "rotation_shortfall": float(np.median(v[inc_]) / np.median(drot)) if np.median(drot) > 0 else float("inf")}
    u = mannwhitneyu(drot[inc_], drot[~inc_], alternative="greater")
    out["drot_in_gt_out_p"] = float(u.pvalue)
    rho, p = spearmanr(drot, v)
    out["spearman_drot_dof"] = {"rho": float(rho), "p_2sided": float(p)}
    return out

res["cells"] = {}
for mname in cells:
    res["cells"][mname] = {
        "m9": P9_MODELS[mname]["m9"],
        "all_c1": cell_stats(f"{mname}_drot", rows)}
    a = res["cells"][mname]["all_c1"]
    logger.info(f"{mname}: med drot={a['med_drot']:.5f} deg "
                f"(in {a['med_drot_in']:.5f} / out {a['med_drot_out']:.5f}) "
                f"max {a['max_drot']:.4f}; shortfall x{a['rotation_shortfall']:.0f}")

# ------------------------------------------------------------------
# Write
# ------------------------------------------------------------------

out = str(RESULTS / "step_b48_mass_ladder.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b48_mass_ladder.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
logger.data_save(out)
logger.data_save(csv_out)
# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

ax = axes[0]
for rec in rows:
    v = ladder_vals(rec)
    ok = np.isfinite(v) & (v > 0)
    col = "crimson" if rec["theta"] < 60 else "0.75"
    ax.plot(LADDER_M[ok], v[ok], "-", color=col, alpha=0.25, lw=0.8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("perturber mass (M_earth)")
ax.set_ylabel("boundary periapsis rotation (deg)")
ax.axhline(np.median(dof[inc]), color="crimson", ls="--", lw=1.2,
           label=f"observed in-cap median {np.median(dof[inc]):.2f} deg")
ax.legend(frameon=False, fontsize=8)
ax.set_title("mass ladder at the BB21 median location", fontsize=10)

ax = axes[1]
mrv = mr[okm]
ax.hist(np.log10(mrv[~inc[okm]]), bins=20, color="0.7", alpha=0.8,
        label="out-of-cap", density=True)
ax.hist(np.log10(mrv[inc[okm]]), bins=20, color="crimson", alpha=0.6,
        label="in-cap", density=True)
for mlim, lab in ((318, "Jupiter"), (6.9, "BB21 med"), (4.4, "SCT25")):
    ax.axvline(math.log10(mlim), color="k", ls=":", lw=1)
    ax.text(math.log10(mlim), ax.get_ylim()[1]*0.9, " " + lab, fontsize=7,
            rotation=90, va="top")
ax.set_xlabel("log10 required mass (M_earth)")
ax.set_ylabel("density")
ax.legend(frameon=False, fontsize=8)
ax.set_title("point mass required to match the anomaly", fontsize=10)

ax = axes[2]
drot6 = np.array([r["drot_6p9"] for r in rows])
ax.scatter(dof[~inc], drot6[~inc], s=11, c="0.6", alpha=0.7, label="out-of-cap")
ax.scatter(dof[inc], drot6[inc], s=11, c="crimson", alpha=0.75, label="in-cap")
ax.plot([0, dof.max()], [0, dof.max()], "k:", lw=1, label="drot = d_of")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("observed rotation $d_{of}$ (deg)")
ax.set_ylabel("P9-injected rotation (deg)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("targeting: the perturber misses the discrepant comets", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b48_mass_ladder.png", dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b48_mass_ladder.png')