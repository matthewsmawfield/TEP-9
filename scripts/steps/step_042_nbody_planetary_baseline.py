"""step_19: N-body planetary baseline -- REBOUND integration of the CODE comets.

step_040 built the planetary baseline from an analytic Jupiter
encounter proxy.  This step replaces the proxy with a genuine
symplectic-family N-body integration:

  For each class-1 CODE comet the osculating elements at the
  perihelion epoch T (TT) are converted to a heliocentric state
  vector and integrated BACKWARD through the full planetary system
  (Sun + planetary-system barycentres 1-9, states from JPL DE440s
  rotated ICRF-equatorial -> ecliptic J2000) with REBOUND/IAS15.
  The integration stops when the comet reaches 255 AU, matching
  the 250 AU asymptote at which the catalogue's original orbit is
  defined.  A MERCURIUS (hybrid-symplectic, close-encounter
  switching) re-integration of the full sample provides an
  integrator-independence check.

  Kick metrics per comet:
    kick_e  = 1/a_osc - 1/a_out   [10^-6 AU^-1]  (energy kick;
              sign convention matches catalogue delta(1/a))
    denc    = minimum comet-planet distance over the integration
              [AU] and the dominant perturber

  Validation: the catalogue's own original-vs-osculating energy
  change cat_kick = 1/a_osc - 1/a_orig is an independent N-body
  measurement of the same planetary perturbation.  The simulated
  and catalogue kicks are compared statistically across the
  sample (deep-encounter comets are chaos-limited, so the
  comparison is reported both full-sample and trimmed).

  Baseline regression mirrors step_040:
    log(d_of) ~ 1 + |kick_e| + q  -> residuals -> cap contrast
    (Mann-Whitney, greater) and continuous Spearman test.

Outputs
-------
results/step_19_nbody_planetary_baseline.json
results/step_19_nbody_kicks.csv   (per-comet kick table)
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_042_nbody_planetary_baseline")
tee_stdout(logger)
logger.header("N-body planetary baseline (REBOUND/IAS15 + DE440s)")
import csv
import hashlib
import json
import math
import re
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
from html.parser import HTMLParser

import rebound
import spiceypy as sp

# ------------------------------------------------------------------
# CODE table parser (keeps the full perihelion-epoch string)
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
# Ephemeris setup: JPL DE440s, planetary-system barycentres
# ------------------------------------------------------------------

SPK = DATA_RAW / "spice" / "de440s.bsp"
if not SPK.exists():
    raise FileNotFoundError(f"JPL ephemeris missing: {SPK}")
sp.furnsh(str(SPK))

SPK_SHA = hashlib.sha256(open(SPK, "rb").read()).hexdigest()

AU_KM   = 149597870.7          # km / AU
DAY_YR  = 365.25               # d / yr (Julian year)
EPS     = math.radians(23.4392911)   # J2000 mean obliquity
RX      = np.array([[1, 0, 0],
                    [0, math.cos(EPS), math.sin(EPS)],
                    [0, -math.sin(EPS), math.cos(EPS)]])   # ICRF-equat -> ecl J2000

# DE440 system mass parameters, GM in km^3 s^-2 (planetary-system
# barycentres; Sun value is the DE440 point-mass GM).
GM_SUN = 1.32712440018e11
GM = {"1": 2.2031868551e4, "2": 3.2485859200e5, "3": 4.0350323562e5,
      "4": 4.2828375814e4, "5": 1.2671276480e8, "6": 3.7940626000e7,
      "7": 5.7945490100e6, "8": 6.8365271006e6, "9": 1.0868657e3}
BODY_NAME = {"1": "MeB", "2": "VB", "3": "EMB", "4": "MaB", "5": "JuB",
             "6": "SaB", "7": "UB", "8": "NB", "9": "PlB"}
PLANET_IDS = list(GM.keys())          # barycentre IDs 1..9

def jd_tt(y, m, d):
    """Julian date (TT scale) for calendar date y-m-d.ddd."""
    if m <= 2: y -= 1; m += 12
    A = y // 100; B = 2 - A + A // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5

def perihelion_et(T):
    """CODE 'T [TT]' epoch -> ephemeris seconds past J2000 (TT~TDB, <2 ms)."""
    m = re.match(r"(\d{4})\s+(\d{1,2})\s+(\d+\.?\d*)", T)
    if not m:
        return None
    return (jd_tt(int(m.group(1)), int(m.group(2)), float(m.group(3))) - 2451545.0) * 86400.0

def body_state(body, et):
    """SSB-barycentric state (J2000 ecliptic, AU & AU/yr) of SPK body."""
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return RX @ np.array(st[:3]) / AU_KM, RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR

# ------------------------------------------------------------------
# REBOUND integration
# ------------------------------------------------------------------

MU = 4 * math.pi ** 2            # G in AU^3 yr^-2 Msun^-1
R_STOP = 255.0                   # AU, heliocentric stop radius
T_MAX  = -20000.0                # yr, integration limit
DT_OUT = 1.0                     # yr, output cadence for min-distance tracking

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
    """Heliocentric ecliptic state vector at periapsis from osculating elements."""
    w_, O_, i_ = map(math.radians, (ro["w"], ro["Om"], ro["i"]))
    ph = perih_dir(w_, O_, i_)
    hh = np.array([math.sin(i_) * math.sin(O_), -math.sin(i_) * math.cos(O_), math.cos(i_)])
    th = np.cross(hh, ph); th /= np.linalg.norm(th)
    rvec = ro["q"] * ph
    vvec = math.sqrt(MU * (1 + ro["e"]) / ro["q"]) * th   # vis-viva at periapsis
    return rvec, vvec

def integrate(ro, et, integrator="ias15"):
    """Backward-integrate one comet to 255 AU.

    Returns dict(aa_out, kick_e, denc, dom, t_years) or None if the
    boundary is not reached.  aa_out and kick_e are barycentric
    1/a values in 10^-6 AU^-1.
    """
    sim, ps, vs = init_sim(et)
    rvec, vvec = state_at_periapsis(ro)
    sim.add(x=rvec[0] + ps[0], y=rvec[1] + ps[1], z=rvec[2] + ps[2],
            vx=vvec[0] + vs[0], vy=vvec[1] + vs[1], vz=vvec[2] + vs[2])
    nc = sim.N - 1
    if integrator == "mercurius":
        sim.integrator = "mercurius"   # hybrid symplectic w/ encounter switching
        sim.dt = 0.1                   # yr; r_crit_hill default (3 Hill radii)
    else:
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
    rc = np.array([p[nc].x, p[nc].y, p[nc].z])
    vc = np.array([p[nc].vx, p[nc].vy, p[nc].vz])
    E = vc.dot(vc) / 2 - MU * mtot / np.linalg.norm(rc)
    aa_out = -2 * E / (MU * mtot) * 1e6          # barycentric 1/a, 10^-6 AU^-1
    j = int(np.argmin(denc))
    return dict(aa_out=float(aa_out), kick_e=float(ro["aa"] - aa_out),
                denc=float(denc[j]), dom=BODY_NAME[PLANET_IDS[j]], t_years=float(t))

# ------------------------------------------------------------------
# Sample: identical selection to step_040 (class-1 CODE, non-Warsaw)
# ------------------------------------------------------------------

osc  = parse_code(str(DATA_RAW / "code" / "code_osculating.html"))
orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))
fut  = parse_code(str(DATA_RAW / "code" / "code_future.html"))
warsaw = {l[5:17].strip() for l in open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")) if len(l) > 115}

sample = []
skipped = {"no_osc": 0, "no_T": 0}
for k, ro in orig.items():
    if k not in fut: continue
    if not (0 < ro["aa"] < 100 and ro["cls"] in ("1a", "1a+", "1b")) or k in warsaw:
        continue
    if k not in osc:
        skipped["no_osc"] += 1; continue
    et = perihelion_et(osc[k]["T"])
    if et is None:
        skipped["no_T"] += 1; continue
    sample.append((k, ro, osc[k], et))

logger.info(f"sample: {len(sample)} class-1 CODE comets "
            f"(skipped {skipped['no_osc']} w/o osculating, {skipped['no_T']} w/o epoch)")

# ------------------------------------------------------------------
# Integrate
# ------------------------------------------------------------------

def _run_comet(k, ro, oo, et):
    r = integrate(oo, et, "ias15")
    rm = integrate(oo, et, "mercurius")
    if r is None:
        return None
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    return dict(desig=k, cls=ro["cls"], q=ro["q"], i=ro["i"], T=oo["T"],
        aa_osc=oo["aa"], aa_orig=ro["aa"], cat_kick=oo["aa"] - ro["aa"],
        aa_out=r["aa_out"], kick_e=r["kick_e"], denc=r["denc"], dom=r["dom"],
        t_years=r["t_years"],
        aa_out_merc=(rm["aa_out"] if rm else float("nan")),
        theta=sep(-po, TNO), d_of=sep(-po, -pf))


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
for (k, ro, oo, et), rec in zip(sample, recs):
    if rec is None:
        failed.append(k); logger.warning(f"{k}: boundary not reached")
    elif isinstance(rec, tuple) and rec[0] == "__error__":
        failed.append(rec[1]); logger.warning(f"{rec[1]}: integration failed ({rec[2]})")
    else:
        rows.append(rec)

logger.info(f"integrated {len(rows)} comets; {len(failed)} failed")

# ------------------------------------------------------------------
# Validation: simulated vs catalogue energy kick
# ------------------------------------------------------------------

cat = np.array([r["cat_kick"] for r in rows])
nb  = np.array([r["kick_e"] for r in rows])
diff = nb - cat
rho_v, p_v = spearmanr(nb, cat)
trim = np.abs(cat) < 400            # outside the deep-encounter chaotic regime
rho_t, p_t = spearmanr(nb[trim], cat[trim])
merc = np.array([r["aa_out_merc"] for r in rows])
ok_m = np.isfinite(merc)
rho_m, p_m = spearmanr(np.array([r["aa_out"] for r in rows])[ok_m], merc[ok_m])

validation = {
    "spearman_nb_vs_catalogue": {"rho": float(rho_v), "p": float(p_v)},
    "spearman_nb_vs_catalogue_trimmed": {"rho": float(rho_t), "p": float(p_t),
        "n": int(trim.sum()), "trim": "|cat_kick| < 400e-6 AU^-1"},
    "median_signed_diff": float(np.median(diff)),
    "median_abs_diff": float(np.median(np.abs(diff))),
    "abs_diff_iqr": [float(np.percentile(np.abs(diff), 25)),
                     float(np.percentile(np.abs(diff), 75))],
    "frac_within_50e-6": float(np.mean(np.abs(diff) < 50)),
    "mercurius_crosscheck": {"n": int(ok_m.sum()),
        "spearman_aa_out": {"rho": float(rho_m), "p": float(p_m)},
        "median_abs_aa_out_diff": float(np.median(np.abs(
            np.array([r["aa_out"] for r in rows])[ok_m] - merc[ok_m])))},
}

# ------------------------------------------------------------------
# Baseline regression on the N-body kick (mirrors step_040)
# ------------------------------------------------------------------

def run(sub):
    th = np.array([r["theta"] for r in sub]); v = np.array([r["d_of"] for r in sub])
    K = np.abs(np.array([r["kick_e"] for r in sub]))
    D = np.array([r["denc"] for r in sub])
    Q = np.array([r["q"] for r in sub])
    inc = th < 60
    out = {"n": len(sub), "n_in": int(inc.sum())}
    for key, x in [("kick_e", K), ("denc", D), ("q", Q), ("i", np.array([r["i"] for r in sub]))]:
        u = mannwhitneyu(x[inc], x[~inc], alternative="two-sided")
        out[f"{key}_p_2s"] = float(u.pvalue)
    rho, p = spearmanr(K, v); out["spearman_dof_kick_e"] = {"rho": float(rho), "p": float(p)}
    rho, p = spearmanr(D, v); out["spearman_dof_denc"] = {"rho": float(rho), "p": float(p)}
    y = np.log(v)
    X = np.column_stack([np.ones(len(y)), K, Q])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    u = mannwhitneyu(resid[inc], resid[~inc], alternative="greater")
    rho2, p2 = spearmanr(th, resid)
    out["residual_cap"] = {"med_in": float(np.median(resid[inc])),
                           "med_out": float(np.median(resid[~inc])), "p": float(u.pvalue)}
    out["residual_continuous"] = {"rho": float(rho2), "p_2sided": float(p2)}
    # robustness: two-covariate baseline (energy kick + closest approach)
    X2 = np.column_stack([np.ones(len(y)), K, D, Q])
    coef2, *_ = np.linalg.lstsq(X2, y, rcond=None)
    resid2 = y - X2 @ coef2
    u2 = mannwhitneyu(resid2[inc], resid2[~inc], alternative="greater")
    rho3, p3 = spearmanr(th, resid2)
    out["residual_cap_2cov"] = {"p": float(u2.pvalue)}
    out["residual_continuous_2cov"] = {"rho": float(rho3), "p_2sided": float(p3)}
    return out

res = {
    "method": "REBOUND N-body: osculating elements at CODE perihelion epoch, "
              "backward integration to 255 AU under Sun + planetary-system "
              "barycentres (DE440s), IAS15 primary + MERCURIUS cross-check",
    "n_all_c1": len(rows), "n_failed": len(failed), "failed": failed,
    "ephemeris": {"file": "data/raw/spice/de440s.bsp", "sha256": SPK_SHA,
                  "frame": "SSB barycentric ICRF -> ecliptic J2000",
                  "planets": "system barycentres 1-9 + Sun"},
    "integrator": {"primary": "ias15", "crosscheck": "mercurius",
                   "boundary_AU": R_STOP, "G": "4*pi^2 (AU, yr, Msun)"},
    "validation": validation,
    "matched": run([r for r in rows if r["q"] < 3.1]),
    "all_c1":  run(rows),
}

out = str(RESULTS / "step_19_nbody_planetary_baseline.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_19_nbody_kicks.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)
logger.data_save(csv_out)
# ------------------------------------------------------------------
# Figure: periapsis-direction discrepancy vs measured energy kick
# ------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

th  = np.array([r["theta"]  for r in rows])
dof = np.array([r["d_of"]   for r in rows])
kk  = np.abs(np.array([r["kick_e"] for r in rows]))
inc = th < 60

fig, ax = plt.subplots(figsize=(7.2, 4.8))
ax.scatter(kk[~inc], dof[~inc], s=22, c="0.55", alpha=0.75,
           label=f"outside 60° cap ($n={int((~inc).sum())}$)")
ax.scatter(kk[inc], dof[inc], s=26, c="crimson", alpha=0.85,
           label=f"inside 60° cap ($n={int(inc.sum())}$)")
ax.axhline(np.median(dof[inc]), color="crimson", lw=1, ls="--", alpha=0.6)
ax.axhline(np.median(dof[~inc]), color="0.4", lw=1, ls="--", alpha=0.6)
ax.set_xscale("log")
ax.set_xlabel(r"N-body energy kick $|\Delta(1/a)|$  ($10^{-6}$ AU$^{-1}$)")
ax.set_ylabel(r"periapsis-direction discrepancy $d_{\rm of}$ (deg)")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_19_rotation_vs_energy.png", dpi=300)
logger.data_save(FIG / "supplementary" / "step_19_rotation_vs_energy.png")