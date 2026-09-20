"""step_063: Bidirectional boundary validation + rotation-per-kick residual.

Every Phase-3 step measures the catalogue's orig->fut periapsis-direction
discrepancy d_of.  step_042 verified that a real N-body integration
reproduces the catalogue's energy kicks (rho=0.98), and step_062 showed
the published Planet Nine cannot inject the observed rotation.  This
step measures what STANDARD dynamics predicts for the rotation itself,
and -- critically -- evaluates the boundary solutions in the same
BARYCENTRIC frame the catalogue uses for its original/future orbits.

For each class-1 CODE comet (identical selection to step_040/042/062)
the osculating elements at the catalogued perihelion epoch are
converted to a state and integrated through Sun + planetary-system
barycentres (DE440s, REBOUND/IAS15) TWICE:

  backward until barycentric |r| = 250 AU  -> analogue of the catalogue's
                                             "original" solution
  forward  until barycentric |r| = 250 AU  -> analogue of the "future"
                                             solution

Outputs
-------
results/step_b28_bidirectional_rotation.json
results/step_b28_bidirectional_rotation.csv
results/figures/supplementary/step_b28_rotation_residual.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
from scripts.utils.coordinates import angular_separation
logger = StepLogger("step_063_bidirectional_rotation")
tee_stdout(logger)
logger.header("Bidirectional boundary validation (REBOUND/IAS15 + DE440s)")

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
# CODE table parser (identical to step_042/062)
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
    return angular_separation(a, b)

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO = lv(34, -13)

# ------------------------------------------------------------------
# Ephemeris setup (identical to step_042/062)
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
# REBOUND integration -- both legs, barycentric boundary elements
# ------------------------------------------------------------------

MU = 4 * math.pi ** 2
R_STOP = 250.0              # barycentric |r|, matching the catalogue sphere
T_MAX  = 20000.0            # |yr| integration limit, either direction
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
    mu = MU * mtot
    r = np.linalg.norm(r_rel)
    h = np.cross(r_rel, v_rel)
    evec = (np.cross(v_rel, h) / mu) - r_rel / r
    en = np.linalg.norm(evec)
    phat = evec / en if en > 1e-12 else r_rel / r
    E = v_rel.dot(v_rel) / 2 - mu / r
    return phat, -2 * E / mu * 1e6

def integrate_leg(ro, et, direction):
    """Integrate one comet to the 250 AU sphere.  direction = -1 backward
    (pre-perihelion asymptote, analogue of the catalogue's "original"
    solution), +1 forward ("future").

    Boundary elements are BARYCENTRIC: the comet state is evaluated
    relative to the Sun + planet-system barycentre, matching the
    catalogue's original/future convention.

    Returns dict(phat, aa, denc, t_years) or None if not reached."""
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
    # barycentric osculating elements at the boundary
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

# ------------------------------------------------------------------
# Sample: identical selection to step_040/042/062
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
# Integrate both legs
# ------------------------------------------------------------------

def _run_comet(k, ro, oo, et):
    rb = integrate_leg(oo, et, -1)
    rf = integrate_leg(oo, et, +1)
    if rb is None or rf is None:
        return None
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    drot_cat = sep(-po, -pf)
    drot_sim = sep(rb["phat"], rf["phat"])
    return dict(desig=k, cls=ro["cls"], q=ro["q"], i=ro["i"],
        theta=sep(-po, TNO),
        drot_cat=drot_cat, drot_sim=drot_sim,
        d_back=sep(rb["phat"], po),      # sim vs catalogue, original leg
        d_fwd=sep(rf["phat"], pf),       # sim vs catalogue, future leg
        daa_cat=fut[k]["aa"] - ro["aa"],
        daa_sim=rf["aa"] - rb["aa"],
        denc=min(rb["denc"], rf["denc"]),
        aa_back=rb["aa"], aa_fwd=rf["aa"], aa_osc=oo["aa"],
        t_back=rb["t_years"], t_fwd=rf["t_years"])


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
        failed.append(rec[1]); logger.warning(f"{rec[1]}: {rec[2]}")
    else:
        rows.append(rec)

logger.info(f"integrated {len(rows)} comets x2 legs; {len(failed)} failed")

# ------------------------------------------------------------------
# Statistics
# ------------------------------------------------------------------

def resid_cap(v, th, X_extra):
    """Model-generated log rotation regressed on supplied covariates.

    This is a residual of the standard-dynamics prediction, not an
    observed-minus-model residual.
    """
    inc = th < 60
    # Report residual amplitudes in dex (the unit used by downstream tables).
    y = np.log10(v)
    X = np.column_stack([np.ones(len(y))] + X_extra)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    u = mannwhitneyu(resid[inc], resid[~inc], alternative="greater")
    rho, p = spearmanr(th, resid)
    return {"p": float(u.pvalue), "rho": float(rho), "p_2sided": float(p),
            "med_resid_in": float(np.median(resid[inc])),
            "med_resid_out": float(np.median(resid[~inc]))}

def run(sub):
    th    = np.array([r["theta"]    for r in sub])
    vc    = np.array([r["drot_cat"] for r in sub])
    vs    = np.array([r["drot_sim"] for r in sub])
    leg   = np.array([r["d_back"]   for r in sub] + [r["d_fwd"] for r in sub])
    daa_c = np.array([r["daa_cat"]  for r in sub])
    daa_s = np.array([r["daa_sim"]  for r in sub])
    K     = np.abs(np.array([r["daa_sim"] for r in sub]))
    Q     = np.array([r["q"]    for r in sub])
    I     = np.array([r["i"]    for r in sub])
    D     = np.array([r["denc"] for r in sub])
    inc   = th < 60
    out = {"n": len(sub), "n_in": int(inc.sum()),
           # --- validation of the independent integration vs the catalogue ---
           "leg_offset_med":  float(np.median(leg)),
           "leg_offset_p95":  float(np.percentile(leg, 95)),
           "leg_offset_max":  float(leg.max()),
           "frac_legs_within_0p01": float(np.mean(leg < 0.01)),
           "frac_legs_within_0p1":  float(np.mean(leg < 0.1)),
           "med_drot_cat_in":  float(np.median(vc[inc])),
           "med_drot_cat_out": float(np.median(vc[~inc])),
           "med_drot_sim_in":  float(np.median(vs[inc])),
           "med_drot_sim_out": float(np.median(vs[~inc])),
           "med_abs_drot_diff": float(np.median(np.abs(vs - vc))),
           "frac_comets_drot_within_0p02": float(np.mean(np.abs(vs - vc) < 0.02))}
    rho, p = spearmanr(vs, vc)
    out["drot_sim_vs_cat"] = {"rho": float(rho), "p_2sided": float(p)}
    rho, p = spearmanr(daa_s, daa_c)
    out["energy_validation"] = {"rho": float(rho), "p_2sided": float(p),
        "med_abs_diff": float(np.median(np.abs(daa_s - daa_c)))}
    # --- does the modelled rotation itself carry a cap preference? ---
    u = mannwhitneyu(vs[inc], vs[~inc], alternative="two-sided")
    out["drot_sim_in_vs_out_p"] = float(u.pvalue)
    rho, p = spearmanr(th, vs)
    out["drot_sim_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
    # --- the anomaly instrument: rotation beyond what the energy kick and
    #     encounter geometry predict (rotation-per-kick residual) ---
    out["rot_per_kick_resid_kick"] = resid_cap(vs, th, [np.log10(K + 1.0)])
    out["rot_per_kick_resid_kick_denc"] = resid_cap(vs, th, [np.log10(K + 1.0), np.log10(D)])
    out["rot_per_kick_resid_full"] = resid_cap(vs, th,
        [np.log10(K + 1.0), np.log10(D), Q, I])
    return out

res = {
    "method": "REBOUND/IAS15 bidirectional: osculating elements at CODE "
              "perihelion epoch integrated to the +/-250 AU barycentric "
              "sphere under Sun + planetary-system barycentres (DE440s). "
              "Boundary solutions evaluated relative to the system "
              "barycentre, matching the catalogue's original/future "
              "convention.",
    "finding": "drot_sim reproduces drot_cat per comet -- the catalogue's "
               "orig->fut rotation is reproduced by standard planetary "
               "dynamics. The cap statistic below is therefore a population "
               "contrast in model-predicted rotation per encounter budget; it "
               "is not, by itself, an independent observed-minus-dynamics "
               "residual or a TEP detection.",
    "n_all_c1": len(rows), "n_failed": len(failed), "failed": failed,
    "boundary_AU_barycentric": R_STOP,
    "ephemeris": "data/raw/spice/de440s.bsp (provenance pinned)",
    "matched": run([r for r in rows if r["q"] < 3.1]),
    "all_c1":  run(rows),
}

out = str(RESULTS / "step_b28_bidirectional_rotation.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b28_bidirectional_rotation.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

# ------------------------------------------------------------------
# Figure: independent integration vs catalogue boundary solutions
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

th  = np.array([r["theta"]    for r in rows])
vc  = np.array([r["drot_cat"] for r in rows])
vs  = np.array([r["drot_sim"] for r in rows])
leg = np.array([r["d_back"]   for r in rows] + [r["d_fwd"] for r in rows])
inc = th < 60

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axes[0]
ax.scatter(vs[~inc], vc[~inc], s=22, c="0.55", alpha=0.75,
           label=f"outside 60° cap ($n={int((~inc).sum())}$)")
ax.scatter(vs[inc], vc[inc], s=26, c="crimson", alpha=0.85,
           label=f"inside 60° cap ($n={int(inc.sum())}$)")
lim = max(vc.max(), vs.max()) * 1.15
lo = min(vc.min(), vs.min()) * 0.8
ax.plot([lo, lim], [lo, lim], "k:", lw=1, alpha=0.5)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lo, lim); ax.set_ylim(lo, lim)
ax.set_xlabel(r"modelled rotation $d_{\rm rot,sim}$ (deg)")
ax.set_ylabel(r"catalogued discrepancy $d_{\rm of}$ (deg)")
ax.legend(frameon=False, fontsize=9, loc="upper left")

ax = axes[1]
ax.hist(np.log10(leg), bins=30, color="steelblue", alpha=0.8)
ax.axvline(np.log10(0.01), color="k", ls=":", lw=1)
ax.set_xlabel(r"per-leg boundary offset $\log_{10}$(deg)")
ax.set_ylabel("legs")

fig.tight_layout()
FIG = RESULTS / "figures" / "supplementary"; FIG.mkdir(parents=True, exist_ok=True)
fig.savefig(FIG / "step_b28_rotation_residual.png", dpi=300)

for tag, d in [("matched", res["matched"]), ("all_c1", res["all_c1"])]:
    logger.info(f"{tag}: leg offset med/p95 {d['leg_offset_med']:.4f}/{d['leg_offset_p95']:.3f} deg | "
                f"sim-vs-cat rho={d['drot_sim_vs_cat']['rho']:.4f} | "
                f"med|drot diff|={d['med_abs_drot_diff']:.4f} | "
                f"rot-per-kick resid (kick+denc+q+i) cap p={d['rot_per_kick_resid_full']['p']:.4f} "
                f"rho={d['rot_per_kick_resid_full']['rho']:.3f}")
logger.data_save(out)
logger.data_save(csv_out)
logger.data_save(FIG / "step_b28_rotation_residual.png")