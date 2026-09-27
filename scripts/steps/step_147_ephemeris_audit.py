"""Step 147: ephemeris-version audit (step_b111).

The bidirectional boundary integration of step_063 runs under DE440s;
the CODE osculating solutions were fitted under an earlier ephemeris
generation (DE405/DE430 era).  A sceptic can ask whether the in-cap
rotation gap is a planetary-ephemeris version artefact -- the catalogue
fit and the independent propagation disagreeing on the planets, not on
the comet.

T1  Per-leg ephemeris displacement: each comet's osculating elements at
    perihelion are integrated to the +/-250 AU barycentric sphere under
    Sun + planetary-system barycentres, once with DE440s and once with
    DE430, and the boundary periapsis directions are compared leg by
    leg (median / p95 / max separation in degrees).

T2  Gap stability: the kick-and-geometry-regressed in-cap residual gap
    (the same four-covariate construction as the declared-anomaly
    instrument) is recomputed under each ephemeris; the gap and its
    MWU p-value are reported for both.

Outputs: results/step_b111_ephemeris_audit.json/.csv and
results/figures/supplementary/step_b111_ephemeris_audit.png.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))

import json
import math
import re
import numpy as np
from html.parser import HTMLParser
from scipy.stats import mannwhitneyu, spearmanr
import rebound
import spiceypy as sp

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (DATA_RAW, RESULTS, tee_stdout,
                                       lv, sep)

logger = StepLogger("step_147_ephemeris_audit")
tee_stdout(logger)

SEED = 20260919
KERNELS = ["de440s", "de430"]
CAP = 60.0

# ------------------------------------------------------------------
# CODE table parser (identical to step_063)
# ------------------------------------------------------------------


class TP(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.cur = []
        self.buf = ""
        self.in_td = False

    def handle_starttag(self, t, a):
        if t == "tr":
            self.cur = []
        elif t == "td":
            self.in_td = True
            self.buf = ""

    def handle_endtag(self, t):
        if t == "td":
            self.in_td = False
            self.cur.append(self.buf.strip())
        elif t == "tr" and self.cur:
            self.rows.append(self.cur)

    def handle_data(self, d):
        if self.in_td:
            self.buf += d


def parse_code(path):
    p = TP()
    p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14:
            continue
        try:
            out[r[0].strip()] = dict(
                desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                T=r[7].strip(), q=float(r[8]), e=float(r[9]),
                w=float(r[10]), Om=float(r[11]), i=float(r[12]),
                aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out


def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = (np.cos(om), np.sin(om), np.cos(Om),
                             np.sin(Om), np.cos(inc), np.sin(inc))
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci, so * si])


TNO = lv(34, -13)

# ------------------------------------------------------------------
# Ephemeris setup (identical to step_063)
# ------------------------------------------------------------------

AU_KM = 149597870.7
DAY_YR = 365.25
EPS = math.radians(23.4392911)
RX = np.array([[1, 0, 0],
               [0, math.cos(EPS), math.sin(EPS)],
               [0, -math.sin(EPS), math.cos(EPS)]])

GM_SUN = 1.32712440018e11
GM = {"1": 2.2031868551e4, "2": 3.2485859200e5, "3": 4.0350323562e5,
      "4": 4.2828375814e4, "5": 1.2671276480e8, "6": 3.7940626000e7,
      "7": 5.7945490100e6, "8": 6.8365271006e6, "9": 1.0868657e3}
PLANET_IDS = list(GM.keys())
MU = 4 * math.pi ** 2
R_STOP = 250.0
T_MAX = 20000.0
DT_OUT = 1.0


def jd_tt(y, m, d):
    a = (14 - m) // 12
    y2 = y + 4800 - a
    m2 = m + 12 * a - 3
    jdn = d + (153 * m2 + 2) // 5 + 365 * y2 + y2 // 4 - y2 // 100 \
        + y2 // 400 - 32045
    return jdn - 0.5 + (int(d) - d)


def perihelion_et(T):
    m = re.match(r"(\d{4})\s+(\d{1,2})\s+(\d+\.?\d*)", T)
    if not m:
        return None
    return (jd_tt(int(m.group(1)), int(m.group(2)), float(m.group(3)))
            - 2451545.0) * 86400.0


def body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return (RX @ np.array(st[:3]) / AU_KM,
            RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR)


def init_sim(et):
    sim = rebound.Simulation()
    sim.G = MU
    ps, vs = body_state("10", et)
    sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2],
            m=1.0)
    for b in PLANET_IDS:
        pp, vv = body_state(b, et)
        sim.add(x=pp[0], y=pp[1], z=pp[2],
                vx=vv[0], vy=vv[1], vz=vv[2], m=GM[b] / GM_SUN)
    return sim, ps, vs


def state_at_periapsis(ro):
    w_, O_, i_ = map(math.radians, (ro["w"], ro["Om"], ro["i"]))
    ph = perih_dir(w_, O_, i_)
    hh = np.array([math.sin(i_) * math.sin(O_),
                   -math.sin(i_) * math.cos(O_), math.cos(i_)])
    th = np.cross(hh, ph)
    th /= np.linalg.norm(th)
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
    sim.exit_min_distance = 0.001  # collision scale: bound IAS15 against step collapse
    denc = np.full(sim.N - 1, np.inf)
    t = direction * DT_OUT
    while abs(t) < T_MAX:
        sim.integrate(t, exact_finish_time=0)
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y,
                          p[nc].z - p[0].z])
        for j in range(1, sim.N - 1):
            d = math.sqrt((p[nc].x - p[j].x) ** 2
                          + (p[nc].y - p[j].y) ** 2
                          + (p[nc].z - p[j].z) ** 2)
            if d < denc[j - 1]:
                denc[j - 1] = d
        if np.linalg.norm(r_rel) >= R_STOP:
            break
        t += direction * DT_OUT
    else:
        return None
    mtot = sum(pp.m for pp in sim.particles)
    rb = np.zeros(3)
    vb = np.zeros(3)
    for pp in sim.particles:
        rb += pp.m * np.array([pp.x, pp.y, pp.z])
        vb += pp.m * np.array([pp.vx, pp.vy, pp.vz])
    rb /= mtot
    vb /= mtot
    p = sim.particles
    r_rel = np.array([p[nc].x, p[nc].y, p[nc].z]) - rb
    v_rel = np.array([p[nc].vx, p[nc].vy, p[nc].vz]) - vb
    phat, aa = boundary_orbit(r_rel, v_rel, mtot)
    return dict(phat=phat, aa=aa, denc=float(denc.min()),
                t_years=float(t))


# ------------------------------------------------------------------
# Sample: identical selection to step_063
# ------------------------------------------------------------------

osc = parse_code(str(DATA_RAW / "code" / "code_osculating.html"))
orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))
fut = parse_code(str(DATA_RAW / "code" / "code_future.html"))
warsaw = {l[5:17].strip()
          for l in open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
          if len(l) > 115}

sample = []
for k, ro in orig.items():
    if k not in fut:
        continue
    if not (0 < ro["aa"] < 100
            and ro["cls"] in ("1a", "1a+", "1b")) or k in warsaw:
        continue
    if k not in osc:
        continue
    et = perihelion_et(osc[k]["T"])
    if et is None:
        continue
    sample.append((k, ro, osc[k], et))

logger.info(f"sample: {len(sample)} class-1 CODE comets")

sp.furnsh(str(DATA_RAW / "naif" / "naif0012.tls"))

# ------------------------------------------------------------------
# Integrate both legs under each kernel
# ------------------------------------------------------------------

def _run_comet(k, ro, oo, et):
    rb = integrate_leg(oo, et, -1)
    rf = integrate_leg(oo, et, +1)
    if rb is None or rf is None:
        return None
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                   math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]),
                   math.radians(fut[k]["Om"]),
                   math.radians(fut[k]["i"]))
    return dict(
        phat_back=rb["phat"], phat_fwd=rf["phat"],
        drot_cat=sep(-po, -pf),
        drot_sim=sep(rb["phat"], rf["phat"]),
        daa_sim=rf["aa"] - rb["aa"],
        denc=min(rb["denc"], rf["denc"]),
        theta=sep(-po, TNO), q=ro["q"], i=ro["i"])


def _work_comet(job):
    k, ro, oo, et = job
    try:
        return _run_comet(k, ro, oo, et)
    except Exception as e:
        return ("__error__", k, str(e))


def _init_kern(kern):
    # forked children inherit the parent's BSP fds; concurrent spkezr
    # reads through a shared descriptor corrupt each other -- reopen the
    # kernels so each worker holds its own file handles.
    sp.kclear()
    sp.furnsh(str(DATA_RAW / "naif" / "naif0012.tls"))
    sp.furnsh(str(DATA_RAW / "spice" / f"{kern}.bsp"))


from scripts.utils.parallel import default_workers as _default_workers

if "--workers" in _sys.argv:
    N_WORK = max(1, int(_sys.argv[_sys.argv.index("--workers") + 1]))
else:
    N_WORK = _default_workers()

import multiprocessing as mp
_ctx = mp.get_context("fork") if _sys.platform != "win32" \
    else mp.get_context("spawn")

runs = {}
for kern in KERNELS:
    if N_WORK > 1 and len(sample) > 1:
        with _ctx.Pool(min(N_WORK, len(sample)),
                       initializer=_init_kern, initargs=(kern,)) as pool:
            recs = pool.map(_work_comet, sample)
    else:
        _init_kern(kern)
        recs = [_work_comet(s) for s in sample]
    rows = {}
    failed = []
    for (k, ro, oo, et), rec in zip(sample, recs):
        if rec is None:
            failed.append(k)
            logger.warning(f"{k}: boundary not reached")
        elif isinstance(rec, tuple) and rec[0] == "__error__":
            failed.append(rec[1])
            logger.warning(f"{rec[1]}: {rec[2]}")
        else:
            rows[k] = rec
    runs[kern] = dict(rows=rows, failed=failed)
    logger.info(f"{kern}: {len(rows)} comets x2 legs; "
                f"{len(failed)} failed")

common = sorted(set(runs["de440s"]["rows"]) & set(runs["de430"]["rows"]))
logger.info(f"{len(common)} comets integrated under both kernels")

# ------------------------------------------------------------------
# T1: per-leg boundary-direction displacement between ephemerides
# ------------------------------------------------------------------

off_back = np.array([sep(runs["de440s"]["rows"][k]["phat_back"],
                         runs["de430"]["rows"][k]["phat_back"])
                     for k in common])
off_fwd = np.array([sep(runs["de440s"]["rows"][k]["phat_fwd"],
                        runs["de430"]["rows"][k]["phat_fwd"])
                    for k in common])
off_all = np.concatenate([off_back, off_fwd])

t1 = {
    "n_legs": int(len(off_all)),
    "median_deg_back": float(np.median(off_back)),
    "median_deg_fwd": float(np.median(off_fwd)),
    "median_deg_all": float(np.median(off_all)),
    "p95_deg_all": float(np.percentile(off_all, 95)),
    "max_deg_all": float(off_all.max()),
    "frac_legs_within_1e-4_deg": float(np.mean(off_all < 1e-4)),
    "frac_legs_within_1e-3_deg": float(np.mean(off_all < 1e-3))}
logger.info(f"T1 ephemeris leg displacement: median "
            f"{t1['median_deg_all']:.3g} deg, max "
            f"{t1['max_deg_all']:.3g} deg")

# ------------------------------------------------------------------
# T2: in-cap residual gap under each ephemeris
# ------------------------------------------------------------------

def resid_gap(rows, keys):
    th = np.array([rows[k]["theta"] for k in keys])
    vs = np.array([rows[k]["drot_sim"] for k in keys])
    K = np.abs(np.array([rows[k]["daa_sim"] for k in keys]))
    D = np.array([rows[k]["denc"] for k in keys])
    Q = np.array([rows[k]["q"] for k in keys])
    I = np.array([rows[k]["i"] for k in keys])
    inc = th < CAP
    y = np.log10(np.clip(vs, 1e-9, None))
    X = np.column_stack([np.ones(len(y)), np.log10(K + 1.0),
                         np.log10(np.clip(D, 1e-3, None)), Q, I])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    u = mannwhitneyu(resid[inc], resid[~inc], alternative="greater")
    return {"n": len(keys), "n_in": int(inc.sum()),
            "gap_dex": float(np.median(resid[inc])
                             - np.median(resid[~inc])),
            "med_resid_in": float(np.median(resid[inc])),
            "med_resid_out": float(np.median(resid[~inc])),
            "p": float(u.pvalue)}

t2 = {kern: resid_gap(runs[kern]["rows"], common) for kern in KERNELS}
for kern in KERNELS:
    logger.info(f"T2 {kern}: in-cap residual gap "
                f"{t2[kern]['gap_dex']:+.4f} dex, p={t2[kern]['p']:.4g}")

# ------------------------------------------------------------------
# verdict
# ------------------------------------------------------------------

leg_ok = t1["p95_deg_all"] < 1e-3
gap_diff = abs(t2["de440s"]["gap_dex"] - t2["de430"]["gap_dex"])
both_sig = (t2["de440s"]["p"] < 0.05 and t2["de430"]["p"] < 0.05)
res = {"step": "step_147_ephemeris_audit", "result": "b111",
       "seed": SEED, "kernels": KERNELS,
       "n_common": len(common),
       "n_failed": {k: len(runs[k]["failed"]) for k in KERNELS},
       "T1_ephemeris_displacement": t1,
       "T2_incap_gap_by_kernel": t2,
       "test_summary": {
           "legs_ephemeris_stable": bool(leg_ok),
           "gap_ephemeris_stable": bool(gap_diff < 0.01),
           "gap_significant_both": bool(both_sig)}}

if leg_ok and both_sig:
    res["verdict"] = (
        "EPHEMERIS-INVARIANT: the DE440s/DE430 boundary periapsis "
        f"directions agree to a median of {t1['median_deg_all']:.3g} deg "
        f"per leg (p95 {t1['p95_deg_all']:.3g}, max "
        f"{t1['max_deg_all']:.3g}), and the kick-and-geometry-regressed "
        f"in-cap gap is {t2['de440s']['gap_dex']:+.4f} dex "
        f"(p={t2['de440s']['p']:.3g}) under DE440s versus "
        f"{t2['de430']['gap_dex']:+.4f} dex (p={t2['de430']['p']:.3g}) "
        "under DE430 -- the ephemeris generation does not perturb the "
        "cometary boundary reconstruction or the in-cap excess.")
elif leg_ok:
    res["verdict"] = (
        "EPHEMERIS-STABLE LEGS, GAP MARGINAL: boundary periapsis "
        f"directions agree to a median of {t1['median_deg_all']:.3g} deg "
        f"per leg, but the in-cap gap is {t2['de440s']['gap_dex']:+.4f} "
        f"versus {t2['de430']['gap_dex']:+.4f} dex (p "
        f"{t2['de440s']['p']:.3g} vs {t2['de430']['p']:.3g}) -- inspect "
        "T2 before reading the gap as ephemeris-independent.")
else:
    res["verdict"] = (
        "EPHEMERIS SENSITIVE: per-leg boundary displacement "
        f"p95 {t1['p95_deg_all']:.3g} deg exceeds the 1e-3 bound -- the "
        "reconstruction carries a kernel dependence that must be "
        "accounted before interpreting the in-cap gap.")
res["evidence_status"] = "reduction-control audit"
res["inputs"] = [
    "data/raw/code/code_{osculating,original,future}.html",
    "data/raw/warsaw/warsaw_tablec.dat (exclusion list)",
    "data/raw/spice/de440s.bsp, de430.bsp, naif0012.tls"]
res["caveats"] = [
    "Tests the two planetary-ephemeris generations against each "
    "other on the catalogue's own osculating elements; a common-mode "
    "ephemeris error shared by both kernels is not probed.",
    "The in-cap gap is the standard-dynamics residual instrument; "
    "its ephemeris stability bounds reduction-level contamination, "
    "not the anomaly's interpretation."]

out = RESULTS / "step_b111_ephemeris_audit.json"
with open(out, "w") as _fh:
    json.dump(res, _fh, indent=1)
logger.data_save(out)

import csv
with open(RESULTS / "step_b111_ephemeris_audit.csv", "w",
          newline="") as f:
    w = csv.writer(f)
    w.writerow(["desig", "off_back_deg", "off_fwd_deg"])
    for k in common:
        w.writerow([k,
                    f"{sep(runs['de440s']['rows'][k]['phat_back'], runs['de430']['rows'][k]['phat_back']):.6e}",
                    f"{sep(runs['de440s']['rows'][k]['phat_fwd'], runs['de430']['rows'][k]['phat_fwd']):.6e}"])

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
ax[0].hist(off_back, bins=25, alpha=0.7, label="inbound leg")
ax[0].hist(off_fwd, bins=25, alpha=0.7, label="outbound leg")
ax[0].set_xscale("log")
ax[0].set_xlabel("DE440s vs DE430 boundary periapsis displacement (deg)")
ax[0].set_ylabel("legs")
ax[0].legend()
ax[0].set_title("T1: per-leg ephemeris displacement")
g440 = [runs["de440s"]["rows"][k]["drot_sim"] for k in common]
g430 = [runs["de430"]["rows"][k]["drot_sim"] for k in common]
ax[1].scatter(g440, g430, s=12, alpha=0.6)
lim = [min(g440 + g430) * 0.9, max(g440 + g430) * 1.1]
ax[1].plot(lim, lim, "k--", lw=1)
ax[1].set_xlabel("drot under DE440s (deg)")
ax[1].set_ylabel("drot under DE430 (deg)")
ax[1].set_title("T2: per-comet boundary rotation, kernel vs kernel")
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.tight_layout()
fig.savefig(FIG / "supplementary" / "step_b111_ephemeris_audit.png", dpi=300)
logger.data_save(FIG / "supplementary" / "step_b111_ephemeris_audit.png")

print("TEST SUMMARY:\n" + json.dumps(res["test_summary"], indent=1))
print(f"VERDICT: {res['verdict']}")
