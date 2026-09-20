"""step_066: Secular stability of the real detached-TNO cluster under
the measured planetary system, with and without the published
Planet Nine.

The shepherding hypothesis requires the observed varpi alignment to be
*maintained*: differential apsidal precession under the known planets
otherwise disperses the cluster on a measurable timescale.  This step
integrates the actual catalogued population -- the 44 detached TNOs
(a > 150 AU, q > 30 AU, condition code <= 3; identical selection to
step_010) -- forward through the real planetary system (Sun + giant
planets on DE440s states) and measures the apsidal-precession rates
directly.

Three realizations are compared:

  giants        : Sun + Jupiter, Saturn, Uranus, Neptune (DE440s)
  giants+bb21ml : + Brown & Batygin (2021) maximum-likelihood Planet Nine
                  (m9=5.0 M_E, a9=300 AU, e9=0.15, i9=17 deg,
                  varpi9=254 deg, Om9=108 deg, M9=0)
  giants+bb21med: + BB21 marginalized-median model
                  (m9=6.9 M_E, a9=461 AU, e9=0.30, i9=15.6 deg,
                  varpi9=246.7 deg, Om9=96.9 deg, M9=0)

The TNOs are massless test particles started from their catalogued
osculating elements at a common epoch (median element epoch; along-track
phase is drawn uniformly since SBDB supplies no mean anomaly -- apsidal
precession is phase-independent).

Measured per object over the integration:
  dvarpi/dt   : linear drift of the longitude of perihelion -> the
                differential-precession dispersal rate
  Dvarpi(t)   : varpi relative to P9's varpi -- the libration coordinate
                a confining perturber must hold near 180 deg

Reported:
  sigma(dvarpi/dt) per model -> dispersal time t_disp = sigma_varpi /
  sigma(dvarpi/dt) for the observed cluster spread;
  the P9 restoring signature: correlation of the P9-added drift
  dvarpi_P9 - dvarpi_giants against the offset (Dvarpi - 180 deg);
  R_varpi(t) of the observed cluster under each model.

Outputs
-------
results/step_b31_tno_secular.json
results/step_b31_tno_secular.csv
results/figures/supplementary/step_b31_tno_secular.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_066_tno_secular_stability")
tee_stdout(logger)
logger.header("Secular stability of the detached-TNO cluster (REBOUND/DE440s)")

import csv
import json
import math
import numpy as np
from scipy.stats import spearmanr

import rebound
import spiceypy as sp

# ------------------------------------------------------------------
# Ephemeris (identical convention to steps 042/062/063)
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
GM = {"5": 1.2671276480e8, "6": 3.7940626000e7,
      "7": 5.7945490100e6, "8": 6.8365271006e6}   # J,S,U,N only
GIANT_IDS = list(GM.keys())

def body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return RX @ np.array(st[:3]) / AU_KM, RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR

MU = 4 * math.pi ** 2
M_EARTH_SUN = 3.0035e-6

P9_MODELS = {
    "bb21_ml":  dict(m9=5.0, a9=300.0, e9=0.15, i9=17.0,
                     varpi9=254.0, Om9=108.0),
    "bb21_med": dict(m9=6.9, a9=461.0, e9=0.30, i9=15.6,
                     varpi9=246.7, Om9=96.9),
}

# ------------------------------------------------------------------
# Detached-TNO sample (identical cuts to step_010 primary sample)
# ------------------------------------------------------------------

def load_sbdb(name="sbdb_outer_ss.json"):
    d = json.load(open(DATA_RAW / "sbdb" / name))
    fields = d["fields"]
    return [dict(zip(fields, rec)) for rec in d["data"]]

def fnum(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return float("nan")

rows = load_sbdb()
tnos = []
for r in rows:
    a, q, cc = fnum(r, "a"), fnum(r, "q"), fnum(r, "condition_code")
    if not (a > 150 and q > 30 and 0 <= cc <= 3):
        continue
    tnos.append(dict(name=r["full_name"], a=a, e=fnum(r, "e"),
                     i=fnum(r, "i"), Om=fnum(r, "om"), w=fnum(r, "w"),
                     q=q, epoch=fnum(r, "epoch"),
                     varpi=(fnum(r, "om") + fnum(r, "w")) % 360.0))

logger.info(f"detached sample: {len(tnos)} TNOs (a>150, q>30, cc<=3)")

# common epoch: median element epoch (JD TT) -> ET
T0_JD = float(np.median([t["epoch"] for t in tnos]))
ET0 = (T0_JD - 2451545.0) * 86400.0
logger.info(f"common epoch JD {T0_JD:.1f}")

SEED = 20260918
rng = np.random.default_rng(SEED)

# ------------------------------------------------------------------
# Integration
# ------------------------------------------------------------------

T_END   = 2.0e7          # 20 Myr -- enough for clean dvarpi/dt slopes
DT      = 0.5            # yr (resolves Jupiter at P/24)
SNAP    = 5.0e4          # snapshot cadence (yr)

def build_sim(p9=None, m_arr=None):
    """Sun + giants (+P9) at T0, TNO test particles from elements."""
    sim = rebound.Simulation()
    sim.G = MU
    ps, vs = body_state("10", ET0)
    sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2], m=1.0)
    for b in GIANT_IDS:
        pp, vv = body_state(b, ET0)
        sim.add(x=pp[0], y=pp[1], z=pp[2], vx=vv[0], vy=vv[1], vz=vv[2],
                m=GM[b] / GM_SUN)
    if p9 is not None:
        w9 = p9["varpi9"] - p9["Om9"]
        sim.add(a=p9["a9"], e=p9["e9"], inc=math.radians(p9["i9"]),
                Omega=math.radians(p9["Om9"]), omega=math.radians(w9),
                M=0.0, m=p9["m9"] * M_EARTH_SUN)
        # P9 elements are heliocentric; shift into the barycentric frame
        i9 = sim.N - 1
        sim.particles[i9].x  += ps[0]; sim.particles[i9].y  += ps[1]
        sim.particles[i9].z  += ps[2]
        sim.particles[i9].vx += vs[0]; sim.particles[i9].vy += vs[1]
        sim.particles[i9].vz += vs[2]
    for j, t in enumerate(tnos):
        sim.add(a=t["a"], e=t["e"], inc=math.radians(t["i"]),
                Omega=math.radians(t["Om"]), omega=math.radians(t["w"]),
                M=float(m_arr[j]))
        i = sim.N - 1
        sim.particles[i].x  += ps[0]; sim.particles[i].y  += ps[1]
        sim.particles[i].z  += ps[2]
        sim.particles[i].vx += vs[0]; sim.particles[i].vy += vs[1]
        sim.particles[i].vz += vs[2]
    sim.integrator = "whfast"
    sim.dt = DT
    return sim

def run_model(label, p9, m_arr=None):
    sim = build_sim(p9, m_arr=m_arr)
    n_tno = len(tnos)
    i0 = sim.N - n_tno                     # first TNO particle index
    i_p9 = i0 - 1 if p9 is not None else None
    n_snap = int(T_END / SNAP) + 1
    varpi_t = np.zeros((n_snap, n_tno))
    p9_varpi_t = np.zeros(n_snap)
    ts = np.zeros(n_snap)
    sun = sim.particles[0]
    for k in range(n_snap):
        t = k * SNAP
        sim.integrate(t, exact_finish_time=1)
        ts[k] = t
        for j in range(n_tno):
            orb = sim.particles[i0 + j].orbit(primary=sun)
            varpi_t[k, j] = orb.Omega + orb.omega
        if i_p9 is not None:
            orb9 = sim.particles[i_p9].orbit(primary=sun)
            p9_varpi_t[k] = orb9.Omega + orb9.omega
    logger.info(f"{label}: integrated {T_END/1e6:.0f} Myr, {n_snap} snapshots")
    return ts, varpi_t, p9_varpi_t

# Pre-draw every M variate in the parent's rng, in the order the serial
# run_model calls consumed them, so parallel runs are bit-identical.
_n = len(tnos)
m_g   = np.array([rng.uniform(0, 2 * math.pi) for _ in range(_n)])
m_ml  = np.array([rng.uniform(0, 2 * math.pi) for _ in range(_n)])
m_med = np.array([rng.uniform(0, 2 * math.pi) for _ in range(_n)])

jobs = [("giants",         None,                 m_g),
        ("giants+bb21_ml", P9_MODELS["bb21_ml"],  m_ml),
        ("giants+bb21_med",P9_MODELS["bb21_med"], m_med)]

def _work_model(job):
    label, p9, ma = job
    return label, run_model(label, p9, m_arr=ma)

if len(jobs) > 1:
    import multiprocessing as mp
    ctx = mp.get_context("fork") if _sys.platform != "win32" \
        else mp.get_context("spawn")

    def _init():
        # forked children inherit the parent's BSP fd; concurrent spkezr
        # reads through a shared descriptor corrupt each other -- reopen
        # the kernel so each worker holds its own file handle.
        sp.kclear()
        sp.furnsh(str(SPK))

    with ctx.Pool(len(jobs), initializer=_init) as pool:
        results = dict(pool.map(_work_model, jobs))
else:
    results = dict([_work_model(j) for j in jobs])

ts, vw_g, _        = results["giants"]
_,  vw_ml, p9_ml   = results["giants+bb21_ml"]
_,  vw_med, p9_med = results["giants+bb21_med"]

# ------------------------------------------------------------------
# Analysis: per-object apsidal drift rates
# ------------------------------------------------------------------

def unwrap_deg(x):
    return np.degrees(np.unwrap(np.radians(x)))

def drift(vw):
    """Linear dvarpi/dt per object (deg/Myr) from snapshot series."""
    n_snap, n = vw.shape
    out = np.zeros(n)
    tt = ts / 1e6
    for j in range(n):
        y = unwrap_deg(vw[:, j])
        out[j] = np.polyfit(tt, y, 1)[0]
    return out

def cluster_r(vw_row):
    z = np.exp(1j * vw_row)
    return float(np.abs(z.mean()))

dv_g   = drift(vw_g)
dv_ml  = drift(vw_ml)
dv_med = drift(vw_med)

# in-cap membership: the pre-declared 60 deg cap about (49,-17)
def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO_AXIS = lv(49.0, -17.0)
theta = np.array([sep(perih_dir(math.radians(t["w"]), math.radians(t["Om"]),
                                math.radians(t["i"])), TNO_AXIS)
                  for t in tnos])
incap = theta < 60

# observed cluster spread today (circular std of varpi, in-cap members)
def circ_std(deg):
    z = np.exp(1j * np.radians(deg))
    R = np.abs(z.mean())
    return math.degrees(math.sqrt(-2 * math.log(max(R, 1e-12))))

sig_varpi_obs = circ_std(np.array([t["varpi"] for t in tnos])[incap])
sig_varpi_all = circ_std(np.array([t["varpi"] for t in tnos]))

def model_stats(dv, label):
    sd_all = float(np.std(dv))
    sd_in  = float(np.std(dv[incap]))
    return {
        "sigma_dv_all_deg_Myr": sd_all,
        "sigma_dv_incap_deg_Myr": sd_in,
        "t_disp_all_Myr": float(sig_varpi_all / sd_all) if sd_all > 0 else None,
        "t_disp_incap_Myr": float(sig_varpi_obs / sd_in) if sd_in > 0 else None,
        "med_dv_incap": float(np.median(dv[incap])),
        "med_dv_outcap": float(np.median(dv[~incap])),
    }

# P9 restoring signature: does the P9-added drift push each object
# toward Dvarpi = 180 deg?  delta_dv = dv_P9 - dv_giants; a confining
# perturber gives delta_dv anticorrelated with (Dvarpi0 - 180).
def restoring(dv_p9, p9_v0):
    dv0 = np.array([t["varpi"] for t in tnos])
    # Dvarpi = varpi_tno - varpi_p9 folded to [-180,180]; offset from the
    # anti-aligned libration point is |Dvarpi| - 180.  Restoring drift in
    # the toward-180 coordinate is sign(Dvarpi) * delta_dv.
    Dvarpi = ((dv0 - math.degrees(p9_v0)) + 540) % 360 - 180
    off = np.abs(Dvarpi) - 180.0
    delta = np.sign(Dvarpi) * (dv_p9 - dv_g)
    rho, p = spearmanr(off, delta)
    return {"rho": float(rho), "p_2sided": float(p),
            "med_delta_incap": float(np.median(delta[incap])),
            "med_delta_outcap": float(np.median(delta[~incap]))}

res = {
    "method": "44 catalogued detached TNOs (a>150,q>30,cc<=3) as test "
              "particles under Sun+giants (DE440s, WHFAST dt=0.5 yr) "
              "for 20 Myr; BB21 Planet Nine added in two realizations.  "
              "Per-object apsidal drift dvarpi/dt from unwrapped osculating "
              "varpi; dispersal time = circular-std(varpi)/sigma(dv).",
    "n_tno": len(tnos), "n_incap": int(incap.sum()),
    "cap_deg": 60.0, "axis": "49,-17",
    "T0_JD": T0_JD, "T_end_yr": T_END, "dt_yr": DT, "snap_yr": SNAP,
    "sig_varpi_incap_deg": sig_varpi_obs,
    "sig_varpi_all_deg": sig_varpi_all,
    "models": {
        "giants":        model_stats(dv_g, "giants"),
        "giants_bb21ml": model_stats(dv_ml, "giants+bb21_ml"),
        "giants_bb21med": model_stats(dv_med, "giants+bb21_med"),
    },
    "p9_restoring": {
        "bb21_ml":  restoring(dv_ml, p9_ml[0]),
        "bb21_med": restoring(dv_med, p9_med[0]),
    },
}

# R_varpi(t) of the observed cluster under each model
def R_series(vw):
    return [cluster_r(vw[k]) for k in range(vw.shape[0])]

res["R_varpi_t"] = {
    "ts_Myr": [float(x) for x in (ts / 1e6)],
    "giants": R_series(vw_g),
    "giants_bb21ml": R_series(vw_ml),
    "giants_bb21med": R_series(vw_med),
}

out = str(RESULTS / "step_b31_tno_secular.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b31_tno_secular.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "a", "q", "i", "varpi0", "theta", "in_cap",
                "dv_giants", "dv_bb21ml", "dv_bb21med"])
    for j, t in enumerate(tnos):
        w.writerow([t["name"], t["a"], t["q"], t["i"], t["varpi"],
                    theta[j], bool(incap[j]), dv_g[j], dv_ml[j], dv_med[j]])

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

ax = axes[0]
ax.plot(ts / 1e6, res["R_varpi_t"]["giants"], color="0.4",
        label="giants only")
ax.plot(ts / 1e6, res["R_varpi_t"]["giants_bb21ml"], color="crimson",
        label="+ BB21 ml")
ax.plot(ts / 1e6, res["R_varpi_t"]["giants_bb21med"], color="darkorange",
        label="+ BB21 med")
ax.set_xlabel("time (Myr)")
ax.set_ylabel(r"cluster resultant $R_\varpi$ (all 44)")
ax.set_ylim(0, 1)
ax.legend(frameon=False, fontsize=8)
ax.set_title("cluster coherence vs time", fontsize=10)

ax = axes[1]
for arr, lab, c in [(dv_g, "giants", "0.4"), (dv_ml, "+ml", "crimson"),
                    (dv_med, "+med", "darkorange")]:
    ax.scatter(theta[~incap], arr[~incap], s=14, c=c, alpha=0.4)
    ax.scatter(theta[incap], arr[incap], s=18, c=c, alpha=0.9, label=lab)
ax.axvline(60, color="k", ls=":", lw=1)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel(r"periapsis--axis separation $\theta$ (deg)")
ax.set_ylabel(r"$d\varpi/dt$ (deg/Myr)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("apsidal drift vs axis distance", fontsize=10)

ax = axes[2]
dv0 = np.array([t["varpi"] for t in tnos])
for dv_p9, p9_v0, lab, c in [(dv_ml, p9_ml[0], "ml", "crimson"),
                             (dv_med, p9_med[0], "med", "darkorange")]:
    Dvarpi = ((dv0 - math.degrees(p9_v0)) + 540) % 360 - 180
    off = np.abs(Dvarpi) - 180.0
    delta = np.sign(Dvarpi) * (dv_p9 - dv_g)
    ax.scatter(off[incap], delta[incap], s=18, c=c, alpha=0.9, label=f"{lab} in-cap")
    ax.scatter(off[~incap], delta[~incap], s=14, c=c, alpha=0.35)
ax.axhline(0, color="k", lw=0.5)
ax.axvline(0, color="k", ls=":", lw=1)
ax.set_xlabel(r"$|\Delta\varpi_0| - 180°$ (deg)")
ax.set_ylabel(r"$\delta(d\varpi/dt)$ toward 180° (deg/Myr)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("P9 restoring signature", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b31_tno_secular.png", dpi=300)

for m in ("giants", "giants_bb21ml", "giants_bb21med"):
    d = res["models"][m]
    logger.info(f"{m}: sigma(dv) all/in = {d['sigma_dv_all_deg_Myr']:.3f}/"
                f"{d['sigma_dv_incap_deg_Myr']:.3f} deg/Myr | "
                f"t_disp in = {d['t_disp_incap_Myr']:.0f} Myr")
for m, d in res["p9_restoring"].items():
    logger.info(f"restoring {m}: rho={d['rho']:.3f} p={d['p_2sided']:.4f} | "
                f"med delta in/out = {d['med_delta_incap']:.3f}/{d['med_delta_outcap']:.3f}")
logger.data_save(out)
logger.data_save(csv_out)
logger.data_save(FIG / "supplementary" / "step_b31_tno_secular.png")