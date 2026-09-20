"""step_088: Galactic-tide insertion on the resident population --
maintenance and regeneration tests.

step_070 measured the tide's contribution to the transit channel
(~1e-9 deg, null).  The resident side carries a different question:
the Galactic tide is the standard "no new physics" mechanism invoked
for outer-system apsidal structure, so the decisive test is whether
the tide can *maintain* the observed varpi cluster or *regenerate*
one from an isotropic start.  This step inserts the same epoch-fixed
inertial tide tensor used in step_070 (Kx=+OmG^2, Ky=-OmG^2,
Kz=-4 pi G rho) into the step_066 secular harness -- the 44 detached
TNOs (a>150, q>30, cc<=3) integrated 20 Myr under Sun+giants -- via
Strang-split velocity kicks.

Three realizations:

  tide_observed   : real catalogue varpi + giants + tide
                    -> does the tide hold or disperse R_varpi?
  tide_scrambled  : initial varpi randomized (a,e,i,Om fixed, w reset)
                    -> does the tide *create* varpi concentration, and
                       if so at what direction?
  tide_p9_med     : real varpi + giants + tide + BB21 median P9
                    -> does tide+perturber jointly confine?

Diagnostics per model: sigma(dvarpi/dt) and dispersal time (as 066);
R_varpi(t); the final mean varpi direction; and a restoring-torque
test against the tide's own ecliptic libration directions (the
galactic-plane points l=90/270 deg, the tide's stationary loci).

Outputs
-------
results/step_b53_tide_resident.json
results/step_b53_tide_resident.csv
results/figures/supplementary/step_b53_tide_resident.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (DATA_RAW, RESULTS, ECL2GAL,
                                       GAL2ECL, gv, lb, lv, tee_stdout)
logger = StepLogger("step_088_tide_resident")
tee_stdout(logger)
logger.header("Resident-side galactic-tide insertion (REBOUND/DE440s)")

import csv
import json
import math
import numpy as np
from scipy.stats import spearmanr

import rebound
import spiceypy as sp

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

# identical tide tensor to step_070 (epoch-fixed inertial approximation)
KMS_AUYR = 0.210805
AU_KPC   = 206264.806 * 1000.0
OMG      = 26.0 * KMS_AUYR / AU_KPC
G_AU     = 4 * math.pi ** 2
RHO_G    = 0.1 / (AU_KPC / 1000.0) ** 3
KZ       = -4.0 * math.pi * G_AU * RHO_G
KX, KY   = OMG ** 2, -OMG ** 2
TIDE     = np.diag([KX, KY, KZ])
logger.info(f"tide coefficients: Kx={KX:.3e} Ky={KY:.3e} Kz={KZ:.3e} yr^-2")

P9_MED = dict(m9=6.9, a9=461.0, e9=0.30, i9=15.6,
              varpi9=246.7, Om9=96.9)

# tide stationary directions: the galactic-plane points l=90,270 deg,
# expressed as ecliptic longitudes for context
TIDE_DIRS = [lb(gv(90.0, 0.0)), lb(gv(270.0, 0.0))]
logger.info("tide libration directions (ecliptic lon,lat): "
            + "; ".join(f"({l:.1f},{b:.1f})" for l, b in TIDE_DIRS))

# ------------------------------------------------------------------
# Detached-TNO sample (identical cuts to step_010 / step_066)
# ------------------------------------------------------------------

def load_sbdb(name="sbdb_outer_ss.json"):
    d = json.load(open(DATA_RAW / "sbdb" / name))
    return [dict(zip(d["fields"], rec)) for rec in d["data"]]

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

T0_JD = float(np.median([t["epoch"] for t in tnos]))
ET0 = (T0_JD - 2451545.0) * 86400.0

SEED = 20260918
rng = np.random.default_rng(SEED)

T_END = 2.0e7
DT    = 0.5
SNAP  = 5.0e4
KICK  = 1000.0          # Strang kick cadence (yr) << tide periods ~1e8 yr

def build_sim(scramble_w0=None, m_arr=None, p9=None):
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
        i9 = sim.N - 1
        sim.particles[i9].x += ps[0]; sim.particles[i9].y += ps[1]
        sim.particles[i9].z += ps[2]
        sim.particles[i9].vx += vs[0]; sim.particles[i9].vy += vs[1]
        sim.particles[i9].vz += vs[2]
    for j, t in enumerate(tnos):
        w0 = t["w"]
        if scramble_w0 is not None:
            w0 = (float(scramble_w0[j]) - t["Om"]) % 360.0
        sim.add(a=t["a"], e=t["e"], inc=math.radians(t["i"]),
                Omega=math.radians(t["Om"]), omega=math.radians(w0),
                M=float(m_arr[j]))
        i = sim.N - 1
        sim.particles[i].x += ps[0]; sim.particles[i].y += ps[1]
        sim.particles[i].z += ps[2]
        sim.particles[i].vx += vs[0]; sim.particles[i].vy += vs[1]
        sim.particles[i].vz += vs[2]
    sim.integrator = "whfast"
    sim.dt = DT
    return sim

def kick_all(sim, i0, n_tno, dtk):
    """tidal velocity kick on every TNO: v += a_tide(r_helio) * dtk."""
    p = sim.particles
    sun = p[0]
    for j in range(n_tno):
        pt = p[i0 + j]
        ag = ECL2GAL @ np.array([pt.x - sun.x, pt.y - sun.y, pt.z - sun.z])
        ae = GAL2ECL @ (TIDE @ ag)
        pt.vx += ae[0] * dtk
        pt.vy += ae[1] * dtk
        pt.vz += ae[2] * dtk

def run_model(label, scramble_w0=None, m_arr=None, p9=None):
    sim = build_sim(scramble_w0=scramble_w0, m_arr=m_arr, p9=p9)
    n_tno = len(tnos)
    i0 = sim.N - n_tno
    n_snap = int(T_END / SNAP) + 1
    varpi_t = np.zeros((n_snap, n_tno))
    ts = np.zeros(n_snap)
    sun = sim.particles[0]
    n_kick = int(SNAP / KICK)
    for k in range(n_snap):
        t = k * SNAP
        while sim.t < t:
            kick_all(sim, i0, n_tno, KICK / 2.0)
            sim.integrate(min(sim.t + KICK, t), exact_finish_time=1)
            kick_all(sim, i0, n_tno, KICK / 2.0)
        ts[k] = t
        for j in range(n_tno):
            orb = sim.particles[i0 + j].orbit(primary=sun)
            varpi_t[k, j] = orb.Omega + orb.omega
    logger.info(f"{label}: integrated {T_END/1e6:.0f} Myr "
                f"({n_kick} tide kicks per snapshot)")
    return ts, varpi_t

# Pre-draw every random variate in the parent's rng, in exactly the order
# the serial calls consumed them, so the parallel runs are bit-identical:
# observed: n_tno M-draws; scrambled: interleaved (w0, M) per TNO;
# p9 model: n_tno M-draws.
_n = len(tnos)
m_obs = np.array([rng.uniform(0, 2 * math.pi) for _ in range(_n)])
w_scr = np.empty(_n); m_scr = np.empty(_n)
for _j in range(_n):
    w_scr[_j] = rng.uniform(0.0, 360.0)
    m_scr[_j] = rng.uniform(0, 2 * math.pi)
m_p9 = np.array([rng.uniform(0, 2 * math.pi) for _ in range(_n)])

jobs = [("tide_observed",  None,  m_obs, None),
        ("tide_scrambled", w_scr, m_scr, None),
        ("tide_p9_med",    None,  m_p9,  P9_MED)]

def _work_model(job):
    label, sw, ma, p9 = job
    return label, run_model(label, scramble_w0=sw, m_arr=ma, p9=p9)

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

ts, vw_obs = results["tide_observed"]
_,  vw_scr = results["tide_scrambled"]
_,  vw_p9  = results["tide_p9_med"]

# ------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------

def unwrap_deg(x):
    return np.degrees(np.unwrap(np.radians(x)))

def drift(vw):
    out = np.zeros(vw.shape[1])
    tt = ts / 1e6
    for j in range(vw.shape[1]):
        out[j] = np.polyfit(tt, unwrap_deg(vw[:, j]), 1)[0]
    return out

def cluster_r(vw_row):
    return float(np.abs(np.exp(1j * vw_row).mean()))

def mean_dir(vw_row):
    return float(np.degrees(np.angle(np.exp(1j * vw_row).mean())) % 360.0)

def circ_std(deg):
    z = np.exp(1j * np.radians(deg))
    R = np.abs(z.mean())
    return math.degrees(math.sqrt(-2 * math.log(max(R, 1e-12))))

# cap membership on the real catalogue geometry (pre-declared resident axis)
TNO_AXIS = lv(49.0, -17.0)

def perih_dir_deg(om, Om, inc):
    om, Om, inc = map(math.radians, (om, Om, inc))
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

theta = np.array([sep(perih_dir_deg(t["w"], t["Om"], t["i"]), TNO_AXIS)
                  for t in tnos])
incap = theta < 60

varpi0 = np.array([t["varpi"] for t in tnos])
sig_varpi_obs  = circ_std(varpi0[incap])
sig_varpi_all  = circ_std(varpi0)

dv_obs = drift(vw_obs)
dv_scr = drift(vw_scr)
dv_p9  = drift(vw_p9)

def model_stats(dv):
    sd_all, sd_in = float(np.std(dv)), float(np.std(dv[incap]))
    return {
        "sigma_dv_all_deg_Myr": sd_all,
        "sigma_dv_incap_deg_Myr": sd_in,
        "t_disp_all_Myr": float(sig_varpi_all / sd_all) if sd_all > 0 else None,
        "t_disp_incap_Myr": float(sig_varpi_obs / sd_in) if sd_in > 0 else None,
        "med_dv_incap": float(np.median(dv[incap])),
        "med_dv_outcap": float(np.median(dv[~incap])),
    }

def R_series(vw):
    return [cluster_r(vw[k]) for k in range(vw.shape[0])]

# tide restoring test: does tide-added drift push varpi toward the
# tide's own stationary directions?  Tide ligation centres in ecliptic
# longitude; each object's offset is distance to the nearer centre.
tide_lons = [d[0] for d in TIDE_DIRS]
def nearest_tide_off(varpi_deg):
    d = min(abs((varpi_deg - l + 180) % 360 - 180) for l in tide_lons)
    # signed: positive drift toward the nearer centre is restoring
    l_near = min(tide_lons,
                 key=lambda l: abs((varpi_deg - l + 180) % 360 - 180))
    return -((varpi_deg - l_near + 180) % 360 - 180)  # deg to cover

off = np.array([nearest_tide_off(v) for v in varpi0])

# restoring diagnostic: under the tide, does each object's drift reduce
# |off| (move toward the nearer libration centre)?  Restoring
# confinement gives a negative correlation between |off| and the
# toward-centre drift (larger offsets restored more strongly).
toward = -np.sign(off) * dv_obs
rho_r, p_r = spearmanr(np.abs(off), toward)

res = {
    "method": "44 detached TNOs (a>150,q>30,cc<=3) under Sun+giants + "
              "epoch-fixed galactic tide (step_070 tensor, Strang-split "
              "kicks every 1 kyr), WHFast dt=0.5 yr, 20 Myr.  Models: "
              "real varpi (maintenance), scrambled varpi (regeneration), "
              "real varpi + BB21 median P9 (joint confinement).",
    "n_tno": len(tnos), "n_incap": int(incap.sum()),
    "cap_deg": 60.0, "axis": "49,-17",
    "T0_JD": T0_JD, "T_end_yr": T_END, "dt_yr": DT,
    "snap_yr": SNAP, "kick_yr": KICK,
    "tide_coeffs_yr2": {"Kx": KX, "Ky": KY, "Kz": KZ},
    "tide_libration_ecliptic_lon": tide_lons,
    "sig_varpi_incap_deg": sig_varpi_obs,
    "sig_varpi_all_deg": sig_varpi_all,
    "models": {
        "tide_observed":  model_stats(dv_obs),
        "tide_scrambled": model_stats(dv_scr),
        "tide_p9_med":    model_stats(dv_p9),
    },
    "tide_restoring": {
        "rho_abs_off_vs_toward": float(rho_r), "p_2sided": float(p_r),
        "med_toward_incap": float(np.median(toward[incap])),
        "med_toward_outcap": float(np.median(toward[~incap])),
        "note": "spearman of |offset from nearer tide libration centre| "
                "against drift toward that centre; restoring confinement "
                "requires a significant negative rho (larger offsets "
                "restored more strongly)",
    },
    "R_varpi_t": {
        "ts_Myr": [float(x) for x in (ts / 1e6)],
        "tide_observed":  R_series(vw_obs),
        "tide_scrambled": R_series(vw_scr),
        "tide_p9_med":    R_series(vw_p9),
    },
    "endpoints": {
        "tide_observed_R20":  cluster_r(vw_obs[-1]),
        "tide_observed_meanvarpi20": mean_dir(vw_obs[-1]),
        "tide_scrambled_R20": cluster_r(vw_scr[-1]),
        "tide_scrambled_meanvarpi20": mean_dir(vw_scr[-1]),
        "tide_scrambled_Rmax": max(R_series(vw_scr)),
        "tide_p9_R20": cluster_r(vw_p9[-1]),
        "tide_p9_meanvarpi20": mean_dir(vw_p9[-1]),
    },
}

out = str(RESULTS / "step_b53_tide_resident.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b53_tide_resident.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "a", "q", "i", "varpi0", "theta", "in_cap",
                "dv_tide_obs", "dv_tide_scr", "dv_tide_p9"])
    for j, t in enumerate(tnos):
        w.writerow([t["name"], t["a"], t["q"], t["i"], t["varpi"],
                    theta[j], bool(incap[j]),
                    dv_obs[j], dv_scr[j], dv_p9[j]])

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

ax = axes[0]
ax.plot(ts / 1e6, res["R_varpi_t"]["tide_observed"], color="crimson",
        label="giants+tide (real $\\varpi$)")
ax.plot(ts / 1e6, res["R_varpi_t"]["tide_p9_med"], color="darkorange",
        label="giants+tide+P9")
ax.plot(ts / 1e6, res["R_varpi_t"]["tide_scrambled"], color="0.4", ls="--",
        label="giants+tide (scrambled)")
ax.axhline(0.332, color="k", ls=":", lw=1)
ax.set_xlabel("time (Myr)")
ax.set_ylabel(r"cluster resultant $R_\varpi$ (all 44)")
ax.set_ylim(0, 1)
ax.legend(frameon=False, fontsize=8)
ax.set_title("tide on the resident cluster", fontsize=10)

ax = axes[1]
for arr, lab, c in [(dv_obs, "tide", "crimson"),
                    (dv_p9, "tide+P9", "darkorange"),
                    (dv_scr, "tide scr.", "0.5")]:
    ax.scatter(theta[~incap], arr[~incap], s=14, c=c, alpha=0.4)
    ax.scatter(theta[incap], arr[incap], s=18, c=c, alpha=0.9, label=lab)
ax.axvline(60, color="k", ls=":", lw=1)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel(r"periapsis--axis separation $\theta$ (deg)")
ax.set_ylabel(r"$d\varpi/dt$ (deg/Myr)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("apsidal drift under the tide", fontsize=10)

ax = axes[2]
ax.scatter(np.abs(off)[incap], toward[incap], s=18, c="crimson",
           alpha=0.9, label="in-cap")
ax.scatter(np.abs(off)[~incap], toward[~incap], s=14, c="0.4",
           alpha=0.4, label="out-of-cap")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel(r"$|\Delta\varpi|$ to nearer tide centre (deg)")
ax.set_ylabel(r"drift toward centre (deg/Myr)")
ax.legend(frameon=False, fontsize=8)
ax.set_title(f"tide restoring test ($\\rho$={rho_r:+.2f})", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b53_tide_resident.png", dpi=300)

for m in ("tide_observed", "tide_scrambled", "tide_p9_med"):
    d = res["models"][m]
    logger.info(f"{m}: sigma(dv) all/in = {d['sigma_dv_all_deg_Myr']:.3f}/"
                f"{d['sigma_dv_incap_deg_Myr']:.3f} deg/Myr | "
                f"t_disp in = {d['t_disp_incap_Myr']:.0f} Myr")
e = res["endpoints"]
logger.info(f"endpoints: R20 obs={e['tide_observed_R20']:.3f} "
            f"scr={e['tide_scrambled_R20']:.3f} "
            f"(mean {e['tide_scrambled_meanvarpi20']:.0f} deg) "
            f"p9={e['tide_p9_R20']:.3f}")
logger.info(f"tide restoring: rho={rho_r:+.3f} p={p_r:.4f}")
logger.data_save(out)
logger.data_save(csv_out)
logger.data_save(FIG / "supplementary" / "step_b53_tide_resident.png")