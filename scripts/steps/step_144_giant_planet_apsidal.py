"""Step 144: Giant-planet secular apsidal channel (step_b108).

The resident axis (azimuth 49 deg) sits 4.0 deg from Neptune's mean
longitude of perihelion (~45 deg; step_112 registered the coincidence
at p ~ 0.02 a priori).  A mean-element coincidence is only one slice
of the secular structure: osculating apsidal longitudes circulate and
oscillate under the giant planets' mutual perturbations, and the
time-averaged eccentricity vector -- the forced apsidal direction of
each planet -- is the quantity a fixed external boundary could shape.
If the outer-system field configuration fixes a preferred direction,
the outermost planets (deepest coupling to the boundary region)
should carry the strongest alignment.

This step integrates Sun + Jupiter..Neptune under DE440s/REBOUND
(WHFAST) for 5 Myr and reads each giant's apsidal vector.

T1  Forced apsidal azimuth per giant: arg of the time-averaged
    complex eccentricity vector <e exp(i varpi)>, its separation from
    the resident-axis azimuth, and the forced/free amplitude ratio.
T2  Joint test: the four forced azimuths' circular mean and the
    axis's proximity to it, priced against 20,000 random target
    longitudes (uniform a priori null).
T3  Dwell fraction: the fraction of the integration each giant's
    instantaneous apsidal longitude spends within +/-15 deg of the
    axis azimuth, against the 15/180 uniform expectation.
T4  Present-epoch 3-D geometry: each giant's osculating periapsis
    unit vector separation from the 3-D axis, alongside the
    mean-element coincidences (Neptune 4.0 deg) for the record.
T5  Inter-giant apsidal coherence: pairwise forced-direction
    separations -- the internal apsidal architecture the axis
    coincidence sits inside.

Outputs: results/step_b108_giant_apsidal.json/.csv and
results/figures/step_b108_giant_apsidal.png.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout, lv, sep
logger = StepLogger("step_144_giant_apsidal")
tee_stdout(logger)
logger.header("Giant-planet secular apsidal channel")

import json
import math
import csv
import numpy as np
import rebound
import spiceypy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20260919
rng = np.random.default_rng(SEED)

SPK = DATA_RAW / "spice" / "de440s.bsp"
sp.furnsh(str(SPK))
sp.furnsh(str(DATA_RAW / "naif" / "naif0012.tls"))

AU_KM = 149597870.7
DAY_YR = 365.25
EPS = math.radians(23.4392911)
RX = np.array([[1, 0, 0],
               [0, math.cos(EPS), math.sin(EPS)],
               [0, -math.sin(EPS), math.cos(EPS)]])

MU = 4 * math.pi ** 2                      # G in AU^3 M_sun^-1 yr^-2
GM_SUN = 1.32712440018e11
GM = {"5": 1.2671276480e8, "6": 3.7940626000e7,
      "7": 5.7945490100e6, "8": 6.8365271006e6}   # J S U N km^3/s^2
NAMES = {"5": "Jupiter", "6": "Saturn", "7": "Uranus", "8": "Neptune"}
T0_JD = 2461200.5
ET0 = (T0_JD - 2451545.0) * 86400.0
AXIS_LON = 49.0
AXIS = lv(49.0, -17.0)
N_AXES_NULL = 20000

T_END = 5.0e6       # covers the longest giant secular periods
DT = 0.5
SNAP = 1.0e3        # snapshot cadence yr


def body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return (RX @ np.array(st[:3]) / AU_KM,
            RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR)


sim = rebound.Simulation()
sim.G = MU
ps, vs = body_state("10", ET0)
sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2], m=1.0)
for b in GM:
    pp, vv = body_state(b, ET0)
    sim.add(x=pp[0], y=pp[1], z=pp[2],
            vx=vv[0], vy=vv[1], vz=vv[2], m=GM[b] / GM_SUN)
sim.integrator = "whfast"
sim.dt = DT

n_snap = int(T_END / SNAP) + 1
exo = {b: np.zeros(n_snap, dtype=complex) for b in GM}
peri3d = {}
ts = np.arange(n_snap) * SNAP

print(f"integrating {T_END/1e6:.0f} Myr, {n_snap} snapshots...")
for k in range(n_snap):
    t = ts[k]
    if k:
        sim.integrate(t * DT / DT, exact_finish_time=1)
    for i, b in enumerate(GM):
        p = sim.particles[i + 1]
        orb = p.orbit(primary=sim.particles[0])
        e = orb.e
        varpi = (orb.Omega + orb.omega) % (2 * math.pi)
        exo[b][k] = e * np.exp(1j * varpi)
        if k == n_snap - 1:
            # 3-D periapsis unit vector at final epoch
            from scripts.utils.tep9_common import perih_dir
            peri3d[b] = perih_dir(orb.omega, orb.Omega, orb.inc)
print("integration done.")

res = {"t_end_yr": T_END, "dt_yr": DT, "snap_yr": SNAP, "seed": SEED,
       "axis_lon_deg": AXIS_LON, "epoch_jd": T0_JD}

# ---- T1: forced apsidal azimuth per giant -----------------------------------
res["T1_forced_apsidal"] = {}
for b in GM:
    z = exo[b]
    zm = z.mean()
    az = math.degrees(math.atan2(zm.imag, zm.real)) % 360.0
    dsep = abs(((az - AXIS_LON + 180) % 360) - 180)
    # free amplitude: scatter about the mean vector
    free_amp = np.abs(z - zm).mean()
    res["T1_forced_apsidal"][NAMES[b]] = {
        "forced_azimuth_deg": az,
        "sep_from_axis_deg": float(dsep),
        "forced_amp": float(abs(zm)),
        "free_amp": float(free_amp),
        "forced_over_free": float(abs(zm) / max(free_amp, 1e-12)),
        "e_mean": float(np.abs(z).mean())}
    print(f"{NAMES[b]:8s} forced apsidal az={az:7.2f} deg "
          f"sep={dsep:6.2f} deg  e_forced={abs(zm):.4f} "
          f"e_free={free_amp:.4f}")

# ---- T2: joint alignment vs axis ---------------------------------------------
azs = np.array([res["T1_forced_apsidal"][NAMES[b]]["forced_azimuth_deg"]
                for b in GM])
# weight each planet's forced direction by its forced amplitude
w = np.array([res["T1_forced_apsidal"][NAMES[b]]["forced_amp"] for b in GM])
zw = (w * np.exp(1j * np.radians(azs))).sum()
mean_dir = math.degrees(math.atan2(zw.imag, zw.real)) % 360
R = abs(zw) / w.sum()
sep_axis = abs(((mean_dir - AXIS_LON + 180) % 360) - 180)
null_dir = rng.uniform(0, 360, N_AXES_NULL)
frac_closer = float((np.minimum(np.abs(null_dir - mean_dir),
                                360 - np.abs(null_dir - mean_dir))
                     <= sep_axis).mean())
res["T2_joint"] = {
    "weighted_mean_forced_azimuth_deg": float(mean_dir),
    "weighted_resultant": float(R),
    "sep_from_axis_deg": float(sep_axis),
    "uniform_p_axis_this_close": frac_closer,
    "unweighted_mean_azimuth_deg": float(
        math.degrees(math.atan2(np.sin(np.radians(azs)).mean(),
                               np.cos(np.radians(azs)).mean())) % 360),
    "note": "forced directions are secular-eigenmode driven and "
            "internally coherent by construction; the measured "
            "quantity is whether the common forced direction lands "
            "on the resident axis"}
print(f"T2 weighted forced azimuth={mean_dir:.2f} (R={R:.3f}), "
      f"sep from axis={sep_axis:.2f}, uniform p={frac_closer:.4f}")

# ---- T3: dwell fraction within +/-15 deg --------------------------------------
res["T3_dwell"] = {}
halfw = 15.0
for b in GM:
    az_t = np.degrees(np.angle(exo[b])) % 360.0
    d = np.minimum(np.abs(az_t - AXIS_LON), 360 - np.abs(az_t - AXIS_LON))
    frac = float((d <= halfw).mean())
    res["T3_dwell"][NAMES[b]] = {
        "frac_within_15deg": frac,
        "uniform_expectation": halfw / 180.0,
        "ratio": frac / (halfw / 180.0)}
    print(f"{NAMES[b]:8s} dwell within 15 deg of axis: {frac:.3f} "
          f"(uniform {halfw/180:.3f})")

# ---- T4: present-epoch 3-D periapsis geometry ---------------------------------
res["T4_present_3d"] = {}
for b in GM:
    s = sep(peri3d[b], AXIS)
    res["T4_present_3d"][NAMES[b]] = {
        "periapsis_sep_from_axis_deg": float(s)}
    print(f"{NAMES[b]:8s} 3-D periapsis sep from axis: {s:.2f} deg")
res["T4_present_3d"]["registered_mean_element_coincidence"] = {
    "neptune_mean_varpi_deg": 44.965,
    "sep_deg": 4.035, "p_a_priori": 0.0224,
    "source": "step_112 resonant control; mean element (Standish 2006)"}

# ---- T5: inter-giant apsidal coherence -----------------------------------------
res["T5_inter_giant"] = {}
for i, b1 in enumerate(GM):
    for b2 in list(GM)[i + 1:]:
        a1 = res["T1_forced_apsidal"][NAMES[b1]]["forced_azimuth_deg"]
        a2 = res["T1_forced_apsidal"][NAMES[b2]]["forced_azimuth_deg"]
        d = abs(((a1 - a2 + 180) % 360) - 180)
        res["T5_inter_giant"][f"{NAMES[b1]}-{NAMES[b2]}"] = {
            "forced_sep_deg": float(d)}

# ---- csv: per-giant apsidal time series (decimated) -----------------------------
with open(RESULTS / "step_b108_giant_apsidal.csv", "w", newline="") as fh:
    wtr = csv.writer(fh)
    wtr.writerow(["t_yr"] + [NAMES[b] + "_varpi_deg" for b in GM])
    for k in range(0, n_snap, 20):
        wtr.writerow([ts[k]] +
                     [float(np.degrees(np.angle(exo[b][k])) % 360)
                      for b in GM])

# ---- verdict --------------------------------------------------------------------
nep = res["T1_forced_apsidal"]["Neptune"]
outer_sep = nep["sep_from_axis_deg"]
verdict = (
    "GIANT-APSIDAL CHANNEL " +
    ("AXIS-ALIGNED" if frac_closer < 0.05 else "QUIET") +
    f": weighted mean forced apsidal azimuth {mean_dir:.1f} deg "
    f"({sep_axis:.1f} deg from the resident axis, uniform "
    f"p={frac_closer:.4f}); Neptune forced azimuth "
    f"{nep['forced_azimuth_deg']:.1f} deg ({outer_sep:.1f} deg off), "
    f"forced/free ratio {nep['forced_over_free']:.2f}.  The mean-element "
    "coincidence (Neptune varpi 4.0 deg, p~0.02) is registered with "
    "its secular context: whether the alignment survives as a fixed "
    "forced direction or is a circulating coincidence is decided here "
    "by the 5 Myr integration.")
res["verdict"] = verdict
res["test_summary"] = {
    "joint_sep_from_axis_deg": float(sep_axis),
    "joint_uniform_p": frac_closer,
    "neptune_forced_sep_deg": float(outer_sep),
    "neptune_dwell_ratio": res["T3_dwell"]["Neptune"]["ratio"],
    "axis_aligned": bool(frac_closer < 0.05),
}

out = RESULTS / "step_b108_giant_apsidal.json"
json.dump(res, open(out, "w"), indent=1)

# ---- figure ---------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
ax = axes[0]
cols = {"5": "orange", "6": "goldenrod", "7": "steelblue", "8": "crimson"}
for b in GM:
    az_t = np.degrees(np.angle(exo[b])) % 360
    ax.plot(ts[::20] / 1e6, az_t[::20], lw=0.4, alpha=0.7,
            color=cols[b], label=NAMES[b])
ax.axhline(AXIS_LON, color="k", ls="--", lw=1, label="resident axis")
ax.set_xlabel("Myr"); ax.set_ylabel("osculating varpi (deg)")
ax.legend(fontsize=7, frameon=False)
ax.set_title("apsidal circulation")

ax = axes[1]
for b in GM:
    z = exo[b]
    zm = z.mean()
    ax.arrow(0, 0, zm.real, zm.imag, head_width=0.004,
             color=cols[b], length_includes_head=True,
             label=NAMES[b])
axa = AXIS_LON * math.pi / 180
ax.arrow(0, 0, 0.06 * math.cos(axa), 0.06 * math.sin(axa),
         head_width=0.006, color="k", ls="--",
         length_includes_head=True, label="axis (scaled)")
ax.set_aspect("equal")
ax.set_xlabel("e cos(varpi)"); ax.set_ylabel("e sin(varpi)")
ax.legend(fontsize=7, frameon=False)
ax.set_title("forced eccentricity vectors")

ax = axes[2]
names = [NAMES[b] for b in GM]
seps = [res["T1_forced_apsidal"][n]["sep_from_axis_deg"] for n in names]
ax.bar(names, seps, color=[cols[b] for b in GM])
ax.axhline(90, color="k", ls=":", lw=1)
ax.set_ylabel("forced azimuth sep from axis (deg)")
ax.set_title("per-giant axis proximity")
fig.tight_layout()
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b108_giant_apsidal.png", dpi=150)

logger.info("verdict: " + res["verdict"])
logger.data_save(out)
logger.data_save(RESULTS / "step_b108_giant_apsidal.csv")
logger.data_save(FIG / "step_b108_giant_apsidal.png")
print(json.dumps(res["test_summary"], indent=1))
print(res["verdict"])
