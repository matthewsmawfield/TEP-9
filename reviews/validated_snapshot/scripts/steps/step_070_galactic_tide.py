"""step_070: Galactic-tide insertion -- the last unsubtracted Newtonian
channel.

The planetary baseline (steps 040/042/063/064) subtracts the measured
kick budget of the Solar system, and the Planet Nine insertion (step
062) shows the published perturber cannot rotate the catalogued
trajectories.  One Newtonian channel remains unevaluated on the real
trajectories: the Galactic tide, the dominant long-term perturber of
Oort-cloud comets.  Its strength at the 250 AU boundary sphere is small
-- the vertical tide coefficient 4 pi G rho ~ 5.6e-15 yr^-2 against
the solar GM/r^2 ~ 6.3e-4 AU yr^-2 -- but the anomaly itself is also
small (~0.1 deg), so the tide's actual contribution is measured here
directly rather than asserted negligible.

Method
------
The identical bidirectional integration of step_063 is rerun with the
Heisler & Tremaine (1986) tidal tensor added as an external potential
via Strang splitting (kick-drift-kick on WHFast).  In the
epoch-fixed approximation -- each leg spans ~600 yr, vastly less than
the ~240 Myr galactic period -- the inertial-frame tensor is

    a_gal = ( +OmG^2 x,  -OmG^2 y,  -4 pi G rho z )

in galactic coordinates (x toward the Galactic centre, z toward the
north Galactic pole), with OmG = A - B = 26 km/s/kpc (A = -B = 13,
delta = 0) and rho = 0.1 Msun/pc^3 (Levison et al. 2001 convention).

Each comet's two legs are integrated with and without the tide; the
per-leg boundary periapsis-direction offset and the change in the
orig->fut rotation are measured on the catalogued trajectories, and
the cap contrast of the tide-induced contribution is tested.

Outputs
-------
results/step_b35_galactic_tide.json
results/step_b35_galactic_tide.csv
results/figures/step_b35_galactic_tide.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (DATA_RAW, RESULTS, ECL2GAL,
                                       GAL2ECL, lv, parse_code,
                                       tee_stdout)
logger = StepLogger("step_070_galactic_tide")
tee_stdout(logger)
logger.header("Galactic-tide insertion control (REBOUND/DE440s)")

import csv
import json
import math
import re
import numpy as np
import rebound
import spiceypy as sp
from html.parser import HTMLParser
from scipy.stats import mannwhitneyu

# ------------------------------------------------------------------
# CODE table parser (identical to step_063) -- keeps the full
# perihelion date string T needed for the integration epoch
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

def parse_code_full(path):
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

# ------------------------------------------------------------------
# Ephemeris setup (identical to step_063)
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
# Galactic tide tensor (epoch-fixed inertial approximation)
# ------------------------------------------------------------------

KMS_AUYR = 0.210805                    # 1 km/s in AU/yr
AU_KPC   = 206264.806 * 1000.0         # AU per kpc
OMG      = 26.0 * KMS_AUYR / AU_KPC    # Omega_G = A-B = 26 km/s/kpc [1/yr]
G_AU     = 4 * math.pi ** 2            # G in AU^3/(yr^2 Msun)
RHO_G    = 0.1 / (AU_KPC / 1000.0) ** 3  # 0.1 Msun/pc^3 in Msun/AU^3
KZ       = -4.0 * math.pi * G_AU * RHO_G   # -4 pi G rho  [1/yr^2]
KX, KY   = OMG ** 2, -OMG ** 2             # epoch-fixed radial terms

TIDE = np.diag([KX, KY, KZ])
logger.info(f"tide coefficients: Kx={KX:.3e} Ky={KY:.3e} "
            f"Kz={KZ:.3e} yr^-2 (z/radial ratio {abs(KZ/KX):.1f})")

MU = 4 * math.pi ** 2
R_STOP = 250.0
DT     = 0.05                          # WHFast step (yr)
T_MAX  = 3000.0

def perih_dir(om, Om, inc):
    co, so = math.cos(om), math.sin(om)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(inc), math.sin(inc)
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci, so * si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

TNO = lv(34.0, -13.0)                  # Warsaw-discovery axis

def init_sim(et):
    sim = rebound.Simulation()
    sim.G = MU
    ps, vs = body_state("10", et)
    sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2], m=1.0)
    for b in PLANET_IDS:
        pp, vv = body_state(b, et)
        sim.add(x=pp[0], y=pp[1], z=pp[2], vx=vv[0], vy=vv[1], vz=vv[2],
                m=GM[b] / GM_SUN)
    return sim

def state_at_periapsis(ro):
    w_, O_, i_ = map(math.radians, (ro["w"], ro["Om"], ro["i"]))
    ph = perih_dir(w_, O_, i_)
    hh = np.array([math.sin(i_) * math.sin(O_),
                   -math.sin(i_) * math.cos(O_), math.cos(i_)])
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

def integrate_leg(ro, et, direction, tide):
    """step_063 leg structure, WHFast + optional Strang-split tide."""
    sim = init_sim(et)
    rvec, vvec = state_at_periapsis(ro)
    ps, vs = body_state("10", et)
    sim.add(x=rvec[0] + ps[0], y=rvec[1] + ps[1], z=rvec[2] + ps[2],
            vx=vvec[0] + vs[0], vy=vvec[1] + vs[1], vz=vvec[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "whfast"
    sim.dt = direction * DT                     # signed step
    STRIDE = 100                                # kick every 5 yr
    t = 0.0
    reached = False
    while abs(t) < T_MAX:
        # KDK over each stride: half kick, STRIDE WHFast steps, half kick
        if tide:
            _kick(sim, nc, sim.dt * STRIDE / 2)
        sim.integrate(sim.t + sim.dt * STRIDE, exact_finish_time=1)
        t += direction * DT * STRIDE
        if tide:
            _kick(sim, nc, sim.dt * STRIDE / 2)
        p = sim.particles
        r_rel = math.sqrt((p[nc].x - p[0].x) ** 2 +
                          (p[nc].y - p[0].y) ** 2 +
                          (p[nc].z - p[0].z) ** 2)
        if r_rel >= R_STOP:
            reached = True
            break
    if not reached:
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
    return dict(phat=phat, aa=aa, t_years=float(t))

def _kick(sim, nc, dtk):
    """tidal velocity kick on the comet: v += a_tide(r_helio) * dtk."""
    p = sim.particles
    rh = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y, p[nc].z - p[0].z])
    ag = ECL2GAL @ rh                       # heliocentric pos, gal frame
    at = TIDE @ ag                          # tidal accel, gal frame
    ae = GAL2ECL @ at                       # back to ecliptic
    p[nc].vx += ae[0] * dtk
    p[nc].vy += ae[1] * dtk
    p[nc].vz += ae[2] * dtk

# ------------------------------------------------------------------
# Sample: identical to step_063/062 -- class-1 CODE comets
# ------------------------------------------------------------------

osc  = parse_code_full(str(DATA_RAW / "code" / "code_osculating.html"))
orig = parse_code_full(str(DATA_RAW / "code" / "code_original.html"))
fut  = parse_code_full(str(DATA_RAW / "code" / "code_future.html"))
warsaw = {l[5:17].strip() for l in
          open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
          if len(l) > 115}

sample = []
for k, ro in orig.items():
    if k not in fut:
        continue
    if not (0 < ro["aa"] < 100 and ro["cls"] in ("1a", "1a+", "1b")) \
            or k in warsaw:
        continue
    if k not in osc:
        continue
    et = perihelion_et(osc[k]["T"])
    sample.append((k, ro, osc[k], et))

logger.info(f"sample: {len(sample)} class-1 CODE comets")

# ------------------------------------------------------------------
# Integrate both legs, tide on/off
# ------------------------------------------------------------------

rows = []
failed = []
for k, ro, oo, et in sample:
    try:
        rb0 = integrate_leg(oo, et, -1, False)
        rf0 = integrate_leg(oo, et, +1, False)
        rbT = integrate_leg(oo, et, -1, True)
        rfT = integrate_leg(oo, et, +1, True)
    except Exception as e:
        failed.append(k); logger.warning(f"{k}: {e}"); continue
    if None in (rb0, rf0, rbT, rfT):
        failed.append(k); logger.warning(f"{k}: boundary not reached")
        continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                   math.radians(ro["i"]))
    dleg_b = sep(rb0["phat"], rbT["phat"])   # tide-induced offset, back leg
    dleg_f = sep(rf0["phat"], rfT["phat"])   # ... forward leg
    drot0 = sep(rb0["phat"], rf0["phat"])
    drotT = sep(rbT["phat"], rfT["phat"])
    rows.append(dict(desig=k, q=ro["q"], i=ro["i"],
                     theta=sep(-po, TNO),
                     dleg_back=dleg_b, dleg_fwd=dleg_f,
                     drot_plain=drot0, drot_tide=drotT,
                     ddrot=drotT - drot0))
    if len(rows) % 20 == 0:
        logger.info(f"  {len(rows)}/{len(sample)}")

logger.info(f"integrated {len(rows)} comets, {len(failed)} failures")

th   = np.array([r["theta"] for r in rows])
inc  = th < 60.0
dlb  = np.array([r["dleg_back"] for r in rows])
dlf  = np.array([r["dleg_fwd"] for r in rows])
dmax = np.maximum(dlb, dlf)
dd   = np.array([r["ddrot"] for r in rows])
dr0  = np.array([r["drot_plain"] for r in rows])
drT  = np.array([r["drot_tide"] for r in rows])

res = {
    "method": "step_063 bidirectional legs rerun as WHFast with the "
              "Heisler-Tremaine tidal tensor added by Strang splitting "
              "(epoch-fixed inertial approximation; A=-B=13 km/s/kpc, "
              "rho=0.1 Msun/pc^3).  Per-leg boundary periapsis offset "
              "and orig->fut rotation change measured tide-on vs off.",
    "tide_coeffs_yr2": {"Kx": KX, "Ky": KY, "Kz": KZ},
    "dt_yr": DT, "n": len(rows), "n_in": int(inc.sum()),
    "n_failed": len(failed), "failed": failed,
    "per_leg_tide_offset_deg": {
        "median": float(np.median(dmax)),
        "p95": float(np.percentile(dmax, 95)),
        "max": float(dmax.max())},
    "rotation_change_deg": {
        "median_abs": float(np.median(np.abs(dd))),
        "p95_abs": float(np.percentile(np.abs(dd), 95)),
        "max_abs": float(np.abs(dd).max()),
        "median_signed_in": float(np.median(dd[inc])),
        "median_signed_out": float(np.median(dd[~inc]))},
    "cap_contrast": {
        "leg_offset_in_gt_out_p": float(
            mannwhitneyu(dmax[inc], dmax[~inc],
                         alternative="greater").pvalue),
        "abs_ddrot_in_gt_out_p": float(
            mannwhitneyu(np.abs(dd[inc]), np.abs(dd[~inc]),
                         alternative="greater").pvalue)},
    "validation": {
        "drot_plain_median": float(np.median(dr0)),
        "drot_tide_median": float(np.median(drT)),
        "observed_incap_median_deg": 0.14},
}

out = str(RESULTS / "step_b35_galactic_tide.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b35_galactic_tide.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axes[0]
ax.scatter(th[~inc], dmax[~inc], s=20, facecolors="none",
           edgecolors="0.55", label=f"outside ($n={int((~inc).sum())}$)")
ax.scatter(th[inc], dmax[inc], s=24, c="crimson",
           label=f"inside ($n={int(inc.sum())}$)")
ax.axvline(60, color="k", ls=":", lw=1)
ax.set_yscale("log")
ax.set_xlabel(r"aphelion--axis separation $\theta$ (deg)")
ax.set_ylabel("tide-induced leg offset (deg)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("per-leg tide contribution vs axis distance", fontsize=10)

ax = axes[1]
bins = np.linspace(-np.percentile(np.abs(dd), 99),
                   np.percentile(np.abs(dd), 99), 40)
ax.hist(dd[~inc], bins=bins, color="0.55", alpha=0.7, density=True,
        label="outside")
ax.hist(dd[inc], bins=bins, color="crimson", alpha=0.6, density=True,
        label="inside")
ax.axvline(0, color="k", ls=":", lw=1)
ax.set_xlabel(r"$\Delta$ rotation from adding the tide (deg)")
ax.set_ylabel("density")
ax.legend(frameon=False, fontsize=8)
ax.set_title("tide contribution to orig$\\to$fut rotation", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b35_galactic_tide.png", dpi=150)

logger.info(f"per-leg tide offset: med {np.median(dmax):.2e} deg, "
            f"max {dmax.max():.2e} deg")
logger.info(f"rotation change: med |ddrot| {np.median(np.abs(dd)):.2e}, "
            f"max {np.abs(dd).max():.2e} deg")
logger.info(f"cap contrast: leg p={res['cap_contrast']['leg_offset_in_gt_out_p']:.3f}, "
            f"ddrot p={res['cap_contrast']['abs_ddrot_in_gt_out_p']:.3f}")
print("wrote", out)
print("wrote", csv_out)
print("wrote", FIG / "step_b35_galactic_tide.png")
