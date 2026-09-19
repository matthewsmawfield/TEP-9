#!/usr/bin/env python3
"""Step 071 -- radial slip profile: where does the anomalous rotation
accumulate?

The step-063 bidirectional legs are extended to a set of boundary shells
(60-400 AU).  On each leg the osculating barycentric periapsis direction
is recorded at every shell crossing, so the inbound/outbound rotation
delta_theta(s) is measured as a function of radius rather than only at
250 AU.  For a discrete phase slip applied at a crossing radius r*,
delta_theta(s) is flat inside r* (baseline only) and offsets by the slip
angle for s > r* -- the kink radius localizes the boundary; for a
distributed field delta_theta(s) grows continuously.  The implied
proper-time profile delta_tau(s) = delta_theta(s) * s^2 / h converts
each shell measurement to time units.

Outputs:
    results/step_b36_radial_slip.json
    results/step_b36_radial_slip.csv      (per comet per shell)
    results/figures/step_b36_radial_slip.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_071_radial_slip_profile")
tee_stdout(logger)
logger.header("Radial slip profile: multi-shell bidirectional legs (REBOUND/DE440s)")

import csv
import json
import math
import re
import numpy as np
import rebound
import spiceypy as sp
from html.parser import HTMLParser
from scipy.stats import spearmanr

# ------------------------------------------------------------------
# CODE table parser (identical to step_063)
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

TNO = lv(34, -13)          # pre-declared transit/Warsaw axis (steps 030-036, 063-070)
CAP = 60.0

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
# Multi-shell integration
# ------------------------------------------------------------------

MU = 4 * math.pi ** 2
SHELLS = [60.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0,
          250.0, 300.0, 350.0, 400.0]
R_MAX  = SHELLS[-1]
T_MAX  = 40000.0
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
    return phat, -2 * E / mu * 1e6, np.linalg.norm(h)

def integrate_leg(ro, et, direction):
    """Integrate one leg; record the barycentric osculating state at the
    first crossing of every shell in SHELLS.  Returns dict shell ->
    (phat, h, t) or None."""
    sim, ps, vs = init_sim(et)
    rvec, vvec = state_at_periapsis(ro)
    sim.add(x=rvec[0] + ps[0], y=rvec[1] + ps[1], z=rvec[2] + ps[2],
            vx=vvec[0] + vs[0], vy=vvec[1] + vs[1], vz=vvec[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    denc = np.full(sim.N - 1, np.inf)
    hits = {}
    nxt = 0                       # index into SHELLS (crossed in order)
    t = direction * DT_OUT
    while abs(t) < T_MAX:
        sim.integrate(t, exact_finish_time=0)
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y, p[nc].z - p[0].z])
        rr = np.linalg.norm(r_rel)
        for j in range(1, sim.N - 1):
            d = math.sqrt((p[nc].x - p[j].x) ** 2 + (p[nc].y - p[j].y) ** 2 +
                          (p[nc].z - p[j].z) ** 2)
            if d < denc[j - 1]: denc[j - 1] = d
        if nxt < len(SHELLS) and rr >= SHELLS[nxt]:
            # barycentric osculating elements at this crossing
            mtot = sum(pp.m for pp in sim.particles)
            rb = np.zeros(3); vb = np.zeros(3)
            for pp in sim.particles:
                rb += pp.m * np.array([pp.x, pp.y, pp.z])
                vb += pp.m * np.array([pp.vx, pp.vy, pp.vz])
            rb /= mtot; vb /= mtot
            r_b = np.array([p[nc].x, p[nc].y, p[nc].z]) - rb
            v_b = np.array([p[nc].vx, p[nc].vy, p[nc].vz]) - vb
            phat, aa, hh = boundary_orbit(r_b, v_b, mtot)
            while nxt < len(SHELLS) and rr >= SHELLS[nxt]:
                hits[SHELLS[nxt]] = (phat, hh, aa, float(t))
                nxt += 1
            if nxt >= len(SHELLS):
                break
        t += direction * DT_OUT
    if len(hits) < len(SHELLS):
        return None
    return dict(hits=hits, denc=float(denc.min()))

# ------------------------------------------------------------------
# Sample: identical selection to step_063
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

logger.info(f"sample: {len(sample)} class-1 CODE comets; shells {SHELLS[0]:.0f}-{SHELLS[-1]:.0f} AU")

# ------------------------------------------------------------------
# Integrate both legs, all shells
# ------------------------------------------------------------------

rows = []
failed = []
for idx, (k, ro, oo, et) in enumerate(sample):
    try:
        rb = integrate_leg(oo, et, -1)
        rf = integrate_leg(oo, et, +1)
    except Exception as e:
        failed.append(k); logger.warning(f"{k}: {e}"); continue
    if rb is None or rf is None:
        failed.append(k); logger.warning(f"{k}: outer shell not reached"); continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    rec = dict(desig=k, q=ro["q"], i=ro["i"],
               theta=sep(-po, TNO),
               denc=min(rb["denc"], rf["denc"]),
               drot_cat=sep(-po, -pf))
    for s in SHELLS:
        pb, hb, ab_, tb = rb["hits"][s]
        pf_, hf, af_, tf = rf["hits"][s]
        rec[f"drot_{int(s)}"] = sep(pb, pf_)
        rec[f"h_{int(s)}"] = 0.5 * (hb + hf)
    rows.append(rec)
    if (idx + 1) % 20 == 0:
        logger.info(f"  {idx+1}/{len(sample)}")

logger.info(f"integrated {len(rows)} comets x {len(SHELLS)} shells, {len(failed)} failures")

# ------------------------------------------------------------------
# Profiles: median rotation and implied slip vs shell radius
# ------------------------------------------------------------------

th = np.array([r["theta"] for r in rows])
inc = th < CAP
rng = np.random.default_rng(20260918)
NB = 20000

def med_boot(vals, n=inc.sum()):
    if len(vals) < 5: return np.array([np.nan, np.nan])
    idx = rng.integers(0, len(vals), (NB, len(vals)))
    bs = np.median(np.asarray(vals)[idx], axis=1)
    return np.percentile(bs, [16, 84])

prof = {"shells": SHELLS, "n_in": int(inc.sum()), "n_out": int((~inc).sum()),
        "median_drot_in": [], "median_drot_out": [],
        "median_dtau_in": [], "median_dtau_out": [],
        "excess_drot": [], "excess_dtau": [], "excess_dtau_ci": []}

for s in SHELLS:
    tag = int(s)
    dr = np.array([r[f"drot_{tag}"] for r in rows])
    hh = np.array([r[f"h_{tag}"] for r in rows])
    dt = np.radians(dr) * s * s / hh          # yr, implied slip at shell s
    prof["median_drot_in"].append(float(np.median(dr[inc])))
    prof["median_drot_out"].append(float(np.median(dr[~inc])))
    prof["median_dtau_in"].append(float(np.median(dt[inc])))
    prof["median_dtau_out"].append(float(np.median(dt[~inc])))
    prof["excess_drot"].append(float(np.median(dr[inc]) - np.median(dr[~inc])))
    ex = dt[inc]                              # excess slip per comet vs out median
    med_out = np.median(dt[~inc])
    prof["excess_dtau"].append(float(np.median(ex) - med_out))
    lo, hi = med_boot(ex - med_out)
    prof["excess_dtau_ci"].append([float(lo), float(hi)])

# ---- onset radius: smallest shell where the in-cap excess exceeds
# ---- its inner-shell baseline by > 1 dex-sigma of the out-of-cap
# ---- scatter; also the implied-slip growth exponent beyond onset.
ex = np.array(prof["excess_drot"])
onset = None
for j, s in enumerate(SHELLS):
    inner = ex[:j]
    if len(inner) >= 3 and ex[j] > np.median(inner) + 2 * (np.std(inner) + 1e-6):
        onset = s; break
prof["onset_shell_au"] = onset

# growth exponent of the in-cap median rotation beyond 100 AU:
# wall -> delta_theta ~ const (b~0); continuous accumulation -> b>0
sh = np.array(SHELLS)
mask = sh >= 100.0
din = np.array(prof["median_drot_in"])[mask]
dout = np.array(prof["median_drot_out"])[mask]
for lab, dd in (("in", din), ("out", dout)):
    rho, p = spearmanr(sh[mask], dd)
    prof[f"drot_vs_shell_{lab}"] = {"rho": float(rho), "p_2sided": float(p)}
    b = np.polyfit(np.log(sh[mask]), np.log(np.maximum(dd, 1e-6)), 1)
    prof[f"drot_loglog_{lab}"] = {"b": float(b[0])}

# implied-slip plateau beyond 150 AU (wall -> dtau grows as s^2 only if
# the slip is measured at the wrong radius; the angle is what is fixed)
tin = np.array(prof["median_dtau_in"])[mask]
rho, p = spearmanr(sh[mask], tin)
prof["dtau_vs_shell_in"] = {"rho": float(rho), "p_2sided": float(p)}

res = {
    "method": "step_063 legs extended to shells 60-400 AU; osculating "
              "barycentric periapsis direction recorded at every shell "
              "crossing on each leg.  delta_theta(s) = separation of "
              "inbound/outbound leg periapsis at shell s; implied slip "
              "dtau(s) = delta_theta(s) * s^2 / h.  A discrete crossing "
              "offsets delta_theta by the slip angle for all s > r* "
              "(kink + plateau); a distributed field grows it "
              "continuously.  Median profiles in/out of the 60-deg cap, "
              "bootstrap 16-84% bands (seed 20260918).",
    "cap_deg": CAP, "axis": "34,-13 (pre-declared transit axis)",
    "seed": 20260918, "n": len(rows), "n_failed": len(failed),
    "failed": failed,
    "profile": prof,
}

# per-comet per-shell CSV
fields = ["desig", "q", "i", "theta", "denc", "drot_cat"] + \
         [f"drot_{int(s)}" for s in SHELLS] + [f"h_{int(s)}" for s in SHELLS]
with open(RESULTS / "step_b36_radial_slip.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
    for r in rows: w.writerow(r)

with open(RESULTS / "step_b36_radial_slip.json", "w") as f:
    json.dump(res, f, indent=1)
print(f"wrote {RESULTS / 'step_b36_radial_slip.json'}")
print(f"wrote {RESULTS / 'step_b36_radial_slip.csv'}")

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

ax = axes[0]
ax.plot(SHELLS, prof["median_drot_in"], "o-", color="crimson",
        label=f"inside cap ($n={int(inc.sum())}$)")
ax.plot(SHELLS, prof["median_drot_out"], "o-", color="0.55",
        label=f"outside ($n={int((~inc).sum())}$)")
if onset:
    ax.axvline(onset, color="k", ls=":", lw=1, label=f"onset {onset:.0f} AU")
ax.set_xlabel("shell radius $s$ (AU)")
ax.set_ylabel(r"median $\delta\theta(s)$ (deg)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("inbound/outbound rotation vs radius", fontsize=10)

ax = axes[1]
ci = np.array(prof["excess_dtau_ci"])
ax.fill_between(SHELLS, ci[:, 0], ci[:, 1], color="crimson", alpha=0.15)
ax.plot(SHELLS, prof["excess_dtau"], "o-", color="crimson")
ax.axhline(0, color="k", ls=":", lw=1)
ax.set_xlabel("shell radius $s$ (AU)")
ax.set_ylabel(r"in-cap excess $\delta\tau(s)$ (yr)")
ax.set_title("implied proper-time excess vs radius", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b36_radial_slip.png", dpi=150)
print(f"wrote {FIG / 'step_b36_radial_slip.png'}")

for j, s in enumerate(SHELLS):
    logger.info(f"s={s:6.0f} AU  drot in/out = {prof['median_drot_in'][j]:.4f}/"
                f"{prof['median_drot_out'][j]:.4f} deg  excess dtau = "
                f"{prof['excess_dtau'][j]:+.2f} yr")
logger.info(f"onset shell: {onset}")
