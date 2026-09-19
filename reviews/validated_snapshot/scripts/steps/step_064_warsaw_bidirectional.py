"""step_064: Warsaw bidirectional boundary validation -- independent
catalogue replication of step_063.

step_063 showed that an independent REBOUND/IAS15 integration of the
class-1 CODE comets reproduces the catalogue's original/future boundary
solutions to a median 0.001 deg per leg, and that the rotation-per-kick
residual retains its cap preference after the full encounter budget is
regressed out.  That validation was confined to a single catalogue.
This step repeats the identical machinery on the Warsaw near-parabolic
catalogue (Krolikowska 2014, VizieR J/A+A/567/A126): a different comet
sample, independently determined three-leg solutions.

For each Warsaw spike comet (0 < original 1/a < 100 x 10^-6 AU^-1,
matching the Phase-3 spike definition) the osculating elements at the
catalogued perihelion epoch are integrated through Sun +
planetary-system barycentres (DE440s, REBOUND/IAS15) twice:

  backward until barycentric |r| = 250 AU  -> analogue of "original"
  forward  until barycentric |r| = 250 AU  -> analogue of "future"

Boundary solutions are evaluated relative to the system barycentre,
matching the catalogue's original/future convention.

Outputs
-------
results/step_b29_warsaw_bidirectional.json
results/step_b29_warsaw_bidirectional.csv
results/figures/step_b29_warsaw_rotation_residual.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_064_warsaw_bidirectional")
tee_stdout(logger)
logger.header("Warsaw bidirectional boundary validation (REBOUND/IAS15 + DE440s)")

import csv
import json
import math
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

import rebound
import spiceypy as sp
from scripts.utils.coordinates import angular_separation

# ------------------------------------------------------------------
# Warsaw table parsers (fixed-width; perihelion time at bytes 28-42)
# ------------------------------------------------------------------

PREF     = {"a": 0, "h": 0, "e": 1, "b": 2}
PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}

def parse_orbit_table(path):
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(
                sample=line[0:2].strip(), com=line[3].strip(),
                desig=line[5:17].strip(),
                tyr=int(line[27:31]), tmo=int(line[31:33]),
                tdy=float(line[33:42]),
                q=float(line[42:56]), e=float(line[56:70]),
                w=float(line[70:82]), Om=float(line[82:94]),
                i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows

def dedup(rows, pref):
    out = {}
    for r in rows:
        k = r["desig"]
        if k not in out or pref.get(r["com"], 9) < pref.get(out[k]["com"], 9):
            out[k] = r
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

def perihelion_et(ro):
    return (jd_tt(ro["tyr"], ro["tmo"], ro["tdy"]) - 2451545.0) * 86400.0

def body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return RX @ np.array(st[:3]) / AU_KM, RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR

# ------------------------------------------------------------------
# REBOUND integration -- both legs, barycentric boundary elements
# ------------------------------------------------------------------

MU = 4 * math.pi ** 2
R_STOP = 250.0
T_MAX  = 20000.0
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
# Sample: Warsaw spike (0 < orig 1/a < 100 x 1e-6), all three legs
# ------------------------------------------------------------------

orig_r = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
osc_r  = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tableb.dat"))
fut_r  = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tabled.dat"))

orig = dedup(orig_r, PREF)
fut  = dedup(fut_r, PREF)
osc  = dedup(osc_r, PREF_OSC)

sample = []
for k, ro in orig.items():
    if not (0 < ro["aa"] < 100):
        continue
    if k not in osc or k not in fut:
        continue
    sample.append((k, ro, osc[k], fut[k]))

logger.info(f"sample: {len(sample)} Warsaw spike comets")

# ------------------------------------------------------------------
# Integrate both legs
# ------------------------------------------------------------------

rows = []
failed = []
for k, ro, oo, fo in sample:
    try:
        et = perihelion_et(oo)
        rb = integrate_leg(oo, et, -1)
        rf = integrate_leg(oo, et, +1)
    except Exception as e:
        failed.append(k); logger.warning(f"{k}: {e}"); continue
    if rb is None or rf is None:
        failed.append(k); logger.warning(f"{k}: boundary not reached"); continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fo["w"]), math.radians(fo["Om"]), math.radians(fo["i"]))
    ps = perih_dir(math.radians(oo["w"]), math.radians(oo["Om"]), math.radians(oo["i"]))
    drot_cat = sep(-po, -pf)
    drot_sim = sep(rb["phat"], rf["phat"])
    rows.append(dict(desig=k, com=oo["com"], q=ro["q"], i=ro["i"],
        theta=sep(-po, TNO),
        drot_cat=drot_cat, drot_sim=drot_sim,
        # Warsaw inbound-leg instrument (step_032): osc -> orig
        # direction discrepancy; the channel the Warsaw anomaly uses
        ddir_cat=sep(ps, po), ddir_sim=sep(ps, rb["phat"]),
        ddir_f_cat=sep(ps, pf), ddir_f_sim=sep(ps, rf["phat"]),
        d_back=sep(rb["phat"], po),
        d_fwd=sep(rf["phat"], pf),
        daa_cat=fo["aa"] - ro["aa"],
        daa_sim=rf["aa"] - rb["aa"],
        denc=min(rb["denc"], rf["denc"]),
        aa_back=rb["aa"], aa_fwd=rf["aa"], aa_osc=oo["aa"],
        t_back=rb["t_years"], t_fwd=rf["t_years"]))

logger.info(f"integrated {len(rows)} comets x2 legs; {len(failed)} failed")

# ------------------------------------------------------------------
# Statistics (identical structure to step_063)
# ------------------------------------------------------------------

def resid_cap(v, th, X_extra):
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
           "leg_offset_med":  float(np.median(leg)),
           "leg_offset_p95":  float(np.percentile(leg, 95)),
           "leg_offset_max":  float(leg.max()),
           "frac_legs_within_0p01": float(np.mean(leg < 0.01)),
           "frac_legs_within_0p1":  float(np.mean(leg < 0.1)),
           "med_drot_cat_in":  float(np.median(vc[inc])) if inc.any() else None,
           "med_drot_cat_out": float(np.median(vc[~inc])) if (~inc).any() else None,
           "med_drot_sim_in":  float(np.median(vs[inc])) if inc.any() else None,
           "med_drot_sim_out": float(np.median(vs[~inc])) if (~inc).any() else None,
           "med_abs_drot_diff": float(np.median(np.abs(vs - vc))),
           "frac_comets_drot_within_0p02": float(np.mean(np.abs(vs - vc) < 0.02))}
    rho, p = spearmanr(vs, vc)
    out["drot_sim_vs_cat"] = {"rho": float(rho), "p_2sided": float(p)}
    rho, p = spearmanr(daa_s, daa_c)
    out["energy_validation"] = {"rho": float(rho), "p_2sided": float(p),
        "med_abs_diff": float(np.median(np.abs(daa_s - daa_c)))}
    if inc.any() and (~inc).any():
        u = mannwhitneyu(vs[inc], vs[~inc], alternative="two-sided")
        out["drot_sim_in_vs_out_p"] = float(u.pvalue)
        u = mannwhitneyu(vc[inc], vc[~inc], alternative="two-sided")
        out["drot_cat_in_vs_out_p"] = float(u.pvalue)
    rho, p = spearmanr(th, vs)
    out["drot_sim_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
    rho, p = spearmanr(th, vc)
    out["drot_cat_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
    out["rot_per_kick_resid_kick"] = resid_cap(vs, th, [np.log10(K + 1.0)])
    out["rot_per_kick_resid_kick_denc"] = resid_cap(vs, th, [np.log10(K + 1.0), np.log10(D)])
    out["rot_per_kick_resid_full"] = resid_cap(vs, th,
        [np.log10(K + 1.0), np.log10(D), Q, I])
    # ---- Warsaw inbound-leg channel (the step_032 instrument) ----
    di_c = np.array([r["ddir_cat"] for r in sub])
    di_s = np.array([r["ddir_sim"] for r in sub])
    df_c = np.array([r["ddir_f_cat"] for r in sub])
    df_s = np.array([r["ddir_f_sim"] for r in sub])
    out["inbound_leg"] = {
        "med_ddir_cat_in":  float(np.median(di_c[inc])) if inc.any() else None,
        "med_ddir_cat_out": float(np.median(di_c[~inc])) if (~inc).any() else None,
        "med_ddir_sim_in":  float(np.median(di_s[inc])) if inc.any() else None,
        "med_ddir_sim_out": float(np.median(di_s[~inc])) if (~inc).any() else None,
    }
    if inc.any() and (~inc).any():
        u = mannwhitneyu(di_c[inc], di_c[~inc], alternative="two-sided")
        out["inbound_leg"]["ddir_cat_in_vs_out_p"] = float(u.pvalue)
        u = mannwhitneyu(di_s[inc], di_s[~inc], alternative="two-sided")
        out["inbound_leg"]["ddir_sim_in_vs_out_p"] = float(u.pvalue)
        u = mannwhitneyu(df_c[inc], df_c[~inc], alternative="two-sided")
        out["inbound_leg"]["outbound_cat_in_vs_out_p"] = float(u.pvalue)
        u = mannwhitneyu(df_s[inc], df_s[~inc], alternative="two-sided")
        out["inbound_leg"]["outbound_sim_in_vs_out_p"] = float(u.pvalue)
    rho, p = spearmanr(di_s, di_c)
    out["inbound_leg"]["ddir_sim_vs_cat"] = {"rho": float(rho), "p_2sided": float(p)}
    rho, p = spearmanr(th, di_c)
    out["inbound_leg"]["ddir_cat_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
    out["inbound_leg"]["resid_full"] = resid_cap(di_s, th,
        [np.log10(K + 1.0), np.log10(D), Q, I])
    return out

res = {
    "method": "REBOUND/IAS15 bidirectional on the Warsaw spike sample: "
              "osculating elements at catalogued perihelion epoch "
              "integrated to the +/-250 AU barycentric sphere under Sun + "
              "planetary-system barycentres (DE440s).  Boundary solutions "
              "barycentric, matching the catalogue convention.  "
              "Independent-catalogue replication of step_063 (CODE).",
    "n_spike": len(rows), "n_failed": len(failed), "failed": failed,
    "boundary_AU_barycentric": R_STOP,
    "ephemeris": "data/raw/spice/de440s.bsp (provenance pinned)",
    "matched": run([r for r in rows if r["q"] < 3.1]),
    "all_spike":  run(rows),
}

out = str(RESULTS / "step_b29_warsaw_bidirectional.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b29_warsaw_bidirectional.csv")
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
ax.set_title("Warsaw: independent integration vs catalogue", fontsize=10)

ax = axes[1]
ax.hist(np.log10(leg), bins=30, color="steelblue", alpha=0.8)
ax.axvline(np.log10(0.01), color="k", ls=":", lw=1)
ax.set_xlabel(r"per-leg boundary offset $\log_{10}$(deg)")
ax.set_ylabel("legs")
ax.set_title("sim vs catalogue, per boundary leg", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b29_warsaw_rotation_residual.png", dpi=150)

for tag, d in [("matched", res["matched"]), ("all_spike", res["all_spike"])]:
    logger.info(f"{tag}: n={d['n']} in-cap={d['n_in']} | leg offset med/p95 "
                f"{d['leg_offset_med']:.4f}/{d['leg_offset_p95']:.3f} deg | "
                f"sim-vs-cat rho={d['drot_sim_vs_cat']['rho']:.4f} | "
                f"rot-per-kick resid cap p={d['rot_per_kick_resid_full']['p']:.4f} "
                f"rho={d['rot_per_kick_resid_full']['rho']:.3f}")
print("wrote", out)
print("wrote", csv_out)
print("wrote", FIG / "step_b29_warsaw_rotation_residual.png")
