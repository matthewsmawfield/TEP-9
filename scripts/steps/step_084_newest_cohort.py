#!/usr/bin/env python3
"""Step 084 -- newest detached-cohort audit.

The SBDB snapshot now carries the 2025-2026 discovery cohort and the
fourth sednoid 2023 KQ14 ('Ammonite', Chen et al. 2025, Nat Astron 9,
1438), whose longitude of perihelion falls opposite the clustered
sednoids and has been widely reported as weakening the Planet Nine
case.  This step measures, on the real catalogue entries:

1. where Ammonite actually sits relative to the measured axis, and
   what removing or adding it does to the resident-sample statistics
   (a single-object jackknife on the cc <= 3 detached sample);
2. the newest detached cohort -- every a > 150, q > 30 object whose
   designation year is 2025 or later -- scored against the registered
   P1/P2 predictions of step_054: the in-cap fraction versus the
   footprint-only baseline (~0.42) and the mixture prediction (0.59),
   and the cohort mean varpi against the frozen 95 per cent interval
   [17, 88] and the footprint-only mean (~10.5 deg);
3. the same footprint model as steps 014/054 (designation half-month
   -> opposition longitude) applied to the new cohort itself, so the
   baseline is computed on the cohort's own discovery seasons.

The 2025-26 cohort orbits are provisional (median condition code ~7-8,
single- or few-opposition arcs); the audit reports them as a
diagnostic read, with the secure-orbit subset (cc <= 3) broken out
separately.  Ammonite is already inside the primary N = 44 sample via
its 2005 precovery, so its jackknife is retrospective -- but it is the
object whose addition was claimed to dissolve the clustering, which
makes the with/without comparison the informative quantity.

Outputs
-------
results/step_b49_newest_cohort.json
results/step_b49_newest_cohort.csv   (per-object scoring)
results/figures/supplementary/step_b49_newest_cohort.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_084_newest_cohort")
tee_stdout(logger)
logger.header("Newest detached-cohort audit")

import csv
import json
import math
import re
import numpy as np
from scipy.stats import binomtest, norm

rng = np.random.default_rng(20260918)
N_MC = 20000

AXIS = 49.0
CAP = 60.0

HALF = {"A": (1, 8), "B": (1, 23), "C": (2, 8), "D": (2, 22),
        "E": (3, 8), "F": (3, 23), "G": (4, 8), "H": (4, 23),
        "J": (5, 8), "K": (5, 23), "L": (6, 23), "M": (6, 23),
        "N": (7, 8), "O": (7, 23), "P": (8, 8), "Q": (8, 23),
        "R": (9, 8), "S": (9, 23), "T": (10, 8), "U": (10, 23),
        "V": (11, 8), "W": (11, 23), "X": (12, 8), "Y": (12, 23)}
DESIG = re.compile(r"\((\d{4})\s*([A-Z])([A-Z]?\d*)\)")

def sun_ecl_lon(year, month, day):
    from astropy.time import Time
    from astropy.coordinates import get_sun
    t = Time(f"{year:04d}-{month:02d}-{day:02d}T00:00:00",
             format="isot", scale="utc")
    return float(get_sun(t).geocentrictrueecliptic.lon.deg) % 360

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

def perih_dir(w, Om, inc):
    co, so, cO, sO, ci, si = np.cos(w), np.sin(w), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def circ_R(a):
    return abs(np.exp(1j * np.asarray(a)).mean())

TNO = lv(AXIS, -17.0)

# ------------------------------------------------------------------
# Load SBDB and build the detached table
# ------------------------------------------------------------------

d = json.loads((DATA_RAW / "sbdb" / "sbdb_outer_ss.json").read_text())
objs = []
for rec in d["data"]:
    o = dict(zip(d["fields"], rec))
    try:
        a, q = float(o["a"]), float(o["q"])
        om, w, i = float(o["om"]), float(o["w"]), float(o["i"])
        cc = int(o["condition_code"] or 9)
    except (TypeError, ValueError):
        continue
    if not (a > 150 and q > 30):
        continue
    name = str(o["full_name"]).strip()
    varpi = (om + w) % 360
    ph = perih_dir(math.radians(w), math.radians(om), math.radians(i))
    lam = np.nan; dyear = np.nan
    m = DESIG.search(name)
    if m and m.group(2) in HALF:
        yr, half = int(m.group(1)), m.group(2)
        mo, dy = HALF[half]
        lam = (sun_ecl_lon(yr, mo, dy) + 180) % 360
        dyear = yr
    objs.append(dict(name=name, a=a, q=q, cc=cc, om=om, w=w, i=i,
                     varpi=varpi, theta=sep(ph, TNO),
                     theta_aph=sep(-ph, TNO),
                     lam_opp=lam, dyear=dyear,
                     first_obs=str(o.get("first_obs") or "")))

logger.info(f"detached a>150 q>30: {len(objs)}")

def varpi_cap_frac(sub, center=AXIS):
    return float((np.abs((np.array([o["varpi"] for o in sub])
                          - center + 180) % 360 - 180) < CAP).mean())

def aph_cap_frac(sub):
    return float((np.array([o["theta"] for o in sub]) < CAP).mean())

def footprint_baseline(sub):
    """P1 baseline: each object's own opposition longitude smeared by
    the OSSOS-measured coupling width (50 deg)."""
    d_ax = np.abs((np.array([o["lam_opp"] for o in sub])
                   - AXIS + 180) % 360 - 180)
    ok = np.isfinite(d_ax)
    if ok.sum() == 0:
        return float("nan")
    return float(np.mean(norm.cdf((CAP - d_ax[ok]) / 50)
                         + norm.cdf((CAP + d_ax[ok]) / 50) - 1))

def conditioned_R(sub):
    """varpi Rayleigh R vs the step_010 conditioned null: hold (a,e,i),
    draw om,w uniform."""
    vp = np.deg2rad(np.array([o["varpi"] for o in sub]))
    R = circ_R(vp)
    null = np.empty(N_MC)
    for j in range(N_MC):
        null[j] = circ_R(rng.uniform(0, 2*np.pi, len(vp)))
    return R, float((int((null >= R).sum()) + 1) / (N_MC + 1))

def resultant3(sub):
    vecs = np.array([perih_dir(math.radians(o["w"]), math.radians(o["om"]),
                             math.radians(o["i"])) for o in sub])
    return float(np.linalg.norm(vecs.mean(axis=0)))

secure = [o for o in objs if o["cc"] <= 3]
newest = [o for o in objs if np.isfinite(o["dyear"]) and o["dyear"] >= 2025]
ammon  = [o for o in objs if "KQ14" in o["name"]]

res = {"method":"newest-cohort audit on the live SBDB snapshot: the "
       "2025+ designation cohort scored against the step_054 registered "
       "predictions on its own designation-derived footprint, plus the "
       "single-object jackknife of 2023 KQ14 inside the secure sample.",
       "axis_deg":AXIS, "axis_beta":-17.0, "cap_deg":CAP,
       "n_detached_all_cc":len(objs), "n_secure_cc3":len(secure),
       "seed":20260918}

# ------------------------------------------------------------------
# 1. Ammonite: position and jackknife
# ------------------------------------------------------------------

if ammon:
    a0 = ammon[0]
    with_a = [o for o in secure]
    without_a = [o for o in secure if "KQ14" not in o["name"]]
    Rw, pw = conditioned_R(with_a)
    Ro, po = conditioned_R(without_a)
    res["ammonite"] = {
        "name": a0["name"], "a": a0["a"], "q": a0["q"], "i": a0["i"],
        "varpi": a0["varpi"], "theta_to_axis": a0["theta"],
        "theta_aphelion_to_axis": a0["theta_aph"],
        "in_cap": a0["theta"] < CAP, "cc": a0["cc"],
        "first_obs": a0["first_obs"],
        "note": "Chen et al. 2025, Nat Astron (FOSSIL II): varpi "
                "anti-aligned with the clustered sednoids; already "
                "inside the secure sample via 2005 precovery",
        "R_varpi_with": Rw, "p_with": pw,
        "R_varpi_without": Ro, "p_without": po,
        "R3_with": resultant3(with_a), "R3_without": resultant3(without_a),
        "n_with": len(with_a), "n_without": len(without_a)}
    logger.info(f"Ammonite: varpi={a0['varpi']:.1f}, theta={a0['theta']:.1f} deg "
                f"({'in' if a0['theta']<CAP else 'out'}-cap); "
                f"R {Ro:.3f}->{Rw:.3f} (p {po:.4f}->{pw:.4f})")

# ------------------------------------------------------------------
# 2. Newest cohort vs registered predictions
# ------------------------------------------------------------------

newest.sort(key=lambda o: o["name"])
n = len(newest)
n_cap_varpi = int(sum(abs((o["varpi"] - AXIS + 180) % 360 - 180) < CAP
                      for o in newest))
n_cap_aph = int(sum(o["theta"] < CAP for o in newest))
frac_v = n_cap_varpi / n if n else float("nan")
frac_a = n_cap_aph / n if n else float("nan")
base = footprint_baseline(newest)
lam = np.array([o["lam_opp"] for o in newest])
okl = np.isfinite(lam)
lam_R = circ_R(np.deg2rad(lam[okl])) if okl.sum() else float("nan")
lam_mean = float(np.rad2deg(np.angle(np.exp(1j*np.deg2rad(lam[okl])).mean())) % 360) \
    if okl.sum() else float("nan")
mv = float(np.rad2deg(np.angle(np.exp(1j*np.deg2rad(
    np.array([o["varpi"] for o in newest]))).mean())) % 360) if n else float("nan")
ccs = [o["cc"] for o in newest]

res["newest_cohort"] = {
    "definition": "a>150 & q>30, designation year >= 2025, any condition code",
    "n": n,
    "cc_median": float(np.median(ccs)) if n else float("nan"),
    "cc_range": [int(min(ccs)), int(max(ccs))] if n else None,
    "n_cc3": int(sum(1 for c in ccs if c <= 3)),
    "varpi_cap": {"n_in": n_cap_varpi, "frac": frac_v,
                  "p_vs_footprint_baseline": float(
                      binomtest(n_cap_varpi, n, base,
                                alternative="greater").pvalue)
                  if n and np.isfinite(base) else float("nan"),
                  "p_vs_prediction_0p59": float(
                      binomtest(n_cap_varpi, n, 0.59).pvalue) if n else float("nan")},
    "peri_cap_3d": {"n_in": n_cap_aph, "frac": frac_a,
                    "note": "perihelion direction within 60 deg of the "
                            "axis -- the resident-signature convention"},
    "footprint": {"baseline_frac": base,
                  "lam_opp_R": float(lam_R),
                  "lam_opp_mean_deg": lam_mean,
                  "n_with_lam": int(okl.sum())},
    "mean_varpi_deg": mv,
    "P2_check": {"predicted_mean": AXIS, "ci_95": [17.0, 88.0],
                 "footprint_mean_baseline": 10.5,
                 "cohort_mean_in_CI": bool(
                     min(abs((mv - AXIS + 180) % 360 - 180),
                         abs((mv - AXIS + 180) % 360 - 180)) < 71.0) if n else None},
    "caveat": "2025-26 designations carry provisional orbits (median "
              "cc ~ 7-8); a diagnostic read, not a confirmation"}

logger.info(f"newest cohort n={n}: varpi-cap {n_cap_varpi}/{n}={frac_v:.2f} "
            f"vs baseline {base:.2f}; aph-cap {n_cap_aph}/{n}={frac_a:.2f}; "
            f"mean varpi {mv:.1f}; footprint R {lam_R:.2f} mean {lam_mean:.0f}")

# secure-orbit subset of the newest cohort (if any)
sec_new = [o for o in newest if o["cc"] <= 3]
res["newest_secure"] = {"n": len(sec_new),
    "objects": [o["name"] for o in sec_new]}

# ------------------------------------------------------------------
# 3. Footprint of the newest cohort alone
# ------------------------------------------------------------------

res["note"] = ("The cohort's own footprint baseline (its discovery "
               "opposition longitudes smeared by the measured coupling) "
               "is the honest comparison set by step_054; the mixture "
               "prediction 0.59 is the registered TEP expectation; "
               "uniform would give 0.33.")

# ------------------------------------------------------------------
# Write
# ------------------------------------------------------------------

with open(RESULTS/"step_b49_newest_cohort.json","w") as f:
    json.dump(res, f, indent=1, default=float)
with open(RESULTS/"step_b49_newest_cohort.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["name","a","q","i","cc","varpi",
                                      "theta","theta_aph","lam_opp","dyear","first_obs"])
    w.writeheader()
    for o in newest + ammon:
        w.writerow({k: o[k] for k in w.fieldnames})
logger.data_save(RESULTS/'step_b49_newest_cohort.json')
logger.data_save(RESULTS/'step_b49_newest_cohort.csv')
# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

ax = axes[0]
vp_all = [o["varpi"] for o in secure]
ax.hist(vp_all, bins=np.arange(0, 361, 20), color="0.75",
        label=f"secure detached (n={len(secure)})")
vp_new = [o["varpi"] for o in newest]
ax.hist(vp_new, bins=np.arange(0, 361, 20), color="crimson", alpha=0.65,
        label=f"2025+ cohort (n={n})")
if ammon:
    ax.axvline(ammon[0]["varpi"], color="darkred", lw=1.8,
               label=f"2023 KQ14 ({ammon[0]['varpi']:.0f} deg)")
ax.axvspan(AXIS-CAP, AXIS+CAP, color="crimson", alpha=0.08)
ax.axvline(AXIS, color="k", ls="--", lw=1, label="axis 49 deg")
ax.set_xlabel("longitude of perihelion $\\varpi$ (deg)")
ax.set_ylabel("count"); ax.legend(frameon=False, fontsize=8)
ax.set_title("newest cohort vs the cluster", fontsize=10)

ax = axes[1]
names = ["footprint\nbaseline", "uniform", "observed\n2025+ cohort",
         "P1\nprediction"]
vals = [base, 2*CAP/360, frac_v, 0.59]
ax.bar(names, vals, color=["steelblue","0.6","crimson","0.35"])
for i, v in enumerate(vals):
    ax.text(i, v+0.01, f"{v:.2f}", ha="center", fontsize=9)
ax.set_ylabel("fraction within 60 deg of axis")
ax.set_ylim(0, 1)
ax.set_title("P1 score on the newest cohort", fontsize=10)

ax = axes[2]
th_new = [o["theta"] for o in newest]
th_sec = [o["theta"] for o in secure]
ax.hist(th_sec, bins=np.arange(0, 181, 15), color="0.75",
        label="secure detached")
ax.hist(th_new, bins=np.arange(0, 181, 15), color="crimson", alpha=0.65,
        label="2025+ cohort")
ax.axvline(CAP, color="k", ls="--", lw=1, label="cap 60 deg")
if ammon:
    ax.axvline(ammon[0]["theta"], color="darkred", lw=1.8,
               label=f"KQ14 ({ammon[0]['theta']:.0f} deg)")
ax.set_xlabel("perihelion angle to axis $\\theta$ (deg)")
ax.set_ylabel("count"); ax.legend(frameon=False, fontsize=8)
ax.set_title("3-D perihelion direction (resident channel)", fontsize=10)

fig.tight_layout()
FIG = RESULTS/"figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG/"supplementary" / "step_b49_newest_cohort.png", dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b49_newest_cohort.png')