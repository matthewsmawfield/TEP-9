#!/usr/bin/env python3
"""Step 085 -- DES independent resident-cohort audit.

Every resident-signature test so far runs on the JPL SBDB catalogue:
one database, one orbit-determination pipeline, discoveries drawn
overwhelmingly from northern surveys.  The Dark Energy Survey TNO
release (Bernardinelli et al. 2022, ApJS 258, 41; data release
bernardinelli/des_tno_catalog) is the largest single-survey TNO
catalogue with uniform selection -- 814 objects discovered in a
contiguous 5000 deg^2 of southern sky with orbits fitted
independently of SBDB.  Its detached subset is therefore the only
public independent resident cohort on which the registered axis can
be scored.

This step:

1. downloads the DES orbit table to data/raw/des/ with the same
   provenance record used by steps 001-004;
2. selects the detached cohort (a > 150, q > 30) -- the 16 objects
   Bernardinelli et al. themselves reported as "consistent with the
   null hypothesis of azimuthal isotropy" (their null is the 2-D
   azimuthal distribution; the pre-declared directional cap is a
   different, more powerful statistic);
3. scores the cohort against the registered P1 axis: in-cap fraction
   versus uniform (1/3), versus the step_054 footprint baseline
   (~0.42), versus the registered mixture prediction (0.59), and
   versus the cohort's own designation-derived footprint baseline
   (the same half-month opposition-longitude model as steps 014/054/
   084) -- the honest comparison, since DES opposition fields sit in
   the axis sector;
4. cross-matches the cohort against the SBDB secure sample to (a)
   measure the sample overlap -- how independent the DES cohort
   actually is in membership -- and (b) measure the resident-side
   inter-catalogue floor: the per-object |Delta varpi| between the
   DES and SBDB orbit determinations of the same bodies, the
   resident analogue of the step_079 comet cross-fitter audit.

Outputs
-------
results/step_b50_des_resident.json
results/step_b50_des_resident.csv   (per-object scoring)
results/figures/step_b50_des_resident.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_085_des_resident")
tee_stdout(logger)
logger.header("DES independent resident-cohort audit")

import csv
import hashlib
import json
import math
import re
import urllib.request
from datetime import datetime, timezone
import numpy as np
from scipy.stats import binomtest, norm, rayleigh
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import get_sun

rng = np.random.default_rng(20260918)
N_MC = 20000

AXIS, BETA, CAP = 49.0, -17.0, 60.0

# ------------------------------------------------------------------
# 1. Download DES orbit catalogue with provenance
# ------------------------------------------------------------------

DES_URL = ("https://raw.githubusercontent.com/bernardinelli/"
           "des_tno_catalog/main/y6_des_tnos_color.fits")
DES_DIR = DATA_RAW / "des"
DES_DIR.mkdir(exist_ok=True)
FIT = DES_DIR / "y6_des_tnos_color.fits"

req = urllib.request.Request(DES_URL, headers={"User-Agent": "TEP9-pipeline"})
with urllib.request.urlopen(req, timeout=120) as r:
    blob = r.read()
FIT.write_bytes(blob)
sha = hashlib.sha256(blob).hexdigest()
prov = {"step": "step_085_des_resident",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "files": {"y6_des_tnos_color.fits": {
            "url": DES_URL,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "bytes": len(blob), "sha256": sha,
            "reference": "Bernardinelli et al. 2022, ApJS, 258, 41 "
                         "(DES six-year TNO catalogue, 814 objects)"}}}
json.dump(prov, open(DES_DIR / "provenance.json", "w"), indent=1)
logger.info(f"DES catalogue: {len(blob)} bytes, sha256 {sha[:16]}...")

tab = fits.open(str(FIT))[1].data
logger.info(f"DES objects: {len(tab)}")

# ------------------------------------------------------------------
# 2. Helpers (same conventions as steps 010/054/084)
# ------------------------------------------------------------------

HALF = {"A": (1, 8), "B": (1, 23), "C": (2, 8), "D": (2, 22),
        "E": (3, 8), "F": (3, 23), "G": (4, 8), "H": (4, 23),
        "J": (5, 8), "K": (5, 23), "L": (6, 8), "M": (6, 23),
        "N": (7, 8), "O": (7, 23), "P": (8, 8), "Q": (8, 23),
        "R": (9, 8), "S": (9, 23), "T": (10, 8), "U": (10, 23),
        "V": (11, 8), "W": (11, 23), "X": (12, 8), "Y": (12, 23)}
DESIG = re.compile(r"^(\d{4})\s*([A-Z])([A-Z]?\d*)")

def sun_ecl_lon(year, month, day):
    t = Time(f"{year:04d}-{month:02d}-{day:02d}T00:00:00",
             format="isot", scale="utc")
    return float(get_sun(t).geocentrictrueecliptic.lon.deg) % 360

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l), math.cos(b) * math.sin(l),
                     math.sin(b)])

def perih_dir(w, Om, inc):
    co, so, cO, sO, ci, si = (np.cos(w), np.sin(w), np.cos(Om),
                             np.sin(Om), np.cos(inc), np.sin(inc))
    return np.array([cO * co - sO * so * ci, sO * co + cO * so * ci, so * si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def circ_R(a):
    return abs(np.exp(1j * np.asarray(a)).mean())

def footprint_baseline(lams):
    d_ax = np.abs((np.asarray(lams) - AXIS + 180) % 360 - 180)
    ok = np.isfinite(d_ax)
    if ok.sum() == 0:
        return float("nan")
    return float(np.mean(norm.cdf((CAP - d_ax[ok]) / 50)
                         + norm.cdf((CAP + d_ax[ok]) / 50) - 1))

TNO = lv(AXIS, BETA)

def desig_lam(name):
    m = DESIG.search(str(name).strip())
    if not m or m.group(2) not in HALF:
        return float("nan")
    mo, dy = HALF[m.group(2)]
    return float((sun_ecl_lon(int(m.group(1)), mo, dy) + 180) % 360)

# ------------------------------------------------------------------
# 3. DES detached cohort
# ------------------------------------------------------------------

sel = (tab["a"] > 150) & (tab["q"] > 30)
d_objs = []
for r in tab[sel]:
    varpi = float((r["lan"] + r["aop"]) % 360)
    ph = perih_dir(math.radians(float(r["aop"])),
                   math.radians(float(r["lan"])),
                   math.radians(float(r["i"])))
    name = str(r["MPC"]).strip()
    d_objs.append(dict(name=name, a=float(r["a"]), q=float(r["q"]),
                       e=float(r["e"]), i=float(r["i"]),
                       varpi=varpi, theta=sep(ph, TNO),
                       cls=str(r["Class"]).strip(),
                       lam_opp=desig_lam(name)))

n = len(d_objs)
n_cap_v = int(sum(abs((o["varpi"] - AXIS + 180) % 360 - 180) < CAP
                  for o in d_objs))
n_cap_3d = int(sum(o["theta"] < CAP for o in d_objs))
base = footprint_baseline([o["lam_opp"] for o in d_objs])
vps = np.array([o["varpi"] for o in d_objs])
R_varpi = circ_R(np.deg2rad(vps))
mean_varpi = float(np.rad2deg(np.angle(
    np.exp(1j * np.deg2rad(vps)).mean())) % 360)

# Rayleigh p for the cohort resultant
z = n * R_varpi ** 2
p_ray = float(rayleigh.sf(np.sqrt(2 * z)))

res = {"method": "DES (Bernardinelli et al. 2022, ApJS 258, 41) "
       "detached cohort scored against the pre-declared axis; the "
       "cohort's own designation-derived footprint baseline is the "
       "honest comparison -- DES opposition fields sit in the axis "
       "sector, so an elevated in-cap rate is partly geometric.",
       "axis_deg": AXIS, "axis_beta": BETA, "cap_deg": CAP,
       "catalogue": {"n_total": int(len(tab)), "detached_n": n,
                     "file": "data/raw/des/y6_des_tnos_color.fits",
                     "sha256": sha,
                     "survey_area_deg2": 5000},
       "detached": {
           "n": n,
           "varpi_cap": {"n_in": n_cap_v, "frac": n_cap_v / n,
                         "p_vs_uniform_1_3": float(
                             binomtest(n_cap_v, n, 1 / 3).pvalue),
                         "p_vs_baseline_0p42": float(
                             binomtest(n_cap_v, n, 0.42).pvalue),
                         "p_vs_prediction_0p59": float(
                             binomtest(n_cap_v, n, 0.59).pvalue),
                         "p_vs_own_footprint": float(
                             binomtest(n_cap_v, n, base).pvalue)
                         if np.isfinite(base) else float("nan")},
           "peri_cap_3d": {"n_in": n_cap_3d, "frac": n_cap_3d / n},
           "footprint": {"baseline_frac": base,
                         "lam_opp": [o["lam_opp"] for o in d_objs]},
           "R_varpi": float(R_varpi), "mean_varpi_deg": mean_varpi,
           "p_rayleigh": p_ray,
           "note": "Bernardinelli et al. report azimuthal isotropy "
                   "for this cohort; the pre-declared directional cap "
                   "is a different statistic -- membership against a "
                   "fixed registered direction, not uniformity of the "
                   "2-D distribution."}}

# ------------------------------------------------------------------
# 4. Overlap with the SBDB secure sample + resident cross-fitter floor
# ------------------------------------------------------------------

d = json.loads((DATA_RAW / "sbdb" / "sbdb_outer_ss.json").read_text())
sb = {}
for rec in d["data"]:
    o = dict(zip(d["fields"], rec))
    try:
        a, q = float(o["a"]), float(o["q"])
        om, w = float(o["om"]), float(o["w"])
        cc = int(o["condition_code"] or 9)
    except (TypeError, ValueError):
        continue
    if a > 150 and q > 30 and cc <= 3:
        nm = str(o["full_name"]).strip()
        m = re.search(r"\((\d{4}\s*[A-Z]+\d*)\)", nm)
        key = m.group(1).replace(" ", "") if m else nm
        sb[key] = dict(name=nm, varpi=(om + w) % 360, cc=cc)

matches = []
for o in d_objs:
    key = o["name"].replace(" ", "")
    if key in sb:
        dv = abs((o["varpi"] - sb[key]["varpi"] + 180) % 360 - 180)
        o["sbdb_match"] = sb[key]["name"]
        o["sbdb_varpi"] = sb[key]["varpi"]
        o["dvarpi_des_sbdb"] = float(dv)
        matches.append(o)
    else:
        o["sbdb_match"] = ""
        o["sbdb_varpi"] = float("nan")
        o["dvarpi_des_sbdb"] = float("nan")

dvv = np.array([o["dvarpi_des_sbdb"] for o in matches])
cap_des = [o["name"] for o in d_objs if abs((o["varpi"] - AXIS + 180) % 360 - 180) < CAP]
cap_sb = [o["name"] for o in d_objs
          if np.isfinite(o["sbdb_varpi"])
          and abs((o["sbdb_varpi"] - AXIS + 180) % 360 - 180) < CAP]

res["sbdb_overlap"] = {
    "n_matched": len(matches),
    "matched_names": [o["name"] for o in matches],
    "n_des_only": n - len(matches),
    "des_only_names": [o["name"] for o in d_objs if not o["sbdb_match"]],
    "dvarpi_median_deg": float(np.median(dvv)) if len(dvv) else float("nan"),
    "dvarpi_max_deg": float(np.max(dvv)) if len(dvv) else float("nan"),
    "cap_membership_agreement": sorted(set(cap_des) & set(cap_sb)),
    "note": "per-object |dvarpi| between DES and SBDB orbit "
            "determinations -- the resident-side inter-catalogue "
            "floor, analogous to the step_079 comet fitter audit"}

logger.info(f"DES detached n={n}: varpi-cap {n_cap_v}/{n} "
            f"({n_cap_v / n:.2f}) vs uniform 0.33, baseline "
            f"{base:.2f}, registered 0.59")
logger.info(f"SBDB overlap {len(matches)}/{n}; dvarpi floor "
            f"{res['sbdb_overlap']['dvarpi_median_deg']:.2f} deg")

# ------------------------------------------------------------------
# 5. Write
# ------------------------------------------------------------------

out = str(RESULTS / "step_b50_des_resident.json")
json.dump(res, open(out, "w"), indent=1, default=float)
csv_out = str(RESULTS / "step_b50_des_resident.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(d_objs[0].keys()))
    w.writeheader(); w.writerows(d_objs)
print("wrote", out)
print("wrote", csv_out)

# ------------------------------------------------------------------
# 6. Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

ax = axes[0]
dsec = json.loads((DATA_RAW / "sbdb" / "sbdb_outer_ss.json").read_text())
vp_all = []
for rec in dsec["data"]:
    o = dict(zip(dsec["fields"], rec))
    try:
        if float(o["a"]) > 150 and float(o["q"]) > 30 \
                and int(o["condition_code"] or 9) <= 3:
            vp_all.append((float(o["om"]) + float(o["w"])) % 360)
    except (TypeError, ValueError):
        continue
ax.hist(vp_all, bins=np.arange(0, 361, 20), color="0.75",
        label=f"SBDB secure (n={len(vp_all)})")
ax.hist(vps, bins=np.arange(0, 361, 20), color="teal", alpha=0.7,
        label=f"DES detached (n={n})")
ax.axvspan(AXIS - CAP, AXIS + CAP, color="crimson", alpha=0.08)
ax.axvline(AXIS, color="k", ls="--", lw=1, label="axis 49 deg")
ax.set_xlabel("longitude of perihelion $\\varpi$ (deg)")
ax.set_ylabel("count"); ax.legend(frameon=False, fontsize=8)
ax.set_title("DES cohort vs the registered axis", fontsize=10)

ax = axes[1]
names = ["uniform", "SBDB-sample\nbaseline", "own footprint\nbaseline",
         "registered\nP1", "observed\nDES"]
vals = [1 / 3, 0.42, base, 0.59, n_cap_v / n]
ax.bar(names, vals, color=["0.6", "steelblue", "teal", "0.35", "crimson"])
for i, v in enumerate(vals):
    ax.text(i, v + 0.015, f"{v:.2f}", ha="center", fontsize=9)
ax.set_ylabel("fraction within 60 deg of axis")
ax.set_ylim(0, 1)
ax.set_title("P1 score on the DES cohort", fontsize=10)

ax = axes[2]
lam = np.array([o["lam_opp"] for o in d_objs])
ok = np.isfinite(lam)
ax.scatter(np.deg2rad(vps[ok]), np.deg2rad(lam[ok]), s=40,
           c="teal", zorder=3)
for o in d_objs:
    if np.isfinite(o["lam_opp"]):
        ax.annotate(o["name"].split()[-1],
                    (np.deg2rad(o["varpi"]), np.deg2rad(o["lam_opp"])),
                    fontsize=6, alpha=0.7)
ax.axvline(np.deg2rad(AXIS), color="k", ls="--", lw=1)
ax.set_xlim(0, 2 * np.pi); ax.set_ylim(0, 2 * np.pi)
ax.set_xticks(np.deg2rad([0, 90, 180, 270, 360]))
ax.set_xticklabels(["0", "90", "180", "270", "360"])
ax.set_yticks(np.deg2rad([0, 90, 180, 270, 360]))
ax.set_yticklabels(["0", "90", "180", "270", "360"])
ax.set_xlabel("$\\varpi$ (deg)"); ax.set_ylabel("opposition longitude (deg)")
ax.set_title("perihelion vs discovery longitude", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b50_des_resident.png", dpi=150)
print(f"wrote {FIG / 'step_b50_des_resident.png'}")
