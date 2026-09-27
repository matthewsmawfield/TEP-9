"""step_152: discovery-geometry-conditioned null for the arrival anisotropy (b116).

Step 149 (b113) found that the class-1 arrival directions themselves
over-concentrate inside the declared 60-deg cap about the resident
axis: 48/131 = 36.6% against the 25% uniform-sky fraction (binomial
p = 0.0020; the 43/131 = 32.8% figure is the comet-axis cap, kept in
b113 T6 for both axes).  That null assumes isotropic arrivals.  A referee's first objection is that cometary
perihelion directions are discovery-biased toward the ecliptic band, so
a cap sitting near the ecliptic collects excess arrivals by
construction.  This step replaces the uniform-sky null with the
conditioned nulls that objection requires.

Tests
  T1  ecliptic-latitude-conditioned null.  Each arrival's ecliptic
      latitude is kept fixed and its longitude is redrawn uniformly
      (N_DRAW draws over the whole sample).  This is the exact null of
      "arrivals uniform in longitude given their observed latitude
      distribution" -- the first-order discovery-pointing model.  The
      conditioned p is the draw fraction with in-cap fraction >=
      observed.  Also reports the null median in-cap fraction, which
      prices how much of the excess the ecliptic band alone explains.
  T2  Galactic-latitude-conditioned null.  Same construction in the
      Galactic frame (the axis/cap converted to Galactic coordinates),
      controlling for concentration toward the Galactic plane rather
      than the ecliptic.
  T3  axis-longitude specificity.  The cap is held at the axis's own
      ecliptic latitude and scanned through longitude: the in-cap
      fraction is recomputed for every 5-deg longitude placement.
      The observed placement's rank among the 72 positions tests
      whether the excess is specific to the axis longitude or generic
      to the axis's latitude band.

Subsets: matched (q < 3.1 AU) and all class-1, mirroring step 149.

Inputs
  results/step_b28_bidirectional_rotation.csv
  data/raw/code/code_original.html

Outputs
  results/step_b116_arrival_latitude_null.json
  results/step_b116_arrival_latitude_null.csv
  results/figures/supplementary/step_b116_arrival_latitude_null.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (
    DATA_RAW, RESULTS, lv, lb, sep, perih_dir, parse_code, tee_stdout)
logger = StepLogger("step_152_arrival_latitude_null")
tee_stdout(logger)
logger.header("Arrival-anisotropy conditioned nulls -- "
              "ecliptic/galactic latitude shuffles and axis-longitude scan")

import csv
import json
import math
import numpy as np
from scipy.stats import binomtest

SEED = 20261030
CAP = 60.0
N_DRAW = 50000
LON_STEP = 5.0
rng = np.random.default_rng(SEED)

TNO = lv(49.0, -17.0)          # resident cluster axis
CAX = lv(34.0, -13.0)          # declared comet transit axis

# ecliptic -> galactic rotation, astropy-free, via J2000 basis vectors
_EPS = math.radians(23.4392911)

def _radec_vec(ra, dec):
    ra, dec = math.radians(ra), math.radians(dec)
    return np.array([math.cos(ra) * math.cos(dec),
                     math.sin(ra) * math.cos(dec), math.sin(dec)])

_GX = _radec_vec(266.405, -28.936)          # galactic centre (l = 0)
_GZ = _radec_vec(192.859508, 27.128336)     # north galactic pole
_GY = np.cross(_GZ, _GX); _GY /= np.linalg.norm(_GY)
_GX = np.cross(_GY, _GZ)
_EQ2GAL = np.stack([_GX, _GY, _GZ])         # equatorial -> galactic

def _ecl2eq(v):
    return np.array([v[0],
                     math.cos(_EPS) * v[1] - math.sin(_EPS) * v[2],
                     math.sin(_EPS) * v[1] + math.cos(_EPS) * v[2]])

def ecl2gal(v):
    """Ecliptic unit vector -> galactic unit vector (J2000)."""
    return _EQ2GAL @ _ecl2eq(v)

TNO_G = ecl2gal(TNO)
CAX_G = ecl2gal(CAX)

# ------------------------------------------------------------------
# data -- identical join to step_149
# ------------------------------------------------------------------
orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))
rows = []
for r in csv.DictReader(open(RESULTS / "step_b28_bidirectional_rotation.csv")):
    k = r["desig"]
    if k not in orig:
        continue
    ro = orig[k]
    u = -perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                   math.radians(ro["i"]))
    rows.append(dict(desig=k, q=float(r["q"]), u=u, u_g=ecl2gal(u)))
logger.info(f"joined {len(rows)} comets (step_b28 x CODE original)")

subsets = {"matched": [r for r in rows if r["q"] < 3.1],
           "all_c1": rows}


def lonlat_rot(u, dlam):
    """Rotate unit vector u about the z axis by dlam radians."""
    c, s = math.cos(dlam), math.sin(dlam)
    return np.array([c * u[0] - s * u[1], s * u[0] + c * u[1], u[2]])


def conditioned_null(U, ax, n_draw):
    """Fraction inside the cap when each direction's longitude is
    shuffled independently at fixed latitude (per-draw)."""
    n = len(U)
    cosc = math.cos(math.radians(CAP))
    f_obs = float((U @ ax >= cosc).mean())
    n_in = int((U @ ax >= cosc).sum())
    sky_frac = float((1 - cosc) / 2)
    z = ax
    draws = np.empty(n_draw)
    for i in range(n_draw):
        dl = rng.uniform(0, 2 * math.pi, n)
        c, s = np.cos(dl), np.sin(dl)
        Ur = np.stack([c * U[:, 0] - s * U[:, 1],
                       s * U[:, 0] + c * U[:, 1], U[:, 2]], axis=1)
        draws[i] = (Ur @ z >= cosc).mean()
    p = float((1 + (draws >= f_obs).sum()) / (n_draw + 1))
    return dict(n=n, n_in=n_in, f_obs=f_obs,
                cap_sky_frac=sky_frac,
                p_binom=float(binomtest(n_in, n, sky_frac,
                                        alternative="greater").pvalue),
                null_med=float(np.median(draws)),
                null_p95=float(np.quantile(draws, 0.95)),
                null_p99=float(np.quantile(draws, 0.99)),
                p_cond=p)


def axis_lon_scan(U, ax):
    """Cap at the axis's own latitude, scanned through longitude."""
    lam0, b0 = lb(ax)
    fr = []
    for lam in np.arange(0.0, 360.0, LON_STEP):
        z = lv(lam, b0)
        fr.append(float((U @ z >= math.cos(math.radians(CAP))).mean()))
    fr = np.array(fr)
    f_obs = float((U @ ax >= math.cos(math.radians(CAP))).mean())
    rank = int((fr >= f_obs).sum())
    return dict(cap_lat=b0, f_obs=f_obs,
                scan_med=float(np.median(fr)), scan_max=float(fr.max()),
                scan_argmax=float(np.arange(0.0, 360.0, LON_STEP)[fr.argmax()]),
                rank_ge=int(rank), n_pos=len(fr),
                frac_rank=rank / len(fr))


out = {"step": "step_152_arrival_latitude_null", "seed": SEED,
       "n_draw": N_DRAW, "results": {}}
csvs = []
for stag, sub in subsets.items():
    U = np.array([r["u"] for r in sub])
    UG = np.array([r["u_g"] for r in sub])
    res = {}
    for aname, ax, axg in [("tno", TNO, TNO_G), ("cax", CAX, CAX_G)]:
        t1 = conditioned_null(U, ax, N_DRAW)
        t2 = conditioned_null(UG, axg, N_DRAW)
        t3 = axis_lon_scan(U, ax)
        res[f"T1_ecl_lat_{aname}"] = t1
        res[f"T2_gal_lat_{aname}"] = t2
        res[f"T3_axis_lon_{aname}"] = t3
        logger.info(f"{stag}/{aname}: in-cap {t1['f_obs']:.3f} "
                    f"(ecl-lat null med {t1['null_med']:.3f}, "
                    f"p_cond {t1['p_cond']:.4f}; gal-lat p_cond "
                    f"{t2['p_cond']:.4f}; lon rank {t3['rank_ge']}/"
                    f"{t3['n_pos']})")
        csvs.append(dict(subset=stag, axis=aname, **{
            f"ecl_{k}": v for k, v in t1.items()}, **{
            f"gal_{k}": v for k, v in t2.items()}, **{
            f"lon_{k}": v for k, v in t3.items()}))
    out["results"][stag] = res

out["verdict"] = ("arrival anisotropy survives ecliptic-band conditioning"
                  if out["results"]["all_c1"]["T1_ecl_lat_tno"]["p_cond"] < 0.05
                  else "arrival anisotropy partly attributable to band "
                       "concentration -- see conditioned nulls")
out["caveats"] = [
    "longitude shuffle preserves latitude exactly; the second-order "
    "bias (opposition-pointing longitude structure) is not modeled -- "
    "it requires per-comet discovery-epoch opposition geometry",
    "Galactic-frame test controls for plane-of-galaxy concentration "
    "only; correlated longitude structure is preserved by neither null",
]

with open(RESULTS / "step_b116_arrival_latitude_null.json", "w") as f:
    json.dump(out, f, indent=1)
logger.add_output_file(RESULTS / "step_b116_arrival_latitude_null.json")
with open(RESULTS / "step_b116_arrival_latitude_null.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=sorted({k for x in csvs for k in x}))
    w.writeheader()
    for x in csvs:
        w.writerow(x)
logger.add_output_file(RESULTS / "step_b116_arrival_latitude_null.csv")

# figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
sub = subsets["all_c1"]
U = np.array([r["u"] for r in sub])
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
for i, (aname, axv) in enumerate([("tno", TNO), ("cax", CAX)]):
    n = len(U)
    cosc = math.cos(math.radians(CAP))
    draws = np.empty(20000)
    for j in range(20000):
        dl = rng.uniform(0, 2 * math.pi, n)
        c, s = np.cos(dl), np.sin(dl)
        Ur = np.stack([c * U[:, 0] - s * U[:, 1],
                       s * U[:, 0] + c * U[:, 1], U[:, 2]], axis=1)
        draws[j] = (Ur @ axv >= cosc).mean()
    f_obs = (U @ axv >= cosc).mean()
    ax[i].hist(draws, bins=60, color="0.7")
    ax[i].axvline(f_obs, color="r", lw=1.5)
    ax[i].set_xlabel("in-cap fraction (ecliptic-latitude null)")
    ax[i].set_title(f"{aname} axis, all class-1 (n={n})")
t3 = out["results"]["all_c1"]["T3_axis_lon_tno"]
lams = np.arange(0.0, 360.0, LON_STEP)
fr = [(U @ lv(la, t3["cap_lat"]) >= math.cos(math.radians(CAP))).mean()
      for la in lams]
ax[2].plot(lams, fr, lw=1)
ax[2].axhline(t3["f_obs"], color="r", lw=1.5)
ax[2].set_xlabel("cap longitude at axis latitude (deg)")
ax[2].set_ylabel("in-cap fraction")
ax[2].set_title("axis-longitude specificity (TNO cap)")
fig.tight_layout()
fig.savefig(FIG / "supplementary" / "step_b116_arrival_latitude_null.png", dpi=300)
logger.add_output_file(FIG / "supplementary" / "step_b116_arrival_latitude_null.png")
logger.info("verdict: " + out["verdict"])
print(json.dumps(out["results"], indent=1)[:1500])
print(f"VERDICT: {out['verdict']}")
