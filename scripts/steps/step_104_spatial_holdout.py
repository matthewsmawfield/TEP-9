#!/usr/bin/env python3
"""Step 104 -- spatial holdout and axis localization of the slip map.

Step 102 validated the bipolar slip map under random five-fold
partitions and cross-catalogue transfer.  Random folds, however,
leave a spatially neighbouring twin of every held-out comet inside
the training set, so they cannot distinguish a genuine sky field
from a pattern carried by a few clustered objects.  This step
closes that loophole and turns the validation into a localization
measurement:

  T1  azimuthal spatial holdout (pooled, n = 229).  Per-comet
      aphelion unit vectors are reconstructed from the catalogues'
      original-orbit elements (the same construction step 063 uses
      for theta) and the sky is split in half by a great circle
      through the axis; the map dtau = a + b*cos(2*theta) is fitted
      on one half and used to predict the signed slip of the other.
      Four dividing planes (0, 45, 90, 135 deg offsets) are scored
      in both directions.  If the map is real spatial structure it
      must transfer across the split; if it is driven by one
      clustered region it cannot.

  T2  band extrapolation (pooled).  The map is fitted excluding the
      mirror cap (theta > 120 deg) and extrapolated into it: the
      bipolar form predicts the lobe's sign to reappear unseen.
      The reverse -- fit on theta > 60 deg, predict the primary
      cap -- is scored identically.  An extrapolation test of
      bipolarity that random folds cannot perform.

  T3  axis-displacement skill curve (pooled).  The assumed axis is
      displaced by delta = 0..90 deg along eight azimuthal
      directions; theta is recomputed per comet from the
      reconstructed aphelion vectors and the five-fold
      out-of-sample skill re-measured.  A real boundary makes the
      skill peak at delta = 0 and decay with displacement; the
      decay scale measures how sharply the comet cohort alone
      localizes the axis.

  T4  fold-partition stability (pooled).  The five-fold CV of
      step 102 is repeated over 200 random partitions; the
      distribution of out-of-sample rho reports the stability of
      the measured skill against the fold draw.

The model is the pre-declared two-parameter form of step 089 /
step 102; nothing is re-tuned.  All randomness uses fixed seeds.

Inputs
------
results/step_b30_proper_time_slip.csv   (229 comets: desig, cohort,
                                         theta, dtau_unexplained)
data/raw/code/code_original.html        (CODE original elements)
data/raw/warsaw/warsaw_tablec.dat       (Warsaw original elements)

Outputs
-------
results/step_b68_spatial_holdout.json
results/step_b68_spatial_holdout.csv    (per-comet azimuth + band flags)
results/figures/supplementary/step_b68_spatial_holdout.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, TNO_AXIS, tee_stdout
logger = StepLogger("step_104_spatial_holdout")
tee_stdout(logger)
logger.header("Spatial holdout and axis localization of the slip map")

import csv
import json
import math
import re
from html.parser import HTMLParser
import numpy as np
from scipy.stats import spearmanr, binomtest

SEED = 20260919
K_FOLD = 5
N_PARTITIONS = 200
DISPLACEMENTS = np.arange(0.0, 91.0, 3.0)          # deg
N_DIRS = 8
rng = np.random.default_rng(SEED)

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l),
                     math.cos(b) * math.sin(l), math.sin(b)])


# The stored theta column is measured against the cap-declaration
# axis (34, -13), the extreme-subsample TNO axis used to pre-declare
# the comet cap (steps 030-038, 063-064); the detached-sample axis
# (49, -17) is the second pre-declared anchor (step 061).
AXIS = lv(34.0, -13.0)
AXIS_DET = np.asarray(TNO_AXIS, dtype=float)
AXIS_DET = AXIS_DET / np.linalg.norm(AXIS_DET)


def unit_basis(a):
    """Orthonormal basis (e1, e2) spanning the plane perp to a."""
    e1 = np.cross([0.0, 0.0, 1.0], a)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    return e1, e2


E1, E2 = unit_basis(AXIS)

# ------------------------------------------------------------------
# Catalogue element parsers (same conventions as steps 063/064)
# ------------------------------------------------------------------


class TP(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.cur, self.buf = [], None, ""
        self.in_td = False

    def handle_starttag(self, t, a):
        if t == "tr":
            self.cur = []
        elif t in ("td", "th") and self.cur is not None:
            self.in_td, self.buf = True, ""

    def handle_endtag(self, t):
        if t in ("td", "th") and self.cur is not None:
            self.cur.append(re.sub(r"\s+", " ", self.buf).strip())
            self.in_td = False
        elif t == "tr" and self.cur:
            self.rows.append(self.cur)
            self.cur = None

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
                w=float(r[10]), Om=float(r[11]), i=float(r[12]))
        except (ValueError, IndexError):
            continue
    return out


def parse_warsaw(path):
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(
                com=line[3].strip(), desig=line[5:17].strip(),
                w=float(line[70:82]), Om=float(line[82:94]),
                i=float(line[94:106])))
        except ValueError:
            continue
    return rows


PREF = {"a": 0, "h": 0, "e": 1, "b": 2}


def dedup(rows, pref):
    out = {}
    for r in rows:
        k = r["desig"]
        if k not in out or pref.get(r["com"], 9) < pref.get(out[k]["com"], 9):
            out[k] = r
    return out


def perih_dir(om, Om, inc):
    om, Om, inc = map(math.radians, (om, Om, inc))
    co, so = math.cos(om), math.sin(om)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(inc), math.sin(inc)
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci,
                     so * si])


# ------------------------------------------------------------------
# Load slip table and reconstruct aphelion unit vectors
# ------------------------------------------------------------------

rows = []
with open(RESULTS / "step_b30_proper_time_slip.csv") as f:
    for r in csv.DictReader(f):
        try:
            rows.append(dict(
                desig=r["desig"], cohort=r["cohort"],
                theta=float(r["theta"]),
                dtau=float(r["dtau_unexplained"])))
        except (ValueError, KeyError):
            continue
n_all = len(rows)

code_el = parse_code(str(DATA_RAW / "code" / "code_original.html"))
war_el = dedup(parse_warsaw(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")), PREF)
logger.info(f"element tables: code {len(code_el)}, warsaw {len(war_el)}")

U = np.full((n_all, 3), np.nan)
n_missing = 0
for j, r in enumerate(rows):
    el = code_el.get(r["desig"]) if r["cohort"] == "code" \
        else war_el.get(r["desig"])
    if el is None:
        n_missing += 1
        continue
    U[j] = -perih_dir(el["w"], el["Om"], el["i"])   # aphelion = -perihelion
have_vec = ~np.isnan(U[:, 0])
logger.info(f"aphelion vectors reconstructed: {int(have_vec.sum())}/{n_all}"
            + (f" ({n_missing} unmatched)" if n_missing else ""))

theta = np.array([r["theta"] for r in rows])
y = np.array([r["dtau"] for r in rows])
cohort = np.array([r["cohort"] for r in rows])

# sanity: recomputed theta must reproduce the stored catalogue values
th_check = np.degrees(np.arccos(np.clip(U[have_vec] @ AXIS, -1, 1)))
med_err = float(np.median(np.abs(th_check - theta[have_vec])))
logger.metric("theta_recompute_med_err_deg", med_err)
if med_err > 1.0:
    raise RuntimeError(
        f"reconstructed aphelion vectors do not reproduce theta "
        f"(median error {med_err:.2f} deg); aborting")

# azimuth of each aphelion around the axis
proj = U - np.outer(U @ AXIS, AXIS)
azim = np.degrees(np.arctan2(proj @ E2, proj @ E1)) % 360.0
c2 = np.cos(2.0 * np.radians(theta))


def fit_map(feat, yy):
    A = np.vstack([np.ones_like(feat), feat]).T
    c, *_ = np.linalg.lstsq(A, yy, rcond=None)
    return c


def skill(tr_idx, te_idx, feat, yy):
    """Fit a+b*feat on tr, predict te, return (rho, sign metrics)."""
    c = fit_map(feat[tr_idx], yy[tr_idx])
    pr = c[0] + c[1] * feat[te_idx]
    rho, p = spearmanr(pr, yy[te_idx])
    n_ok = int(np.sum((pr > 0) == (yy[te_idx] > 0)))
    n = int(len(te_idx))
    return dict(rho=float(rho), rho_p=float(p), n=n,
                sign_correct=n_ok, sign_frac=n_ok / n,
                sign_binom_p=float(binomtest(n_ok, n, 0.5,
                                             alternative="greater").pvalue))


def oob_cv(feat, yy, k=K_FOLD, rng_local=None):
    """k-fold out-of-sample Spearman (no permutation null)."""
    r_ = rng_local if rng_local is not None else rng
    folds = np.empty(len(yy), dtype=int)
    idx = np.arange(len(yy))
    r_.shuffle(idx)
    for j, part in enumerate(np.array_split(idx, k)):
        folds[part] = j
    pred = np.full(len(yy), np.nan)
    for kk in range(k):
        tr, te = folds != kk, folds == kk
        c = fit_map(feat[tr], yy[tr])
        pred[te] = c[0] + c[1] * feat[te]
    return float(spearmanr(pred, yy)[0])


# ------------------------------------------------------------------
# T1: azimuthal spatial holdout (pooled, vector-bearing rows)
# ------------------------------------------------------------------

idx_v = np.where(have_vec)[0]
t1 = []
for off in (0.0, 45.0, 90.0, 135.0):
    half = ((azim + off) % 360.0) >= 180.0
    for direction, tr_m in (("A->B", ~half), ("B->A", half)):
        tr = idx_v[tr_m[idx_v]]
        te = idx_v[~tr_m[idx_v]]
        if len(tr) < 20 or len(te) < 20:
            continue
        s = skill(tr, te, c2, y)
        s.update(offset_deg=off, direction=direction,
                 n_train=int(len(tr)), n_test=int(len(te)))
        t1.append(s)
        logger.metric(f"azhold_{int(off)}_{direction.replace('->','_')}",
                      f"rho={s['rho']:+.3f} (p={s['rho_p']:.4f}), "
                      f"sign {s['sign_correct']}/{s['n']} "
                      f"(p={s['sign_binom_p']:.3f})")

rho_t1 = np.array([s["rho"] for s in t1])

# ------------------------------------------------------------------
# T2: band extrapolation (pooled)
# ------------------------------------------------------------------

t2 = {}
for name, tr_m, te_m in (
        ("fit_le120_predict_mirror", theta <= 120.0, theta > 120.0),
        ("fit_ge60_predict_cap", theta >= 60.0, theta < 60.0)):
    tr, te = np.where(tr_m & have_vec)[0], np.where(te_m & have_vec)[0]
    s = skill(tr, te, c2, y)
    c_ = fit_map(c2[tr], y[tr])
    s.update(n_train=int(len(tr)), n_test=int(len(te)),
             obs_median=float(np.median(y[te])),
             pred_median=float(np.median(c_[0] + c_[1] * c2[te])))
    t2[name] = s
    logger.metric(name, f"rho={s['rho']:+.3f} (p={s['rho_p']:.4f}), "
                  f"sign {s['sign_correct']}/{s['n']}, "
                  f"pred med {s['pred_median']:+.2f} vs "
                  f"obs {s['obs_median']:+.2f} yr")

# ------------------------------------------------------------------
# T3: axis-displacement skill curve (pooled)
# ------------------------------------------------------------------

dirs = [E1 * np.cos(2 * np.pi * j / N_DIRS)
        + E2 * np.sin(2 * np.pi * j / N_DIRS)
        for j in range(N_DIRS)]

curve = np.full((N_DIRS, len(DISPLACEMENTS)), np.nan)
for j, d in enumerate(dirs):
    for m, delta in enumerate(DISPLACEMENTS):
        ax = AXIS * np.cos(np.radians(delta)) \
            + d * np.sin(np.radians(delta))
        ax /= np.linalg.norm(ax)
        th_d = np.degrees(np.arccos(
            np.clip(U[have_vec] @ ax, -1, 1)))
        feat = np.full(n_all, np.nan)
        feat[have_vec] = np.cos(2.0 * np.radians(th_d))
        rng_loc = np.random.default_rng(SEED + 7)
        curve[j, m] = oob_cv(feat[have_vec], y[have_vec],
                             rng_local=rng_loc)

curve_mean = np.nanmean(curve, axis=0)
curve_lo = np.nanmin(curve, axis=0)
curve_hi = np.nanmax(curve, axis=0)
rho0 = float(curve_mean[0])
izero = np.where(curve_mean <= 0)[0]
delta_zero = float(DISPLACEMENTS[izero[0]]) if len(izero) else 91.0
frac_dirs_neg_at60 = float(np.mean(curve[:, int(60 / 3)] <= 0))
# localization: displacement at which mean skill halves
ihalf = np.where(curve_mean <= rho0 / 2)[0]
delta_half = float(DISPLACEMENTS[ihalf[0]]) if len(ihalf) else 91.0

# second pre-declared anchor: skill with theta measured against the
# detached-sample axis (49, -17)
th_det = np.degrees(np.arccos(np.clip(U[have_vec] @ AXIS_DET, -1, 1)))
feat_det = np.cos(2.0 * np.radians(th_det))
rng_loc = np.random.default_rng(SEED + 7)
rho_det = oob_cv(feat_det, y[have_vec], rng_local=rng_loc)
sep_axes = float(np.degrees(np.arccos(np.clip(AXIS @ AXIS_DET, -1, 1))))
logger.metric("disp_skill_at0", f"{rho0:+.3f}")
logger.metric("skill_detached_axis", f"{rho_det:+.3f} "
              f"(axis {sep_axes:.1f} deg away)")
logger.metric("disp_zero_crossing_deg", delta_zero)
logger.metric("disp_half_skill_deg", delta_half)
logger.metric("disp_frac_dirs_neg_60deg", frac_dirs_neg_at60)

# ------------------------------------------------------------------
# T4: fold-partition stability (pooled)
# ------------------------------------------------------------------

rho_dist = np.empty(N_PARTITIONS)
for i in range(N_PARTITIONS):
    rng_loc = np.random.default_rng(SEED + 1000 + i)
    rho_dist[i] = oob_cv(c2[have_vec], y[have_vec], rng_local=rng_loc)
q16, q50, q84 = np.percentile(rho_dist, [16, 50, 84])
frac_pos = float(np.mean(rho_dist > 0))
logger.metric("cv_stability",
              f"median {q50:+.3f} [{q16:+.3f},{q84:+.3f}], "
              f"frac>0 {frac_pos:.3f} over {N_PARTITIONS} partitions")

# ------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------

verdict = dict(
    spatial_transfer=("positive" if np.median(rho_t1) > 0 else "null"),
    mirror_extrapolation=(
        "bipolar sign recovered out-of-sample"
        if t2["fit_le120_predict_mirror"]["sign_binom_p"] < 0.1
        else "not resolved"),
    localization=f"CV skill at zero displacement {rho0:+.3f}; "
                 f"half-skill near {delta_half:.0f} deg, zero "
                 f"crossing near {delta_zero:.0f} deg",
    reading=(
        f"Azimuthal half-sky holdout: median rho {np.median(rho_t1):+.3f} "
        f"over {len(t1)} split/direction tests (range "
        f"{rho_t1.min():+.3f}..{rho_t1.max():+.3f}).  Fitted without "
        f"the mirror cap, the map predicts the held-out lobe at "
        f"sign accuracy "
        f"{t2['fit_le120_predict_mirror']['sign_frac']:.2f} "
        f"(p={t2['fit_le120_predict_mirror']['sign_binom_p']:.3f}); "
        f"the mirror cap's observed median "
        f"{t2['fit_le120_predict_mirror']['obs_median']:+.2f} yr "
        f"against predicted "
        f"{t2['fit_le120_predict_mirror']['pred_median']:+.2f} yr.  "
        f"Displacing the axis, the pooled cohort's out-of-sample "
        f"skill is {rho0:+.3f} at zero displacement, halves near "
        f"{delta_half:.0f} deg and crosses zero near "
        f"{delta_zero:.0f} deg; the detached-sample axis "
        f"({sep_axes:.0f} deg away) returns {rho_det:+.3f}.  "
        f"Over {N_PARTITIONS} random "
        f"partitions the CV skill is positive in "
        f"{100*frac_pos:.0f} per cent of draws (median {q50:+.3f})."))

# ------------------------------------------------------------------
# Write
# ------------------------------------------------------------------

res = dict(
    step="step_104_spatial_holdout",
    description=("Spatial holdout and axis-localization validation of "
                 "the bipolar slip map: azimuthal half-sky holdout, "
                 "mirror-cap and primary-cap band extrapolation, an "
                 "axis-displacement out-of-sample skill curve, and "
                 "fold-partition stability.  Aphelion unit vectors "
                 "are reconstructed from the catalogues' "
                 "original-orbit elements (step-063 convention); "
                 "model dtau = a + b*cos(2theta), pre-declared "
                 "(step 089)."),
    inputs=["results/step_b30_proper_time_slip.csv",
            "data/raw/code/code_original.html",
            "data/raw/warsaw/warsaw_tablec.dat"],
    seed=SEED,
    n=n_all, n_vectors=int(have_vec.sum()),
    theta_recompute_med_err_deg=med_err,
    T1_azimuthal_holdout=dict(tests=t1, median_rho=float(np.median(rho_t1))),
    T2_band_extrapolation=t2,
    T3_axis_displacement=dict(
        displacements_deg=DISPLACEMENTS.tolist(), n_dirs=N_DIRS,
        curve_mean=curve_mean.tolist(), curve_min=curve_lo.tolist(),
        curve_max=curve_hi.tolist(), skill_at_zero=rho0,
        half_skill_deg=delta_half, zero_crossing_deg=delta_zero,
        frac_dirs_le0_at_60deg=frac_dirs_neg_at60,
        detached_axis=dict(lon=49.0, lat=-17.0,
                           sep_from_cap_axis_deg=sep_axes,
                           cv_skill=rho_det)),
    T4_partition_stability=dict(
        n_partitions=N_PARTITIONS,
        rho_median=float(q50), rho_p16=float(q16), rho_p84=float(q84),
        frac_positive=frac_pos),
    verdict=verdict,
    caveats=[
        "The azimuthal split tests field-vs-cluster structure; it "
        "does not make the two halves statistically independent "
        "since comets near the dividing plane are spatially close "
        "across it.",
        "The displacement curve uses the same two-parameter map; it "
        "measures where the comet cohort places the axis under the "
        "pre-declared model, not a free axis search.",
        "Mirror-cap extrapolation is a single-region holdout; its "
        "sign test has modest power at n ~ 30."])

out = RESULTS / "step_b68_spatial_holdout.json"
json.dump(res, open(out, "w"), indent=1, default=float)
logger.data_save(out)

csv_out = RESULTS / "step_b68_spatial_holdout.csv"
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["desig", "cohort", "theta", "dtau_unexplained",
                "azimuth_deg", "band"])
    for i, r in enumerate(rows):
        az = float(azim[i]) if have_vec[i] else ""
        band = ("cap" if r["theta"] < 60 else
                "mid" if r["theta"] <= 120 else "mirror")
        w.writerow([r["desig"], r["cohort"], r["theta"], r["dtau"],
                    az, band])
logger.data_save(csv_out)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.4))

ax = axes[0]
sc = ax.scatter(azim[have_vec], np.clip(y[have_vec], -60, 60),
                s=14, c=theta[have_vec], cmap="coolwarm",
                alpha=0.65, edgecolor="none")
ax.axhline(0, color="0.5", lw=0.8)
ax.set_xlabel("aphelion azimuth about the axis (deg)")
ax.set_ylabel("unexplained slip $\\delta\\tau$ (yr)")
ax.set_title(f"azimuthal halves: median $\\rho$ "
             f"{np.median(rho_t1):+.3f}", fontsize=10)
cb = fig.colorbar(sc, ax=ax, pad=0.01)
cb.set_label("$\\theta$ (deg)", fontsize=8)

ax = axes[1]
th_g = np.linspace(0, 180, 400)
m = (theta <= 120.0) & have_vec
c_ = fit_map(c2[m], y[m])
ax.scatter(theta[m], np.clip(y[m], -60, 60), s=12, c="0.5",
           alpha=0.45, label="fitted ($\\theta\\leq120$)")
ax.scatter(theta[~m & have_vec],
           np.clip(y[~m & have_vec], -60, 60), s=22, c="crimson",
           alpha=0.75, label="held-out mirror cap")
ax.plot(th_g, c_[0] + c_[1] * np.cos(2 * np.radians(th_g)), "k-",
        lw=1.5, label="extrapolated map")
ax.axvspan(120, 180, color="crimson", alpha=0.05)
ax.axhline(0, color="0.5", lw=0.7, ls=":")
ax.set_xlabel("transit angle $\\theta$ (deg)")
ax.set_ylabel("$\\delta\\tau$ (yr)")
s2 = t2["fit_le120_predict_mirror"]
ax.set_title(f"mirror-cap extrapolation: sign "
             f"{s2['sign_correct']}/{s2['n']} "
             f"(p={s2['sign_binom_p']:.3f})", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="lower left")

ax = axes[2]
ax.fill_between(DISPLACEMENTS, curve_lo, curve_hi,
                color="steelblue", alpha=0.25, label="8 azimuths")
ax.plot(DISPLACEMENTS, curve_mean, "o-", color="steelblue",
        ms=3, lw=1.4, label="mean CV skill")
ax.axhline(0, color="0.5", lw=0.8)
ax.axvline(0, color="crimson", lw=1.0, ls="--",
           label="measured axis")
ax.set_xlabel("axis displacement $\\delta$ (deg)")
ax.set_ylabel("out-of-sample $\\rho$ (pooled)")
ax.set_title(f"axis localization: half-skill "
             f"{delta_half:.0f} deg, zero {delta_zero:.0f} deg",
             fontsize=10)
ax.legend(frameon=False, fontsize=8)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b68_spatial_holdout.png", dpi=300)
logger.data_save(FIG / "supplementary" / "step_b68_spatial_holdout.png")
logger.success("Spatial holdout and axis localization complete")
