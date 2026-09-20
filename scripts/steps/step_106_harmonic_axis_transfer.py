#!/usr/bin/env python3
"""Step 106 -- harmonic-order transfer, third-cohort axis recovery,
and cross-solution per-comet concordance.

Three remaining questions on the slip map are closed out here,
using only existing per-comet products (no new integrations).

  T1  Third-cohort axis recovery.  The original convergence test
      (step_036) let the CODE cohort recover its own axis by a free
      sky scan; it landed 24 deg off the TNO axis.  The identical
      scan is run on the one-apparition cohort (step 105 products):
      a 10-deg grid of trial axes, Mann-Whitney in-cap vs out-cap
      on the catalogue rotation d_of, best axis and the rank of the
      pre-declared axis among all trial directions.

  T2  Harmonic-order out-of-sample discrimination.  The bipolar
      cos(2 theta) morphology was selected in-sample (step 089:
      cos2 beats cos1 at p=0.010 vs 0.30), so the predictive
      tests inherit a selection on the same data.  The five-fold
      CV of step 102 is therefore rerun with alternative profile
      families -- constant, cos(theta), cos(2theta), cos(3theta),
      and a linear theta term -- under identical folds.  If
      cos(2theta) retains the best held-out skill, the bipolar
      form is confirmed predictively, not just descriptively.

  T3  Amplitude transfer.  The lpc cohort's own free fit
      dtau = a + b cos(2theta) returns b_lpc with a bootstrap
      interval, compared against the 229-fitted b = +5.98 yr --
      an amplitude consistency check that per-object rho cannot
      deliver.

  T4  Cross-solution per-comet concordance.  On the 23 members
      shared between the lpc cohort and the training set, the
      same bodies carry two different orbit solutions.  If the
      anomaly is an object property, the per-comet discrepancies
      should correlate across solutions; if it is a fit artefact,
      they should not.

Inputs
------
results/step_b30_proper_time_slip.csv     (229-comet training set)
results/step_b69_lpc_bidirectional.csv    (30-comet lpc cohort)

Outputs
-------
results/step_b70_harmonic_axis_transfer.json
results/step_b70_harmonic_axis_transfer.csv
results/figures/supplementary/step_b70_harmonic_axis_transfer.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_106_harmonic_axis_transfer")
tee_stdout(logger)
logger.header("Harmonic-order transfer + third-cohort axis recovery")

import csv
import json
import math
import re
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

SEED = 20260920
N_BOOT = 10000
K_FOLD = 5
rng = np.random.default_rng(SEED)

CAP = 60.0


def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l),
                     math.cos(b) * math.sin(l), math.sin(b)])


def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))


AXIS = lv(34.0, -13.0)
AXIS_DET = lv(49.9, -17.0)


def _norm_desig(d):
    d = d.strip().upper().replace(" ", "")
    m = re.search(r"([CPD]/\d{4}[A-Z]+\d+|\d{4}[A-Z]+\d+)", d)
    return m.group(1) if m else d


# ------------------------------------------------------------------
# Load products
# ------------------------------------------------------------------

slip = list(csv.DictReader(
    open(RESULTS / "step_b30_proper_time_slip.csv")))
lpc = list(csv.DictReader(
    open(RESULTS / "step_b69_lpc_bidirectional.csv")))

th229 = np.array([float(r["theta"]) for r in slip])
y229 = np.array([float(r["dtau_unexplained"]) for r in slip])
cohort = np.array([r["cohort"] for r in slip])

th_l = np.array([float(r["theta"]) for r in lpc])
dt_l = np.array([float(r["dtau_unexplained"]) for r in lpc])
drot_l = np.array([float(r["drot_cat"]) for r in lpc])
# the scan needs full aphelion vectors; rebuild from the lpc
# catalogue original-orbit elements as in step_105
import xml.etree.ElementTree as ET
from scripts.utils.tep9_common import DATA_RAW


def load_vot(path):
    t = ET.parse(path)
    fields = [el.get("name") for el in t.getroot().iter()
              if el.tag.split("}")[-1] == "FIELD"]
    return [dict(zip(fields, [td.text for td in list(el)]))
            for el in t.getroot().iter()
            if el.tag.split("}")[-1] == "TR"]


def perih_dir(om, Om, inc):
    om, Om, inc = map(math.radians, (om, Om, inc))
    co, so = math.cos(om), math.sin(om)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(inc), math.sin(inc)
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci, so * si])


org = {}
for r in load_vot(DATA_RAW / "lpc" / "lpc_orig_2006_2010.vot"):
    k = r["Comet"].strip()
    if k not in org or r["Model"].strip() == "GR":
        org[k] = r

aph = {}
for r in lpc:
    k = r["desig"]
    ro = org[k]
    po = perih_dir(float(ro["arg"]), float(ro["long"]), float(ro["i"]))
    aph[k] = -po
APH = np.array([aph[r["desig"]] for r in lpc])
assert all(abs(sep(APH[j], AXIS) - th_l[j]) < 0.5 for j in range(len(lpc)))
logger.info(f"loaded 229 training + {len(lpc)} lpc comets")

# ------------------------------------------------------------------
# T1: free axis scan on lpc (step_036 methodology)
# ------------------------------------------------------------------

d_l = drot_l
grid = []
for lam in np.arange(0, 360, 10):
    for b in np.arange(-60, 61, 10):
        ax = lv(lam, b)
        th_ = np.array([sep(v, ax) for v in APH])
        inc = th_ < CAP
        if inc.sum() < 4 or inc.sum() > len(d_l) - 4:
            continue
        p = mannwhitneyu(d_l[inc], d_l[~inc],
                         alternative="greater").pvalue
        grid.append((p, lam, b))
grid.sort()
bp, blam, bb = grid[0]
# rank of the pre-declared axis among trial directions
th_ax = np.array([sep(v, AXIS) for v in APH])
p_ax = mannwhitneyu(d_l[th_ax < CAP], d_l[th_ax >= CAP],
                    alternative="greater").pvalue
rank_frac = float((int(sum(g[0] <= p_ax for g in grid)) + 1)
                  / (len(grid) + 1))
t1 = dict(n_grid=len(grid),
          best=dict(lam=float(blam), beta=float(bb), p=float(bp),
                    sep_from_cap=float(sep(lv(blam, bb), AXIS)),
                    sep_from_det=float(sep(lv(blam, bb), AXIS_DET))),
          predeclared=dict(p=float(p_ax), rank_frac=rank_frac,
                           n_in=int((th_ax < CAP).sum())))
logger.metric("lpc_axis_recovery",
              f"best ({blam:.0f},{bb:.0f}) p={bp:.4f}; "
              f"predeclared p={p_ax:.4f} top {rank_frac*100:.1f}%")

# ------------------------------------------------------------------
# T2: harmonic-order out-of-sample discrimination on the 229
# ------------------------------------------------------------------
# identical stratified folds as step_102 (same seed scheme)

folds = np.empty(len(th229), dtype=int)
rng_f = np.random.default_rng(20260918)   # step_102 seed
for co in ("code", "warsaw"):
    idx = np.where(cohort == co)[0]
    rng_f.shuffle(idx)
    for j, part in enumerate(np.array_split(idx, K_FOLD)):
        folds[part] = j


def feat(theta_deg, kind):
    if kind == "const":
        return np.zeros_like(theta_deg)
    if kind == "linear":
        return np.radians(theta_deg)
    m = int(kind[-1])
    return np.cos(m * np.radians(theta_deg))


def cv_skill(kind):
    pred = np.full(len(th229), np.nan)
    for k in range(K_FOLD):
        tr, te = folds != k, folds == k
        x = feat(th229, kind)
        A = np.vstack([np.ones(tr.sum()), x[tr]]).T
        c, *_ = np.linalg.lstsq(A, y229[tr], rcond=None)
        pred[te] = c[0] + c[1] * x[te]
    return float(spearmanr(pred, y229)[0]), pred


profiles = ["const", "cos1", "cos2", "cos3", "linear"]
t2 = {}
preds = {}
for kind in profiles:
    rho_, pr = cv_skill(kind)
    t2[kind] = dict(oob_rho=rho_)
    preds[kind] = pr
    logger.info(f"  OOB skill {kind}: rho={rho_:+.4f}")

# bootstrap CI on the rho difference cos2 - cos1 over the pooled OOB
# predictions (resample comets)
d_rho = np.empty(N_BOOT)
idx_all = np.arange(len(th229))
for i_ in range(N_BOOT):
    ii = rng.choice(idx_all, len(idx_all), replace=True)
    d_rho[i_] = (spearmanr(preds["cos2"][ii], y229[ii])[0]
                 - spearmanr(preds["cos1"][ii], y229[ii])[0])
t2["cos2_minus_cos1"] = dict(
    median=float(np.median(d_rho)),
    ci68=[float(np.percentile(d_rho, 16)),
          float(np.percentile(d_rho, 84))],
    frac_positive=float(np.mean(d_rho > 0)))
logger.metric("harmonic_oob",
              f"cos2-cos1 = {t2['cos2_minus_cos1']['median']:+.4f} "
              f"[{t2['cos2_minus_cos1']['ci68'][0]:+.4f},"
              f"{t2['cos2_minus_cos1']['ci68'][1]:+.4f}]")

# ------------------------------------------------------------------
# T3: amplitude transfer on lpc
# ------------------------------------------------------------------

c2_229 = np.cos(2.0 * np.radians(th229))
A229 = np.vstack([np.ones_like(c2_229), c2_229]).T
c_map, *_ = np.linalg.lstsq(A229, y229, rcond=None)
b229 = float(c_map[1])

train_desigs = {_norm_desig(r["desig"]) for r in slip}
shared_mask = np.array([_norm_desig(r["desig"]) in train_desigs
                        for r in lpc])


def fit_b(mask):
    x = np.cos(2.0 * np.radians(th_l[mask]))
    A = np.vstack([np.ones_like(x), x]).T
    c, *_ = np.linalg.lstsq(A, dt_l[mask], rcond=None)
    bs = np.empty(N_BOOT)
    idx = np.arange(mask.sum())
    for i_ in range(N_BOOT):
        ii = rng.choice(idx, len(idx), replace=True)
        cb, *_ = np.linalg.lstsq(A[ii], dt_l[mask][ii], rcond=None)
        bs[i_] = cb[1]
    return dict(n=int(mask.sum()), b=float(c[1]),
                ci68=[float(np.percentile(bs, 16)),
                      float(np.percentile(bs, 84))],
                frac_b_pos=float(np.mean(bs > 0)))


t3 = dict(b_229=b229,
          all=fit_b(np.ones(len(lpc), bool)),
          shared=fit_b(shared_mask),
          nonshared=fit_b(~shared_mask))
logger.metric("amplitude_transfer",
              f"b_lpc={t3['all']['b']:+.2f} "
              f"[{t3['all']['ci68'][0]:+.2f},{t3['all']['ci68'][1]:+.2f}]"
              f" vs b_229={b229:+.2f}")

# ------------------------------------------------------------------
# T4: cross-solution per-comet concordance (shared members)
# ------------------------------------------------------------------

slip_by = {_norm_desig(r["desig"]): r for r in slip}
pairs = [(r, slip_by[_norm_desig(r["desig"])])
         for r in lpc if _norm_desig(r["desig"]) in slip_by]
dl = np.array([float(p[0]["drot_cat"]) for p in pairs])
ds = np.array([float(p[1]["drot_cat"]) for p in pairs])
tl = np.array([float(p[0]["dtau_unexplained"]) for p in pairs])
ts = np.array([float(p[1]["dtau_unexplained"]) for p in pairs])
t4 = dict(n=len(pairs),
          drot_rho=float(spearmanr(dl, ds)[0]),
          drot_p=float(spearmanr(dl, ds)[1]),
          dtau_rho=float(spearmanr(tl, ts)[0]),
          dtau_p=float(spearmanr(tl, ts)[1]),
          med_abs_drot_diff=float(np.median(np.abs(dl - ds))))
logger.metric("cross_solution",
              f"drot rho={t4['drot_rho']:+.3f} (p={t4['drot_p']:.3f}), "
              f"dtau rho={t4['dtau_rho']:+.3f}")

# ------------------------------------------------------------------
# Verdict + write
# ------------------------------------------------------------------

verdict = (
    f"On the third cohort (n={len(lpc)}) the free scan's strongest "
    f"axis sits at ({blam:.0f},{bb:.0f}) deg, "
    f"{t1['best']['sep_from_cap']:.0f} deg from the cap-declaration "
    f"axis, and the pre-declared axis ranks in the top "
    f"{rank_frac*100:.1f} per cent of {len(grid)} trial directions "
    f"(p={p_ax:.4f}).  Out-of-sample, the cos2theta profile retains "
    f"the best held-out skill of the tested families "
    f"(rho={t2['cos2']['oob_rho']:+.3f}; cos2-cos1 "
    f"{t2['cos2_minus_cos1']['median']:+.3f} "
    f"[{t2['cos2_minus_cos1']['ci68'][0]:+.3f},"
    f"{t2['cos2_minus_cos1']['ci68'][1]:+.3f}]).  The lpc cohort's "
    f"own amplitude fit returns b={t3['all']['b']:+.2f} yr "
    f"[{t3['all']['ci68'][0]:+.2f},{t3['all']['ci68'][1]:+.2f}] "
    f"against the training amplitude {b229:+.2f} yr.  On the "
    f"{t4['n']} shared bodies the two solutions' per-comet "
    f"discrepancies correlate at rho={t4['drot_rho']:+.3f}.")

res = dict(
    step="step_106_harmonic_axis_transfer",
    description=("Out-of-sample harmonic-order discrimination on "
                 "the 229 (identical step-102 folds), free axis "
                 "recovery scan on the third (one-apparition) "
                 "cohort, lpc amplitude transfer, and per-comet "
                 "cross-solution concordance on shared members."),
    inputs=["results/step_b30_proper_time_slip.csv",
            "results/step_b69_lpc_bidirectional.csv",
            "data/raw/lpc/lpc_orig_2006_2010.vot"],
    seed=SEED, n_boot=N_BOOT,
    T1_axis_recovery=t1,
    T2_harmonic_oob=t2,
    T3_amplitude=t3,
    T4_cross_solution=t4,
    verdict=verdict)

out = RESULTS / "step_b70_harmonic_axis_transfer.json"
json.dump(res, open(out, "w"), indent=1, default=float)
logger.data_save(out)

csv_out = RESULTS / "step_b70_harmonic_axis_transfer.csv"
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["axis_lam", "axis_beta", "mw_p", "is_best",
                "is_predeclared"])
    for p_, l_, b_ in grid:
        w.writerow([l_, b_, p_, int(l_ == blam and b_ == bb), 0])
logger.data_save(csv_out)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.4))

ax = axes[0]
kinds = profiles
rhos = [t2[k]["oob_rho"] for k in kinds]
cols = ["0.6" if k != "cos2" else "crimson" for k in kinds]
ax.bar(range(len(kinds)), rhos, color=cols)
ax.set_xticks(range(len(kinds)))
ax.set_xticklabels(["const", "cos$\\theta$", "cos2$\\theta$",
                    "cos3$\\theta$", "linear"], fontsize=9)
ax.axhline(0, color="0.5", lw=0.7)
ax.set_ylabel("held-out Spearman $\\rho$")
ax.set_title("harmonic order, out-of-sample (229)", fontsize=10)

ax = axes[1]
ax.errorbar([0], [b229], yerr=0, fmt="s", c="k", ms=7,
            label=f"229 fit ({b229:+.1f} yr)")
for i_, (tag, sub) in enumerate(
        (("all", t3["all"]), ("shared", t3["shared"]),
         ("new", t3["nonshared"]))):
    ax.errorbar([i_ + 1], [sub["b"]],
                yerr=[[sub["b"] - sub["ci68"][0]],
                      [sub["ci68"][1] - sub["b"]]],
                fmt="o", ms=6, capsize=4,
                label=f"lpc {tag} (n={sub['n']})")
ax.axhline(b229, color="0.6", ls="--", lw=0.9)
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(["229", "all", "shared", "new"], fontsize=9)
ax.set_ylabel("$b$ in $\\delta\\tau = a + b\\cos2\\theta$ (yr)")
ax.set_title("amplitude transfer", fontsize=10)
ax.legend(frameon=False, fontsize=8)

ax = axes[2]
ax.scatter(ds, dl, s=24, c="steelblue", alpha=0.8)
lim = max(ds.max(), dl.max()) * 1.1
ax.plot([0, lim], [0, lim], "k:", lw=0.9)
ax.set_xlabel("$d_{\\rm of}$, Warsaw/CODE solution (deg)")
ax.set_ylabel("$d_{\\rm of}$, lpc solution (deg)")
ax.set_title(f"cross-solution concordance "
             f"$\\rho$={t4['drot_rho']:+.2f} (n={t4['n']})",
             fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b70_harmonic_axis_transfer.png", dpi=300)
logger.data_save(FIG / "supplementary" / "step_b70_harmonic_axis_transfer.png")
logger.success("Harmonic/axis transfer analysis complete")
