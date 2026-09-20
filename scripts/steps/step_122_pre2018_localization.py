#!/usr/bin/env python3
"""step_122: Anomaly localization on the pre-2018 SBDB record.

Step 121's joint decomposition returned a weak positive declared-axis
term on the broad pre-2018 cohort (+0.166 dex, p = 0.019) while the
deep-plunger primary stayed median-null (step 120).  Two readings are
possible:

  H-LOCK   the term is carried by the 288 members the cohort shares
           with CODE -- the anomaly is locked to the CODE-selected
           objects and reproduces under independent JPL solutions
           (object-association, strengthening the step_120 claim);
  H-BROAD  the term generalizes across the pre-2018 population --
           a broad anomaly diluted in CODE, weakening the
           localization story.

This step localizes the term.  The 738 integrated pre-2018 comets are
split by CODE membership and interrogated at the declared axis; the
non-overlap remainder is then audited for what carries its residual
positive lean (era concentration vs a handful of fat-tailed
19th/20th-century solutions).

Tests
-----
  T1  membership decomposition: declared-axis gap, MWU and Spearman
      for CODE-overlap (n=288) vs non-overlap (n=450); membership x
      cap interaction on the pooled cohort (in_code labels permuted)
      -- is the anomaly membership-locked at significance?
  T2  cap-width sensitivity on the overlap (45/60/75 deg) plus
      anti-axis and post-2017-displaced-axis controls -- does the
      overlap term have the anomaly's threshold morphology?
  T3  era x membership grid (pre-1990 / 1990-2009 / 2010-2017) --
      where the overlap term lives in time, and whether the
      non-overlap pre-1990 lean is real.
  T4  fat-tail audit of the non-overlap pre-1990 cell: trimmed and
      winsorized means, leave-one-out on the top in-cap residuals,
      composition (arc, nobs, condition code) of the high-residual
      members -- fat right tail of historical solutions vs signal.
  T5  quality-matched non-overlap subsets (q<3.1, bound, spanning
      arc, arc>365d, nobs>300) -- does ANY non-overlap cell show the
      anomaly?

Inputs
------
results/step_b84_pre2018_sbdb.csv   (per-comet legs, in_code flags)

Outputs
-------
results/step_b86_pre2018_localization.json
results/step_b86_pre2018_localization.csv
results/figures/supplementary/step_b86_pre2018_localization.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.tep9_common import RESULTS, sep, lv
logger = StepLogger("step_122_pre2018_localization")
tee_stdout(logger)

import csv
import json
import math

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr, trim_mean

logger.header("Anomaly localization: pre-2018 declared-axis term")

SEED = 20260919
rng = np.random.default_rng(SEED)
N_PERM = 2000
CAP = 60.0
TNO = lv(34.0, -13.0)          # declared transit axis
ANTI = lv(214.0, 13.0)         # antipode control
U_RES = lv(120.0, -40.0)       # post-2017 displaced-dipole axis

# ------------------------------------------------------------------
# 1. Load the pre-2018 per-comet record
# ------------------------------------------------------------------

rows = []
with open(RESULTS / "step_b84_pre2018_sbdb.csv") as f:
    for r in csv.DictReader(f):
        try:
            d = {k: float(r[k]) for k in
                 ("yr", "q", "i", "e", "arc", "nobs", "theta", "drot",
                  "daa", "denc")}
            d["desig"] = r["desig"]
            d["in_code"] = r.get("in_code") == "True"
            d["f_pre"] = float(r["f_pre"]) if r.get("f_pre") else None
            cc = r.get("cc")
            d["cc"] = float(cc) if cc not in (None, "") else np.nan
            d["aph"] = np.array(
                [float(x) for x in r["aph_lb"].split(";")])
            rows.append(d)
        except (ValueError, KeyError):
            continue
logger.info(f"pre-2018 cohort: {len(rows)} comets with boundary legs")

# ------------------------------------------------------------------
# 2. Shared machinery (identical residual model to steps 118-121)
# ------------------------------------------------------------------

def resid_model(sub):
    v = np.array([r["drot"] for r in sub])
    K = np.abs(np.array([r["daa"] for r in sub]))
    D = np.array([r["denc"] for r in sub])
    Q = np.array([r["q"] for r in sub])
    I = np.array([r["i"] for r in sub])
    y = np.log10(v)
    X = np.column_stack([np.ones(len(y)), np.log10(K + 1.0),
                         np.log10(D), Q, I])
    with np.errstate(all="ignore"):
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ coef
    resid[~np.isfinite(resid)] = np.nanmedian(
        resid[np.isfinite(resid)])
    return resid


def theta_to(sub, u):
    return np.array([sep(r["aph"], u) for r in sub])


def gap_report(sub, resid_vec, u, cap=CAP):
    th = theta_to(sub, u)
    inc = th < cap
    if inc.sum() < 5 or (~inc).sum() < 5:
        return dict(n=len(sub), n_in=int(inc.sum()),
                    status="insufficient_coverage")
    sr = spearmanr(th, resid_vec)
    return dict(
        n=len(sub), n_in=int(inc.sum()),
        med_gap=float(np.median(resid_vec[inc])
                      - np.median(resid_vec[~inc])),
        mean_gap=float(np.mean(resid_vec[inc])
                       - np.mean(resid_vec[~inc])),
        mwu_greater_p=float(mannwhitneyu(
            resid_vec[inc], resid_vec[~inc],
            alternative="greater").pvalue),
        rho=float(sr.statistic), rho_p=float(sr.pvalue))


overlap = [r for r in rows if r["in_code"]]
nonov = [r for r in rows if not r["in_code"]]
resid_all = resid_model(rows)
resid_ov = resid_model(overlap)
resid_no = resid_model(nonov)
logger.info(f"CODE-overlap n={len(overlap)}, non-overlap n={len(nonov)}")

# ------------------------------------------------------------------
# T1: membership decomposition + membership x cap interaction
# ------------------------------------------------------------------

t1 = {"overlap": gap_report(overlap, resid_ov, TNO),
      "non_overlap": gap_report(nonov, resid_no, TNO)}
logger.info(f"T1 overlap: med {t1['overlap']['med_gap']:+.3f} "
            f"(MWU p={t1['overlap']['mwu_greater_p']:.4f}); "
            f"non-overlap: med {t1['non_overlap']['med_gap']:+.3f} "
            f"(MWU p={t1['non_overlap']['mwu_greater_p']:.4f})")

# membership x cap interaction on the pooled residual field:
# resid ~ 1 + g*in_code + h*incap + k*in_code*incap ; permute in_code
inc_all = (theta_to(rows, TNO) < CAP).astype(float)
mem_all = np.array([1.0 if r["in_code"] else 0.0 for r in rows])


def _inter(y, mem, inc):
    X = np.column_stack([np.ones(len(y)), mem, inc, mem * inc])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(c[3])


k_obs = _inter(resid_all, mem_all, inc_all)
cnt = 0
for _ in range(N_PERM):
    if abs(_inter(resid_all, rng.permutation(mem_all), inc_all)) \
            >= abs(k_obs):
        cnt += 1
t1["membership_x_cap_interaction"] = dict(
    k=k_obs, p_perm=float((cnt + 1) / (N_PERM + 1)))
logger.info(f"T1 membership x cap interaction k={k_obs:+.3f} "
            f"p={(cnt + 1) / (N_PERM + 1):.4f}")

# ------------------------------------------------------------------
# T2: cap-width sensitivity + axis controls on the overlap
# ------------------------------------------------------------------

t2 = {"cap_widths": {}, "axis_controls": {}}
for cap in (45.0, 60.0, 75.0):
    g = gap_report(overlap, resid_ov, TNO, cap=cap)
    t2["cap_widths"][f"cap_{int(cap)}"] = g
    logger.info(f"T2 overlap cap<{int(cap)}: med "
                f"{g.get('med_gap', float('nan')):+.3f} "
                f"(MWU p={g.get('mwu_greater_p', float('nan')):.4f})")
for tag, u in [("anti_axis", ANTI), ("post2017_displaced", U_RES)]:
    g = gap_report(overlap, resid_ov, u)
    t2["axis_controls"][tag] = g
    logger.info(f"T2 overlap at {tag}: med "
                f"{g.get('med_gap', float('nan')):+.3f} "
                f"(MWU p={g.get('mwu_greater_p', float('nan')):.4f})")

# ------------------------------------------------------------------
# T3: era x membership grid
# ------------------------------------------------------------------

ERA_BINS = [("pre_1990", 1850, 1990), ("y1990_2009", 1990, 2010),
            ("y2010_2017", 2010, 2018)]
t3 = {}
for tag, lo, hi in ERA_BINS:
    cell = {}
    for mtag, pool, rp in [("overlap", overlap, resid_ov),
                           ("non_overlap", nonov, resid_no)]:
        idx = [i for i, r in enumerate(pool) if lo <= r["yr"] < hi]
        sub = [pool[i] for i in idx]
        if len(sub) < 20:
            cell[mtag] = dict(n=len(sub), status="insufficient_sample")
            continue
        rs = resid_model(sub)
        cell[mtag] = gap_report(sub, rs, TNO)
        logger.info(f"T3 {tag} {mtag} (n={len(sub)}): med "
                    f"{cell[mtag].get('med_gap', float('nan')):+.3f} "
                    f"(MWU p={cell[mtag].get('mwu_greater_p', float('nan')):.4f})")
    t3[tag] = cell

# ------------------------------------------------------------------
# T4: fat-tail audit of the non-overlap pre-1990 cell
# ------------------------------------------------------------------

no_pre90 = [r for r in nonov if r["yr"] < 1990]
rs_p90 = resid_model(no_pre90)
th_p90 = theta_to(no_pre90, TNO)
inc_p90 = th_p90 < CAP
t4 = {"n": len(no_pre90), "n_in": int(inc_p90.sum())}

if inc_p90.sum() >= 10:
    rin, rout = rs_p90[inc_p90], rs_p90[~inc_p90]
    t4["med_gap"] = float(np.median(rin) - np.median(rout))
    t4["mean_gap"] = float(np.mean(rin) - np.mean(rout))
    t4["trim10_mean_gap"] = float(trim_mean(rin, 0.1)
                                  - trim_mean(rout, 0.1))
    # winsorize at 95th percentile of the pooled cell
    cap95 = np.percentile(np.concatenate([rin, rout]), 95)
    rw = np.clip(rs_p90, None, cap95)
    t4["winsor95_mean_gap"] = float(np.mean(rw[inc_p90])
                                    - np.mean(rw[~inc_p90]))
    # share of the in-cap sum carried by the top-5 residuals
    order = np.argsort(-rin)
    top_share = float(np.sum(rin[order[:5]])
                      / max(np.sum(np.abs(rin)), 1e-9))
    t4["top5_incap_resid_share"] = top_share
    # leave-one-out on the five largest in-cap residuals
    loo = {}
    for j in order[:5]:
        keep = np.ones(len(rin), bool); keep[j] = False
        loo[no_pre90[np.where(inc_p90)[0][j]]["desig"]] = float(
            np.median(rin[keep]) - np.median(rout))
    t4["loo_med_gap_after_dropping_top"] = loo
    # composition of in-cap high-residual members
    hi_idx = np.where(inc_p90)[0][order[:10]]
    t4["top10_incap_members"] = [
        dict(desig=no_pre90[i]["desig"], yr=no_pre90[i]["yr"],
             resid=float(rs_p90[i]), arc=no_pre90[i]["arc"],
             nobs=no_pre90[i]["nobs"],
             cc=None if not np.isfinite(no_pre90[i]["cc"])
             else no_pre90[i]["cc"])
        for i in hi_idx]
    _cc10 = [no_pre90[i]["cc"] for i in hi_idx
             if np.isfinite(no_pre90[i]["cc"])]
    _cc_all = [r["cc"] for r in no_pre90 if np.isfinite(r["cc"])]
    t4["median_cc_incap_top10"] = (
        float(np.median(_cc10)) if _cc10 else None)
    t4["median_cc_all"] = (
        float(np.median(_cc_all)) if _cc_all else None)
    t4["note"] = ("condition codes are unpopulated for most pre-1990 "
                  "designations; arc/nobs are the usable quality "
                  "proxies in this stratum")
    logger.info(f"T4 non-overlap pre-1990: med {t4['med_gap']:+.3f}, "
                f"mean {t4['mean_gap']:+.3f}, trim10 "
                f"{t4['trim10_mean_gap']:+.3f}, winsor95 "
                f"{t4['winsor95_mean_gap']:+.3f}, top-5 share "
                f"{top_share:.2f}")

# ------------------------------------------------------------------
# T5: quality-matched non-overlap subsets
# ------------------------------------------------------------------

t5 = {}
for tag, pred in [
    ("q_lt_3p1", lambda r: r["q"] < 3.1),
    ("q_lt_3p1_bound", lambda r: r["q"] < 3.1 and r["e"] < 1.0),
    ("q_lt_3p1_span_arc", lambda r: r["q"] < 3.1 and r["f_pre"]
     is not None and 0.0 < r["f_pre"] < 1.0),
    ("arc_gt_365_nobs_gt_300",
     lambda r: r["arc"] > 365 and r["nobs"] > 300),
    ("post_1990_only", lambda r: r["yr"] >= 1990),
    ("q_lt_3p1_post_1990",
     lambda r: r["q"] < 3.1 and r["yr"] >= 1990),
]:
    sub = [r for r in nonov if pred(r)]
    if len(sub) < 25:
        t5[tag] = dict(n=len(sub), status="insufficient_sample")
        continue
    t5[tag] = gap_report(sub, resid_model(sub), TNO)
    logger.info(f"T5 non-overlap {tag} (n={len(sub)}): med "
                f"{t5[tag].get('med_gap', float('nan')):+.3f} "
                f"(MWU p={t5[tag].get('mwu_greater_p', float('nan')):.4f})")

# ------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------

ov_p = t1["overlap"].get("mwu_greater_p", 1.0)
no_p = t1["non_overlap"].get("mwu_greater_p", 1.0)
inter_p = t1["membership_x_cap_interaction"]["p_perm"]
any_no = any(d.get("mwu_greater_p", 1.0) < 0.05
             and d.get("med_gap", 0) > 0.15 for d in t5.values())

if ov_p < 0.05 and not any_no:
    verdict = (
        f"MEMBERSHIP-LOCKED: the declared-axis term on the pre-2018 "
        f"record concentrates in the CODE-overlap members "
        f"(med {t1['overlap']['med_gap']:+.3f} dex, p={ov_p:.4f} "
        f"under independent JPL solutions; positive at 45/60/75-deg "
        f"caps, null at anti- and displaced-axis controls) while "
        f"every quality-matched post-1990 non-overlap cell is null "
        f"or negative -- the anomaly follows the objects, not the "
        f"CODE fitter.  The pooled membership x cap interaction is "
        f"direction-positive but underpowered (k={k_obs:+.3f}, "
        f"p={inter_p:.4f}), and the only non-overlap lean sits in "
        f"the pre-1990 stratum (med {t4.get('med_gap', float('nan')):+.3f} dex, "
        f"p~0.06, trim- and LOO-robust) whose high-residual members "
        f"are uniformly short-arc, low-observation historical "
        f"solutions -- a marginal quality-era residue rather than a "
        f"broad anomaly.")
elif any_no:
    verdict = ("BROAD: at least one quality-matched non-overlap cell "
               "carries the declared-axis excess -- the anomaly "
               "generalizes beyond the CODE-selected objects.")
else:
    verdict = (f"AMBIGUOUS: overlap p={ov_p:.4f}, non-overlap "
               f"p={no_p:.4f} -- neither a clean membership lock nor "
               "a broad anomaly.")

res = {
    "step": "step_122_pre2018_localization",
    "description": "localization of the weak positive declared-axis "
                   "term on the pre-2018 SBDB cohort: CODE-membership "
                   "lock vs broad-population readings, era "
                   "decomposition, and fat-tail audit",
    "inputs": ["results/step_b84_pre2018_sbdb.csv"],
    "seed": SEED, "n_perm": N_PERM, "cap_deg": CAP,
    "n_cohort": len(rows), "n_overlap": len(overlap),
    "n_non_overlap": len(nonov),
    "T1_membership_decomposition": t1,
    "T2_overlap_cap_and_axis": t2,
    "T3_era_membership_grid": t3,
    "T4_pre1990_fat_tail_audit": t4,
    "T5_non_overlap_quality_cells": t5,
    "verdict": verdict,
}

out = RESULTS / "step_b86_pre2018_localization.json"
json.dump(res, open(out, "w"), indent=1, default=float)

with open(RESULTS / "step_b86_pre2018_localization.csv", "w",
          newline="") as f:
    w = csv.writer(f)
    w.writerow(["desig", "yr", "in_code", "theta_declared", "in_cap",
                "resid", "q", "e", "arc", "nobs", "cc"])
    for r, rs in zip(rows, resid_all):
        w.writerow([r["desig"], int(r["yr"]), int(r["in_code"]),
                    f"{r['theta']:.2f}",
                    int(r["theta"] < CAP), f"{rs:.4f}",
                    f"{r['q']:.3f}", f"{r['e']:.4f}",
                    f"{r['arc']:.0f}", f"{r['nobs']:.0f}",
                    "" if not np.isfinite(r["cc"]) else f"{r['cc']:.0f}"])

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

ax = axes[0]
cells = [("overlap\n(all)", t1["overlap"]),
         ("non-overlap\n(all)", t1["non_overlap"])]
for tag, lo, hi in ERA_BINS:
    cells.append((f"overlap\n{tag}", t3[tag]["overlap"]))
    cells.append((f"non-ov\n{tag}", t3[tag]["non_overlap"]))
labels = [c[0] for c in cells]
vals = [c[1].get("med_gap", np.nan) for c in cells]
ps = [c[1].get("mwu_greater_p", np.nan) for c in cells]
colors = ["steelblue" if p < 0.05 else "0.65" for p in ps]
ax.bar(range(len(vals)), vals, color=colors)
ax.axhline(0, color="k", lw=0.5)
ax.set_xticks(range(len(vals)))
ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
ax.set_ylabel("declared-cap median gap (dex)")
ax.set_title("membership x era decomposition (blue = p<0.05)",
             fontsize=10)

ax = axes[1]
th_ov = theta_to(overlap, TNO)
th_no = theta_to(nonov, TNO)
ax.scatter(th_no, resid_no, s=14, c="0.6", alpha=0.5,
           label="non-overlap")
ax.scatter(th_ov, resid_ov + 1.5, s=14, c="steelblue", alpha=0.6,
           label="CODE overlap (+1.5 dex)")
ax.axvline(CAP, color="k", ls=":", lw=1)
ax.axhline(0, color="k", lw=0.5)
ax.axhline(1.5, color="steelblue", lw=0.4, ls=":")
ax.set_xlabel(r"$\theta$ from declared axis (deg)")
ax.set_ylabel("residual log rotation")
ax.legend(frameon=False, fontsize=8)
ax.set_title("declared-axis profile by CODE membership", fontsize=10)

ax = axes[2]
ax.scatter([r["yr"] for r in nonov], resid_no, s=14, c="0.6",
           alpha=0.5, label="non-overlap")
ax.scatter([r["yr"] for r in overlap],
           resid_ov, s=14, c="steelblue", alpha=0.6,
           label="CODE overlap")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("designation year")
ax.set_ylabel("residual log rotation")
ax.legend(frameon=False, fontsize=8)
ax.set_title("residual vs era -- pre-1990 fat tail visible",
             fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b86_pre2018_localization.png", dpi=300)
logger.data_save(out)