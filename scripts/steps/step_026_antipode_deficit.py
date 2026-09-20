#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b6: the anti-cluster side of the sky
=============================================================

Every published discussion asks only whether the cluster side is
excess.  The discriminating question is what the OPPOSITE side
looks like, because the two candidate mechanisms predict
different tails:

  * point-mass shepherding: confined orbits pile onto the
    perturber-facing direction; the rest of the circle should be
    ~uniform background -- a bump, not a dip.

  * fixed field boundary (TEP): perihelia are phase-locked by a
    directed structure; objects pointing into the far side are
    systematically disfavoured -- a bump AND a compensating
    deficit, i.e. the whole distribution is repelled from the
    anti-axis, not just attracted to the axis.

Inspection of the detached sample suggests the anti-sector is
sparse: of 44 objects, only ~4-5 sit at varpi in 200-290 deg.

C1  anti-sector deficit: count within +-60 deg of the anti-axis
    (varpi ~ 229 deg) vs the conditioned-null distribution of
    the same count -- lower-tail p.

C2  shape discrimination: fit the varpi distribution with
    (a) uniform + von-Mises bump at the cluster axis,
    (b) uniform + VM bump + VM dip (negative component) at the
    anti-axis,
    (c) a single wrapped "repulsion" model p(th) ~ 1 - a*cos(
    th - mu_anti).
    Compare log-likelihoods / AIC -- is a dip component
    required?

C3  who lives on the far side: compare a, q, i, discovery year,
    and orbit quality of the anti-sector objects vs the
    clustered population -- is the anti-cluster a different
    population (e.g. low-q plungers or recent discoveries)?

Outputs: results/step_b6_antipode_deficit.json,
         figures/supplementary/step_b6_antipode_deficit.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_026_antipode_deficit")
tee_stdout(logger)
logger.header("Anti-axis deficit")

from pathlib import Path
import json
import re
import numpy as np
from scipy.stats import vonmises, kstest
from scipy.optimize import minimize

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260919)
N_MC = 20000

AXIS = 49.1          # cluster axis, deg
ANTI = (AXIS + 180) % 360   # 229.1


def load_sbdb(fname="sbdb_outer_ss.json"):
    d = json.loads((DATA / fname).read_text())
    objs = []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            for k in ("a", "e", "i", "om", "w", "q"):
                o[k] = float(o[k])
            o["cc"] = int(o["condition_code"]) \
                if o["condition_code"] is not None else 9
            o["arc"] = float(o["data_arc"]) if o["data_arc"] else 0
        except (TypeError, ValueError):
            continue
        o["name"] = str(o["full_name"]).strip()
        try:
            o["first_year"] = float(str(o["first_obs"])[:4])
        except (TypeError, ValueError):
            o["first_year"] = np.nan
        objs.append(o)
    return objs


def circ_R(a):
    return abs(np.exp(1j * np.asarray(a)).mean())


def dphi(a, b):
    return (a - b + np.pi) % (2 * np.pi) - np.pi


def loglik(theta, pdf):
    p = pdf(theta)
    p = np.clip(p, 1e-12, None)
    return float(np.log(p).sum())


def main():
    objs = load_sbdb()
    det = [o for o in objs if o["a"] > 150 and o["q"] > 30
           and o["cc"] <= 3]
    n = len(det)
    om = np.deg2rad([o["om"] for o in det])
    w = np.deg2rad([o["w"] for o in det])
    vp = (om + w) % (2 * np.pi)
    i_rad = np.deg2rad([o["i"] for o in det])

    ax_dir = np.deg2rad(AXIS)
    anti = np.deg2rad(ANTI)

    out = {"sample": "a>150, q>30, cc<=3", "N": n,
           "axis_deg": AXIS, "anti_deg": ANTI}

    # ---------------- C1 anti-sector deficit ----------------
    cap = np.abs(dphi(vp, anti)) < np.pi / 3     # +-60 deg
    n_anti = int(cap.sum())
    # conditioned null: fix i, uniform om,w -> distribution of the
    # anti-cap count
    om_mc = rng.uniform(0, 2 * np.pi, (N_MC, n))
    w_mc = rng.uniform(0, 2 * np.pi, (N_MC, n))
    vp_mc = (om_mc + w_mc) % (2 * np.pi)
    n_mc = (np.abs(dphi(vp_mc, anti)) < np.pi / 3).sum(axis=1)
    p_low = float((int((n_mc <= n_anti).sum()) + 1) / (N_MC + 1))
    # symmetric check: the cluster-side count for context
    cap_c = np.abs(dphi(vp, ax_dir)) < np.pi / 3
    n_cl = int(cap_c.sum())
    n_mc_c = (np.abs(dphi(vp_mc, ax_dir)) < np.pi / 3).sum(axis=1)
    p_high = float((int((n_mc_c >= n_cl).sum()) + 1)
                   / (N_MC + 1))
    out["C1_anti_deficit"] = {
        "n_in_anti_cap": n_anti,
        "null_median": round(float(np.median(n_mc)), 1),
        "null_p10": round(float(np.percentile(n_mc, 10)), 1),
        "p_lower": p_low,
        "n_in_cluster_cap": n_cl,
        "p_upper_cluster": p_high,
        "note": "lower-tail test of the count within +-60 deg of "
                "the anti-axis vs conditioned null"}
    print(f"C1 anti-cap N={n_anti} (null median "
          f"{np.median(n_mc):.0f}) p_low={p_low:.4f} | "
          f"cluster-cap N={n_cl} p_high={p_high:.4f}")

    # ---------------- C2 shape discrimination ----------------
    def pdf_mix_cluster(th, f, kap):
        return (1 - f) / (2 * np.pi) + f * vonmises.pdf(
            th, kap, loc=ax_dir)

    def pdf_mix_bump_dip(th, f1, k1, f2, k2):
        return ((1 - f1 - f2) / (2 * np.pi)
                + f1 * vonmises.pdf(th, k1, loc=ax_dir)
                - f2 * vonmises.pdf(th, k2, loc=anti))

    def pdf_repel(th, a, mu):
        # first-harmonic repulsion from mu (the anti-axis)
        return (1 - a * np.cos(th - mu)) / (2 * np.pi)

    # fit (a) uniform+VM bump
    def neg_a(par):
        f, k = par
        if not (0 <= f <= 1) or k <= 0:
            return 1e9
        return -loglik(vp, lambda t: pdf_mix_cluster(t, f, k))

    ra = minimize(neg_a, x0=[0.4, 1.0], method="Nelder-Mead")
    fa, ka = ra.x
    ll_a = -ra.fun

    # fit (b) uniform + bump + dip
    def neg_b(par):
        f1, k1, f2, k2 = par
        if not (0 <= f1) or not (0 <= f2) or f1 + f2 > 0.9 \
                or k1 <= 0 or k2 <= 0:
            return 1e9
        pdf = lambda t: pdf_mix_bump_dip(t, f1, k1, f2, k2)
        if np.any(pdf(t_g) < 0):
            return 1e9
        return -loglik(vp, pdf)

    t_g = np.linspace(0, 2 * np.pi, 721)
    rb = minimize(neg_b, x0=[0.4, 1.0, 0.15, 1.0],
                  method="Nelder-Mead",
                  options={"maxiter": 20000})
    ll_b = -rb.fun

    # fit (c) repulsion p ~ 1 - a cos(th - mu_anti), a in [0,1]
    def neg_c(par):
        a = par[0]
        if not (0 <= a <= 1):
            return 1e9
        return -loglik(vp, lambda t: pdf_repel(t, a, anti))

    rc = minimize(neg_c, x0=[0.3], method="Nelder-Mead")
    ac = rc.x[0]
    ll_c = -rc.fun
    # repulsion-direction free fit too
    def neg_c2(par):
        a, mu = par
        if not (0 <= a <= 1):
            return 1e9
        return -loglik(vp, lambda t: pdf_repel(t, a, mu))

    rc2 = minimize(neg_c2, x0=[0.3, anti], method="Nelder-Mead")
    ll_c2 = -rc2.fun
    mu_c2 = float(rc2.x[1]) % (2 * np.pi)

    # AIC (k params: a=2, b=4, c=1, c2=2)
    aic = {"uniform_plus_bump": 2 * 2 - 2 * ll_a,
           "uniform_plus_bump_plus_dip": 2 * 4 - 2 * ll_b,
           "repulsion_fixed_axis": 2 * 1 - 2 * ll_c,
           "repulsion_free_axis": 2 * 2 - 2 * ll_c2}
    out["C2_shape"] = {
        "models": {
            "uniform_plus_bump": {
                "f": round(float(fa), 3),
                "kappa": round(float(ka), 2),
                "loglik": round(ll_a, 2)},
            "uniform_plus_bump_plus_dip": {
                "f_bump": round(float(rb.x[0]), 3),
                "k_bump": round(float(rb.x[1]), 2),
                "f_dip": round(float(rb.x[2]), 3),
                "k_dip": round(float(rb.x[3]), 2),
                "loglik": round(ll_b, 2)},
            "repulsion_fixed_axis": {
                "a": round(float(ac), 3),
                "loglik": round(ll_c, 2)},
            "repulsion_free_axis": {
                "a": round(float(rc2.x[0]), 3),
                "mu_deg": round(float(np.rad2deg(mu_c2)), 1),
                "loglik": round(ll_c2, 2)}},
        "AIC": {k: round(v, 2) for k, v in aic.items()},
        "best_AIC": min(aic, key=aic.get),
        "note": "AIC prefers the model that needs a dip/repulsion "
                "component opposite the axis, if any"}
    print(f"C2 AIC: " + ", ".join(f"{k}={v:.1f}"
                                  for k, v in aic.items())
          + f" | bump+dip params {np.round(rb.x, 2)}")

    # ---------------- C3 who lives on the far side ----------------
    far = [o for o, c in zip(det, cap) if c]
    near = [o for o, c in zip(det, cap_c) if c]
    def summ(s, k):
        v = [o[k] for o in s]
        return {"N": len(v),
                "median": round(float(np.median(v)), 1),
                "lo": round(float(np.min(v)), 1),
                "hi": round(float(np.max(v)), 1)}
    out["C3_far_side"] = {
        "anti_cap_objects": [{"name": o["name"],
                              "a": o["a"], "q": o["q"],
                              "i": o["i"],
                              "varpi": round((o["om"] + o["w"])
                                             % 360, 1),
                              "first_year": o["first_year"],
                              "cc": o["cc"],
                              "arc": o["arc"]}
                             for o in far],
        "compare": {
            "a": {"anti": summ(far, "a"),
                  "cluster": summ(near, "a")},
            "i": {"anti": summ(far, "i"),
                  "cluster": summ(near, "i")},
            "q": {"anti": summ(far, "q"),
                  "cluster": summ(near, "q")},
            "first_year": {
                "anti": summ(far, "first_year"),
                "cluster": summ(near, "first_year")}},
        "note": "are anti-axis objects a different population?"}
    print(f"C3 anti-cap objects: "
          f"{[(o['name'], round((o['om']+o['w'])%360)) for o in far]}")

    RES.mkdir(exist_ok=True)
    (RES / "step_b6_antipode_deficit.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.4))

    th = np.linspace(0, 360, 361)
    ax[0].hist(np.rad2deg(vp), bins=np.arange(0, 361, 20),
               density=True, color="steelblue", alpha=0.8,
               edgecolor="k", lw=0.4, label="detached varpi")
    tt = np.deg2rad(th)
    ax[0].plot(th, pdf_mix_cluster(tt, fa, ka), "r-", lw=1.6,
               label=f"U+VM bump (AIC {aic['uniform_plus_bump']:.0f})")
    ax[0].plot(th, pdf_mix_bump_dip(tt, *rb.x), "g-", lw=1.6,
               label=f"U+bump+dip (AIC "
                     f"{aic['uniform_plus_bump_plus_dip']:.0f})")
    ax[0].plot(th, pdf_repel(tt, ac, anti), "m--", lw=1.4,
               label=f"repulsion (AIC "
                     f"{aic['repulsion_fixed_axis']:.0f})")
    ax[0].axvline(AXIS, color="r", ls=":", lw=1)
    ax[0].axvline(ANTI, color="k", ls=":", lw=1)
    ax[0].set(xlabel="varpi (deg)", ylabel="density",
              title="C2: bump vs dip shape")
    ax[0].legend(fontsize=6)

    ax[1].hist(n_mc, bins=np.arange(-0.5, 20), density=True,
               color="0.7", edgecolor="k", lw=0.3,
               label="null anti-cap count")
    ax[1].axvline(n_anti, color="crimson", lw=1.8,
                  label=f"observed = {n_anti}")
    ax[1].set(xlabel="N within +-60 deg of anti-axis",
              ylabel="density",
              title=f"C1: anti-sector deficit\np_low={p_low:.3f}")
    ax[1].legend(fontsize=8)

    far_v = np.rad2deg([o["varpi"] for o in
                        out["C3_far_side"]["anti_cap_objects"]])
    ax[2].scatter([o["a"] for o in det],
                  [(o["om"] + o["w"]) % 360 for o in det],
                  s=30, c="steelblue", edgecolor="k", lw=0.4)
    ax[2].axhspan(ANTI - 60, ANTI + 60, color="0.85", alpha=0.5,
                  label="anti-sector")
    ax[2].axhline(AXIS, color="r", ls="--", lw=1,
                  label="cluster axis")
    ax[2].set_xscale("log")
    ax[2].set(xlabel="a (AU)", ylabel="varpi (deg)",
              title="C3: who lives on the far side")
    ax[2].legend(fontsize=7)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b6_antipode_deficit.png", dpi=300)
    logger.data_save(RESULTS / "step_b6_antipode_deficit.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b6_antipode_deficit.png")


if __name__ == "__main__":
    main()
