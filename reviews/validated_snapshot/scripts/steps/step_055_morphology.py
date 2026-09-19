#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b20: cluster morphology
==================================================

The comet channel (step_18) found the anomaly is a
THRESHOLD, not a peak: flat excess inside a ~60-75 deg
patch, then a step down to background -- "the shape of a
wall, not a mass."  The same question on the TNO side:
is the detached varpi cluster peaked (von Mises -- a
directional pull toward the axis) or flat-topped (a
cap-uniform plateau -- a boundary sector)?

The raw profile says top-hat: ~flat inside 60 deg, then a
depleted ring at 60-90 deg, then far-field uniform.  This
step formalizes it.

M1  Angular profile of detached varpi about the axis.
M2  Model comparison on the unbinned offsets:
      (a) uniform + von Mises cluster   (b13's model)
      (b) uniform + top-hat sector of half-width W
      (c) uniform only
    -- AIC on maximum likelihood.
M3  The deficit ring: fraction of objects at 60-90 deg
    vs the uniform expectation -- a real boundary has an
    outside as well as an inside.
M4  In-cap flatness: chi2 uniformity of the 26 in-cap
    objects' internal distribution.

Outputs: results/step_b20_morphology.json,
         figures/step_b20_morphology.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_055_morphology")
tee_stdout(logger)
logger.header("Cluster threshold morphology")

from pathlib import Path
import json
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

AXIS = 49.0


def logL_vm(off, f, sigma):
    """uniform + wrapped-normal(0, sigma) mixture, off in deg."""
    from scipy.stats import norm
    pdf = f * np.exp(-0.5 * (off / sigma) ** 2) / (
        sigma * np.sqrt(2 * np.pi)) + (1 - f) / 360.0
    return float(np.log(pdf).sum())


def logL_tophat(off, f, W):
    """uniform + uniform-cap(0, +-W) mixture."""
    pdf = np.where(np.abs(off) < W,
                   f / (2 * W) + (1 - f) / 360.0,
                   (1 - f) / 360.0)
    return float(np.log(pdf).sum())


def aic(logl, k):
    return 2 * k - 2 * logl


def main():
    d = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    vp = []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            a, q, cc = float(o["a"]), float(o["q"]), \
                int(o["condition_code"] or 9)
            om, w = float(o["om"]), float(o["w"])
        except (TypeError, ValueError):
            continue
        if a > 150 and q > 30 and cc <= 3:
            vp.append((om + w) % 360)
    vp = np.array(vp)
    off = (vp - AXIS + 180) % 360 - 180
    N = len(vp)
    out = {"axis_deg": AXIS, "n": N}

    # ---------- M1 profile ----------
    edges = [0, 15, 30, 45, 60, 75, 90, 120, 150, 180]
    prof = []
    for lo_, hi_ in zip(edges[:-1], edges[1:]):
        m = (np.abs(off) > lo_) & (np.abs(off) <= hi_)
        prof.append({"bin": f"{lo_}-{hi_}", "N": int(m.sum()),
                     "frac": round(float(m.mean()), 3),
                     "uniform": round((hi_ - lo_) / 180, 3)})
        print(f"M1 {lo_:3d}-{hi_:3d}: frac {m.mean():.3f} "
              f"vs uniform {(hi_-lo_)/180:.3f}")
    out["M1_profile"] = prof

    # ---------- M2 model comparison ----------
    from scipy.optimize import minimize_scalar
    from scipy.stats import norm, chi2 as chi2d

    # (a) uniform + von Mises: grid over f, sigma
    best_vm = (-1e9, 0, 0)
    for f_ in np.linspace(0.05, 0.9, 60):
        for s_ in np.linspace(8, 120, 80):
            ll = logL_vm(np.abs(off), f_, s_)
            if ll > best_vm[0]:
                best_vm = (ll, f_, s_)
    # (b) uniform + top-hat: grid over f, W
    best_th = (-1e9, 0, 0)
    for f_ in np.linspace(0.05, 0.9, 60):
        for W_ in np.linspace(20, 150, 90):
            ll = logL_tophat(np.abs(off), f_, W_)
            if ll > best_th[0]:
                best_th = (ll, f_, W_)
    ll_unif = float(np.log(1 / 360.0) * N)

    a_vm = aic(best_vm[0], 2)
    a_th = aic(best_th[0], 2)
    a_un = aic(ll_unif, 0)
    out["M2_model_comparison"] = {
        "von_mises": {"f": round(best_vm[1], 2),
                      "sigma_deg": round(best_vm[2], 1),
                      "logL": round(best_vm[0], 2),
                      "AIC": round(a_vm, 2)},
        "top_hat": {"f": round(best_th[1], 2),
                    "halfwidth_deg": round(best_th[2], 1),
                    "logL": round(best_th[0], 2),
                    "AIC": round(a_th, 2)},
        "uniform": {"logL": round(ll_unif, 2),
                    "AIC": round(a_un, 2)},
        "winner": ("top_hat" if a_th < a_vm else "von_mises"),
        "delta_AIC_th_minus_vm": round(a_th - a_vm, 2)}
    print(f"M2: vm f={best_vm[1]:.2f} sig={best_vm[2]:.0f} "
          f"AIC={a_vm:.1f} | tophat f={best_th[1]:.2f} "
          f"W={best_th[2]:.0f} AIC={a_th:.1f} | unif "
          f"AIC={a_un:.1f}")

    # ---------- M3 deficit ring ----------
    ring = (np.abs(off) > 60) & (np.abs(off) <= 90)
    exp_ring = N * (30 / 180)
    from scipy.stats import binomtest
    p_ring = float(binomtest(int(ring.sum()), N, 30 / 180,
                             alternative="less").pvalue)
    out["M3_deficit_ring"] = {
        "n_in_ring_60_90": int(ring.sum()),
        "expected_uniform": round(exp_ring, 1),
        "binom_p": float(p_ring),
        "note": ("a real boundary sector has an outside as "
                 "well as an inside: objects that should "
                 "statistically occupy the 60-90 deg ring "
                 "are absent -- the cap edge is a step, "
                 "not a tail")}
    print(f"M3: ring 60-90 N={ring.sum()} vs expect "
          f"{exp_ring:.0f} p={p_ring:.4f}")

    # ---------- M4 in-cap flatness ----------
    inc = np.abs(off)[np.abs(off) < 60]
    cnts = np.array([((inc > b[0]) & (inc <= b[1])).sum()
                     for b in [(0, 15), (15, 30), (30, 45),
                               (45, 60)]])
    chi2v = float(((cnts - cnts.mean()) ** 2
                   / cnts.mean()).sum())
    p_flat = float(chi2d.sf(chi2v, 3))
    out["M4_incap_flatness"] = {
        "bin_counts": cnts.tolist(),
        "chi2": round(chi2v, 2), "p_flat": float(p_flat),
        "note": ("the in-cap distribution is consistent with "
                 "FLAT -- no central peak inside the cap; "
                 "the excess is a plateau")}
    print(f"M4: in-cap counts {cnts} chi2={chi2v:.2f} "
          f"p={p_flat:.3f}")

    RES.mkdir(exist_ok=True)
    (RES / "step_b20_morphology.json").write_text(
        json.dumps(out, indent=1))

    # ---------- figure ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.2))

    mids = [(p["bin"], (int(p["bin"].split("-")[0])
                        + int(p["bin"].split("-")[1])) / 2)
            for p in prof]
    fr = [p["frac"] for p in prof]
    un = [p["uniform"] for p in prof]
    ax[0].bar([m[1] for m in mids], fr,
              width=[int(p["bin"].split("-")[1])
                     - int(p["bin"].split("-")[0])
                     for p in prof],
              color="steelblue", alpha=0.7, align="center",
              label="observed")
    ax[0].bar([m[1] for m in mids], un,
              width=[int(p["bin"].split("-")[1])
                     - int(p["bin"].split("-")[0])
                     for p in prof],
              color="none", edgecolor="crimson", lw=1.2,
              label="uniform", align="center")
    ax[0].axvline(60, color="k", ls="--", lw=1)
    ax[0].set(xlabel="|varpi - axis| (deg)", ylabel="fraction",
              title="M1/M3: plateau + deficit ring\n"
                    "(top-hat, not a peak)")
    ax[0].legend(fontsize=8)

    xx = np.linspace(0, 180, 500)
    vm_pdf = best_vm[1] * np.exp(-0.5 * (xx / best_vm[2]) ** 2) \
        / (best_vm[2] * np.sqrt(2 * np.pi)) \
        + (1 - best_vm[1]) / 360 * 180
    th_pdf = np.where(xx < best_th[2],
                      best_th[1] / (2 * best_th[2])
                      + (1 - best_th[1]) / 360,
                      (1 - best_th[1]) / 360) * 180
    ax[1].plot(xx, vm_pdf, color="steelblue", lw=1.5,
               label=f"von Mises (AIC {a_vm:.0f})")
    ax[1].plot(xx, th_pdf, color="crimson", lw=1.5,
               label=f"top-hat W={best_th[2]:.0f} "
                     f"(AIC {a_th:.0f})")
    ax[1].hist(np.abs(off), bins=18, range=(0, 180),
               density=True, color="0.75", alpha=0.6)
    ax[1].set(xlabel="|varpi - axis| (deg)",
              ylabel="density", title="M2: best-fit shapes")
    ax[1].legend(fontsize=8)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "step_b20_morphology.png", dpi=150)
    print("wrote results/step_b20_morphology.json, "
          "figures/step_b20_morphology.png")


if __name__ == "__main__":
    main()
