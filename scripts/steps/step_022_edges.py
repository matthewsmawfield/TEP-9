#!/usr/bin/env python3
"""
TEP / Planet-9 -- step 04: two real edges in the outer solar system
===================================================================

Under TEP a field domain boundary is literally a wall in the clock
structure: dynamics change character across it.  The solar system
contains two sharp, unexplained radial edges that are candidate
boundary signatures:

E1  The Kuiper cliff.  The classical Kuiper belt's radial density
    ends abruptly near 47-50 AU (Trujillo & Brown 2001; Allen,
    Bernstein & Malhotra 2001; Petit+ 2011).  No accepted
    explanation: a smooth truncation event is too soft, a hidden
    planet is post-hoc.  Here the semimajor-axis distribution of
    cold classical objects (i < 5 deg) is extracted from the full
    SBDB catalog and the edge location and sharpness are measured
    with a logistic-edge fit, separately for cold and hot classicals.

E2  The detachment boundary.  Objects with a > 150 AU separate into
    Neptune-coupled (q <~ 38 AU) and detached (q >~ 40 AU) sets; the
    intermediate perihelion range is anomalously empty -- the
    "perihelion gap" that marks where perturbative access ends.  A
    field boundary selects by crossing: the gap is the wall's
    inner-edge signature.  The q-distribution edge is measured the
    same way.

Also computed: the aphelion-distance distribution of the detached
population, locating the outer edge of the affected zone.

Outputs: results/step_04_edges.json,
         figures/supplementary/step_04_edges.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_022_edges")
tee_stdout(logger)
logger.header("Radial edge structure")

from pathlib import Path
import json
import numpy as np
from scipy.optimize import curve_fit

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260917)


def load_objects(fname):
    d = json.loads((DATA / fname).read_text())
    rows, fields = d["data"], d["fields"]
    objs = []
    for r in rows:
        o = dict(zip(fields, r))
        for k in ("a", "e", "i", "om", "w", "q", "ad", "H"):
            try:
                o[k] = float(o[k]) if o[k] not in (None, "") else np.nan
            except (TypeError, ValueError):
                o[k] = np.nan
        objs.append(o)
    return objs


def edge_fit(centers, counts):
    """Fit counts ~ A/(1+exp((x-x0)/w)) + B; return x0, w."""
    def model(x, A, x0, w, B):
        return A / (1 + np.exp((x - x0) / w)) + B
    p0 = [counts.max(), centers[np.argmax(np.abs(np.gradient(counts)))],
          1.0, 0.0]
    try:
        p, c = curve_fit(model, centers, counts, p0=p0,
                         bounds=([0, centers.min(), 0.01, 0],
                                 [np.inf, centers.max(), 10, np.inf]),
                         maxfev=20000)
        err = np.sqrt(np.diag(c))
        return {"A": float(p[0]), "edge_AU": float(p[1]),
                "width_AU": float(p[2]), "B": float(p[3]),
                "edge_err": float(err[1]), "width_err": float(err[2]),
                "converged": True}
    except Exception as e:
        return {"converged": False, "err": str(e)}


def main():
    objs = load_objects("sbdb_outer_ss.json")
    out = {}

    # ---------------- E1 Kuiper cliff ----------------
    for label, i_max in [("cold_i<5", 5), ("hot_i>=5", np.inf)]:
        if label.startswith("cold"):
            pop = [o for o in objs if np.isfinite(o["a"])
                   and 38 < o["a"] < 62 and np.isfinite(o["i"])
                   and o["i"] < 5 and np.isfinite(o["q"]) and o["q"] > 36]
        else:
            pop = [o for o in objs if np.isfinite(o["a"])
                   and 38 < o["a"] < 62 and np.isfinite(o["i"])
                   and o["i"] >= 5 and np.isfinite(o["q"]) and o["q"] > 36]
        a = np.array([o["a"] for o in pop])
        edges = np.arange(38, 62.5, 0.5)
        h, _ = np.histogram(a, bins=edges)
        c = 0.5 * (edges[:-1] + edges[1:])
        fit = edge_fit(c, h)
        out[label] = {"N": len(a), "fit": fit,
                      "bin_centers": c.tolist(),
                      "counts": h.tolist()}
        if fit.get("converged"):
            print(f"E1 {label}: N={len(a)}, edge={fit['edge_AU']:.2f}"
                  f"+-{fit['edge_err']:.2f} AU, "
                  f"width={fit['width_AU']:.2f}+-{fit['width_err']:.2f} AU")
        else:
            print(f"E1 {label}: N={len(a)}, fit failed: {fit}")

    # cold-vs-hot edge ratio: does the hot belt continue past the cliff?
    # (published result: hot objects are found past 50 AU, cold are not)

    # ---------------- E2 detachment structure ----------------
    # The q distribution for a>150 AU is bimodal: a detached-belt
    # pile at 34-45 AU, a "perihelion desert" near 48-62 AU, then the
    # sednoid group (Sedna, VP113, TG387, KG163) at 64-82 AU.
    det = [o for o in objs if np.isfinite(o["a"]) and o["a"] > 150
           and np.isfinite(o["q"]) and o["q"] > 28
           and o["condition_code"] is not None
           and int(o["condition_code"]) <= 3]
    q = np.array([o["q"] for o in det])
    edges = np.arange(30, 90, 2.0)
    h, _ = np.histogram(q, bins=edges)
    c = 0.5 * (edges[:-1] + edges[1:])
    out["q_histogram"] = {"N": len(q), "bin_centers": c.tolist(),
                          "counts": h.tolist()}

    # inner edge: non-parametric lower boundary of the detached
    # population -- 10th percentile of q with bootstrap 68% interval
    boot = np.percentile(rng.choice(q, (10000, len(q))), 10, axis=1)
    q10 = float(np.percentile(q, 10))
    lo, hi = np.percentile(boot, [16, 84])
    out["detachment_inner_edge"] = {
        "q10_AU": round(q10, 2),
        "ci68": [round(float(lo), 2), round(float(hi), 2)],
        "min_q_AU": round(float(q.min()), 2)}
    print(f"E2 inner edge: q10={q10:.1f} AU "
          f"CI68 [{lo:.1f},{hi:.1f}], min q={q.min():.1f}")

    # perihelion-desert deficit: observed count in 48-62 AU vs the
    # expectation from a KDE-smoothed version of the observed
    # distribution itself (tests whether the dip is deeper than any
    # smooth density allows)
    from scipy.stats import gaussian_kde, poisson
    kde = gaussian_kde(q, bw_method=0.25)
    xs = np.linspace(30, 90, 6000)
    dens = kde(xs)
    dens /= np.trapezoid(dens, xs)
    p_gap_range = float(np.trapezoid(dens[(xs >= 48) & (xs < 62)],
                                     xs[(xs >= 48) & (xs < 62)]))
    lam = len(q) * p_gap_range
    n_obs = int(((q >= 48) & (q < 62)).sum())
    p_def = float(poisson.cdf(n_obs, lam))
    out["perihelion_desert"] = {
        "range_AU": [48, 62], "n_observed": n_obs,
        "expected_smooth": round(float(lam), 2),
        "poisson_p_deficit": p_def,
        "objects_in_desert": [o["full_name"].strip() for o in det
                              if 48 <= o["q"] < 62]}
    print(f"E2 perihelion desert 48-62 AU: observed {n_obs}, "
          f"smooth-expected {lam:.1f}, p={p_def:.2e}")

    # sednoid group census
    sed = [o for o in det if o["q"] > 62]
    out["sednoid_group"] = [{"name": o["full_name"].strip(),
                             "q": o["q"], "a": o["a"]} for o in sed]
    print(f"E2 sednoid group q>62: {[o['full_name'].strip() for o in sed]}")

    # aphelion distribution of detached population (outer edge)
    ad = np.array([o["ad"] for o in det])
    out["detached_aphelion"] = {"median_AU": float(np.median(ad)),
                                "p16": float(np.percentile(ad, 16)),
                                "p84": float(np.percentile(ad, 84)),
                                "max_AU": float(ad.max())}
    print(f"detached aphelia: median {np.median(ad):.0f} AU, "
          f"[{np.percentile(ad,16):.0f},{np.percentile(ad,84):.0f}]")

    RES.mkdir(exist_ok=True)
    (RES / "step_04_edges.json").write_text(json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(12.5, 4))
    for label, color in [("cold_i<5", "crimson"), ("hot_i>=5", "steelblue")]:
        d = out[label]
        ax[0].step(d["bin_centers"], d["counts"], where="mid",
                   color=color, label=f"{label} (N={d['N']})")
        if d["fit"].get("converged"):
            x = np.linspace(38, 62, 300)
            ax[0].plot(x, d["fit"]["A"] /
                       (1 + np.exp((x - d["fit"]["edge_AU"]) /
                                   d["fit"]["width_AU"])) + d["fit"]["B"],
                       color=color, ls=":", lw=1)
            ax[0].axvline(d["fit"]["edge_AU"], color=color, ls=":",
                          lw=0.8)
    ax[0].set_xlabel("a (AU)"); ax[0].set_ylabel("N per 0.5 AU")
    ax[0].legend(fontsize=8)
    ax[0].set_title("E1: Kuiper cliff, classical belt", fontsize=9)

    d = out["q_histogram"]
    ax[1].step(d["bin_centers"], d["counts"], where="mid",
               color="steelblue")
    ax[1].axvspan(48, 62, color="r", alpha=0.12,
                  label="perihelion desert")
    ax[1].axvspan(62, 90, color="g", alpha=0.08,
                  label="sednoid group")
    ie = out["detachment_inner_edge"]
    ax[1].axvline(ie["q10_AU"], color="k", ls=":", lw=0.8,
                  label=f"q10 {ie['q10_AU']:.1f} AU")
    ax[1].legend(fontsize=8)
    ax[1].set_xlabel("q (AU)"); ax[1].set_ylabel("N per 2 AU")
    ax[1].set_title("E2: detachment structure, a>150 AU", fontsize=9)

    ad = np.array([o["ad"] for o in det])
    ax[2].hist(ad, bins=30, range=(150, 3000), color="steelblue",
               edgecolor="k", lw=0.4)
    ax[2].set_xlabel("aphelion (AU)"); ax[2].set_ylabel("N")
    ax[2].set_title("detached-population aphelia", fontsize=9)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_04_edges.png", dpi=300)
    logger.data_save(RESULTS / "step_04_edges.json")
    logger.data_save(RESULTS / "figures/supplementary/step_04_edges.png")


if __name__ == "__main__":
    main()
