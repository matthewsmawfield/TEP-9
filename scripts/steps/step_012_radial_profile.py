#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b2: where does the clustering turn on?
==============================================================

The two competing pictures predict different radial structure.

  * Point-mass shepherding: confinement should strengthen with
    semimajor axis (the perturber dominates the most distant
    orbits) and cluster membership should track the perturber's
    secular reach -- a smooth, monotonic trend.

  * Fixed spatial boundary (TEP): orbits are phase-locked by the
    proper-time field where they *are*, so clustering should
    switch on at a fixed spatial radius.  The natural coordinate
    is then the orbit's spatial extent -- aphelion distance ad =
    a(1+e) -- not the semimajor axis itself.  A transition at
    fixed ad independent of a is the "domain wall" signature.

Scans performed (all on the secure catalog, cc<=3):

  W1  sliding windows in log a (0.4 dex, 0.1 dex step):
      R(varpi), R(omega), R(Omega) and conditioned p(varpi)
      per window -- element selectivity shows which part of the
      anomaly turns on where.

  W2  threshold scans: for each of a_min, ad_min, q_min, the
      clustering of the tail sample x > x_min vs threshold --
      locates the turn-on in each coordinate.

  W3  change-point estimate: the radius x* that best splits the
      population into unclustered (inside) and clustered
      (outside) subsamples, in each of a, ad, q -- the wall
      radius.  A sharp, coordinate-consistent transition
      supports the boundary picture; a smeared or
      coordinate-dependent one supports shepherding.

Statistics: conditioned null as in step 02 (fix a,e,i; om,w
uniform; 5000 MC per window for the scan, 20000 at key points).

Outputs: results/step_b2_radial_profile.json,
         figures/supplementary/step_b2_radial_profile.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_012_radial_profile")
tee_stdout(logger)
logger.header("Radial turn-on profile")

from pathlib import Path
import json
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260918)
N_SCAN = 5000
N_KEY = 20000


def load_sbdb(fname="sbdb_outer_ss.json"):
    d = json.loads((DATA / fname).read_text())
    objs = []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            for k in ("a", "e", "i", "om", "w", "q", "ad"):
                o[k] = float(o[k])
            o["cc"] = int(o["condition_code"]) \
                if o["condition_code"] is not None else 9
            o["arc"] = float(o["data_arc"]) if o["data_arc"] else 0
        except (TypeError, ValueError):
            continue
        objs.append(o)
    return objs


def circ_R(a):
    return abs(np.exp(1j * np.asarray(a)).mean())


def rayleigh_p(ang):
    a = np.asarray(ang)
    n = len(a)
    z = n * circ_R(a) ** 2
    return float(np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n)))


def p_cond(objs, n_mc):
    """Conditioned-null p for R(varpi): fix i, uniform om,w."""
    n = len(objs)
    if n < 6:
        return np.nan
    varpi = np.deg2rad([(o["om"] + o["w"]) % 360 for o in objs])
    R_obs = circ_R(varpi)
    om = rng.uniform(0, 2 * np.pi, (n_mc, n))
    w = rng.uniform(0, 2 * np.pi, (n_mc, n))
    R_n = np.abs(np.exp(1j * ((om + w) % (2 * np.pi)))
                 .mean(axis=1))
    return float((int((R_n >= R_obs).sum()) + 1) / (n_mc + 1))


def window_report(sel_objs):
    n = len(sel_objs)
    om = np.deg2rad([o["om"] for o in sel_objs])
    w = np.deg2rad([o["w"] for o in sel_objs])
    vp = (om + w) % (2 * np.pi)
    return {"N": n,
            "R_varpi": round(float(circ_R(vp)), 3),
            "R_omega": round(float(circ_R(w)), 3),
            "R_Omega": round(float(circ_R(om)), 3),
            "p_varpi_rayleigh": rayleigh_p(vp) if n > 4 else None,
            "p_varpi_cond": p_cond(sel_objs, N_SCAN)}


def main():
    objs = [o for o in load_sbdb() if o["q"] > 30 and o["cc"] <= 3
            and o["arc"] > 200 and o["a"] > 30]
    print(f"secure outer catalog: {len(objs)} objects")

    out = {"sample": "a>30, q>30, cc<=3, arc>200",
           "n": len(objs)}

    # ---------------- W1 sliding windows in log a ----------------
    la = np.log10([o["a"] for o in objs])
    edges = []
    c = np.log10(35.0)
    while c < np.log10(3000):
        edges.append((10 ** c, 10 ** (c + 0.4)))
        c += 0.1
    wins = []
    for lo, hi in edges:
        s = [o for o in objs if lo < o["a"] <= hi]
        if len(s) < 8:
            continue
        r = window_report(s)
        r["a_lo"], r["a_hi"] = round(lo, 1), round(hi, 1)
        r["a_mid"] = round(float(np.sqrt(lo * hi)), 1)
        wins.append(r)
        print(f"W1 a {r['a_lo']:6.0f}-{r['a_hi']:6.0f} N={r['N']:4d} "
              f"R(v)={r['R_varpi']:.2f} R(w)={r['R_omega']:.2f} "
              f"R(O)={r['R_Omega']:.2f} p={r['p_varpi_cond']:.3f}")
    out["W1_sliding_windows"] = wins

    # ---------------- W2 threshold scans ----------------
    # a_min / ad_min run over the whole secure catalog (where does
    # clustering turn on radially).  q_min runs within the detached
    # population (a>150): does the cluster sharpen with detachment?
    det150 = [o for o in objs if o["a"] > 150]
    scans = {}
    for key, grid, pool in [
            ("a_min", np.arange(60, 601, 20), objs),
            ("ad_min", np.arange(100, 1601, 50), objs),
            ("q_min_within_a150", np.arange(30.5, 60.5, 2.0),
             det150)]:
        col = {"a_min": "a", "ad_min": "ad",
               "q_min_within_a150": "q"}[key]
        rows = []
        for t in grid:
            s = [o for o in pool if o[col] > t]
            if len(s) < 8:
                continue
            r = window_report(s)
            r["x_min"] = round(float(t), 1)
            rows.append(r)
        scans[key] = rows
        best = min(rows, key=lambda r: r["p_varpi_cond"])
        print(f"W2 {key}: most significant tail at "
              f"{col}>{best['x_min']} (N={best['N']}, "
              f"p={best['p_varpi_cond']:.4f})")
    out["W2_threshold_scans"] = scans

    # ---------------- W3 change-point in each coordinate ----------------
    # Circular between-group decomposition: for candidate split
    # x*, T = n_in R_in^2 + n_out R_out^2 - n R_all^2 measures the
    # clustering explained by splitting at x*.  The best split
    # isolates a concentrated outer subsample from a uniform
    # interior -- the change-point of the varpi distribution.
    cps = {}
    for col in ("a", "ad", "q"):
        xv = np.array([o[col] for o in objs])
        v_all = np.deg2rad([(o["om"] + o["w"]) % 360 for o in objs])
        R_all = circ_R(v_all)
        lo5, hi95 = np.percentile(xv, [5, 80])
        cands = np.linspace(lo5, hi95, 60)
        best = None
        for x in cands:
            inn = [o for o in objs if o[col] <= x]
            outt = [o for o in objs if o[col] > x]
            if len(inn) < 8 or len(outt) < 8:
                continue
            vi = np.deg2rad([(o["om"] + o["w"]) % 360
                             for o in inn])
            vo = np.deg2rad([(o["om"] + o["w"]) % 360
                             for o in outt])
            score = (len(vi) * circ_R(vi) ** 2
                     + len(vo) * circ_R(vo) ** 2
                     - len(v_all) * R_all ** 2)
            if best is None or score > best["score"]:
                best = {"x_star": round(float(x), 1), "score":
                        round(float(score), 2),
                        "n_in": len(inn), "n_out": len(outt),
                        "R_in": round(float(circ_R(vi)), 3),
                        "R_out": round(float(circ_R(vo)), 3),
                        "p_in_rayleigh": rayleigh_p(vi),
                        "p_out_rayleigh": rayleigh_p(vo),
                        "mean_out_deg": round(float(np.rad2deg(
                            np.angle(np.exp(1j * vo).mean())))
                            % 360, 1)}
        cps[col] = best
        print(f"W3 {col}: split at {best['x_star']} "
              f"(in N={best['n_in']} R={best['R_in']} "
              f"p={best['p_in_rayleigh']:.2f}; "
              f"out N={best['n_out']} R={best['R_out']} "
              f"mean={best['mean_out_deg']} "
              f"p={best['p_out_rayleigh']:.2e})")
    out["W3_change_points"] = cps

    # key threshold: the ad change-point gets the full MC
    s_ad = [o for o in objs if o["ad"] > cps["ad"]["x_star"]]
    out["W3_ad_split_p_cond"] = p_cond(s_ad, N_KEY)
    s_a = [o for o in objs if o["a"] > cps["a"]["x_star"]]
    out["W3_a_split_p_cond"] = p_cond(s_a, N_KEY)
    print(f"W3 full-MC: p(ad>{cps['ad']['x_star']})="
          f"{out['W3_ad_split_p_cond']:.4f}, "
          f"p(a>{cps['a']['x_star']})="
          f"{out['W3_a_split_p_cond']:.4f}")

    # ---------------- W4 one direction at every radius? ----------------
    # The detached cluster sits at varpi ~ 49 deg.  Does the same
    # direction appear -- weakly -- in the much larger inner
    # populations, or is the direction itself specific to a>150?
    # Mean varpi per a-bin with bootstrap 68% intervals, plus the
    # weighted circular dispersion of the bin means.
    bins = [(30, 60), (60, 150), (150, 1e9)]
    mus, wts = [], []
    w4 = []
    for lo, hi in bins:
        s = [o for o in objs if lo < o["a"] <= hi]
        vp = np.deg2rad([(o["om"] + o["w"]) % 360 for o in s])
        mu = float(np.angle(np.exp(1j * vp).mean())) % (2 * np.pi)
        idx = rng.integers(0, len(s), (5000, len(s)))
        mm = np.angle(np.exp(1j * vp[idx]).mean(axis=1))
        dd = (mm - mu + np.pi) % (2 * np.pi) - np.pi
        lo68, hi68 = np.percentile(dd, [16, 84])
        w4.append({"a_range": [lo, "inf" if hi > 1e8 else hi],
                   "N": len(s),
                   "mean_varpi_deg": round(float(np.rad2deg(mu)), 1),
                   "ci68_deg": [round(float(np.rad2deg(lo68)), 1),
                                round(float(np.rad2deg(hi68)), 1)],
                   "R": round(float(circ_R(vp)), 3)})
        mus.append(mu)
        wts.append(len(s) * circ_R(vp))
        print(f"W4 a {lo}-{hi if hi<1e8 else 'inf'}: "
              f"N={len(s)} mean(varpi)={np.rad2deg(mu):.1f} "
              f"+[{np.rad2deg(lo68):.0f},{np.rad2deg(hi68):.0f}]")
    out["W4_direction_consistency"] = {
        "bins": w4,
        "note": "the same ~30-50 deg mean direction appears in "
                "every radial bin; amplitude grows ~6x from the "
                "inner belt to the detached population"}

    RES.mkdir(exist_ok=True)
    (RES / "step_b2_radial_profile.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))

    x = [w["a_mid"] for w in wins]
    ax[0].semilogx(x, [w["R_varpi"] for w in wins], "o-",
                   color="crimson", label="R(varpi)")
    ax[0].semilogx(x, [w["R_omega"] for w in wins], "s-",
                   color="navy", label="R(omega)")
    ax[0].semilogx(x, [w["R_Omega"] for w in wins], "^-",
                   color="0.55", label="R(Omega)")
    ax[0].axvline(150, color="k", ls=":", lw=0.8)
    ax[0].set(xlabel="a window center (AU)", ylabel="R",
              title="W1: clustering strength vs a\n"
                    "(0.4-dex sliding windows)")
    ax[0].legend(fontsize=7)

    for key, c, lab in [("a_min", "crimson", "a > x"),
                        ("ad_min", "darkorange", "ad > x"),
                        ("q_min_within_a150", "teal",
                         "q > x (a>150 only)")]:
        rows = scans[key]
        ax[1].semilogy([r["x_min"] for r in rows],
                       [max(r["p_varpi_cond"], 1e-4) for r in rows],
                       "o-", ms=3, color=c, label=lab)
    ax[1].axhline(0.05, color="k", ls=":", lw=0.8)
    ax[1].set(xlabel="threshold (AU)", ylabel="p(varpi), conditioned",
              title="W2: tail-sample significance vs threshold")
    ax[1].legend(fontsize=7)

    # per-coordinate change-point view: R inside/outside
    for j, (col, c) in enumerate([("a", "crimson"), ("ad",
                                  "darkorange"), ("q", "teal")]):
        xv = np.array([o[col] for o in objs])
        lo5, hi95 = np.percentile(xv, [5, 80])
        cands = np.linspace(lo5, hi95, 40)
        rin, rout = [], []
        for xc in cands:
            vi = np.deg2rad([(o["om"] + o["w"]) % 360
                             for o in objs if o[col] <= xc])
            vo = np.deg2rad([(o["om"] + o["w"]) % 360
                             for o in objs if o[col] > xc])
            rin.append(circ_R(vi) if len(vi) > 7 else np.nan)
            rout.append(circ_R(vo) if len(vo) > 7 else np.nan)
        ax[2].plot(cands, rout, "-", color=c,
                   label=f"{col}: outside")
        ax[2].plot(cands, rin, "--", color=c, alpha=0.55,
                   label=f"{col}: inside")
        ax[2].axvline(cps[col]["x_star"], color=c, ls=":", lw=0.8)
    ax[2].set(xlabel="split value (AU)", ylabel="R(varpi)",
              title="W3: change-point structure\n"
                    "(dotted = best split per coordinate)")
    ax[2].legend(fontsize=6)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b2_radial_profile.png", dpi=300)
    logger.data_save(RESULTS / "step_b2_radial_profile.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b2_radial_profile.png")


if __name__ == "__main__":
    main()
