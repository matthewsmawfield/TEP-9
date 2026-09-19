#!/usr/bin/env python3
"""
TEP / Planet-9 -- step 03: field-wall vs point-mass discriminators
==================================================================

The two candidate explanations of the extreme-TNO clustering make
different structural predictions for the same catalog:

  Point mass (Planet 9): a moving perturber shepherds orbits.  It
  couples to all orbital elements through Newtonian dynamics;
  confinement should strengthen with semimajor axis; surviving
  objects are expected to occupy mean-motion resonances with the
  perturber (resonant substructure in a); non-aligned objects are
  transitional, scattered smoothly.

  Temporal field boundary (TEP): a spatially fixed lapse gradient.
  Its gradient has one direction in inertial space, so the anomaly
  concentrates in the clock-like elements (perihelion longitudes,
  nodal phases) rather than the geometric elements; selection is by
  *crossing*, so membership should have a sharp edge in the variable
  that controls crossing (perihelion/aphelion geometry), not a smooth
  gradient in a; semimajor axes should be smooth (no resonances);
  every population that crosses the boundary shares the axis.

Tests implemented on the detached sample (a>150 AU, q>30 AU,
condition_code<=3):

D1  Element selectivity: conditioned-MC significance separately for
    node om, argument w, longitude varpi, perihelion vector phat and
    folded pole nhat.

D2  Clustered-mixture structure: uniform + von Mises fit to varpi;
    clustered fraction, center, concentration kappa.

D3  Membership predictor: is alignment with the mean axis controlled
    by a (shepherding expectation) or by the crossing variables
    q, aphelion ad?  Spearman tests + logistic AUC.

D4  Resonant substructure: for a trial perturber axis a9 in
    300-900 AU, does the observed a-set sit preferentially at
    mean-motion resonances (period ratios near small-integer p/q)?
    Statistic: mean distance of (a/a9)^{3/2} to the nearest small
    rational, calibrated by Monte-Carlo.

D5  Sharp-edge test: does the transition between aligned and
    unaligned objects occur sharply at a boundary in perihelion
    direction space (wall), or smoothly (shepherd)?

D6  Cross-population: Centaurs (class CEN) -- the Neptune-coupled
    population -- should show no axis; outer comets are checked for a
    preferred perihelion axis for comparison with the TNO axis.

Outputs: results/step_03_discriminators.json,
         figures/step_03_discriminators.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_011_discriminators")
tee_stdout(logger)
logger.header("Field-wall vs point-mass discriminators")

from pathlib import Path
import json
import numpy as np
from scipy import stats as st
from scipy.optimize import minimize
from scipy.special import i0, i1

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260917)
N_MC = 20000


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


def perih_vec(om, w, i):
    return np.stack([
        np.cos(om) * np.cos(w) - np.sin(om) * np.sin(w) * np.cos(i),
        np.sin(om) * np.cos(w) + np.cos(om) * np.sin(w) * np.cos(i),
        np.sin(w) * np.sin(i)], axis=-1)


def pole_vec(om, i):
    n = np.stack([np.sin(i) * np.sin(om),
                  -np.sin(i) * np.cos(om),
                  np.cos(i)], axis=-1)
    n[n[..., 2] < 0] *= -1
    return n


def sel(objs, a_min, q_min, cc_max=None, a_max=np.inf):
    return [o for o in objs
            if np.isfinite(o["a"]) and np.isfinite(o["q"])
            and a_min < o["a"] <= a_max and o["q"] > q_min
            and (cc_max is None or (o["condition_code"] is not None
                 and int(o["condition_code"]) <= cc_max))]


def main():
    objs = load_objects("sbdb_outer_ss.json")
    cen = load_objects("sbdb_centaurs.json")
    com = load_objects("sbdb_comets_ext.json")
    s = sel(objs, 150, 30, cc_max=3)
    out = {"sample": "a>150, q>30, cc<=3", "N": len(s)}

    om = np.deg2rad(np.array([o["om"] for o in s]))
    w = np.deg2rad(np.array([o["w"] for o in s]))
    i = np.deg2rad(np.array([o["i"] for o in s]))
    varpi = (om + w) % (2 * np.pi)
    mean_v = np.angle(np.exp(1j * varpi).mean())

    # ---------- D1 element selectivity ----------
    i_b = np.broadcast_to(i, (N_MC, len(s)))
    om_mc = rng.uniform(0, 2 * np.pi, (N_MC, len(s)))
    w_mc = rng.uniform(0, 2 * np.pi, (N_MC, len(s)))
    varpi_mc = (om_mc + w_mc) % (2 * np.pi)

    null_R = {
        "om": abs(np.exp(1j * om_mc).mean(axis=1)),
        "w": abs(np.exp(1j * w_mc).mean(axis=1)),
        "varpi": abs(np.exp(1j * varpi_mc).mean(axis=1)),
        "phat": np.linalg.norm(perih_vec(om_mc, w_mc, i_b).mean(axis=1),
                              axis=1),
        "nhat": np.linalg.norm(pole_vec(om_mc, i_b).mean(axis=1),
                               axis=1),
    }
    obs_R = {
        "om": abs(np.exp(1j * om).mean()),
        "w": abs(np.exp(1j * w).mean()),
        "varpi": abs(np.exp(1j * varpi).mean()),
        "phat": np.linalg.norm(perih_vec(om, w, i).mean(axis=0)),
        "nhat": np.linalg.norm(pole_vec(om, i).mean(axis=0)),
    }
    d1 = {}
    for k in obs_R:
        d1[k] = {"R_obs": round(float(obs_R[k]), 4),
                 "R_null_mean": round(float(null_R[k].mean()), 4),
                 "p_conditioned": float(
                     (int((null_R[k] >= obs_R[k]).sum()) + 1)
                     / (null_R[k].size + 1))}
        print(f"D1 {k:6s} R={obs_R[k]:.4f} null={null_R[k].mean():.4f} "
              f"p={d1[k]['p_conditioned']:.2e}")
    out["D1_element_selectivity"] = d1

    # ---------- D2 von Mises + uniform mixture ----------
    def nll(par):
        mu, kap, f = par
        if kap < 0 or kap > 200 or f < 0 or f > 1:
            return 1e12
        vm = np.exp(kap * np.cos(varpi - mu)) / (2 * np.pi * i0(kap))
        return -np.sum(np.log(f * vm + (1 - f) / (2 * np.pi)))

    best = None
    for mu0 in np.linspace(0, 2 * np.pi, 25):
        for k0 in (0.5, 1, 2, 5):
            r = minimize(nll, [mu0, k0, 0.6],
                         bounds=[(0, 2 * np.pi), (0, 200), (0, 1)])
            if best is None or r.fun < best.fun:
                best = r
    mu, kap, f = best.x
    out["D2_mixture"] = {"mu_deg": round(float(np.rad2deg(mu) % 360), 1),
                         "kappa": round(float(kap), 2),
                         "sigma_deg": round(float(np.rad2deg(
                             np.sqrt(-2 * np.log(i1(kap) / i0(kap))))
                             if kap > 0.05 else 180), 1),
                         "f_clustered": round(float(f), 2),
                         "nll": round(float(best.fun), 2)}
    print(f"D2 mixture: mu={out['D2_mixture']['mu_deg']} deg, "
          f"kappa={kap:.2f}, f={f:.2f}")

    # ---------- D3 membership predictor ----------
    align = np.cos(varpi - mean_v)
    a = np.array([o["a"] for o in s])
    q = np.array([o["q"] for o in s])
    ad = np.array([o["ad"] for o in s])
    ii = np.array([o["i"] for o in s])
    d3 = {}
    for name, x in [("a", a), ("q", q), ("ad", ad), ("i", ii),
                    ("log_a", np.log10(a))]:
        rho, p = st.spearmanr(align, x)
        d3[name] = {"spearman_rho": round(float(rho), 3),
                    "p": float(p)}
        print(f"D3 align vs {name:6s} rho={rho:+.3f} p={p:.3f}")
    # logistic comparison: AUC for each variable predicting "aligned"
    aligned = align > np.cos(np.deg2rad(60))  # within +-60 deg of mean
    for name, x in [("a", a), ("q", q), ("ad", ad)]:
        xs = (x - x.mean()) / x.std()
        A = np.vstack([np.ones_like(xs), xs]).T
        try:
            beta, *_ = np.linalg.lstsq(
                np.hstack([A, np.ones((len(xs), 1))]),
                aligned.astype(float), rcond=None)
        except Exception:
            beta = None
        # Mann-Whitney AUC: P(x_aligned > x_unaligned)
        xa, xu = x[aligned], x[~aligned]
        auc = float(st.mannwhitneyu(xa, xu).statistic /
                    (len(xa) * len(xu))) if len(xu) else np.nan
        d3[name]["AUC_aligned_vs_not"] = round(auc, 3)
        print(f"   AUC({name})={auc:.3f}  "
              f"(n_aligned={aligned.sum()}, n_not={(~aligned).sum()})")
    out["D3_membership"] = d3

    # ---------- D4 resonant substructure ----------
    # period ratio r = (a/a9)^{3/2}; resonant if r near small p/q
    res = np.array(sorted({p / qq for qq in range(1, 9)
                           for p in range(qq + 1, qq + 9)}))
    res = res[(res > 1) & (res < 20)]

    def res_stat(a_vals, a9):
        r = (a_vals / a9) ** 1.5
        idx = np.clip(np.searchsorted(res, r), 1, len(res) - 1)
        d = np.minimum(np.abs(r - res[idx]),
                       np.abs(r - res[idx - 1]))
        return d.mean()

    a_det = a
    # null: resample the observed a-set (preserves the marginal);
    # shepherding claims resonant concentration ON TOP of it
    a_null = np.array([rng.choice(a_det, len(a_det), replace=True)
                       for _ in range(2000)])
    a9_grid = np.linspace(200, 1500, 261)
    T_obs = np.array([res_stat(a_det, a9) for a9 in a9_grid])
    T_null = np.array([res_stat(a_null[j], a9)
                       for j in range(2000) for a9 in a9_grid])
    T_null = T_null.reshape(2000, len(a9_grid))
    z = (T_null.mean(axis=0) - T_obs) / T_null.std(axis=0)
    jbest = int(np.argmax(z))
    # global significance: fraction of null realizations whose own
    # max-z exceeds observed max-z
    znull = (T_null.mean(axis=0)[None, :] - T_null) / \
        T_null.std(axis=0)[None, :]
    maxz_null = znull.max(axis=1)
    p_global = float((int((maxz_null >= z[jbest]).sum()) + 1)
                     / (maxz_null.size + 1))
    out["D4_resonances"] = {
        "a9_best_AU": float(a9_grid[jbest]),
        "z_max": round(float(z[jbest]), 2),
        "p_global_vs_loguniform": p_global,
        "resonance_set": "p/q, 1<=q<=8, p<q+9, 1<r<20",
        "note": "shepherding predicts a values concentrated at "
                "mean-motion resonances; a smooth a-set (low z) is "
                "the field-wall expectation"}
    print(f"D4 resonances: best a9={a9_grid[jbest]:.0f} AU, "
          f"z={z[jbest]:.2f}, global p={p_global:.3f}")

    # ---------- D5 sharp edge in varpi membership ----------
    # wall predicts a discontinuity; test largest gap between
    # consecutive sorted varpi (log-spaced gap significance)
    sv = np.sort(varpi)
    gaps = np.diff(np.r_[sv, sv[0] + 2 * np.pi])
    jg = int(np.argmax(gaps))
    gap_deg = float(np.rad2deg(gaps[jg]))
    # MC: expected max gap for a clustered+uniform mixture? use the
    # fitted mixture: draw n angles, max gap distribution
    vm_mc = []
    for _ in range(4000):
        u = rng.uniform(0, 1, len(s))
        samp = np.empty(len(s))
        mask = u < f
        samp[mask] = rng.vonmises(mu, kap, mask.sum()) % (2 * np.pi)
        samp[~mask] = rng.uniform(0, 2 * np.pi, (~mask).sum())
        sm = np.sort(samp)
        vm_mc.append(np.rad2deg(
            np.diff(np.r_[sm, sm[0] + 2 * np.pi]).max()))
    p_gap = float(((np.array(vm_mc) >= gap_deg).sum() + 1)
                  / (len(vm_mc) + 1))
    out["D5_max_gap"] = {"max_gap_deg": round(gap_deg, 1),
                         "gap_edges_deg": [round(float(np.rad2deg(
                             sv[jg])) % 360, 1),
                             round(float(np.rad2deg(
                                 sv[(jg + 1) % len(sv)])) % 360, 1)],
                         "p_vs_fitted_mixture": p_gap}
    print(f"D5 max gap in varpi: {gap_deg:.0f} deg between "
          f"{out['D5_max_gap']['gap_edges_deg']}, p={p_gap:.3f}")

    # ---------- D6 cross-population ----------
    d6 = {}
    # Centaurs with a>30 AU (Neptune-coupled control)
    cen30 = [o for o in cen if np.isfinite(o["a"]) and o["a"] > 30
             and np.isfinite(o["om"]) and np.isfinite(o["w"])]
    if cen30:
        vc = np.deg2rad([(o["om"] + o["w"]) % 360 for o in cen30])
        R_c = abs(np.exp(1j * vc).mean())
        p_c = float(np.exp(-len(vc) * R_c**2))
        d6["centaurs_a>30"] = {"N": len(vc), "R": round(float(R_c), 3),
                               "p_rayleigh_approx": p_c,
                               "varpi_mean_deg": round(float(np.rad2deg(
                                   np.angle(np.exp(1j * vc).mean()))
                                   % 360), 1)}
        print(f"D6 Centaurs a>30: N={len(vc)} R={R_c:.3f} p={p_c:.3f}")
    # outer comets (q>=5): perihelion-direction anisotropy on sphere,
    # total and per dynamical class
    co = [o for o in com if np.isfinite(o["om"]) and np.isfinite(o["w"])
          and np.isfinite(o["i"])]

    def sphere_rep(objs):
        pc = perih_vec(np.deg2rad([o["om"] for o in objs]),
                       np.deg2rad([o["w"] for o in objs]),
                       np.deg2rad([o["i"] for o in objs]))
        R_p = np.linalg.norm(pc.mean(axis=0))
        p_p = float(np.exp(-3 * len(objs) * R_p**2 / 2))
        mv_c = pc.mean(axis=0) / np.linalg.norm(pc.mean(axis=0))
        return {"N": len(objs), "R": round(float(R_p), 3),
                "p_sphere_approx": p_p,
                "mean_lam_deg": round(float(np.rad2deg(np.arctan2(
                    mv_c[1], mv_c[0])) % 360), 1),
                "mean_beta_deg": round(float(np.rad2deg(
                    np.arcsin(mv_c[2]))), 1)}

    if co:
        d6["comets_q>=5"] = sphere_rep(co)
        print(f"D6 comets q>=5: N={len(co)} "
              f"R={d6['comets_q>=5']['R']} "
              f"p={d6['comets_q>=5']['p_sphere_approx']:.1e}")
        d6["comets_by_class"] = {}
        for cl in sorted(set(o["class"] for o in co)):
            sub = [o for o in co if o["class"] == cl]
            if len(sub) > 3:
                d6["comets_by_class"][cl] = sphere_rep(sub)
                r = d6["comets_by_class"][cl]
                print(f"   {cl}: N={r['N']} R={r['R']} "
                      f"p={r['p_sphere_approx']:.1e} -> "
                      f"({r['mean_lam_deg']},{r['mean_beta_deg']})")
        d6["comets_note"] = (
            "q>=5 comets: JFc perihelion longitudes are known to be "
            "non-uniform from capture/discovery geometry; HYP are "
            "hyperbolic (Oort/inbound) objects. Discovery bias is "
            "not modeled here -- axes are descriptive.")
    out["D6_cross_population"] = d6

    # ---------- D7 harmonic decomposition ----------
    # directional (m=1: a vector gradient -- TEP lapse gradient or a
    # point perturber) vs axial (m=2: quadrupolar, e.g. galactic
    # tide) vs higher-order power in the varpi distribution
    d7 = {}
    n = len(varpi)
    for m in (1, 2, 3):
        R_m = abs(np.exp(1j * m * varpi).mean())
        z_m = n * R_m**2
        p_m = float(np.exp(-z_m) *
                    (1 + (2 * z_m - z_m**2) / (4 * n)
                     - (24 * z_m - 76 * z_m**2 + 9 * z_m**3)
                     / (24 * n**2)))
        mu_m = float(np.rad2deg(
            np.angle(np.exp(1j * m * varpi).mean()) / m) % 360)
        d7[f"m{m}"] = {"R": round(float(R_m), 4), "z": round(z_m, 2),
                       "p_rayleigh": float(p_m),
                       "mean_dir_deg": round(mu_m, 1)}
        print(f"D7 m={m} R={R_m:.3f} z={z_m:.2f} p={p_m:.2e} "
              f"dir={mu_m:.0f}")
    out["D7_harmonic_decomp"] = d7
    out["D7_note"] = ("m=1 dominance => directed (vector) anomaly: "
                      "consistent with a lapse gradient OR a point "
                      "mass; an axial/quadrupolar field (e.g. "
                      "galactic tide) would peak at m=2")

    # per-object alignment table
    out["object_alignment"] = sorted(
        [{"name": o["full_name"].strip(),
          "varpi": round(float(np.rad2deg(v)) % 360, 1),
          "align_cos": round(float(np.cos(v - mean_v)), 3),
          "a": o["a"], "q": o["q"], "ad": o["ad"], "i": o["i"]}
         for o, v in zip(s, varpi)],
        key=lambda r: -r["align_cos"])

    RES.mkdir(exist_ok=True)
    (RES / "step_03_discriminators.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(10.5, 8))
    ax = ax.ravel()

    ax[0].hist(np.rad2deg(varpi), bins=36, range=(0, 360),
               color="steelblue", edgecolor="k", lw=0.4)
    th = np.linspace(0, 2 * np.pi, 720)
    pdf = f * np.exp(kap * np.cos(th - mu)) / (2 * np.pi * i0(kap)) \
        + (1 - f) / (2 * np.pi)
    ax[0].plot(np.rad2deg(th), pdf * len(s) * 10, "r-", lw=2,
               label=f"fit: f={f:.2f}, $\\kappa$={kap:.1f}")
    ax[0].axvline(np.rad2deg(mean_v) % 360, color="r", ls=":")
    ax[0].set_xlabel("$\\varpi$ (deg)")
    ax[0].legend(fontsize=8)
    ax[0].set_title("D2: uniform + von Mises mixture", fontsize=9)

    ax[1].scatter(a, align, s=22, c="steelblue", edgecolor="k", lw=0.3)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("a (AU)")
    ax[1].set_ylabel("cos($\\varpi - \\bar{\\varpi}$)")
    ax[1].set_title(f"D3: alignment vs a "
                    f"(rho={d3['a']['spearman_rho']:+.2f}, "
                    f"p={d3['a']['p']:.3f})", fontsize=9)

    ax[2].plot(a9_grid, z, "k-", lw=1.5)
    ax[2].axhline(0, color="0.6", lw=0.5)
    ax[2].set_xlabel("trial perturber a$_9$ (AU)")
    ax[2].set_ylabel("resonance-excess z")
    ax[2].set_title(f"D4: resonant substructure scan "
                    f"(max z={z[jbest]:.1f}, p={p_global:.2f})",
                    fontsize=9)

    qq = np.array([o["q"] for o in objs
                   if np.isfinite(o["a"]) and o["a"] > 150
                   and np.isfinite(o["q"])
                   and o["condition_code"] is not None
                   and int(o["condition_code"]) <= 3])
    ax[3].hist(qq, bins=40, range=(25, 90), color="steelblue",
               edgecolor="k", lw=0.4)
    ax[3].axvline(30, color="r", ls=":")
    ax[3].set_xlabel("q (AU)")
    ax[3].set_title("perihelion distribution, a>150 AU (cc<=3)",
                    fontsize=9)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "step_03_discriminators.png", dpi=150)
    print("\nwrote results/step_03_discriminators.json, "
          "figures/step_03_discriminators.png")


if __name__ == "__main__":
    main()
