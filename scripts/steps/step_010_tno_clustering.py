#!/usr/bin/env python3
"""
TEP / Planet-9 -- step 02: the perihelion-clustering anomaly
============================================================

Reproduces, from the full JPL SBDB catalog, the anomalous clustering
of distant trans-Neptunian object orbits that motivates the Planet-9
hypothesis (Trujillo & Sheppard 2014; Batygin & Brown 2016;
Brown & Batygin 2021).

Under the standard interpretation the clustering is shepherding by an
unseen point mass.  Under TEP the clustered perihelion longitudes are
orbital clocks phase-locked about a fixed inertial direction -- the
signature of a spatially structured proper-time field, not a body.

Method -- a *conditioned* null.  Each orbit is decomposed into
geometric elements (a, e, i: shape and plane tilt) and clock elements
(node om, argument of perihelion w, and varpi = om + w).  The null
hypothesis holds every object's observed (a, e, i) fixed and draws om,
w uniform on the circle.  Clustering is therefore tested in the
angular/clock sector only, conditioned on the observed geometric
sector -- precisely where a temporal-field boundary acts, and free of
the inclination-marginal bias that afflicts naive pole-vector tests.

Sample definitions follow the literature convention:
  detached:    a > 150 AU, q > 30 AU
  extreme:     a > 250 AU, q > 30 AU
  sednoid:     a > 150 AU, q > 50 AU
Orbit quality: primary results use condition_code <= 3 (secure,
multi-opposition orbits), as is standard in this literature; the
unfiltered catalog is shown for contrast.

Statistics per sample: Rayleigh R of varpi; spherical resultant R of
the perihelion unit vectors phat; spherical R of the folded orbit
poles nhat; each calibrated against 20000 conditioned Monte-Carlo
realizations.  Bootstrap 68% intervals on the mean direction.

Outputs: results/step_02_clustering.json,
         figures/step_02_clustering.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_010_tno_clustering")
tee_stdout(logger)
logger.header("Detached-TNO orientation clustering")

from pathlib import Path
import json
import numpy as np

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
    """Unit vector to the perihelion point; angles in radians."""
    return np.stack([
        np.cos(om) * np.cos(w) - np.sin(om) * np.sin(w) * np.cos(i),
        np.sin(om) * np.cos(w) + np.cos(om) * np.sin(w) * np.cos(i),
        np.sin(w) * np.sin(i)], axis=-1)


def pole_vec(om, i):
    """Orbit normal; folded so n_z >= 0 (planes as axes)."""
    n = np.stack([np.sin(i) * np.sin(om),
                  -np.sin(i) * np.cos(om),
                  np.cos(i)], axis=-1)
    n[n[..., 2] < 0] *= -1
    return n


def stats_of(om, w, i):
    """The three clustering statistics for one realization."""
    varpi = (om + w) % (2 * np.pi)
    R_varpi = abs(np.exp(1j * varpi).mean())
    ph = perih_vec(om, w, i)
    R_phat = np.linalg.norm(ph.mean(axis=0))
    nh = pole_vec(om, i)
    R_nhat = np.linalg.norm(nh.mean(axis=0))
    return R_varpi, R_phat, R_nhat


def conditioned_mc(obs, stat_obs):
    """Null: fix (a,e,i); draw om,w uniform.  Return MC p-values."""
    i_rad = np.deg2rad(np.array([o["i"] for o in obs]))
    n = len(obs)
    om_mc = rng.uniform(0, 2 * np.pi, (N_MC, n))
    w_mc = rng.uniform(0, 2 * np.pi, (N_MC, n))
    i_b = np.broadcast_to(i_rad, (N_MC, n))

    varpi = (om_mc + w_mc) % (2 * np.pi)
    Rv = abs(np.exp(1j * varpi).mean(axis=1))

    ph = perih_vec(om_mc, w_mc, i_b)
    Rp = np.linalg.norm(ph.mean(axis=1), axis=1)

    nh = pole_vec(om_mc, i_b)
    Rn = np.linalg.norm(nh.mean(axis=1), axis=1)

    return [(float((int((Rv >= stat_obs[0]).sum()) + 1)
                  / (N_MC + 1)), float(Rv.mean())),
            (float((int((Rp >= stat_obs[1]).sum()) + 1)
                  / (N_MC + 1)), float(Rp.mean())),
            (float((int((Rn >= stat_obs[2]).sum()) + 1)
                  / (N_MC + 1)), float(Rn.mean()))]


def boot_ci(ang, n_boot=20000):
    th = np.asarray(ang)
    n = th.size
    idx = rng.integers(0, n, (n_boot, n))
    m = np.angle(np.exp(1j * th[idx]).mean(axis=1))
    c = np.angle(np.exp(1j * th).mean())
    d = np.rad2deg((m - c + np.pi) % (2 * np.pi) - np.pi)
    lo, hi = np.percentile(d, [16, 84])
    return float(lo), float(hi)


def sample_report(objs, label):
    om = np.deg2rad(np.array([o["om"] for o in objs]))
    w = np.deg2rad(np.array([o["w"] for o in objs]))
    i = np.deg2rad(np.array([o["i"] for o in objs]))
    varpi = (om + w) % (2 * np.pi)

    s_obs = stats_of(om, w, i)
    (pv, mv), (pp, mp), (pn, mn) = conditioned_mc(objs, s_obs)

    mean_v = float(np.rad2deg(np.angle(np.exp(1j * varpi).mean())) % 360)
    lo, hi = boot_ci(varpi)

    ph = perih_vec(om, w, i)
    mvec = ph.mean(axis=0)
    mvec /= np.linalg.norm(mvec)
    lam = float(np.rad2deg(np.arctan2(mvec[1], mvec[0])) % 360)
    bet = float(np.rad2deg(np.arcsin(mvec[2])))

    nh = pole_vec(om, i)
    nvec = nh.mean(axis=0)
    nvec /= np.linalg.norm(nvec)
    lam_n = float(np.rad2deg(np.arctan2(nvec[1], nvec[0])) % 360)
    bet_n = float(np.rad2deg(np.arcsin(nvec[2])))

    rep = {
        "label": label, "N": len(objs),
        "varpi_deg": sorted(round(float(np.rad2deg(x)), 2)
                            for x in varpi),
        "varpi_mean_deg": round(mean_v, 1),
        "varpi_mean_ci68": [round(mean_v + lo, 1),
                            round(mean_v + hi, 1)],
        "R_varpi": round(s_obs[0], 4),
        "R_varpi_null_mean": round(mv, 4),
        "p_varpi_conditioned": pv,
        "R_phat": round(s_obs[1], 4),
        "R_phat_null_mean": round(mp, 4),
        "p_phat_conditioned": pp,
        "phat_ecliptic_lam_deg": round(lam, 1),
        "phat_ecliptic_beta_deg": round(bet, 1),
        "R_nhat_folded": round(s_obs[2], 4),
        "R_nhat_null_mean": round(mn, 4),
        "p_nhat_conditioned": pn,
        "nhat_ecliptic_lam_deg": round(lam_n, 1),
        "nhat_ecliptic_beta_deg": round(bet_n, 1),
    }
    return rep


def main():
    objs = load_objects("sbdb_outer_ss.json")
    print(f"catalog: {len(objs)} objects with a >= 30 AU")

    def sel(a_min, q_min, cc_max=None):
        return [o for o in objs
                if np.isfinite(o["a"]) and np.isfinite(o["q"])
                and o["a"] > a_min and o["q"] > q_min
                and (cc_max is None
                     or (o["condition_code"] is not None
                         and int(o["condition_code"]) <= cc_max))]

    out = {"catalog": "JPL SBDB sbdb_query.api, asteroids a>=30 AU",
           "n_catalog": len(objs),
           "null_model": "fix (a,e,i); om,w uniform; 20000 realizations",
           "samples": {}}

    for label, am, qm, cc in [
            ("detached_a150_q30_cc3", 150, 30, 3),
            ("extreme_a250_q30_cc3", 250, 30, 3),
            ("sednoid_a150_q50_cc3", 150, 50, 3),
            ("detached_a150_q30_allcc", 150, 30, None),
            ("extreme_a150_q40_cc3", 150, 40, 3)]:
        s = sel(am, qm, cc)
        rep = sample_report(s, label)
        out["samples"][label] = rep
        print(f"\n== {label}: N={rep['N']}")
        print(f"   varpi mean={rep['varpi_mean_deg']} "
              f"CI68 {rep['varpi_mean_ci68']} "
              f"R={rep['R_varpi']} (null {rep['R_varpi_null_mean']}) "
              f"p={rep['p_varpi_conditioned']:.2e}")
        print(f"   phat R={rep['R_phat']} (null {rep['R_phat_null_mean']}) "
              f"p={rep['p_phat_conditioned']:.2e} "
              f"-> lam={rep['phat_ecliptic_lam_deg']} "
              f"beta={rep['phat_ecliptic_beta_deg']}")
        print(f"   nhat R={rep['R_nhat_folded']} "
              f"(null {rep['R_nhat_null_mean']}) "
              f"p={rep['p_nhat_conditioned']:.2e} "
              f"-> lam={rep['nhat_ecliptic_lam_deg']} "
              f"beta={rep['nhat_ecliptic_beta_deg']}")

    # per-object table for the primary sample
    s = sel(150, 30, 3)
    table = sorted([{"name": o["full_name"].strip(), "a": o["a"],
                     "q": o["q"], "i": o["i"], "om": o["om"], "w": o["w"],
                     "varpi": round((o["om"] + o["w"]) % 360, 1),
                     "cc": o["condition_code"],
                     "arc_days": o["data_arc"]}
                    for o in s], key=lambda r: r["varpi"])
    out["detached_table"] = table

    RES.mkdir(exist_ok=True)
    (RES / "step_02_clustering.json").write_text(json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rep = out["samples"]["detached_a150_q30_cc3"]
    varpi = np.deg2rad(np.array(rep["varpi_deg"]))
    mv = np.deg2rad(rep["varpi_mean_deg"])

    fig = plt.figure(figsize=(13, 4.4))
    ax1 = fig.add_subplot(131, projection="polar")
    ax1.hist(varpi, bins=18, range=(0, 2 * np.pi), color="steelblue",
             alpha=0.85, edgecolor="k", lw=0.4)
    ax1.plot([mv, mv], [0, ax1.get_ylim()[1]], "r-", lw=2)
    ax1.set_theta_zero_location("N")
    ax1.set_title("$\\varpi$ distribution, detached sample\n"
                  f"a>150, q>30, cc<=3 (N={rep['N']}), "
                  f"p={rep['p_varpi_conditioned']:.1e}", fontsize=9)

    ax2 = fig.add_subplot(132)
    a_all = np.array([o["a"] for o in objs if np.isfinite(o["a"])])
    q_all = np.array([o["q"] for o in objs if np.isfinite(o["q"])])
    m = a_all < 5000
    ax2.scatter(a_all[m], q_all[m], s=1, c="0.75", rasterized=True)
    sa = np.array([o["a"] for o in s])
    sq = np.array([o["q"] for o in s])
    sv = np.deg2rad(np.array([t["varpi"] for t in table]))
    ax2.scatter(sa, sq, s=26, c=np.cos(sv - mv), cmap="RdBu_r",
                vmin=-1, vmax=1, edgecolor="k", lw=0.3, zorder=5)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("a (AU)"); ax2.set_ylabel("q (AU)")
    ax2.set_title("detached objects colored by\n"
                  "alignment with mean $\\varpi$", fontsize=9)

    ax3 = fig.add_subplot(133)
    for lb, c in [("detached_a150_q30_cc3", "crimson"),
                  ("extreme_a250_q30_cc3", "navy"),
                  ("detached_a150_q30_allcc", "0.6")]:
        v = np.array(out["samples"][lb]["varpi_deg"])
        ax3.hist(v, bins=36, range=(0, 360), histtype="step",
                 color=c, lw=1.4, label=lb + f" (N={out['samples'][lb]['N']})")
    ax3.axvline(rep["varpi_mean_deg"], color="r", lw=2)
    ax3.set_xlabel("$\\varpi$ (deg)")
    ax3.legend(fontsize=7)
    ax3.set_title("$\\varpi$ histograms by sample", fontsize=9)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "step_02_clustering.png", dpi=150)
    print("\nwrote results/step_02_clustering.json, "
          "figures/step_02_clustering.png")


if __name__ == "__main__":
    main()
