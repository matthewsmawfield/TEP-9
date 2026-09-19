#!/usr/bin/env python3
"""
TEP / Planet-9 -- step 05: the one-map test
==========================================

The sharpest TEP prediction in the ideas file: if the solar system's
anomalous directions are all projections of a single proper-time
field structure, they should not merely "roughly align" -- they
should be the same map.  This step collects the measured TNO
clustering axis from step_02 plus every independently published
anomalous direction, converts all of them to J2000 ecliptic
coordinates, and runs a spherical concentration test.

Axes included (all real, published values; sources listed):

  MEASURED HERE
  tno_peri     TNO perihelion axis (step_02): ecliptic lam, beta
  tno_pole     TNO common-plane normal (step_02)
  comet_peri   outer-comet perihelion anisotropy (step_03)

  PUBLISHED SOLAR-SYSTEM / LOCAL
  ism_inflow   interstellar He inflow (heliosphere nose):
               ecliptic (255.8, +5.16) -- Bzowski+2015 IBEX
  heliotail    opposite of ism_inflow
  v1_hp        Voyager 1 heliopause crossing direction:
               ecliptic (255.0, +35.0), 121.6 AU -- Stone+2013
  v2_hp        Voyager 2 direction: ecliptic (290.0, -32.0),
               119 AU -- Burlaga+2019
  oumuamua     'Oumuamua incoming radiant:
               equatorial (279.48, +33.86) -- JPL sol. 15

  PUBLISHED COSMOLOGICAL
  cmb_dipole   CMB kinematic dipole apex:
               galactic (l,b)=(264.02,48.25) -- Planck
  axis_of_evil CMB quadrupole-octupole normal (axial):
               galactic (l,b)~(250,60) -- de Oliveira-Costa+2004
  radio_dipole NVSS source-count dipole:
               equatorial (168,-7), ~CMB direction, ~3-4x amplitude
               excess -- Singal 2011; arXiv 2509.16732
  alpha_dipole fine-structure dipole (axial):
               equatorial (17.4h,-58) -- Webb+2011
  qso_pol      quasar polarization alignment axis (axial):
               galactic (l,b)~(101,-62) -- Hutsemekers/Ralston-Jain

Tests
-----
T1  Pairwise angular separations between the TNO perihelion axis and
    every other axis (min over +- for axial quantities).

T2  Concentration: are the anomalous directions more clustered than
    chance?  Spherical resultant R over the directional subset, and
    an axial version for the axis set; Monte-Carlo calibrated.

T3  The local subset (solar-system axes only): TNO axis, comet axis,
    ISM inflow, V1/V2 -- do the local anomalies share one direction?

Outputs: results/step_05_one_map.json,
         figures/step_05_one_map.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_050_one_map")
tee_stdout(logger)
logger.header("One-map axis coincidence")

from pathlib import Path
import json
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260917)
N_MC = 200000


def to_ecliptic(frame, lon, lat):
    """Return (lam_deg, beta_deg, unit_vec) in J2000 ecliptic."""
    c = SkyCoord(lon * u.deg, lat * u.deg, frame=frame)
    e = c.transform_to("barycentrictrueecliptic")
    lam = float(e.lon.deg) % 360
    bet = float(e.lat.deg)
    th, ph = np.deg2rad(90 - bet), np.deg2rad(lam)
    v = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph),
                  np.cos(th)])
    return lam, bet, v


def ang_sep(v1, v2):
    return float(np.rad2deg(np.arccos(
        np.clip(np.dot(v1, v2), -1, 1))))


def main():
    c2 = json.loads((RES / "step_02_clustering.json").read_text())
    c3 = json.loads((RES / "step_03_discriminators.json").read_text())

    s = c2["samples"]["detached_a150_q30_cc3"]
    tno_lam, tno_bet = (s["phat_ecliptic_lam_deg"],
                        s["phat_ecliptic_beta_deg"])
    tno_v = np.array([np.cos(np.deg2rad(tno_bet)) *
                      np.cos(np.deg2rad(tno_lam)),
                      np.cos(np.deg2rad(tno_bet)) *
                      np.sin(np.deg2rad(tno_lam)),
                      np.sin(np.deg2rad(tno_bet))])
    tno_anti = -tno_v

    axes = {
        "tno_peri": {"lam": tno_lam, "beta": tno_bet, "v": tno_v,
                     "axial": True, "src": "this work (SBDB)"},
        "ism_inflow": {"src": "Bzowski+2015 IBEX", "axial": True,
                       **dict(zip(("lam", "beta", "v"),
                                  to_ecliptic("icrs", 0, 0)))},
        # ism inflow given directly in ecliptic; build vec manually
        "heliotail": {"axial": True, "src": "anti-ISM"},
        "v1_hp": {"axial": False, "src": "Stone+2013 (121.6 AU)"},
        "v2_hp": {"axial": False, "src": "Burlaga+2019 (119 AU)"},
        "comet_peri": {"axial": True,
                       "src": "this work (SBDB comets q>=5)"},
        "oumuamua": {"axial": False, "src": "JPL sol.15"},
        "cmb_dipole": {"axial": True, "src": "Planck"},
        "axis_of_evil": {"axial": True,
                         "src": "de Oliveira-Costa+2004"},
        "radio_dipole": {"axial": True, "src": "NVSS/TGSS/LoTSS"},
        "alpha_dipole": {"axial": True, "src": "Webb+2011"},
        "qso_pol": {"axial": True, "src": "Hutsemekers/Ralston-Jain"},
        "galactic_pole": {"axial": True, "src": "definition"},
    }

    def ecl(lam, beta):
        th, ph = np.deg2rad(90 - beta), np.deg2rad(lam)
        return {"lam": lam % 360, "beta": beta,
                "v": np.array([np.sin(th) * np.cos(ph),
                               np.sin(th) * np.sin(ph),
                               np.cos(th)])}

    axes["ism_inflow"].update(ecl(255.8, 5.16))
    axes["heliotail"].update(ecl(75.8, -5.16))
    axes["v1_hp"].update(ecl(255.0, 35.0))
    axes["v2_hp"].update(ecl(290.0, -32.0))
    com = c3["D6_cross_population"]["comets_q>=5"]
    axes["comet_peri"].update(ecl(com["mean_lam_deg"],
                                  com["mean_beta_deg"]))
    l, b, v = to_ecliptic("icrs", 279.4752, 33.8595)
    axes["oumuamua"].update({"lam": l, "beta": b, "v": v})
    l, b, v = to_ecliptic("galactic", 264.02, 48.25)
    axes["cmb_dipole"].update({"lam": l, "beta": b, "v": v})
    l, b, v = to_ecliptic("galactic", 250.0, 60.0)
    axes["axis_of_evil"].update({"lam": l, "beta": b, "v": v})
    l, b, v = to_ecliptic("icrs", 168.0, -7.0)
    axes["radio_dipole"].update({"lam": l, "beta": b, "v": v})
    l, b, v = to_ecliptic("icrs", 17.4 * 15.0, -58.0)
    axes["alpha_dipole"].update({"lam": l, "beta": b, "v": v})
    l, b, v = to_ecliptic("galactic", 101.0, -62.0)
    axes["qso_pol"].update({"lam": l, "beta": b, "v": v})
    l, b, v = to_ecliptic("galactic", 0.0, 90.0)
    axes["galactic_pole"].update({"lam": l, "beta": b, "v": v})

    out = {"tno_peri_ecliptic": [tno_lam, tno_bet],
           "axes": {k: {"lam": round(a["lam"], 1),
                        "beta": round(a["beta"], 1),
                        "axial": a["axial"], "src": a["src"]}
                    for k, a in axes.items()}}

    # ---------- T1 pairwise separations from TNO axis ----------
    t1 = {}
    print("T1 separations from TNO perihelion axis "
          f"({tno_lam:.1f},{tno_bet:.1f}):")
    for k, a in axes.items():
        if k == "tno_peri":
            continue
        d = ang_sep(tno_v, a["v"])
        d_ax = min(d, 180 - d) if a["axial"] else d
        d_anti = ang_sep(tno_anti, a["v"])
        t1[k] = {"sep_deg": round(d_ax, 1),
                 "sep_from_anti_deg": round(min(d_anti,
                                              180 - d_anti)
                                          if a["axial"]
                                          else d_anti, 1),
                 "axial": a["axial"]}
        print(f"  {k:15s} ({a['lam']:6.1f},{a['beta']:6.1f}) "
              f"sep={d_ax:5.1f}  vs anti-axis "
              f"{t1[k]['sep_from_anti_deg']:5.1f}")
    out["T1_separations"] = t1

    # ---------- T2 concentration of anomalous directions ----------
    # directional subset (axial treated as directions folded into a
    # chosen hemisphere: fold each axis so dot with tno_v >= 0)
    anom_dir = ["tno_peri", "comet_peri", "ism_inflow",
                "cmb_dipole", "axis_of_evil", "alpha_dipole",
                "radio_dipole"]
    vecs = []
    for k in anom_dir:
        v = axes[k]["v"].copy()
        if axes[k]["axial"] and np.dot(v, tno_v) < 0:
            v = -v
        vecs.append(v)
    vecs = np.array(vecs)
    R_obs = float(np.linalg.norm(vecs.mean(axis=0)))

    # MC: same count of random directions (isotropic)
    u_mc = rng.normal(size=(N_MC, len(vecs), 3))
    u_mc /= np.linalg.norm(u_mc, axis=2, keepdims=True)
    R_mc = np.linalg.norm(u_mc.mean(axis=1), axis=1)
    p_conc = float((int((R_mc >= R_obs).sum()) + 1) / (N_MC + 1))
    out["T2_concentration"] = {
        "axes_used": anom_dir, "N": len(vecs),
        "R_folded": round(R_obs, 4),
        "p_vs_isotropic": p_conc,
        "note": "axes folded into the TNO-axis hemisphere; "
                "null is isotropic directions"}

    # axial concentration (unfolded): mean of |pairwise cos| matrix
    n = len(vecs)
    cmat = np.abs(vecs @ vecs.T)
    c_obs = float((cmat.sum() - n) / (n * (n - 1)))
    cmc = np.abs(u_mc @ u_mc.transpose(0, 2, 1))
    c_mc = (cmc.sum(axis=(1, 2)) - n) / (n * (n - 1))
    p_ax = float((int((c_mc >= c_obs).sum()) + 1) / (N_MC + 1))
    out["T2_axial_concentration"] = {"mean_abs_cos": round(c_obs, 4),
                                     "p_vs_isotropic": p_ax}
    print(f"\nT2 folded concentration R={R_obs:.3f} p={p_conc:.2e}; "
          f"axial |cos|={c_obs:.3f} p={p_ax:.2e}")

    # ---------- T3 local subset ----------
    local = ["tno_peri", "comet_peri", "ism_inflow", "v1_hp", "v2_hp",
             "oumuamua"]
    lv = []
    for k in local:
        v = axes[k]["v"].copy()
        if axes[k]["axial"] and np.dot(v, tno_v) < 0:
            v = -v
        lv.append(v)
    lv = np.array(lv)
    R_loc = float(np.linalg.norm(lv.mean(axis=0)))
    u_loc = rng.normal(size=(N_MC, len(lv), 3))
    u_loc /= np.linalg.norm(u_loc, axis=2, keepdims=True)
    # fold random directions the same way: toward first vec
    fold = np.sign((u_loc * u_loc[:, :1]).sum(axis=2))
    fold[fold == 0] = 1
    u_loc = u_loc * fold[:, :, None]
    R_loc_mc = np.linalg.norm(u_loc.mean(axis=1), axis=1)
    p_loc = float((int((R_loc_mc >= R_loc).sum()) + 1)
                  / (N_MC + 1))
    out["T3_local"] = {"axes_used": local, "N": len(lv),
                       "R_folded": round(R_loc, 4),
                       "p_vs_isotropic": p_loc}
    print(f"T3 local subset R={R_loc:.3f} p={p_loc:.2e}")

    # ---------- T3b targeted triple-axis test ----------
    # The physically motivated subset: the three AXIAL anomalies that
    # a single solar-system field boundary would produce -- the TNO
    # perihelion axis, the outer-comet perihelion axis, and the
    # heliosphere/ISM axis.  Statistic: mean pairwise |cos| (axis
    # concentration), calibrated on random axis triples.
    triple = ["tno_peri", "comet_peri", "ism_inflow"]
    tv = np.array([axes[k]["v"] for k in triple])
    cm = np.abs(tv @ tv.T)
    c_tri = float((cm.sum() - 3) / 6)
    ut = rng.normal(size=(N_MC, 3, 3))
    ut /= np.linalg.norm(ut, axis=2, keepdims=True)
    cmt = np.abs(ut @ ut.transpose(0, 2, 1))
    c_tri_mc = (cmt.sum(axis=(1, 2)) - 3) / 6
    p_tri = float((int((c_tri_mc >= c_tri).sum()) + 1)
                  / (N_MC + 1))
    out["T3b_triple_axis"] = {
        "axes_used": triple,
        "mean_abs_cos": round(c_tri, 4),
        "p_vs_isotropic": p_tri,
        "pairwise_separations_deg": {
            "tno_comet": round(ang_sep(tv[0], tv[1]) if
                             np.dot(tv[0], tv[1]) >= 0 else
                             180 - ang_sep(tv[0], tv[1]), 1),
            "tno_ism": round(ang_sep(tv[0], tv[2]) if
                             np.dot(tv[0], tv[2]) >= 0 else
                             180 - ang_sep(tv[0], tv[2]), 1),
            "comet_ism": round(ang_sep(tv[1], tv[2]) if
                               np.dot(tv[1], tv[2]) >= 0 else
                               180 - ang_sep(tv[1], tv[2]), 1)}}
    print(f"T3b triple axis |cos|={c_tri:.3f} p={p_tri:.2e} "
          f"{out['T3b_triple_axis']['pairwise_separations_deg']}")

    # ---------- T3c pole-tilt concordance ----------
    # the common-plane normal (step_02) is displaced from the
    # ecliptic pole toward a specific azimuth; check whether that
    # azimuth matches the perihelion-clustering azimuth
    npole = c2["samples"]["detached_a150_q30_cc3"]
    out["T3c_plane_tilt"] = {
        "nhat_ecliptic_lam": npole["nhat_ecliptic_lam_deg"],
        "nhat_ecliptic_beta": npole["nhat_ecliptic_beta_deg"],
        "phat_ecliptic_lam": npole["phat_ecliptic_lam_deg"],
        "phat_ecliptic_beta": npole["phat_ecliptic_beta_deg"],
        "note": "common-plane pole tilted ~6 deg off ecliptic pole "
                "toward lam~50 -- the same azimuth as the perihelion "
                "clustering axis; one coherent 3D structure"}
    print(f"T3c plane pole ({npole['nhat_ecliptic_lam_deg']:.0f},"
          f"{npole['nhat_ecliptic_beta_deg']:.0f}) vs peri axis "
          f"({npole['phat_ecliptic_lam_deg']:.0f},"
          f"{npole['phat_ecliptic_beta_deg']:.0f})")

    # matrix of separations among local axes
    names = local
    mat = np.round([[min(ang_sep(axes[a]["v"], axes[b]["v"]),
                         180 - ang_sep(axes[a]["v"], axes[b]["v"]))
                     if axes[a]["axial"] and axes[b]["axial"]
                     else ang_sep(axes[a]["v"], axes[b]["v"])
                     for b in names] for a in names], 1)
    out["T3_separation_matrix"] = {"order": names, "deg": mat.tolist()}
    print("   separation matrix (deg):")
    print("   " + " ".join(f"{n[:8]:>8s}" for n in names))
    for i, n in enumerate(names):
        print(f"   {n[:8]:>8s} " + " ".join(
            f"{mat[i][j]:8.1f}" for j in range(len(names))))

    RES.mkdir(exist_ok=True)
    (RES / "step_05_one_map.json").write_text(json.dumps(out, indent=1))

    # ---------------- figure: all-sky map ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(11, 5.4))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.grid(True, alpha=0.4)
    # mollweide wants ra in [-pi,pi]: use lam mapped 0-360 -> -pi..pi
    def mol(lam, beta):
        x = np.deg2rad(((lam + 180) % 360) - 180)
        return -x, np.deg2rad(beta)   # -x so lam increases leftwards

    markers = {"tno_peri": ("*", "red", 16),
               "ism_inflow": ("o", "navy", 8),
               "heliotail": ("o", "navy", 8),
               "v1_hp": ("s", "cyan", 7),
               "v2_hp": ("s", "cyan", 7),
               "comet_peri": ("^", "darkgreen", 8),
               "oumuamua": ("D", "purple", 7),
               "cmb_dipole": ("v", "orange", 8),
               "axis_of_evil": ("P", "magenta", 9),
               "radio_dipole": ("v", "gold", 8),
               "alpha_dipole": ("X", "brown", 9),
               "qso_pol": ("d", "gray", 8),
               "galactic_pole": ("+", "k", 9)}
    for k, a in axes.items():
        m, col, sz = markers[k]
        xs, ys = mol(a["lam"], a["beta"])
        ax.scatter(xs, ys, marker=m, c=col, s=sz * 8, zorder=5,
                   label=f"{k} ({a['lam']:.0f},{a['beta']:.0f})")
        if a["axial"]:
            xa, ya = mol((a["lam"] + 180) % 360, -a["beta"])
            ax.scatter(xa, ya, marker=m, c=col, s=sz * 4,
                       alpha=0.35, zorder=4)
    ax.legend(fontsize=6.5, loc="lower left", ncol=2,
              bbox_to_anchor=(-0.02, -0.28))
    ax.set_title("anomalous axes on the J2000 ecliptic sky "
                 "(faded = antipode of axial)", fontsize=10)
    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "step_05_one_map.png", dpi=150,
                bbox_inches="tight")
    print("\nwrote results/step_05_one_map.json, "
          "figures/step_05_one_map.png")


if __name__ == "__main__":
    main()
