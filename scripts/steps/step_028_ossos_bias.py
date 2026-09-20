#!/usr/bin/env python3
"""
TEP / Planet-9 -- step 06: bias control with the OSSOS ensemble
================================================================

The standard objection to the extreme-TNO clustering is discovery
bias: an eccentric object is brightest near perihelion, so a survey
pointing at ecliptic longitude lam_d preferentially detects objects
whose perihelion directions lie near lam_d.  The characterized OSSOS
ensemble (Bannister+ 2018, ApJS 236, 18; VizieR J/ApJS/236/18) is
the only TNO sample with fully documented discovery circumstances --
per-object discovery RA/Dec, distance, magnitude, and survey block.

Three real-data tests on the 31 detached objects (cl='det'):

B1  Pointing coupling.  Circular correlation between discovery
    ecliptic longitude lam_d and longitude of perihelion varpi.
    Under the "found near perihelion" bias, varpi - lam_d ~ 0.

B2  Is OSSOS's detached sample direction different from the global
    (SBDB) cluster axis?  OSSOS's deep blocks concentrated near
    lam ~ 210-240 deg -- the anti-cluster side.  A uniform
    population filtered through that footprint produces an apparent
    concentration on the anti side: explicit demonstration of the
    bias mechanism, and evidence the global cluster is not an
    OSSOS artifact.

B3  Discovery geometry of the clustered group: were objects in the
    clustered direction (varpi ~ 50 deg) found only by blocks
    pointing at ~50 deg (suspicious), or across diverse pointings?

B4  Discovery distance: are clustered objects systematically nearer
    perihelion at discovery (the bias channel's expected marker)?

Outputs: results/step_06_ossos_bias.json,
         figures/supplementary/step_06_ossos_bias.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_028_ossos_bias")
tee_stdout(logger)
logger.header("OSSOS pointing coupling")

from pathlib import Path
import json
import numpy as np
from astropy.io.votable import parse
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260917)


def ecl_lon(ra_deg, dec_deg):
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    e = c.transform_to("barycentrictrueecliptic")
    return float(e.lon.deg) % 360


def circ_R(ang_rad):
    return abs(np.exp(1j * np.asarray(ang_rad)).mean())


def rayleigh_p(ang_rad):
    a = np.asarray(ang_rad)
    n = len(a)
    z = n * circ_R(a) ** 2
    return float(np.exp(-z) * (1 + (2 * z - z**2) / (4 * n)))


def main():
    t = parse((DATA_RAW / "ossos" / "ossos_t3char.vot")).get_first_table().to_table()
    cl = np.array([str(x) for x in t["cl"]])
    a = np.array(t["a"], dtype=float)
    e = np.array(t["e"], dtype=float)
    Om = np.array(t["Omega"], dtype=float)
    w = np.array(t["omega"], dtype=float)
    ra = np.array(t["RAJ2000"], dtype=float)
    de = np.array(t["DEJ2000"], dtype=float)
    dist = np.array(t["Dist"], dtype=float)
    q = a * (1 - e)
    varpi = np.deg2rad((Om + w) % 360)
    lam_d = np.deg2rad([ecl_lon(r, d) for r, d in zip(ra, de)])

    out = {"catalog": "VizieR J/ApJS/236/18 OSSOS t3char",
           "n_total": int(len(t))}

    det = cl == "det"
    vd, ld = varpi[det], lam_d[det]
    out["n_detached"] = int(det.sum())

    # ---------- B1 pointing coupling ----------
    dphi = (vd - ld + np.pi) % (2 * np.pi) - np.pi
    R_c = circ_R(dphi)
    # null: shuffle discovery longitudes among objects
    R_null = np.array([circ_R(vd - rng.permutation(ld))
                       for _ in range(20000)])
    p_c = float((int((R_null >= R_c).sum()) + 1)
                / (R_null.size + 1))
    out["B1_pointing_coupling"] = {
        "mean_abs_dphi_deg": round(float(np.rad2deg(
            np.mean(np.abs(dphi)))), 1),
        "R_coupling": round(float(R_c), 3),
        "R_null_mean": round(float(R_null.mean()), 3),
        "p_vs_shuffled": p_c,
        "note": "R~1 => objects found at their own perihelion "
                "direction (full pointing bias); R~null => "
                "discovery longitude independent of varpi"}
    print(f"B1 |varpi-lam_d| mean={out['B1_pointing_coupling']['mean_abs_dphi_deg']} deg, "
          f"R={R_c:.3f} (null {R_null.mean():.3f}) p={p_c:.2e}")

    # ---------- B2 OSSOS direction vs global cluster ----------
    mu_d = float(np.rad2deg(np.angle(np.exp(1j * vd).mean())) % 360)
    out["B2_ossos_axis"] = {
        "varpi_mean_deg": round(mu_d, 1),
        "R": round(float(circ_R(vd)), 3),
        "p_rayleigh": rayleigh_p(vd),
        "sbdb_cluster_axis_deg": 49.9,
        "note": "OSSOS deep blocks concentrated near lam~210-240, "
                "the anti-cluster side"}
    print(f"B2 OSSOS detached varpi mean={mu_d:.0f} deg "
          f"R={out['B2_ossos_axis']['R']} "
          f"p={out['B2_ossos_axis']['p_rayleigh']:.3f}")

    # footprint: discovery-longitude distribution of all objects
    out["B2_footprint"] = {
        "lam_d_hist_edges": np.arange(0, 361, 30).tolist(),
        "lam_d_hist": np.histogram(np.rad2deg(lam_d),
                                 bins=np.arange(0, 361, 30))[0].tolist()}

    # ---------- B3 clustered objects' discovery geometry ----------
    mu = np.deg2rad(49.9)
    cl_mask = np.abs((vd - mu + np.pi) % (2 * np.pi) - np.pi) < \
        np.deg2rad(60)
    out["B3_clustered_discoveries"] = {
        "n_within60": int(cl_mask.sum()),
        "lam_d_deg": sorted(round(float(np.rad2deg(x)), 1)
                            for x in ld[cl_mask]),
        "lam_d_spread_deg": round(float(np.rad2deg(
            np.max(ld[cl_mask]) - np.min(ld[cl_mask]))), 1)
        if cl_mask.sum() else None}
    print(f"B3 clustered (varpi within 60 of 50): "
          f"{cl_mask.sum()} objects discovered at lam_d "
          f"{out['B3_clustered_discoveries']['lam_d_deg']}")

    # ---------- B4 discovery distance ----------
    dd = dist[det]
    qd = q[det]
    r_over_q = dd / qd   # ~1 near perihelion, ~e-fold far otherwise
    out["B4_distance"] = {
        "r_over_q_median": round(float(np.median(r_over_q)), 2),
        "r_over_q_p84": round(float(np.percentile(r_over_q, 84)), 2),
        "note": "r/q~1 means found near perihelion; >>1 means "
                "found far from perihelion (bias-immune direction)"}
    print(f"B4 discovery r/q median={np.median(r_over_q):.2f}")

    # ---------- B5 nodal-bias audit across the full catalog ----------
    # Omega (node) non-uniformity is a known ecliptic-discovery
    # artifact: objects are found near their nodes, so Omega tracks
    # the survey footprint in EVERY population.  The decisive
    # comparison: if pointing bias alone produced the a>150 varpi
    # cluster, the same footprint should sculpt varpi in the much
    # larger lower-a samples.  It does not -- varpi is uniform for
    # a<150 while the a>150 cluster stands.
    sbdb = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    rows, flds = sbdb["data"], sbdb["fields"]
    sob = []
    for r in rows:
        o = dict(zip(flds, r))
        try:
            for k in ("a", "e", "i", "om", "w", "q"):
                o[k] = float(o[k])
            o["cc"] = int(o["condition_code"]) \
                if o["condition_code"] is not None else 9
            o["arc"] = float(o["data_arc"]) if o["data_arc"] else 0
        except (TypeError, ValueError):
            continue
        sob.append(o)
    b5 = {}
    for lab, lo, hi, qm in [("a30_40", 30, 40, 30),
                            ("a40_60", 40, 60, 30),
                            ("a60_150", 60, 150, 30),
                            ("a150_detached", 150, 1e9, 30)]:
        ss = [o for o in sob if lo < o["a"] <= hi and o["q"] > qm
              and o["cc"] <= 3 and o["arc"] > 200]
        if len(ss) < 10:
            continue
        omv = [o["om"] for o in ss]
        vv = [(o["om"] + o["w"]) % 360 for o in ss]
        b5[lab] = {"N": len(ss),
                   "p_Omega": rayleigh_p(np.deg2rad(omv)),
                   "p_varpi": rayleigh_p(np.deg2rad(vv))}
        print(f"B5 {lab:16s} N={len(ss):5d}  p(Om)={b5[lab]['p_Omega']:.2e}"
              f"  p(varpi)={b5[lab]['p_varpi']:.2e}")
    out["B5_nodal_audit"] = b5
    out["B5_note"] = ("Omega is non-uniform at EVERY a (nodal "
                      "discovery bias); varpi is uniform for a<150 "
                      "in the same catalog under the same footprint "
                      "-- the a>150 varpi cluster is "
                      "population-specific, not a generic artifact")

    # per-object table
    ids = np.array([str(x) for x in t["ID"]])
    out["detached_table"] = [
        {"id": ids[i], "a": round(float(a[det][i]), 1),
         "q": round(float(qd[i]), 1),
         "varpi": round(float(np.rad2deg(vd[i])) % 360, 1),
         "lam_d": round(float(np.rad2deg(ld[i])), 1),
         "dphi": round(float(np.rad2deg(dphi[i])), 1),
         "dist_AU": round(float(dd[i]), 1),
         "r_over_q": round(float(r_over_q[i]), 2)}
        for i in range(det.sum())]

    RES.mkdir(exist_ok=True)
    (RES / "step_06_ossos_bias.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    ax[0].scatter(np.rad2deg(ld), np.rad2deg(vd), s=30,
                  c="steelblue", edgecolor="k", lw=0.4)
    ax[0].plot([0, 360], [0, 360], "r:", lw=1, label="varpi = lam_d")
    ax[0].axhline(49.9, color="k", ls="--", lw=0.8,
                  label="global cluster axis")
    ax[0].set(xlabel="discovery ecliptic longitude (deg)",
              ylabel="varpi (deg)",
              title="B1: pointing coupling, OSSOS detached")
    ax[0].legend(fontsize=7)

    ax[1].hist(np.rad2deg(lam_d), bins=np.arange(0, 361, 30),
               color="steelblue", edgecolor="k", lw=0.4)
    ax[1].axvline(49.9, color="r", ls="--", lw=1,
                  label="cluster axis")
    ax[1].axvline(229.9, color="g", ls="--", lw=1,
                  label="anti-axis")
    ax[1].set(xlabel="discovery ecliptic longitude (deg)",
              ylabel="N (all OSSOS)", title="B2: OSSOS footprint")
    ax[1].legend(fontsize=7)

    ax[2].hist(np.rad2deg(vd), bins=np.arange(0, 361, 30),
               color="steelblue", edgecolor="k", lw=0.4,
               label=f"OSSOS detached (N={det.sum()})")
    ax[2].axvline(49.9, color="r", ls="--", lw=1)
    ax[2].axvline(229.9, color="g", ls="--", lw=1)
    ax[2].set(xlabel="varpi (deg)", ylabel="N",
              title="B2: OSSOS detached varpi")
    ax[2].legend(fontsize=7)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_06_ossos_bias.png", dpi=300)
    logger.data_save(RESULTS / "step_06_ossos_bias.json")
    logger.data_save(RESULTS / "figures/supplementary/step_06_ossos_bias.png")


if __name__ == "__main__":
    main()
