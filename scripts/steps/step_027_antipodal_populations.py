#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b12: the antipodal-population test
=========================================================

Step_03 found the outer-comet perihelion anisotropy near
lam ~ 219 deg and the detached-TNO cluster near lam ~ 49 deg --
the two ends of ONE axis.  A directed field structure predicts
exactly this bipolar pattern on different dynamical
populations; a point mass does not.  But comet perihelion
statistics are discovery-biased, so the JFC sample gets the
full bias machinery developed for the TNOs:

B1  Reconstruct the JFC discovery footprint (provisional
    designation -> lam_opp) and measure the realized coupling
    Delta = varpi - lam_opp.

B2  The pairing test: if JFC varpi clustering were pure
    footprint smear, the observed (lam_opp, Delta) pairing is
    a typical shuffle of the marginals.  Under intrinsic
    clustering the pairing compensates -- same machinery as
    b8/D2.

B3  Ceiling: R(varpi) vs R(lam_opp) -- can the realized
    footprint reach the observed concentration at all?

B4  Era and quality stability: designation-era splits,
    and cross-checks in the CTc and HYP classes.

Outputs: results/step_b12_antipodal_populations.json,
         figures/supplementary/step_b12_antipodal_populations.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_027_antipodal_populations")
tee_stdout(logger)
logger.header("JFC bipolar population")

from pathlib import Path
import json
import re
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260924)
N_MC = 20000
AXIS, ANTI = 49.0, 229.0

HALF = {"A": (1, 8), "B": (1, 23), "C": (2, 8), "D": (2, 22),
        "E": (3, 8), "F": (3, 23), "G": (4, 8), "H": (4, 23),
        "J": (5, 8), "K": (5, 23), "L": (6, 8), "M": (6, 23),
        "N": (7, 8), "O": (7, 23), "P": (8, 8), "Q": (8, 23),
        "R": (9, 8), "S": (9, 23), "T": (10, 8), "U": (10, 23),
        "V": (11, 8), "W": (11, 23), "X": (12, 8), "Y": (12, 23)}
DESIG = re.compile(r"(\d{4})\s*([A-Z])([A-Z]?\d*)")


def sun_ecl_lon(year, month, day):
    from astropy.time import Time
    from astropy.coordinates import get_sun
    t = Time(f"{year:04d}-{month:02d}-{day:02d}T00:00:00",
             format="isot", scale="utc")
    return float(get_sun(t).geocentrictrueecliptic.lon.deg) % 360


def circ_R(th):
    return float(np.abs(np.exp(1j * np.asarray(th)).mean()))


def circ_mu(th):
    return float(np.rad2deg(
        np.angle(np.exp(1j * np.asarray(th)).mean())) % 360)


def load():
    d = json.loads(((DATA_RAW / "sbdb" / "sbdb_comets_ext.json")).read_text())
    pops = {}
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        cl = o.get("class")
        try:
            om, w = float(o["om"]), float(o["w"])
            q, a = float(o["q"]), float(o["a"])
        except (TypeError, ValueError):
            continue
        m = DESIG.search(str(o["full_name"]))
        yr = None
        lam_opp = None
        if m and m.group(2) in HALF:
            yr = int(m.group(1))
            mo, dy = HALF[m.group(2)]
            lam_opp = (sun_ecl_lon(yr, mo, dy) + 180) % 360
        pops.setdefault(cl, []).append(
            {"varpi": (om + w) % 360, "lam_opp": lam_opp,
             "yr": yr, "q": q, "a": a})
    return pops


def analyze(vp, lo, label, out):
    """B1-B3 for one population with parsed footprints."""
    vp, lo = np.deg2rad(vp), np.deg2rad(lo)
    n = len(vp)
    D = (vp - lo + np.pi) % (2 * np.pi) - np.pi
    R_vp, R_lo, R_D = circ_R(vp), circ_R(lo), circ_R(D)
    mu_vp, mu_lo = circ_mu(vp), circ_mu(lo)

    # pairing/shuffle null
    idx = rng.integers(0, n, (N_MC, n))
    vp_sh = (np.broadcast_to(lo, (N_MC, n)) + D[idx]) % (2 * np.pi)
    R_sh = np.abs(np.exp(1j * vp_sh).mean(axis=1))
    p_pair = float((int((R_sh >= R_vp).sum()) + 1) / (N_MC + 1))

    res = {"N": n,
           "R_varpi": round(R_vp, 3),
           "mu_varpi": round(mu_vp, 1),
           "R_lam_opp": round(R_lo, 3),
           "mu_lam_opp": round(mu_lo, 1),
           "R_Delta": round(R_D, 3),
           "R_indep_conv": round(R_lo * R_D, 3),
           "ceiling_R_lam_opp": round(R_lo, 3),
           "shuffle_null_R": round(float(R_sh.mean()), 3),
           "p_pairing": float(p_pair),
           "sep_to_anti_axis_deg": round(
               float(np.abs((mu_vp - ANTI + 180) % 360 - 180)), 1)}
    print(f"{label}: N={n} R(varpi)={R_vp:.3f}@{mu_vp:.0f} "
          f"R(lo)={R_lo:.3f}@{mu_lo:.0f} R(D)={R_D:.2f} "
          f"p_pair={p_pair:.4f} sep_to_229="
          f"{res['sep_to_anti_axis_deg']}")
    out[label] = res
    return res


def main():
    pops = load()
    out = {"axis_deg": AXIS, "anti_axis_deg": ANTI}

    # ---------------- JFC full analysis ----------------
    jfc = [r for r in pops.get("JFc", []) if r["lam_opp"] is not None]
    vp = np.array([r["varpi"] for r in jfc])
    lo = np.array([r["lam_opp"] for r in jfc])
    yr = np.array([r["yr"] for r in jfc])
    analyze(vp, lo, "JFc_all", out)

    # era splits
    for lab, m in [("JFc_pre2000", yr < 2000),
                   ("JFc_post2000", yr >= 2000)]:
        analyze(vp[m], lo[m], lab, out)

    # ---------------- CTc / HYP cross-checks ----------------
    for cl in ["CTc", "HYP", "COM"]:
        sub = [r for r in pops.get(cl, [])
               if r["lam_opp"] is not None]
        if len(sub) < 8:
            continue
        analyze(np.array([r["varpi"] for r in sub]),
                np.array([r["lam_opp"] for r in sub]), cl, out)

    out["note"] = ("the JFC cluster direction (~220 deg) sits "
                   "~10 deg off the anti-axis (229 deg).  "
                   "Whether the pairing test finds compensation "
                   "beyond the footprint marginals determines "
                   "if this is a bipolar population datum or a "
                   "pointing artifact; era stability is the "
                   "second check")

    RES.mkdir(exist_ok=True)
    (RES / "step_b12_antipodal_populations.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))

    ax[0].hist(np.rad2deg(np.deg2rad(lo)), bins=np.arange(0, 361, 20),
               color="steelblue", alpha=0.7, label="lam_opp footprint")
    ax[0].hist(np.rad2deg(np.deg2rad(vp)), bins=np.arange(0, 361, 20),
               histtype="step", color="crimson", lw=1.6, label="varpi")
    ax[0].axvline(ANTI, color="purple", ls=":", lw=1.2,
                  label="anti-axis 229")
    ax[0].axvline(AXIS, color="r", ls="--", lw=1)
    ax[0].set(xlabel="ecliptic longitude (deg)", ylabel="N",
              title="B1: JFC footprint vs varpi")
    ax[0].legend(fontsize=7)

    jres = out["JFc_all"]
    D = (np.deg2rad(vp) - np.deg2rad(lo) + np.pi) % (2 * np.pi) - np.pi
    idx = rng.integers(0, len(vp), (5000, len(vp)))
    vp_sh = (np.broadcast_to(np.deg2rad(lo), (5000, len(vp)))
             + D[idx]) % (2 * np.pi)
    R_sh = np.abs(np.exp(1j * vp_sh).mean(axis=1))
    ax[1].hist(R_sh, bins=40, density=True, color="0.6",
               label="shuffle null")
    ax[1].axvline(jres["R_varpi"], color="r", lw=1.8,
                  label=f"observed {jres['R_varpi']:.2f}")
    ax[1].axvline(jres["R_indep_conv"], color="purple", ls=":",
                  label=f"indep conv {jres['R_indep_conv']:.2f}")
    ax[1].set(xlabel="R(varpi)", ylabel="density",
              title=f"B2: JFC pairing test\n"
                    f"p={jres['p_pairing']:.4f}")
    ax[1].legend(fontsize=7)

    labels = [k for k in out if isinstance(out[k], dict)
              and "mu_varpi" in out[k]]
    mus = [out[k]["mu_varpi"] for k in labels]
    Rs = [out[k]["R_varpi"] for k in labels]
    ax[2].scatter(mus, Rs, s=50, c="steelblue", edgecolor="k")
    for k, x, y in zip(labels, mus, Rs):
        ax[2].annotate(k.replace("JFc_", ""), (x, y), fontsize=7,
                       xytext=(4, 4), textcoords="offset points")
    ax[2].axvline(ANTI, color="purple", ls=":", label="anti-axis")
    ax[2].axvline(AXIS, color="r", ls="--", label="axis")
    ax[2].set(xlabel="mean varpi (deg)", ylabel="R(varpi)",
              xlim=(0, 360),
              title="B4: all comet classes\nvs both axes")
    ax[2].legend(fontsize=7)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b12_antipodal_populations.png",
                dpi=300)
    logger.data_save(RESULTS / "step_b12_antipodal_populations.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b12_antipodal_populations.png")


if __name__ == "__main__":
    main()
