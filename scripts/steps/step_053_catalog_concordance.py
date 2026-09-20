#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b18: catalog concordance
==================================================

Every clustering result so far uses JPL Small-Body Database
orbit solutions.  An independent orbit pipeline exists: the
Minor Planet Center's MPCORB.DAT -- the same objects, fitted
by a different group with different software and different
arc weighting.  If the perihelion cluster were an orbit-
solver systematic (epoch handling, planetary-ephemeris
version, outlier rejection), it need not reproduce on the
MPC solutions.

C1  Match the 44 detached objects to MPCORB by readable
    designation; element-by-element concordance (a, q, om,
    w, varpi).
C2  Clustering on the MPC elements alone: R(varpi), mean
    direction, Rayleigh p -- does the anomaly survive the
    change of orbit solver?
C3  Any systematic direction offset JPL-vs-MPC.

Data: data/sbdb_outer_ss.json (JPL), data/MPCORB.DAT.gz.

Outputs: results/step_b18_catalog_concordance.json,
         figures/supplementary/step_b18_catalog_concordance.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_053_catalog_concordance")
tee_stdout(logger)
logger.header("MPC catalogue concordance")

from pathlib import Path
import gzip
import json
import re
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW


def cmean(th):
    z = np.exp(1j * np.deg2rad(th)).mean()
    return float(np.abs(z)), float(np.rad2deg(np.angle(z)) % 360)


def rayleigh_p(R, N):
    return float(np.exp(-N * R * R))


def key_of(name):
    """Extract matchable key from an SBDB full_name.

    Numbered objects appear either as '541132 Leleakuhonua
    (2015 TG387)' (leading number) or '(541132)' (paren
    number); unnumbered as '2015 TG387'."""
    name = str(name)
    m = re.match(r"\s*(\d+)\s", name)
    if m:
        return ("num", int(m.group(1)))
    m = re.search(r"\((\d+)\)", name)
    if m:
        return ("num", int(m.group(1)))
    m = re.search(r"(\d{4})\s*([A-Z]{2}\d*)", name)
    if m:
        return ("prov", m.group(0).replace(" ", ""))
    return None


def main():
    # ---------- JPL detached sample ----------
    d = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    jpl = {}
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            a, q, cc = float(o["a"]), float(o["q"]), \
                int(o["condition_code"] or 9)
            om, w = float(o["om"]), float(o["w"])
        except (TypeError, ValueError):
            continue
        if a > 150 and q > 30 and cc <= 3:
            k = key_of(o["full_name"])
            if k:
                jpl[k] = {"name": o["full_name"], "a": a, "q": q,
                          "om": om, "w": w,
                          "varpi": (om + w) % 360}
    print(f"JPL detached sample: {len(jpl)}")

    # ---------- MPCORB: collect outer objects ----------
    mpc = {}
    with gzip.open((DATA_RAW / "mpc" / "MPCORB.DAT.gz"), "rt") as f:
        for line in f:
            if line.startswith("-") or len(line) < 160:
                continue
            try:
                # after the 5-char packed epoch, the next 7
                # whitespace-separated fields are
                # M, w, Om, i, e, n, a
                tk = line[25:103].split()
                M_, w, om, i_, e, n_, a = (float(t) for t in tk[:7])
            except (ValueError, IndexError):
                continue
            q = a * (1 - e)
            if not (a > 150 and q > 30):
                continue
            readable = line[166:194].strip()
            m = re.match(r"\((\d+)\)", readable)
            if m:
                k = ("num", int(m.group(1)))
            else:
                k = ("prov", readable.replace(" ", ""))
            mpc[k] = {"name": readable, "a": a, "q": q,
                      "om": om % 360, "w": w % 360,
                      "varpi": (om + w) % 360}
    print(f"MPCORB outer detached (a>150,q>30): {len(mpc)}")

    # ---------- match ----------
    both = sorted(set(jpl) & set(mpc))
    out = {"n_jpl": len(jpl), "n_mpc_outer": len(mpc),
           "n_matched": len(both)}
    print(f"matched: {len(both)}")

    if len(both) < 20:
        print("too few matches -- aborting")
        return

    ja = np.array([jpl[k]["a"] for k in both])
    ma = np.array([mpc[k]["a"] for k in both])
    jq = np.array([jpl[k]["q"] for k in both])
    mq = np.array([mpc[k]["q"] for k in both])
    jv = np.array([jpl[k]["varpi"] for k in both])
    mv = np.array([mpc[k]["varpi"] for k in both])
    dvp = (jv - mv + 180) % 360 - 180

    # ---------- C1 concordance ----------
    out["C1_concordance"] = {
        "med_abs_d_varpi_deg": round(float(np.median(np.abs(dvp))), 4),
        "p90_abs_d_varpi_deg": round(
            float(np.percentile(np.abs(dvp), 90)), 4),
        "med_da_au": round(float(np.median(ma - ja)), 3),
        "med_abs_da_au": round(float(np.median(np.abs(ma - ja))), 3),
        "med_dq_au": round(float(np.median(np.abs(mq - jq))), 3),
        "note": ("element concordance between JPL and MPC "
                 "orbit solutions for the same objects")}
    print(f"C1: med|d varpi|={np.median(np.abs(dvp)):.2f} "
          f"p90={np.percentile(np.abs(dvp), 90):.1f} "
          f"med|da|={np.median(np.abs(ma-ja)):.3f} AU")

    # ---------- C2 clustering on MPC elements ----------
    Rj, muj = cmean(jv)
    Rm, mum = cmean(mv)
    out["C2_clustering_mpc"] = {
        "jpl": {"R": round(Rj, 3), "mu": round(muj, 1),
                "rayleigh_p": f"{rayleigh_p(Rj, len(jv)):.2e}"},
        "mpc": {"R": round(Rm, 3), "mu": round(mum, 1),
                "rayleigh_p": f"{rayleigh_p(Rm, len(mv)):.2e}"},
        "axis_offset_deg": round(
            float(np.abs((mum - muj + 180) % 360 - 180)), 1),
        "note": ("the anomaly recomputed on MPC orbit "
                 "solutions for the identical object list.  "
                 "Caveat: element-for-element identity to "
                 "printed precision indicates the catalogs "
                 "share fitted solutions (JPL ingests MPC "
                 "fits) rather than being independent fits -- "
                 "this excludes solver/parse systematics on "
                 "our side but is not a fully independent "
                 "orbit determination")}
    print(f"C2: JPL R={Rj:.3f}@{muj:.0f} | MPC R={Rm:.3f}@"
          f"{mum:.0f} (offset {out['C2_clustering_mpc']['axis_offset_deg']} deg)")

    # ---------- C3 systematic offset ----------
    Rz, muz = cmean(dvp)
    out["C3_systematic_offset"] = {
        "mean_d_varpi_deg": round(muz, 2),
        "R_of_offsets": round(Rz, 3),
        "note": ("if orbit-solver systematics displaced varpi "
                 "coherently, the offset distribution would "
                 "be narrow and nonzero")}
    print(f"C3: mean d varpi = {muz:.2f} deg (R={Rz:.3f})")

    RES.mkdir(exist_ok=True)
    (RES / "step_b18_catalog_concordance.json").write_text(
        json.dumps(out, indent=1))

    # ---------- figure ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
    ax[0].hist(np.abs(dvp), bins=30, color="steelblue")
    ax[0].set(xlabel="|varpi_JPL - varpi_MPC| (deg)", ylabel="N",
              title="C1: element concordance")
    ax[0].axvline(np.median(np.abs(dvp)), color="r", ls="--",
                  label=f"med {np.median(np.abs(dvp)):.1f}")
    ax[0].legend(fontsize=8)

    ax[1].hist(jv, bins=np.arange(0, 361, 20), color="steelblue",
               alpha=0.6, label="JPL")
    ax[1].hist(mv, bins=np.arange(0, 361, 20), histtype="step",
               color="crimson", lw=1.8, label="MPC")
    ax[1].axvline(49, color="r", ls="--", lw=1)
    ax[1].set(xlabel="varpi (deg)", ylabel="N",
              title="C2: cluster on two catalogs")
    ax[1].legend(fontsize=8)

    ax[2].scatter(jv, mv, s=14, color="steelblue")
    ax[2].plot([0, 360], [0, 360], "r--", lw=1)
    ax[2].set(xlabel="varpi JPL (deg)", ylabel="varpi MPC (deg)",
              title="C3: per-object agreement")
    ax[2].set_xlim(0, 360); ax[2].set_ylim(0, 360)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b18_catalog_concordance.png",
                dpi=300)
    logger.data_save(RESULTS / "step_b18_catalog_concordance.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b18_catalog_concordance.png")


if __name__ == "__main__":
    main()
