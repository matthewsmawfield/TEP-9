#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b16: the decorrelation radius
======================================================

Where does the intrinsic axis take over from the footprint?
For each semimajor-axis bin, count the fraction of objects
within +-60 deg of (a) the detached-cluster axis lam = 49 deg
and (b) the bin's OWN discovery-footprint mean direction.
Where pointing dominates, the footprint cap wins; where the
intrinsic structure dominates, the axis cap wins.  The
crossover radius is an empirical measurement of where the
boundary structure begins to organize the clock sector.

R1  Axis-score vs footprint-score in sliding a windows.

R2  Resonant control: Plutinos (3:2, a ~ 39.4 AU) have
    perihelion longitudes pinned by resonance geometry, not
    by any external axis or footprint.  Their axis-alignment
    fraction is the natural control for the inner belt's
    weak offset (b11).

R3  The sednoid subsample (a>150, q>50): the cleanest objects
    -- deepest boundary-residents -- for comparison.

Outputs: results/step_b16_decorrelation_radius.json,
         figures/supplementary/step_b16_decorrelation_radius.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_017_decorrelation_radius")
tee_stdout(logger)
logger.header("Axis-footprint decorrelation radius")

from pathlib import Path
import json
import re
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260927)
AXIS = 49.0

HALF = {"A": (1, 8), "B": (1, 23), "C": (2, 8), "D": (2, 22),
        "E": (3, 8), "F": (3, 23), "G": (4, 8), "H": (4, 23),
        "J": (5, 8), "K": (5, 23), "L": (6, 8), "M": (6, 23),
        "N": (7, 8), "O": (7, 23), "P": (8, 8), "Q": (8, 23),
        "R": (9, 8), "S": (9, 23), "T": (10, 8), "U": (10, 23),
        "V": (11, 8), "W": (11, 23), "X": (12, 8), "Y": (12, 23)}
DESIG = re.compile(r"\((\d{4})\s*([A-Z])([A-Z]?\d*)\)")


def sun_ecl_lon(year, month, day):
    from astropy.time import Time
    from astropy.coordinates import get_sun
    t = Time(f"{year:04d}-{month:02d}-{day:02d}T00:00:00",
             format="isot", scale="utc")
    return float(get_sun(t).geocentrictrueecliptic.lon.deg) % 360


def cap_frac(th, center, cap=60.0):
    return float((np.abs((th - center + 180) % 360 - 180)
                  < cap).mean())


def main():
    d = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    rows = []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            a, q, cc = float(o["a"]), float(o["q"]), \
                int(o["condition_code"] or 9)
            om, w = float(o["om"]), float(o["w"])
        except (TypeError, ValueError):
            continue
        if not (a > 30 and q > 30 and cc <= 3):
            continue
        lam_opp = np.nan
        m = DESIG.search(str(o["full_name"]))
        if m and m.group(2) in HALF:
            yr, half = int(m.group(1)), m.group(2)
            mo, dy = HALF[half]
            lam_opp = (sun_ecl_lon(yr, mo, dy) + 180) % 360
        rows.append({"a": a, "q": q, "varpi": (om + w) % 360,
                     "lam_opp": lam_opp})
    a = np.array([r["a"] for r in rows])
    q = np.array([r["q"] for r in rows])
    vp = np.array([r["varpi"] for r in rows])
    lo = np.array([r["lam_opp"] for r in rows])
    out = {"axis_deg": AXIS, "cap_deg": 60, "n": len(rows)}

    # ---------------- R1 axis vs footprint score ----------------
    # per-object discriminator: is varpi closer to the axis
    # or to the object's OWN discovery longitude?  Under pure
    # pointing the latter wins; under intrinsic structure the
    # former.  (A bin-mean cap comparison fails because the
    # detached footprint mean ~10 deg sits inside the axis
    # cap's overlap range.)
    bins = [(30, 40), (40, 50), (50, 70), (70, 100),
            (100, 150), (150, 300), (300, np.inf)]
    prof = []
    for lo_, hi_ in bins:
        m = (a > lo_) & (a <= hi_) & np.isfinite(lo)
        if m.sum() < 10:
            continue
        d_ax = np.abs((vp[m] - AXIS + 180) % 360 - 180)
        d_fp = np.abs((vp[m] - lo[m] + 180) % 360 - 180)
        f_ax = float((d_ax < d_fp).mean())
        prof.append({"bin": f"{lo_:.0f}-{hi_:.0f}",
                     "a_mid": round(float(np.median(a[m])), 1),
                     "N": int(m.sum()),
                     "frac_closer_to_axis": round(f_ax, 3),
                     "frac_closer_to_fp": round(1 - f_ax, 3),
                     "med_d_axis": round(float(np.median(d_ax)), 0),
                     "med_d_fp": round(float(np.median(d_fp)), 0)})
        print(f"R1 a {lo_:.0f}-{hi_:.0f}: N={m.sum():4d} "
              f"closer-to-axis {f_ax:.2f} "
              f"(med|vp-49|={np.median(d_ax):.0f} vs "
              f"med|vp-lam_opp|={np.median(d_fp):.0f})")
    out["R1_axis_vs_footprint"] = {
        "rows": prof,
        "note": ("fraction of objects whose varpi is closer to "
                 "the axis than to their own discovery "
                 "longitude.  The mid-belt dips to ~0.23 -- "
                 "footprint coupling genuinely dominates "
                 "there (med|varpi-lam_opp| ~ 31-38 deg).  "
                 "Beyond 150 AU the fraction rises to "
                 "~0.40-0.42 even though detached objects have "
                 "the TIGHTEST coupling of all (med 24-29 "
                 "deg) -- the axis beats the footprint more "
                 "often exactly where pointing is strongest.  "
                 "The turnover at ~150 AU is the decorrelation "
                 "signature: perihelion directions decouple "
                 "from discovery geometry at the boundary")}

    # ---------------- R2 resonant control ----------------
    # Plutinos: 3:2 at a~39.4; Twotinos: 2:1 at a~47.8
    plut = (np.abs(a - 39.4) < 0.6) & (q > 30)
    two = (np.abs(a - 47.8) < 0.7) & (q > 30)
    inner_nr = (a > 35) & (a < 42) & ~plut & (q > 30)
    out["R2_resonant_control"] = {
        "plutino": {"N": int(plut.sum()),
                    "axis_frac": round(cap_frac(vp[plut], AXIS), 3)},
        "twotino": {"N": int(two.sum()),
                    "axis_frac": round(cap_frac(vp[two], AXIS), 3)},
        "inner_nonresonant": {
            "N": int(inner_nr.sum()),
            "axis_frac": round(cap_frac(vp[inner_nr], AXIS), 3)},
        "note": ("resonance-locked perihelia are pinned by "
                 "Neptune, not by footprint or external axis -- "
                 "their axis-fraction is the control for the "
                 "inner belt's weak offset")}
    for k in ("plutino", "twotino", "inner_nonresonant"):
        print(f"R2 {k}: N={out['R2_resonant_control'][k]['N']} "
              f"axis_frac={out['R2_resonant_control'][k]['axis_frac']}")

    # ---------------- R3 sednoids ----------------
    sed = (a > 150) & (q > 50)
    det_ = (a > 150) & (q > 30)
    out["R3_sednoids"] = {
        "N": int(sed.sum()),
        "R_varpi": round(float(np.abs(
            np.exp(1j * np.deg2rad(vp[sed])).mean())), 3),
        "mean_varpi_deg": round(float(np.rad2deg(np.angle(
            np.exp(1j * np.deg2rad(vp[sed])).mean())) % 360), 1),
        "axis_frac": round(cap_frac(vp[sed], AXIS), 3),
        "detached_axis_frac": round(cap_frac(vp[det_], AXIS), 3),
        "note": ("the deepest boundary-residents: their "
                 "clustering strength and direction vs the "
                 "wider detached sample")}
    print(f"R3 sednoids N={sed.sum()} "
          f"R={out['R3_sednoids']['R_varpi']} "
          f"mu={out['R3_sednoids']['mean_varpi_deg']} "
          f"axis_frac={out['R3_sednoids']['axis_frac']}")

    RES.mkdir(exist_ok=True)
    (RES / "step_b16_decorrelation_radius.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))

    mids = [p["a_mid"] for p in prof]
    ax[0].plot(mids, [p["frac_closer_to_axis"] for p in prof],
               "o-", color="crimson",
               label="frac |varpi-49| < |varpi-lam_opp|")
    ax[0].axhline(0.5, color="0.5", ls=":", label="tie")
    ax[0].set(xscale="log", xlabel="a bin median (AU)",
              ylabel="fraction closer to axis",
              title="R1: decorrelation profile\n(axis beats "
                    "footprint where curve rises)")
    ax[0].legend(fontsize=7)

    ks = ["plutino", "twotino", "inner_nonresonant"]
    ax[1].bar(range(3),
              [out["R2_resonant_control"][k]["axis_frac"]
               for k in ks],
              color=["steelblue", "steelblue", "crimson"])
    ax[1].axhline(1 / 3, color="0.5", ls=":")
    ax[1].set_xticks(range(3))
    ax[1].set_xticklabels(["Plutino\n(3:2)", "Twotino\n(2:1)",
                           "inner\nnon-res"], fontsize=8)
    ax[1].set(ylabel="frac within 60 of axis",
              title="R2: resonant control")

    ax[2].hist(vp[det_], bins=np.arange(0, 361, 20),
               color="steelblue", alpha=0.6, label="detached")
    ax[2].hist(vp[sed], bins=np.arange(0, 361, 20),
               histtype="step", color="crimson", lw=1.8,
               label="sednoid q>50")
    ax[2].axvline(AXIS, color="r", ls="--", lw=1)
    ax[2].set(xlabel="varpi (deg)", ylabel="N",
              title="R3: sednoid subsample")
    ax[2].legend(fontsize=7)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b16_decorrelation_radius.png",
                dpi=300)
    logger.data_save(RESULTS / "step_b16_decorrelation_radius.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b16_decorrelation_radius.png")


if __name__ == "__main__":
    main()
