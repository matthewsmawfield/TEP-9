#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b19: prospective predictions
======================================================

Everything so far is retrospective.  A rigorous program
should also freeze what the structure predicts BEFORE the
next discoveries arrive.  This step computes dated,
quantitative predictions from the current catalogs for
three future data releases:

P1  The next detached-class discoveries (a>150, q>30).
    From the b13 mixture fit (cluster fraction f, von Mises
    width sigma): predicted fraction of NEW detached
    objects within a 60 deg cap of the axis -- compared with
    the footprint-only baseline (what pure pointing would
    predict given the measured detached discovery
    footprint).

P2  The axis direction on new data: predicted mean varpi
    of future discoveries = 49 deg with the bootstrap CI
    from b13; the footprint-only baseline predicts the
    footprint mean (~10 deg).

P3  The perturber cross-check: under the conventional
    anti-aligned-confinement reading, the cluster's
    existence would place a shepherding perturber's own
    perihelion near ~229 deg.  That is where the shallow
    JFCs cluster (b17) -- recorded as the discriminator
    between the two hypotheses: future surveys covering
    that region either find the mass or they do not.

P4  The comet channel: fresh long-period comets'
    reconstruction discrepancies (the numbered-pipeline
    observable) should concentrate in the ~60 deg patch
    around the axis -- the pre-declared cap and expected
    excess are frozen here.

Each prediction is written to the JSON with the sample it
applies to, the numerical value, and the competing-model
baseline -- a dated, falsifiable record.

Outputs: results/step_b19_prospective_predictions.json,
         figures/supplementary/step_b19_prospective_predictions.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_054_prospective_predictions")
tee_stdout(logger)
logger.header("Dated falsifiable predictions")

from pathlib import Path
import json
import re
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260927)
AXIS = 49.0
ANTI = 229.0
CAP = 60.0

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


def cap_frac(th, center, cap=CAP):
    return float((np.abs((th - center + 180) % 360 - 180)
                  < cap).mean())


def vm_cap_frac(sigma_deg, cap=CAP):
    """fraction of a wrapped-normal(0,sigma) within +-cap."""
    from scipy.stats import norm
    return float(norm.cdf(cap / sigma_deg)
                 - norm.cdf(-cap / sigma_deg))


def main():
    d = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    vp, lo = [], []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            a, q, cc = float(o["a"]), float(o["q"]), \
                int(o["condition_code"] or 9)
            om, w = float(o["om"]), float(o["w"])
        except (TypeError, ValueError):
            continue
        if not (a > 150 and q > 30 and cc <= 3):
            continue
        lam = np.nan
        m = DESIG.search(str(o["full_name"]))
        if m and m.group(2) in HALF:
            yr, half = int(m.group(1)), m.group(2)
            mo, dy = HALF[half]
            lam = (sun_ecl_lon(yr, mo, dy) + 180) % 360
        vp.append((om + w) % 360)
        lo.append(lam)
    vp = np.array(vp); lo = np.array(lo)
    out = {"axis_deg": AXIS, "cap_deg": CAP, "n_detached": len(vp)}

    # ---------- P1 predicted in-cap fraction ----------
    f, sig = 0.43, 34.0          # b13 mixture fit
    pred = f * vm_cap_frac(sig) + (1 - f) * (2 * CAP / 360)
    obs = cap_frac(vp, AXIS)
    # footprint-only baseline: each object's own lam_opp
    # smeared by a gaussian of the OSSOS-measured coupling
    # width (~50 deg, conservative -- the measured near-
    # perihelion median |vp-lam_d| is ~44 deg); the fraction
    # of that smear falling inside the cap.
    from scipy.stats import norm
    d_ax = np.abs((lo - AXIS + 180) % 360 - 180)
    ok = np.isfinite(d_ax)
    base = float(np.mean(
        norm.cdf((CAP - d_ax[ok]) / 50)
        + norm.cdf((CAP + d_ax[ok]) / 50) - 1))
    out["P1_next_detached"] = {
        "sample": "next a>150, q>30 discoveries (any cc)",
        "predicted_frac_in_cap": round(pred, 2),
        "observed_current": round(obs, 2),
        "footprint_only_baseline": round(base, 2),
        "falsifier": ("if the next ~20 detached discoveries "
                      "land in-cap at the baseline rate "
                      "(~%.0f%%) or below, the intrinsic-axis "
                      "interpretation loses its strongest "
                      "prospective support" % (base * 100)),
        "note": ("mixture prediction = f*P_vm(cap) + "
                 "(1-f)*uniform; baseline = each object's "
                 "footprint position smeared by the "
                 "OSSOS-measured coupling at 100%")}
    print(f"P1: predicted {pred:.2f} vs baseline {base:.2f} "
          f"(observed {obs:.2f})")

    # ---------- P2 direction ----------
    n_mc = 4000
    boot = np.empty(n_mc)
    for k in range(n_mc):
        s = rng.choice(vp, size=len(vp), replace=True)
        boot[k] = np.rad2deg(np.angle(
            np.exp(1j * np.deg2rad(s)).mean())) % 360
    lo_ci = np.percentile(np.abs((boot - AXIS + 180) % 360
                                 - 180) * np.sign(
        (boot - AXIS + 180) % 360 - 180), [2.5, 97.5]) + AXIS
    out["P2_direction"] = {
        "predicted_mean_varpi_deg": AXIS,
        "ci_95_deg": [round(float(lo_ci[0]), 0),
                      round(float(lo_ci[1]), 0)],
        "footprint_only_baseline_deg": 10.5,
        "falsifier": ("a future detached subsample whose mean "
                      "varpi lands outside the CI while its "
                      "footprint stays near ~10 deg would "
                      "contradict the axis being a fixed "
                      "inertial structure")}
    print(f"P2: direction 49 CI {lo_ci}")

    # ---------- P3 perturber cross-check ----------
    out["P3_perturber_crosscheck"] = {
        "if_point_mass": ("anti-aligned confinement places the "
                          "perturber's perihelion near ~229 deg "
                          "-- coincident with the shallow-JFC "
                          "concentration found in b17"),
        "discriminator": ("deep surveys covering lam ~ 229 deg, "
                          "beta ~ +17 deg either find the mass "
                          "or they do not; under the field "
                          "interpretation that region is the "
                          "outflow end of the boundary, under "
                          "the mass interpretation it is the "
                          "perturber's own apsidal direction"),
        "note": ("the two hypotheses make different "
                 "predictions for that sky region: a detectable "
                 "moving point source vs nothing")}

    # ---------- P4 comet channel ----------
    # compute the real CODE-measured rates rather than quoting text
    from scripts.utils.tep9_common import parse_code
    import math
    orig = parse_code(DATA_RAW / "code" / "code_original.html")
    fut = parse_code(DATA_RAW / "code" / "code_future.html")
    warsaw = {l[5:17].strip()
              for l in open(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
              if len(l) > 115}

    def _pdir(w, Om, i):
        w, Om, i = map(math.radians, (w, Om, i))
        return np.array([math.cos(Om) * math.cos(w)
                         - math.sin(Om) * math.sin(w) * math.cos(i),
                         math.sin(Om) * math.cos(w)
                         + math.cos(Om) * math.sin(w) * math.cos(i),
                         math.sin(w) * math.sin(i)])

    def _sep(a, b):
        return math.degrees(
            math.acos(float(np.clip(np.dot(a, b), -1, 1))))

    _l, _b = math.radians(34.0), math.radians(-13.0)
    _axis = np.array([math.cos(_b) * math.cos(_l),
                      math.cos(_b) * math.sin(_l), math.sin(_b)])
    din, dout = [], []
    for k, r in orig.items():
        if k not in fut or k in warsaw:
            continue
        if not (0 < r["aa"] < 100 and r["q"] < 3.1
                and r["cls"] in ("1a", "1a+", "1b")):
            continue
        aph = -_pdir(r["w"], r["Om"], r["i"])
        pf = -_pdir(fut[k]["w"], fut[k]["Om"], fut[k]["i"])
        (din if _sep(aph, _axis) < CAP else dout).append(
            _sep(aph, pf))
    din, dout = np.array(din), np.array(dout)
    thr = 0.20  # deg -- predeclared large-discrepancy threshold
    out["P4_comet_channel"] = {
        "sample": "fresh long-period comets entering catalog",
        "predeclared_cap_deg": 60,
        "predeclared_axis": "lam ~ 34, beta ~ -13 (comet-side)",
        "n_code_only_class1": int(len(din) + len(dout)),
        "n_in": int(len(din)), "n_out": int(len(dout)),
        "frac_cap_of_sample": float(len(din) / (len(din) + len(dout))),
        "rate_gt_0p2_in": float((din > thr).mean()),
        "rate_gt_0p2_out": float((dout > thr).mean()),
        "med_in_deg": float(np.median(din)),
        "med_out_deg": float(np.median(dout)),
        "expected": ("reconstruction discrepancies "
                     "concentrated in the cap: the "
                     "CODE-measured rate of d(orig->fut) > 0.2 "
                     "deg is "
                     f"{(din > thr).mean() * 100:.0f}% in-cap vs "
                     f"{(dout > thr).mean() * 100:.0f}% out-cap, "
                     "energy channel flat"),
        "falsifier": ("a fresh comet cohort with flat "
                      "rotation-channel vs the cap would "
                      "undercut the clock-channel claim")}
    RES.mkdir(exist_ok=True)
    (RES / "step_b19_prospective_predictions.json").write_text(
        json.dumps(out, indent=1))

    # ---------- figure ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(9, 4.2))
    ax[0].bar(["footprint\nbaseline", "mixture\nprediction",
               "observed\n(current)"],
              [base, pred, obs],
              color=["steelblue", "crimson", "0.4"])
    ax[0].axhline(2 * CAP / 360, color="0.5", ls=":",
                  label="uniform")
    ax[0].set(ylabel="fraction within 60 deg of axis",
              title="P1: next detached\ndiscoveries")
    ax[0].legend(fontsize=8)
    ax[0].set_ylim(0, 1)

    ax[1].hist(boot, bins=40, color="steelblue", alpha=0.7,
               label="bootstrap mean varpi")
    ax[1].axvline(AXIS, color="r", lw=1.5, label="axis 49")
    ax[1].axvline(10.5, color="purple", ls=":", lw=1.5,
                  label="footprint mean")
    ax[1].set(xlabel="mean varpi (deg)", ylabel="count",
              title="P2: direction prediction")
    ax[1].legend(fontsize=8)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b19_prospective_predictions.png",
                dpi=300)
    logger.data_save(RESULTS / "step_b19_prospective_predictions.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b19_prospective_predictions.png")


if __name__ == "__main__":
    main()
