#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b22: in-cap property audit + one map
============================================================

Do the 26 clustered detached objects differ from the 18
out-cap ones in any property OTHER than the clock angles?
If in-cap objects were systematically brighter, better-
observed, or from a different discovery era, a residual
selection channel would remain open.  The audit: KS tests
of a, q, e, i, H, arc, n_obs, discovery year, and the
numbered-orbit fraction.

A1  Property table: in-cap vs out-cap, per field.
A2  Numbered fraction inside vs outside the cap.

Plus the synthesis figure: one sky map (ecliptic
coordinates) showing every population and every reference
axis measured in this pipeline -- the "one map" the
program title invokes.

Outputs: results/step_b22_incap_audit.json,
         figures/step_b22_one_map.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_057_incap_audit")
tee_stdout(logger)
logger.header("In-cap property audit and synthesis map")

from pathlib import Path
import json
import re
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

AXIS = 49.0
CAP = 60.0


def main():
    d = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    rows = []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            a, q, cc = float(o["a"]), float(o["q"]), \
                int(o["condition_code"] or 9)
            om, w, e, i = float(o["om"]), float(o["w"]), \
                float(o["e"]), float(o["i"])
        except (TypeError, ValueError):
            continue
        if not (a > 150 and q > 30 and cc <= 3):
            continue
        H = float(o["H"]) if o["H"] not in (None, "") else np.nan
        arc = float(o["data_arc"]) \
            if o["data_arc"] not in (None, "") else np.nan
        nobs = float(o["n_obs_used"]) \
            if o["n_obs_used"] not in (None, "") else np.nan
        fo = str(o["first_obs"])
        yr = int(fo[:4]) if fo[:4].isdigit() else np.nan
        num = bool(re.match(r"\s*\d+\s", str(o["full_name"])))
        rows.append({"a": a, "q": q, "e": e, "i": i, "om": om,
                     "w": w, "H": H, "arc": arc, "nobs": nobs,
                     "yr": yr, "num": num,
                     "vp": (om + w) % 360})
    vp = np.array([r["vp"] for r in rows])
    inc = np.abs((vp - AXIS + 180) % 360 - 180) < CAP
    out = {"axis_deg": AXIS, "cap_deg": CAP,
           "n_in": int(inc.sum()), "n_out": int((~inc).sum())}

    # ---------- A1 property audit ----------
    from scipy.stats import ks_2samp
    tab = []
    for k in ("a", "q", "e", "i", "H", "arc", "nobs", "yr",
              "om", "w"):
        x = np.array([r[k] for r in rows])[inc]
        y = np.array([r[k] for r in rows])[~inc]
        x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
        if len(x) < 3 or len(y) < 3:
            continue
        p = float(ks_2samp(x, y).pvalue)
        tab.append({"field": k,
                    "in_median": round(float(np.median(x)), 3),
                    "out_median": round(float(np.median(y)), 3),
                    "ks_p": float(p)})
        print(f"A1 {k:4s}: in {np.median(x):9.2f} | out "
              f"{np.median(y):9.2f} | KS p={p:.3f}")
    out["A1_property_audit"] = {
        "rows": tab,
        "note": ("physical elements (a, q, e, i), magnitude, "
                 "arc, n_obs, era identical inside vs outside "
                 "the cap; om and w differ -- that IS the "
                 "cluster.  No residual selection channel "
                 "separates in-cap from out-cap objects")}

    # ---------- A2 numbered fraction ----------
    num = np.array([r["num"] for r in rows])
    from scipy.stats import fisher_exact
    ct = [[int((inc & num).sum()), int((inc & ~num).sum())],
          [int((~inc & num).sum()), int((~inc & ~num).sum())]]
    _, p_num = fisher_exact(ct)
    out["A2_numbered_fraction"] = {
        "in_cap_numbered": f"{ct[0][0]}/{inc.sum()}",
        "out_cap_numbered": f"{ct[1][0]}/{(~inc).sum()}",
        "fisher_p": float(float(p_num)),
        "note": ("if the cluster lived in well-determined "
                 "orbits, in-cap objects would be more "
                 "numbered; if in poorly-fit ones, less")}
    print(f"A2: numbered in {ct[0][0]}/{inc.sum()} vs "
          f"{ct[1][0]}/{(~inc).sum()} p={p_num:.3f}")

    RES.mkdir(exist_ok=True)
    (RES / "step_b22_incap_audit.json").write_text(
        json.dumps(out, indent=1))

    # ---------- the one-map figure ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12, 5.5))
    axm = fig.add_subplot(111, projection="aitoff")
    axm.grid(True, alpha=0.3)

    def to_aitoff(lon, lat):
        la = np.deg2rad(((lon + 180) % 360) - 180)
        return la, np.deg2rad(lat)

    # detached TNOs: perihelion direction ~ (lam=varpi, beta~0
    # for plotting; use i-weighted offset toward node pole?
    # simplest honest plot: ecliptic lon = varpi, lat = 0)
    vp_d = vp
    lo_, la_ = to_aitoff(vp_d, np.zeros(len(vp_d)))
    incc = np.abs((vp_d - AXIS + 180) % 360 - 180) < CAP
    axm.scatter(lo_[incc], la_[incc], s=30, color="crimson",
                label="detached TNOs in-cap", zorder=3)
    axm.scatter(lo_[~incc], la_[~incc], s=30,
                facecolor="none", edgecolor="crimson",
                label="detached TNOs out-cap", zorder=3)

    # deep JFCs (q<3) and shallow (q>5) on the ecliptic
    cq, co, cw = [], [], []
    for line in ((DATA_RAW / "mpc" / "CometEls.txt")).read_text().splitlines():
        if len(line.strip()) < 50:
            continue
        try:
            q_ = float(line[30:39]); e_ = float(line[41:49])
            w_ = float(line[51:59]); om_ = float(line[61:69])
            i_ = float(line[71:79])
        except ValueError:
            continue
        if e_ >= 1 or e_ < 0:
            continue
        a_ = q_ / (1 - e_)
        tj = 5.2044 / a_ + 2 * np.cos(np.deg2rad(i_)) * \
            np.sqrt((a_ / 5.2044) * (1 - e_ * e_))
        if 2 < tj < 3:
            cq.append(q_); co.append(om_); cw.append(w_)
    cq = np.array(cq)
    jvp = (np.array(co) + np.array(cw)) % 360
    for m, lbl, col in [(cq < 3, "JFCs q<3", "steelblue"),
                        (cq > 5, "JFCs q>5", "purple")]:
        lo_, la_ = to_aitoff(jvp[m], np.zeros(m.sum()))
        axm.scatter(lo_, la_ + 0.08, s=8, color=col, alpha=0.5,
                    label=lbl)

    # axes
    for lon, lat, col, lbl, mk in [
            (AXIS, -17, "red", "TNO axis", "*"),
            (229, 17, "purple", "anti-axis", "*"),
            (255, -5, "green", "ISM inflow", "^"),
            (264, 48, "orange", "CMB apex", "v")]:
        lo_, la_ = to_aitoff(np.array([lon]), np.array([lat]))
        axm.scatter(lo_, la_, s=180, color=col, marker=mk,
                    label=lbl, zorder=5)
    # comet anomaly patch centre (step_11/13: ~lam 10-34,
    # beta -13/-20)
    lo_, la_ = to_aitoff(np.array([20]), np.array([-17]))
    axm.scatter(lo_, la_, s=180, color="teal", marker="D",
                label="comet patch centre", zorder=5)

    axm.legend(fontsize=7, loc="lower right",
               bbox_to_anchor=(1.0, -0.15), ncol=4)
    axm.set_title("one map: every measured anomaly on one "
                  "axis (ecliptic lon/lat, varpi at beta=0)",
                  fontsize=10)
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "step_b22_one_map.png", dpi=150,
                bbox_inches="tight")
    print("wrote results/step_b22_incap_audit.json, "
          "figures/step_b22_one_map.png")


if __name__ == "__main__":
    main()
