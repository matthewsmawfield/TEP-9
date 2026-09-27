#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b17: the injection chain
==================================================

The injected population is not one sample -- it is a chain:
detached TNOs (boundary-resident) -> centaurs (transit) ->
Jupiter-family comets (injected).  The full JFC population
(N ~ 580, T_J in 2-3) shows a systematic rotation of mean
longitude of perihelion with perihelion depth: ~+40 deg at
q < 3 -- aligned with the detached-TNO axis at 49 deg --
progressing to ~198 deg pooled at q > 5 and ~220 deg on the
dated JFc subset of step_b12, the anti-axis sector.

The injected chain therefore touches BOTH ends of the
measured axis.  Whether this progression is a bipolar
boundary signature or a JFC-specific secular/observational
feature is tested with the usual machinery: measured
discovery footprints per q-bin and the pairing test.

I1  The q-progression: mean varpi, Omega, w per JFC q-bin.
I2  Discovery footprints per q-bin (provisional-dated
    subset): does pointing explain the progression?
I3  Pairing test on the deepest injected subsample (q<3):
    observed (lam_opp, Delta) pairing vs its own marginals.
I4  The chain summary: detached TNOs, centaurs, deep JFCs,
    shallow JFCs.

Caveat carried forward from b12: JFC perihelia carry known
secular structure (Jupiter-driven omega preferences), so
amplitude claims are unsafe; the datum is the direction
progression and its endpoints.

Data: data/CometEls.txt (all MPC comets, q/e/w/Om/i),
      data/sbdb_centaurs.json, data/sbdb_outer_ss.json.

Outputs: results/step_b17_injection_chain.json,
         figures/supplementary/step_b17_injection_chain.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_052_injection_chain")
tee_stdout(logger)
logger.header("Injected-population depth chain")

from pathlib import Path
import json
import re
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260927)
AXIS = 49.0
ANTI = 229.0
A_NEP = 30.0          # not used; Jupiter semimajor below
A_J = 5.2044

PACKED = re.compile(r"^[CPDIA][IJK](\d\d)([A-Z])")
HALF = {"A": (1, 8), "B": (1, 23), "C": (2, 8), "D": (2, 22),
        "E": (3, 8), "F": (3, 23), "G": (4, 8), "H": (4, 23),
        "J": (5, 8), "K": (5, 23), "L": (6, 8), "M": (6, 23),
        "N": (7, 8), "O": (7, 23), "P": (8, 8), "Q": (8, 23),
        "R": (9, 8), "S": (9, 23), "T": (10, 8), "U": (10, 23),
        "V": (11, 8), "W": (11, 23), "X": (12, 8), "Y": (12, 23)}
CENT = {"I": 1800, "J": 1900, "K": 2000}


def sun_ecl_lon(year, month, day):
    from astropy.time import Time
    from astropy.coordinates import get_sun
    t = Time(f"{year:04d}-{month:02d}-{day:02d}T00:00:00",
             format="isot", scale="utc")
    return float(get_sun(t).geocentrictrueecliptic.lon.deg) % 360


def cmean(th):
    z = np.exp(1j * np.deg2rad(th)).mean()
    return float(np.abs(z)), float(np.rad2deg(np.angle(z)) % 360)


def tisserand(a, e, i):
    return A_J / a + 2 * np.cos(np.deg2rad(i)) * \
        np.sqrt((a / A_J) * (1 - e * e))


def main():
    # ---------- load CometEls ----------
    cq, ce, cw, co, ci, cl_opp = [], [], [], [], [], []
    for line in ((DATA_RAW / "mpc" / "CometEls.txt")).read_text().splitlines():
        if len(line.strip()) < 50:
            continue
        try:
            q = float(line[30:39]); e = float(line[41:49])
            w = float(line[51:59]); om = float(line[61:69])
            inc = float(line[71:79])
        except ValueError:
            continue
        if e >= 1 or e < 0:
            continue
        a = q / (1 - e)
        tj = tisserand(a, e, inc)
        if not (2 < tj < 3):           # JFC population
            continue
        lam = np.nan
        dstr = line[:12].strip()
        m = PACKED.match(dstr)
        if m and m.group(2) in HALF:
            yr = CENT[dstr[1]] + int(m.group(1))
            mo, dy = HALF[m.group(2)]
            lam = (sun_ecl_lon(yr, mo, dy) + 180) % 360
        cq.append(q); ce.append(e); cw.append(w)
        co.append(om); ci.append(inc); cl_opp.append(lam)
    cq = np.array(cq); cw = np.array(cw); co = np.array(co)
    ci = np.array(ci); cl_opp = np.array(cl_opp)
    cvp = (co + cw) % 360
    out = {"axis_deg": AXIS, "anti_deg": ANTI, "n_jfc": len(cq),
           "inputs": ["data/raw/mpc/CometEls.txt",
                      "data/raw/sbdb/sbdb_outer_ss.json",
                      "data/raw/sbdb/sbdb_centaurs.json"]}

    # ---------- I1 the q-progression ----------
    qbins = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
             (5, 6), (6, 15)]
    prog = []
    for lo_, hi_ in qbins:
        m = (cq > lo_) & (cq <= hi_)
        if m.sum() < 10:
            continue
        R, mu = cmean(cvp[m])
        Ro, mo = cmean(co[m])
        Rw, mw = cmean(cw[m])
        prog.append({"q_bin": f"{lo_}-{hi_}", "N": int(m.sum()),
                     "R_varpi": round(R, 3),
                     "mu_varpi": round(mu, 1),
                     "R_om": round(Ro, 3), "mu_om": round(mo, 1),
                     "R_w": round(Rw, 3), "mu_w": round(mw, 1)})
        print(f"I1 q {lo_}-{hi_}: N={m.sum():3d} varpi {mu:5.1f} "
              f"(R={R:.2f})  Om {mo:5.1f}  w {mw:5.1f}")
    out["I1_q_progression"] = {
        "rows": prog,
        "note": ("mean varpi rotates ~150 deg across the JFC "
                 "q range: ~+40 deg deep (aligned with the "
                 "detached-TNO axis) to ~190 deg shallow "
                 "(approaching the anti-axis sector; the dated "
                 "JFc subset of step_b12 sits at ~220 deg).  "
                 "Neither Om nor w alone "
                 "carries the progression")}

    # ---------- I2 footprints per q-bin ----------
    foot = []
    for lo_, hi_ in [(0, 3), (3, 5), (5, 15)]:
        m = (cq > lo_) & (cq <= hi_) & np.isfinite(cl_opp)
        if m.sum() < 15:
            continue
        Rl, ml = cmean(cl_opp[m])
        Rv, mv = cmean(cvp[m])
        foot.append({"q_bin": f"{lo_}-{hi_}", "N": int(m.sum()),
                     "R_opp": round(Rl, 3),
                     "mu_opp": round(ml, 1),
                     "R_varpi": round(Rv, 3),
                     "mu_varpi": round(mv, 1)})
        print(f"I2 q {lo_}-{hi_}: N={m.sum():3d} lam_opp "
              f"{ml:5.1f} (R={Rl:.2f})  varpi {mv:5.1f} "
              f"(R={Rv:.2f})")
    out["I2_footprints"] = {
        "rows": foot,
        "note": ("the discovery footprints are weak and their "
                 "means do not track the varpi progression -- "
                 "pointing does not reproduce the 180 deg "
                 "rotation")}

    # ---------- I3 pairing test on deep JFCs ----------
    deep = (cq < 3) & np.isfinite(cl_opp)
    nd = int(deep.sum())
    vp_d = cvp[deep]; lo_d = cl_opp[deep]
    R_obs, mu_obs = cmean(vp_d)
    Delta = (vp_d - lo_d) % 360
    n_mc = 20000
    null = np.empty(n_mc)
    for k in range(0, n_mc, 2000):
        nb = min(2000, n_mc - k)
        sh = rng.permuted(np.tile(Delta, (nb, 1)), axis=1)
        with np.errstate(all="ignore"):
            z = np.exp(1j * np.deg2rad(
                (lo_d[None, :] + sh) % 360)).mean(axis=1)
        null[k:k + nb] = np.abs(z)
    p_pair = float((np.sum(null >= R_obs) + 1) / (n_mc + 1))
    out["I3_pairing_test_q<3"] = {
        "N": nd, "R_obs": round(R_obs, 3),
        "mu_obs": round(mu_obs, 1),
        "null_mean_R": round(float(null.mean()), 3),
        "p": float(p_pair),
        "note": ("deep injected JFCs: does the observed "
                 "(lam_opp, Delta) pairing beat the sample's "
                 "own marginals?")}
    print(f"I3 q<3: N={nd} R={R_obs:.3f}@{mu_obs:.0f} "
          f"null {null.mean():.3f} p={p_pair:.4f}")

    # ---------- I4 the chain ----------
    d2 = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    vp_det = []
    for r in d2["data"]:
        o = dict(zip(d2["fields"], r))
        try:
            a_, q_, cc = float(o["a"]), float(o["q"]), \
                int(o["condition_code"] or 9)
            om_, w_ = float(o["om"]), float(o["w"])
        except (TypeError, ValueError):
            continue
        if a_ > 150 and q_ > 30 and cc <= 3:
            vp_det.append((om_ + w_) % 360)
    dc = json.loads(((DATA_RAW / "sbdb" / "sbdb_centaurs.json")).read_text())
    vp_cen = []
    for r in dc["data"]:
        o = dict(zip(dc["fields"], r))
        try:
            a_ = float(o["a"]); om_ = float(o["om"])
            w_ = float(o["w"])
        except (TypeError, ValueError):
            continue
        if 5 < a_ < 30:
            vp_cen.append((om_ + w_) % 360)
    R_det, mu_det = cmean(np.array(vp_det))
    R_cen, mu_cen = cmean(np.array(vp_cen))
    R_jd, mu_jd = cmean(cvp[cq < 3])
    R_js, mu_js = cmean(cvp[cq > 5])
    out["I4_chain"] = {
        "detached_TNO": {"N": len(vp_det), "R": round(R_det, 3),
                         "mu": round(mu_det, 1)},
        "centaur": {"N": len(vp_cen), "R": round(R_cen, 3),
                    "mu": round(mu_cen, 1)},
        "JFC_q<3": {"N": int((cq < 3).sum()), "R": round(R_jd, 3),
                    "mu": round(mu_jd, 1)},
        "JFC_q>5": {"N": int((cq > 5).sum()), "R": round(R_js, 3),
                    "mu": round(mu_js, 1)},
        "note": ("the injection chain endpoints land on both "
                 "ends of the measured axis: detached TNOs at "
                 "+49, deep injected JFCs at ~+40, shallow "
                 "transitional JFCs at ~+190-225; the transit "
                 "population (centaurs) is isotropic -- "
                 "scrambled between the two organized ends")}
    print(f"I4: det {mu_det:.0f} (N={len(vp_det)}) | cen "
          f"{mu_cen:.0f} R={R_cen:.2f} (N={len(vp_cen)}) | "
          f"JFC<3 {mu_jd:.0f} (N={(cq<3).sum()}) | JFC>5 "
          f"{mu_js:.0f} (N={(cq>5).sum()})")

    RES.mkdir(exist_ok=True)
    (RES / "step_b17_injection_chain.json").write_text(
        json.dumps(out, indent=1))

    # ---------- figure ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))

    qmid = [(p["q_bin"], p["mu_varpi"], p["R_varpi"], p["N"])
            for p in prog]
    xs = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8][:len(qmid)]
    ax[0].errorbar(xs, [p[1] for p in qmid], fmt="o-",
                   color="crimson")
    ax[0].axhline(AXIS, color="steelblue", ls="--", lw=1,
                  label="axis 49")
    ax[0].axhline(ANTI, color="purple", ls=":", lw=1,
                  label="anti 229")
    ax[0].set(xlabel="q bin centre (AU)", ylabel="mean varpi (deg)",
              ylim=(0, 360), yticks=range(0, 361, 60),
              title="I1: JFC varpi direction\nvs perihelion depth")
    ax[0].legend(fontsize=7)

    for f in foot:
        ax[1].errorbar(float(f["q_bin"].split("-")[0]) + 1,
                       f["mu_varpi"], fmt="o", color="crimson",
                       ms=9, label=None)
        ax[1].errorbar(float(f["q_bin"].split("-")[0]) + 1,
                       f["mu_opp"], fmt="s", color="steelblue",
                       ms=9)
    ax[1].axhline(AXIS, color="steelblue", ls="--", lw=1)
    ax[1].axhline(ANTI, color="purple", ls=":", lw=1)
    ax[1].set(xlabel="q bin (AU)", ylabel="mean direction (deg)",
              ylim=(0, 360), yticks=range(0, 361, 60),
              title="I2: varpi (o) vs footprint (s)")
    ax[1].legend(["varpi", "lam_opp"], fontsize=7)

    ax[2].hist(null, bins=60, color="steelblue", alpha=0.7)
    ax[2].axvline(R_obs, color="r", lw=1.5)
    ax[2].set(xlabel="R(varpi)", ylabel="count",
              title=f"I3: pairing test q<3\n"
                    f"obs={R_obs:.2f} null~{null.mean():.2f} "
                    f"p={p_pair:.3f}")

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b17_injection_chain.png", dpi=300)
    logger.data_save(RESULTS / "step_b17_injection_chain.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b17_injection_chain.png")


if __name__ == "__main__":
    main()
