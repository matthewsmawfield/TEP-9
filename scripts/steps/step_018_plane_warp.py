#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b7: is the disk warped with radius?
===========================================================

The mean plane of a population is read from its mean orbit pole.
A coherent offset of the mean pole from the ecliptic pole is the
signature of a shared tilted plane; an offset whose direction or
magnitude changes with semimajor axis is a WARP -- the prediction
of the inclination-instability picture (Madigan & McCourt) and
also of a radially-structured field boundary, but NOT of simple
shepherding by a single perturber's plane.

The subtlety: the pole azimuth equals Omega - 90 deg, so the
measured footprint (nodal discovery bias) imprints a coherent
pole offset at EVERY radius even for an unwarped disk.  The
calibration must therefore be against the empirical null --
control-pool (Om, w) resampling (step b1), which carries the
footprint -- not against isotropy.

T1  mean pole per a-bin: offset direction (ecliptic longitude
    of the mean pole) and magnitude (tilt of the mean plane),
    with the empirical-null band -- is the detached mean pole
    offset larger than the footprint alone predicts?

T2  warp profile: offset magnitude and direction in sliding
    a-windows -- a coherent drift of the tilt direction with
    radius is the warp signature.

T3  detached vs classical-belt planes: angular separation of the
    two populations' mean poles, bootstrap significance.

Outputs: results/step_b7_plane_warp.json,
         figures/supplementary/step_b7_plane_warp.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_018_plane_warp")
tee_stdout(logger)
logger.header("Coherent mean-plane warp")

from pathlib import Path
import json
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260919)
N_MC = 5000
N_BOOT = 5000


def load_sbdb(fname="sbdb_outer_ss.json"):
    d = json.loads((DATA / fname).read_text())
    objs = []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            for k in ("a", "e", "i", "om", "w", "q"):
                o[k] = float(o[k])
            o["cc"] = int(o["condition_code"]) \
                if o["condition_code"] is not None else 9
            o["arc"] = float(o["data_arc"]) if o["data_arc"] else 0
        except (TypeError, ValueError):
            continue
        objs.append(o)
    return objs


def pole_vec(om, i):
    n = np.stack([np.sin(i) * np.sin(om),
                  -np.sin(i) * np.cos(om),
                  np.cos(i)], axis=-1)
    n[n[..., 2] < 0] *= -1
    return n


def mean_pole(om, i):
    """Direction of the mean (folded) pole: offset from the
    ecliptic pole (deg), offset azimuth (deg), raw mean vector.
    The mean vector is normalised before the offset is read, so
    the tilt measures the COHERENT component only -- a ring of
    poles from random nodes gives ~0 offset regardless of i."""
    n = pole_vec(np.asarray(om), np.asarray(i))
    m = n.mean(axis=0)
    mn = m / np.linalg.norm(m)
    off = float(np.rad2deg(np.arccos(np.clip(mn[2], -1, 1))))
    azm = float(np.rad2deg(np.arctan2(mn[1], mn[0])) % 360)
    return off, azm, m


def main():
    objs = [o for o in load_sbdb()
            if o["q"] > 30 and o["cc"] <= 3 and o["arc"] > 200
            and o["a"] > 30]
    ctl = [o for o in objs if 30 < o["a"] <= 150]
    om_c = np.deg2rad([o["om"] for o in ctl])
    w_c = np.deg2rad([o["w"] for o in ctl])
    pool = np.stack([om_c, w_c], axis=1)
    print(f"secure catalog N={len(objs)}, control pool {len(ctl)}")

    out = {"sample": "a>30, q>30, cc<=3, arc>200",
           "caveat": "the empirical null assumes the detached "
                     "population shares the generic (30-150 AU) "
                     "discovery footprint; if the true footprint "
                     "narrows with a -- plausible, since distant "
                     "objects are found by fewer, deeper surveys "
                     "-- part of the tilt excess at 60-150 AU "
                     "could be deeper selection rather than a "
                     "physical warp.  The radial GROWTH of the "
                     "coherent tilt (T2) is the datum either "
                     "interpretation must explain."}

    # ---------------- T1 mean pole per a-bin vs empirical null ----------------
    bins = [(30, 60), (60, 150), (150, 1e9)]
    t1 = []
    for lo, hi in bins:
        s = [o for o in objs if lo < o["a"] <= hi]
        om = np.deg2rad([o["om"] for o in s])
        i = np.deg2rad([o["i"] for o in s])
        off, azm, _ = mean_pole(om, i)
        # empirical null: same N, i fixed, (om,w) from control
        n = len(s)
        idx = rng.integers(0, len(ctl), (N_MC, n))
        offs = np.empty(N_MC)
        for k in range(N_MC):
            offs[k], _, _ = mean_pole(pool[idx[k], 0],
                                      np.broadcast_to(i, n))
        t1.append({"a_range": [lo, "inf" if hi > 1e8 else hi],
                   "N": n,
                   "tilt_deg": round(off, 2),
                   "azimuth_deg": round(azm, 1),
                   "null_tilt_median": round(
                       float(np.median(offs)), 2),
                   "null_tilt_p95": round(
                       float(np.percentile(offs, 95)), 2),
                   "p_tilt": float((int((offs >= off).sum()) + 1)
                                  / (offs.size + 1))})
        print(f"T1 a {lo}-{hi if hi<1e8 else 'inf'}: N={n} "
              f"tilt={off:.2f} az={azm:.0f} "
              f"(null {np.median(offs):.2f}, "
              f"p={t1[-1]['p_tilt']:.4f})")
    out["T1_pole_offset_per_bin"] = t1

    # ---------------- T2 sliding-window warp profile ----------------
    wins = []
    c = np.log10(35.0)
    while c < np.log10(2000):
        alo, ahi = 10 ** c, 10 ** (c + 0.4)
        s = [o for o in objs if alo < o["a"] <= ahi]
        if len(s) >= 10:
            om = np.deg2rad([o["om"] for o in s])
            i = np.deg2rad([o["i"] for o in s])
            off, azm, _ = mean_pole(om, i)
            # empirical-null tilt band
            n = len(s)
            idx = rng.integers(0, len(ctl), (1000, n))
            no = []
            for k in range(1000):
                o2, _, _ = mean_pole(pool[idx[k], 0],
                                     np.broadcast_to(i, n))
                no.append(o2)
            wins.append({"a_mid": round(float(np.sqrt(alo * ahi)),
                                        1),
                         "N": n, "tilt_deg": round(off, 2),
                         "azimuth_deg": round(azm, 1),
                         "null_p95_tilt": round(
                             float(np.percentile(no, 95)), 2)})
        c += 0.1
    out["T2_warp_profile"] = wins
    print("T2 window profile:",
          [(w["a_mid"], w["tilt_deg"], w["azimuth_deg"])
           for w in wins])

    # ---------------- T3 detached vs classical plane ----------------
    s_cl = [o for o in objs if 40 < o["a"] <= 47 and o["i"] < 10]
    s_de = [o for o in objs if o["a"] > 150]
    om1 = np.deg2rad([o["om"] for o in s_cl])
    i1 = np.deg2rad([o["i"] for o in s_cl])
    om2 = np.deg2rad([o["om"] for o in s_de])
    i2 = np.deg2rad([o["i"] for o in s_de])
    _, _, m1 = mean_pole(om1, i1)
    _, _, m2 = mean_pole(om2, i2)
    m1n, m2n = m1 / np.linalg.norm(m1), m2 / np.linalg.norm(m2)
    sep = float(np.rad2deg(np.arccos(np.clip(
        m1n @ m2n, -1, 1))))
    # bootstrap the separation
    seps = []
    for _ in range(N_BOOT):
        j1 = rng.integers(0, len(s_cl), len(s_cl))
        j2 = rng.integers(0, len(s_de), len(s_de))
        _, _, b1 = mean_pole(om1[j1], i1[j1])
        _, _, b2 = mean_pole(om2[j2], i2[j2])
        b1 /= np.linalg.norm(b1)
        b2 /= np.linalg.norm(b2)
        seps.append(np.rad2deg(np.arccos(np.clip(b1 @ b2, -1, 1))))
    # null: same samples but om randomized (i fixed) -> is the
    # observed separation beyond footprint scatter?
    sep0 = []
    for _ in range(2000):
        _, _, b1 = mean_pole(rng.uniform(0, 2 * np.pi, len(s_cl)),
                             i1)
        _, _, b2 = mean_pole(rng.uniform(0, 2 * np.pi, len(s_de)),
                             i2)
        b1 /= np.linalg.norm(b1)
        b2 /= np.linalg.norm(b2)
        sep0.append(np.rad2deg(np.arccos(np.clip(b1 @ b2, -1, 1))))
    out["T3_detached_vs_classical_plane"] = {
        "n_classical": len(s_cl), "n_detached": len(s_de),
        "separation_deg": round(sep, 2),
        "boot_ci68_deg": [round(float(np.percentile(seps, 16)), 2),
                          round(float(np.percentile(seps, 84)), 2)],
        "uniform_om_null_median_deg": round(
            float(np.median(sep0)), 2),
        "p_vs_uniform": float(((np.array(sep0) >= sep).sum() + 1)
                              / (len(sep0) + 1)),
        "note": "angular distance between the mean poles of the "
                "classical belt (40-47 AU, i<10) and the detached "
                "population -- the physical warp between planes"}
    print(f"T3 planes separated by {sep:.2f} deg "
          f"(CI68 {out['T3_detached_vs_classical_plane']['boot_ci68_deg']})")

    RES.mkdir(exist_ok=True)
    (RES / "step_b7_plane_warp.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.4))

    xm = [w["a_mid"] for w in wins]
    ax[0].semilogx(xm, [w["tilt_deg"] for w in wins], "o-",
                   color="crimson", label="observed mean-pole tilt")
    ax[0].semilogx(xm, [w["null_p95_tilt"] for w in wins], "--",
                   color="0.55", label="footprint null 95th")
    ax[0].axvline(150, color="k", ls=":", lw=0.8)
    ax[0].set(xlabel="a window center (AU)",
              ylabel="mean-pole offset from ecliptic pole (deg)",
              title="T2: warp profile\n(coherent plane tilt vs a)")
    ax[0].legend(fontsize=7)

    ax[1].semilogx(xm, [w["azimuth_deg"] for w in wins], "s-",
                   color="navy")
    ax[1].axhline(49.1, color="r", ls="--", lw=1,
                  label="cluster axis 49 deg")
    ax[1].set(xlabel="a window center (AU)",
              ylabel="mean-pole azimuth (deg)",
              title="T2: tilt direction vs a\n(= mean(Om)-90)")
    ax[1].legend(fontsize=7)

    for j, b in enumerate(t1):
        lab = f"a {b['a_range'][0]}-{b['a_range'][1]}"
        ax[2].bar(j, b["tilt_deg"], color="steelblue",
                  edgecolor="k", lw=0.4,
                  label=lab + f" (N={b['N']})")
        ax[2].plot([j - 0.3, j + 0.3],
                   [b["null_tilt_p95"]] * 2, "r--", lw=1.2)
    ax[2].set_xticks(range(len(t1)))
    ax[2].set_xticklabels([f"a {b['a_range'][0]}-" +
                           str(b['a_range'][1]) for b in t1],
                          fontsize=8)
    ax[2].set(ylabel="mean-pole tilt (deg)",
              title="T1: coherent tilt per bin\n(red = footprint "
                    "null 95th)")
    ax[2].legend(fontsize=7)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b7_plane_warp.png", dpi=300)
    logger.data_save(RESULTS / "step_b7_plane_warp.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b7_plane_warp.png")


if __name__ == "__main__":
    main()
