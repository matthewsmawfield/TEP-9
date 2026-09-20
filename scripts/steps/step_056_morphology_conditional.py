#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b21: morphology vs footprint
======================================================

step_b20 found the detached cluster is a plateau, not a
peak: flat inside the 60 deg cap, depleted ring at 60-90,
uniform far field -- top-hat beats von Mises on AIC.

But is the morphology itself footprint-explicable?  The
ring region (varpi ~ 109-139 and ~319-349 deg) sits far
from the detached discovery-footprint centre (~10 deg), so
pointing could plausibly produce BOTH the flat interior
(smeared pile-up near the footprint) AND the ring deficit
(objects rarely discovered with perihelia there).

The test: generate conditional-null realizations of the
detached varpi sample -- each object draws (Om,w) from
controls discovered within +-30 deg of ITS OWN opposition
longitude (the b8 machinery, which inherits the measured
footprint and footprint->elements map at full strength) --
and ask what morphology the null produces.

K1  Expected |varpi - axis| profile under the conditional
    null vs observed.
K2  Null distribution of (AIC_tophat - AIC_vm): does the
    footprint also prefer the top-hat?  If the null prefers
    von Mises/uniform while the data prefer the top-hat,
    the morphology is an independent datum; if the null
    produces top-hats too, it is not.
K3  How often does the null produce a 60-90 deg ring as
    depleted as observed (<= 4 objects)?

Outputs: results/step_b21_morphology_conditional.json,
         figures/supplementary/step_b21_morphology_conditional.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_056_morphology_conditional")
tee_stdout(logger)
logger.header("Morphology vs footprint null")

from pathlib import Path
import json
import re
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260927)
AXIS = 49.0
WIN = 30.0
N_MC = 4000

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


def logL_vm(off, f, sigma):
    pdf = f * np.exp(-0.5 * (off / sigma) ** 2) / (
        sigma * np.sqrt(2 * np.pi)) + (1 - f) / 360.0
    return float(np.log(np.clip(pdf, 1e-12, None)).sum())


def logL_tophat(off, f, W):
    pdf = np.where(np.abs(off) < W,
                   f / (2 * W) + (1 - f) / 360.0,
                   (1 - f) / 360.0)
    return float(np.log(np.clip(pdf, 1e-12, None)).sum())


def fit_delta_aic(off):
    """AIC(tophat) - AIC(vm) on unbinned offsets."""
    best_vm, best_th = -1e9, -1e9
    for f_ in np.linspace(0.05, 0.9, 40):
        for s_ in np.linspace(8, 120, 40):
            best_vm = max(best_vm, logL_vm(off, f_, s_))
        for W_ in np.linspace(20, 150, 40):
            best_th = max(best_th, logL_tophat(off, f_, W_))
    return (2 * 2 - 2 * best_th) - (2 * 2 - 2 * best_vm)


def main():
    d = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    rows = []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            a, q = float(o["a"]), float(o["q"])
            cc = int(o["condition_code"] or 9)
            om, w = float(o["om"]), float(o["w"])
        except (TypeError, ValueError):
            continue
        if not (a > 30 and q > 30 and cc <= 3):
            continue
        m = DESIG.search(str(o["full_name"]))
        if not m or m.group(2) not in HALF:
            continue
        yr, half = int(m.group(1)), m.group(2)
        mo, dy = HALF[half]
        lam = (sun_ecl_lon(yr, mo, dy) + 180) % 360
        rows.append({"a": a, "lam": lam, "om": om, "w": w,
                     "vp": (om + w) % 360})
    a = np.array([r["a"] for r in rows])
    lo = np.array([r["lam"] for r in rows])
    vp = np.array([r["vp"] for r in rows])
    om = np.array([r["om"] for r in rows])
    w = np.array([r["w"] for r in rows])

    det = a > 150
    ctl = ~det
    off_obs = np.abs((vp[det] - AXIS + 180) % 360 - 180)
    N = int(det.sum())
    out = {"axis_deg": AXIS, "n": N, "n_mc": N_MC,
           "window_deg": WIN}

    # ---------- conditional-null realizations ----------
    loc = np.deg2rad(lo[ctl])
    omc = np.deg2rad(om[ctl])
    wc = np.deg2rad(w[ctl])
    lod = np.deg2rad(lo[det])
    pools = []
    for x in lod:
        dd = np.abs((loc - x + np.pi) % (2 * np.pi) - np.pi)
        pools.append(np.where(dd < np.deg2rad(WIN))[0])

    # generate MC realizations of the detached varpi sample
    draw = np.empty((N_MC, N))
    for j, p in enumerate(pools):
        sel = rng.choice(p, size=N_MC)
        draw[:, j] = (omc[sel] + wc[sel]) % (2 * np.pi)
    null_vp = np.rad2deg(draw)
    null_off = np.abs((null_vp - AXIS + 180) % 360 - 180)

    # ---------- K1 profile ----------
    edges = [0, 15, 30, 45, 60, 75, 90, 120, 150, 180]
    obs_frac = [float(((off_obs > l_) & (off_obs <= h_)).mean())
                for l_, h_ in zip(edges[:-1], edges[1:])]
    nul_frac = np.array(
        [[float(((o > l_) & (o <= h_)).mean())
          for o in null_off]
         for l_, h_ in zip(edges[:-1], edges[1:])])
    prof = []
    for i, (l_, h_) in enumerate(zip(edges[:-1], edges[1:])):
        p = float((int((nul_frac[i] <= obs_frac[i]).sum()) + 1)
                  / (nul_frac[i].size + 1)) \
            if obs_frac[i] > nul_frac[i].mean() else \
            float((int((nul_frac[i] >= obs_frac[i]).sum()) + 1)
                  / (nul_frac[i].size + 1))
        prof.append({"bin": f"{l_}-{h_}",
                     "obs": round(obs_frac[i], 3),
                     "null_mean": round(float(nul_frac[i].mean()), 3),
                     "null_p16": round(float(np.percentile(
                         nul_frac[i], 16)), 3),
                     "null_p84": round(float(np.percentile(
                         nul_frac[i], 84)), 3),
                     "tail_p": float(p)})
        print(f"K1 {l_:3d}-{h_:3d}: obs {obs_frac[i]:.3f} vs "
              f"null {nul_frac[i].mean():.3f} "
              f"[{np.percentile(nul_frac[i],16):.3f},"
              f"{np.percentile(nul_frac[i],84):.3f}]")
    out["K1_profile_vs_null"] = prof

    # ---------- K2 morphology preference ----------
    da_obs = fit_delta_aic(off_obs)
    da_null = np.array([fit_delta_aic(o) for o in
                        null_off[:800]])
    out["K2_morphology_null"] = {
        "observed_deltaAIC_th_minus_vm": round(da_obs, 2),
        "null_median": round(float(np.median(da_null)), 2),
        "null_p16_p84": [round(float(np.percentile(da_null, 16)), 2),
                         round(float(np.percentile(da_null, 84)), 2)],
        "frac_null_prefers_tophat": round(
            float((da_null < 0).mean()), 3),
        "note": ("negative delta = top-hat preferred.  If the "
                 "conditional null also prefers top-hat, the "
                 "morphology is footprint-consistent; if the "
                 "null prefers von Mises, the plateau is an "
                 "independent datum")}
    print(f"K2: obs dAIC {da_obs:.1f} vs null med "
          f"{np.median(da_null):.1f} "
          f"[{np.percentile(da_null,16):.1f},"
          f"{np.percentile(da_null,84):.1f}] "
          f"frac tophat {(da_null<0).mean():.2f}")

    # ---------- K3 ring deficit ----------
    ring_obs = int(((off_obs > 60) & (off_obs <= 90)).sum())
    ring_null = ((null_off > 60) & (null_off <= 90)).sum(axis=1)
    p_ring = float((int((ring_null <= ring_obs).sum()) + 1)
                   / (ring_null.size + 1))
    out["K3_ring_deficit"] = {
        "observed_n": ring_obs,
        "null_median": float(np.median(ring_null)),
        "p": float(p_ring),
        "note": ("does the measured footprint produce the "
                 "60-90 deg deficit ring on its own?")}
    print(f"K3: ring obs {ring_obs} vs null med "
          f"{np.median(ring_null):.1f} p={p_ring:.4f}")

    RES.mkdir(exist_ok=True)
    (RES / "step_b21_morphology_conditional.json").write_text(
        json.dumps(out, indent=1))

    # ---------- figure ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))

    mids = [(int(p["bin"].split("-")[0])
             + int(p["bin"].split("-")[1])) / 2 for p in prof]
    wid = [int(p["bin"].split("-")[1])
           - int(p["bin"].split("-")[0]) for p in prof]
    ax[0].bar(mids, [p["obs"] for p in prof], width=wid,
              color="crimson", alpha=0.7, label="observed")
    ax[0].bar(mids, [p["null_mean"] for p in prof], width=wid,
              color="none", edgecolor="steelblue", lw=1.3,
              label="conditional null")
    ax[0].axvline(60, color="k", ls="--", lw=1)
    ax[0].set(xlabel="|varpi - axis| (deg)", ylabel="fraction",
              title="K1: profile vs footprint\nconditional null")
    ax[0].legend(fontsize=8)

    ax[1].hist(da_null, bins=40, color="steelblue", alpha=0.7)
    ax[1].axvline(da_obs, color="r", lw=1.5,
                  label=f"obs {da_obs:.1f}")
    ax[1].axvline(0, color="k", ls=":", lw=1)
    ax[1].set(xlabel="AIC(tophat) - AIC(von Mises)",
              ylabel="count",
              title="K2: morphology under\nthe footprint null")
    ax[1].legend(fontsize=8)

    ax[2].hist(ring_null, bins=np.arange(-0.5, 15, 1),
               color="steelblue", alpha=0.7)
    ax[2].axvline(ring_obs, color="r", lw=1.5,
                  label=f"obs {ring_obs}")
    ax[2].set(xlabel="objects in 60-90 deg ring",
              ylabel="count",
              title="K3: the deficit ring")
    ax[2].legend(fontsize=8)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b21_morphology_conditional.png",
                dpi=300)
    logger.data_save(RESULTS / "step_b21_morphology_conditional.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b21_morphology_conditional.png")


if __name__ == "__main__":
    main()
