#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b14: environment checks
==============================================

Two external-structure checks that the bias machinery has not
yet covered:

E1  The galactic tide.  The known external field acting on
    outer solar-system orbits is the galactic tidal field --
    quadrupolar, with preferred directions set by the galactic
    plane, not a pointed direction.  If the detached cluster
    were tidal, perihelia should organize relative to the
    galactic equator/poles and the signature should be
    axis-symmetric (m=2).  Test: galactic latitude of the
    perihelion vectors vs isotropic and conditioned nulls;
    plus the axis's own galactic position.

E2  The omega decomposition.  varpi = Om + w.  If the varpi
    cluster were intrinsic and Om merely footprint-biased,
    mean w should land near (mean varpi - mean Om) ~ 269 deg.
    The observed mean w tells whether the w-clustering
    (Trujillo-Sheppard element) is a trivial consequence or
    carries independent structure.

Outputs: results/step_b14_environment_checks.json,
         figures/step_b14_environment_checks.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_021_environment_checks")
tee_stdout(logger)
logger.header("Galactic-tide and omega decomposition")

from pathlib import Path
import json
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260926)
N_MC = 20000


def load():
    d = json.loads(((DATA_RAW / "sbdb" / "sbdb_outer_ss.json")).read_text())
    det = []
    for r in d["data"]:
        o = dict(zip(d["fields"], r))
        try:
            a, q, cc = float(o["a"]), float(o["q"]), \
                int(o["condition_code"] or 9)
            om, w, i = float(o["om"]), float(o["w"]), float(o["i"])
        except (TypeError, ValueError):
            continue
        if a > 150 and q > 30 and cc <= 3:
            det.append({"om": om, "w": w, "i": i})
    return det


def main():
    det = load()
    n = len(det)
    om = np.deg2rad([o["om"] for o in det])
    w = np.deg2rad([o["w"] for o in det])
    i = np.deg2rad([o["i"] for o in det])
    out = {"sample": "a>150, q>30, cc<=3", "N": n}

    # ---------------- E1 galactic tide ----------------
    pv = np.stack([np.cos(w) * np.cos(om)
                   - np.sin(w) * np.cos(i) * np.sin(om),
                   np.cos(w) * np.sin(om)
                   + np.sin(w) * np.cos(i) * np.cos(om),
                   np.sin(w) * np.sin(i)], axis=1)
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    lam = np.arctan2(pv[:, 1], pv[:, 0])
    bet = np.arcsin(np.clip(pv[:, 2], -1, 1))
    sc = SkyCoord(lon=np.rad2deg(lam) * u.deg,
                  lat=np.rad2deg(bet) * u.deg,
                  frame="geocentrictrueecliptic")
    gb = np.deg2rad(sc.galactic.b.deg)
    sinb = np.abs(np.sin(gb))

    # isotropic + conditioned null for mean |sin b|
    u_mc = rng.uniform(-1, 1, (N_MC, n))
    ph = rng.uniform(0, 2 * np.pi, (N_MC, n))
    rxy = np.sqrt(1 - u_mc ** 2)
    nul = np.stack([rxy * np.cos(ph), rxy * np.sin(ph), u_mc],
                   axis=2)
    # galactic pole in ecliptic coords: l=0,b=90 -> convert
    pole_gal = SkyCoord(l=0 * u.deg, b=90 * u.deg, frame="galactic")
    pole_ecl = pole_gal.geocentrictrueecliptic
    pl = np.deg2rad(float(pole_ecl.lon.deg))
    pb = np.deg2rad(float(pole_ecl.lat.deg))
    gp = np.array([np.cos(pb) * np.cos(pl),
                   np.cos(pb) * np.sin(pl), np.sin(pb)])
    sinb_mc = np.abs(nul @ gp).mean(axis=1)
    p_e1 = float((int((sinb_mc >= sinb.mean()).sum()) + 1)
                 / (N_MC + 1))

    out["E1_galactic_tide"] = {
        "mean_abs_sin_b_gal": round(float(sinb.mean()), 3),
        "isotropic_expectation": 0.5,
        "p_isotropic": float(p_e1),
        "median_abs_b_deg": round(
            float(np.rad2deg(np.median(np.abs(gb)))), 1),
        "axis_galactic": {"l": 182.0, "b": -44.0},
        "m1_vs_m2": {"m1_p": 0.0074, "m2_p": 0.081,
                     "source": "step_02/03 harmonic decomposition"},
        "note": ("perihelia sit at mid galactic latitudes "
                 "(mean|sin b| 0.61, p=0.006 vs isotropic) -- "
                 "but this is a restatement of the cluster "
                 "sitting at b=-44, not independent tidal "
                 "evidence.  The discriminators: no "
                 "concentration toward the galactic equator "
                 "band the tidal quadrupole would produce, "
                 "the axis is off both plane and pole, and "
                 "the signature is m=1 where a tide is m=2.  "
                 "The galactic tide is disfavoured as the "
                 "organizing field")}
    print(f"E1 mean|sin b|={sinb.mean():.3f} p={p_e1:.4f}")

    # ---------------- E2 omega decomposition ----------------
    m_om = np.rad2deg(np.angle(np.exp(1j * om).mean())) % 360
    m_w = np.rad2deg(np.angle(np.exp(1j * w).mean())) % 360
    m_vp = np.rad2deg(np.angle(np.exp(1j * (om + w)).mean())) % 360
    pred_w = (m_vp - m_om) % 360
    d_w = float((m_w - pred_w + 180) % 360 - 180)

    # distribution of w - (vp - om) per-object is trivially 0;
    # the informative quantity is whether the CIRCULAR means
    # decompose: compute the residual direction consistency
    resid = (w - (m_vp - m_om) * np.pi / 180 + np.pi) % (2 * np.pi) - np.pi
    out["E2_omega_decomposition"] = {
        "mean_Om_deg": round(float(m_om), 1),
        "mean_w_deg": round(float(m_w), 1),
        "mean_varpi_deg": round(float(m_vp), 1),
        "predicted_w_deg": round(float(pred_w), 1),
        "offset_deg": round(d_w, 1),
        "note": ("if the varpi cluster were intrinsic and Om "
                 "purely footprint, mean w would sit near "
                 "mean(varpi)-mean(Om) = 269 deg; observed is "
                 "301 deg -- a ~33 deg discrepancy, so the w "
                 "clustering carries structure not reducible "
                 "to footprint-Om plus clustered-varpi "
                 "(consistent with the b1 finding that R(w) "
                 "survives the joint empirical null)")}
    print(f"E2 mean w={m_w:.0f} vs predicted {pred_w:.0f} "
          f"(offset {d_w:+.0f})")

    RES.mkdir(exist_ok=True)
    (RES / "step_b14_environment_checks.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(9, 4.2))

    ax[0].hist(np.rad2deg(np.abs(gb)), bins=18, color="steelblue",
               edgecolor="k", lw=0.3, density=True)
    ax[0].axvline(30, color="0.5", ls=":", label="isotropic median")
    ax[0].axvline(np.rad2deg(np.median(np.abs(gb))), color="r",
                  ls="--", label=f"observed {np.rad2deg(np.median(np.abs(gb))):.0f}")
    ax[0].set(xlabel="|galactic latitude| of perihelion (deg)",
              ylabel="density",
              title="E1: tide check -- no\nequator concentration")
    ax[0].legend(fontsize=7)

    th = np.deg2rad([0, m_om, m_w, m_vp, pred_w])
    labels = ["", f"Om {m_om:.0f}", f"w {m_w:.0f}",
              f"varpi {m_vp:.0f}", f"w pred {pred_w:.0f}"]
    cols = ["k", "steelblue", "crimson", "purple", "gray"]
    for t, l, c in zip(th[1:], labels[1:], cols[1:]):
        ax[1].arrow(0, 0, np.cos(t), np.sin(t),
                    head_width=0.05, color=c, lw=1.6)
        ax[1].text(1.15 * np.cos(t), 1.15 * np.sin(t), l,
                   fontsize=8, ha="center", va="center")
    ax[1].set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5),
              aspect="equal",
              title="E2: mean directions -- w is not\n"
                    "varpi - Om")
    ax[1].axis("off")

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "step_b14_environment_checks.png", dpi=150)
    print("wrote results/step_b14_environment_checks.json, "
          "figures/step_b14_environment_checks.png")


if __name__ == "__main__":
    main()
