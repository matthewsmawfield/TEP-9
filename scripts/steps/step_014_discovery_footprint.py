#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b5: the sample's own discovery footprint
=================================================================

The last open loophole in the bias discussion (step b1) is that
the detached objects may have been discovered by surveys whose
combined footprint was narrower than the generic control
footprint -- the "targeted pointing" channel.  This step
reconstructs the *actual* discovery footprint of the 44 detached
objects from real metadata rather than assuming one.

Every provisional designation encodes its own discovery epoch:
(YYYY Lx) -> the letter gives the half-month of discovery
(A = Jan 1-15, B = Jan 16-31, ... Y = Dec 16-31; I skipped).
An object found in a given half-month was, with high
probability, found near opposition: its discovery ecliptic
longitude lies within ~+-60-90 deg of the opposition longitude
lam_opp = lam_sun + 180 on that date.  The distribution of
lam_opp across the sample IS the true discovery footprint --
measured, not modeled.

F1  The footprint itself: R(lam_opp) and its distribution.
    If the extreme-TNO discoveries came from a narrow seasonal
    window, R(lam_opp) is large and the pointing channel can in
    principle reach the observed R(varpi) = 0.33.

F2  Realized coupling: varpi - lam_opp per object vs the OSSOS
    benchmark (mean |dphi| ~ 51 deg).

F3  The exact ceiling test: under the strongest possible
    pointing model every object's perihelion direction equals
    its discovery longitude, so R(varpi)_max = R(lam_opp) --
    a hard bound computed on the real footprint, not a proxy.
    Then the realistic version: varpi_i ~ VM(lam_opp_i, kappa)
    for a fraction f of objects, uniform otherwise, with
    kappa calibrated to the OSSOS |dphi| distribution --
    does the TRUE footprint + measured coupling reproduce
    the cluster?

F4  Comparison against the generic control footprint
    (control-pool Omega marginal): is the detached discovery
    footprint actually narrower than the generic one?

Caveats: lam_opp is a proxy for discovery longitude (objects
can be found up to ~90 deg off opposition, smearing the true
footprint further -- which only widens it, making the ceiling
argument conservative in the direction that matters).
first_observation precoveries are not discovery events; the
designation date is the discovery epoch.

Outputs: results/step_b5_discovery_footprint.json,
         figures/step_b5_discovery_footprint.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_014_discovery_footprint")
tee_stdout(logger)
logger.header("Per-object discovery footprint")

from pathlib import Path
import json
import re
import numpy as np
from scipy.stats import vonmises

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260919)
N_MC = 20000

# half-month letter -> (month, mid-day)
HALF = {"A": (1, 8), "B": (1, 23), "C": (2, 8), "D": (2, 22),
        "E": (3, 8), "F": (3, 23), "G": (4, 8), "H": (4, 23),
        "J": (5, 8), "K": (5, 23), "L": (6, 8), "M": (6, 23),
        "N": (7, 8), "O": (7, 23), "P": (8, 8), "Q": (8, 23),
        "R": (9, 8), "S": (9, 23), "T": (10, 8), "U": (10, 23),
        "V": (11, 8), "W": (11, 23), "X": (12, 8), "Y": (12, 23)}

DESIG = re.compile(r"\((\d{4})\s*([A-Z])([A-Z]?\d*)\)")


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
        except (TypeError, ValueError):
            continue
        o["name"] = str(o["full_name"]).strip()
        objs.append(o)
    return objs


def sun_ecl_lon(year, month, day):
    """Apparent geocentric ecliptic longitude of the Sun (deg)."""
    from astropy.time import Time
    from astropy.coordinates import get_sun
    t = Time(f"{year:04d}-{month:02d}-{day:02d}T00:00:00",
             format="isot", scale="utc")
    return float(get_sun(t).geocentrictrueecliptic.lon.deg) % 360


def circ_R(a):
    return abs(np.exp(1j * np.asarray(a)).mean())


def main():
    objs = load_sbdb()
    det = [o for o in objs if o["a"] > 150 and o["q"] > 30
           and o["cc"] <= 3]
    ctl = [o for o in objs if 30 < o["a"] <= 150 and o["q"] > 30
           and o["cc"] <= 3]

    rows = []
    fails = []
    for o in det:
        m = DESIG.search(o["name"])
        if not m:
            fails.append(o["name"])
            continue
        yr, half = int(m.group(1)), m.group(2)
        if half not in HALF:
            fails.append(o["name"])
            continue
        mo, dy = HALF[half]
        lam_sun = sun_ecl_lon(yr, mo, dy)
        lam_opp = (lam_sun + 180.0) % 360
        varpi = (o["om"] + o["w"]) % 360
        rows.append({"name": o["name"],
                     "desig": f"{yr} {half}{m.group(3)}",
                     "disc_date": f"{yr}-{mo:02d}-{dy:02d}",
                     "lam_opp": round(lam_opp, 1),
                     "varpi": round(varpi, 1),
                     "om": round(o["om"], 1),
                     "w": round(o["w"], 1)})
    print(f"parsed discovery epochs for {len(rows)}/{len(det)} "
          f"objects; unparseable: {fails}")

    vp = np.deg2rad([r["varpi"] for r in rows])
    lo = np.deg2rad([r["lam_opp"] for r in rows])
    om_c = np.deg2rad([o["om"] for o in ctl])

    out = {"sample": "a>150, q>30, cc<=3",
           "n_parsed": len(rows), "unparseable": fails,
           "objects": rows}

    # ---------------- F1 the true footprint ----------------
    R_lo = circ_R(lo)
    lam_mean = float(np.rad2deg(
        np.angle(np.exp(1j * lo).mean())) % 360)
    R_vp = circ_R(vp)
    out["F1_footprint"] = {
        "R_lam_opp": round(float(R_lo), 3),
        "mean_lam_opp_deg": round(lam_mean, 1),
        "R_varpi": round(float(R_vp), 3),
        "mean_varpi_deg": round(float(np.rad2deg(
            np.angle(np.exp(1j * vp).mean())) % 360), 1),
        "note": "R(lam_opp) is the concentration of the sample's "
                "ACTUAL discovery-opposition footprint"}
    print(f"F1 R(lam_opp)={R_lo:.3f} mean={lam_mean:.0f} | "
          f"R(varpi)={R_vp:.3f}")

    # ---------------- F2 realized coupling ----------------
    dphi = (vp - lo + np.pi) % (2 * np.pi) - np.pi
    out["F2_coupling"] = {
        "mean_abs_dphi_deg": round(
            float(np.rad2deg(np.abs(dphi).mean())), 1),
        "frac_within_60": round(
            float((np.abs(dphi) < np.pi / 3).mean()), 3),
        "frac_within_90": round(
            float((np.abs(dphi) < np.pi / 2).mean()), 3),
        "ossos_benchmark_abs_dphi_deg": 50.8,
        "note": "compare to OSSOS measured mean |varpi-lam_d| "
                "= 50.8 deg"}
    print(f"F2 mean|varpi-lam_opp|="
          f"{np.rad2deg(np.abs(dphi).mean()):.0f} deg, "
          f"frac<60: {(np.abs(dphi) < np.pi/3).mean():.2f}")

    # ---------------- F3 ceiling + realistic footprint model ----------------
    out["F3_ceiling"] = {
        "max_pointing_R_varpi": round(float(R_lo), 3),
        "observed_R_varpi": round(float(R_vp), 3),
        "note": "even if varpi_i = lam_opp_i exactly (perfect "
                "perihelion-proximity + opposition discovery), "
                "the realized footprint gives R = R(lam_opp)"}

    # realistic: fraction f at VM(lam_opp_i, kappa), rest uniform
    fs = [0.0, 0.25, 0.5, 0.75, 1.0]
    kap = 2.4   # OSSOS-calibrated |dphi| spread
    sweep = []
    for f in fs:
        nb = int(round(f * len(rows)))
        draw = np.empty((N_MC, len(rows)))
        # biased part: VM about each object's own lam_opp
        if nb:
            draw[:, :nb] = (np.broadcast_to(lo[:nb], (N_MC, nb))
                            + vonmises.rvs(kap, size=(N_MC, nb),
                                           random_state=rng)) \
                % (2 * np.pi)
        draw[:, nb:] = rng.uniform(0, 2 * np.pi,
                                   (N_MC, len(rows) - nb))
        R_m = np.abs(np.exp(1j * draw).mean(axis=1))
        sweep.append({"f_coupled": f,
                      "R_mean": round(float(R_m.mean()), 3),
                      "R_p95": round(float(np.percentile(R_m, 95)),
                                     3),
                      "p": float((int((R_m >= R_vp).sum()) + 1)
                                 / (N_MC + 1))})
        print(f"F3 f={f:.2f}: <R>={R_m.mean():.3f} "
              f"p95={np.percentile(R_m, 95):.3f} "
              f"p={sweep[-1]['p']:.4f}")
    out["F3_realistic_sweep"] = {
        "kappa": kap, "rows": sweep,
        "note": "varpi_i ~ VM(lam_opp_i, kappa) for fraction f; "
                "uniform otherwise.  lam_opp_i is each object's "
                "OWN measured discovery longitude."}

    # ---------------- F4 vs generic footprint ----------------
    out["F4_vs_generic"] = {
        "R_control_Omega": round(float(circ_R(om_c)), 3),
        "R_detached_lam_opp": round(float(R_lo), 3),
        "note": "is the detached discovery footprint narrower "
                "than the generic control footprint?"}
    print(f"F4 R(control Om)={circ_R(om_c):.3f} vs "
          f"R(detached lam_opp)={R_lo:.3f}")

    # ---------------- F5 direction mismatch ----------------
    # Under maximal pointing the model reproduces the footprint's
    # own mean direction (~11 deg), not the observed 49 deg.
    # Quantify: in the f=1 realizations, how often does the mean
    # direction land near the observed axis?
    mu_obs = np.deg2rad(49.1)
    nb = len(rows)
    draw1 = (np.broadcast_to(lo, (N_MC, nb))
             + vonmises.rvs(kap, size=(N_MC, nb),
                            random_state=rng)) % (2 * np.pi)
    mu_m = np.angle(np.exp(1j * draw1).mean(axis=1))
    dmu = (mu_m - mu_obs + np.pi) % (2 * np.pi) - np.pi
    out["F5_direction_mismatch"] = {
        "footprint_mean_deg": round(lam_mean, 1),
        "observed_mean_deg": 49.1,
        "offset_deg": round(float(np.rad2deg(
            (mu_obs - np.deg2rad(lam_mean) + np.pi)
            % (2 * np.pi) - np.pi)), 1),
        "frac_within30_of_obs": round(
            float((np.abs(dmu) < np.pi / 6).mean()), 4),
        "note": "pure pointing predicts mean(varpi) ~ mean("
                "lam_opp) = 11 deg; observed is 49 deg -- a ~38 "
                "deg systematic offset the footprint does not "
                "explain"}
    print(f"F5 footprint mean={lam_mean:.0f} vs observed 49 "
          f"(offset {out['F5_direction_mismatch']['offset_deg']}); "
          f"f=1 realizations within 30 deg of 49: "
          f"{out['F5_direction_mismatch']['frac_within30_of_obs']}")

    # ---------------- F6 clustered objects' discoveries ----------------
    mu = np.deg2rad(49.1)
    cl = np.abs((vp - mu + np.pi) % (2 * np.pi) - np.pi) < np.pi / 3
    out["F6_clustered_discoveries"] = {
        "n_clustered": int(cl.sum()),
        "their_lam_opp_deg": sorted(
            round(float(np.rad2deg(x)), 1) for x in lo[cl]),
        "frac_clustered_found_near_axis": round(float(
            (np.abs((lo[cl] - mu + np.pi) % (2 * np.pi) - np.pi)
             < np.pi / 3).mean()), 3) if cl.sum() else None,
        "note": "if the clustered objects were all discovered "
                "near lam~49, pointing could account for them; "
                "a broad discovery spread cannot"}
    print(f"F6 {cl.sum()} clustered objects; their lam_opp: "
          f"{out['F6_clustered_discoveries']['their_lam_opp_deg']}")

    RES.mkdir(exist_ok=True)
    (RES / "step_b5_discovery_footprint.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.4))

    ax[0].hist(np.rad2deg(lo), bins=np.arange(0, 361, 20),
               color="steelblue", edgecolor="k", lw=0.4,
               label="lam_opp at discovery")
    ax[0].hist(np.rad2deg(vp), bins=np.arange(0, 361, 20),
               histtype="step", color="crimson", lw=1.6,
               label="varpi")
    ax[0].axvline(49.1, color="r", ls="--", lw=1)
    ax[0].set(xlabel="ecliptic longitude (deg)", ylabel="N",
              title="F1: discovery footprint vs\nperihelion "
                    "directions")
    ax[0].legend(fontsize=7)

    ax[1].scatter(np.rad2deg(lo), np.rad2deg(vp), s=36,
                  c="steelblue", edgecolor="k", lw=0.4)
    ax[1].plot([0, 360], [0, 360], "r:", lw=1,
               label="varpi = lam_opp")
    ax[1].axhline(49.1, color="k", ls="--", lw=0.8,
                  label="cluster axis")
    ax[1].set(xlabel="lam_opp at discovery (deg)",
              ylabel="varpi (deg)",
              title="F2: realized coupling\n(dash = OSSOS "
                    "benchmark)")
    ax[1].legend(fontsize=7)

    fr = [r["f_coupled"] for r in sweep]
    ax[2].plot(fr, [r["R_mean"] for r in sweep], "o-",
               color="navy", label="null mean")
    ax[2].plot(fr, [r["R_p95"] for r in sweep], "o--",
               color="navy", alpha=0.5, label="null 95th")
    ax[2].axhline(R_vp, color="r", lw=1.6,
                  label=f"observed R(varpi)={R_vp:.2f}")
    ax[2].axhline(R_lo, color="purple", ls=":", lw=1.2,
                  label=f"R(lam_opp)={R_lo:.2f} = hard ceiling")
    ax[2].set(xlabel="fraction coupled to own lam_opp",
              ylabel="R(varpi)",
              title="F3: can the TRUE footprint\nfake the "
                    "cluster?")
    ax[2].legend(fontsize=7)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "step_b5_discovery_footprint.png", dpi=150)
    print("wrote results/step_b5_discovery_footprint.json, "
          "figures/step_b5_discovery_footprint.png")


if __name__ == "__main__":
    main()
