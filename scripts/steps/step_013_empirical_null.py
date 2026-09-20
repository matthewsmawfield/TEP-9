#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b1: empirical, bias-informed null model
================================================================

Step 02 tested the detached-TNO clustering against a conditioned
null in which om and w are drawn uniform on the circle.  That null
is unrealistic: real surveys point at specific ecliptic longitudes,
objects are found near opposition and near their nodes, and
eccentric objects are found preferentially near perihelion (step 06
measured this directly in OSSOS: mean |varpi - lam_d| = 51 deg,
R_coupling = 0.513).

This step replaces the uniform angular null with an *empirical*
one built from the catalog itself, and then asks how much
perihelion-proximity coupling would be required to manufacture the
observed signal.

E1  Control-resample null.  The ~1200-object control population
    (30 < a <= 150 AU, q > 30, cc <= 3, arc > 200 d) carries the
    real angular discovery filter: its (om, w) joint distribution
    encodes nodal bias, opposition geometry, and ecliptic-pointing
    structure.  For each MC realization the detached sample keeps
    its observed (a, e, i) and draws (om, w) pairs -- jointly, to
    preserve their correlation -- from the control pool.  Any
    clustering the filter itself produces is now inside the null.

E2  Perihelion-proximity augmentation sweep.  The control pool
    under-represents one channel: detached objects, being distant
    and faint, are found much closer to perihelion than the
    controls.  An object found near perihelion near the ecliptic
    has w ~ 0 or 180 deg, so varpi ~ om and the perihelion
    direction inherits the footprint longitude.  A fraction f of
    the null objects is therefore drawn from a "found near
    perihelion" component: om from the control footprint,
    w ~ von Mises(0 or 180, kappa).  The sweep in f asks: at what
    coupling strength does the bias reproduce the observed
    R(varpi) = 0.33?  The OSSOS coupling (step 06: effective
    biased fraction ~ 0.4) marks the empirically calibrated point.

E3  Direct OSSOS calibration.  The OSSOS detached objects' own
    (om, w, varpi, lam_d) joint is used to estimate the biased
    fraction and its angular spread, giving an empirical anchor
    for the E2 sweep rather than an assumed one.

Outputs: results/step_b1_empirical_null.json,
         figures/supplementary/step_b1_empirical_null.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_013_empirical_null")
tee_stdout(logger)
logger.header("Empirical discovery-filter null")

from pathlib import Path
import json
import numpy as np
from scipy.stats import vonmises

ROOT = PROJECT_ROOT
RES, FIG, DATA = RESULTS, RESULTS / "figures", DATA_RAW

rng = np.random.default_rng(20260918)
N_MC = 20000


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


def perih_vec(om, w, i):
    return np.stack([
        np.cos(om) * np.cos(w) - np.sin(om) * np.sin(w) * np.cos(i),
        np.sin(om) * np.cos(w) + np.cos(om) * np.sin(w) * np.cos(i),
        np.sin(w) * np.sin(i)], axis=-1)


def pole_vec(om, i):
    n = np.stack([np.sin(i) * np.sin(om),
                  -np.sin(i) * np.cos(om),
                  np.cos(i)], axis=-1)
    n[n[..., 2] < 0] *= -1
    return n


def circ_R(a):
    return abs(np.exp(1j * np.asarray(a)).mean())


def all_stats(om, w, i):
    varpi = (om + w) % (2 * np.pi)
    return (circ_R(varpi), circ_R(w), circ_R(om),
            np.linalg.norm(perih_vec(om, w, i).mean(axis=0)),
            np.linalg.norm(pole_vec(om, i).mean(axis=0)))


def main():
    objs = load_sbdb()
    det = [o for o in objs if o["a"] > 150 and o["q"] > 30
           and o["cc"] <= 3]
    ctl = [o for o in objs if 30 < o["a"] <= 150 and o["q"] > 30
           and o["cc"] <= 3 and o["arc"] > 200]
    n = len(det)
    print(f"detached N={n}, control pool N={len(ctl)}")

    om_d = np.deg2rad([o["om"] for o in det])
    w_d = np.deg2rad([o["w"] for o in det])
    i_d = np.deg2rad([o["i"] for o in det])
    s_obs = all_stats(om_d, w_d, i_d)
    names = ["R_varpi", "R_omega", "R_Omega", "R_phat", "R_nhat"]
    print("observed:", {k: round(v, 3) for k, v in zip(names, s_obs)})

    # control pool's own angular structure (diagnostic)
    om_c = np.deg2rad([o["om"] for o in ctl])
    w_c = np.deg2rad([o["w"] for o in ctl])
    out = {"sample": {"n_detached": n, "n_control": len(ctl)},
           "control_angular_structure": {
               "R_Omega": round(float(circ_R(om_c)), 3),
               "R_omega": round(float(circ_R(w_c)), 3),
               "R_varpi": round(float(circ_R((om_c + w_c)
                                             % (2 * np.pi))), 3),
               "note": "nonzero R_Omega is the nodal/footprint "
                       "bias the empirical null inherits"},
           "observed": dict(zip(names, [round(float(v), 4)
                                        for v in s_obs]))}
    print("control pool:", out["control_angular_structure"])

    i_mc = np.broadcast_to(i_d, (N_MC, n))

    # ---------------- E1 control-resample null ----------------
    pool = np.stack([om_c, w_c], axis=1)
    idx = rng.integers(0, len(ctl), (N_MC, n))
    om_mc = pool[idx, 0]
    w_mc = pool[idx, 1]

    varpi = (om_mc + w_mc) % (2 * np.pi)
    Rv = np.abs(np.exp(1j * varpi).mean(axis=1))
    Rw = np.abs(np.exp(1j * w_mc).mean(axis=1))
    Ro = np.abs(np.exp(1j * om_mc).mean(axis=1))
    Rp = np.linalg.norm(perih_vec(om_mc, w_mc, i_mc).mean(axis=1),
                        axis=1)
    Rn = np.linalg.norm(pole_vec(om_mc, i_mc).mean(axis=1), axis=1)
    nulls = {"R_varpi": Rv, "R_omega": Rw, "R_Omega": Ro,
             "R_phat": Rp, "R_nhat": Rn}

    out["E1_control_resample_null"] = {}
    for k, arr in nulls.items():
        p = float((int((arr >= s_obs[names.index(k)]).sum()) + 1)
                  / (N_MC + 1))
        out["E1_control_resample_null"][k] = {
            "obs": round(float(s_obs[names.index(k)]), 4),
            "null_mean": round(float(arr.mean()), 4),
            "null_p95": round(float(np.percentile(arr, 95)), 4),
            "p": p}
        print(f"E1 {k:9s} obs={s_obs[names.index(k)]:.3f} "
              f"null={arr.mean():.3f} p={p:.4f}")

    # ---------------- E2 perihelion-proximity sweep ----------------
    # Biased component: om ~ footprint (control marginal),
    # w ~ VM(0 or 180, kappa).  f = fraction of sample in that
    # component.  OSSOS anchor: effective biased fraction ~ 0.4,
    # spread ~ |dphi| ~ 51 deg => kappa ~ 1.5-3.
    om_pool = om_c
    fs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    kappas = {"vm3": 3.0, "vm1p5": 1.5}
    sweep = {}
    for kn, kap in kappas.items():
        rows = []
        for f in fs:
            nb = int(round(f * n))
            # biased draws
            om_b = om_pool[rng.integers(0, len(ctl), (N_MC, nb))] \
                if nb else np.zeros((N_MC, 0))
            branch = rng.integers(0, 2, (N_MC, nb)) * np.pi
            w_b = (branch + vonmises.rvs(kap, size=(N_MC, nb),
                                         random_state=rng)) \
                % (2 * np.pi) if nb else np.zeros((N_MC, 0))
            # unbiased draws
            nu = n - nb
            ui = rng.integers(0, len(ctl), (N_MC, nu))
            om_u = pool[ui, 0]
            w_u = pool[ui, 1]
            om_m = np.concatenate([om_b, om_u], axis=1)
            w_m = np.concatenate([w_b, w_u], axis=1)
            vv = (om_m + w_m) % (2 * np.pi)
            Rv_m = np.abs(np.exp(1j * vv).mean(axis=1))
            Rw_m = np.abs(np.exp(1j * w_m).mean(axis=1))
            rows.append({
                "f_biased": f,
                "p_varpi": float((int((Rv_m >= s_obs[0]).sum()) + 1)
                                / (N_MC + 1)),
                "R_varpi_null_mean": round(float(Rv_m.mean()), 3),
                "R_varpi_null_p95": round(
                    float(np.percentile(Rv_m, 95)), 3),
                "p_omega": float((int((Rw_m >= s_obs[1]).sum()) + 1)
                                / (N_MC + 1)),
                "R_omega_null_mean": round(float(Rw_m.mean()), 3)})
            print(f"E2 {kn} f={f:.1f}: <R(varpi)>={Rv_m.mean():.3f} "
                  f"p={rows[-1]['p_varpi']:.4f} | "
                  f"<R(w)>={Rw_m.mean():.3f} "
                  f"p={rows[-1]['p_omega']:.4f}")
        sweep[kn] = rows

    # E2b -- the maximal pointing channel: the biased object's
    # perihelion longitude equals its discovery longitude
    # (varpi = lam_d + eps) regardless of node geometry.  lam_d
    # drawn from the control footprint; this is the strongest
    # clustering pure pointing can produce -- its ceiling is the
    # footprint's own concentration, R(lam_d) ~ 0.2.
    rows = []
    for f in fs:
        nb = int(round(f * n))
        ld = om_pool[rng.integers(0, len(ctl), (N_MC, nb))] \
            if nb else np.zeros((N_MC, 0))
        vp_b = (ld + vonmises.rvs(3.0, size=(N_MC, nb),
                                  random_state=rng)) \
            % (2 * np.pi) if nb else np.zeros((N_MC, 0))
        om_b = om_pool[rng.integers(0, len(ctl), (N_MC, nb))] \
            if nb else np.zeros((N_MC, 0))
        w_b = (vp_b - om_b) % (2 * np.pi)
        nu = n - nb
        ui = rng.integers(0, len(ctl), (N_MC, nu))
        om_m = np.concatenate([om_b, pool[ui, 0]], axis=1)
        w_m = np.concatenate([w_b, pool[ui, 1]], axis=1)
        vv = (om_m + w_m) % (2 * np.pi)
        Rv_m = np.abs(np.exp(1j * vv).mean(axis=1))
        Rw_m = np.abs(np.exp(1j * w_m).mean(axis=1))
        rows.append({
            "f_biased": f,
            "p_varpi": float((int((Rv_m >= s_obs[0]).sum()) + 1)
                            / (N_MC + 1)),
            "R_varpi_null_mean": round(float(Rv_m.mean()), 3),
            "R_varpi_null_p95": round(
                float(np.percentile(Rv_m, 95)), 3),
            "p_omega": float((int((Rw_m >= s_obs[1]).sum()) + 1)
                            / (N_MC + 1)),
            "R_omega_null_mean": round(float(Rw_m.mean()), 3)})
        print(f"E2b direct f={f:.1f}: <R(varpi)>={Rv_m.mean():.3f} "
              f"p={rows[-1]['p_varpi']:.4f}")
    sweep["direct_varpi_equals_lamd"] = rows
    out["E2_perihelion_proximity_sweep"] = sweep
    out["E2_note"] = ("f_biased = fraction of objects drawn from "
                      "the 'found near perihelion near the ecliptic'"
                      " component (om~footprint, w~VM(0/180,kappa));"
                      " OSSOS calibration suggests f~0.4")

    # ---------------- E3 OSSOS calibration of f ----------------
    try:
        from astropy.io.votable import parse
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        t = parse((DATA_RAW / "ossos" / "ossos_t3char.vot")).get_first_table() \
            .to_table()
        cl = np.array([str(x) for x in t["cl"]])
        det_o = cl == "det"
        Om_o = np.deg2rad(np.array(t["Omega"], float)[det_o])
        w_o = np.deg2rad(np.array(t["omega"], float)[det_o])
        ra = np.array(t["RAJ2000"], float)[det_o]
        de = np.array(t["DEJ2000"], float)[det_o]
        lam = np.deg2rad([float(
            SkyCoord(ra=r * u.deg, dec=d * u.deg, frame="icrs")
            .transform_to("barycentrictrueecliptic").lon.deg) % 360
            for r, d in zip(ra, de)])
        vp = (Om_o + w_o) % (2 * np.pi)
        dphi = (vp - lam + np.pi) % (2 * np.pi) - np.pi
        # effective biased fraction: objects with |dphi| < 90
        f_eff = float((np.abs(dphi) < np.pi / 2).mean())
        # spread of the coupled component
        k_eff = float(vonmises.fit(dphi[np.abs(dphi) < np.pi / 2],
                                 fscale=1)[0]) \
            if (np.abs(dphi) < np.pi / 2).sum() > 3 else np.nan
        w_clus = float(np.mean([circ_R(w_o % (2 * np.pi)),
                                circ_R((w_o + np.pi)
                                       % (2 * np.pi))]))
        out["E3_ossos_calibration"] = {
            "n_detached_ossos": int(det_o.sum()),
            "f_coupled_absdphi_lt90": round(f_eff, 3),
            "kappa_dphi": round(k_eff, 2),
            "mean_abs_dphi_deg": round(
                float(np.rad2deg(np.abs(dphi).mean())), 1),
            "R_omega_ossos": round(float(circ_R(w_o)), 3),
            "R_omega_folded_ossos": round(
                float(circ_R(2 * w_o)), 3),
            "note": "f_coupled is the fraction of OSSOS detached "
                    "objects whose discovery longitude lay within "
                    "90 deg of varpi -- the empirical anchor for f"}
        print(f"E3 OSSOS: f_coupled={f_eff:.2f} "
              f"kappa={k_eff:.2f} R(w)={circ_R(w_o):.3f} "
              f"R(2w)={circ_R(2*w_o):.3f}")
    except Exception as e:
        out["E3_ossos_calibration"] = {"error": str(e)}
        print("E3 failed:", e)

    RES.mkdir(exist_ok=True)
    (RES / "step_b1_empirical_null.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))

    ax[0].hist(np.rad2deg(om_c), bins=np.arange(0, 361, 15),
               color="0.7", label=f"Omega, control (N={len(ctl)})")
    ax[0].hist(np.rad2deg(om_d), bins=np.arange(0, 361, 15),
               histtype="step", color="crimson", lw=1.6,
               label=f"Omega, detached (N={n})")
    ax[0].set(xlabel="deg", ylabel="N",
              title="footprint: node distribution\n(control = empirical null)")
    ax[0].legend(fontsize=7)

    for kn, c in [("vm3", "navy"), ("vm1p5", "teal"),
                  ("direct_varpi_equals_lamd", "purple")]:
        fr = [r["f_biased"] for r in sweep[kn]]
        rm = [r["R_varpi_null_mean"] for r in sweep[kn]]
        r95 = [r["R_varpi_null_p95"] for r in sweep[kn]]
        ax[1].plot(fr, rm, "o-", color=c,
                   label=f"null mean, {kn.replace('_',' ')}")
        ax[1].plot(fr, r95, "o--", color=c, alpha=0.4)
    ax[1].axhline(s_obs[0], color="r", lw=1.6,
                  label=f"observed R(varpi)={s_obs[0]:.2f}")
    ax[1].axvline(0.4, color="k", ls=":", lw=1,
                  label="OSSOS coupling ~0.4")
    ax[1].set(xlabel="fraction biased (found near perihelion)",
              ylabel="R(varpi)",
              title="E2: bias needed to fake the cluster")
    ax[1].legend(fontsize=7)

    for kn, c in [("vm3", "navy"), ("vm1p5", "teal"),
                  ("direct_varpi_equals_lamd", "purple")]:
        fr = [r["f_biased"] for r in sweep[kn]]
        rm = [r["R_omega_null_mean"] for r in sweep[kn]]
        ax[2].plot(fr, rm, "o-", color=c,
                   label=f"null, {kn.replace('_',' ')}")
    ax[2].axhline(s_obs[1], color="r", lw=1.6,
                  label=f"observed R(w)={s_obs[1]:.2f}")
    ax[2].axvline(0.4, color="k", ls=":", lw=1)
    ax[2].set(xlabel="fraction biased", ylabel="R(omega)",
              title="E2: same for argument of perihelion")
    ax[2].legend(fontsize=7)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b1_empirical_null.png", dpi=300)
    logger.data_save(RESULTS / "step_b1_empirical_null.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b1_empirical_null.png")


if __name__ == "__main__":
    main()
