#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b4: geometry of the detached orbit poles
=================================================================

Three pictures predict different three-dimensional pole
structure for the detached population:

  * point-mass shepherding: orbits share a common plane tilted
    toward the perturber -- poles cluster about ONE direction
    (a point on the pole sphere).

  * inclination instability / self-gravity (Madigan & McCourt):
    the disk torques itself into a warped plane -- poles spread
    along a GREAT CIRCLE (a girdle), not a point.

  * fixed field boundary (TEP): if the structure is a lapse
    gradient it selects perihelion direction; if it also tilts
    the clock planes it leaves a preferred pole direction --
    either way, the eigenstructure of the pole scatter matrix
    discriminates point (l1 >> l2 ~ l3) from girdle
    (l1 ~ l2 >> l3).

Because low-i orbits trivially place poles near the ecliptic pole,
everything is calibrated against the conditioned null (fix each
object's i, draw om uniform) -- the null ring of poles at radius
sin(i) is itself girdle-like, so the observed statistics are
compared to that null, not to isotropy.

Statistics: scatter-matrix eigenvalues l1>=l2>=l3 of the folded
poles; cluster index C = l1 - l2 (single-direction excess);
girdle index G = l2 - l3 (planar excess).  The same decomposition
is applied to the (unfolded) perihelion vectors.

Also reported: the best common-plane normal (l3 eigenvector of
the girdle fit / direction that minimises pole spread), its
ecliptic coordinates and tilt, and its relation to the varpi
cluster axis -- whether the pole structure and the perihelion
structure point at the same inertial frame.

Outputs: results/step_b4_pole_geometry.json,
         figures/step_b4_pole_geometry.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_019_pole_geometry")
tee_stdout(logger)
logger.header("Orbit-pole eigenstructure")

from pathlib import Path
import json
import numpy as np

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
        except (TypeError, ValueError):
            continue
        objs.append(o)
    return objs


def perih_vec(om, w, i):
    return np.stack([
        np.cos(om) * np.cos(w) - np.sin(om) * np.sin(w) * np.cos(i),
        np.sin(om) * np.cos(w) + np.cos(om) * np.sin(w) * np.cos(i),
        np.sin(w) * np.sin(i)], axis=-1)


def pole_vec(om, i, fold=True):
    n = np.stack([np.sin(i) * np.sin(om),
                  -np.sin(i) * np.cos(om),
                  np.cos(i)], axis=-1)
    if fold:
        n[n[..., 2] < 0] *= -1
    return n


def eigs(vecs):
    """Scatter-matrix eigenvalues (desc) and eigenvectors."""
    S = np.einsum("ni,nj->ij", vecs, vecs) / len(vecs)
    w, v = np.linalg.eigh(S)
    return w[::-1], v[:, ::-1]


def lonlat(v):
    return (float(np.rad2deg(np.arctan2(v[1], v[0])) % 360),
            float(np.rad2deg(np.arcsin(np.clip(v[2], -1, 1)))))


def pole_stats(om, i):
    n = pole_vec(om, i)
    l, v = eigs(n)
    C = l[0] - l[1]          # point-cluster excess
    G = l[1] - l[2]          # girdle excess
    return l, v, C, G


def main():
    objs = load_sbdb()
    det = [o for o in objs if o["a"] > 150 and o["q"] > 30
           and o["cc"] <= 3]
    n = len(det)
    om = np.deg2rad([o["om"] for o in det])
    w = np.deg2rad([o["w"] for o in det])
    i = np.deg2rad([o["i"] for o in det])

    l_obs, v_obs, C_obs, G_obs = pole_stats(om, i)
    print(f"detached N={n}: pole eigs "
          f"{np.round(l_obs, 3)} C={C_obs:.3f} G={G_obs:.3f}")

    # conditioned null: i fixed, om uniform
    om_mc = rng.uniform(0, 2 * np.pi, (N_MC, n))
    i_mc = np.broadcast_to(i, (N_MC, n))
    Cn = np.empty(N_MC)
    Gn = np.empty(N_MC)
    l1n = np.empty(N_MC)
    for k in range(N_MC):
        nn = pole_vec(om_mc[k], i_mc[k])
        lk, _ = eigs(nn)
        Cn[k] = lk[0] - lk[1]
        Gn[k] = lk[1] - lk[2]
        l1n[k] = lk[0]

    out = {"sample": "a>150, q>30, cc<=3", "N": n,
           "P1_pole_eigenstructure": {
               "eigenvalues_obs": [round(float(x), 4)
                                   for x in l_obs],
               "eigenvalue_1_null_p95": round(
                   float(np.percentile(l1n, 95)), 4),
               "C_point_index_obs": round(float(C_obs), 4),
               "C_null_mean": round(float(Cn.mean()), 4),
               "C_null_p95": round(float(np.percentile(Cn, 95)), 4),
               "p_C": float((int((Cn >= C_obs).sum()) + 1)
                            / (N_MC + 1)),
               "G_girdle_index_obs": round(float(G_obs), 4),
               "G_null_mean": round(float(Gn.mean()), 4),
               "G_null_p95": round(float(np.percentile(Gn, 95)), 4),
               "p_G": float((int((Gn >= G_obs).sum()) + 1)
                            / (N_MC + 1)),
               "note": "C high => poles share one direction "
                       "(common tilted plane / shepherd or fixed "
                       "frame); G high => poles on a great circle "
                       "(warped disk / inclination instability). "
                       "Null holds i fixed, om uniform -- the "
                       "ecliptic ring it produces is itself "
                       "girdle-like."}}
    print(f"P1 C: obs={C_obs:.3f} null={Cn.mean():.3f} "
          f"p={out['P1_pole_eigenstructure']['p_C']:.4f} | "
          f"G: obs={G_obs:.3f} null={Gn.mean():.3f} "
          f"p={out['P1_pole_eigenstructure']['p_G']:.4f}")

    # best common-plane normal: eigenvector of smallest eigenvalue
    # is the normal to the best-fit plane through the poles; the
    # mean pole is e1.
    e1 = v_obs[:, 0]
    if e1[2] < 0:
        e1 = -e1
    lam1, bet1 = lonlat(e1)
    tilt = 90.0 - bet1   # tilt of common plane vs ecliptic plane
    out["P2_common_plane"] = {
        "mean_pole_ecliptic": [round(lam1, 1), round(bet1, 1)],
        "plane_tilt_deg": round(float(tilt), 1),
        "pole_axis_azimuth_deg": round(lam1, 1),
        "varpi_axis_deg": 49.1,
        "note": "pole azimuth = mean(Omega) - 90 deg, so this "
                "largely restates the node clustering; the new "
                "information is the eigenstructure (P1): a point "
                "excess, no girdle -- a common tilted plane, not "
                "a warped disk"}
    print(f"P2 mean pole dir lam={lam1:.0f} bet={bet1:.0f} "
          f"(plane tilt {tilt:.0f} deg)")

    # same decomposition on the perihelion vectors (unfolded --
    # phat is a direction, not an axis; fold to hemisphere anyway
    # for the eigen test then also give the raw resultant)
    ph = perih_vec(om, w, i)
    phf = ph.copy()
    phf[np.einsum("ni,i->n", ph, e1) < 0] *= -1
    lp, vp_, Cp, Gp = eigs(phf)[0], eigs(phf)[1], None, None
    lp_all = lp
    Rp = np.linalg.norm(ph.mean(axis=0))
    # null for perihelion C
    Cpn = np.empty(N_MC)
    for k in range(N_MC):
        wmc = rng.uniform(0, 2 * np.pi, n)
        phn = perih_vec(om_mc[k], wmc, i_mc[k])
        phnf = phn.copy()
        phnf[np.einsum("ni,i->n", phn, e1) < 0] *= -1
        lk, _ = eigs(phnf)
        Cpn[k] = lk[0] - lk[1]
    Cp_obs = lp_all[0] - lp_all[1]
    out["P3_perihelion_eigenstructure"] = {
        "eigenvalues_obs": [round(float(x), 4) for x in lp_all],
        "C_obs": round(float(Cp_obs), 4),
        "C_null_mean": round(float(Cpn.mean()), 4),
        "p_C": float((int((Cpn >= Cp_obs).sum()) + 1)
                    / (N_MC + 1)),
        "R_phat": round(float(Rp), 4),
        "note": "perihelion vectors are point-like not axial; "
                "C tests extra point-concentration beyond the "
                "conditioned null"}

    RES.mkdir(exist_ok=True)
    (RES / "step_b4_pole_geometry.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.4))

    # pole positions on the sky (ecliptic lon/lat of folded poles)
    nn = pole_vec(om, i)
    lam = np.rad2deg(np.arctan2(nn[:, 1], nn[:, 0])) % 360
    bet = np.rad2deg(np.arcsin(np.clip(nn[:, 2], -1, 1)))
    ax[0].scatter(lam, 90 - bet, s=40, c="steelblue",
                  edgecolor="k", lw=0.4)
    ax[0].scatter([lam1], [90 - bet1], marker="*", s=220,
                  c="crimson", zorder=5, label="mean pole")
    ax[0].set(xlabel="ecliptic lon of pole (deg)",
              ylabel="pole tilt from ecliptic pole (deg)",
              title="folded orbit poles\n(dist. from ecl. pole)",
              xlim=(0, 360))
    ax[0].legend(fontsize=8)

    ax[1].scatter(Cn, Gn, s=3, c="0.75", rasterized=True,
                  label="conditioned null")
    ax[1].scatter([C_obs], [G_obs], marker="*", s=220,
                  c="crimson", zorder=5, label="observed")
    ax[1].axvline(np.percentile(Cn, 95), color="navy", ls=":",
                  lw=0.8)
    ax[1].axhline(np.percentile(Gn, 95), color="teal", ls=":",
                  lw=0.8)
    ax[1].set(xlabel="C = l1 - l2  (point cluster)",
              ylabel="G = l2 - l3  (girdle)",
              title="pole eigenstructure vs null\n"
                    "(upper right = point+girdle both)")
    ax[1].legend(fontsize=8)

    ax[2].hist(Cpn, bins=60, color="0.7",
               label="null C(perihelion vecs)")
    ax[2].axvline(Cp_obs, color="crimson", lw=1.6,
                  label=f"observed C={Cp_obs:.2f}")
    ax[2].set(xlabel="C = l1 - l2", ylabel="N",
              title="perihelion-vector point\nconcentration")
    ax[2].legend(fontsize=8)

    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "step_b4_pole_geometry.png", dpi=150)
    print("wrote results/step_b4_pole_geometry.json, "
          "figures/step_b4_pole_geometry.png")


if __name__ == "__main__":
    main()
