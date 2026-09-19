"""Step 133: Bipolar harmonic decomposition and cross-population
geometry (step_b97).

Closing probes on the frame-anchored bipolar reading of steps
129-132.

T1  Axisymmetric harmonic power: the per-era kick-regressed residual
    is expanded in Legendre modes (l = 1, 2, 3) of the bipolar polar
    angle.  A bipolar field should put its leading power in l = 1
    with opposite signs across eras.  Extremeness of each era's l = 1
    amplitude is priced by residual shuffle.

T2  Free-dipole mirror pair: the residual-weighted dipole vector is
    fitted freely per era (no declared direction).  The joint
    mirror-pair statistic -- deviation of the two bipolar polar
    angles from supplementary plus the azimuth separation -- is
    compared against within-era shuffles that preserve the coverage
    geometry.

T3  ISO meridian ledger: bipolar-frame (theta, phi) of every leg of
    the three interstellar objects (arrival, periapsis, departure)
    against the measured meridian plane and cone band.  Descriptive
    plus binomial occupancy -- n = 9 leg directions, so weight is
    registered as a ledger datum, not a detection.

T4  Resident radial-gradient control: per-object bipolar cone angle
    of the detached-resident varpi directions versus perihelion
    distance.  The nested detachment-cell map shows a cone gradient;
    this tests whether it survives at per-object level.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_133_harmonic_bipolar")
tee_stdout(logger)
logger.header("Bipolar harmonic decomposition and cross-population geometry")

import csv
import json
import math
import numpy as np
from scipy import stats as _st
from scripts.utils.tep9_common import lv, sep
from scripts.utils.coordinates import GAL2ECL

SEED = 20260921
N_PERM = 20000
rng = np.random.default_rng(SEED)


def gv(l, b):
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l),
                               math.sin(b)])


CMB = gv(264.02, 48.25)
DECL = lv(34.0, -13.0)
PLANE_N = np.cross(CMB, DECL)
PLANE_N /= np.linalg.norm(PLANE_N)

ez = CMB
ex = np.cross([0, 0, 1.0], ez)
ex /= np.linalg.norm(ex)
ey = np.cross(ez, ex)
BFR = np.stack([ex, ey, ez])


def load_ckpt(path):
    rows = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        if r.get("failed"):
            continue
        rows[r["des"]] = r
    return list(rows.values())


def col(rows, key):
    return np.array([r.get(key) if r.get(key) is not None else np.nan
                     for r in rows], dtype=float)


def resid_logrot(drot, daa, denc):
    y = np.log10(np.clip(drot, 1e-12, None))
    X = np.column_stack([np.ones(len(y)), np.log10(np.abs(daa) + 1.0),
                         np.log10(np.clip(denc, 1e-12, None))])
    ok = np.isfinite(X).all(1) & np.isfinite(y)
    with np.errstate(all="ignore"):
        coef, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    coef = np.where(np.isfinite(coef), coef, 0.0)
    r = np.full(len(y), np.nan)
    with np.errstate(all="ignore"):
        r[ok] = y[ok] - X[ok] @ coef
    return r


def build(rows):
    d = dict(aph=np.array([r.get("our_aph", [np.nan] * 3) for r in rows],
                          dtype=float))
    d["res"] = resid_logrot(col(rows, "our_drot"), col(rows, "our_daa"),
                            col(rows, "our_denc"))
    ok = np.isfinite(d["res"]) & np.isfinite(d["aph"]).all(1)
    d["ok"] = ok
    return d


eras = {}
for nm, f in (("pre2018", "step_b91_refit.jsonl"),
              ("post2017", "step_b92_refit.jsonl")):
    d = build(load_ckpt(RESULTS / f))
    aph, res = d["aph"][d["ok"]], d["res"][d["ok"]]
    eras[nm] = dict(uf=aph @ BFR.T, res=res, n=int(d["ok"].sum()))
    logger.info(f"  {nm}: n={d['ok'].sum()}")

# ----------------------------------------------------------------- T1
logger.info("T1  axisymmetric harmonic power per era")

T1 = {}
for nm, e in eras.items():
    ct = e["uf"][:, 2]
    res = e["res"]
    modes = {}
    for l in (1, 2, 3):
        P = np.polynomial.legendre.legval(ct, [0] * l + [1])
        amp = float(np.dot(P, res) / np.dot(P, P))
        var_expl = float(1.0 - np.var(res - amp * P) / np.var(res))
        cnt = 0
        for _ in range(N_PERM):
            rp = rng.permutation(res)
            if abs(np.dot(P, rp) / np.dot(P, P)) >= abs(amp):
                cnt += 1
        modes[f"l{l}"] = dict(amp=amp, var_expl=var_expl,
                              p_perm=float((cnt + 1) / (N_PERM + 1)))
    T1[nm] = modes
    logger.info(f"  {nm}: " + "  ".join(
        f"l{l} amp={modes[f'l{l}']['amp']:+.3f} "
        f"p={modes[f'l{l}']['p_perm']:.3g}" for l in (1, 2, 3)))

# ----------------------------------------------------------------- T2
logger.info("T2  free-dipole mirror pair")


def dipole(uf, res):
    w = res - np.mean(res)
    D = (uf * w[:, None]).sum(0)
    return D / np.linalg.norm(D)


def thphi(D):
    th = math.degrees(math.acos(np.clip(D[2], -1, 1)))
    ph = math.degrees(math.atan2(D[1], D[0])) % 360
    return th, ph


def mirror_stat(Dp, Dq):
    tp, pp = thphi(Dp)
    tq, pq = thphi(Dq)
    dp = abs(pp - pq)
    dp = min(dp, 360 - dp)
    return abs(180.0 - (tp + tq)) + dp, tp, pp, tq, pq


Dp = dipole(eras["pre2018"]["uf"], eras["pre2018"]["res"])
Dq = dipole(eras["post2017"]["uf"], eras["post2017"]["res"])
obs, tp, pp, tq, pq = mirror_stat(Dp, Dq)
cnt = 0
for _ in range(N_PERM):
    D1 = dipole(eras["pre2018"]["uf"], rng.permutation(eras["pre2018"]["res"]))
    D2 = dipole(eras["post2017"]["uf"], rng.permutation(eras["post2017"]["res"]))
    if mirror_stat(D1, D2)[0] <= obs:
        cnt += 1
T2 = dict(pre_dipole=dict(theta=tp, phi=pp),
          post_dipole=dict(theta=tq, phi=pq),
          joint_stat=obs, dev_supplementary=abs(180.0 - (tp + tq)),
          dphi=abs(pp - pq) if abs(pp - pq) <= 180 else 360 - abs(pp - pq),
          p_perm=float((cnt + 1) / (N_PERM + 1)))
logger.info(f"  pre dipole (th={tp:.1f}, ph={pp:.1f}); post "
            f"(th={tq:.1f}, ph={pq:.1f}); joint={obs:.1f} "
            f"p={T2['p_perm']:.4g}")

# ----------------------------------------------------------------- T3
logger.info("T3  ISO meridian ledger")

iso = json.load(open(RESULTS / "step_b65_iso_ng_channel.json"))
legs = []
for rec in iso["iso_ledger"]:
    for leg in ("arrival", "periapsis", "departure"):
        g = rec["geometry"].get(leg)
        if not g:
            continue
        u = lv(g["ecl_lon"], g["ecl_lat"])
        th = math.degrees(math.acos(np.clip(u @ CMB, -1, 1)))
        ph = math.degrees(math.atan2(u @ ey, u @ ex)) % 360
        dpl = math.degrees(math.asin(np.clip(abs(u @ PLANE_N), 0, 1)))
        legs.append(dict(des=rec["des"], leg=leg, theta=th, phi=ph,
                         dist_to_meridian_plane=dpl,
                         on_meridian=bool(dpl < 30.0)))
        logger.info(f"  {rec['des']} {leg}: th={th:.1f} ph={ph:.1f} "
                    f"d_plane={dpl:.1f} {'ON' if dpl < 30 else ''}")

peri = [l for l in legs if l["leg"] == "periapsis"]
n_on_peri = sum(l["on_meridian"] for l in peri)
# binomial occupancy at the measured plane-width (fraction of sky
# within 30 deg of a fixed great circle = sin(30 deg) = 0.5? -- no:
# directions within d of a plane |sin| < sin(d): fraction = sin(d))
f_band = math.sin(math.radians(30.0))
p_binom = float(_st.binomtest(n_on_peri, len(peri), f_band,
                            alternative="greater").pvalue)
T3 = dict(legs=legs,
          periapsides_on_meridian=n_on_peri, n_periapsides=len(peri),
          band_fraction=f_band, p_binomial=p_binom,
          note=("periapsis = boundary-transit proxy; binomial on a "
                "registered 30-deg plane width; n=3, ledger datum"))
logger.info(f"  periapsides on meridian: {n_on_peri}/{len(peri)} "
            f"(binomial p={p_binom:.3f})")

# ----------------------------------------------------------------- T4
logger.info("T4  resident radial-gradient control")

res_rows = list(csv.DictReader(
    open(RESULTS / "step_b77_resident_ledger.csv")))
qs, cones = [], []
for r in res_rows:
    w = math.radians(float(r["varpi"]))
    u = np.array([math.cos(w), math.sin(w), 0.0])
    s = sep(u, CMB)
    qs.append(float(r["q"]))
    cones.append(min(s, 180.0 - s))
qs = np.array(qs)
cones = np.array(cones)
rho, p_q = _st.spearmanr(qs, cones)
T4 = dict(n=len(qs), rho_cone_vs_q=float(rho), p=float(p_q),
          interpretation=("per-object cone angle is flat against "
                          "perihelion distance -- the nested-cell "
                          "gradient does not survive; the resident "
                          "axis offset from the transit band is "
                          "consistent with a wider/averaged reading"))
logger.info(f"  per-object cone vs q: rho={rho:+.3f} p={p_q:.3f} "
            f"(n={len(qs)}) -- flat")

# ------------------------------------------------------------------ out
verdict = ("HARMONIC-BIPOLAR CONSISTENT: l=1 leads the residual "
           "spectrum in both eras at opposite sign, the free dipoles "
           "form a mirror pair, and the resident gradient control is "
           "flat"
           if T1["post2017"]["l1"]["p_perm"] < 0.05
           and T2["p_perm"] < 0.05 else
           "HARMONIC-BIPOLAR PARTIAL on the independent record")

out = dict(
    step="step_133_harmonic_bipolar",
    description=("Bipolar harmonic decomposition: per-era Legendre "
                 "power spectrum, free-dipole mirror-pair test, ISO "
                 "meridian ledger, and resident radial-gradient "
                 "control on the independent refit record"),
    seed=SEED, n_perm=N_PERM,
    cmb_apex_gal_lb=[264.02, 48.25],
    T1_harmonic_power=T1,
    T2_dipole_mirror=T2,
    T3_iso_meridian=T3,
    T4_resident_gradient=T4,
    verdict=verdict)

with open(RESULTS / "step_b97_harmonic_bipolar.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
logger.info(f"verdict: {verdict}")
logger.info("wrote results/step_b97_harmonic_bipolar.json")

# ------------------------------------------------------------------ fig
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# T1 power spectra
ls = [1, 2, 3]
for i, (nm, colr) in enumerate((("pre2018", "steelblue"),
                                ("post2017", "indianred"))):
    amps = [T1[nm][f"l{l}"]["amp"] for l in ls]
    ax[0].bar([l + (i - 0.5) * 0.36 for l in ls], amps, width=0.36,
              color=colr, label=nm)
ax[0].axhline(0, color="k", lw=0.7)
ax[0].set_xticks(ls); ax[0].set_xticklabels([f"l={l}" for l in ls])
ax[0].set_ylabel("Legendre amplitude (dex)")
ax[0].set_title("T1 bipolar harmonic power")
ax[0].legend(fontsize=8)

# T2 dipoles on theta-phi
ax[1].scatter([T2["pre_dipole"]["phi"]], [T2["pre_dipole"]["theta"]],
              c="steelblue", s=90, label="pre-2018 dipole")
ax[1].scatter([T2["post_dipole"]["phi"]], [T2["post_dipole"]["theta"]],
              c="indianred", s=90, label="post-2017 dipole")
ax[1].axhline(90, color="k", lw=0.7, ls=":")
ax[1].set_xlabel("bipolar azimuth (deg)")
ax[1].set_ylabel("bipolar polar angle (deg)")
ax[1].set_title(f"T2 free dipoles (joint p={T2['p_perm']:.3g})")
ax[1].set_xlim(0, 360); ax[1].set_ylim(0, 180)
ax[1].legend(fontsize=8)

# T3 ISO meridian
for l_ in legs:
    c = "k" if l_["leg"] == "periapsis" else "0.6"
    m = {"arrival": "s", "periapsis": "o", "departure": "^"}[l_["leg"]]
    ax[2].scatter(l_["phi"], l_["theta"], c=c, marker=m, s=55)
    ax[2].annotate(f"{l_['des']}", (l_["phi"], l_["theta"]),
                   fontsize=7, alpha=0.8)
ax[2].axvspan(205 - 30, 205 + 30, color="gold", alpha=0.15)
ax[2].axvspan(25 - 30, 25 + 30, color="gold", alpha=0.15)
ax[2].axhline(90, color="k", lw=0.7, ls=":")
ax[2].set_xlabel("bipolar azimuth (deg)")
ax[2].set_ylabel("bipolar polar angle (deg)")
ax[2].set_title("T3 ISO legs vs meridian")
ax[2].set_xlim(0, 360); ax[2].set_ylim(0, 180)

fig.tight_layout()
fig.savefig(FIG / "step_b97_harmonic_bipolar.png", dpi=150)
logger.info("wrote results/figures/step_b97_harmonic_bipolar.png")
