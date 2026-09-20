"""Step 132: Bipolar meridian and polarity-mirror coherence (step_b96).

Follow-on to steps 129-131.  If the two eras' transit structures are
mirror lobes of one bipolar shell, the displaced axis should sit near
the bipolar-equator mirror of the declared axis, and all measured
anomaly directions should share a common meridian about the frame
axis.  Both relations are tested on measured directions and on the
independent per-comet residual record; every null is a permutation or
axis-shuffle that preserves the data.

T1  Polarity-mirror pair: separation between the equator-mirror of
    the declared axis (reflection through the bipolar-equatorial
    plane, the polarity-flip map) and the displaced axis.  Reported
    against the fixed-target isotropic cap probability and against a
    20,000-axis shuffle (fraction of random frame axes under which the
    mirror pair is at least as tight).

T2  Meridian coplanarity: azimuthal coherence (circular resultant R)
    of the three independently measured axes -- declared transit cap,
    detached-resident axis, displaced structure -- about the CMB
    bipolar axis.  Nulls: three isotropic directions about the fixed
    axis, and the same three axes about random bipolar axes.

T3  Meridian-resident residual: within each era on the independent
    record, the residual of comets on the era-correct bipolar side
    within 30 deg of the meridian plane (defined by the CMB axis and
    the declared axis -- both externally specified) versus the rest.
    Mann-Whitney and label permutation.

T4  Cross-channel note: detached-resident axis azimuth included in
    the coplanarity set (T2); ISO asymptotes reported descriptively.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.statistics import monte_carlo_tail, monte_carlo_p
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_132_bipolar_meridian")
tee_stdout(logger)
logger.header("Bipolar meridian and polarity-mirror coherence")

import json
import math
import numpy as np
from scipy import stats as _st
from scripts.utils.tep9_common import lv, sep
from scripts.utils.coordinates import GAL2ECL

SEED = 20260920
N_AX = 20000
N_PERM = 20000
PLANE_DEG = 30.0
rng = np.random.default_rng(SEED)


def gv(l, b):
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l),
                               math.sin(b)])


CMB = gv(264.02, 48.25)
DECL = lv(34.0, -13.0)
DISP = lv(120.0, -40.0)
DET = lv(49.9, -17.0)
AXES = dict(declared=DECL, displaced=DISP, detached=DET)


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


# ----------------------------------------------------------------- T1
logger.info("T1  polarity-mirror pair")


def eq_mirror(u, a):
    return u - 2.0 * (u @ a) * a


obs_mirror_sep = float(sep(eq_mirror(DECL, CMB), DISP))
p_fixed = (1.0 - math.cos(math.radians(obs_mirror_sep))) / 2.0

ax_sh = rng.normal(size=(N_AX, 3))
ax_sh /= np.linalg.norm(ax_sh, axis=1)[:, None]
with np.errstate(all="ignore"):
    proj = ax_sh @ DECL
    mir = 2.0 * proj[:, None] * ax_sh - DECL
    sh_sep = np.degrees(np.arccos(np.clip(mir @ DISP, -1, 1)))
p_shuffle = monte_carlo_tail(sh_sep, obs_mirror_sep, "less")
T1 = dict(mirror_sep_deg=obs_mirror_sep, p_fixed_target=p_fixed,
          p_axis_shuffle=p_shuffle, n_axes=N_AX)
logger.info(f"  mirror(declared) vs displaced: {obs_mirror_sep:.2f} deg  "
            f"p_fixed={p_fixed:.4f}  p_shuffle={p_shuffle:.4f}")

# ----------------------------------------------------------------- T2
logger.info("T2  meridian coplanarity of the measured axes")


def azimuths(us, ax):
    x = np.cross([0, 0, 1.0], ax)
    if np.linalg.norm(x) < 1e-8:
        return np.full(len(us), np.nan)
    x /= np.linalg.norm(x)
    y = np.cross(ax, x)
    return np.degrees(np.arctan2([u @ y for u in us],
                                 [u @ x for u in us])) % 360.0


def circ_R(phis_deg):
    th = np.deg2rad(np.asarray(phis_deg, dtype=float))
    return float(np.abs(np.exp(1j * th).mean()))


phi_obs = azimuths(list(AXES.values()), CMB)
R_obs = circ_R(phi_obs)

# null 1: three isotropic directions about the fixed CMB axis
u3 = rng.normal(size=(N_PERM, 3, 3))
u3 /= np.linalg.norm(u3, axis=2)[:, :, None]
R_null1 = np.array([circ_R(azimuths(list(uu), CMB)) for uu in u3])
p_iso = monte_carlo_tail(R_null1, R_obs)

# null 2: the same three axes about random bipolar axes
R_null2 = np.array([circ_R(azimuths(list(AXES.values()), a))
                    for a in ax_sh])
p_axis = float(np.nanmean(R_null2 >= R_obs))
T2 = dict(azimuths_deg={k: float(v) for k, v in zip(AXES, phi_obs)},
          circular_R=R_obs, p_isotropic_dirs=p_iso,
          p_axis_shuffle=p_axis)
logger.info(f"  azimuths {np.round(phi_obs,1)}  R={R_obs:.3f}  "
            f"p_iso={p_iso:.4f}  p_axis={p_axis:.4f}")

# ----------------------------------------------------------------- T3
logger.info("T3  meridian-resident residual on the independent record")

PLANE_N = np.cross(CMB, DECL)
PLANE_N /= np.linalg.norm(PLANE_N)

pre = build(load_ckpt(RESULTS / "step_b91_refit.jsonl"))
post = build(load_ckpt(RESULTS / "step_b92_refit.jsonl"))
logger.info(f"  loaded: pre {pre['ok'].sum()}, post {post['ok'].sum()}")

T3 = {}
for nm, d, sgn in (("pre2018", pre, -1.0), ("post2017", post, 1.0)):
    aph, res = d["aph"][d["ok"]], d["res"][d["ok"]]
    with np.errstate(all="ignore"):
        cos = aph @ CMB
        dplane = np.degrees(np.arcsin(np.clip(np.abs(aph @ PLANE_N),
                                              0, 1)))
    cell = (np.sign(cos) == sgn) & (dplane < PLANE_DEG)
    rest = ~cell
    med_in, med_out = float(np.median(res[cell])), float(np.median(res[rest]))
    mw = _st.mannwhitneyu(res[cell], res[rest], alternative="greater")
    gap = med_in - med_out
    cnt = 0
    for _ in range(N_PERM):
        rp = rng.permutation(res)
        if np.median(rp[cell]) - np.median(rp[rest]) >= gap:
            cnt += 1
    T3[nm] = dict(side="apex" if sgn > 0 else "antapex",
                  plane_deg=PLANE_DEG, n_cell=int(cell.sum()),
                  med_in=med_in, med_out=med_out, gap=gap,
                  mw_p=float(mw.pvalue),
                  p_perm=float((cnt + 1) / (N_PERM + 1)))
    logger.info(f"  {nm} {T3[nm]['side']} near-meridian: {med_in:+.3f} "
                f"vs {med_out:+.3f} (n={cell.sum()})  "
                f"perm p={T3[nm]['p_perm']:.4g}")

# ----------------------------------------------------------------- out
verdict = ("MERIDIAN-MIRROR COHERENT: the era structures form a "
           "polarity-mirror pair about the bipolar equator and the "
           "measured axes share one meridian about the CMB axis"
           if (p_shuffle < 0.05 and p_iso < 0.05)
           else "MERIDIAN-MIRROR PARTIAL on the independent record")

out = dict(
    step="step_132_bipolar_meridian",
    description=("Bipolar meridian and polarity-mirror coherence: "
                 "equator-mirror pair of the era axes, azimuthal "
                 "coplanarity of the measured axes about the CMB "
                 "bipolar axis, and meridian-resident residual cells "
                 "on the independent refit record"),
    seed=SEED, n_axes=N_AX, n_perm=N_PERM,
    cmb_apex_gal_lb=[264.02, 48.25],
    axes_used={k: "measured" for k in AXES},
    T1_mirror_pair=T1,
    T2_meridian_coplanarity=T2,
    T3_meridian_resident=T3,
    verdict=verdict)

with open(RESULTS / "step_b96_bipolar_meridian.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
logger.info(f"verdict: {verdict}")
logger.data_save(RESULTS / "step_b96_bipolar_meridian.json")

# ------------------------------------------------------------------ fig
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

ax[0].hist(sh_sep, bins=60, color="0.6")
ax[0].axvline(obs_mirror_sep, color="indianred", lw=2,
              label=f"CMB axis ({obs_mirror_sep:.1f} deg)")
ax[0].set_xlabel("mirror(declared)-to-displaced separation (deg)")
ax[0].set_ylabel("random bipolar axes")
ax[0].set_title("T1 polarity-mirror pair: axis shuffle")
ax[0].legend(fontsize=8)

ax[1].hist(R_null1, bins=60, color="0.6", label="3 isotropic dirs")
ax[1].axvline(R_obs, color="steelblue", lw=2,
              label=f"measured axes (R={R_obs:.3f})")
ax[1].set_xlabel("circular resultant R about bipolar axis")
ax[1].set_ylabel("draws")
ax[1].set_title("T2 meridian coplanarity")
ax[1].legend(fontsize=8)

for i, nm in enumerate(("pre2018", "post2017")):
    t = T3[nm]
    ax[2].bar(i - 0.18, t["med_in"], width=0.36,
              color="steelblue" if i == 0 else "indianred",
              label="era-side near-meridian" if i == 0 else None)
    ax[2].bar(i + 0.18, t["med_out"], width=0.36, color="0.6",
              label="rest of sky" if i == 0 else None)
    ax[2].text(i, max(t["med_in"], t["med_out"]) + 0.02,
               f"p={t['p_perm']:.3g}", ha="center", fontsize=8)
ax[2].axhline(0, color="k", lw=0.7)
ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["pre-2018", "post-2017"])
ax[2].set_ylabel("median residual (dex)")
ax[2].set_title("T3 era-side near-meridian cells")
ax[2].legend(fontsize=8)

fig.tight_layout()
fig.savefig(FIG / "supplementary" / "step_b96_bipolar_meridian.png", dpi=300)
logger.data_save(RESULTS / "figures/supplementary/step_b96_bipolar_meridian.png")
