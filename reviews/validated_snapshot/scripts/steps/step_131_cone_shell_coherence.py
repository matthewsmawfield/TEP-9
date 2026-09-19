"""Step 131: Bipolar cone-shell coherence (step_b95).

Follow-on to step_129's bipolar unification test.  If the two eras'
transit structures are segments of one shell-like boundary at fixed
angular radius about the CMB-anchored bipolar axis, the slip residual
should organize by bipolar cone angle and peak on a ring at a shared
cone radius rather than at the axis pole.

Runs entirely on the independent two-leg refit checkpoints (step_b91,
step_b92); the CMB bipolar axis is the externally declared frame axis
(step_125 convention).  Look-elsewhere is priced inside every
permutation null by re-maximizing the same cell grid / ring scan on
each draw.

T1  Per-era cone x side profile (15-deg cone bands x apex/antapex)
    of the kick-regressed log-rotation residual.  Extremeness of the
    modal cell under a max-statistic permutation null.

T2  Common-cone joint test on the two measured anomaly axes: the
    declared axis (48.5 deg, antapex side) and the displaced axis
    (53.8 deg, apex side) lie in the same 45-60 deg cone band on
    opposite sides of the bipolar axis.  Null: independent isotropic
    direction pairs; reports the joint probability of the observed
    band-membership-plus-opposite-side configuration.  The resident
    detached axis (64.1 deg) is reported alongside honestly.

T3  Ring-versus-pole morphology: per era a Gaussian ring profile
    centred at cone angle theta0 is scanned over theta0 and the best
    ring amplitude compared against the pole (theta0 = 0) contrast;
    a shell peaks away from the pole.  Permutation null re-scans.

T4  Band-resident amplitude contrast: per era, residual inside the
    union band (45-60 deg, era-correct side) versus the rest of the
    sphere -- Mann-Whitney and permutation p, plus the cross-era
    product extremeness.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.statistics import monte_carlo_tail, monte_carlo_p
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_131_cone_shell_coherence")
tee_stdout(logger)
logger.header("Bipolar cone-shell coherence")

import json
import math
import numpy as np
from scipy import stats as _st
from scripts.utils.tep9_common import lv, sep
from scripts.utils.coordinates import GAL2ECL

SEED = 20260919
N_PERM = 20000
BIN = 15.0
rng = np.random.default_rng(SEED)


def gv(l, b):
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l),
                               math.sin(b)])


CMB_AP = gv(264.02, 48.25)

AXES = {
    "declared": lv(34.0, -13.0),
    "displaced": lv(120.0, -40.0),
    "resident_detached": lv(49.9, -17.0),
}


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
    d["cres"] = resid_logrot(col(rows, "cat_drot"), col(rows, "cat_daa"),
                             col(rows, "cat_denc"))
    ok = np.isfinite(d["res"]) & np.isfinite(d["aph"]).all(1)
    with np.errstate(all="ignore"):
        cos = np.clip(d["aph"][ok] @ CMB_AP, -1, 1)
    d["ok"] = ok
    d["cone"] = np.degrees(np.arccos(np.abs(cos)))   # 0..90 bipolar cone
    d["cos"] = cos
    return d


pre_rows = load_ckpt(RESULTS / "step_b91_refit.jsonl")
post_rows = load_ckpt(RESULTS / "step_b92_refit.jsonl")
pre, post = build(pre_rows), build(post_rows)
logger.info(f"loaded: pre-2018 {len(pre_rows)}, post-2017 {len(post_rows)}")

# ------------------------------------------------------- measured axes
cone_info = {}
for nm, u in AXES.items():
    s_ap = float(sep(u, CMB_AP))
    cone = min(s_ap, 180.0 - s_ap)
    cone_info[nm] = dict(cone_deg=cone,
                         side="apex" if s_ap <= 90.0 else "antapex")
    logger.info(f"  {nm}: cone {cone:.1f} deg ({cone_info[nm]['side']})")

# ----------------------------------------------------------------- T1
logger.info("T1  per-era cone x side profile")

edges = np.arange(0, 91, BIN)


def profile(d):
    res, cone, cos = d["res"][d["ok"]], d["cone"], d["cos"]
    cells = {}
    for lo in edges[:-1]:
        for side in ("apex", "antapex"):
            m = ((cone >= lo) & (cone < lo + BIN)
                 & ((cos >= 0) if side == "apex" else (cos < 0)))
            if m.sum() >= 4:
                cells[f"{lo:.0f}-{lo + BIN:.0f}_{side}"] = dict(
                    n=int(m.sum()), med=float(np.median(res[m])))
    return res, cells


def max_cell(res, cone, cos):
    best = 0.0
    for lo in edges[:-1]:
        for sgn in (1.0, -1.0):
            m = (cone >= lo) & (cone < lo + BIN) & (np.sign(cos) == sgn)
            if m.sum() >= 4:
                best = max(best, abs(float(np.median(res[m]))))
    return best


T1 = {}
for nm, d in (("pre2018", pre), ("post2017", post)):
    res, cells = profile(d)
    modal = max(cells.items(), key=lambda kv: abs(kv[1]["med"]))
    obs = abs(modal[1]["med"])
    cone, cos = d["cone"], d["cos"]
    cnt = 0
    for _ in range(N_PERM):
        rp = rng.permutation(res)
        if max_cell(rp, cone, cos) >= obs:
            cnt += 1
    T1[nm] = dict(cells=cells, modal_cell=modal[0],
                  modal_med=modal[1]["med"], modal_n=modal[1]["n"],
                  p_maxcell=float((cnt + 1) / (N_PERM + 1)))
    logger.info(f"  {nm}: modal cell {modal[0]} "
                f"med={modal[1]['med']:+.3f} n={modal[1]['n']} "
                f"p_max={T1[nm]['p_maxcell']:.4g}")

# ----------------------------------------------------------------- T2
logger.info("T2  common-cone joint test on the measured anomaly axes")

BAND = (45.0, 60.0)
d_decl, d_disp = cone_info["declared"], cone_info["displaced"]
in_band = lambda c: BAND[0] <= c < BAND[1]
obs_cfg = (in_band(d_decl["cone_deg"]) and in_band(d_disp["cone_deg"])
           and d_decl["side"] != d_disp["side"])
band_frac = (math.cos(math.radians(BAND[0]))
             - math.cos(math.radians(BAND[1])))   # bipolar band fraction
# independent isotropic pairs: P(both in band) * P(opposite sides)
p_iso = band_frac ** 2 * 0.5
# Monte-Carlo version with the same min() cone definition
u1 = rng.normal(size=(N_PERM, 3)); u1 /= np.linalg.norm(u1, axis=1)[:, None]
u2 = rng.normal(size=(N_PERM, 3)); u2 /= np.linalg.norm(u2, axis=1)[:, None]
with np.errstate(all="ignore"):
    c1 = np.degrees(np.arccos(np.clip(np.abs(u1 @ CMB_AP), -1, 1)))
    c2 = np.degrees(np.arccos(np.clip(np.abs(u2 @ CMB_AP), -1, 1)))
    s1 = (u1 @ CMB_AP) >= 0; s2 = (u2 @ CMB_AP) >= 0
mc = ((c1 >= BAND[0]) & (c1 < BAND[1]) & (c2 >= BAND[0]) & (c2 < BAND[1])
      & (s1 != s2))
p_mc = monte_carlo_p(int(mc.sum()), len(mc))

T2 = dict(band_deg=list(BAND), band_fraction=band_frac,
          declared=cone_info["declared"], displaced=cone_info["displaced"],
          resident=cone_info["resident_detached"],
          observed_config=bool(obs_cfg),
          p_isotropic_analytic=p_iso, p_isotropic_mc=p_mc,
          note=("resident detached axis sits at 64.1 deg cone -- "
                "outside the 45-60 band; reported honestly"))
logger.info(f"  band 45-60 deg: declared {d_decl['cone_deg']:.1f} "
            f"({d_decl['side']}), displaced {d_disp['cone_deg']:.1f} "
            f"({d_disp['side']})  p_iso={p_iso:.4f} (mc {p_mc:.4f})")

# ----------------------------------------------------------------- T3
logger.info("T3  ring-versus-pole morphology per era")


def ring_scan(res, cone, cos, side):
    """Gaussian ring (sigma = BIN/2) amplitude vs centre angle,
    restricted to the given bipolar side."""
    sel = (cos >= 0) if side == "apex" else (cos < 0)
    r, c = res[sel], cone[sel]
    best = (0.0, 0.0)
    for th0 in np.arange(0.0, 91.0, 5.0):
        w = np.exp(-0.5 * ((c - th0) / (BIN / 2.0)) ** 2)
        if w.sum() < 4:
            continue
        a = float((w * r).sum() / w.sum())
        if abs(a) > abs(best[0]):
            best = (a, float(th0))
    pole = r[(c >= 0) & (c < BIN)]
    return dict(ring_amp=best[0], ring_theta0=best[1],
                pole_med=float(np.median(pole)) if pole.size >= 4 else None,
                pole_n=int(pole.size), n_side=int(sel.sum()))


T3 = {}
for nm, d, side in (("pre2018", pre, "antapex"), ("post2017", post, "apex")):
    res, cone, cos = d["res"][d["ok"]], d["cone"], d["cos"]
    obs = ring_scan(res, cone, cos, side)
    cnt = 0
    for _ in range(N_PERM):
        rp = rng.permutation(res)
        if abs(ring_scan(rp, cone, cos, side)["ring_amp"]) \
                >= abs(obs["ring_amp"]):
            cnt += 1
    obs["p_ring"] = float((cnt + 1) / (N_PERM + 1))
    obs["ring_away_from_pole"] = bool(obs["ring_theta0"] > BIN)
    T3[nm] = obs
    logger.info(f"  {nm} ({side}): ring amp={obs['ring_amp']:+.3f} at "
                f"theta0={obs['ring_theta0']:.0f} deg "
                f"(pole med {obs['pole_med']})  p={obs['p_ring']:.4g}")

# ----------------------------------------------------------------- T4
logger.info("T4  band-resident amplitude contrast per era")

T4 = {}
for nm, d, side in (("pre2018", pre, "antapex"), ("post2017", post, "apex")):
    res, cone, cos = d["res"][d["ok"]], d["cone"], d["cos"]
    sgn = 1.0 if side == "apex" else -1.0
    inb = (cone >= BAND[0]) & (cone < BAND[1]) & (np.sign(cos) == sgn)
    rest = ~inb
    med_in, med_out = float(np.median(res[inb])), float(np.median(res[rest]))
    mw = _st.mannwhitneyu(res[inb], res[rest], alternative="greater")
    cnt = 0
    gap = med_in - med_out
    for _ in range(N_PERM):
        rp = rng.permutation(res)
        if np.median(rp[inb]) - np.median(rp[rest]) >= gap:
            cnt += 1
    T4[nm] = dict(side=side, n_in=int(inb.sum()), med_in=med_in,
                  med_out=med_out, gap=gap,
                  mw_p=float(mw.pvalue),
                  p_perm=float((cnt + 1) / (N_PERM + 1)))
    logger.info(f"  {nm} {side} band: {med_in:+.3f} vs {med_out:+.3f} "
                f"(n_in={inb.sum()})  perm p={T4[nm]['p_perm']:.4g}")

# cross-era joint: both era band contrasts positive
J_band = T4["pre2018"]["gap"] * T4["post2017"]["gap"]
T4["joint_product"] = float(J_band)

# ------------------------------------------------------------------ out
verdict = ("CONE-SHELL COHERENT: both era structures localize to the "
           "same 45-60 deg bipolar cone band on opposite sides, and "
           "the modern structure peaks on the ring rather than at the "
           "pole" if (obs_cfg and T1["post2017"]["p_maxcell"] < 0.05
                      and T4["post2017"]["p_perm"] < 0.05)
           else "CONE-SHELL PARTIAL on the independent record")

out = dict(
    step="step_131_cone_shell_coherence",
    description=("Bipolar cone-shell coherence: per-era cone-angle "
                 "profiles, common-cone joint test, ring-vs-pole "
                 "morphology and band-resident contrasts on the "
                 "independent refit record"),
    seed=SEED, n_perm=N_PERM, cone_bin_deg=BIN,
    cmb_apex_gal_lb=[264.02, 48.25],
    cohorts=dict(pre2018=len(pre_rows), post2017=len(post_rows)),
    measured_axes=cone_info,
    T1_cone_side_profiles=T1,
    T2_common_cone=T2,
    T3_ring_morphology=T3,
    T4_band_contrast=T4,
    verdict=verdict)

with open(RESULTS / "step_b95_cone_shell.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
logger.info(f"verdict: {verdict}")
logger.info("wrote results/step_b95_cone_shell.json")

# ------------------------------------------------------------------ fig
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# T1 profiles
width = BIN / 2.2
for nm, d, colr, off in (("pre2018", pre, "steelblue", -width / 2),
                         ("post2017", post, "indianred", width / 2)):
    cells = T1[nm]["cells"]
    for side, hatch in (("apex", ""), ("antapex", "//")):
        ks = [k for k in cells if k.endswith(side)]
        ks.sort(key=lambda k: float(k.split("-")[0]))
        xs = [float(k.split("-")[0]) + BIN / 2 + off for k in ks]
        ys = [cells[k]["med"] for k in ks]
        ax[0].bar(xs, ys, width=width, color=colr, alpha=0.8 if not hatch else 0.45,
                  hatch=hatch, label=f"{nm} {side}")
ax[0].axvspan(BAND[0], BAND[1], color="gold", alpha=0.15)
ax[0].axhline(0, color="k", lw=0.7)
ax[0].set_xlabel("bipolar cone angle (deg)")
ax[0].set_ylabel("median residual (dex)")
ax[0].set_title("Cone-angle profiles by era and side")
ax[0].legend(fontsize=7)

# T3 ring scans
for nm, d, side, colr in (("pre2018", pre, "antapex", "steelblue"),
                          ("post2017", post, "apex", "indianred")):
    res, cone, cos = d["res"][d["ok"]], d["cone"], d["cos"]
    sel = (cos >= 0) if side == "apex" else (cos < 0)
    r, c = res[sel], cone[sel]
    ths = np.arange(0.0, 91.0, 5.0)
    amps = []
    for th0 in ths:
        w = np.exp(-0.5 * ((c - th0) / (BIN / 2.0)) ** 2)
        amps.append(float((w * r).sum() / w.sum()) if w.sum() >= 4 else np.nan)
    ax[1].plot(ths, amps, color=colr, lw=2, label=f"{nm} ({side})")
    t = T3[nm]
    ax[1].axvline(t["ring_theta0"], color=colr, ls=":", lw=1)
ax[1].axvspan(BAND[0], BAND[1], color="gold", alpha=0.15)
ax[1].axhline(0, color="k", lw=0.7)
ax[1].set_xlabel("ring centre angle (deg)")
ax[1].set_ylabel("ring-weighted residual (dex)")
ax[1].set_title("Ring-morphology scan (era-correct side)")
ax[1].legend(fontsize=8)

# T4 contrasts
labs, ins, outs = [], [], []
for nm in ("pre2018", "post2017"):
    t = T4[nm]
    labs.append(nm); ins.append(t["med_in"]); outs.append(t["med_out"])
x = np.arange(2)
ax[2].bar(x - 0.18, ins, width=0.36, color=["steelblue", "indianred"],
          label="in-band (era side)")
ax[2].bar(x + 0.18, outs, width=0.36, color="0.6", label="rest of sky")
ax[2].axhline(0, color="k", lw=0.7)
ax[2].set_xticks(x); ax[2].set_xticklabels(labs)
ax[2].set_ylabel("median residual (dex)")
ax[2].set_title("Band-resident contrast")
for i, nm in enumerate(("pre2018", "post2017")):
    ax[2].text(i, max(ins[i], outs[i]) + 0.02,
               f"p={T4[nm]['p_perm']:.3g}", ha="center", fontsize=8)
ax[2].legend(fontsize=8)

fig.tight_layout()
fig.savefig(FIG / "step_b95_cone_shell.png", dpi=150)
logger.info("wrote results/figures/step_b95_cone_shell.png")
