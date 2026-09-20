"""Step 129: Frame-anchored bipolar unification test (step_b93).

Registered test of the single-field reading of the two-era transit
record, run entirely on the independent two-leg refit checkpoints
(step_b91 pre-2018, step_b92 post-2017) so that no Warsaw- or
JPL-lineage solution enters the observable.

The bipolar hypothesis: both eras' slip residuals are the same m=1
structure anchored on one axis, evaluated at opposite signed polarity
by era.  The CMB dipole apex is used as the externally declared frame
axis (the scalar-field rest frame already used in the GNSS, flyby and
lunar-ranging channels); no axis is selected from the data for the
registered tests.

T1  Per-era CMB-frame dipole of the kick-regressed log-rotation
    residual (covariate model identical to step_125) on our own fits,
    and identically on the catalogue channels carried in the same
    checkpoints for a lineage cross-check.  Permutation nulls.

T2  Joint opposite-polarity statistic J = -b_pre * b_post (positive
    when the two eras carry opposite-signed CMB dipoles).  Null: the
    residual labels are permuted within each era and J recomputed;
    p_J is the fraction of draws reaching J_obs.  One-sided era
    p-values are Fisher-combined as a second joint measure.

T3  Direction shuffle: J(u) evaluated on 500 random bipolar axes --
    the extremeness of the CMB axis among arbitrary axes on the
    independent record.

T4  Bipolar-cap audit: 60-deg caps about the CMB apex and antapex;
    per era the residual, per-leg deviations and leg-disagreement are
    contrasted apex-cap vs antapex-cap (characterization of which
    lobe and which leg carries each era's structure).

T5  Controls: post-2017 NG-clean subset, n_obs median strata, and
    the pre-2018 CODE-overlap carrier subset.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_129_bipolar_unification")
tee_stdout(logger)
logger.header("Frame-anchored bipolar unification test")

import json
import math
import numpy as np
from scripts.utils.tep9_common import lv, sep
from scripts.utils.coordinates import GAL2ECL

SEED = 20260919
N_PERM = 20000
N_DIRS = 500
CAP = 60.0
rng = np.random.default_rng(SEED)


def gv(l, b):
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l),
                               math.sin(b)])


CMB_AP = gv(264.02, 48.25)          # CMB dipole apex (step_125 convention)
CMB_AN = -CMB_AP


def load_ckpt(path):
    rows = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        if r.get("failed"):
            continue
        rows[r["des"]] = r            # last-wins dedupe
    return list(rows.values())


def col(rows, key):
    return np.array([r.get(key, np.nan) if r.get(key) is not None
                     else np.nan for r in rows], dtype=float)


def aph_arr(rows):
    return np.array([r.get("our_aph", [np.nan] * 3) for r in rows],
                    dtype=float)


def resid_logrot(drot, daa, denc):
    """Kick-regressed log-rotation residual, identical covariate
    model to step_125 (_resid_logrot)."""
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


def dipole_b(res, aph, u):
    ok = np.isfinite(res) & np.isfinite(aph).all(1)
    with np.errstate(all="ignore"):
        cos = np.clip(aph[ok] @ u, -1, 1)
    Xm = np.column_stack([np.ones(int(ok.sum())), cos])
    cc, *_ = np.linalg.lstsq(Xm, res[ok], rcond=None)
    return float(cc[1])


def dipole_b_se(res, aph, u):
    """OLS standard error of the dipole slope b."""
    ok = np.isfinite(res) & np.isfinite(aph).all(1)
    with np.errstate(all="ignore"):
        cos = np.clip(aph[ok] @ u, -1, 1)
    y = res[ok]
    Xm = np.column_stack([np.ones(len(y)), cos])
    cc, *_ = np.linalg.lstsq(Xm, y, rcond=None)
    resid = y - Xm @ cc
    dof = max(len(y) - 2, 1)
    s2 = float(resid @ resid) / dof
    sxx = float(((cos - cos.mean()) ** 2).sum())
    return float(np.sqrt(s2 / sxx)) if sxx > 0 else np.nan


def dipole_perm(res, aph, u, n_perm=N_PERM):
    b = dipole_b(res, aph, u)
    ok = np.isfinite(res) & np.isfinite(aph).all(1)
    rok, aok = res[ok], aph[ok]
    cnt_g = cnt_l = 0
    for _ in range(n_perm):
        bp = dipole_b(rng.permutation(rok), aok, u)
        if bp >= b:
            cnt_g += 1
        if abs(bp) >= abs(b):
            cnt_l += 1
    return dict(b=b, p_greater=float((cnt_g + 1) / (n_perm + 1)),
                p_less=float((n_perm - cnt_g + 1) / (n_perm + 1)),
                p_2sided=float((cnt_l + 1) / (n_perm + 1)),
                n=int(ok.sum()))


pre_rows = load_ckpt(RESULTS / "step_b91_refit.jsonl")
post_rows = load_ckpt(RESULTS / "step_b92_refit.jsonl")
logger.info(f"loaded: pre-2018 {len(pre_rows)}, post-2017 {len(post_rows)}")


def build(rows):
    d = dict(
        aph=aph_arr(rows),
        drot=col(rows, "our_drot"), daa=col(rows, "our_daa"),
        denc=col(rows, "our_denc"), ddirf=col(rows, "our_ddirf"),
        cdrot=col(rows, "cat_drot"), cdaa=col(rows, "cat_daa"),
        cdenc=col(rows, "cat_denc"),
        nobs=col(rows, "n_obs"))
    d["res"] = resid_logrot(d["drot"], d["daa"], d["denc"])
    d["cres"] = resid_logrot(d["cdrot"], d["cdaa"], d["cdenc"])
    return d


pre, post = build(pre_rows), build(post_rows)

# ---------------------------------------------------------------- T1
logger.info("T1  per-era CMB-frame dipoles (independent + catalogue)")

T1 = {}
for nm, d in (("pre2018", pre), ("post2017", post)):
    T1[nm] = dict(independent=dipole_perm(d["res"], d["aph"], CMB_AP),
                  catalogue=dipole_perm(d["cres"], d["aph"], CMB_AP))
    for rec in ("independent", "catalogue"):
        c = T1[nm][rec]
        logger.info(f"  {nm:<9} {rec:<11} b={c['b']:+.4f} dex  "
                    f"p_2sided={c['p_2sided']:.4g}  n={c['n']}")

# ---------------------------------------------------------------- T2
logger.info("T2  joint opposite-polarity statistic J = -b_pre*b_post")

b_pre = T1["pre2018"]["independent"]["b"]
b_post = T1["post2017"]["independent"]["b"]
J_obs = -b_pre * b_post

okp = np.isfinite(pre["res"]) & np.isfinite(pre["aph"]).all(1)
okq = np.isfinite(post["res"]) & np.isfinite(post["aph"]).all(1)
rp, ap = pre["res"][okp], pre["aph"][okp]
rq, aq = post["res"][okq], post["aph"][okq]

cnt = 0
for _ in range(N_PERM):
    bp = dipole_b(rng.permutation(rp), ap, CMB_AP)
    bq = dipole_b(rng.permutation(rq), aq, CMB_AP)
    if -bp * bq >= J_obs:
        cnt += 1
p_J = float((cnt + 1) / (N_PERM + 1))

# The sign pattern was seen in these data. Use two-sided component tests;
# the primary opposite-polarity product above already includes either polarity.
p_pre_dir = T1["pre2018"]["independent"]["p_2sided"]
p_post_dir = T1["post2017"]["independent"]["p_2sided"]
fisher = -2.0 * (math.log(max(p_pre_dir, 1e-12)) +
                 math.log(max(p_post_dir, 1e-12)))
from scipy import stats as _st
p_fisher = float(_st.chi2.sf(fisher, 4))

T2 = dict(b_pre=b_pre, b_post=b_post, J=J_obs, p_J=p_J,
          p_pre_two_sided=p_pre_dir, p_post_two_sided=p_post_dir,
          interpretation="Exploratory sign pattern; joint product allows either opposite polarity; component Fisher uses two-sided tests.",
          p_fisher=p_fisher,
          signs_opposite=bool(b_pre * b_post < 0))
logger.info(f"  b_pre={b_pre:+.4f}  b_post={b_post:+.4f}  "
            f"J={J_obs:.5f}  p_J={p_J:.4g}  p_fisher={p_fisher:.4g}")

# ---------------------------------------------------------------- T3
logger.info(f"T3  direction shuffle: J(u) on {N_DIRS} random bipolar axes")

Js = np.empty(N_DIRS)
for k in range(N_DIRS):
    v = rng.normal(size=3)
    Js[k] = -dipole_b(rp, ap, v / np.linalg.norm(v)) * \
        dipole_b(rq, aq, v / np.linalg.norm(v))
rank = int((Js >= J_obs).sum())
T3 = dict(n_dir=N_DIRS, J_cmb=J_obs,
          frac_axes_ge_cmb=float(rank / N_DIRS),
          J_max=float(Js.max()), J_med=float(np.median(Js)))
logger.info(f"  J(CMB)={J_obs:.5f}  frac random axes >= CMB: "
            f"{rank}/{N_DIRS} ({rank / N_DIRS:.3f})")

# ---------------------------------------------------------------- T4
logger.info("T4  bipolar-cap audit (60-deg caps about CMB axis)")

def cap_audit(d):
    ok = np.isfinite(d["res"]) & np.isfinite(d["aph"]).all(1)
    with np.errstate(all="ignore"):
        cos_ap = d["aph"][ok] @ CMB_AP
    res = d["res"][ok]
    out = {}
    for nm, sel in (("apex_cap", cos_ap >= math.cos(math.radians(CAP))),
                    ("antapex_cap", cos_ap <= -math.cos(math.radians(CAP))),
                    ("mid", np.abs(cos_ap) < math.cos(math.radians(CAP)))):
        out[nm] = dict(n=int(sel.sum()),
                       med_res=float(np.median(res[sel]))
                       if sel.any() else None,
                       se_med_res=float(1.2533 * np.std(res[sel], ddof=1)
                                        / np.sqrt(sel.sum()))
                       if sel.sum() > 1 else None)
    return out

T4 = {}
for nm, d, rows in (("pre2018", pre, pre_rows), ("post2017", post, post_rows)):
    cell = cap_audit(d)
    ok = np.isfinite(d["res"]) & np.isfinite(d["aph"]).all(1)
    with np.errstate(all="ignore"):
        cos_ap = d["aph"][ok] @ CMB_AP
    sel_ap = cos_ap >= math.cos(math.radians(CAP))
    sel_an = cos_ap <= -math.cos(math.radians(CAP))
    for ch in ("our_d_in_leg", "our_d_out_leg", "our_ddirf"):
        v = col(rows, ch)[ok]
        for lab, sel in (("apex", sel_ap), ("antapex", sel_an)):
            vv = v[sel]
            vv = vv[np.isfinite(vv)]
            cell[f"{ch}_{lab}"] = dict(
                n=int(vv.size),
                med=float(np.median(vv)) if vv.size else None)
    T4[nm] = cell
    logger.info(f"  {nm:<9} apex-cap res={cell['apex_cap']['med_res']} "
                f"(n={cell['apex_cap']['n']})  "
                f"antapex res={cell['antapex_cap']['med_res']} "
                f"(n={cell['antapex_cap']['n']})")

# ---------------------------------------------------------------- T5
logger.info("T5  controls")

T5 = {}
# NG-clean post-2017
ng = np.array([r.get("ng") is True for r in post_rows])
if ng.any():
    T5["post2017_ng_clean"] = dipole_perm(post["res"][~ng],
                                          post["aph"][~ng], CMB_AP)
    c = T5["post2017_ng_clean"]
    logger.info(f"  NG-clean post-2017: b={c['b']:+.4f} "
                f"p={c['p_2sided']:.4g} n={c['n']}")
else:
    T5["post2017_ng_clean"] = None
    logger.info("  NG flag absent -- subset skipped")

# CODE-overlap pre-2018 subset
incode = np.array([bool(r.get("in_code")) for r in pre_rows])
if incode.any():
    T5["pre2018_code_overlap"] = dipole_perm(pre["res"][incode],
                                           pre["aph"][incode], CMB_AP)
    c = T5["pre2018_code_overlap"]
    logger.info(f"  CODE-overlap pre-2018: b={c['b']:+.4f} "
                f"p={c['p_2sided']:.4g} n={c['n']}")

# n_obs median strata
for nm, d in (("pre2018", pre), ("post2017", post)):
    ok = np.isfinite(d["nobs"]) & np.isfinite(d["res"])
    med = np.nanmedian(d["nobs"][ok])
    sel = ok & (d["nobs"] >= med)
    T5[f"{nm}_nobs_hi"] = dipole_perm(d["res"][sel], d["aph"][sel],
                                     CMB_AP)
    c = T5[f"{nm}_nobs_hi"]
    logger.info(f"  {nm} n_obs>=med ({med:.0f}): b={c['b']:+.4f} "
                f"p={c['p_2sided']:.4g} n={c['n']}")

# ---------------------------------------------------------------- out
verdict = ("BIPOLAR-CONSISTENT: opposite-signed CMB-frame dipoles on "
           "the independent record (J extremal), same-axis reading "
           "supported" if (T2["signs_opposite"] and p_J < 0.05)
           else "BIPOLAR-NOT-ESTABLISHED on the independent record")

out = dict(
    step="step_129_bipolar_unification",
    description=("Frame-anchored bipolar unification test "
                           "on the independent two-leg refit record"),
    seed=SEED, n_perm=N_PERM, n_dirs=N_DIRS, cap_deg=CAP,
    cmb_apex_gal_lb=[264.02, 48.25],
    cohorts=dict(pre2018=len(pre_rows), post2017=len(post_rows)),
    T1_cmb_dipoles=T1, T2_joint_opposite_polarity=T2,
    T3_direction_shuffle=T3, T4_bipolar_cap_audit=T4,
    T5_controls=T5, verdict=verdict)

with open(RESULTS / "step_b93_bipolar_unification.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
logger.info(f"verdict: {verdict}")
logger.data_save(RESULTS / "step_b93_bipolar_unification.json")

# ---------------------------------------------------------------- fig
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

ax[0].hist(Js, bins=40, color="#566573", alpha=0.8)
ax[0].axvline(J_obs, color="#b43b4e", lw=2,
              label=f"J(CMB) = {J_obs:.4f}\nrank {rank}/{N_DIRS}")
ax[0].set_xlabel("J(u) = $-b_{pre}(u)\\,b_{post}(u)$")
ax[0].set_ylabel("random bipolar axes")
ax[0].legend()

# T1: signed CMB-frame dipole amplitude by era with OLS standard errors.
xs = ["pre-2018\nindep.", "pre-2018\ncat.", "post-2017\nindep.",
      "post-2017\ncat."]
bs = [T1["pre2018"]["independent"]["b"], T1["pre2018"]["catalogue"]["b"],
      T1["post2017"]["independent"]["b"], T1["post2017"]["catalogue"]["b"]]
ses = [dipole_b_se(pre["res"], pre["aph"], CMB_AP),
       dipole_b_se(pre["cres"], pre["aph"], CMB_AP),
       dipole_b_se(post["res"], post["aph"], CMB_AP),
       dipole_b_se(post["cres"], post["aph"], CMB_AP)]
ax[1].bar(range(4), bs, yerr=ses, capsize=4,
          error_kw=dict(ecolor="#17202A", elinewidth=1.2),
          color=["#1C2E4A", "#84a3aa", "#b43b4e", "#d98c96"])
ax[1].axhline(0, color="k", lw=0.7)
ax[1].set_xticks(range(4)); ax[1].set_xticklabels(xs)
ax[1].set_ylabel("$b_{\\cos\\theta_{CMB}}$ (dex)")

labels, vals, ns, ses4 = [], [], [], []
for nm in ("pre2018", "post2017"):
    for cap_ in ("apex_cap", "antapex_cap", "mid"):
        c = T4[nm][cap_]
        if c["med_res"] is not None:
            labels.append(f"{nm.replace('2018', '-18').replace('2017', '-17')}"
                          f"\n{cap_.replace('_cap', '')}")
            vals.append(c["med_res"])
            ns.append(c["n"]); ses4.append(c["se_med_res"])
cap_cols = {"apex": "#b43b4e", "antapex": "#1C2E4A",
            "mid": "#566573"}
bar_cols = [cap_cols[l.split("\n")[1]] for l in labels]
bars = ax[2].bar(range(len(vals)), vals, yerr=ses4, capsize=3,
                 error_kw=dict(ecolor="#17202A", elinewidth=1.0),
                 color=bar_cols)
for i, (rect, n) in enumerate(zip(bars, ns)):
    se = ses4[i] or 0.0
    h = rect.get_height()
    ax[2].text(rect.get_x() + rect.get_width() / 2,
               h + se + 0.012 if h >= 0 else h - se - 0.012,
               f"n={n}", ha="center",
               va="bottom" if h >= 0 else "top",
               fontsize=8, color="#17202A")
ax[2].axhline(0, color="k", lw=0.7)
ax[2].set_xticks(range(len(vals)))
ax[2].set_xticklabels(labels)
ax[2].set_ylabel("median residual (dex)")

fig.tight_layout()
fig.savefig(FIG / "step_b93_bipolar_unification.png", dpi=300)
logger.data_save(RESULTS / "figures/step_b93_bipolar_unification.png")
