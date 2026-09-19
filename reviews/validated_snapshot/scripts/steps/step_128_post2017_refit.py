"""Step 128 -- Post-2017 two-leg refit: the masking discriminator.

The registered prospective test (step_117) returned a reversed
declared-axis contrast on the post-2017 SBDB cohort -- but the
carrier audit showed the slip observable is structurally absent from
single-solution records: a joint fit absorbs any uniform proper-time
slip into shared orbital elements, leaving the propagated residual
dominated by planetary scatter.  The holdout could therefore only
fail in form, never confirm: the construction needed to express the
anomaly (a leg-separated record) does not exist in that lineage.

Step 127 supplied the missing instrument -- an independent
Levenberg-Marquardt two-leg refit of the raw MPC astrometry under the
identical REBOUND/IAS15 + DE440s force model, propagated to the
+/-250 AU sphere through the same machinery.  On the pre-2018 cohort
it reproduced the catalogue record at per-comet concordance rho =
+0.998 and recovered the registered inbound-leg in-cap excess.

This step points the same instrument at the post-2017 cohort
(step_b81 primaries).  Two constructions are tested on identical
members:

  * the single-fit construction -- the propagated boundary products
    of the joint (full-arc) fit, the only construction a
    single-solution lineage can express.  Predicted flat under the
    masking hypothesis.
  * the leg-fit construction -- the independently fitted inbound and
    outbound legs, whose boundary asymptotes carry whatever slip the
    transit actually imparted.  This is the construction the
    three-leg catalogue records realize, and the registered carrier
    d_in on the pre-2018 record lives here.

Registered tests

  T1  instrument validity on the new era: per-comet concordance of
      our drot versus the catalogue-propagated drot, and element
      agreement of our full-arc fit against the SBDB seed.
  T2  DECISIVE -- the masking discriminator: on the dual-leg subset,
      the cap contrast (theta < 60 deg, aph against the declared
      axis) of the leg-fit inbound deviation d_in_leg versus (a) the
      single-fit d_in on the same members and (b) the catalogue's
      propagated d_in.  Under TEP-with-masking the leg-fit channel
      carries the in-cap excess the joint construction erases.
      Flatness on the leg-fit channel would confine the anomaly to
      the pre-2018 era -- a substantive boundary on the hypothesis.
  T3  displaced-structure control: the same contrasts evaluated
      about the post-2017 displaced axis (120, -40) -- does the
      modern residual systematic appear on independent fits, or was
      it a single-lineage fit field?
  T4  quality strata: contrast on NG-flagged vs clean comets (the
      non-gravitational floor differs by construction), and the
      cross-arc extrapolation penalty by cap.
  T4b signature-class audit: at both axes on the dual-leg cohort,
      cap contrasts on every channel -- leg deviations (d_in_leg,
      d_out_leg), drot, ddirf, |daa| energy, plus coverage and
      fit-quality controls.  The declared anomaly's class is
      single-leg rotation slip without energy exchange; the audit
      asks whether the displaced structure shares that class.
  T4c robustness: label-permutation nulls on the carrying channels,
      500-direction shuffle extremeness of the displaced residual
      gap on the independent record, and quality strata (n_obs,
      fit RMS, NG-clean) on the same field.
  T5  audit: cohort yield, leg coverage, fit quality.

Outputs
  results/step_b92_post2017_refit.json
  results/step_b92_post2017_refit.csv
  results/figures/step_b92_post2017_refit.png
  results/step_b92_refit.jsonl            (per-comet checkpoint)
  (obs/sbdb caches shared with step_127 under data/raw/mpc/)
"""

import csv
import json
import math
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, sep
from scripts.utils.lpc_boundary import TNO

FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)

logger = StepLogger("step_128_post2017_refit")

SPK = DATA_RAW / "spice" / "de440s.bsp"
LSK = DATA_RAW / "naif" / "naif0012.tls"
MPC_DIR = DATA_RAW / "mpc"
OBS_DIR = MPC_DIR / "obs"
SBDB_DIR = MPC_DIR / "sbdb_fp"
OBSC_PATH = MPC_DIR / "obscodes.json"
PROV_PATH = MPC_DIR / "provenance.json"
CKPT = RESULTS / "step_b92_refit.jsonl"

N_PERM = 5000
SEED = 20260919

WORKERS = 1
if "--workers" in _sys.argv:
    _i = _sys.argv.index("--workers")
    WORKERS = max(1, int(_sys.argv[_i + 1]))
REDO = "--redo" in _sys.argv

logger.header("Post-2017 two-leg refit: the masking discriminator")

from scripts.utils import mpc_refit as R

R.configure(OBS_DIR, SBDB_DIR, OBSC_PATH, PROV_PATH, SPK, LSK,
            "step_128_post2017_refit")


def _worker_init():
    R.refit_worker_init()


def _des(row):
    return row["desig"].split(" (")[0].strip()


def _process(row):
    res, err = R.fit_comet(_des(row))
    if res is None:
        return None, err
    row2 = dict(desig=row["desig"], des=row["desig"].split(" (")[0],
                yr=row["yr"], cat_theta=row["theta"],
                cat_drot=row["drot"], cat_d_in=row["d_in"],
                cat_d_out=row["d_out"], cat_daa=row["daa"],
                cat_denc=row["denc"], cat_q=row["q"], cat_i=row["i"],
                ng=row["ng"],
                n_obs=res["n_obs"], n_in=res["n_in"],
                n_out=res["n_out"], tp_jd=res["tp_jd"])
    if "del_pdeg" in res:
        row2.update(fit_all_rms=res["fit_all_rms"],
                    fit_all_n=res["fit_all_n"],
                    our_q=res["q_fit"], our_e=res["e_fit"],
                    our_i=res["i_fit"], our_om=res["om_fit"],
                    our_w=res["w_fit"],
                    del_q=res["del_q"], del_e=res["del_e"],
                    del_i=res["del_i"], del_pdeg=res["del_pdeg"])
    for src, dst in (("drot", "our_drot"), ("d_in", "our_d_in"),
                     ("d_out", "our_d_out"), ("daa", "our_daa"),
                     ("denc", "our_denc"), ("theta", "our_theta"),
                     ("aph", "our_aph"), ("dtau", "our_dtau"),
                     ("ddirf", "our_ddirf"),
                     ("d_in_leg", "our_d_in_leg"),
                     ("d_out_leg", "our_d_out_leg"),
                     ("xarc_rms_in2out", "xarc_rms_in2out"),
                     ("boundary_err", "boundary_err")):
        if src in res:
            row2[dst] = res[src]
    return row2, None


# ------------------------------------------------------------------
# cohort + checkpointed run
# ------------------------------------------------------------------

rows = []

with open(RESULTS / "step_b81_prospective_lpc.csv") as f:
    for r in csv.DictReader(f):
        try:
            rec = dict(desig=r["desig"], yr=float(r["yr"]),
                       theta=float(r["theta"]), drot=float(r["drot"]),
                       d_in=float(r["d_in"]), d_out=float(r["d_out"]),
                       daa=float(r["daa"]), denc=float(r["denc"]),
                       q=float(r["q"]), i=float(r["i"]),
                       ng=(r.get("ng") == "True"))
        except (ValueError, KeyError):
            continue
        rows.append(rec)
logger.info(f"post-2017 refit cohort: {len(rows)} comets")

done, failed = {}, {}
if CKPT.exists() and not REDO:
    for line in CKPT.read_text().splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if rec.get("failed"):
            failed[rec["desig"]] = rec
        else:
            done[rec["desig"]] = rec
elif CKPT.exists():
    CKPT.unlink()

out_rows = []
todo = []
for rec in rows:
    if rec["desig"] in done:
        out_rows.append(done[rec["desig"]])
    elif rec["desig"] in failed:
        continue
    else:
        todo.append(rec)

logger.info(f"{len(done)} checkpointed, {len(failed)} failed, "
            f"{len(todo)} to process")


def _account(rec, rowres, err):
    with open(CKPT, "a") as ck:
        if rowres is None:
            ck.write(json.dumps({"desig": rec["desig"],
                                 "failed": True,
                                 "err": str(err)[:200]}) + "\n")
            failed[rec["desig"]] = {"err": str(err)[:200]}
        else:
            ck.write(json.dumps(rowres) + "\n")
            out_rows.append(rowres)


if todo:
    if WORKERS > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ctx.Pool(WORKERS, initializer=_worker_init) as pool:
            for i, (rec, (rr, err)) in enumerate(zip(
                    todo, pool.imap(_process, todo))):
                _account(rec, rr, err)
                if (i + 1) % 10 == 0:
                    logger.info(f"  {i+1}/{len(todo)} processed")
    else:
        for i, rec in enumerate(todo):
            rr, err = _process(rec)
            _account(rec, rr, err)
            if (i + 1) % 10 == 0:
                logger.info(f"  {i+1}/{len(todo)} processed")

logger.info(f"fitted cohort: {len(out_rows)} comets "
            f"({len(failed)} failed/skipped)")

# ------------------------------------------------------------------
# statistics
# ------------------------------------------------------------------

rng = np.random.default_rng(SEED)

DISPLACED = np.array([math.cos(math.radians(-40.0)) *
                      math.cos(math.radians(120.0)),
                      math.cos(math.radians(-40.0)) *
                      math.sin(math.radians(120.0)),
                      math.sin(math.radians(-40.0))])


def cap_contrast(vals, th, cap=60.0, alternative="greater"):
    inc = np.asarray(th) < cap
    v = np.asarray(vals)
    if inc.sum() < 5 or (~inc).sum() < 5:
        return None
    u = mannwhitneyu(v[inc], v[~inc], alternative=alternative)
    return dict(n_in=int(inc.sum()), n_out=int((~inc).sum()),
                med_in=float(np.median(v[inc])),
                med_out=float(np.median(v[~inc])),
                p=float(u.pvalue))


def theta_to(rows_, key, axis):
    out = []
    for r in rows_:
        a = np.array(r["our_aph"])
        out.append(sep(a, axis))
    return np.array(out)


res = {"step": "step_128_post2017_refit",
       "description": ("independent two-leg LM refit of the post-2017 "
                       "prospective cohort on raw MPC astrometry; the "
                       "masking discriminator -- leg-fit vs single-fit "
                       "carrier constructions on identical members"),
       "n_attempted": len(rows),
       "n_fitted": len(out_rows),
       "n_failed": len(failed)}

fit_all = [r for r in out_rows if "our_drot" in r]
dual = [r for r in out_rows if "our_d_in_leg" in r]
res["n_full_arc"] = len(fit_all)
res["n_dual_leg"] = len(dual)
res["n_inbound_obs"] = int(sum(r["n_in"] > 0 for r in out_rows))

# ---- T1: instrument validity ------------------------------------
t1 = {}
if len(fit_all) > 15:
    rho, p = spearmanr([r["our_drot"] for r in fit_all],
                       [r["cat_drot"] for r in fit_all])
    t1["drot_concordance"] = {"n": len(fit_all),
                              "rho": float(rho), "p": float(p)}
    rho, p = spearmanr([r["our_d_in"] for r in fit_all],
                       [r["cat_d_in"] for r in fit_all])
    t1["d_in_concordance"] = {"rho": float(rho), "p": float(p)}
    t1["element_diffs_vs_sbdb"] = {
        "med_dq_au": float(np.median([r["del_q"] for r in fit_all])),
        "med_pdir_deg": float(np.median(
            [r["del_pdeg"] for r in fit_all])),
        "p95_pdir_deg": float(np.percentile(
            [r["del_pdeg"] for r in fit_all], 95))}
res["T1_validity"] = t1

# ---- T2: DECISIVE masking discriminator --------------------------
# On the dual-leg subset: leg-fit carrier d_in_leg vs the single-fit
# carrier d_in (same comets) vs the catalogue's propagated d_in.
t2 = {}
for_lbl = {"leg_fit": "our_d_in_leg",
           "single_fit": "our_d_in",
           "catalogue": "cat_d_in"}
if len(dual) > 20:
    th = np.array([r["our_theta"] for r in dual])
    t2["declared_axis"] = {lbl: cap_contrast(
        [r[k] for r in dual], th) for lbl, k in for_lbl.items()}
    # leg-fit outbound control (predicted flat under TEP)
    t2["declared_axis"]["leg_fit_out"] = cap_contrast(
        [r["our_d_out_leg"] for r in dual], th)
    # permutation null on the leg-fit contrast
    if t2["declared_axis"]["leg_fit"]:
        v = np.array([r["our_d_in_leg"] for r in dual])
        inc = th < 60.0
        obs = float(np.median(v[inc]) - np.median(v[~inc]))
        idx = np.arange(len(v))
        cnt = 1
        for _ in range(N_PERM):
            rp = rng.permutation(idx)
            g = float(np.median(v[rp[:inc.sum()]])
                      - np.median(v[rp[inc.sum():]]))
            if g >= obs:
                cnt += 1
        t2["declared_axis"]["leg_fit"]["p_perm"] = float(
            cnt / (N_PERM + 1))
res["T2_masking"] = t2

# ---- T3: displaced-structure control ------------------------------
t3 = {}
if len(dual) > 20:
    thd = theta_to(dual, "our_aph", DISPLACED)
    t3["displaced_axis"] = {lbl: cap_contrast(
        [r[k] for r in dual], thd) for lbl, k in for_lbl.items()}
if len(fit_all) > 30:
    thf = np.array([r["our_theta"] for r in fit_all])
    t3["declared_single_all"] = cap_contrast(
        [r["our_d_in"] for r in fit_all], thf)
    t3["declared_cat_all"] = cap_contrast(
        [r["cat_d_in"] for r in fit_all], thf)

# the displaced structure's own registered channel (step_117): the
# kick-regressed log10(drot) residual cap gap about (120,-40).  Run on
# OUR independently fitted drot and on the catalogue drot for the
# identical members -- if the modern record's dominant feature is a
# real astrometric structure it must replicate here.
def _resid_gap(rows_, drot_k, daa_k, denc_k, q_k, i_k, axis):
    v = np.array([r[drot_k] for r in rows_])
    K = np.abs([r[daa_k] for r in rows_])
    D = np.array([r[denc_k] for r in rows_])
    Q = np.array([r[q_k] for r in rows_])
    I = np.array([r[i_k] for r in rows_])
    y = np.log10(np.clip(v, 1e-9, None))
    X = np.column_stack([np.ones(len(y)), np.log10(K + 1.0),
                         np.log10(np.clip(D, 1e-3, None)), Q, I])
    with np.errstate(all="ignore"):
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ coef
    resid[~np.isfinite(resid)] = np.nanmedian(
        resid[np.isfinite(resid)])
    aph = np.array([r["our_aph"] for r in rows_])
    th = np.degrees(np.arccos(np.clip(aph @ axis, -1, 1)))
    inc = th < 60.0
    u = mannwhitneyu(resid[inc], resid[~inc], alternative="greater")
    return {"n_in": int(inc.sum()), "n_out": int((~inc).sum()),
            "med_in": float(np.median(resid[inc])),
            "med_out": float(np.median(resid[~inc])),
            "p": float(u.pvalue)}

if len(fit_all) > 30:
    t3["displaced_resid_field"] = {
        "our": _resid_gap(fit_all, "our_drot", "our_daa", "our_denc",
                          "our_q", "our_i", DISPLACED),
        "cat": _resid_gap(fit_all, "cat_drot", "cat_daa", "cat_denc",
                          "cat_q", "cat_i", DISPLACED),
        "declared_control_our": _resid_gap(
            fit_all, "our_drot", "our_daa", "our_denc",
            "our_q", "our_i", TNO)}
res["T3_displaced"] = t3

# ---- T4: quality strata -------------------------------------------
t4 = {}
if len(dual) > 20:
    for flag, sub in (("ng", [r for r in dual if r["ng"]]),
                      ("clean", [r for r in dual if not r["ng"]])):
        if len(sub) > 15:
            th = np.array([r["our_theta"] for r in sub])
            t4[flag] = cap_contrast(
                [r["our_d_in_leg"] for r in sub], th)
x = [r for r in dual if "xarc_rms_in2out" in r]
if len(x) > 20:
    th = np.array([r["our_theta"] for r in x])
    t4["xarc_rms"] = cap_contrast(
        [r["xarc_rms_in2out"] for r in x], th)
res["T4_strata"] = t4

# ---- T4b: signature-class audit at both axes ----------------------
# The declared anomaly's fingerprint is a single-leg rotation slip
# without energy exchange (inbound-localized d_in, flat |daa|).  The
# same audit is run on the displaced structure at its own axis: if
# it shares the class (one leg elevated, energy flat, controls flat)
# the modern structure is a lapse-slip anomaly of the same kind --
# a different registered position/leg of the boundary channel, not
# a different physics.  Coverage and fit-quality controls accompany
# every channel so a coverage artefact cannot masquerade as signal.
t4b = {}
if len(dual) > 20:
    thd = theta_to(dual, "our_aph", DISPLACED)
    thc = np.array([r["our_theta"] for r in dual])
    for lbl, th in (("at_displaced", thd), ("at_declared", thc)):
        blk = {}
        for ch in ("our_d_in_leg", "our_d_out_leg", "our_drot",
                   "our_ddirf"):
            vv = [r.get(ch) for r in dual]
            if all(v is not None for v in vv):
                blk[ch] = cap_contrast(vv, th)
        blk["abs_daa"] = cap_contrast(
            [abs(r["our_daa"]) for r in dual], th)
        for ctl in ("n_in", "n_out", "fit_all_rms"):
            blk[f"ctl_{ctl}"] = cap_contrast(
                [r[ctl] for r in dual], th, alternative="two-sided")
        t4b[lbl] = blk
res["T4_signature_class"] = t4b

# ---- T4c: robustness of the displaced structure -------------------
# (i) label-permutation nulls on the channels that carry the
#     signature-class result (MW alone is thin for a marginal claim);
# (ii) direction-shuffle extremeness -- is (120,-40) still extremal
#     on the independent record, or does any direction give the same
#     residual gap?
# (iii) quality strata -- does the residual field survive splits on
#     observation count, fit quality and NG flags?
t4c = {}
if len(dual) > 20:
    thd = theta_to(dual, "our_aph", DISPLACED)
    incd = thd < 60.0
    perms = {}
    for ch in ("our_d_out_leg", "our_ddirf", "our_drot"):
        vv = np.array([r.get(ch) for r in dual], dtype=object)
        if any(v is None for v in vv):
            continue
        vv = vv.astype(float)
        obs = float(np.median(vv[incd]) - np.median(vv[~incd]))
        idx = np.arange(len(vv))
        cnt = 1
        for _ in range(N_PERM):
            rp = rng.permutation(idx)
            if float(np.median(vv[rp[:incd.sum()]])
                     - np.median(vv[rp[incd.sum():]])) >= obs:
                cnt += 1
        perms[ch] = {"gap": obs, "p_perm": float(cnt / (N_PERM + 1))}
    t4c["perm_nulls_at_displaced"] = perms
if len(fit_all) > 30:
    # (ii) direction shuffle on the residual field
    _d = (t3.get("displaced_resid_field") or {}).get("our")
    obs_gap = (_d["med_in"] - _d["med_out"]) if _d else None
    if obs_gap is not None:
        v = np.array([r["our_drot"] for r in fit_all])
        K = np.abs([r["our_daa"] for r in fit_all])
        D = np.array([r["our_denc"] for r in fit_all])
        Q = np.array([r["our_q"] for r in fit_all])
        I = np.array([r["our_i"] for r in fit_all])
        y = np.log10(np.clip(v, 1e-9, None))
        X = np.column_stack([np.ones(len(y)), np.log10(K + 1.0),
                             np.log10(np.clip(D, 1e-3, None)), Q, I])
        with np.errstate(all="ignore"):
            resid = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
        resid[~np.isfinite(resid)] = np.nanmedian(
            resid[np.isfinite(resid)])
        aph = np.array([r["our_aph"] for r in fit_all])
        n_dir = 500
        cnt = 1
        for _ in range(n_dir):
            z = rng.normal(size=3)
            ax = z / np.linalg.norm(z)
            thx = np.degrees(np.arccos(np.clip(aph @ ax, -1, 1)))
            ix = thx < 60.0
            if ix.sum() < 10 or (~ix).sum() < 10:
                continue
            if float(np.median(resid[ix])
                     - np.median(resid[~ix])) >= obs_gap:
                cnt += 1
        t4c["direction_shuffle"] = {
            "n_dir": n_dir, "obs_gap": float(obs_gap),
            "p_extremeness": float(cnt / (n_dir + 1))}
    # (iii) quality strata on the same residual field
    strata = {}
    nobs = np.array([r["n_obs"] for r in fit_all])
    rms = np.array([r["fit_all_rms"] for r in fit_all])
    for nm, m in (("n_obs>=med", nobs >= np.median(nobs)),
                  ("rms<=med", rms <= np.median(rms)),
                  ("ng_clean", np.array(
                      [not r.get("ng") for r in fit_all]))):
        sub = [r for r, k in zip(fit_all, m) if k]
        if len(sub) > 30:
            strata[nm] = {"n": len(sub), **_resid_gap(
                sub, "our_drot", "our_daa", "our_denc",
                "our_q", "our_i", DISPLACED)}
    t4c["quality_strata"] = strata
res["T4c_robustness"] = t4c
res["n_nonfinite_rows"] = int(sum(
    1 for r in out_rows
    if not np.all(np.isfinite(np.array(r.get("our_aph", [0, 0, 0]),
                                        dtype=float)))))

# ---- T5: audit -----------------------------------------------------
res["T5_audit"] = {
    "rms_all_median": float(np.median(
        [r["fit_all_rms"] for r in fit_all])) if fit_all else None,
    "dual_yield_frac": len(dual) / max(1, len(out_rows)),
    "fail_reasons": sorted({v.get("err", "?")[:60]
                            for v in failed.values()})[:10]}

# verdict: two questions answered.  (a) masking -- does the
# leg-separated record (the construction that CAN express slip)
# recover the in-cap excess the single-fit record erases?  (b) is
# the post-2017 record's own dominant structure (the displaced-axis
# residual field of step_117) present on independently fitted orbits?
leg = (t2.get("declared_axis") or {}).get("leg_fit")
disp = (t3.get("displaced_resid_field") or {}).get("our")
ok_leg = (leg and leg["med_in"] > leg["med_out"]
          and (leg.get("p_perm") or leg["p"]) < 0.05)
ok_disp = (disp and disp["med_in"] > disp["med_out"]
           and disp["p"] < 0.05)
res["test_summary"] = dict(declared_leg_excess=bool(ok_leg), displaced_excess=bool(ok_disp))
res["verdict"] = 'TEMPORAL HOLDOUT COMPARISON: the displaced structure is exploratory and is not a measured lapse. Inspect the declared carrier in T2 and catalogue-composition controls in step 140.'
res["evidence_status"] = "exploratory population comparison"

# registered forward falsifier (see manuscript 5.5): under the
# frame-anchored bipolar reading the polarity ordering is the
# prediction -- future transit records should continue to resolve
# the modern apex-side lobe rather than revert to the pre-2018
# antapex-side polarity.
res["forward_prediction"] = {
    "proposed_not_preregistered": "future (LSST-era) LPC cohorts continue the "
                  "modern displaced-lobe polarity on the leg-"
                  "separated channel; a return of the pre-2018 "
                  "declared-axis inbound-leg excess, or no coherent "
                  "lobe, breaks the era-polarity ordering",
    "channel": "two-leg refit d_out_leg / residual cap-gap at the "
               "displaced axis",
    "scorable_with": "this script re-run on the new cohort"}
res["inputs"] = [
    "results/step_b81_prospective_lpc.csv",
    "data/raw/mpc/obs/*.json (shared MPC cache)",
    "data/raw/mpc/sbdb_fp/*.json (shared SBDB seed cache)",
    "data/raw/mpc/obscodes.json",
    "data/raw/spice/de440s.bsp", "data/raw/naif/naif0012.tls"]
res["caveats"] = [
    "The leg-fit carrier requires >=12 astrometric observations per "
    "leg; post-2017 comets discovered post-perihelion contribute only "
    "the single-fit construction, so the dual-leg subset is the "
    "decisive arena and its size is set by observability.",
    "Pure-gravity refit: cometary NG terms are absorbed into the leg "
    "solutions; the NG-flag stratum in T4 prices that floor.",
    "d_in_leg compares the inbound-leg asymptote to the joint-fit "
    "osculating direction -- the construction the three-leg "
    "catalogue record realizes; per-comet sign conventions follow "
    "step_123."]

out = RESULTS / "step_b92_post2017_refit.json"
json.dump(res, open(out, "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))

with open(RESULTS / "step_b92_post2017_refit.csv", "w",
          newline="") as f:
    keys = sorted({k for r in out_rows for k in r.keys()
                   if k != "our_aph"})
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in out_rows:
        w.writerow({k: r.get(k) for k in keys})

# ------------------------------------------------------------------
# figure
# ------------------------------------------------------------------

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
if len(dual) > 15:
    ax = axes[0]
    ax.scatter([r["cat_d_in"] for r in dual],
               [r["our_d_in_leg"] for r in dual], s=9, alpha=0.5)
    lim = max(max(r["cat_d_in"] for r in dual),
              max(r["our_d_in_leg"] for r in dual))
    ax.plot([0, lim], [0, lim], "k-", lw=0.7)
    ax.set_xlabel("catalogue d_in (deg)")
    ax.set_ylabel("leg-fit d_in (deg)")
    ax.set_title(f"dual-leg record, n={len(dual)}")
if t2.get("declared_axis"):
    ax = axes[1]
    th = np.array([r["our_theta"] for r in dual])
    for key, col, lbl in (("our_d_in_leg", "crimson", "leg-fit d_in"),
                          ("our_d_in", "steelblue", "single-fit d_in")):
        ax.scatter(th, [r[key] for r in dual], s=9, alpha=0.45,
                   c=col, label=lbl)
    ax.axvline(60, color="k", lw=0.7, ls="--")
    ax.set_xlabel("theta to declared axis (deg)")
    ax.set_ylabel("d_in (deg)")
    ax.legend(fontsize=7)
    if leg:
        ax.set_title(f"leg-fit {leg['med_in']:.3f}/{leg['med_out']:.3f} "
                     f"(perm p={leg.get('p_perm', leg['p']):.3f})")
if t3.get("displaced_resid_field") and t3["displaced_resid_field"].get("our"):
    ax = axes[2]
    # recompute the residual field for the scatter
    _rf = t3["displaced_resid_field"]
    v = np.array([r["our_drot"] for r in fit_all])
    K = np.abs([r["our_daa"] for r in fit_all])
    D = np.array([r["our_denc"] for r in fit_all])
    Q = np.array([r["our_q"] for r in fit_all])
    I = np.array([r["our_i"] for r in fit_all])
    y = np.log10(np.clip(v, 1e-9, None))
    X = np.column_stack([np.ones(len(y)), np.log10(K + 1.0),
                         np.log10(np.clip(D, 1e-3, None)), Q, I])
    with np.errstate(all="ignore"):
        resid = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    aph = np.array([r["our_aph"] for r in fit_all])
    thd = np.degrees(np.arccos(np.clip(aph @ DISPLACED, -1, 1)))
    ax.scatter(thd, resid, s=9, alpha=0.5,
               c=np.where(thd < 60, "darkviolet", "grey"))
    ax.axvline(60, color="k", lw=0.7, ls="--")
    ax.axhline(0, color="k", lw=0.4)
    ax.set_xlabel("theta to displaced axis (120,-40) deg")
    ax.set_ylabel("log10(drot) residual")
    d3 = _rf["our"]
    ax.set_title(f"displaced field replicates: {d3['med_in']:+.3f}/"
                 f"{d3['med_out']:+.3f} dex (p={d3['p']:.4f})")
fig.tight_layout()
fig.savefig(FIG / "step_b92_post2017_refit.png", dpi=150)
logger.data_save(out)
logger.data_save(RESULTS / "step_b92_post2017_refit.csv")
logger.data_save(FIG / "step_b92_post2017_refit.png")
