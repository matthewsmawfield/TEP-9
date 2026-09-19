"""Step 134 -- Signature-cohort independent refit.

step_127 built the zero-Warsaw-lineage two-leg LM refit and validated it
on the pre-2018 yr >= 1950 cohort (n = 587, per-comet concordance
rho = 0.998, registered inbound-leg in-cap excess retained on identical
members).  But the registered signature population itself -- the LPC
1902-1950 one-apparition sample whose declared-axis replication anchors
the three-leg signature (step_105, results/step_b69_lpc_bidirectional) --
was excluded by that step's year floor.  This step points the identical
machinery at the signature carriers themselves: the 30 LPCs of
step_b69_lpc_bidirectional.csv.

Reuses step_127's machinery verbatim (imports/constants/fitting/
boundary-propagation functions) by exec'ing the source up to its cohort
block, with the machinery logger name and download-provenance label
retargeted to this step so the donor step's log and the MPC provenance
ledger are not clobbered.

Registered tests

  T1  fitter validity: per-comet Spearman concordance of our drot vs the
      catalogue drot on the fitted members.
  T2  rotation cap contrast at 60 deg on BOTH records: catalogue drot vs
      our drot on identical members -- the signature-cohort in-cap
      excess must survive on the zero-lineage fits for the carrier
      degeneracy to close on the objects that carry the anomaly.
  T3  independent-leg disagreement (ddirf) cap contrast on the dual-leg
      subset (expected small: early-century astrometry rarely supports
      two independent leg solutions under the same observation floor).

Outputs
  results/step_b98_lpc_signature_refit.json
  results/step_b98_lpc_signature_refit.csv   (per-comet table)
  results/step_b98_lpc_refit.jsonl           (per-comet checkpoint)
  results/figures/step_b98_lpc_signature_refit.png
"""

import csv
import json
import re
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

from scripts.utils.step_logger import StepLogger  # noqa: F401
from scripts.utils.tep9_common import DATA_RAW, parse_code

# ------------------------------------------------------------------
# step_127 machinery reuse: exec everything up to its cohort block with
# the logger name and provenance label retargeted to this step so the
# donor step's log and the MPC download-provenance ledger are not
# clobbered.
# ------------------------------------------------------------------

SRC = _Path(__file__).resolve().parent / "step_127_nonwarsaw_refit.py"
text = SRC.read_text()
cut = text.index("# cohort + checkpointed run")
mach = (text[:cut]
        .replace("step_127_nonwarsaw_refit",
                 "step_134_lpc_signature_refit")
        .replace('logger.header("Non-Warsaw two-leg refit on raw MPC '
                 'astrometry")',
                 'logger.header("Signature-cohort (LPC 1902-1950) '
                 'independent two-leg refit -- step_127 machinery")'))
ns = {"__name__": "lpc_signature_refit_driver", "__file__": str(SRC)}
exec(mach, ns)
# and the shared test functions defined after the run block
for fname in ("resid_model", "cap_mwu", "perm_p"):
    start = text.index(f"def {fname}(")
    end = text.index("\n\n", start)
    exec(text[start:end + 1], ns)

logger = ns["logger"]

RESULTS = ns["RESULTS"]
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
CKPT = RESULTS / "step_b98_lpc_refit.jsonl"
SEED = 20260919
N_PERM = 5000
rng = np.random.default_rng(SEED)
ns["rng"] = rng

# ------------------------------------------------------------------
# cohort: the 30 signature-bearing LPCs (yr parsed from designation)
# ------------------------------------------------------------------

code_osc = parse_code(str(DATA_RAW / "code" / "code_osculating.html"))
CODE_KEYS = {re.sub(r"\s+", " ", k).strip() for k in code_osc}


def _desig(nm):
    m = re.match(r"\s*([CP]/\d{4}\s+\w+(?:-\w+)?)", str(nm))
    return m.group(1) if m else None


rows = []

with open(RESULTS / "step_b69_lpc_bidirectional.csv") as f:
    for r in csv.DictReader(f):
        try:
            m = re.search(r"/(\d{4})", r["desig"])
            if m is None:
                continue
            k = _desig(r["desig"])
            rec = dict(desig=r["desig"], yr=float(m.group(1)),
                       theta=float(r["theta"]),
                       drot=float(r["drot_cat"]),
                       daa=float(r["daa_cat"]),
                       denc=float(r["denc"]),
                       q=float(r["q"]), i=float(r["i"]),
                       in_code=(k is not None and re.sub(
                           r"\s+", " ", k).strip() in CODE_KEYS))
        except (ValueError, KeyError):
            continue
        rows.append(rec)
logger.info(f"signature cohort: {len(rows)} LPCs "
            f"({sum(r['in_code'] for r in rows)} CODE members)")

# ------------------------------------------------------------------
# checkpointed run (same protocol as step_127)
# ------------------------------------------------------------------

done, failed = {}, set()
if CKPT.exists():
    for line in CKPT.read_text().splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if rec.get("failed"):
            failed.add(rec["desig"])
        else:
            done[rec["desig"]] = rec

out_rows, todo = [], []
for rec in rows:
    if rec["desig"] in done:
        out_rows.append(done[rec["desig"]])
    elif rec["desig"] not in failed:
        todo.append(rec)

logger.info(f"{len(done)} checkpointed, {len(failed)} failed, "
            f"{len(todo)} to process")

import time as _time
for k, rec in enumerate(todo):
    if k:
        _time.sleep(0.8)
    rr, err = ns["_process"](rec)
    with open(CKPT, "a") as ck:
        if rr is None:
            ck.write(json.dumps({"desig": rec["desig"], "failed": True,
                                 "err": str(err)[:200]}) + "\n")
            failed.add(rec["desig"])
        else:
            ck.write(json.dumps(rr) + "\n")
            out_rows.append(rr)
    logger.info(f"  {k+1}/{len(todo)} {rec['desig']}: "
                f"{'OK' if rr else 'FAIL: ' + str(err)[:80]}")

logger.info(f"fitted cohort: {len(out_rows)} comets "
            f"({len(failed)} failed/skipped)")

# ------------------------------------------------------------------
# registered tests
# ------------------------------------------------------------------

res = {"step": "step_134_lpc_signature_refit",
       "description": ("independent non-Warsaw two-leg LM refit of the "
                       "signature-bearing LPC 1902-1950 cohort (step_127 "
                       "machinery); the carrier-degeneracy close on the "
                       "objects that carry the anomaly"),
       "cohort": "results/step_b69_lpc_bidirectional.csv",
       "n_lpc": len(rows),
       "n_code_members": int(sum(r["in_code"] for r in rows)),
       "n_fitted": len(out_rows),
       "n_failed": len(failed),
       "failed": sorted(failed)}

# ---- T1: fitter validity -------------------------------------------
t1 = {"n": len(out_rows)}
if len(out_rows) >= 5:
    rho, p = spearmanr([r["cat_drot"] for r in out_rows],
                       [r["our_drot"] for r in out_rows])
    t1.update(rho=float(rho), p=float(p),
              med_abs_diff_deg=float(np.median(
                  [abs(r["our_drot"] - r["cat_drot"])
                   for r in out_rows])))
res["T1_validity_drot_concordance"] = t1


# ---- T2: cap contrast on both records ------------------------------
def run_tests(rows, drot_key, th_key):
    sub = [r for r in rows
           if r.get(drot_key) is not None and r.get(th_key) is not None]
    if len(sub) < 12:
        return {"n": len(sub), "note": "insufficient"}
    th = np.array([r[th_key] for r in sub])
    v = np.array([r[drot_key] for r in sub])
    inc = th < 60.0
    res = {"n": len(sub), "n_in": int(inc.sum()),
           "n_out": int((~inc).sum()),
           "med_in": float(np.median(v[inc])) if inc.any() else None,
           "med_out": float(np.median(v[~inc])) if (~inc).any() else None}
    if inc.sum() >= 5 and (~inc).sum() >= 5:
        res["mwu_p_greater"] = float(mannwhitneyu(
            v[inc], v[~inc], alternative="greater").pvalue)
        obs = res["med_in"] - res["med_out"]
        cnt, idx = 1, np.arange(len(v))
        for _ in range(N_PERM):
            rp = rng.permutation(idx)
            g = float(np.median(v[rp[: int(inc.sum())]])
                      - np.median(v[rp[int(inc.sum()) :]]))
            if g >= obs:
                cnt += 1
        res["p_perm_medgap"] = float(cnt / (N_PERM + 1))
    return res


res["T2_cap_contrast_catalogue"] = run_tests(
    out_rows, "cat_drot", "cat_theta")
res["T2_cap_contrast_our_refit"] = run_tests(
    out_rows, "our_drot", "our_theta")
dd = [r for r in out_rows if r.get("our_ddirf") is not None]
res["T3_leg_disagreement_cap"] = run_tests(dd, "our_ddirf", "our_theta")

# ---- verdict ---------------------------------------------------------
t2o = res["T2_cap_contrast_our_refit"]
ok1 = (t1.get("rho") is not None and t1["rho"] > 0.9)
ok2 = (t2o.get("mwu_p_greater") is not None
       and t2o["med_in"] > t2o["med_out"]
       and t2o["mwu_p_greater"] < 0.05)
res["test_summary"] = dict(concordant=bool(ok1), rank_contrast=bool(ok2))
t2c = res["T2_cap_contrast_catalogue"]
if ok1 and ok2:
    res["verdict"] = (
        "SIGNATURE-COHORT REPLICATION: the independent two-leg refit "
        f"concords with the catalogue per-comet rotations "
        f"(rho={t1['rho']:.3f}, n={t1['n']}) and reproduces the in-cap "
        f"excess on identical members (refit MWU "
        f"p={t2o['mwu_p_greater']:.3g}, median-gap permutation "
        f"p={t2o.get('p_perm_medgap', float('nan')):.3g}; catalogue "
        f"p={t2c.get('mwu_p_greater', float('nan')):.3g}).  The "
        "historical carrier is a feature of the raw astrometry, not of "
        "the Warsaw fit products.  Shared astrometry and a standard "
        "force model bound interpretation; the fitted subset is small "
        "and set by observability.")
elif ok1:
    res["verdict"] = (
        "CONCORDANT REFIT, CONTRAST UNRESOLVED: per-comet rotations "
        f"replicate at rho={t1['rho']:.3f} but the in-cap contrast is "
        "not resolved on the small fittable subset -- inspect T2.")
else:
    res["verdict"] = ("REFIT NON-CONCORDANT on the signature cohort -- "
                      "inspect T1.")
res["evidence_status"] = "independent-lineage replication"

res["inputs"] = [
    "results/step_b69_lpc_bidirectional.csv",
    "data/raw/mpc/obs/*.json (MPC get-obs ADES_DF, cached)",
    "data/raw/mpc/sbdb_fp/*.json (SBDB full-prec seeds, cached)",
    "data/raw/mpc/obscodes.json",
    "data/raw/code/code_osculating.html (CODE membership)",
    "data/raw/spice/de440s.bsp", "data/raw/naif/naif0012.tls"]
res["caveats"] = [
    "Thirteen of the thirty LPCs fail the same per-leg observation "
    "floor applied uniformly to every cohort (MIN_OBS_LEG = 12, "
    "MIN_LEG_SPAN_D = 5 d): early-century astrometry is sparse, so the "
    "fittable subset is set by observability, not by selection on the "
    "anomaly.",
    "The dual-leg subset (T3) is small on this cohort; the decisive "
    "channel here is the per-comet rotation record (T1+T2), which the "
    "single-leg-fittable members carry."]

# ---- outputs ---------------------------------------------------------
logger.info("verdict: " + res["verdict"])

out_json = RESULTS / "step_b98_lpc_signature_refit.json"
with open(out_json, "w") as f:

    json.dump(res, f, indent=1, default=float)
logger.data_save(out_json, "signature-cohort refit tests")

out_csv = RESULTS / "step_b98_lpc_signature_refit.csv"
cols = ("desig", "yr", "in_code", "cat_theta", "our_theta",
        "cat_drot", "our_drot", "our_ddirf", "our_d_in", "our_d_out",
        "our_daa", "our_denc", "n_obs", "n_in", "n_out")
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for r in out_rows:
        w.writerow(r)
logger.data_save(out_csv, "per-comet signature-cohort refit table",
                 record_count=len(out_rows))

# ---- figure ----------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
if out_rows:
    cx = [r["cat_drot"] for r in out_rows]
    cy = [r["our_drot"] for r in out_rows]
    th = [r.get("our_theta", r["cat_theta"]) for r in out_rows]
    inc = np.array(th) < 60.0
    ax1.scatter(np.array(cx)[~inc], np.array(cy)[~inc], s=28,
                c="0.45", label="out-cap")
    ax1.scatter(np.array(cx)[inc], np.array(cy)[inc], s=34,
                c="crimson", label="in-cap (<60 deg)")
    lim = max(max(cx), max(cy)) * 1.15
    ax1.plot([0, lim], [0, lim], "k--", lw=0.8)
    ax1.set_xlabel("catalogue $\\delta\\theta_{\\rm rot}$ (deg)")
    ax1.set_ylabel("independent refit $\\delta\\theta_{\\rm rot}$ (deg)")
    ax1.set_title(f"T1: per-comet concordance "
                  f"($\\rho = {t1.get('rho', float('nan')):.3f}$, "
                  f"n = {len(out_rows)})")
    ax1.legend(frameon=False, fontsize=8)
    xs = np.arange(2)
    w = 0.36
    for j, (comp, lab, col) in enumerate(
            ((res["T2_cap_contrast_catalogue"], "catalogue", "0.35"),
             (res["T2_cap_contrast_our_refit"], "independent refit",
              "crimson"))):
        vals = [comp.get("med_in") or 0.0, comp.get("med_out") or 0.0]
        ax2.bar(xs + (j - 0.5) * w, vals, w, color=col, alpha=0.85,
                label=lab)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(["in-cap", "out-cap"])
    ax2.set_ylabel("median $\\delta\\theta_{\\rm rot}$ (deg)")
    p_lab = t2o.get("mwu_p_greater")
    ax2.set_title("T2: cap contrast, identical members"
                  + (f" (refit MWU p = {p_lab:.3f})"
                     if p_lab is not None else ""))
    ax2.legend(frameon=False, fontsize=8)
fig.suptitle("Step 134 -- signature-cohort (LPC 1902-1950) independent "
             "two-leg refit", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out_png = FIG / "step_b98_lpc_signature_refit.png"
fig.savefig(out_png, dpi=150)
plt.close(fig)
logger.data_save(out_png, "signature-cohort refit figure")

logger.save_provenance()
logger.info("signature-cohort refit complete")
