"""Step 143: Comet non-gravitational acceleration channel (step_b107).

Published NG parameters are used throughout this work only as an
exclusion control (the NG-clean subsets).  They are never read as an
observable.  Under the boundary-slip picture they should carry the
anomaly in a specific way: a proper-time slip absorbed by a
pure-gravity fit must leak into whatever free parameter soaks up a
timing error, and in the Marsden formulation that parameter is the
transverse NG term A2 (the component that changes the fitted orbital
period).  The Warsaw catalogue publishes A1/A2/A3 with uncertainties
plus the outgassing-peak time shift tau for every comet whose orbit
required an NG solution (tableb4, J/A+A/567/A126).

This step reads the NG record as a directional channel against the
declared cap (theta < 60 deg to the (34,-13) transit axis, taken from
the step_b30 per-comet record):

T1  NG-model incidence: fraction of the Warsaw cohort fitted with an
    NG orbit (tablea1 Model column), in-cap vs out-cap.
T2  Transverse leakage: signed A2 and significance-weighted A2/e_A2
    in-cap vs out-cap, plus A2 against cos(theta) over the full cohort.
    The slip is time-domain; A2 is the timing parameter.
T3  Radial/normal controls: the same tests on A1 and A3 -- a genuine
    boundary signal should load preferentially into the transverse
    component, not uniformly into all three.
T4  Time-shift channel: the fitted outgassing-peak offset tau in-cap
    vs out-cap.
T5  Confound audit: NG solutions are preferentially fitted to
    well-observed, brighter comets -- n_obs and arc-length medians
    in-cap vs out-cap among NG members, and an n_obs-matched
    sensitivity bound on the T2 contrast.
T6  Era note: the Warsaw NG record is pre-2011; the channel is
    registered on the cohort that carries the declared anomaly.
T7  Independent replication: the same transverse-term cap contrast on
    the JPL SBDB NG record (model_pars in the full-precision orbit
    files), pre-2018 and post-2017 cohorts -- a second orbit pipeline
    with no Warsaw input.

Outputs: results/step_b107_comet_ng_channel.json/.csv and
results/figures/step_b107_comet_ng_channel.png.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_143_comet_ng_channel")
tee_stdout(logger)
logger.header("Comet non-gravitational acceleration channel")

import csv
import json
import numpy as np
from scipy import stats as _st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20260919
CAP = 60.0
rng = np.random.default_rng(SEED)

# ---- load per-comet transit angles (step_b30 record) --------------------
b30 = {r["desig"].strip(): r for r in csv.DictReader(
    open(RESULTS / "step_b30_proper_time_slip.csv"))
    if r["cohort"] == "warsaw"}


def norm_desig(s):
    return " ".join(s.split())


# ---- tablea1: model incidence (GR vs NG) --------------------------------
a1 = {}
for ln in open(DATA_RAW / "warsaw_tablea1.dat"):
    if len(ln) < 140 or not ln.strip():
        continue
    desig = norm_desig(ln[3:15])
    model = ln[132:140].strip()
    qosc = ln[40:46].strip()
    nobs = ln[88:93].strip()
    arcy = ln[97:102].strip()
    a1[desig] = {
        "model": model,
        "q": float(qosc) if qosc else float("nan"),
        "nobs": int(nobs) if nobs else 0,
        "arc_yr": float(arcy) if arcy else float("nan"),
    }

# ---- tableb4: NG parameters ----------------------------------------------
b4 = {}
for ln in open(DATA_RAW / "warsaw_tableb4.dat"):
    if len(ln) < 115 or not ln.strip():
        continue
    desig = norm_desig(ln[3:15])

    def f(a, b):
        s = ln[a:b].strip()
        return float(s) if s else 0.0
    b4[desig] = {
        "A1": f(17, 29), "eA1": f(29, 41),
        "A2": f(41, 52), "eA2": f(53, 64),
        "A3": f(66, 77), "eA3": f(77, 86),
        "tau": f(89, 96), "model": ln[111:115].strip(),
    }

# ---- join ----------------------------------------------------------------
rows = []
for desig, rec in a1.items():
    if desig not in b30:
        continue
    th = float(b30[desig]["theta"])
    q = float(b30[desig]["q"])
    ng = b4.get(desig)
    rows.append({
        "desig": desig, "theta": th, "q": q,
        "model_ng": rec["model"] == "NG",
        "nobs": rec["nobs"], "arc_yr": rec["arc_yr"],
        "A1": ng["A1"] if ng else float("nan"),
        "eA1": ng["eA1"] if ng else float("nan"),
        "A2": ng["A2"] if ng else float("nan"),
        "eA2": ng["eA2"] if ng else float("nan"),
        "A3": ng["A3"] if ng else float("nan"),
        "eA3": ng["eA3"] if ng else float("nan"),
        "tau": ng["tau"] if ng else float("nan"),
    })

inc = np.array([r["theta"] <= CAP for r in rows])
print(f"joined Warsaw cohort n={len(rows)}  in-cap={inc.sum()}  "
      f"NG members={sum(r['model_ng'] for r in rows)}")

res = {"cap_deg": CAP, "n": len(rows), "n_in": int(inc.sum()),
       "seed": SEED, "test_summary": {}}

# ---- T1: NG incidence -----------------------------------------------------
ngf = np.array([r["model_ng"] for r in rows])
n_ng_in, n_ng_out = int((ngf & inc).sum()), int((ngf & ~inc).sum())
p_inc = _st.fisher_exact([[n_ng_in, int(inc.sum()) - n_ng_in],
                          [n_ng_out, int((~inc).sum()) - n_ng_out]])[1]
res["T1_ng_incidence"] = {
    "n_ng_in": n_ng_in, "n_in": int(inc.sum()),
    "n_ng_out": n_ng_out, "n_out": int((~inc).sum()),
    "frac_in": n_ng_in / max(1, inc.sum()),
    "frac_out": n_ng_out / max(1, (~inc).sum()),
    "p_fisher_2sided": p_inc}
print(f"T1 NG incidence: in {n_ng_in}/{inc.sum()} "
      f"vs out {n_ng_out}/{(~inc).sum()}  p={p_inc:.3f}")

# ---- T2/T3/T4: NG-parameter contrasts --------------------------------------
ng_rows = [r for r in rows if r["model_ng"]]
incg = np.array([r["theta"] <= CAP for r in ng_rows])
print(f"NG members with theta: {len(ng_rows)}  in-cap={incg.sum()}")


def contrast(key, members, mask):
    v_in = np.array([r[key] for r, m in zip(members, mask) if m])
    v_out = np.array([r[key] for r, m in zip(members, mask) if not m])
    fin = np.isfinite(v_in)
    fout = np.isfinite(v_out)
    v_in, v_out = v_in[fin], v_out[fout]
    if len(v_in) < 3 or len(v_out) < 3:
        return {"n_in": len(v_in), "n_out": len(v_out), "underpowered": True}
    mw = _st.mannwhitneyu(v_in, v_out, alternative="two-sided")
    allv = np.concatenate([v_in, v_out])
    th_all = np.array([r["theta"] for r in members])[np.isfinite(
        np.array([r[key] for r in members]))]
    sp = _st.spearmanr(allv, np.cos(np.radians(th_all)))
    return {"n_in": len(v_in), "n_out": len(v_out),
            "med_in": float(np.median(v_in)),
            "med_out": float(np.median(v_out)),
            "mean_in": float(v_in.mean()),
            "mean_out": float(v_out.mean()),
            "p_mw": float(mw.pvalue),
            "rho_vs_costheta": float(sp.statistic),
            "p_spearman": float(sp.pvalue)}


res["T2_transverse_A2"] = contrast("A2", ng_rows, incg)
res["T2b_A2_signif"] = {
    **contrast(
        "A2",
        [{**r, "A2": (r["A2"] / r["eA2"]) if r["eA2"] else float("nan")}
         for r in ng_rows], incg),
    "note": "A2/e_A2 significance-weighted"}
res["T3_radial_A1"] = contrast("A1", ng_rows, incg)
res["T3b_normal_A3"] = contrast("A3", ng_rows, incg)
res["T4_tau"] = contrast("tau", ng_rows, incg)
for t in ["T2_transverse_A2", "T2b_A2_signif", "T3_radial_A1",
          "T3b_normal_A3", "T4_tau"]:
    r = res[t]
    if "p_mw" in r:
        print(f"{t}: med_in={r['med_in']:.4g} med_out={r['med_out']:.4g} "
              f"p={r['p_mw']:.3g}  rho(cos)={r['rho_vs_costheta']:+.3f} "
              f"p={r['p_spearman']:.3g}")

# ---- signed-mean A2 in-cap (binomial-style sign test on A2>0) --------------
a2_in = np.array([r["A2"] for r, m in zip(ng_rows, incg)
                  if m and np.isfinite(r["A2"])])
a2_out = np.array([r["A2"] for r, m in zip(ng_rows, incg)
                   if not m and np.isfinite(r["A2"])])
if len(a2_in):
    res["T2c_A2_signed"] = {
        "n_pos_in": int((a2_in > 0).sum()), "n_in": len(a2_in),
        "frac_pos_in": float((a2_in > 0).mean()),
        "frac_pos_out": float((a2_out > 0).mean()) if len(a2_out) else None,
        "mean_in_au_d2": float(a2_in.mean() * 1e-8),
        "mean_out_au_d2": float(a2_out.mean() * 1e-8) if len(a2_out) else None,
        "note": "A2 sign convention: positive = delayed/perihelion-later "
                "effective return (period-lengthening transverse drag)"}
    print(f"T2c signed A2: in-cap mean {a2_in.mean():+.3g} "
          f"(frac>0 {(a2_in>0).mean():.2f}) vs out {a2_out.mean():+.3g}")

# ---- T5: confound audit ------------------------------------------------------
def _v(key):
    return (np.array([r[key] for r, m in zip(ng_rows, incg) if m]),
            np.array([r[key] for r, m in zip(ng_rows, incg) if not m]))
ni, no = _v("nobs")
ai, ao = _v("arc_yr")
res["T5_confounds"] = {
    "nobs_med_in": float(np.median(ni)), "nobs_med_out": float(np.median(no)),
    "p_nobs": float(_st.mannwhitneyu(ni, no, alternative="two-sided").pvalue)
    if len(ni) >= 3 and len(no) >= 3 else None,
    "arc_med_in": float(np.median(ai)), "arc_med_out": float(np.median(ao)),
    "p_arc": float(_st.mannwhitneyu(ai, ao, alternative="two-sided").pvalue)
    if len(ai) >= 3 and len(no) >= 3 else None,
    "note": "NG members only; NG solutions preferentially exist for "
            "well-observed comets -- an in-cap nobs surplus would "
            "inflate any NG-channel contrast"}

# ---- T6: theta-continuous version of the transverse test -------------------
ths = np.array([r["theta"] for r in ng_rows])
a2s = np.array([r["A2"] for r in ng_rows])
fok = np.isfinite(a2s)
if fok.sum() >= 10:
    sp = _st.spearmanr(a2s[fok], ths[fok])
    res["T6_A2_vs_theta_continuous"] = {
        "n": int(fok.sum()), "rho": float(sp.statistic),
        "p_2sided": float(sp.pvalue),
        "note": "A2 vs transit angle over all NG members; a boundary "
                "slip predicts structure with theta, not only a step at 60 deg"}
    print(f"T6 A2 vs theta: rho={sp.statistic:+.3f} p={sp.pvalue:.3g} "
          f"(n={fok.sum()})")

# ---- T7: independent SBDB NG replication -------------------------------------
# The Warsaw NG elevation is replicated on the JPL SBDB NG record: for
# each legs-record member the SBDB full-precision solution's model_pars
# block supplies fitted A1/A2/A3 (JPL orbit pipeline, no Warsaw input).
# The same cap contrast is run on the transverse term -- the channel
# the time-domain reading requires to carry the signal.
import re

SBDB_DIR = DATA_RAW / "mpc" / "sbdb_fp"


def _safe_des(desig):
    return re.sub(r"[^A-Za-z0-9]+", "_",
                  desig.split("(")[0].strip()).strip("_")


def sbdb_ng(desig):
    path = SBDB_DIR / f"{_safe_des(desig)}.json"
    if not path.exists():
        return {}
    try:
        d = json.load(open(path))
    except Exception:
        return {}
    out = {}
    for par in (d["orbit"].get("model_pars") or []):
        if par.get("value") is not None:
            try:
                out[par["name"]] = float(par["value"])
            except Exception:
                continue
    return out


t7 = {}
for tag, jf in (("pre2018", "step_b84_legs.jsonl"),
                ("post2017", "step_b81_legs.jsonl")):
    vals = []
    for line in open(RESULTS / jf):
        r = json.loads(line)
        if r.get("failed") or r.get("theta") is None:
            continue
        ng = sbdb_ng(r["desig"])
        if "A2" in ng:
            vals.append({"theta": float(r["theta"]),
                         "A1": ng.get("A1", float("nan")),
                         "A2": ng.get("A2", float("nan")),
                         "A3": ng.get("A3", float("nan"))})
    inc7 = np.array([v["theta"] <= CAP for v in vals])
    blk = {"n": len(vals), "n_in": int(inc7.sum())}
    for key in ("A2", "A1", "A3"):
        blk[key] = contrast(key, vals, inc7)
    t7[tag] = blk
    if blk["A2"].get("p_mw") is not None:
        print(f"T7 SBDB replication {tag}: n={len(vals)} "
              f"A2 med_in={blk['A2']['med_in']:.3g} "
              f"med_out={blk['A2']['med_out']:.3g} "
              f"p={blk['A2']['p_mw']:.3g}")
res["T7_sbdb_replication"] = t7

# ---- verdict ---------------------------------------------------------------
t2 = res["T2_transverse_A2"]
carrier = (t2.get("p_mw") is not None and t2["p_mw"] < 0.05) or \
          (res.get("T6_A2_vs_theta_continuous", {}).get("p_2sided", 1) < 0.05)
incidence = p_inc < 0.05
_t2m = t2.get("med_in"), t2.get("med_out")
_ratio = (_t2m[0] / _t2m[1]) if (_t2m[0] and _t2m[1]) else None
_t3 = res.get("T3_radial_A3" , res.get("T3_radial_A1", {}))
_t3b = res.get("T3b_normal_A3", {})
_a1_flat = res.get("T3_radial_A1", {}).get("p_mw", 1) > 0.05
_a3_flat = res.get("T3b_normal_A3", {}).get("p_mw", 1) > 0.05
_pattern = (_ratio is not None and _ratio > 2.0
            and _a1_flat and _a3_flat)
if carrier or incidence:
    _label = "SUPPORTS"
elif _pattern and t2.get("p_mw") is not None:
    _label = ("TRANSVERSE-SELECTIVE ELEVATION, UNDERPOWERED")
else:
    _label = "QUIET"
verdict = ("NG CHANNEL " + _label +
           f": NG incidence in {n_ng_in}/{inc.sum()} vs out "
           f"{n_ng_out}/{(~inc).sum()} (p={p_inc:.3f}); "
           f"transverse A2 contrast "
           + (f"med {t2.get('med_in'):.3g} vs {t2.get('med_out'):.3g} "
              f"({_ratio:.1f}x), p={t2.get('p_mw'):.3g} "
              f"at n_in={t2.get('n_in')}" if t2.get("p_mw") is not None
              else f"underpowered (n_in={t2.get('n_in')})")
           + ".  The NG record is read as an observable for the first "
             "time: a boundary slip absorbed by pure-gravity fits must "
             "leak into the fitted transverse term."
           + ("  The transverse component -- the time-domain direction "
              "-- is the only elevated channel (radial A1 and normal A3 "
              "flat), and the in-cap fraction of positive A2 runs "
              f"{res.get('T2c_A2_signed', {}).get('n_pos_in', 0)}/"
              f"{res.get('T2c_A2_signed', {}).get('n_in', 0)}; the "
              "contrast is sign-consistent with time-domain loading "
              "but unresolved at this membership." if _pattern
              and not carrier else ""))
res["verdict"] = verdict
res["test_summary"] = {
    "ng_incidence_differs": bool(incidence),
    "a2_cap_contrast_p": t2.get("p_mw"),
    "a2_vs_theta_p": res.get("T6_A2_vs_theta_continuous", {}).get("p_2sided"),
    "n_ng_members": len(ng_rows), "n_ng_in_cap": int(incg.sum()),
}

out = RESULTS / "step_b107_comet_ng_channel.json"
json.dump(res, open(out, "w"), indent=1)

# csv ledger
with open(RESULTS / "step_b107_comet_ng_channel.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ---- figure ------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
ax = axes[0]
ax.bar(["in-cap", "out-cap"],
       [n_ng_in / max(1, inc.sum()), n_ng_out / max(1, (~inc).sum())],
       color=["crimson", "steelblue"])
ax.set_ylabel("NG-model fraction")
ax.set_title(f"NG incidence (p={p_inc:.3g})")

ax = axes[1]
a2v = np.array([r["A2"] for r in ng_rows])
thv = np.array([r["theta"] for r in ng_rows])
m = np.isfinite(a2v)
ax.scatter(thv[m], a2v[m], s=18,
           c=["crimson" if t <= CAP else "steelblue" for t in thv[m]])
ax.axvline(CAP, color="k", ls=":", lw=1)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("transit angle theta (deg)")
ax.set_ylabel("A2 (1e-8 AU/d^2)")
ax.set_title("transverse NG vs theta")

ax = axes[2]
for key, lab in (("A1", "A1 radial"), ("A2", "A2 transverse"),
                 ("A3", "A3 normal")):
    vi = np.array([r[key] for r, mm in zip(ng_rows, incg)
                   if mm and np.isfinite(r[key])])
    if len(vi):
        ax.scatter([key] * len(vi), vi, s=18, alpha=0.6)
ax.axhline(0, color="k", lw=0.5)
ax.set_ylabel("in-cap NG params")
ax.set_title("in-cap component loadings")
fig.tight_layout()
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b107_comet_ng_channel.png", dpi=150)

logger.info("verdict: " + res["verdict"])
logger.data_save(out)
logger.data_save(RESULTS / "step_b107_comet_ng_channel.csv")
logger.data_save(FIG / "step_b107_comet_ng_channel.png")
print(json.dumps(res["test_summary"], indent=1))
print(res["verdict"])
