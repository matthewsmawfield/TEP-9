#!/usr/bin/env python3
"""step_120: Era x lineage factorial -- the pre-2018 SBDB LPC cohort.

The prospective programme has produced a 2x2 design with one cell
missing:

  CODE solutions x pre-2018 objects   -> declared-axis anomaly
                                         (steps 030-106)
  SBDB solutions x post-2017 objects  -> displaced unipolar dipole,
                                         cohort-private (117-119)
  SBDB solutions x pre-2018 objects   -> ??? (this step)

The reading of the prospective failure hinges on the missing cell:
  - if the pre-2018 SBDB cohort recovers the DECLARED axis, the
    displaced dipole is carried by the modern objects' data (short
    arcs, modern discovery geometry), and the CODE anomaly stands
    unchallenged on the historical record;
  - if it recovers the post-2017 DISPLACED dipole instead, the
    structure is a property of the SBDB fit record irrespective of
    object era -- a far stronger systematic, and one that turns the
    question back onto the CODE anomaly itself.

A second discriminator is built into the SBDB record: every
pre-2018 comet carries a last_obs date, so the cohort splits by FIT
VINTAGE -- solutions whose arcs end before 2018 (historical record)
versus solutions updated on modern data.  If the displaced structure
tracks fit vintage rather than object era it appears in the
modern-refit members of the pre-2018 cohort.

Cohort: identical cuts to step_117 -- C/ designations, primary
objects, 0.95 < e < 1.5, arc >= 30 d, nobs >= 20 -- except
designation year < 2018.  The instrument is the shared
scripts/utils/lpc_boundary.py bidirectional integrator (identical to
step_063/064/065/117): REBOUND/IAS15, DE440s, Sun plus planetary
barycentres, 250 AU barycentric boundary.

Tests
-----
  T1  declared-axis instruments on the deep-plunger primary (q<3.1):
      cap contrast, residual cap test, Spearman vs theta -- does
      SBDB pre-2018 reproduce the CODE anomaly?
  T2  free axis scan (10-deg grid, MWU on drot and on the residual
      field) -- does the cohort's own axis land near the declared
      axis or the post-2017 displaced direction?
  T3  displaced-axis gap: cap-gap at the post-2017 recovered axes
      (120,-40) residual / (160,-10) raw -- is the displaced dipole
      present in the pre-2018 record?
  T4  fit-vintage split: last_obs >= 2018 vs < 2018 -- does the
      displaced structure track fit vintage or object era?
  T5  longitude harmonic of the residual field -- phase against the
      post-2017 (182 deg) and CODE (26 deg) references.
  T6  mechanism-split replication: bound vs hyperbolic,
      perihelion-spanning vs one-sided arcs (the step_119
      localization) on the pre-2018 cohort.
  T7  CODE-overlap subset: members also in the CODE class-1 table --
      SBDB-solution behaviour on bodies whose CODE solutions carry
      the anomaly.
  T8  quality-stratified declared-axis test (the step_118 strata) --
      does the in-cap residual concentrate in the best-determined
      members, as in CODE?

Inputs
------
data/raw/sbdb/sbdb_comets_all.json
data/raw/spice/de440s.bsp
data/raw/code/code_osculating.html          (overlap membership)

Outputs
-------
results/step_b84_legs.jsonl               (integration checkpoint)
results/step_b84_pre2018_sbdb.json
results/step_b84_pre2018_sbdb.csv
results/figures/step_b84_pre2018_sbdb.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (
    DATA_RAW, RESULTS, sep, lv, lb, parse_code)
from scripts.utils.lpc_boundary import (
    R_STOP, TNO, process_comet, worker_init as _lpc_worker_init)
logger = StepLogger("step_120_pre2018_sbdb_factorial")

import csv
import datetime as _dt
import json
import math
import re

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

import spiceypy as sp

logger.header("Era x lineage factorial: pre-2018 SBDB cohort")

SEED = 20260919
rng = np.random.default_rng(SEED)
N_PERM = 20000
CAP = 60.0
U_RES = lv(120.0, -40.0)       # post-2017 residual-field recovered axis
U_RAW = lv(160.0, -10.0)       # post-2017 raw-drot recovered axis

STATS_ONLY = "--stats-only" in _sys.argv

WORKERS = 1
for _i, _a in enumerate(_sys.argv):
    if _a == "--workers" and _i + 1 < len(_sys.argv):
        WORKERS = max(1, int(_sys.argv[_i + 1]))

SPK = DATA_RAW / "spice" / "de440s.bsp"
sp.furnsh(str(SPK))


def _worker_init():
    _lpc_worker_init(SPK)


def _process_comet(rec):
    return process_comet(rec, axis=TNO)


# ------------------------------------------------------------------
# 1. Cohort selection (identical cuts to step_117, year < 2018)
# ------------------------------------------------------------------

sbdb = json.load(open(DATA_RAW / "sbdb" / "sbdb_comets_all.json"))
fields = sbdb["fields"]
rows = [dict(zip(fields, rec)) for rec in sbdb["data"]]


def fnum(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return float("nan")


def desig_year(name):
    m = re.match(r"\s*[CP]/(\d{4})", name)
    return int(m.group(1)) if m else None


def is_fragment(name):
    return bool(re.match(r"\s*C/\d{4}\s+\S+-\w", name))


def _jd(datestr):
    d = _dt.date.fromisoformat(str(datestr)[:10])
    return d.toordinal() + 1721424.5


cohort_all = []
excluded = {"year": 0, "frag": 0, "e": 0, "q": 0, "arc": 0, "iso": 0}
for r in rows:
    name = str(r["full_name"]).strip()
    if not name.startswith("C/"):
        continue
    yr = desig_year(name)
    if yr is None or yr >= 2018:
        excluded["year"] += 1
        continue
    if is_fragment(name):
        excluded["frag"] += 1
        continue
    e = fnum(r, "e")
    if not np.isfinite(e):
        excluded["e"] += 1
        continue
    if e >= 1.5:
        excluded["iso"] += 1
        continue
    q = fnum(r, "q")
    arc = fnum(r, "data_arc")
    nobs = fnum(r, "n_obs_used")
    if e < 0.95:
        continue
    if not (np.isfinite(q) and 0.1 <= q):
        excluded["q"] += 1
        continue
    if not (np.isfinite(arc) and arc >= 30 and np.isfinite(nobs)
            and nobs >= 20):
        excluded["arc"] += 1
        continue
    rec = dict(name=name, yr=yr, q=q, e=e, i=fnum(r, "i"),
               om=fnum(r, "om"), w=fnum(r, "w"), a=fnum(r, "a"),
               epoch=fnum(r, "epoch"), tp=fnum(r, "tp"),
               cc=r.get("condition_code"), arc=arc, nobs=nobs,
               cls=r.get("class"),
               ng=any(np.isfinite(fnum(r, k)) and fnum(r, k) != 0.0
                      for k in ("A1", "A2", "A3", "DT")))
    if not (np.isfinite(rec["epoch"]) and np.isfinite(rec["tp"])):
        excluded["e"] += 1
        continue
    try:
        rec["last_jd"] = _jd(r["last_obs"])
        rec["first_jd"] = _jd(r["first_obs"])
    except (TypeError, ValueError):
        rec["last_jd"] = np.nan
    rec["modern_fit"] = (np.isfinite(rec.get("last_jd", np.nan))
                         and rec["last_jd"] >= _jd("2018-01-01"))
    cohort_all.append(rec)

cohort = [r for r in cohort_all if r["q"] < 3.1]
logger.info(f"pre-2018 comets: near-parabolic usable {len(cohort_all)}; "
            f"deep-plunger primary (q<3.1) {len(cohort)}; "
            f"modern-fit (last_obs>=2018) "
            f"{sum(r['modern_fit'] for r in cohort_all)}; "
            f"excluded {excluded}")

# CODE overlap set
code_osc = parse_code(str(DATA_RAW / "code" / "code_osculating.html"))


def _desig(nm):
    m = re.match(r"\s*([CP]/\d{4}\s+\w+(?:-\w+)?)", str(nm))
    return m.group(1) if m else None


CODE_KEYS = {re.sub(r"\s+", " ", k).strip() for k in code_osc}
for r in cohort_all:
    k = _desig(r["name"])
    r["in_code"] = k is not None and re.sub(
        r"\s+", " ", k).strip() in CODE_KEYS

# ------------------------------------------------------------------
# 2. Integration with checkpointing (identical to step_117)
# ------------------------------------------------------------------

CKPT = RESULTS / "step_b84_legs.jsonl"


def _load_checkpoint():
    done, failed_names, hard = {}, set(), set()
    if CKPT.exists():
        for line in CKPT.read_text().splitlines():
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if row.get("failed"):
                failed_names.add(row["desig"])
                if not row.get("skipped"):
                    hard.add(row["desig"])
            else:
                done[row["desig"]] = row
    return done, failed_names, hard


def run_cohort(recs, tag, aph_limit=None):
    done, failed_names, hard = _load_checkpoint()
    out_rows, failed = [], []
    todo = []
    for rec in recs:
        if rec["name"] in done:
            out_rows.append(done[rec["name"]])
        elif rec["name"] in failed_names:
            failed.append(rec["name"])
        else:
            todo.append(rec)
    if not out_rows and not todo and any(n in hard for n in failed):
        names = {r["name"] for r in recs}
        keep, quarantined = [], []
        for line in CKPT.read_text().splitlines():
            try:
                desig = json.loads(line).get("desig")
            except json.JSONDecodeError:
                keep.append(line)
                continue
            (quarantined if desig in names else keep).append(line)
        suspect = CKPT.with_suffix(f".{tag}.suspect")
        suspect.write_text("\n".join(quarantined) + "\n")
        CKPT.write_text("\n".join(keep) + "\n")
        logger.warning(f"{tag}: checkpoint holds {len(failed)} failures "
                       f"and 0 successes -- quarantined to "
                       f"{suspect.name}; re-attempting")
        todo, failed = list(recs), []
    if todo:
        logger.info(f"{tag}: {len(todo)} comets to integrate "
                    f"({len(done) + len(failed)} checkpointed)")
    n_done = 0

    def _checkpoint(row, name, err=None, skipped=False):
        with open(CKPT, "a") as ck:
            if row is None:
                ck.write(json.dumps({"desig": name, "failed": True,
                                     "skipped": skipped,
                                     "err": str(err)[:200]})
                         + "\n")
            else:
                ck.write(json.dumps(row) + "\n")

    def _account(rec, row, err, skipped=False):
        nonlocal n_done
        if row is None:
            failed.append(rec["name"])
            logger.warning(f"{rec['name']}: {err}")
        else:
            out_rows.append(row)
        _checkpoint(row, rec["name"], err, skipped)
        n_done += 1
        if n_done % 25 == 0:
            logger.info(f"{tag}: {n_done}/{len(todo)} integrated")

    def _preclassify(rec):
        if aph_limit is not None and rec["e"] < 1:
            aph = rec["q"] * (1 + rec["e"]) / (1 - rec["e"])
            if aph < aph_limit:
                return (f"osculating aphelion {aph:.0f} AU < "
                        f"{aph_limit:.0f} AU -- boundary unreachable")
        return None

    if todo and WORKERS > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        todo_int = [r for r in todo if _preclassify(r) is None]
        for rec in todo:
            err = _preclassify(rec)
            if err:
                _account(rec, None, err, skipped=True)
        with ctx.Pool(WORKERS, initializer=_worker_init) as pool:
            for rec, (row, err) in zip(
                    todo_int, pool.imap(_process_comet, todo_int)):
                _account(rec, row, err)
    else:
        for rec in todo:
            err = _preclassify(rec)
            if err:
                _account(rec, None, err, skipped=True)
                continue
            row, err = _process_comet(rec)
            _account(rec, row, err)
    logger.info(f"{tag}: integrated {len(out_rows)} comets x2 legs; "
                f"{len(failed)} failed")
    return out_rows, failed


if STATS_ONLY:
    prows = []
    with open(RESULTS / "step_b84_pre2018_sbdb.csv") as _f:
        for r in csv.DictReader(_f):
            for k in ("yr", "q", "i", "e", "arc", "nobs", "theta", "drot",
                      "d_in", "d_out", "daa", "denc", "aa_back", "aa_fwd",
                      "t_back", "t_fwd", "dtau", "dtau_in", "dtau_out"):
                r[k] = float(r[k])
            r["ng"] = r.get("ng") in ("True", True)
            r["aph_lb"] = [float(x) for x in r["aph_lb"].split(";")]
            r["modern_fit"] = r.get("modern_fit") in ("True", True)
            r["in_code"] = r.get("in_code") in ("True", True)
            r.pop("dtau_unexplained", None)
            prows.append(r)
    logger.info(f"stats-only mode: {len(prows)} rows loaded from CSV")
else:
    logger.info("integrating pre-2018 cohort ...")
    prows, pfailed = run_cohort(cohort_all, "pre2018",
                              aph_limit=0.9 * R_STOP)
    if not prows:
        raise RuntimeError(
            "pre-2018 cohort produced 0 integrated comets -- refusing "
            "to overwrite results with an empty analysis")
    _mf = {r["name"]: r["modern_fit"] for r in cohort_all}
    _ic = {r["name"]: r["in_code"] for r in cohort_all}
    for r in prows:
        r["modern_fit"] = _mf.get(r["desig"], False)
        r["in_code"] = _ic.get(r["desig"], False)

# ------------------------------------------------------------------
# 3. Instruments (identical to steps 117/118)
# ------------------------------------------------------------------

def resid_model(sub):
    v = np.array([r["drot"] for r in sub])
    K = np.abs(np.array([r["daa"] for r in sub]))
    D = np.array([r["denc"] for r in sub])
    Q = np.array([r["q"] for r in sub])
    I = np.array([r["i"] for r in sub])
    y = np.log10(v)
    X = np.column_stack([np.ones(len(y)), np.log10(K + 1.0),
                         np.log10(D), Q, I])
    with np.errstate(all="ignore"):
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ coef
    return resid


def aph_u(r):
    return np.asarray(r["aph_lb"] if isinstance(r["aph_lb"], list)
                    else [float(x) for x in r["aph_lb"].split(";")])


def gap_at(sub, resid, u):
    th = np.array([sep(aph_u(r), u) for r in sub])
    inc = th < CAP
    if inc.sum() < 4 or (~inc).sum() < 4:
        return dict(n_in=int(inc.sum()), status="insufficient_coverage")
    return dict(
        n_in=int(inc.sum()),
        med_in=float(np.median(resid[inc])),
        med_out=float(np.median(resid[~inc])),
        gap=float(np.median(resid[inc]) - np.median(resid[~inc])),
        mwu_greater_p=float(mannwhitneyu(
            resid[inc], resid[~inc], alternative="greater").pvalue),
        rho=float(spearmanr(th, resid).statistic),
        p_2sided=float(spearmanr(th, resid).pvalue))




def _fg(d):
    g = d.get("gap")
    return f"{g:+.3f}" if g is not None else "n/a"

def run_declared(sub):
    """The step_117 channel set at the pre-declared axis."""
    v = np.array([r["drot"] for r in sub])
    th = np.array([r["theta"] for r in sub])
    inc = th < CAP
    if inc.sum() < 4 or (~inc).sum() < 4 or len(sub) < 10:
        return dict(n=len(sub), n_in=int(inc.sum()),
                    status="insufficient_sample")
    resid = resid_model(sub)
    out = dict(n=len(sub), n_in=int(inc.sum()),
               med_drot_in=float(np.median(v[inc])),
               med_drot_out=float(np.median(v[~inc])),
               drot_in_gt_out_p=float(mannwhitneyu(
                   v[inc], v[~inc], alternative="greater").pvalue),
               rho_vs_theta=float(spearmanr(th, v).statistic),
               p_2sided=float(spearmanr(th, v).pvalue))
    rinc = th < CAP
    out["resid_cap"] = dict(
        med_in=float(np.median(resid[rinc])),
        med_out=float(np.median(resid[~rinc])),
        p=float(mannwhitneyu(resid[rinc], resid[~rinc],
                             alternative="greater").pvalue))
    return out


GRID = [(l, b) for l in np.arange(0, 360, 10)
        for b in np.arange(-60, 61, 10)]
GRID_U = [lv(l, b) for l, b in GRID]


def free_scan(sub, field):
    """Min-p MWU axis over the 10-deg grid (the step_118 instrument)."""
    A = np.stack([aph_u(r) for r in sub])
    G = np.stack(GRID_U)
    dots = np.clip(A @ G.T, -1.0, 1.0)
    TH = np.degrees(np.arccos(dots))
    best = {"p": 2.0}
    n_better_declared = 0
    p_decl = None
    for j, (l, b) in enumerate(GRID):
        inc = TH[:, j] < CAP
        if inc.sum() < 4 or (~inc).sum() < 4:
            continue
        p = mannwhitneyu(field[inc], field[~inc],
                         alternative="greater").pvalue
        if abs(l - 34) < 1 and abs(b + 13) < 1:
            p_decl = p
        if p < best["p"]:
            best = {"p": float(p), "lam": float(l), "beta": float(b),
                    "n_in": int(inc.sum())}
    u_best = lv(best["lam"], best["beta"])
    best["sep_from_declared"] = sep(u_best, TNO)
    best["sep_from_post2017_resid"] = sep(u_best, U_RES)
    return best, p_decl


def lon_harmonic(sub, resid, order=1):
    L = np.radians([lb(aph_u(r))[0] for r in sub])
    c = np.sum(resid * np.exp(-1j * order * L)) / len(resid)
    amp = 2 * abs(c)
    # c = (A/2) e^{-i phi} for a field ~ A cos(order*L - phi), so the
    # phase of the field MAXIMUM is angle(conj(c)) = -angle(c).
    # angle(-conj(c)) would return the antipode (the field minimum).
    ph = float(np.degrees(np.angle(np.conj(c)) / order) % 360)
    cnt = 0
    for _ in range(2000):
        cp = np.sum(rng.permutation(resid)
                    * np.exp(-1j * order * L)) / len(resid)
        cnt += 2 * abs(cp) >= amp
    return dict(amplitude=float(amp), phase_deg=ph,
                p_perm=float((cnt + 1) / 2001))


# ------------------------------------------------------------------
# 4. Analysis
# ------------------------------------------------------------------

for r in prows:
    r["aph"] = aph_u(r)

resid_all = resid_model(prows)
drot_all = np.array([r["drot"] for r in prows])

primary = [r for r in prows if r["q"] < 3.1]
prim_idx = [i for i, r in enumerate(prows) if r["q"] < 3.1]
resid_p = resid_all[prim_idx]

t1 = run_declared(primary)
logger.info(f"T1 primary declared-axis: n={t1['n']} in={t1['n_in']} "
            f"med {t1.get('med_drot_in')}/{t1.get('med_drot_out')} "
            f"resid cap p={t1.get('resid_cap', {}).get('p')}")

best_raw, pdecl_raw = free_scan(primary, np.array(
    [r["drot"] for r in primary]))
best_res, pdecl_res = free_scan(primary, resid_p)
t2 = {"raw_drot": {"best": best_raw, "declared_p": pdecl_raw},
      "resid": {"best": best_res, "declared_p": pdecl_res}}
logger.info(f"T2 free scan raw: ({best_raw['lam']:.0f},"
            f"{best_raw['beta']:.0f}) p={best_raw['p']:.2e} | "
            f"resid: ({best_res['lam']:.0f},{best_res['beta']:.0f}) "
            f"p={best_res['p']:.2e}")

t3 = {
    "at_post2017_resid_axis": gap_at(primary, resid_p, U_RES),
    "at_post2017_raw_axis": gap_at(primary, resid_p, U_RAW),
    "at_declared_axis": gap_at(primary, resid_p, TNO),
}
logger.info(f"T3 displaced-axis gaps: (120,-40) "
            f"{_fg(t3['at_post2017_resid_axis'])}, "
            f"(160,-10) {_fg(t3['at_post2017_raw_axis'])}, "
            f"declared {_fg(t3['at_declared_axis'])}")

# T4 fit-vintage split
t4 = {}
for tag, sub in [("modern_fit_last_obs_ge_2018",
                  [r for r in primary if r.get("modern_fit")]),
                 ("historical_fit_last_obs_lt_2018",
                  [r for r in primary if not r.get("modern_fit")])]:
    if len(sub) < 15:
        t4[tag] = dict(n=len(sub), status="insufficient_sample")
        continue
    rs = resid_model(sub)
    t4[tag] = dict(
        n=len(sub),
        declared=gap_at(sub, rs, TNO),
        post2017_resid_axis=gap_at(sub, rs, U_RES))
    logger.info(f"T4 {tag} (n={len(sub)}): declared gap "
                f"{_fg(t4[tag]['declared'])}, "
                f"displaced gap "
                f"{_fg(t4[tag]['post2017_resid_axis'])}")

t5 = {"resid_lon_harmonic_order1": lon_harmonic(primary, resid_p),
      "note": "post-2017 reference phase 182 deg; CODE 26 deg"}
logger.info(f"T5 pre-2018 resid order-1: amp "
            f"{t5['resid_lon_harmonic_order1']['amplitude']:.3f} "
            f"phase {t5['resid_lon_harmonic_order1']['phase_deg']:.0f} "
            f"p={t5['resid_lon_harmonic_order1']['p_perm']:.4f}")

# T6 mechanism splits (step_119 localization replicated here)
JD2018 = _jd("2018-01-01")
_mj = {r["name"]: (r.get("first_jd"), r.get("last_jd"), r.get("tp"))
       for r in cohort_all}
for r in primary:
    fj, lj, tp = _mj.get(r["desig"], (None, None, None))
    try:
        r["f_pre"] = float(np.clip(
            (float(tp) - fj) / max(lj - fj, 1e-9), 0.0, 1.0))
    except (TypeError, ValueError):
        r["f_pre"] = np.nan
t6 = {}
for tag, sub in [
        ("bound_e_lt_1", [r for r in primary if r["e"] < 1.0]),
        ("hyperbolic_e_ge_1", [r for r in primary if r["e"] >= 1.0]),
        ("arc_spans_perihelion",
         [r for r in primary if 0.0 < r.get("f_pre", 0) < 1.0]),
        ("arc_one_sided",
         [r for r in primary if r.get("f_pre", 0) in (0.0, 1.0)])]:
    if len(sub) < 15:
        t6[tag] = dict(n=len(sub), status="insufficient_sample")
        continue
    rs = resid_model(sub)
    t6[tag] = dict(n=len(sub),
                   declared=gap_at(sub, rs, TNO),
                   post2017_axis=gap_at(sub, rs, U_RES))
    logger.info(f"T6 {tag} (n={len(sub)}): declared "
                f"{_fg(t6[tag]['declared'])}, displaced "
                f"{_fg(t6[tag]['post2017_axis'])}")

# T8 quality-stratified declared-axis test (mirrors the step_118
# strata): does the in-cap residual concentrate in the
# best-determined members, as in CODE?
def _gap_for(sub, u):
    if len(sub) < 15:
        return dict(n=len(sub), status="insufficient_sample")
    return dict(n=len(sub), **gap_at(sub, resid_model(sub), u))


_arcs = sorted(r["arc"] for r in primary)
_arc_q3 = _arcs[int(0.75 * len(_arcs))]
t8 = {
    "arc_gt_730_nobs_gt_500": _gap_for(
        [r for r in primary if r["arc"] > 730 and r["nobs"] > 500], TNO),
    "arc_gt_365_nobs_gt_300": _gap_for(
        [r for r in primary if r["arc"] > 365 and r["nobs"] > 300], TNO),
    "arc_top_quartile": _gap_for(
        [r for r in primary if r["arc"] >= _arc_q3], TNO),
    "nobs_gt_1000": _gap_for(
        [r for r in primary if r["nobs"] > 1000], TNO),
}
for tag, d in t8.items():
    logger.info(f"T8 {tag} (n={d.get('n')}): declared gap "
                f"{d.get('gap')}, p={d.get('mwu_greater_p')}")

# T7 CODE-overlap members
overlap = [r for r in primary if r.get("in_code")]
t7 = {"n_overlap": len(overlap)}
if len(overlap) >= 15:
    rs = resid_model(overlap)
    t7.update(declared=gap_at(overlap, rs, TNO),
              post2017_axis=gap_at(overlap, rs, U_RES))
    logger.info(f"T7 CODE-overlap (n={len(overlap)}): declared "
                f"{_fg(t7['declared'])}, displaced "
                f"{_fg(t7['post2017_axis'])}")
else:
    t7["status"] = "insufficient_sample"

# ------------------------------------------------------------------
# 5. Verdict
# ------------------------------------------------------------------

decl_gap = t3["at_declared_axis"].get("gap", 0)
disp_gap = t3["at_post2017_resid_axis"].get("gap", 0)
decl_p = t1.get("resid_cap", {}).get("p", 1.0)

if decl_gap > 0.15 and decl_p < 0.05:
    verdict = ("the pre-2018 SBDB cohort REPRODUCES the declared-axis "
               "anomaly -- the displaced dipole is carried by the "
               "modern objects' data record, not by the SBDB fitter; "
               "the CODE anomaly stands on the historical record")
elif disp_gap > 0.15 and decl_gap < 0.10:
    verdict = ("the pre-2018 SBDB cohort shows the DISPLACED dipole "
               "rather than the declared-axis anomaly -- the structure "
               "is a property of the SBDB fit record itself, and the "
               "question turns back onto the CODE anomaly's own "
               "lineage dependence")
else:
    verdict = (f"ambiguous: declared-axis gap {decl_gap:+.3f} "
               f"(resid cap p={decl_p:.3f}) vs displaced-axis gap "
               f"{disp_gap:+.3f} -- the pre-2018 SBDB record matches "
               "neither cell cleanly")

res = {
    "step": "step_120_pre2018_sbdb_factorial",
    "description": __doc__.strip().splitlines()[0],
    "inputs": ["data/raw/sbdb/sbdb_comets_all.json",
               "data/raw/spice/de440s.bsp",
               "data/raw/code/code_osculating.html"],
    "seed": SEED, "cap_deg": CAP, "axis_ecl_deg": [34, -13],
    "boundary_AU_barycentric": R_STOP,
    "cohort_counts": {"usable_near_parabolic": len(cohort_all),
                      "deep_plunger_primary": len(cohort),
                      "integrated": len(prows),
                      "modern_fit_last_obs_ge_2018": sum(
                          r["modern_fit"] for r in cohort_all),
                      "code_overlap": sum(
                          r["in_code"] for r in cohort_all),
                      "excluded": excluded},
    "T1_declared_axis": t1,
    "T2_free_scan": t2,
    "T3_displaced_axis": t3,
    "T4_fit_vintage": t4,
    "T5_longitude_harmonic": t5,
    "T6_mechanism_splits": t6,
    "T7_code_overlap": t7,
    "T8_quality_strata": t8,
    "verdict": verdict,
}

out = RESULTS / "step_b84_pre2018_sbdb.json"
json.dump(res, open(out, "w"), indent=1, default=float)

with open(RESULTS / "step_b84_pre2018_sbdb.csv", "w",
          newline="") as f:
    flat = []
    keys = []
    for r in prows:
        d = {k: (";".join(f"{x:.6f}" for x in v)
                 if isinstance(v, (list, np.ndarray)) else v)
             for k, v in r.items() if k != "aph"}
        flat.append(d)
        for k in d:
            if k not in keys:
                keys.append(k)
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    w.writerows(flat)

# ------------------------------------------------------------------
# 6. Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

ax = axes[0]
th = np.array([r["theta"] for r in primary])
rs = resid_model(primary)
ax.scatter(th, rs, s=14, c="0.55", alpha=0.5)
inc = th < CAP
ax.scatter(th[inc], rs[inc], s=16, c="crimson", alpha=0.7)
ax.axvline(CAP, color="k", ls=":", lw=1)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel(r"$\theta$ from declared axis (deg)")
ax.set_ylabel("residual log rotation")
ax.set_title("pre-2018 SBDB cohort at the declared axis", fontsize=10)

ax = axes[1]
thr = np.array([sep(r["aph"], U_RES) for r in primary])
incr = thr < CAP
ax.scatter(thr[~incr], rs[~incr], s=14, c="0.55", alpha=0.5)
ax.scatter(thr[incr], rs[incr], s=16, c="steelblue", alpha=0.7)
ax.axvline(CAP, color="k", ls=":", lw=1)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel(r"$\theta$ from post-2017 axis (120,-40) (deg)")
ax.set_ylabel("residual log rotation")
ax.set_title("pre-2018 cohort at the displaced axis", fontsize=10)

ax = axes[2]
labels, gaps = [], []
for tag, d in [("modern fit\n(last_obs>=2018)",
                t4.get("modern_fit_last_obs_ge_2018")),
               ("historical fit\n(last_obs<2018)",
                t4.get("historical_fit_last_obs_lt_2018")),
               ("bound", t6.get("bound_e_lt_1")),
               ("hyperbolic", t6.get("hyperbolic_e_ge_1")),
               ("spanning arc", t6.get("arc_spans_perihelion")),
               ("one-sided arc", t6.get("arc_one_sided")),
               ("CODE overlap", t7)]:
    if d and d.get("post2017_axis", d.get("post2017_resid_axis", {})):
        g = d.get("post2017_axis") or d.get("post2017_resid_axis") or {}
        if "gap" in g:
            labels.append(tag)
            gaps.append(g["gap"])
ax.bar(range(len(gaps)), gaps, color="steelblue")
ax.axhline(0, color="k", lw=0.5)
ax.set_xticks(range(len(gaps)))
ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
ax.set_ylabel("gap at displaced axis (dex)")
ax.set_title("displaced dipole by subset", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b84_pre2018_sbdb.png", dpi=150)
print("wrote", out)
