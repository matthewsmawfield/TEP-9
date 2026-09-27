"""Step 159 -- Perihelion-time leg-slip channel on the independent
refit record.

The transit anomaly has so far been measured in angle: the inbound-
and outbound-leg boundary asymptotes disagree (ddirf / d_in_leg), and
the catalogue leg rotation converts to an effective proper-time slip
through om_b = h/R^2 (step 065).  A conformal-lapse crossing is, at
root, a *timing* event -- the leg disagreement should therefore also
express directly as a perihelion-epoch discrepancy between orbit
solutions fitted independently to the pre- and post-perihelion arcs.

This step extracts that channel.  Each comet's two independent leg
fits (step_127 machinery on raw MPC astrometry, REBOUND/IAS15 +
DE440s, MAD clip) are re-derived on the cached observation record --
identical inputs, so the fits are deterministic reproductions -- and
each fitted leg state is converted to its osculating perihelion epoch
through a Keplerian two-body reduction anchored at the leg mid-epoch:

    tp_leg = et_leg - M0 / n ,  n = sqrt(mu / |a|^3)

The boundary products (theta, ddirf, d_in_leg, cross-arc rms) are
joined from the step_b91 checkpoint rather than recomputed, since the
timing channel needs only the leg states.  The catalogue-side
companions (dtau, drot) come from step_b84.

Registered tests

  T1  validity: per-comet |dtp_io| (inbound minus outbound osculating
      perihelion epoch) must correlate with the catalogue slip |dtau|
      and with our own leg disagreement ddirf -- both are projections
      of the same non-closure, so concordance certifies the channel.
  T2  cap contrast of |dtp_io| at the declared 60-deg cap on the
      independent theta, MWU + label-shuffle permutation, reported on
      the full dual-leg set and the CODE-overlap seat (the arena the
      registered anomaly was detected on, step_123).
  T3  signed dtp_io polarity: median sign inside the declared cap and
      inside the mirror cap -- the bipolar slip field predicts
      opposite polarities on the two lobes (step_089).
  T4  timing-vs-geometry dissociation: |dtp_io| residualized against
      (arc length, n_obs, leg-span asymmetry) -- a pure
      astrometric-noise floor should flatten the cap contrast.

Outputs
  results/step_b125_tp_slip.json
  results/step_b125_tp_slip.csv
  results/step_b125_tp_slip.jsonl        (per-comet checkpoint)
"""

import csv
import json
import math
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

import matplotlib
matplotlib.use("Agg")
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.tep9_common import (DATA_RAW, RESULTS, sep,
                                       load_jsonl_dedup)

logger = StepLogger("step_159_tp_slip")
tee_stdout(logger)

SPK = DATA_RAW / "spice" / "de440s.bsp"
LSK = DATA_RAW / "naif" / "naif0012.tls"
MPC_DIR = DATA_RAW / "mpc"
OBS_DIR = MPC_DIR / "obs"
SBDB_DIR = MPC_DIR / "sbdb_fp"
OBSC_PATH = MPC_DIR / "obscodes.json"
PROV_PATH = MPC_DIR / "provenance.json"
CKPT = RESULTS / "step_b125_tp_slip.jsonl"

CAP_DEG = 60.0
N_PERM = 5000
SEED = 20260920

from scripts.utils.parallel import cli_workers as _cli_workers
WORKERS = _cli_workers(_sys.argv)
REDO = "--redo" in _sys.argv

logger.header("Perihelion-time leg-slip channel on the independent refit record")

from scripts.utils import mpc_refit as R
from scripts.utils.lpc_boundary import TNO
import spiceypy as sp

R.configure(OBS_DIR, SBDB_DIR, OBSC_PATH, PROV_PATH, SPK, LSK,
            "step_159_tp_slip")


def _worker_init():
    R.refit_worker_init()


def osc_tp_jd(r, v, et):
    """Osculating perihelion epoch (JD TDB) of a heliocentric state.

    r in AU, v in AU/yr, et in s past J2000.  Returns (tp_jd, ecc) or
    (None, ecc) when the orbit is too close to parabolic for the
    mean-motion inversion to be stable."""
    el = sp.oscelt(np.r_[r * R.AU_KM, v * R.AU_KM / (R.DAY_YR * 86400.0)],
                   et, R.GM_SUN)
    rp, ecc, m0 = el[0], el[1], el[5]
    if abs(1.0 - ecc) < 1e-5:
        return None, ecc
    a = rp / (1.0 - ecc)                       # km; a<0 hyperbolic
    n = math.sqrt(R.GM_SUN / abs(a) ** 3)      # rad/s
    tp_et = et - m0 / n
    return tp_et / 86400.0 + 2451545.0, ecc


def _process(rec):
    """Light re-fit of one checkpointed comet: leg states only, no
    boundary propagation.  Returns (row_dict, None) or (None, err)."""
    des = rec["des"]
    try:
        raw = R.get_obs(des)
        obs = R.parse_obs(raw)
        sb = R.get_sbdb(rec.get("sbdb_des") or des)
        el = {e["name"]: float(e["value"])
              for e in sb["orbit"]["elements"] if e.get("value")}
        el["epoch"] = float(sb["orbit"]["epoch"])
        et0 = (el["epoch"] - 2451545.0) * 86400.0
        tp_et = (el["tp"] - 2451545.0) * 86400.0
        st = np.array(sp.conics(
            [el["q"] * R.AU_KM, el["e"], math.radians(el["i"]),
             math.radians(el["om"]), math.radians(el["w"]),
             0.0, tp_et, R.GM_SUN], et0))
        r_seed = st[:3] / R.AU_KM
        v_seed = st[3:] / R.AU_KM * 86400 * R.DAY_YR

        inb = [o for o in obs if o["et"] < tp_et]
        out = [o for o in obs if o["et"] >= tp_et]

        def leg_ok(leg):
            return (len(leg) >= R.MIN_OBS_LEG and
                    (leg[-1]["et"] - leg[0]["et"]) / 86400.0
                    >= R.MIN_LEG_SPAN_D)

        legs = {"in": inb if leg_ok(inb) else None,
                "out": out if leg_ok(out) else None,
                "all": obs if leg_ok(obs) else None}
        if legs["in"] is None or legs["out"] is None:
            return None, "single-leg only"

        def anchor(r_s, v_s, et_s, et_t):
            st6 = R.propagate_states(r_s, v_s, et_s,
                                     np.array([et_t]))[0]
            return st6[:3], st6[3:]

        fits = {}
        todo = ["out", "in"] if et0 >= tp_et else ["in", "out"]
        r_s, v_s, et_s = r_seed, v_seed, et0
        for lab in todo:
            leg = legs[lab]
            mid = 0.5 * (leg[0]["et"] + leg[-1]["et"])
            ra, va = anchor(r_s, v_s, et_s, mid)
            rl, vl, rmsl, kl, _ = R.fit(ra, va, mid, leg)
            fits[lab] = dict(r=rl, v=vl, et=mid, rms=rmsl)
            r_s, v_s, et_s = rl, vl, mid
        leg = legs["all"]
        mid = 0.5 * (leg[0]["et"] + leg[-1]["et"])
        ra, va = anchor(r_seed, v_seed, et0, mid)
        rl, vl, rmsl, kl, _ = R.fit(ra, va, mid, leg)
        fits["all"] = dict(r=rl, v=vl, et=mid, rms=rmsl)
    except Exception as exc:
        return None, f"fit: {exc}"

    row = dict(desig=rec["desig"], des=des,
               tp_seed_jd=el["tp"], n_in=len(inb), n_out=len(out),
               span_in_d=(inb[-1]["et"] - inb[0]["et"]) / 86400.0,
               span_out_d=(out[-1]["et"] - out[0]["et"]) / 86400.0)
    for lab in ("in", "out", "all"):
        tp_jd, ecc = osc_tp_jd(fits[lab]["r"], fits[lab]["v"],
                               fits[lab]["et"])
        row[f"tp_{lab}_jd"] = tp_jd
        row[f"e_{lab}"] = ecc
        row[f"rms_{lab}"] = fits[lab]["rms"]
    if row["tp_in_jd"] and row["tp_out_jd"]:
        row["dtp_io_d"] = row["tp_in_jd"] - row["tp_out_jd"]
        row["dtp_in_cat_d"] = row["tp_in_jd"] - row["tp_seed_jd"]
        row["dtp_out_cat_d"] = row["tp_out_jd"] - row["tp_seed_jd"]
    return row, None


# ------------------------------------------------------------------
# cohort: the step_b91 dual-leg refit record
# ------------------------------------------------------------------

recs = {}
for line in open(CKPT).readlines() if CKPT.exists() and not REDO else []:
    try:
        d = json.loads(line)
    except json.JSONDecodeError:
        continue
    recs[d["desig"]] = d
if REDO and CKPT.exists():
    CKPT.unlink()
    recs = {}

cohort = []
b91 = RESULTS / "step_b91_refit.jsonl"
for d in load_jsonl_dedup(b91, keep=lambda d: not d.get("failed")):
    if "our_ddirf" not in d:
        continue
    cohort.append(d)
logger.info(f"dual-leg cohort from step_b91 checkpoint: {len(cohort)}")

out_rows = []
todo = []
done = {k for k, v in recs.items() if not v.get("failed")}
for d in cohort:
    if d["desig"] in done:
        out_rows.append(recs[d["desig"]])
    else:
        todo.append(d)
logger.info(f"{len(done)} checkpointed, {len(todo)} to process")


def _account(rec, rowres, err):
    with open(CKPT, "a") as ck:
        if rowres is None:
            ck.write(json.dumps({"desig": rec["desig"],
                                 "failed": True,
                                 "err": str(err)[:200]}) + "\n")
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
                if (i + 1) % 25 == 0:
                    logger.info(f"  {i+1}/{len(todo)} processed")
    else:
        for i, rec in enumerate(todo):
            rr, err = _process(rec)
            _account(rec, rr, err)
            if (i + 1) % 25 == 0:
                logger.info(f"  {i+1}/{len(todo)} processed")

logger.info(f"tp-slip record built: {len(out_rows)} comets")

# join refit + catalogue products onto the tp record
for r in out_rows:
    d = next((x for x in cohort if x["desig"] == r["desig"]), None)
    if not d:
        continue
    for k in ("our_theta", "our_ddirf", "our_d_in_leg",
              "our_d_out_leg", "our_aph", "xarc_rms_in2out", "cat_drot",
              "cat_daa", "cat_denc", "cat_q", "cat_i", "in_code",
              "yr"):
        if k in d:
            r[k] = d[k]
b84 = {}
with open(RESULTS / "step_b84_pre2018_sbdb.csv") as f:
    for rr in csv.DictReader(f):
        b84[rr["desig"]] = rr
for r in out_rows:
    rr = b84.get(r["desig"])
    if not rr:
        continue
    for k in ("dtau", "dtau_in", "dtau_out"):
        try:
            r["cat_" + k] = float(rr[k])
        except (ValueError, KeyError):
            pass

raw_have = [r for r in out_rows if "dtp_io_d" in r and
            "our_theta" in r]
# Conditioning cut: the osculating tp = et - M/n inversion is
# unbounded for the near-parabolic leg fits, so a leg whose implied
# perihelion epoch lands outside the observed arc is not measuring
# perihelion timing at all.  A row enters the tests only when both
# leg tp epochs sit within one total-arc length of the seed tp --
# symmetric in cap membership, so the cut cannot bias the geometry.
have = []
for r in raw_have:
    arc = r["span_in_d"] + r["span_out_d"]
    r["dtp_valid"] = (abs(r["dtp_in_cat_d"]) < arc and
                      abs(r["dtp_out_cat_d"]) < arc)
    if r["dtp_valid"]:
        have.append(r)
logger.info(f"{len(raw_have)} comets carry both dtp_io and theta; "
            f"{len(have)} pass the conditioning cut")

rng = np.random.default_rng(SEED)
res = {"step": "step_159_tp_slip",
       "description": ("perihelion-epoch disagreement between "
                       "independently fitted inbound and outbound "
                       "legs on raw MPC astrometry; boundary products "
                       "joined from the step_b91 checkpoint"),
       "cap_deg": CAP_DEG, "n_dual_leg": len(cohort),
       "n_with_dtp": len(raw_have), "n_conditioned": len(have),
       "conditioning": ("both leg osculating tp within one total arc "
                        "of the SBDB seed tp; near-parabolic fits "
                        "outside that window carry an unbounded "
                        "M/n inversion and are excluded from tests")}


def _perm_p(vals, inc, n=N_PERM):
    obs = float(np.median(vals[inc]) - np.median(vals[~inc]))
    cnt, idx = 1, np.arange(len(vals))
    for _ in range(n):
        rp = rng.permutation(idx)
        g = float(np.median(vals[rp[:inc.sum()]]) -
                  np.median(vals[rp[inc.sum():]]))
        if g >= obs:
            cnt += 1
    return obs, cnt / (n + 1)


def cap_test(pool, key="dtp_io_d", tag=""):
    th = np.array([r["our_theta"] for r in pool])
    v = np.abs(np.array([r[key] for r in pool]))
    inc = th < CAP_DEG
    if inc.sum() < 5 or (~inc).sum() < 5:
        return None
    u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
    obs, p = _perm_p(v, inc)
    return {"tag": tag, "n_in": int(inc.sum()),
            "n_out": int((~inc).sum()),
            "med_in_d": float(np.median(v[inc])),
            "med_out_d": float(np.median(v[~inc])),
            "p_mwu": float(u.pvalue), "med_diff_d": obs,
            "p_perm": float(p)}


# ---- T1: channel validity -----------------------------------------
t1 = {}
if len(have) > 30:
    v = np.abs([r["dtp_io_d"] for r in have])
    for src, key in (("catalogue_dtau", "cat_dtau"),
                     ("our_ddirf", "our_ddirf")):
        sub = [r for r in have if key in r and
               np.isfinite(r.get(key, np.nan))]
        if len(sub) > 30:
            rho, p = spearmanr(np.abs([r["dtp_io_d"] for r in sub]),
                               np.abs([r[key] for r in sub]))
            t1[src] = {"n": len(sub), "rho": float(rho),
                       "p": float(p)}
    t1["med_abs_dtp_io_d"] = float(np.median(v))
    t1["med_abs_dtp_io_frac_of_arc"] = float(np.median(
        np.abs([r["dtp_io_d"] for r in have]) /
        np.maximum(1.0, [r["span_in_d"] + r["span_out_d"]
                         for r in have])))
    if raw_have:
        t1["n_excluded_unconditioned"] = len(raw_have) - len(have)
        t1["excluded_cap_frac"] = float(np.mean(
            [r["our_theta"] < CAP_DEG for r in raw_have
             if not r["dtp_valid"]]))
res["T1_validity"] = t1

# ---- T2: cap contrast ----------------------------------------------
t2 = {"full": cap_test(have, tag="full"),
      "code_overlap": cap_test([r for r in have if r.get("in_code")],
                               tag="code_overlap")}
res["T2_cap_contrast"] = t2

# ---- T3: signed polarity vs bipolar lobes ---------------------------
t3 = {}
th = np.array([r["our_theta"] for r in have])
v = np.array([r["dtp_io_d"] for r in have])
inc = th < CAP_DEG
if inc.sum() >= 5:
    t3["declared_cap"] = {"n": int(inc.sum()),
                          "med_signed_d": float(np.median(v[inc])),
                          "frac_positive": float(
                              np.mean(v[inc] > 0))}
    u = mannwhitneyu(v[inc], v[~inc], alternative="two-sided")
    t3["declared_cap"]["p_signed_shift"] = float(u.pvalue)
# mirror cap: the transit channel's bipolar lobe is the antipode of
# the declared axis our_theta is measured against (TNO = (34,-13) in
# lpc_boundary), NOT AXES["anti"] -- that antipode belongs to the
# detached-TNO cluster axis (49,-17) used by the resident/ISO
# channels.  Membership is the aphelion direction inside 60 deg of
# the transit antipode.
ANTI_TNO = -np.asarray(TNO, dtype=float)
aph = np.array([
    sep(np.array(r["our_aph"]), ANTI_TNO) if "our_aph" in r
    else np.nan for r in have])
mir = np.where(np.isfinite(aph), aph < CAP_DEG, False)
if mir.sum() >= 5:
    t3["mirror_cap"] = {"n": int(mir.sum()),
                        "med_signed_d": float(np.median(v[mir])),
                        "frac_positive": float(
                            np.mean(v[mir] > 0))}
res["T3_signed_polarity"] = t3

# ---- T4: timing vs geometry dissociation ----------------------------
t4 = {}
if len(have) > 30:
    y = np.log10(np.clip(np.abs([r["dtp_io_d"] for r in have]),
                         1e-3, None))
    X = np.column_stack([
        np.ones(len(have)),
        np.log10([r["n_in"] + r["n_out"] for r in have]),
        np.log10(np.clip([r["span_in_d"] + r["span_out_d"]
                          for r in have], 1.0, None)),
        np.abs([r["span_in_d"] - r["span_out_d"]
                for r in have])])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    inc = th < CAP_DEG
    if inc.sum() >= 5 and (~inc).sum() >= 5:
        u = mannwhitneyu(resid[inc], resid[~inc],
                         alternative="greater")
        obs, p = _perm_p(resid, inc)
        t4["resid_cap"] = {"med_diff_dex": obs, "p_mwu": float(u.pvalue),
                           "p_perm": float(p)}
    rho, p = spearmanr(np.abs([r["dtp_io_d"] for r in have]),
                       [r["xarc_rms_in2out"] for r in have])
    t4["vs_xarc_rms"] = {"rho": float(rho), "p": float(p)}
res["T4_dissociation"] = t4

res["reading"] = (
    "T1 certifies the channel: the perihelion-epoch disagreement is a "
    "projection of the same inter-leg non-closure the angular record "
    "carries, so it must concord with |dtau| and ddirf.  T2/T3 test "
    "whether the timing slip carries the same cap localization and "
    "bipolar polarity; T4 checks whether it survives residualization "
    "against astrometric-noise covariates.")

with open(RESULTS / "step_b125_tp_slip.json", "w") as f:
    json.dump(res, f, indent=1)
logger.data_save(RESULTS / "step_b125_tp_slip.json")

with open(RESULTS / "step_b125_tp_slip.csv", "w", newline="") as f:
    keys = sorted({k for r in out_rows for k in r})
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in out_rows:
        w.writerow({k: r.get(k) for k in keys})
logger.data_save(RESULTS / "step_b125_tp_slip.csv")

logger.info(json.dumps(res, indent=1))
logger.info("done")
