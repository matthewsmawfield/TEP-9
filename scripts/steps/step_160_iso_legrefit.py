"""Step 160 -- Independent two-leg refit of the interstellar objects.

Step 101 (iso_ng_channel) tests the ISOs through their fitted
non-gravitational ledger against the cap-transit geometry; step 128
built the leg-split machinery that resolves the post-2017 lapse-slip
class on raw MPC astrometry.  This step crosses the two: each ISO's
raw MPC observation set is split at perihelion, each leg is fitted
independently by the step_127 LM differential-correction instrument
(REBOUND/IAS15 + DE440s, one-way light time, topocentric parallax,
MAD clip), and each fitted leg is propagated to the 250 AU
barycentric sphere through the shared integrate_leg machinery.

The product is a zero-catalogue boundary record per ISO: per-leg
asymptote directions and their cap classifications, the inter-leg
non-closure ddirf, the leg-vs-osculating deviations d_in_leg /
d_out_leg, the cross-arc prediction loss, and the perihelion-epoch
disagreement between legs (the step_159 timing channel).

3I/ATLAS is the case of interest: its outbound asymptote sits inside
the declared cap (46.9 deg from the axis, step 101), so its outbound
leg realises a post-2017-style in-cap transit through the machinery
that carries the registered anomaly class.  Its fitted NG parameters
are retained as the volatile-confound context: a real outgassing
acceleration inflates leg disagreement through solver absorption
symmetrically, so the interpretation registers magnitude and
localization jointly, not either alone.

Outputs
  results/step_b126_iso_legrefit.json
  results/step_b126_iso_legrefit.csv
"""

import json
import math
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

import numpy as np

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.tep9_common import (
    DATA_RAW, RESULTS, AXES, perih_dir, sep, lb)

logger = StepLogger("step_160_iso_legrefit")
tee_stdout(logger)

SPK = DATA_RAW / "spice" / "de440s.bsp"
LSK = DATA_RAW / "naif" / "naif0012.tls"
MPC_DIR = DATA_RAW / "mpc"
OBS_DIR = MPC_DIR / "obs"
SBDB_DIR = MPC_DIR / "sbdb_fp"
OBSC_PATH = MPC_DIR / "obscodes.json"
PROV_PATH = MPC_DIR / "provenance.json"

CAP_DEG = 60.0
EDGE_HI = 75.0

logger.header("Independent two-leg refit of the interstellar objects")

from scripts.utils import mpc_refit as R
from scripts.utils.lpc_boundary import integrate_leg
import spiceypy as sp

R.configure(OBS_DIR, SBDB_DIR, OBSC_PATH, PROV_PATH, SPK, LSK,
            "step_160_iso_legrefit")

# MPC observation designators differ from the SBDB designators for
# the ISOs; each entry lists candidates tried in order.
ISOS = [
    {"label": "1I/'Oumuamua", "sbdb": "1I",
     "obs_desigs": ["1I", "A/2017 U1", "1I/'Oumuamua"]},
    {"label": "2I/Borisov", "sbdb": "2I",
     "obs_desigs": ["2I", "C/2019 Q4"]},
    {"label": "3I/ATLAS", "sbdb": "3I",
     "obs_desigs": ["C/2025 N1", "3I", "3I/ATLAS"]},
]


def resolve_obs_des(cands):
    """Return the first designator yielding a usable MPC record."""
    for des in cands:
        try:
            raw = R.get_obs(des)
            obs = R.parse_obs(raw)
        except Exception as exc:
            logger.info(f"  get-obs {des}: {exc}")
            continue
        if len(obs) >= 2 * R.MIN_OBS_LEG:
            return des, len(obs)
        logger.info(f"  get-obs {des}: only {len(obs)} usable obs")
    return None, 0


def osc_tp_jd(r, v, et):
    el = sp.oscelt(np.r_[r * R.AU_KM, v * R.AU_KM / (R.DAY_YR * 86400.0)],
                   et, R.GM_SUN)
    rp, ecc, m0 = el[0], el[1], el[5]
    if abs(1.0 - ecc) < 1e-5:
        return None
    a = rp / (1.0 - ecc)
    n = math.sqrt(R.GM_SUN / abs(a) ** 3)
    return (et - m0 / n) / 86400.0 + 2451545.0


def cap_class(p):
    """Sector membership of a unit direction under the declared caps."""
    sa = sep(p, AXES["tno"])
    sm = sep(p, AXES["anti"])
    in_cap = sa < CAP_DEG
    in_mir = sm < CAP_DEG
    if in_cap:
        sector = "cap"
    elif in_mir:
        sector = "mirror_cap"
    elif sa < EDGE_HI or sm < EDGE_HI:
        sector = "edge"
    else:
        sector = "field"
    lon, lat = lb(p)
    return dict(ecl_lon=round(lon, 2), ecl_lat=round(lat, 2),
                sep_axis_deg=round(sa, 2), sep_antiaxis_deg=round(sm, 2),
                in_cap=bool(in_cap), in_mirror_cap=bool(in_mir),
                sector=sector)


ledger = []
for iso in ISOS:
    logger.info(f"--- {iso['label']} ---")
    obs_des, n_obs = resolve_obs_des(iso["obs_desigs"])
    if obs_des is None:
        ledger.append({"label": iso["label"], "sbdb": iso["sbdb"],
                       "error": "no usable MPC record under any desig"})
        continue
    logger.info(f"  obs desig: {obs_des} ({n_obs} usable obs)")
    res, err = R.fit_comet(obs_des, sbdb_des=iso["sbdb"])
    if res is None:
        ledger.append({"label": iso["label"], "sbdb": iso["sbdb"],
                       "obs_des": obs_des, "error": err})
        logger.info(f"  fit failed: {err}")
        continue

    rec = {"label": iso["label"], "sbdb": iso["sbdb"],
           "obs_des": obs_des, "n_obs": res["n_obs"],
           "n_in": res["n_in"], "n_out": res["n_out"],
           "tp_seed_jd": res["tp_jd"],
           "elements": {k: res["elements"][k] for k in
                        ("q", "e", "i", "om", "w", "tp")
                        if k in res["elements"]}}
    for k in ("fit_all_rms", "fit_all_n", "drot", "d_in", "d_out",
              "daa", "denc", "theta", "dtau", "ddirf", "d_in_leg",
              "d_out_leg", "xarc_rms_in2out", "in_rms", "out_rms",
              "in_nkeep", "out_nkeep", "boundary_err"):
        if k in res:
            rec[k] = res[k]

    # NG context from the SBDB seed (model_pars), where fitted
    sb = R.get_sbdb(iso["sbdb"])
    mp = {p["name"]: {"value": float(p["value"]),
                      "sigma": (float(p["sigma"])
                                if p.get("sigma") else None)}
          for p in sb["orbit"].get("model_pars", [])
          if p.get("value") is not None}
    if mp:
        rec["sbdb_model_pars"] = mp
    rec["soln_date"] = sb["orbit"].get("soln_date")
    rec["orbit_id"] = sb["orbit"].get("orbit_id")

    # per-leg osculating perihelion epochs (step_159 channel)
    fits = res.get("fits", {})
    for lab in ("in", "out", "all"):
        if lab in fits:
            rec[f"tp_{lab}_jd"] = osc_tp_jd(
                np.array(fits[lab]["r"]), np.array(fits[lab]["v"]),
                fits[lab]["et"])
    if rec.get("tp_in_jd") and rec.get("tp_out_jd"):
        rec["dtp_io_d"] = rec["tp_in_jd"] - rec["tp_out_jd"]

    # per-leg boundary periapsis directions and their cap
    # classification -- phat is the osculating periapsis direction of
    # each leg's own solution at the 250 AU sphere, the same quantity
    # the step_101 mirror-lobe classification uses
    try:
        if "in" in fits and "out" in fits:
            fi, fo = fits["in"], fits["out"]
            rbi = integrate_leg(np.array(fi["r"]), np.array(fi["v"]),
                                fi["et"], -1, v_max=40.0)
            rfo = integrate_leg(np.array(fo["r"]), np.array(fo["v"]),
                                fo["et"], +1, v_max=40.0)
            if rbi and rfo:
                rec["inleg_periapsis"] = cap_class(np.array(rbi["phat"]))
                rec["outleg_periapsis"] = cap_class(np.array(rfo["phat"]))
                rec["ddirf_legs"] = sep(np.array(rbi["phat"]),
                                        np.array(rfo["phat"]))
    except Exception as exc:
        rec["leg_boundary_err"] = str(exc)[:150]
    ledger.append(rec)

# ------------------------------------------------------------------
# cohort context: where the ISO products sit inside the post-2017
# dual-leg refit distribution (step_b92)
# ------------------------------------------------------------------

coh = []
b92 = RESULTS / "step_b92_refit.jsonl"
if b92.exists():
    for line in open(b92):
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not d.get("failed"):
            coh.append(d)
logger.info(f"post-2017 dual-leg cohort: {len(coh)} records")


def pctile(key, val):
    v = np.array([c[key] for c in coh
                  if key in c and np.isfinite(c[key])])
    if len(v) < 20 or val is None:
        return None
    return {"cohort_n": int(len(v)),
            "cohort_med": float(np.median(v)),
            "cohort_p95": float(np.percentile(v, 95)),
            "iso_value": float(val),
            "iso_percentile": float(np.mean(v <= val) * 100)}


ctx = {}
for rec in ledger:
    if "error" in rec:
        continue
    ctx[rec["label"]] = {
        "ddirf": pctile("our_ddirf", rec.get("ddirf")),
        "d_in_leg": pctile("our_d_in_leg", rec.get("d_in_leg")),
        "d_out_leg": pctile("our_d_out_leg", rec.get("d_out_leg")),
        "xarc_rms_in2out": pctile("xarc_rms_in2out",
                                  rec.get("xarc_rms_in2out")),
        "fit_all_rms": pctile("fit_all_rms", rec.get("fit_all_rms"))}

out = {"step": "step_160_iso_legrefit",
       "description": ("independent two-leg LM refit of the ISOs on "
                       "raw MPC astrometry; per-leg boundary "
                       "asymptotes, cap classification, inter-leg "
                       "non-closure and perihelion-epoch slip; "
                       "cohort percentiles vs the post-2017 refit "
                       "record"),
       "cap_deg": CAP_DEG, "edge_hi_deg": EDGE_HI,
       "ledger": ledger, "cohort_context": ctx,
       "caveats": [
        "pure-gravity leg fits absorb real NG outgassing into the "
        "fitted states; for the volatile-active ISOs (2I, 3I) the "
        "leg disagreement is therefore a joint lapse+outgassing "
        "channel and is interpreted against the SBDB NG ledger, "
        "not as a clean datum.",
        "the ISO arcs are short (months), so leg fits carry wider "
        "element covariance than the multi-apparition comet cohort; "
        "cohort percentiles are contextual, not calibrated nulls."]}

with open(RESULTS / "step_b126_iso_legrefit.json", "w") as f:
    json.dump(out, f, indent=1, default=str)
logger.data_save(RESULTS / "step_b126_iso_legrefit.json")

import csv
keys = sorted({k for r in ledger for k in r
               if not isinstance(r[k], (dict, list))})
with open(RESULTS / "step_b126_iso_legrefit.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in ledger:
        w.writerow({k: r.get(k) for k in keys})
logger.data_save(RESULTS / "step_b126_iso_legrefit.csv")

logger.info(json.dumps(out, indent=1, default=str))
logger.info("done")
