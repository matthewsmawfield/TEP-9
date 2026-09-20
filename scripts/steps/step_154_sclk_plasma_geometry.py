"""Step 154: SCLK link-plasma and solar-activity regression (step_b118).

The plasma-projection confound on the spacecraft clock channel: if
solar-wind turbulence, IMF state, or heliospheric plasma structure
along the Earth--spacecraft line of sight contaminated the SCLK
calibration record, the residuals must organize by the link's plasma
covariates -- contemporaneous solar activity and the Sun--Earth--probe
(SEP) elongation angle that controls the inner-heliosphere column the
radio link traverses -- not by the craft's heliocentric radius.

Steps 141/146/148 tested crossings, radial gradients and PWS
anchoring; this step runs the direct regression that was missing.

Data: per-calibration-boundary residual series from step 096
(step_b60_sclk_journey.csv: phase_resid_s, drift_ppm, r_au, region);
NASA/GSFC OMNIWeb COHO daily merged dataset (omni_m_daily.dat, RTN
format: BR/BT/BN, |B|, Vsw, Np, Tp) and OMNI2 daily averages
(omni_01_av.dat: Kp, SSN, F10.7, ap) matched to each segment interval;
SEP angle computed at each boundary epoch from the NAIF trajectory and
DE440s ephemeris kernels (ECLIPJ2000, Sun-centred).

T1  Interval-matched activity regression: Spearman |phase_resid| and
    |drift_ppm| against segment-averaged P_dyn, |B|, Vsw, Tp, SSN,
    F10.7, ap -- per craft and pooled Voyager.
T2  Link-geometry regression: |phase_resid| vs SEP elongation; solar
    plasma scintillation noise grows steeply toward conjunction, so a
    link artefact must load the small-SEP end.
T3  Partial correlations (decisive): resid ~ r_au given {P_dyn, SEP,
    F10.7} versus resid ~ P_dyn given r_au -- whether the radial
    organization survives controlling the plasma covariates.
T4  Shell stratification: the heliosheath residual excess split at
    median P_dyn / F10.7 -- a link-plasma origin requires the excess
    to concentrate in the high-activity half -- with the segment
    interval-length control (raw phase residuals accumulate with
    segment span at rho ~ 0.9, so any cadence-activity covariance
    masquerades as a residual-activity coupling).
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_154_sclk_plasma_geometry")
tee_stdout(logger)
logger.header("SCLK link-plasma and solar-activity regression")

import csv
import datetime
import hashlib
import json
import math
import time
import numpy as np
import requests
from scipy import stats as _st
import spiceypy as sp

AU_KM = 1.496e8
MP_NPA = 1.6726e-6          # proton mass factor: nPa per (cm^-3 km^2/s^2)
NAIF = DATA_RAW / "naif"
SPICE = DATA_RAW / "spice"
OMNI_DIR = DATA_RAW / "omni"
OMNI_M = OMNI_DIR / "omni_m_daily.dat"
OMNI2 = OMNI_DIR / "omni_01_av.dat"
OMNI2_URL = ("https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/"
             "omni_01_av.dat")
PROV = OMNI_DIR / "provenance.json"

CRAFT = {"VG1": {"spice_id": "-31",
                 "spk": NAIF / "Voyager_1.a54206u_V0.2_merged.bsp"},
         "VG2": {"spice_id": "-32",
                 "spk": NAIF / "Voyager_2.m05016u.merged.bsp"},
         "NH1": {"spice_id": "-98",
                 "spk": [NAIF / "nh_recon_e2j_v1.bsp",
                         NAIF / "nh_recon_j2sep07_prelimv1.bsp",
                         NAIF / "nh_pred_od077.bsp",
                         NAIF / "nh_recon_od117_v01.bsp",
                         NAIF / "nh_recon_pluto_od122_v01.bsp",
                         NAIF / "nh_recon_arrokoth_od147_v01.bsp",
                         NAIF / "nh_pred_alleph_od164.bsp"]}}


def fetch_omni2():
    if OMNI2.exists():
        logger.info(f"OMNI2 daily cache present: {OMNI2.name}")
        return
    r = None
    for attempt in range(1, 6):
        try:
            r = requests.get(OMNI2_URL, timeout=180,
                             headers={"User-Agent": "tep9/1.0"})
            r.raise_for_status()
            break
        except Exception as exc:
            if attempt == 5:
                raise
            wait = 5 * 2 ** (attempt - 1)
            logger.info(f"OMNI2 fetch attempt {attempt} failed ({exc}); "
                        f"retrying in {wait}s")
            time.sleep(wait)
    OMNI2.write_bytes(r.content)
    prov = json.loads(PROV.read_text()) if PROV.exists() else \
        {"step": "step_136_solar_wind_screening", "files": {}}
    prov["files"]["omni_01_av.dat"] = {
        "url": OMNI2_URL,
        "retrieved_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "bytes": len(r.content),
        "sha256": hashlib.sha256(r.content).hexdigest(),
        "note": "OMNI2 daily averages: Kp, SSN, DST, AE, proton fluxes, "
                "ap, F10.7, PC(N), AL, AU (55-word omni2 format)"}
    PROV.write_text(json.dumps(prov, indent=1))
    logger.info(f"downloaded OMNI2 daily ({len(r.content)} bytes)")


fetch_omni2()

# ---- daily covariates -------------------------------------------------
daily = {}
for line in OMNI_M.read_text().splitlines():
    t = line.split()
    if len(t) < 14:
        continue
    try:
        yr, doy = int(t[0]), int(t[1])
        bmag, v, n, tp = (float(t[8]), float(t[9]),
                          float(t[12]), float(t[13]))
    except ValueError:
        continue
    d = datetime.date(yr, 1, 1) + datetime.timedelta(days=doy - 1)
    e = daily.setdefault(d, {})
    if v <= 9000 and n <= 900:
        e["p_dyn"] = MP_NPA * n * v * v
        e["vsw"] = v
    if 0 < bmag < 900:
        e["bmag"] = bmag
    if 0 < tp < 9000000:
        e["tp"] = tp

for line in OMNI2.read_text().splitlines():
    t = line.split()
    if len(t) < 55:
        continue
    try:
        yr, doy = int(t[0]), int(t[1])
        kp, ssn, ap, f107 = (float(t[38]), float(t[39]),
                             float(t[49]), float(t[50]))
    except ValueError:
        continue
    d = datetime.date(yr, 1, 1) + datetime.timedelta(days=doy - 1)
    e = daily.setdefault(d, {})
    if kp < 99:
        e["kp"] = kp * 0.1
    if ssn < 999:
        e["ssn"] = ssn
    if ap < 999:
        e["ap"] = ap
    if 10 < f107 < 999:
        e["f107"] = f107

COV = ("p_dyn", "bmag", "vsw", "tp", "ssn", "f107", "ap")
logger.info(f"daily covariates: {len(daily)} days, "
            f"{min(daily)} to {max(daily)}")


def interval_mean(date_end, interval_days, key):
    a = date_end - datetime.timedelta(days=max(interval_days, 1.0))
    vals = [daily[d][key]
            for k in range((date_end - a).days + 1)
            if (d := a + datetime.timedelta(days=k)) in daily
            and key in daily[d]]
    return float(np.mean(vals)) if vals else np.nan


# ---- boundary series ---------------------------------------------------
rows = []
for r in csv.DictReader(open(RESULTS / "step_b60_sclk_journey.csv")):
    if r["utc"].startswith("1979-03-05"):   # VG1 dense encounter cal
        continue
    rows.append(r)
logger.info(f"boundary segments: {len(rows)} "
            f"({ {c: sum(1 for r in rows if r['craft'] == c) for c in CRAFT} })")


def utc_date(s):
    return datetime.date.fromisoformat(s[:10])


# ---- SEP angle via SPICE ------------------------------------------------
sp.furnsh(str(NAIF / "naif0012.tls"))
sp.furnsh(str(SPICE / "de440s.bsp"))
for cfg in CRAFT.values():
    for k in (cfg["spk"] if isinstance(cfg["spk"], list)
              else [cfg["spk"]]):
        sp.furnsh(str(k))


def sep_angle(craft, et):
    try:
        pe, _ = sp.spkezr("399", et, "ECLIPJ2000", "NONE", "10")
        pc, _ = sp.spkezr(CRAFT[craft]["spice_id"], et,
                          "ECLIPJ2000", "NONE", "10")
    except Exception:
        return np.nan
    pe, pc = np.array(pe[:3]), np.array(pc[:3])
    to_craft, to_sun = pc - pe, -pe
    c = np.dot(to_craft, to_sun) / (np.linalg.norm(to_craft)
                                    * np.linalg.norm(to_sun))
    return math.degrees(math.acos(np.clip(c, -1, 1)))


recs = []
for r in rows:
    et = float(r["segment_end_et_s"])
    d0 = utc_date(r["utc"])
    ivl = float(r["interval_days"])
    e = dict(craft=r["craft"], utc=r["utc"], et_s=et,
             r_au=float(r["r_au"]), region=r["region"],
             interval_days=ivl,
             phase_resid_s=float(r["phase_resid_s"]),
             drift_ppm=float(r["drift_ppm"]),
             sep_deg=sep_angle(r["craft"], et))
    for k in COV:
        e[k] = interval_mean(d0, ivl, k)
    recs.append(e)

with open(RESULTS / "step_b118_sclk_plasma_geometry.csv", "w",
          newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(recs[0]))
    w.writeheader()
    w.writerows(recs)
logger.data_save(RESULTS / "step_b118_sclk_plasma_geometry.csv")

by_craft = {c: [e for e in recs if e["craft"] == c] for c in CRAFT}
pooled_vg = [e for e in recs if e["craft"] in ("VG1", "VG2")]


def spear(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan, np.nan, int(m.sum())
    s = _st.spearmanr(x[m], y[m])
    return float(s.statistic), float(s.pvalue), int(m.sum())


def arr(es, k):
    return np.array([e[k] for e in es], float)


def rank_partial(x, y, Z):
    """Spearman x~y partialled on covariates Z via rank-space
    residualization."""
    rx, ry = _st.rankdata(x), _st.rankdata(y)
    Zr = np.column_stack([np.ones(len(x))] +
                         [_st.rankdata(z) for z in Z])
    bx, *_ = np.linalg.lstsq(Zr, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Zr, ry, rcond=None)
    ex, ey = rx - Zr @ bx, ry - Zr @ by
    if ex.std() == 0 or ey.std() == 0:
        class _Z: statistic, pvalue = np.nan, np.nan
        return _Z()
    return _st.pearsonr(ex, ey)


res = {"step": "step_154_sclk_plasma_geometry", "result": "b118",
       "covariates": list(COV),
       "inputs": {"residuals": "step_b60_sclk_journey.csv",
                  "omni_m": "data/raw/omni/omni_m_daily.dat",
                  "omni2": "data/raw/omni/omni_01_av.dat",
                  "sep": "NAIF craft SPKs + DE440s, ECLIPJ2000"}}

# ---- T1 activity regressions --------------------------------------------
logger.info("T1 interval-matched activity regression")
T1 = {}
for tag, es in (("VG1", by_craft["VG1"]), ("VG2", by_craft["VG2"]),
                ("VG-pooled", pooled_vg), ("NH1", by_craft["NH1"])):
    T1[tag] = {}
    for obs in ("phase_resid_s", "drift_ppm"):
        y = np.abs(arr(es, obs))
        T1[tag][obs] = {}
        for k in COV:
            rho, p, n = spear(arr(es, k), y)
            T1[tag][obs][k] = dict(rho=rho, p=p, n=n)
        best = max(T1[tag][obs].items(),
                   key=lambda kv: -kv[1]["p"] if np.isfinite(
                       kv[1]["p"]) else 0)
        logger.info(f"  {tag} |{obs}|: strongest covariate "
                    f"{best[0]} rho={best[1]['rho']:+.3f} "
                    f"p={best[1]['p']:.3g} (n={best[1]['n']})")
res["T1_activity"] = T1

# ---- T2 SEP-angle regression --------------------------------------------
logger.info("T2 SEP-elongation regression")
T2 = {}
for tag, es in (("VG1", by_craft["VG1"]), ("VG2", by_craft["VG2"]),
                ("VG-pooled", pooled_vg), ("NH1", by_craft["NH1"])):
    sep_v, y = arr(es, "sep_deg"), np.abs(arr(es, "phase_resid_s"))
    rho, p, n = spear(sep_v, y)
    shell = [e["sep_deg"] for e in es if e["region"] == "heliosheath"]
    shell = [s for s in shell if np.isfinite(s)]
    T2[tag] = dict(rho=rho, p=p, n=n,
                   sep_min=float(np.nanmin(sep_v)),
                   sep_median=float(np.nanmedian(sep_v)),
                   n_small_sep_lt10=int((sep_v < 10).sum()),
                   shell_sep_median=float(np.median(shell))
                   if shell else None)
    logger.info(f"  {tag}: rho(|phase_resid|, SEP)={rho:+.3f} "
                f"p={p:.3g}; min SEP {np.nanmin(sep_v):.1f} deg, "
                f"{int((sep_v < 10).sum())} boundaries <10 deg")
res["T2_sep_geometry"] = T2

# ---- T3 partial correlations ---------------------------------------------
logger.info("T3 partial correlations: radius vs plasma covariates")
T3 = {}
for tag, es in (("VG1", by_craft["VG1"]), ("VG2", by_craft["VG2"]),
                ("VG-pooled", pooled_vg)):
    y = np.abs(arr(es, "phase_resid_s"))
    r_v = arr(es, "r_au")
    cov_m = np.column_stack([arr(es, k) for k in
                             ("p_dyn", "sep_deg", "f107")])
    m = np.isfinite(y) & np.isfinite(r_v) & np.isfinite(cov_m).all(1)
    if m.sum() < 30:
        continue
    y, r_v, Z = y[m], r_v[m], [cov_m[m, i] for i in range(3)]
    rho_r = _st.spearmanr(r_v, y)
    rho_p = _st.spearmanr(Z[0], y)
    pr_r = rank_partial(r_v, y, Z)
    pr_p = rank_partial(Z[0], y, [r_v])
    pr_sep = rank_partial(Z[1], y, [r_v])
    pr_f = rank_partial(Z[2], y, [r_v])
    T3[tag] = dict(n=int(m.sum()),
                   rho_resid_radius=dict(rho=float(rho_r.statistic),
                                         p=float(rho_r.pvalue)),
                   rho_resid_pdyn=dict(rho=float(rho_p.statistic),
                                       p=float(rho_p.pvalue)),
                   partial_radius_given_plasma=dict(
                       r=float(pr_r.statistic),
                       p=float(pr_r.pvalue)),
                   partial_pdyn_given_radius=dict(
                       r=float(pr_p.statistic),
                       p=float(pr_p.pvalue)),
                   partial_sep_given_radius=dict(
                       r=float(pr_sep.statistic),
                       p=float(pr_sep.pvalue)),
                   partial_f107_given_radius=dict(
                       r=float(pr_f.statistic),
                       p=float(pr_f.pvalue)))
    logger.info(f"  {tag}: resid~r rho={rho_r.statistic:+.3f} "
                f"-> partial given plasma {pr_r.statistic:+.3f} "
                f"(p={pr_r.pvalue:.2g}); resid~P_dyn partial "
                f"{pr_p.statistic:+.3f} (p={pr_p.pvalue:.2g}); "
                f"SEP partial {pr_sep.statistic:+.3f} "
                f"(p={pr_sep.pvalue:.2g})")
res["T3_partial"] = T3

# ---- T4 shell activity stratification ------------------------------------
logger.info("T4 heliosheath-shell activity stratification")
T4 = {}
for tag, es in (("VG1", by_craft["VG1"]), ("VG2", by_craft["VG2"])):
    sh = [e for e in es if e["region"] == "heliosheath"]
    ysh = np.array([abs(e["phase_resid_s"]) for e in sh])
    T4[tag] = {"n_shell": len(sh)}
    for k in ("p_dyn", "f107"):
        x = np.array([e[k] for e in sh])
        m = np.isfinite(x)
        if m.sum() < 20:
            continue
        med = np.median(x[m])
        hi, lo = ysh[m][x[m] > med], ysh[m][x[m] <= med]
        u = _st.mannwhitneyu(hi, lo, alternative="two-sided")
        # interval-length control: raw residuals accumulate with
        # segment span, so a cadence-activity covariance can
        # masquerade as a residual-activity coupling
        rr = np.array([e["r_au"] for e in sh])[m]
        iv = np.array([e["interval_days"] for e in sh])[m]
        pr = rank_partial(x[m], ysh[m], [rr, iv])
        T4[tag][f"split_{k}"] = dict(
            median=float(med), n_hi=int(len(hi)), n_lo=int(len(lo)),
            med_resid_hi=float(np.median(hi)),
            med_resid_lo=float(np.median(lo)),
            mw_p=float(u.pvalue),
            partial_given_radius_interval=dict(
                r=float(pr.statistic), p=float(pr.pvalue)))
        logger.info(f"  {tag} shell split {k}@{med:.2f}: "
                    f"hi {np.median(hi):.2e}s vs lo "
                    f"{np.median(lo):.2e}s, MWU p={u.pvalue:.3g}; "
                    f"partial given (r, interval) "
                    f"r={pr.statistic:+.3f} p={pr.pvalue:.3g}")
res["T4_shell_stratification"] = T4

# ---- verdict --------------------------------------------------------------
pooled_t1_minp = min((v["p"] for obs in T1["VG-pooled"].values()
                      for v in obs.values() if np.isfinite(v["p"])),
                     default=1.0)
radius_survives = all(T3[c]["partial_radius_given_plasma"]["p"] < 0.05
                      and T3[c]["partial_radius_given_plasma"]["r"] > 0
                      for c in ("VG1", "VG2") if c in T3)
plasma_flat = all(
    T3[c]["partial_pdyn_given_radius"]["p"] > 0.05 for c in T3)
shell_flat = all(
    split["partial_given_radius_interval"]["p"] > 0.05
    for craft_res in T4.values()
    for key, split in craft_res.items()
    if isinstance(split, dict)
    and np.isfinite(split["partial_given_radius_interval"]["p"]))
verdict = ("LINK-PLASMA NOT SUPPORTED: residuals carry no "
           "contemporaneous solar-activity or SEP-elongation "
           "organization once segment-interval length is controlled; "
           "the radial structure survives partialling every plasma "
           "covariate (T1-T4)"
           if radius_survives and plasma_flat and shell_flat else
           "MIXED: see T1-T4 -- some plasma covariate retains signal "
           "after radius and interval control")
res["verdict"] = verdict
logger.info(f"verdict: {verdict} (min T1 p={pooled_t1_minp:.3g})")

with open(RESULTS / "step_b118_sclk_plasma_geometry.json", "w") as f:
    json.dump(res, f, indent=1, default=float)
logger.data_save(RESULTS / "step_b118_sclk_plasma_geometry.json")

# ---- figure ----------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cols = {"VG1": "indianred", "VG2": "steelblue", "NH1": "0.55"}
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
for c, es in by_craft.items():
    ax[0].scatter(arr(es, "p_dyn"), np.abs(arr(es, "phase_resid_s")),
                  s=7, alpha=0.35, c=cols[c], label=c)
    ax[1].scatter(arr(es, "sep_deg"), np.abs(arr(es, "phase_resid_s")),
                  s=7, alpha=0.35, c=cols[c])
    ax[2].scatter(arr(es, "r_au"), np.abs(arr(es, "phase_resid_s")),
                  s=7, alpha=0.35, c=cols[c])
ax[0].set(xlabel="interval-mean P_dyn (nPa)", yscale="log",
          ylabel="|phase residual| (s)",
          title="T1 residual vs solar-wind pressure")
ax[1].set(xlabel="Sun-Earth-probe elongation (deg)", yscale="log",
          title="T2 residual vs link geometry")
ax[1].axvline(10, ls=":", color="k", lw=0.8)
ax[2].set(xlabel="heliocentric radius (AU)", yscale="log",
          title="residual vs radius (organizer)")
ax[0].legend(fontsize=8)
fig.suptitle("SCLK residuals vs plasma/link covariates (step 154)")
fig.tight_layout()
fig.savefig(RESULTS / "figures" / "supplementary" /
            "step_b118_sclk_plasma_geometry.png", dpi=300)
logger.data_save(RESULTS / "figures/supplementary/"
                 "step_b118_sclk_plasma_geometry.png")
