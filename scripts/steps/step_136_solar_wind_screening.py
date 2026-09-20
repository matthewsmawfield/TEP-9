"""Step 136: Solar-wind screening-epoch test (step_b100).

Tests the plasma-breathing remediation hypothesis for the era
polarity structure: if the proper-time boundary is screened by
local heliospheric plasma density, the comet transit structure
should track the contemporaneous solar-wind dynamic pressure at
each perihelion epoch.  Solar Cycle 24 peaked in 2014 and declined
through the 2018-2021 minimum -- the same window in which the
declared/displaced polarity transition sits (step 130).

Data: NASA/GSFC OMNIWeb COHO daily merged dataset
(omni_m_daily.dat), provenance-pinned into data/raw/omni/.
Dynamic pressure P_dyn = m_p * n * v^2 in nPa, monthly averaged.

T1  Pooled per-comet residual ~ P_dyn Spearman correlation.

T2  Within-era decomposition of the apex-side shell-band
    correlation -- the era-proxy control: if the pooled
    correlation is only the era contrast restated, within-era
    coupling vanishes.

T3  Era P_dyn contrast -- documents the confound directly:
    post-2017 perihelia sample the cycle-24 decline/minimum.

T4  Regime-split axis test: dipole polarity of the residual
    field about the bipolar axis in high- versus low-pressure
    perihelion subsets (median split), the direct screen test.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_136_solar_wind_screening")
tee_stdout(logger)
logger.header("Solar-wind screening-epoch test")

import datetime
import time
import hashlib
import json
import math
import numpy as np
import requests
from scipy import stats as _st
from scripts.utils.tep9_common import lv, sep
from scripts.utils.coordinates import GAL2ECL

OMNI_URL = ("https://spdf.gsfc.nasa.gov/pub/data/omni/"
            "low_res_omni/omni_m_daily.dat")
OMNI_DIR = RESULTS.parent / "data" / "raw" / "omni"
OMNI_FILE = OMNI_DIR / "omni_m_daily.dat"
PROV = OMNI_DIR / "provenance.json"
MP_MEV = 1.6726e-6  # proton mass in nPa per (cm^-3 * km^2/s^2)


def fetch_omni():
    if OMNI_FILE.exists():
        logger.info(f"OMNI cache present: {OMNI_FILE.name}")
        return
    OMNI_DIR.mkdir(parents=True, exist_ok=True)
    r = None
    for attempt in range(1, 6):
        try:
            r = requests.get(OMNI_URL, timeout=180,
                             headers={"User-Agent": "tep9/1.0"})
            r.raise_for_status()
            break
        except Exception as exc:
            if attempt == 5:
                raise
            wait = 5 * 2 ** (attempt - 1)
            logger.info(f"OMNI fetch attempt {attempt} failed ({exc}); "
                        f"retrying in {wait}s")
            time.sleep(wait)
    OMNI_FILE.write_bytes(r.content)
    prov = {"step": "step_136_solar_wind_screening",
            "files": {"omni_m_daily.dat": {
                "url": OMNI_URL,
                "retrieved_utc": datetime.datetime.now(
                    datetime.timezone.utc).isoformat(),
                "bytes": len(r.content),
                "sha256": hashlib.sha256(r.content).hexdigest(),
                "note": "COHO daily merged OMNI: IMF + solar-wind "
                        "plasma, 1963-present"}}}
    PROV.write_text(json.dumps(prov, indent=1))
    logger.info(f"downloaded OMNI daily ({len(r.content)} bytes)")


fetch_omni()

# monthly dynamic pressure from daily speed + density
mm = {}
for line in OMNI_FILE.read_text().splitlines():
    t = line.split()
    if len(t) < 14:
        continue
    try:
        yr, doy = int(t[0]), int(t[1])
        v, n = float(t[9]), float(t[12])
    except ValueError:
        continue
    if v > 9000 or n > 900:
        continue
    d = datetime.date(yr, 1, 1) + datetime.timedelta(days=doy - 1)
    mm.setdefault((d.year, d.month), []).append(MP_MEV * n * v * v)
P = {k: float(np.mean(v)) for k, v in mm.items()}
logger.info(f"OMNI monthly P_dyn: {len(P)} months, "
            f"{min(P)} to {max(P)}")


def pdyn(jd):
    if not np.isfinite(jd):
        return np.nan
    d = datetime.date(2000, 1, 1) + datetime.timedelta(
        days=float(jd) - 2451545.0)
    return P.get((d.year, d.month), np.nan)


def load_ckpt(f):
    out = {}
    for line in open(RESULTS / f):
        r = json.loads(line)
        if r.get("failed"):
            continue
        out[r["des"]] = r
    return list(out.values())


def resid_logrot(rows):
    drot = np.array([r.get("our_drot", np.nan) for r in rows])
    daa = np.array([r.get("our_daa", np.nan) for r in rows])
    denc = np.array([r.get("our_denc", np.nan) for r in rows])
    # signed log rotation: reversed (negative) legs keep their
    # magnitude instead of being floored into a degenerate channel
    y = np.sign(drot) * np.log10(np.clip(np.abs(drot), 1e-12, None))
    X = np.column_stack([np.ones(len(y)), np.log10(np.abs(daa) + 1),
                         np.log10(np.clip(denc, 1e-12, None))])
    ok = np.isfinite(X).all(1) & np.isfinite(y)
    with np.errstate(all="ignore"):
        c, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    c = np.where(np.isfinite(c), c, 0.0)
    r = np.full(len(y), np.nan)
    with np.errstate(all="ignore"):
        r[ok] = y[ok] - X[ok] @ c
    return r


def gv(l, b):
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l),
                               math.sin(b)])


CMB = gv(264.02, 48.25)

eras = {}
for nm, f in (("pre2018", "step_b91_refit.jsonl"),
              ("post2017", "step_b92_refit.jsonl")):
    rows = load_ckpt(f)
    res = resid_logrot(rows)
    aph = np.array([r.get("our_aph", [np.nan] * 3) for r in rows])
    pdy = np.array([pdyn(r.get("tp_jd", np.nan)) for r in rows])
    yrs = np.array([r.get("yr", np.nan) for r in rows])
    ok = np.isfinite(res) & np.isfinite(aph).all(1) & np.isfinite(pdy)
    u = np.nan_to_num(aph[ok], nan=0.0, posinf=0.0, neginf=0.0)
    nrm = np.linalg.norm(u, axis=1, keepdims=True)
    u = np.divide(u, nrm, out=np.zeros_like(u), where=nrm > 0)
    with np.errstate(all="ignore"):
        cone = np.degrees(np.arccos(np.clip(u @ CMB, -1, 1)))
    eras[nm] = dict(res=res[ok], pd=pdy[ok], cone=cone,
                    yrs=yrs[ok], u=u, n=int(ok.sum()))

# ---- T1 pooled ----------------------------------------------------
logger.info("T1 pooled residual ~ P_dyn")
r_all = np.concatenate([eras[k]["res"] for k in eras])
p_all = np.concatenate([eras[k]["pd"] for k in eras])
c_all = np.concatenate([eras[k]["cone"] for k in eras])
m = np.isfinite(r_all)
rho1, p1 = _st.spearmanr(r_all[m], p_all[m])
band = (c_all >= 45) & (c_all < 60)
mb = m & band
rho1b, p1b = _st.spearmanr(r_all[mb], p_all[mb])
logger.info(f"  pooled n={m.sum()}: rho={rho1:+.3f} p={p1:.3g}; "
            f"apex-band n={mb.sum()}: rho={rho1b:+.3f} p={p1b:.3g}")
T1 = dict(pooled=dict(n=int(m.sum()), rho=float(rho1), p=float(p1)),
          apex_band=dict(n=int(mb.sum()), rho=float(rho1b),
                         p=float(p1b)))

# ---- T2 within-era control ----------------------------------------
logger.info("T2 within-era decomposition (era-proxy control)")
T2 = {}
for nm, e in eras.items():
    b = (e["cone"] >= 45) & (e["cone"] < 60)
    rb, pb_ = e["res"][b], e["pd"][b]
    mm2 = np.isfinite(rb)
    rho, p = _st.spearmanr(rb[mm2], pb_[mm2])
    rp, pp = _st.spearmanr(e["res"], e["pd"])
    T2[nm] = dict(apex_band=dict(n=int(mm2.sum()), rho=float(rho),
                                 p=float(p)),
                  pooled=dict(n=int(len(rb)), rho=float(rp),
                              p=float(pp)))
    logger.info(f"  {nm}: apex-band n={mm2.sum()} rho={rho:+.3f} "
                f"p={p:.3g}; pooled rho={rp:+.3f} p={pp:.3g}")

# ---- T3 era P_dyn contrast ----------------------------------------
logger.info("T3 era P_dyn contrast (the confound)")
pa, pb = eras["pre2018"]["pd"], eras["post2017"]["pd"]
u3 = _st.mannwhitneyu(pa, pb, alternative="two-sided")
T3 = dict(pre_median=float(np.median(pa)),
          post_median=float(np.median(pb)),
          mw_p=float(u3.pvalue),
          note="post-2017 perihelia sample the SC24 decline/minimum")
logger.info(f"  pre {np.median(pa):.2f} vs post {np.median(pb):.2f} "
            f"nPa, MWU p={u3.pvalue:.3g}")

# ---- T4 regime-split dipole polarity -------------------------------
logger.info("T4 regime-split axis test (median P_dyn split)")
pooled_p = np.concatenate([eras[k]["pd"] for k in eras])
pmed = float(np.median(pooled_p))
T4 = {"p_split_npa": pmed, "regimes": {}}
u_all = np.concatenate([eras[k]["u"] for k in eras])
with np.errstate(all="ignore"):
    cs = u_all @ CMB
for tag, sel in (("high", p_all[m] > pmed), ("low", p_all[m] <= pmed)):
    # exact least-squares slope of the residual on the dipole
    # component (r is a centred residual, so no intercept term)
    b = float(np.sum(r_all[m][sel] * cs[sel]) /
              np.sum(cs[sel] ** 2))
    T4["regimes"][tag] = dict(n=int(sel.sum()), dipole_dex=b)
    logger.info(f"  {tag}-P_dyn regime: n={sel.sum()}, "
                f"bipolar dipole b={b:+.3f} dex")

# verdict
weak = T1["apex_band"]["p"] < 0.05 and all(
    T2[k]["apex_band"]["p"] > 0.05 for k in T2)
verdict = ("SCREENING NOT SUPPORTED: pooled correlation is an "
           "era-proxy artefact; within-era coupling is flat"
           if weak else
           "SCREENING FLAT: no contemporaneous pressure coupling")
logger.info(f"verdict: {verdict}")

out = dict(step="step_136_solar_wind_screening", result="b100",
           verdict=verdict,
           omni=dict(url=OMNI_URL, n_months=len(P),
                     span=[f"{min(P)}", f"{max(P)}"],
                     p_dyn_median_npa=pmed),
           T1_pooled=T1, T2_within_era=T2, T3_era_contrast=T3,
           T4_regime_split=T4)
with open(RESULTS / "step_b100_solar_wind_screening.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
logger.data_save(RESULTS / "step_b100_solar_wind_screening.json")

FIG = RESULTS / "figures"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# T1 pooled scatter
ax[0].scatter(p_all[m], r_all[m], s=8, alpha=0.35, c="0.4")
ax[0].scatter(p_all[mb], r_all[mb], s=14, alpha=0.6, c="indianred")
ax[0].set_xlabel("perihelion-epoch P_dyn (nPa)")
ax[0].set_ylabel("residual (dex)")
ax[0].set_title(f"T1 pooled: rho={rho1:+.2f} (band {rho1b:+.2f})")

# T3 era pressure distributions
ax[1].hist([pa, pb], bins=30, label=["pre-2018", "post-2017"],
           alpha=0.65)
ax[1].set_xlabel("P_dyn (nPa)")
ax[1].set_title(f"T3 era contrast p={u3.pvalue:.1g}")
ax[1].legend(fontsize=8)

# T4 monthly P_dyn timeline with comet switch marked
yrs_m = sorted(P)
x = [k[0] + k[1] / 12 for k in yrs_m]
ax[2].plot(x, [P[k] for k in yrs_m], lw=0.7, color="0.35")
ax[2].axvline(2018, color="indianred", ls="--", lw=0.9,
              label="comet switch")
ax[2].axvline(2014.3, color="gold", ls=":", lw=0.9, label="SC24 max")
ax[2].axvline(2019.9, color="steelblue", ls=":", lw=0.9,
              label="SC24/25 min")
ax[2].set_xlim(1995, 2027)
ax[2].set_xlabel("year")
ax[2].set_ylabel("monthly P_dyn (nPa)")
ax[2].set_title("T4 pressure timeline vs switch epoch")
ax[2].legend(fontsize=8)

fig.tight_layout()
fig.savefig(FIG / "supplementary" / "step_b100_solar_wind_screening.png", dpi=300)
logger.data_save(RESULTS / "figures/supplementary/step_b100_solar_wind_screening.png")
