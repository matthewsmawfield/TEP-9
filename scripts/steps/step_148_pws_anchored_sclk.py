"""Step 148: PWS-anchored SCLK epoch-lock test (step_b112).

Step 141 found negative signed-oscillator steps in the 3-yr windows
centred on the four canonical boundary crossings (joint negative-tail
Fisher p ~ 0.0024).  The sharpening question is whether the clock
steps are *locked to the wall*: does the strongest negative step near
each crossing sit at the crossing epoch rather than anywhere in the
window?  This step prices the epoch alignment directly.

T1  Epoch lock at the canonical crossings: within +/-2 yr of each of
    V1 TS/HP and V2 TS/HP, scan the windowed signed-step statistic
    (w = 3 yr, as step 141, detrended) and take the strongest
    negative changepoint; the statistic is its offset from the
    canonical crossing epoch.  The null places the same window at
    random target epochs over each record and measures the same
    best-offset -- pricing "how close to target does the strongest
    nearby step typically land".
T2  Joint lock statistic: Fisher combination of the four per-crossing
    offset tail fractions.
T3  Plasma-event coincidence (V1 only): the PWS VLISM density series
    carries its own measured step epochs (step_097 top events).
    Distance from each PWS step to the nearest strong negative SCLK
    changepoint, against the same distance for randomly placed
    epoch sets -- testing whether clock events and plasma events
    share epochs in the VLISM record.
T4  Caveat ledger: cadence gaps near crossings and kernel
    segmentation remain the standing confounds; the channel
    registers the epoch evidence alongside them.

Outputs: results/step_b112_pws_anchored_sclk.json/.csv and
results/figures/supplementary/step_b112_pws_anchored_sclk.png.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_148_pws_anchored_sclk")
tee_stdout(logger)
logger.header("PWS-anchored SCLK epoch-lock test")

import csv
import json
import math
import numpy as np
from scipy import stats as _st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20260919
rng = np.random.default_rng(SEED)
W = 3.0
HALF_WIN = 2.0
N_NULL = 2000

RATES = RESULTS / "step_b60_sclk_rates.csv"
rows = [r for r in csv.DictReader(open(RATES))
        if r["mode"] == "calibration"]
CROSS = {"VG1": [("TS", 2004.96, 94.0), ("HP", 2012.65, 121.6)],
         "VG2": [("TS", 2007.66, 84.0), ("HP", 2018.90, 119.0)]}


def series(craft):
    rs = [r for r in rows if r["craft"] == craft]
    et = np.array([float(r["et_s"]) for r in rs])
    dnu = np.array([float(r["dnu_frac"]) for r in rs]) * 1e6
    return et / 31557600 + 2000.0, dnu


def detrend(yr, dnu, exclude):
    keep = ~exclude
    age = yr[keep] - yr.min() + 0.5
    X = np.column_stack([np.ones(keep.sum()), np.log(age),
                         np.log(age) ** 2])
    c, *_ = np.linalg.lstsq(X, dnu[keep], rcond=None)
    Xf = np.column_stack([np.ones(len(yr)), np.log(yr - yr.min() + 0.5),
                          np.log(yr - yr.min() + 0.5) ** 2])
    return dnu - Xf @ c


def step_stat(yr, dnu, tc, w):
    a = (yr >= tc - w) & (yr < tc)
    b = (yr >= tc) & (yr <= tc + w)
    if a.sum() < 8 or b.sum() < 8:
        return None
    return float(dnu[b].mean() - dnu[a].mean())


def best_neg_offset(yr, dnu, t_target, half=HALF_WIN, w=W):
    """Strongest negative step within +/-half of target -> (offset, step)."""
    ts = np.arange(t_target - half, t_target + half, 0.1)
    best, bt = 0.0, None
    for tc in ts:
        s = step_stat(yr, dnu, tc, w)
        if s is not None and s < best:
            best, bt = s, tc
    if bt is None:
        return None, 0.0
    return abs(bt - t_target), best


res = {"w_yr": W, "half_win_yr": HALF_WIN, "seed": SEED,
       "n_null": N_NULL, "test_summary": {}}

# ---- T1/T2: epoch lock ------------------------------------------------------
res["T1_epoch_lock"] = {}
tail_fracs = []
for craft in ("VG1", "VG2"):
    yr, dnu = series(craft)
    ex = np.zeros(len(yr), bool)
    for _, t_hp, _ in CROSS[craft]:
        ex |= np.abs(yr - t_hp) <= max(HALF_WIN, W)
    rd = detrend(yr, dnu, ex)
    for name, t_hp, r_hp in CROSS[craft]:
        off, step = best_neg_offset(yr, rd, t_hp)
        # null: random target epochs, same record, same procedure
        nt = 0
        t_lo, t_hi = yr.min() + W + HALF_WIN, yr.max() - W - HALF_WIN
        targets = rng.uniform(t_lo, t_hi, N_NULL)
        for tt in targets:
            o2, _ = best_neg_offset(yr, rd, tt)
            if o2 is not None and off is not None and o2 <= off:
                nt += 1
        p = (nt + 1) / (N_NULL + 1)
        tail_fracs.append(p)
        res["T1_epoch_lock"][f"{craft}_{name}"] = {
            "crossing_epoch": t_hp, "radius_au": r_hp,
            "best_neg_step_ppm": float(step),
            "offset_from_crossing_yr": float(off) if off else None,
            "offset_tail_frac": float(p)}
        print(f"{craft} {name}: strongest neg step {step:+.2f} ppm "
              f"at offset {off:.2f} yr (tail {p:.4f})")

fisher = float(_st.combine_pvalues(tail_fracs, method="fisher")[1])
res["T2_joint"] = {"tail_fracs": tail_fracs,
                   "fisher_p": fisher}
print(f"T2 joint epoch-lock Fisher p = {fisher:.4f}")

# ---- T3: V1 PWS plasma-step coincidence ---------------------------------------
pws = json.load(open(RESULTS / "step_b61_pws_channel.json"))
pws_epochs_raw = []
for s in pws["craft"]["VG1"]["top_steps"]:
    y0 = float(s["utc0"][:4]) + (float(s["utc0"][5:7]) - 0.5) / 12.0
    pws_epochs_raw.append(y0)
# merge steps belonging to the same event (contiguous samples split
# into several top_steps entries) -- keep the earliest epoch
pws_epochs = []
for y0 in sorted(pws_epochs_raw):
    if not pws_epochs or y0 - pws_epochs[-1] > 1.0:
        pws_epochs.append(y0)
res["T3_pws_coincidence"] = {"pws_step_epochs": pws_epochs,
                             "n_top_steps_raw": len(pws_epochs_raw)}
if pws_epochs:
    yr, dnu = series("VG1")
    ex = np.zeros(len(yr), bool)
    for _, t_hp, _ in CROSS["VG1"]:
        ex |= np.abs(yr - t_hp) <= max(HALF_WIN, W)
    rd = detrend(yr, dnu, ex)
    # catalogue of strong negative changepoints on the V1 record
    ts = np.arange(yr.min() + W + 0.5, yr.max() - W - 0.5, 0.25)
    cps = np.array([(tc, step_stat(yr, rd, tc, W)) for tc in ts])
    cps = cps[np.isfinite(cps[:, 1])]
    neg = cps[cps[:, 1] < cps[:, 1].mean() - 2 * cps[:, 1].std()]
    neg = neg[neg[:, 1] < 0]
    res["T3_pws_coincidence"]["n_strong_neg_changepoints"] = int(len(neg))
    if len(neg):
        dists = np.array([min(abs(pe - t) for t in neg[:, 0])
                          for pe in pws_epochs])
        # null: same count of random changepoint epochs over the
        # scanned range (not the span of the detected changepoints)
        nd = []
        for _ in range(N_NULL):
            fake = rng.uniform(ts.min(), ts.max(), len(neg))
            nd.append(np.array([min(abs(pe - t) for t in fake)
                                for pe in pws_epochs]).mean())
        obs_mean = float(dists.mean())
        p3 = float((np.array(nd) <= obs_mean).sum() + 1) / (N_NULL + 1)
        res["T3_pws_coincidence"].update({
            "mean_dist_yr": obs_mean,
            "per_step_dist_yr": [float(x) for x in dists],
            "null_mean_dist_yr": float(np.mean(nd)),
            "p_closer_than_random": p3,
            "direction": ("coincident" if p3 < 0.5 else "anticoincident"),
            "note": "p is the fraction of random-changepoint null "
                    "draws at least as close as observed; p~0 means "
                    "the real negative changepoints sit CLOSER to "
                    "the plasma steps than random placements -- "
                    "the epoch-coincidence direction; p~1 means the "
                    "changepoints are farther than random"})
        print(f"T3: PWS steps mean {obs_mean:.2f} yr from nearest "
              f"negative SCLK changepoint (null {np.mean(nd):.2f}, "
              f"p={p3:.4f})")

# ---- verdict ---------------------------------------------------------------------
verdict = (
    "EPOCH-LOCK " +
    ("CONFIRMED" if fisher < 0.05 else "NOT ESTABLISHED") +
    f": strongest negative steps near the four crossings land at "
    f"offsets {[res['T1_epoch_lock'][k]['offset_from_crossing_yr'] for k in res['T1_epoch_lock']]} yr "
    f"from the canonical epochs (joint Fisher p={fisher:.4f}); "
    + (f"V1 PWS plasma steps sit {res['T3_pws_coincidence'].get('mean_dist_yr', float('nan')):.2f} yr "
       f"from the nearest strong negative SCLK changepoints "
       f"(p={res['T3_pws_coincidence'].get('p_closer_than_random', float('nan')):.4f}). " if res["T3_pws_coincidence"].get("mean_dist_yr") else "")
    + "The step amplitude evidence (step 141) and the epoch evidence "
      "(this step) are registered separately; kernel segmentation "
      "remains the standing caveat.")
res["verdict"] = verdict
res["test_summary"] = {
    "joint_epoch_lock_fisher_p": fisher,
    "per_crossing_offsets_yr":
        {k: v["offset_from_crossing_yr"]
         for k, v in res["T1_epoch_lock"].items()},
    "per_crossing_tails":
        {k: v["offset_tail_frac"]
         for k, v in res["T1_epoch_lock"].items()},
    "pws_coincidence_p":
        res["T3_pws_coincidence"].get("p_closer_than_random"),
}

out = RESULTS / "step_b112_pws_anchored_sclk.json"
with open(out, "w") as _fh:
    json.dump(res, _fh, indent=1)

with open(RESULTS / "step_b112_pws_anchored_sclk.csv", "w",
          newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["crossing", "epoch", "best_neg_step_ppm",
                "offset_yr", "tail_frac"])
    for k, v in res["T1_epoch_lock"].items():
        w.writerow([k, v["crossing_epoch"], v["best_neg_step_ppm"],
                    v["offset_from_crossing_yr"], v["offset_tail_frac"]])

fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
ax = axes[0]
ks = list(res["T1_epoch_lock"])
offs = [res["T1_epoch_lock"][k]["offset_from_crossing_yr"] for k in ks]
tails = [res["T1_epoch_lock"][k]["offset_tail_frac"] for k in ks]
ax.bar(ks, offs, color=["crimson" if t < 0.05 else "steelblue"
                        for t in tails])
ax.set_ylabel("offset of strongest neg step (yr)")
ax.set_title(f"epoch lock (joint p={fisher:.3f})")
for i, t in enumerate(tails):
    ax.text(i, offs[i] + 0.05, f"{t:.3f}", ha="center", fontsize=7)

ax = axes[1]
if res["T3_pws_coincidence"].get("per_step_dist_yr"):
    ax.bar(range(len(res["T3_pws_coincidence"]["per_step_dist_yr"])),
           res["T3_pws_coincidence"]["per_step_dist_yr"],
           color="seagreen")
    ax.axhline(res["T3_pws_coincidence"]["null_mean_dist_yr"],
               color="k", ls=":", lw=1, label="random-epoch mean")
    ax.set_xlabel("V1 PWS step #")
    ax.set_ylabel("dist to nearest neg SCLK cp (yr)")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(
        f"PWS-SCLK coincidence (p={res['T3_pws_coincidence']['p_closer_than_random']:.3f})")
fig.tight_layout()
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b112_pws_anchored_sclk.png", dpi=300)

logger.info("verdict: " + res["verdict"])
logger.data_save(out)
logger.data_save(RESULTS / "step_b112_pws_anchored_sclk.csv")
logger.data_save(FIG / "supplementary" / "step_b112_pws_anchored_sclk.png")
print("TEST SUMMARY:\n" + json.dumps(res["test_summary"], indent=1))
print(f"VERDICT: {res['verdict']}")
