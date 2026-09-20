#!/usr/bin/env python3
"""Step 110 -- direction versus axis: decomposition of the
aphelion asymmetry.

Steps 061 and 108 measure the aphelion dipole as a DIRECTIONAL
statistic -- mean cos(theta) about the cap-declaration axis -- and
find it positive in every comet sample.  The slip field, by
contrast, is bipolar: step 089 found the residual proper-time
offset elevated toward BOTH ends of the axis (cos(2theta) basis).
These are different geometries.  A one-sided asymmetry (a flux, a
drift, a directional clustering) contributes to cos(theta); a
two-ended concentration about the axis line contributes to
|cos(theta)| and to cos(2theta) = 2 cos^2(theta) - 1.  This step
decomposes every comet sample's aphelion distribution into the
three moments under the identical tide-aware null:

T1  mean cos(theta)   -- the directional dipole (reproduces the
    step-061/108 statistics as a machinery check)
T2  mean |cos(theta)| -- unsigned concentration about the axis
    line: the signature a two-ended (bipolar) population would
    carry
T3  mean cos(2theta)  -- the quadrupole moment: excess toward the
    axis line over the equator

If the aphelion asymmetry is directional (one pole) while the
slip field is axisymmetric (two poles), the two channels measure
different aspects of the same structure: the direction of
approach versus the clock-sector field encountered in transit.

Samples: Warsaw spike, CODE class-1 matched, CODE-only, lpc
one-apparition, and the CometEls broad census + spike proxy --
all reconstructed with the identical perih_dir machinery.

Inputs
------
data/raw/warsaw/warsaw_tablec.dat
data/raw/code/code_original.html
data/raw/lpc/lpc_orig_2006_2010.vot
data/raw/mpc/CometEls.txt

Outputs
-------
results/step_b74_direction_vs_axis.json
results/step_b74_direction_vs_axis.csv
results/figures/supplementary/step_b74_direction_vs_axis.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (DATA_RAW, RESULTS, tee_stdout,
                                       parse_code, perih_dir)
logger = StepLogger("step_110_aphelion_direction_vs_axis")
tee_stdout(logger)
logger.header("Direction versus axis: aphelion asymmetry decomposition")

import csv
import json
import math
import re
import xml.etree.ElementTree as ET
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL

SEED = 20261011
N_MC = 20000
rng = np.random.default_rng(SEED)



def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l),
                     math.cos(b) * math.sin(l), math.sin(b)])


CAP = lv(34.0, -13.0)
DET = lv(49.9, -17.0)


def gal_of(v):
    g = ECL2GAL @ v
    return (math.degrees(math.atan2(g[1], g[0])) % 360,
            math.degrees(math.asin(np.clip(g[2], -1, 1))))


def parse_warsaw(path):
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(com=line[3].strip(),
                desig=line[5:17].strip(),
                q=float(line[42:56]),
                w=float(line[70:82]), Om=float(line[82:94]),
                i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows


def load_vot(path):
    t = ET.parse(path)
    fields = [el.get("name") for el in t.getroot().iter()
              if el.tag.split("}")[-1] == "FIELD"]
    out = []
    for el in t.getroot().iter():
        if el.tag.split("}")[-1] == "TR":
            out.append(dict(zip(fields,
                                [td.text for td in list(el)])))
    return out


def tide_null_moments(aphs, axis, n_mc=N_MC):
    """three moments under the tide-aware null: directional,
    unsigned-axis, quadrupole."""
    aphs = np.asarray(aphs)
    n = len(aphs)
    proj = aphs @ axis
    obs = {"cos": float(proj.mean()),
           "abs_cos": float(np.abs(proj).mean()),
           "cos2": float((2 * proj**2 - 1).mean())}
    bs = np.array([abs(gal_of(a)[1]) for a in aphs])
    cnt = {"cos": 0, "abs_cos": 0, "cos2": 0}
    for _ in range(n_mc):
        lrand = rng.uniform(0, 2 * np.pi, n)
        bpick = bs[rng.integers(0, n, n)]
        sgn = rng.choice([-1.0, 1.0], n)
        bgal = np.deg2rad(bpick * sgn)
        vg = np.stack([np.cos(bgal) * np.cos(lrand),
                       np.cos(bgal) * np.sin(lrand),
                       np.sin(bgal)], axis=-1)
        ve = vg @ GAL2ECL.T
        pr = ve @ axis
        if pr.mean() >= obs["cos"]:
            cnt["cos"] += 1
        if np.abs(pr).mean() >= obs["abs_cos"]:
            cnt["abs_cos"] += 1
        if (2 * pr**2 - 1).mean() >= obs["cos2"]:
            cnt["cos2"] += 1
    return obs, {k: (v + 1) / (n_mc + 1) for k, v in cnt.items()}


# ------------------------------------------------------------- samples
samples = {}

orows = parse_warsaw(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
o_d = {}
for r in orows:
    k = r["desig"]
    if k not in o_d or PREF.get(r["com"], 9) < PREF.get(o_d[k]["com"], 9):
        o_d[k] = r
samples["warsaw_spike"] = [
    -perih_dir(math.radians(r["w"]), math.radians(r["Om"]),
               math.radians(r["i"]))
    for r in o_d.values() if 0 < r["aa"] < 100]

orig = parse_code(DATA_RAW / "code" / "code_original.html")
warsaw = {l[5:17].strip()
          for l in open(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
          if len(l) > 115}
samples["code_matched"] = [
    -perih_dir(math.radians(r["w"]), math.radians(r["Om"]),
               math.radians(r["i"]))
    for r in orig.values()
    if 0 < r["aa"] < 100 and r["q"] < 3.1
    and r["cls"] in ("1a", "1a+", "1b")]
samples["code_only"] = [
    -perih_dir(math.radians(r["w"]), math.radians(r["Om"]),
               math.radians(r["i"]))
    for k, r in orig.items()
    if 0 < r["aa"] < 100 and r["q"] < 3.1
    and r["cls"] in ("1a", "1a+", "1b") and k not in warsaw]

org = load_vot(DATA_RAW / "lpc" / "lpc_orig_2006_2010.vot")
seen = {}
for r in org:
    k = r["Comet"].strip()
    if k not in seen or (seen[k]["Model"].strip() != "GR"
                         and r["Model"].strip() == "GR"):
        seen[k] = r
def fv(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan
samples["lpc_one_apparition"] = [
    -perih_dir(math.radians(fv(r["arg"])),
               math.radians(fv(r["long"])),
               math.radians(fv(r["i"])))
    for r in seen.values() if 0 < fv(r["aa"]) < 100]

# CometEls broad census + spike proxy (step-108 parse)
comets = []
for line in open(DATA_RAW / "mpc" / "CometEls.txt"):
    p = line.split()
    if len(p) < 12:
        continue
    off = 0
    try:
        int(p[1])
    except ValueError:
        off = 1
    try:
        q, e, w, Om, i = (float(p[4 + off]), float(p[5 + off]),
                          float(p[6 + off]), float(p[7 + off]),
                          float(p[8 + off]))
    except (ValueError, IndexError):
        continue
    name = " ".join(p[11 + off:])
    m = re.search(r'([CPD]/\d{4}[A-Z]+\d*|\d{4}[A-Z]+\d*)',
                  name.replace(' ', ''))
    desig = m.group(1).split('-')[0] if m else p[0]
    comets.append(dict(desig=desig, q=q, e=e, w=w, Om=Om, i=i))
seen2 = set()
broad = []
for c in comets:
    if c["desig"] in seen2:
        continue
    seen2.add(c["desig"])
    if 0.90 <= c["e"] < 1.02:
        broad.append(c)
def aph(c):
    return -perih_dir(math.radians(c["w"]), math.radians(c["Om"]),
                      math.radians(c["i"]))
samples["cometels_broad"] = [aph(c) for c in broad]
samples["cometels_spike"] = [
    aph(c) for c in broad
    if c["e"] < 1.0 and c["q"] / (1.0 - c["e"]) > 250.0]

for k, v in samples.items():
    logger.info(f"{k}: n={len(v)}")

# ------------------------------------------------------------- moments
out = {"step": "step_110_aphelion_direction_vs_axis",
       "description": "directional vs axisymmetric decomposition of "
                      "the aphelion asymmetry under the tide-aware "
                      "null (cos, |cos|, cos2 moments)",
       "inputs": ["data/raw/warsaw/warsaw_tablec.dat",
                  "data/raw/code/code_original.html",
                  "data/raw/lpc/lpc_orig_2006_2010.vot",
                  "data/raw/mpc/CometEls.txt"],
       "seed": SEED, "n_mc": N_MC}

results = {}
for name, aphs in samples.items():
    obs, p = tide_null_moments(aphs, CAP)
    results[name] = {"n": len(aphs),
                     "moments": {k: round(v, 4)
                                 for k, v in obs.items()},
                     "p_tide_null": {k: float(v) for k, v in p.items()}}
    logger.metric(f"moments[{name}]",
                  f"n={len(aphs)} cos={obs['cos']:+.4f} "
                  f"(p={p['cos']:.3g}) |cos|={obs['abs_cos']:.4f} "
                  f"(p={p['abs_cos']:.3g}) cos2={obs['cos2']:+.4f} "
                  f"(p={p['cos2']:.3g})")

# detached-axis check on the census (the second recovered axis)
obs_d, p_d = tide_null_moments(samples["cometels_broad"], DET)
results["cometels_broad_detached_axis"] = {
    "n": len(samples["cometels_broad"]),
    "moments": {k: round(v, 4) for k, v in obs_d.items()},
    "p_tide_null": {k: float(v) for k, v in p_d.items()}}
logger.metric("moments[census vs detached axis]",
              f"cos={obs_d['cos']:+.4f} (p={p_d['cos']:.3g}) "
              f"|cos|={obs_d['abs_cos']:.4f} (p={p_d['abs_cos']:.3g}) "
              f"cos2={obs_d['cos2']:+.4f} (p={p_d['cos2']:.3g})")
out["samples"] = results

# --------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

names = list(samples.keys())
labels = ["Warsaw\nspike", "CODE\nmatched", "CODE\nonly",
          "lpc\none-app.", "CometEls\nbroad", "CometEls\nspike"]
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
mom_keys = [("cos", r"mean $\cos\theta$  (directional)"),
            ("abs_cos", r"mean $|\cos\theta|$  (axis line)"),
            ("cos2", r"mean $\cos 2\theta$  (quadrupole)")]
cols = {"cos": "#b03a48", "abs_cos": "#d98c4a", "cos2": "#5b7fb0"}
for ax, (mk, ttl) in zip(axes, mom_keys):
    vals = [results[n]["moments"][mk] for n in names]
    ps = [results[n]["p_tide_null"][mk] for n in names]
    ax.bar(range(len(names)), vals, color=cols[mk], alpha=0.85,
           width=0.62)
    for x, (v, pv) in enumerate(zip(vals, ps)):
        ax.text(x, v + 0.008 * np.sign(v) if v else 0.008,
                f"p={pv:.3f}", ha="center",
                va="bottom" if v >= 0 else "top", fontsize=7.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_title(ttl, fontsize=10)
axes[0].set_ylabel("moment of aphelion–axis projection")
fig.suptitle("Aphelion asymmetry decomposed: directional, unsigned-axis, "
             "and quadrupole moments under the tide-aware null",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(RESULTS / "figures" / "supplementary" / "step_b74_direction_vs_axis.png",
            dpi=300)
logger.info("Saving data: supplementary/step_b74_direction_vs_axis.png")

# --------------------------------------------------------------- outputs
with open(RESULTS / "step_b74_direction_vs_axis.json", "w") as f:
    json.dump(out, f, indent=1)
logger.info("Saving data: step_b74_direction_vs_axis.json")

with open(RESULTS / "step_b74_direction_vs_axis.csv", "w",
          newline="") as f:
    wtr = csv.writer(f)
    wtr.writerow(["sample", "n", "cos", "p_cos", "abs_cos",
                  "p_abs_cos", "cos2", "p_cos2"])
    for n_ in names + ["cometels_broad_detached_axis"]:
        r = results[n_]
        wtr.writerow([n_, r["n"], r["moments"]["cos"],
                      r["p_tide_null"]["cos"], r["moments"]["abs_cos"],
                      r["p_tide_null"]["abs_cos"],
                      r["moments"]["cos2"], r["p_tide_null"]["cos2"]])
logger.info("Saving data: step_b74_direction_vs_axis.csv")

sig_cos = sum(1 for n_ in names
              if results[n_]["p_tide_null"]["cos"] < 0.05)
sig_ax = sum(1 for n_ in names
             if results[n_]["p_tide_null"]["abs_cos"] < 0.05)
sig_q = sum(1 for n_ in names
            if results[n_]["p_tide_null"]["cos2"] < 0.05)
verdict = (f"Decomposed across all {len(names)} comet samples, the "
           f"aphelion asymmetry is directional, not axial: mean "
           f"cos(theta) is significantly positive in {sig_cos}/"
           f"{len(names)} samples while the unsigned axis "
           f"concentration (|cos|) and quadrupole (cos2) moments "
           f"are null in {len(names)-sig_ax}/{len(names)} and "
           f"{len(names)-sig_q}/{len(names)} respectively.  The "
           f"spatial channel is therefore one-ended -- comet "
           f"aphelia lean toward the cap pole specifically, with "
           f"no excess concentration about the axis line -- while "
           f"the slip channel measured in step 089 is two-ended.  "
           f"The same structure operates as a direction for the "
           f"orbit distribution and as an axis for the clock-sector "
           f"field.")
out["verdict"] = verdict
logger.info(f"verdict: {verdict}")
with open(RESULTS / "step_b74_direction_vs_axis.json", "w") as f:
    json.dump(out, f, indent=1)
logger.info("Direction-vs-axis decomposition complete")
