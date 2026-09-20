#!/usr/bin/env python3
"""
TEP-9 step 061 -- aphelion-sky channel: tide calibration and the
m=1 dipole
=====================================================================

The transit analyses (steps 030-060) test the RECONSTRUCTION channel:
how badly a standard-dynamics clock describes the inbound trajectory.
This step tests the complementary SPATIAL channel on the same comets:
where their original-orbit aphelia actually point on the sky.

Two measurements, both strictly real-data:

  T1  galactic-tide calibration.  The dominant known anisotropy of
      Oort-spike aphelia is the galactic tide, which is axisymmetric
      about the galactic pole and concentrates aphelia toward the
      galactic plane.  The pipeline's directional machinery must
      recover that known structure to be trusted on the unknown one.
      Test: distribution of |b_gal| of the aphelia vs isotropic
      (two-sample KS against a uniform-sphere draw, and the
      mean |b_gal|).

  T2  the m=1 dipole toward the boundary axis.  A tide cannot
      produce an m=1 asymmetry along a fixed ecliptic direction --
      it is axisymmetric about the galactic pole.  The measured
      dipole component along the pre-declared axis is therefore a
      clean observable:
          d_par = < aph_i . axis >
      calibrated against a TIDE-AWARE Monte-Carlo null that resamples
      galactic longitudes uniform while retaining the observed
      |b_gal| distribution (random hemisphere sign).  This is the
      spatial analogue of the clock-channel test: does the sky
      itself lean toward the axis?

Samples: Warsaw spike (orig aphelia, 0<1/a<100e-6), CODE class-1
matched, and the lpc one-apparition cohort -- three increasingly
independent draws of the same population.

Axes: both TNO axes plus the antipode and galactic-pole control.

Outputs: results/step_b26_aphelion_dipole.json,
         figures/supplementary/step_b26_aphelion_dipole.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_061_aphelion_dipole")
tee_stdout(logger)
logger.header("Aphelion-sky dipole and tide calibration")

import json
import math
import re
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
import xml.etree.ElementTree as ET
from scipy.stats import kstest
from html.parser import HTMLParser

RNG = np.random.default_rng(20260918)
N_MC = 20000



def perih_dir(om, Om, inc):
    co, so = np.cos(om), np.sin(om)
    cO, sO, ci, si = np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci, so * si])


def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))


def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l),
                     math.cos(b) * math.sin(l), math.sin(b)])


def gv(l, b):
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l), math.sin(b)])


AXES = {
    "tno_extreme_34_-13": lv(34.0, -13.0),
    "tno_detached_50_-17": lv(49.9, -17.0),
    "anti_extreme": lv(214.0, 13.0),
    "gal_pole": gv(0, 90),
    "gal_antictr_180_0": gv(180, 0),
    "gal_ctr_0_0": gv(0, 0),
}


class TP(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.cur, self.buf, self.in_td = [], [], "", False

    def handle_starttag(self, t, a):
        if t == "tr":
            self.cur = []
        elif t == "td":
            self.in_td, self.buf = True, ""

    def handle_endtag(self, t):
        if t == "td":
            self.in_td = False
            self.cur.append(self.buf.strip())
        elif t == "tr" and self.cur:
            self.rows.append(self.cur)

    def handle_data(self, d):
        if self.in_td:
            self.buf += d


def parse_code(path):
    p = TP()
    p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14:
            continue
        try:
            out[r[0].strip()] = dict(
                desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out


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
            d = dict(zip(fields, [td.text for td in list(el)]))
            out.append(d)
    return out


def gal_of(v):
    g = ECL2GAL @ v
    return (math.degrees(math.atan2(g[1], g[0])) % 360,
            math.degrees(math.asin(np.clip(g[2], -1, 1))))


def tide_aware_null(aphs, axis, n_mc=N_MC):
    """MC null: preserve observed |b_gal|; draw galactic longitude
    uniform and random hemisphere sign -> the tide's axisymmetric
    sky.  Return p for mean aph.axis >= observed."""
    n = len(aphs)
    obs = float(np.mean([np.dot(a, axis) for a in aphs]))
    bs = np.array([abs(gal_of(a)[1]) for a in aphs])
    cnt = 0
    for _ in range(n_mc):
        lrand = RNG.uniform(0, 2 * np.pi, n)
        bpick = bs[RNG.integers(0, n, n)]
        sgn = RNG.choice([-1.0, 1.0], n)
        bgal = np.deg2rad(bpick * sgn)
        vg = np.stack([np.cos(bgal) * np.cos(lrand),
                       np.cos(bgal) * np.sin(lrand),
                       np.sin(bgal)], axis=-1)
        ve = vg @ GAL2ECL.T  # rotate back to ecliptic
        m = float(np.mean(ve @ axis))
        if m >= obs:
            cnt += 1
    return obs, (cnt + 1) / (n_mc + 1)


def dipole_test(aphs, label):
    out = {"label": label, "n": len(aphs)}
    # T1 tide calibration: |b_gal| distribution vs isotropic
    bs = np.array([abs(gal_of(a)[1]) for a in aphs])
    u = RNG.uniform(-1, 1, 200000)
    iso_b = np.rad2deg(np.arcsin(np.abs(u)))
    ks = kstest(bs, iso_b)
    out["T1_tide_calibration"] = {
        "mean_abs_bgal": float(np.mean(bs)),
        "isotropic_mean": float(np.mean(iso_b)),
        "ks_p_vs_isotropic": float(ks.pvalue),
        "note": ("the galactic tide is axisymmetric about the gal "
                 "pole; the observed |b_gal| profile is the known "
                 "structure the machinery must recover")}
    # T2 dipole along candidate axes
    t2 = {}
    for name, ax in AXES.items():
        obs, p = tide_aware_null(aphs, ax)
        t2[name] = {"d_par": round(obs, 4), "p_tide_null": float(p)}
    out["T2_dipole_tide_null"] = t2
    print(f"  {label}: n={len(aphs)} mean|b|={np.mean(bs):.1f} "
          f"(iso {np.mean(iso_b):.1f}) KS p={ks.pvalue:.3f}")
    for name, t in t2.items():
        print(f"    {name:24s} d_par={t['d_par']:+.4f} p={t['p_tide_null']}")
    return out


def main():
    res = {}

    # Warsaw spike
    orows = parse_warsaw(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
    PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
    o_d = {}
    for r in orows:
        k = r["desig"]
        if k not in o_d or PREF.get(r["com"], 9) < PREF.get(o_d[k]["com"], 9):
            o_d[k] = r
    waph = [-perih_dir(math.radians(r["w"]), math.radians(r["Om"]),
                       math.radians(r["i"]))
            for r in o_d.values() if 0 < r["aa"] < 100]
    res["warsaw_spike"] = dipole_test(waph, "warsaw spike")

    # CODE class-1
    orig = parse_code(DATA_RAW / "code" / "code_original.html")
    warsaw = {l[5:17].strip()
              for l in open(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
              if len(l) > 115}
    caph = [-perih_dir(math.radians(r["w"]), math.radians(r["Om"]),
                       math.radians(r["i"]))
            for r in orig.values()
            if 0 < r["aa"] < 100 and r["q"] < 3.1
            and r["cls"] in ("1a", "1a+", "1b")]
    c_only = [-perih_dir(math.radians(r["w"]), math.radians(r["Om"]),
                         math.radians(r["i"]))
              for k, r in orig.items()
              if 0 < r["aa"] < 100 and r["q"] < 3.1
              and r["cls"] in ("1a", "1a+", "1b") and k not in warsaw]
    res["code_matched"] = dipole_test(caph, "code class-1 matched")
    res["code_only"] = dipole_test(c_only, "code-only")

    # lpc one-apparition cohort (original orbits)
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
    laph = [-perih_dir(math.radians(fv(r["arg"])),
                       math.radians(fv(r["long"])),
                       math.radians(fv(r["i"])))
            for r in seen.values()
            if 0 < fv(r["aa"]) < 100]
    res["lpc_one_apparition"] = dipole_test(laph, "lpc one-apparition")

    (RESULTS / "step_b26_aphelion_dipole.json").write_text(
        json.dumps(res, indent=1, default=float))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    bs = np.array([abs(gal_of(a)[1]) for a in waph])
    u = RNG.uniform(-1, 1, 200000)
    iso_b = np.rad2deg(np.arcsin(np.abs(u)))
    ax[0].hist(bs, bins=18, density=True, alpha=0.6, color="crimson",
               label="Warsaw spike aphelia")
    ax[0].hist(iso_b, bins=60, density=True, histtype="step",
               color="k", label="isotropic")
    ax[0].set(xlabel="|b_gal| of aphelion (deg)", ylabel="density",
              title="T1: tide calibration")
    ax[0].legend(fontsize=7)
    labels = list(res.keys())
    x = np.arange(len(labels))
    for j, axn in enumerate(("tno_extreme_34_-13", "anti_extreme")):
        vals = [res[l]["T2_dipole_tide_null"][axn]["d_par"]
                for l in labels]
        ax[1].bar(x + j * 0.35 - 0.18, vals, width=0.32,
                  color="crimson" if j == 0 else "0.6",
                  label=axn)
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xticks(x)
    ax[1].set_xticklabels([l.replace("_", "\n") for l in labels],
                          fontsize=7)
    ax[1].set(ylabel="d_par = <aph . axis>",
              title="T2: m=1 dipole along axis")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(RESULTS / "figures" / "supplementary" / "step_b26_aphelion_dipole.png",
                dpi=300)
    logger.data_save(RESULTS / "step_b26_aphelion_dipole.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b26_aphelion_dipole.png")


if __name__ == "__main__":
    main()
