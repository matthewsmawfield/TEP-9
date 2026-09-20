#!/usr/bin/env python3
"""
TEP-9 step 059 -- axis-consistency test
========================================

Two distinct TNO-derived directions appear in the pipeline:

  * the detached-sample axis (a>150, q>30, cc<=3):
      (lam, beta) = (49.9, -17.0) deg   [step_02]
  * the extreme-sample axis (a>250, q>30, cc<=3):
      (lam, beta) = (34.0, -12.6) deg   [step_02]

The comet-side cap tests (steps 030-041) were pre-declared on the
extreme axis (34,-13).  This step re-runs the headline cap contrast
and the continuous statistic under BOTH TNO axes, the CODE-recovered
axis, their antipodes, and non-TNO reference directions -- on the
independent CODE-only sample.  If the anomaly is a real structure in
the sky rather than a cap-placement artefact, it must survive under
either TNO axis and remain absent under the control directions.

Tests (CODE-only matched, N=54):
  C1  in-cap vs out-cap d_of contrast (MWU, greater) under each axis
  C2  continuous Spearman correlation d_of vs angular distance
  C3  same two tests on the Warsaw spike sample for cross-checking

Outputs: results/step_b24_axis_consistency.json,
         figures/supplementary/step_b24_axis_consistency.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_059_axis_consistency")
tee_stdout(logger)
logger.header("Axis-consistency test")

import json
import math
import re
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
from scipy.stats import mannwhitneyu, spearmanr
from html.parser import HTMLParser



def perih_dir(om, Om, inc):
    co, so = np.cos(om), np.sin(om)
    cO, sO, ci, si = np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci,
                     so * si])


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
    "tno_detached_50_-17": lv(49.9, -17.0),
    "tno_extreme_34_-13": lv(34.0, -13.0),
    "code_recovered_10_-20": lv(10.0, -20.0),
    "anti_detached": lv(229.9, 17.0),
    "anti_extreme": lv(214.0, 13.0),
    "ism_inflow": lv(255.8, 5.16),
    "gal_pole": gv(0, 90),
    "ecl_pole": np.array([0, 0, 1.0]),
}
GROUPS = {
    "tno_detached_50_-17": "TNO axis",
    "tno_extreme_34_-13": "TNO axis",
    "code_recovered_10_-20": "TNO axis",
    "anti_detached": "antipode",
    "anti_extreme": "antipode",
    "ism_inflow": "control",
    "gal_pole": "control",
    "ecl_pole": "control",
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
                desig=r[0].strip(), model=r[1].strip(),
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
                q=float(line[42:56]), e=float(line[56:70]),
                w=float(line[70:82]), Om=float(line[82:94]),
                i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows


def cap_test(sub, axis, capdeg=60):
    th = np.array([sep(r["aph"], axis) for r in sub])
    v = np.array([r["d_of"] for r in sub])
    inc = th < capdeg
    if inc.sum() < 4 or (~inc).sum() < 4:
        return None
    u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
    rho, pr = spearmanr(th, v)
    return {"n_in": int(inc.sum()), "n_out": int((~inc).sum()),
            "med_in": float(np.median(v[inc])),
            "med_out": float(np.median(v[~inc])),
            "cap_p": float(u.pvalue),
            "spearman_rho": float(rho), "spearman_p": float(pr)}


def main():
    # ---------- CODE-only matched sample ----------
    orig = parse_code(DATA_RAW / "code" / "code_original.html")
    osc = parse_code(DATA_RAW / "code" / "code_osculating.html")
    fut = parse_code(DATA_RAW / "code" / "code_future.html")
    warsaw = {l[5:17].strip()
              for l in open(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
              if len(l) > 115}
    rows = []
    for k, ro in orig.items():
        if k not in osc or k not in fut:
            continue
        po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                       math.radians(ro["i"]))
        pf = perih_dir(math.radians(fut[k]["w"]),
                       math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
        rows.append(dict(desig=k, warsaw=k in warsaw, aa=ro["aa"],
                         q=ro["q"], cls=ro["cls"], aph=-po,
                         d_of=sep(-po, -pf)))
    sp = [r for r in rows if 0 < r["aa"] < 100 and r["q"] < 3.1
          and r["cls"] in ("1a", "1a+", "1b")]
    co = [r for r in sp if not r["warsaw"]]
    print(f"CODE-only matched: {len(co)}")

    res = {"n_code_only": len(co), "cap_deg": 60,
           "code_only": {}}
    for name, ax in AXES.items():
        t = cap_test(co, ax)
        if t:
            t["group"] = GROUPS[name]
            res["code_only"][name] = t
            print(f"  {name:26s} cap p={t['cap_p']:.4f} "
                  f"rho={t['spearman_rho']:+.3f} p={t['spearman_p']:.4f}")

    # ---------- Warsaw spike cross-check ----------
    orows = parse_warsaw(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
    frows = parse_warsaw(DATA_RAW / "warsaw" / "warsaw_tabled.dat")
    PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
    PREF_FUT = {"i": 0, "l": 0, "j": 2, "k": 2}
    o_d, f_d = {}, {}
    for r in orows:
        k = r["desig"]
        if k not in o_d or PREF.get(r["com"], 9) < PREF.get(o_d[k]["com"], 9):
            o_d[k] = r
    for r in frows:
        k = r["desig"]
        if k not in f_d or PREF_FUT.get(r["com"], 9) < PREF_FUT.get(f_d[k]["com"], 9):
            f_d[k] = r
    wrows = []
    for k, ro in o_d.items():
        if k not in f_d:
            continue
        po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                       math.radians(ro["i"]))
        pf = perih_dir(math.radians(f_d[k]["w"]),
                       math.radians(f_d[k]["Om"]), math.radians(f_d[k]["i"]))
        wrows.append(dict(desig=k, aa=ro["aa"], aph=-po,
                          d_of=sep(-po, -pf)))
    wsp = [r for r in wrows if 0 < r["aa"] < 100]
    res["warsaw_spike"] = {"n": len(wsp)}
    for name in ("tno_detached_50_-17", "tno_extreme_34_-13",
                 "anti_extreme", "gal_pole"):
        t = cap_test(wsp, AXES[name])
        if t:
            res["warsaw_spike"][name] = t
            print(f"  W {name:26s} cap p={t['cap_p']:.4f} "
                  f"rho={t['spearman_rho']:+.3f}")

    # axis separations (for the record)
    res["axis_separations_deg"] = {
        "detached_vs_extreme": round(sep(AXES["tno_detached_50_-17"],
                                         AXES["tno_extreme_34_-13"]), 1),
        "detached_vs_code_recovered": round(
            sep(AXES["tno_detached_50_-17"],
                AXES["code_recovered_10_-20"]), 1),
        "extreme_vs_code_recovered": round(
            sep(AXES["tno_extreme_34_-13"],
                AXES["code_recovered_10_-20"]), 1)}

    (RESULTS / "step_b24_axis_consistency.json").write_text(
        json.dumps(res, indent=1, default=float))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4.6))
    names = list(res["code_only"].keys())
    ps = [res["code_only"][n]["cap_p"] for n in names]
    cols = {"TNO axis": "crimson", "antipode": "0.6",
            "control": "steelblue"}
    x = np.arange(len(names))
    ax.bar(x, [-math.log10(p) for p in ps],
           color=[cols[res["code_only"][n]["group"]] for n in names])
    ax.axhline(-math.log10(0.05), color="k", ls="--", lw=1,
               label="p=0.05")
    ax.set_xticks(x)
    ax.set_xticklabels([n.replace("_", "\n") for n in names],
                       fontsize=6.5)
    ax.set(ylabel="-log10(cap p)", title="CODE-only d_of cap contrast "
           "under candidate axes (N=54)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=g)
                       for g, c in cols.items()], fontsize=7)
    fig.tight_layout()
    fig.savefig(RESULTS / "figures" / "supplementary" / "step_b24_axis_consistency.png",
                dpi=300)
    logger.data_save(RESULTS / "step_b24_axis_consistency.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b24_axis_consistency.png")


if __name__ == "__main__":
    main()
