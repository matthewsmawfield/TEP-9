#!/usr/bin/env python3
"""Step 075 -- signed direction of the slip: lag or advance?

A proper-time lapse is not a random torque: a comet that loses proper
time during the inside passage arrives at each observed position late
in orbital phase, and a standard-clock reconstruction absorbs the
delay as a systematic rotation of the fitted orbit -- always in the
same signed sense relative to the comet's own direction of motion.
A stray gravitational perturbation, by contrast, rotates the orbit
either way depending on encounter geometry.

For every comet the signed rotation from the original-leg to the
future-leg periapsis direction is measured about the comet's own
orbit pole:  positive = periapsis advanced along the direction of
orbital motion.  CODE uses the orig->fut pair (its anomaly channel);
Warsaw uses osc->orig (the catalogue's inbound channel, where its
anomaly lives).  Sign coherence is tested in and out of the cap,
pooled and per cohort, and split by prograde/retrograde orientation
to separate a universal phase sign from a sky-fixed direction kick.

Outputs:
    results/step_b40_signed_slip.json
    results/step_b40_signed_slip.csv
    results/figures/supplementary/step_b40_signed_slip.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout, parse_code
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_075_signed_slip")
tee_stdout(logger)
logger.header("Signed slip direction: lag vs advance")

import csv
import json
import math
import re
import numpy as np
from scipy.stats import binomtest, spearmanr
from html.parser import HTMLParser

CAP  = 60.0
SEED = 20260918
rng  = np.random.default_rng(SEED)

def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def pole(Om, inc):
    sO, cO, si, ci = math.sin(Om), math.cos(Om), math.sin(inc), math.cos(inc)
    return np.array([si*sO, -si*cO, ci])

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

TNO = lv(34, -13)

def signed_rot(p_from, p_to, h):
    """signed angle p_from -> p_to about pole h (deg); + = along motion."""
    x = np.cross(p_from, p_to)
    s = np.dot(x, h)
    c = np.dot(p_from, p_to)
    return math.degrees(math.atan2(s, c))

# ------------------------------------------------------------------
# CODE: full parser for orig + fut (three-leg signed rotation)
# ------------------------------------------------------------------

class TP(HTMLParser):
    def __init__(self):
        super().__init__(); self.rows = []; self.cur = []; self.buf = ""; self.in_td = False
    def handle_starttag(self, t, a):
        if t == "tr": self.cur = []
        elif t == "td": self.in_td = True; self.buf = ""
    def handle_endtag(self, t):
        if t == "td": self.in_td = False; self.cur.append(self.buf.strip())
        elif t == "tr" and self.cur: self.rows.append(self.cur)
    def handle_data(self, d):
        if self.in_td: self.buf += d

def parse_code_full(path):
    p = TP(); p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14: continue
        try:
            out[r[0].strip()] = dict(desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                T=r[7].strip(),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out

orig = parse_code_full(str(DATA_RAW / "code" / "code_original.html"))
fut  = parse_code_full(str(DATA_RAW / "code" / "code_future.html"))
osc  = parse_code_full(str(DATA_RAW / "code" / "code_osculating.html"))
warsaw_ids = {l[5:17].strip() for l in
              open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
              if len(l) > 115}

code_rows = []
for k, ro in orig.items():
    if k not in fut: continue
    if not (0 < ro["aa"] < 100 and ro["cls"] in ("1a", "1a+", "1b")) or k in warsaw_ids:
        continue
    if k not in osc: continue
    w_, O_, i_ = map(math.radians, (ro["w"], ro["Om"], ro["i"]))
    po = perih_dir(w_, O_, i_)
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]),
                   math.radians(fut[k]["i"]))
    h = pole(O_, i_)
    code_rows.append(dict(desig=k, cohort="code", q=ro["q"],
        theta=sep(-po, TNO),
        srot=signed_rot(po, pf, h),
        prograde=bool(h[2] > 0)))

# ------------------------------------------------------------------
# Warsaw: three legs live in tableb (osculating), tablec (original),
# tabled (future); signed osc -> orig rotation is the inbound channel
# ------------------------------------------------------------------

PREF     = {"a": 0, "h": 0, "e": 1, "b": 2}
PREF_FUT = {"i": 0, "l": 0, "j": 2, "k": 2}
PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}

def parse_orbit_table(path):
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(desig=line[5:17].strip(),
                             com=line[3].strip(),
                             q=float(line[42:56]),
                             w=float(line[70:82]), Om=float(line[82:94]),
                             i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows

def dedup(rows, pref):
    out = {}
    for r in rows:
        k = r["desig"]
        if k not in out or pref.get(r["com"], 9) < pref.get(out[k]["com"], 9):
            out[k] = r
    return out

war_osc  = dedup(parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tableb.dat")), PREF_OSC)
war_orig = dedup(parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")), PREF)
war_fut  = dedup(parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tabled.dat")), PREF_FUT)

war_rows = []
for d, b_leg in war_orig.items():
    a_leg = war_osc.get(d)
    if not a_leg or d not in war_fut:
        continue
    if not (0 < b_leg["aa"] < 100):
        continue
    w_, O_, i_ = map(math.radians, (b_leg["w"], b_leg["Om"], b_leg["i"]))
    po = perih_dir(w_, O_, i_)
    pa = perih_dir(math.radians(a_leg["w"]), math.radians(a_leg["Om"]),
                   math.radians(a_leg["i"]))
    h = pole(O_, i_)
    war_rows.append(dict(desig=d, cohort="warsaw", q=b_leg["q"],
        theta=sep(-po, TNO),
        srot=signed_rot(pa, po, h),          # osc -> orig (inbound leg)
        prograde=bool(h[2] > 0)))

rows = code_rows + war_rows
logger.info(f"signed rotations: code {len(code_rows)}, warsaw {len(war_rows)}")

# ------------------------------------------------------------------
# Sign statistics
# ------------------------------------------------------------------

def sign_stats(sub):
    th = np.array([r["theta"] for r in sub]); inc = th < CAP
    sr = np.array([r["srot"] for r in sub])
    out = {"n": int(len(sub)), "n_in": int(inc.sum())}
    for tag, m in (("in", inc), ("out", ~inc)):
        s = sr[m]
        if len(s) < 5:
            out[tag] = None; continue
        npos = int((s > 0).sum())
        bt = binomtest(npos, len(s), 0.5)
        out[tag] = {"n": int(len(s)), "n_pos": npos,
                    "frac_pos": npos / len(s),
                    "median_signed_deg": float(np.median(s)),
                    "median_abs_deg": float(np.median(np.abs(s))),
                    "p_binom": float(bt.pvalue)}
    return out

res = {"method": "signed periapsis rotation about each comet's own "
                 "orbit pole; + = advanced along direction of motion.  "
                 "CODE uses orig->fut (its anomaly channel), Warsaw "
                 "osc->orig (its inbound channel).  A lapse produces "
                 "one signed sense relative to motion; a random torque "
                 "does not.",
       "cap_deg": CAP, "axis": "34,-13", "seed": SEED,
       "sign_convention": "+ = periapsis advanced along orbital motion"}

for co in ("code", "warsaw", "pooled"):
    sub = [r for r in rows if co == "pooled" or r["cohort"] == co]
    res[f"sign_{co}"] = sign_stats(sub)
    for tag in ("in", "out"):
        s = res[f"sign_{co}"].get(tag)
        if s:
            logger.info(f"{co}/{tag}: n={s['n']} pos={s['n_pos']} "
                        f"({s['frac_pos']:.2f}) p={s['p_binom']:.4f} "
                        f"med={s['median_signed_deg']:+.4f} deg")

# orientation split on the pooled in-cap set: prograde vs retrograde
inc_rows = [r for r in rows if r["theta"] < CAP]
for lab, pr in (("prograde", True), ("retrograde", False)):
    s = np.array([r["srot"] for r in inc_rows if r["prograde"] == pr])
    if len(s) >= 5:
        npos = int((s > 0).sum())
        bt = binomtest(npos, len(s), 0.5)
        res[f"incap_{lab}"] = {"n": int(len(s)), "n_pos": npos,
                               "frac_pos": npos / len(s),
                               "median_signed_deg": float(np.median(s)),
                               "p_binom": float(bt.pvalue)}
        logger.info(f"in-cap {lab}: n={len(s)} pos={npos} "
                    f"med={np.median(s):+.4f} p={bt.pvalue:.4f}")

# ------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------

with open(RESULTS / "step_b40_signed_slip.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["desig", "cohort", "q", "theta",
                                      "srot", "prograde"])
    w.writeheader()
    for r in rows: w.writerow(r)

with open(RESULTS / "step_b40_signed_slip.json", "w") as f:
    json.dump(res, f, indent=1)
logger.data_save(RESULTS / 'step_b40_signed_slip.json')
logger.data_save(RESULTS / 'step_b40_signed_slip.csv')
# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7.2, 4.4))
bins = np.linspace(-1.5, 1.5, 60)
for co, mkm in (("code", "code"), ("warsaw", "warsaw")):
    sub = [r for r in rows if r["cohort"] == co]
    sr = np.array([r["srot"] for r in sub])
    th = np.array([r["theta"] for r in sub])
    ax.hist(sr[th >= CAP], bins=bins, histtype="step", color="0.55",
            lw=1.2, label=f"{co} outside")
    ax.hist(sr[th < CAP], bins=bins, histtype="step", color="crimson",
            lw=1.5, label=f"{co} inside")
ax.axvline(0, color="k", ls=":", lw=1)
ax.set_xlabel(r"signed rotation about orbit pole (deg; $+$ = along motion)")
ax.set_ylabel("count")
ax.legend(frameon=False, fontsize=8)
ax.set_title("slip direction: coherent sign = phase lag/advance, "
             "mixed = random torque", fontsize=9)
fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b40_signed_slip.png", dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b40_signed_slip.png')