"""step_14: blind validation -- the tests that close the last loopholes.

The remaining critiques of the step_11-13 confirmation:
  (a) the axis seed came from the TNO side
  (b) the q<3.1/class-1 regime was inherited from Warsaw
  (c) the 60 deg cap was data-chosen

This battery removes each:

  B1  continuous statistic, no cap: Spearman rho between angular
      distance to the axis and d_of on ALL class-1 CODE-only
      Oort-spike comets (no q cut at all) -- and on the matched
      subset.  A real localized structure must decay monotonically
      with angle; the cap becomes irrelevant.
  B2  permutation null: shuffle aphelion directions among the same
      comets (preserves every marginal) -> empirical p for the
      cap contrast, no Gaussian assumptions.
  B3  look-elsewhere count: of the full-sky grid scan, how many
      independent directions reach p <= the observed?  The honest
      trials factor.
  B4  split-half blind validation: derive the best axis on a random
      half of the data, apply the cap-60 d_of test to the held-out
      half, repeat 500x.  If the structure is real, held-out p is
      small AND the recovered axes cluster at the TNO direction --
      axis and significance derived on data that never met.
  B5  reverse cross-validation: evaluate the Warsaw discovery
      statistic (d_so, full-arc) at the CODE-recovered axis
      (lam 10, beta -20) -- symmetric confirmation: each catalogue
      predicts the other's anomaly.
  B6  all-class-1 no-q-cut cap test: the matched-regime criterion
      dropped entirely -- does the signal survive losing the
      Warsaw-inherited cut?
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_037_blind_validation")
tee_stdout(logger)
logger.header("Blind validation battery")
import json, math, os, re
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
from html.parser import HTMLParser

rng = np.random.default_rng(42)

def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO = lv(34, -13); CODEAX = lv(10, -20)

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

def parse_code(path):
    p = TP(); p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14: continue
        try:
            out[r[0].strip()] = dict(desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out

orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))
fut  = parse_code(str(DATA_RAW / "code" / "code_future.html"))
warsaw = {l[5:17].strip() for l in open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")) if len(l) > 115}

allco, mat = [], []
for k, ro in orig.items():
    if k not in fut: continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    if 0 < ro["aa"] < 100 and ro["cls"] in ("1a", "1a+", "1b") and k not in warsaw:
        r = dict(desig=k, q=ro["q"], aph=-po, d_of=sep(-po, -pf))
        allco.append(r)
        if ro["q"] < 3.1: mat.append(r)
print(f"CODE-only class-1 spike: {len(allco)} total, {len(mat)} deep plungers (q<3.1)")

res = {"n_all_c1": len(allco), "n_matched": len(mat)}

# ---- B1 continuous statistic (no cap)
for lab, sub, ax in [("all_c1", allco, TNO), ("matched", mat, TNO), ("matched_codeaxis", mat, CODEAX)]:
    th = np.array([sep(r["aph"], ax) for r in sub]); v = np.array([r["d_of"] for r in sub])
    rho, p = spearmanr(th, v)
    res.setdefault("B1_spearman", {})[lab] = {"rho": float(rho), "p_2sided": float(p), "n": len(sub)}

# ---- B2 permutation null on matched CODE-only, cap 60
th = np.array([sep(r["aph"], TNO) for r in mat]); v = np.array([r["d_of"] for r in mat])
inc = th < 60
obs = np.median(v[inc]) - np.median(v[~inc])
cnt = 0
N = 20000
for _ in range(N):
    vs = rng.permutation(v)
    if np.median(vs[inc]) - np.median(vs[~inc]) >= obs: cnt += 1
res["B2_perm"] = {"obs_meddiff": float(obs), "p_perm": float((cnt + 1) / (N + 1)), "n_perm": N}

# ---- B3 look-elsewhere: count grid axes beating the observed p
hits = []
for lam in np.arange(0, 360, 5):
    for b in np.arange(-75, 76, 5):
        ax = lv(lam, b)
        t2 = np.array([sep(r["aph"], ax) for r in mat])
        i2 = t2 < 60
        if i2.sum() < 4 or (~i2).sum() < 4: continue
        p2 = mannwhitneyu(v[i2], v[~i2], alternative="greater").pvalue
        hits.append((p2, lam, b))
hits.sort()
n_grid = len(hits); n_better = sum(1 for p2, _, _ in hits if p2 <= 0.0029505468290258216)
res["B3_look_elsewhere"] = {"n_grid": n_grid, "n_axes_p_le_obs": n_better,
    "fraction": n_better / n_grid, "best": hits[:5],
    "note": "nearby grid points are correlated; the effective number of "
            "independent sky directions is ~n_grid/(cap solid-angle factor) "
            "-- the raw fraction is the honest upper bound on trials"}

# ---- B4 split-half blind validation
loo_ax = []; hp = []
for rep in range(500):
    idx = rng.permutation(len(mat)); A = idx[:27]; B = idx[27:]
    subA = [mat[i] for i in A]; subB = [mat[i] for i in B]
    best = (1.0, None, None)
    for lam in np.arange(0, 360, 15):
        for b in np.arange(-60, 61, 15):
            ax = lv(lam, b)
            t2 = np.array([sep(r["aph"], ax) for r in subA])
            v2 = np.array([r["d_of"] for r in subA]); i2 = t2 < 60
            if i2.sum() < 3 or (~i2).sum() < 3: continue
            p2 = mannwhitneyu(v2[i2], v2[~i2], alternative="greater").pvalue
            if p2 < best[0]: best = (p2, lam, b)
    ax = lv(best[1], best[2])
    loo_ax.append(sep(ax, TNO))
    t2 = np.array([sep(r["aph"], ax) for r in subB]); v2 = np.array([r["d_of"] for r in subB])
    i2 = t2 < 60
    hp.append(mannwhitneyu(v2[i2], v2[~i2], alternative="greater").pvalue
              if i2.sum() >= 3 and (~i2).sum() >= 3 else 1.0)
hp = np.array(hp)
res["B4_split_half"] = {
    "heldout_p_median": float(np.median(hp)),
    "heldout_p_lt_0.05_frac": float((hp < 0.05).mean()),
    "heldout_p_lt_0.10_frac": float((hp < 0.10).mean()),
    "recovered_axis_to_tno_med_deg": float(np.median(loo_ax)),
    "recovered_within_45deg_frac": float((np.array(loo_ax) < 45).mean()),
    "note": "axis derived on half A, tested on held-out half B, 500 reps"}

# ---- B5 reverse: Warsaw d_so at the CODE-recovered axis
def parse_orbit_table(path):
    rows = []
    for line in open(path):
        if len(line) < 115: continue
        try:
            rows.append(dict(sample=line[0:2].strip(), com=line[3].strip(),
                desig=line[5:17].strip(), q=float(line[42:56]),
                w=float(line[70:82]), Om=float(line[82:94]),
                i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows

def parse_a1(path):
    out = {}
    for line in open(path):
        if len(line) < 160: continue
        d = line[3:15].strip()
        if d: out.setdefault(d, []).append(dict(datat=line[122:132].strip()))
    return out

PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
def dedup(rows):
    out = {}
    for r in rows:
        k = r["desig"]
        if k not in out or PREF.get(r["com"], 9) < PREF.get(out[k]["com"], 9): out[k] = r
    return out

wor = dedup(parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")))
wosc_rows = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tableb.dat"))
PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}
wosc = {}
for r in wosc_rows:
    k = r["desig"]
    if k not in wosc or PREF_OSC.get(r["com"], 9) < PREF_OSC.get(wosc[k]["com"], 9): wosc[k] = r
meta = parse_a1(str(DATA_RAW / "warsaw" / "warsaw_tablea1.dat"))

wsub = []
for k, ro in wor.items():
    if k not in wosc: continue
    if not (0 < ro["aa"] < 100): continue
    if "full" not in str(meta.get(k, [{}])[0].get("datat", "")): continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    ps = perih_dir(math.radians(wosc[k]["w"]), math.radians(wosc[k]["Om"]), math.radians(wosc[k]["i"]))
    wsub.append(dict(desig=k, aph=-po, d_so=sep(-po, -ps)))

for lab, ax in [("tno", TNO), ("code_recovered", CODEAX)]:
    t2 = np.array([sep(r["aph"], ax) for r in wsub]); v2 = np.array([r["d_so"] for r in wsub])
    i2 = t2 < 60
    if i2.sum() >= 4 and (~i2).sum() >= 4:
        u = mannwhitneyu(v2[i2], v2[~i2], alternative="greater")
        res.setdefault("B5_reverse_warsaw", {})[lab] = {
            "n": len(wsub), "n_in": int(i2.sum()),
            "med_in": float(np.median(v2[i2])), "med_out": float(np.median(v2[~i2])),
            "p": float(u.pvalue)}

# ---- B6 no-q-cut, all class-1 CODE-only, cap 60 at TNO axis
t2 = np.array([sep(r["aph"], TNO) for r in allco]); v2 = np.array([r["d_of"] for r in allco])
i2 = t2 < 60
u = mannwhitneyu(v2[i2], v2[~i2], alternative="greater")
res["B6_no_qcut"] = {"n": len(allco), "n_in": int(i2.sum()),
    "med_in": float(np.median(v2[i2])), "med_out": float(np.median(v2[~i2])), "p": float(u.pvalue)}

out = str(RESULTS / "step_14_blind_validation.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)