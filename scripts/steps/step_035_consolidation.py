"""step_12: consolidate the CODE confirmation into its strongest form.

The step_11 result (CODE-only, matched regime): d_of excess at p=0.003,
axis-unique among 8, energy channel flat, rotation-per-kick doubled.

This step stress-tests and extends it on the SAME independent data:

  S1  cap-radius sensitivity -- a real boundary should peak at a
      characteristic radius and dilute beyond it, not grow with cap
  S2  NG-model mix -- do in-cap comets need non-gravitational solutions
      at a different rate than out-cap?  (unexplained momentum)
  S3  element decomposition of the d_of excess -- dOm, dw, di, de
      separately: a domain wall rotates the PLANE (Om, i) not the
      in-plane elements
  S4  era split -- pre-2006 vs post-2006 comets (survey-generation
      control)
  S5  anti-cap deficit -- repulsion signature: underpopulation of
      aphelia opposite the axis
  S6  shallow-plunger control -- q>3.1 comets (do NOT cross the inner
      planetary zone deeply) should show a weaker/null effect if the
      anomaly requires deep traversal
  S7  combined significance -- Fisher combination of the independent
      statistics (Warsaw d_so discovery, CODE d_of confirmation,
      CODE ratio_of), reported with the testing context
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_035_consolidation")
tee_stdout(logger)
logger.header("Confirmation consolidation")
import json, math, os, re
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
from scipy.stats import mannwhitneyu, fisher_exact, chi2
from html.parser import HTMLParser



def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO = lv(34, -13)
ANTI = -TNO

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

def parse(path):
    p = TP(); p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14: continue
        try:
            dates = re.findall(r"(\d{4})\s+(\d{2})\s+(\d{2})", r[5])
            out[r[0].strip()] = dict(desig=r[0].strip(), model=r[1].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                tyr=int(r[7].split()[0]) if r[7].split() else 0,
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out

orig = parse(str(DATA_RAW / "code" / "code_original.html"))
fut  = parse(str(DATA_RAW / "code" / "code_future.html"))
warsaw = {l[5:17].strip() for l in open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")) if len(l) > 115}

# CODE model nomenclature: the two-character solution code alone is
# ambiguous ('d5', 'bn', 'm1' are NG; 'b5', 'da' are GR).  The orbit
# class is carried by each row's abbr tooltip, which prefixes the
# class ('GR - gravitational orbit' vs 'NS/NT/NC/CT/NI - non-grav').
_MODEL_CLASS = {}
for _t, _c in re.findall(
        r"<abbr title=[\"']([^\"'>]+)[\"'][^>]*>([^<]+)</abbr>",
        open(str(DATA_RAW / "code" / "code_original.html"),
             encoding="utf-8", errors="replace").read()):
    _MODEL_CLASS[_c.strip()] = _t.strip().split(" ")[0]
def _is_ng(model):
    return _MODEL_CLASS.get(model, "GR") != "GR"

rows = []
for k, ro in orig.items():
    if k not in fut: continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    rows.append(dict(desig=k, warsaw=k in warsaw, aa=ro["aa"], q=ro["q"], i=ro["i"],
        cls=ro["cls"], model=ro["model"], tyr=ro["tyr"], aph=-po,
        d_of=sep(-po, -pf), daa_of=abs(fut[k]["aa"] - ro["aa"]),
        dOm=abs((fut[k]["Om"] - ro["Om"] + 180) % 360 - 180),
        dw=abs((fut[k]["w"] - ro["w"] + 180) % 360 - 180),
        di=abs(fut[k]["i"] - ro["i"]), de=abs(fut[k]["e"] - ro["e"]),
        theta=sep(-po, TNO), theta_anti=sep(-po, ANTI)))

sp  = [r for r in rows if 0 < r["aa"] < 100 and r["q"] < 3.1 and r["cls"] in ("1a", "1a+", "1b")]
co  = [r for r in sp if not r["warsaw"]]
sh  = [r for r in rows if 0 < r["aa"] < 100 and r["q"] >= 3.1 and r["cls"] in ("1a", "1a+", "1b") and not r["warsaw"]]
res = {"n_co": len(co), "n_shallow_control": len(sh)}

# S1 cap-radius sensitivity on CODE-only d_of
s1 = {}
for cap in (30, 40, 50, 60, 75, 90):
    t = np.array([r["theta"] for r in co]); v = np.array([r["d_of"] for r in co])
    inc = t < cap
    if inc.sum() >= 4 and (~inc).sum() >= 4:
        u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
        s1[str(cap)] = {"n_in": int(inc.sum()), "med_in": float(np.median(v[inc])),
                        "med_out": float(np.median(v[~inc])), "p": float(u.pvalue)}
res["S1_cap_sensitivity"] = s1

# S2 NG-model need in/out cap
th = np.array([r["theta"] for r in co])
ng = np.array([1 if _is_ng(r["model"]) else 0 for r in co])
inc60 = th < 60
tab = [[int(ng[inc60].sum()), int(inc60.sum() - ng[inc60].sum())],
       [int(ng[~inc60].sum()), int((~inc60).sum() - ng[~inc60].sum())]]
orr, pf = fisher_exact(tab)
res["S2_ng_need_60"] = {"in_cap_ng/tot": tab[0], "out_cap_ng/tot": tab[1],
                       "odds_ratio": float(orr), "p_2sided": float(pf)}

# S3 element decomposition of d_of excess
s3 = {}
for key in ("dOm", "dw", "di", "de"):
    v = np.array([r[key] for r in co])
    u = mannwhitneyu(v[inc60], v[~inc60], alternative="greater")
    s3[key] = {"med_in": float(np.median(v[inc60])), "med_out": float(np.median(v[~inc60])),
               "p": float(u.pvalue)}
res["S3_elements"] = s3

# S4 era split
s4 = {}
for lab, sub in [("pre2006", [r for r in co if r["tyr"] < 2006]),
                 ("post2006", [r for r in co if r["tyr"] >= 2006])]:
    if len(sub) >= 10:
        t = np.array([r["theta"] for r in sub]); v = np.array([r["d_of"] for r in sub])
        inc = t < 60
        if inc.sum() >= 3 and (~inc).sum() >= 3:
            u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
            s4[lab] = {"n": len(sub), "n_in": int(inc.sum()),
                       "med_in": float(np.median(v[inc])), "med_out": float(np.median(v[~inc])),
                       "p": float(u.pvalue)}
res["S4_era"] = s4

# S5 anti-cap deficit: fraction of aphelia within 60 of anti-axis
fa = np.array([r["theta_anti"] for r in co]) < 60
res["S5_anticap"] = {"n_in_anticap": int(fa.sum()), "n": len(co),
    "frac": float(fa.mean()), "expect_iso": float((1 - math.cos(math.radians(60))) / 2)}

# S6 shallow-plunger control (q>=3.1, CODE-only, class-1)
if len(sh) >= 10:
    t = np.array([r["theta"] for r in sh]); v = np.array([r["d_of"] for r in sh])
    inc = t < 60
    if inc.sum() >= 4 and (~inc).sum() >= 4:
        u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
        res["S6_shallow_control"] = {"n": len(sh), "n_in": int(inc.sum()),
            "med_in": float(np.median(v[inc])), "med_out": float(np.median(v[~inc])),
            "p": float(u.pvalue)}

# S8 all-quality-class stability: drop the class-1 restriction entirely
s8 = {}
allc = [r for r in rows if 0 < r["aa"] < 100 and r["q"] < 3.1]
for lab, sub in [("code_only", [r for r in allc if not r["warsaw"]]),
                 ("all_matched", allc)]:
    t = np.array([r["theta"] for r in sub]); v = np.array([r["d_of"] for r in sub])
    inc = t < 60
    if inc.sum() >= 4 and (~inc).sum() >= 4:
        u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
        s8[lab] = {"n": len(sub), "n_in": int(inc.sum()),
                   "med_in": float(np.median(v[inc])),
                   "med_out": float(np.median(v[~inc])),
                   "p": float(u.pvalue)}
res["S8_all_class"] = s8

# S7 combined significance (Fisher) across independent statistics
#   - Warsaw discovery: d_so full-arc p (step_09 axis_scan_ddir_full.tno)
#   - CODE-only: d_of p (step_11 code_only.d_of)
# CODE-only channels share data -> treat as one (min p, then a second
#   independent stat is not strictly independent; report conservative)
# Values are harvested from the live result files so the combination
# tracks the current pipeline outputs, not stale manuscript literals.
def _harvest(src, *path):
    d = json.loads((RESULTS / src).read_text())
    for k in path:
        d = d[k]
    return float(d)

p_warsaw = _harvest("step_09_nested_field.json",
                    "axis_scan_ddir_full", "tno", "p")
p_code = _harvest("step_11_code_confirm.json",
                  "code_only", "d_of", "p")
ps = [p_warsaw, p_code]
fstat = -2 * sum(math.log(p) for p in ps)
res["S7_combined"] = {
    "inputs": {"warsaw_d_so_fullarc": p_warsaw,
               "code_only_d_of": p_code},
    "fisher_chi2": float(fstat), "dof": 2 * len(ps),
    "p_combined": float(chi2.sf(fstat, 2 * len(ps))),
    "note": "conservative: Warsaw discovery + CODE-only confirmation are "
            "independent samples; CODE-internal secondary stats not added"}

out = str(RESULTS / "step_12_consolidate.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)