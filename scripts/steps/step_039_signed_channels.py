"""step_16: fresh signed channels on the independent CODE data.

Everything so far uses magnitudes (d_of).  TEP's boundary makes
SIGNED predictions never yet tested on this data:

  R1  repulsion sign: if the boundary deflects trajectories away
      from the axis, in-cap comets' future aphelia should move
      AWAY from it -- sign(theta_fut - theta_orig) > 0 in-cap,
      tested vs out-cap and vs a binomial null.  This is the
      comet-side analog of the anti-cap deficit and the b12
      JFC bipolar population.
  R2  rotation-axis coherence: the axis of each comet's
      orig->fut aphelion rotation.  A localized structure
      twists crossing trajectories about a shared direction;
      random planetary kicks do not.  Mean resultant length of
      in-cap rotation axes vs permutation.
  R3  doubly-matched control: each in-cap comet matched to the
      out-cap comet nearest in (q, |daa_of|) -- does d_of still
      differ when kick magnitude AND perihelion depth are
      held fixed?
  R4  joint patch fit: scan the axis on the COMBINED Warsaw+CODE
      independent data for the d-type statistic each catalogue
      supports, and report the joint best axis vs the TNO axis.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_039_signed_channels")
tee_stdout(logger)
logger.header("Signed channels and matched control")
import json, math, os, re
import numpy as np
from scipy.stats import mannwhitneyu, binomtest
from html.parser import HTMLParser

rng = np.random.default_rng(11)

def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO = lv(34, -13)

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

mat = []
for k, ro in orig.items():
    if k not in fut: continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    if 0 < ro["aa"] < 100 and ro["q"] < 3.1 and ro["cls"] in ("1a", "1a+", "1b") and k not in warsaw:
        a_o, a_f = -po, -pf
        mat.append(dict(desig=k, q=ro["q"], aph=a_o, aph_o=a_o, aph_f=a_f,
            theta_o=sep(a_o, TNO), theta_f=sep(a_f, TNO),
            d_of=sep(a_o, a_f), daa=abs(fut[k]["aa"] - ro["aa"])))

th = np.array([r["theta_o"] for r in mat]); v = np.array([r["d_of"] for r in mat])
inc = th < 60
res = {"n": len(mat), "n_in": int(inc.sum())}

# ---- R1 repulsion sign: theta_f - theta_o
dth = np.array([r["theta_f"] - r["theta_o"] for r in mat])
sgn_in = dth[inc]; sgn_out = dth[~inc]
res["R1_repulsion"] = {
    "in_cap_pos": int((sgn_in > 0).sum()), "in_cap_n": int(len(sgn_in)),
    "in_cap_med_dtheta": float(np.median(sgn_in)),
    "out_cap_pos": int((sgn_out > 0).sum()), "out_cap_n": int(len(sgn_out)),
    "out_cap_med_dtheta": float(np.median(sgn_out)),
    "p_binom_in": float(binomtest(int((sgn_in > 0).sum()), int(len(sgn_in)), 0.5,
                                 alternative="greater").pvalue),
    "p_binom_in_2sided": float(binomtest(int((sgn_in > 0).sum()), int(len(sgn_in)), 0.5).pvalue),
    "p_mw_in_vs_out": float(mannwhitneyu(sgn_in, sgn_out, alternative="greater").pvalue)}

# ---- R2 rotation-axis coherence in-cap
# rotation axis of each orig->fut rotation; sign the axis so it
# points roughly with the TNO-axis hemisphere for fold-invariance
rax = []
for r in mat:
    c = np.cross(r["aph_o"], r["aph_f"])
    n = np.linalg.norm(c)
    if n < 1e-12: rax.append(None); continue
    u = c / n
    if np.dot(u, TNO) < 0: u = -u   # fold: axis unoriented
    rax.append(u)
rax_i = np.array([rax[i] for i, x in enumerate(rax) if x is not None and inc[i]])
rax_o = np.array([rax[i] for i, x in enumerate(rax) if x is not None and not inc[i]])
def Rbar(A):
    return np.linalg.norm(A.mean(axis=0)) if len(A) else 0.0
Ri, Ro = Rbar(rax_i), Rbar(rax_o)
# permutation: shuffle in/out labels
cnt = 0; N = 20000
allr = np.array([x for x in rax if x is not None])
lab = np.array([inc[i] for i, x in enumerate(rax) if x is not None])
for _ in range(N):
    p = rng.permutation(lab)
    if Rbar(allr[p]) >= Ri: cnt += 1
res["R2_rotaxis"] = {
    "R_in": float(Ri), "R_out": float(Ro),
    "n_in_axes": int(len(rax_i)),
    "p_perm_in": float((cnt + 1) / (N + 1)),
    "mean_axis_sep_from_tno": float(sep(rax_i.mean(axis=0) / np.linalg.norm(rax_i.mean(axis=0)), TNO)) if len(rax_i) else None}

# ---- R3 doubly-matched control (q + |daa|)
qin = np.array([r["q"] for r in mat]); daa = np.array([r["daa"] for r in mat])
matched_pairs = []
out_idx = np.where(~inc)[0]
used = set()
for i in np.where(inc)[0]:
    cand = [j for j in out_idx if j not in used]
    if not cand: break
    j = min(cand, key=lambda j: abs(qin[j] - qin[i]) / 1.0 + abs(daa[j] - daa[i]) / 400.0)
    used.add(j); matched_pairs.append((i, j))
if len(matched_pairs) >= 8:
    din = np.array([v[i] for i, _ in matched_pairs])
    dout = np.array([v[j] for _, j in matched_pairs])
    u = mannwhitneyu(din, dout, alternative="greater")
    res["R3_matched"] = {"n_pairs": len(matched_pairs),
        "med_in": float(np.median(din)), "med_out": float(np.median(dout)),
        "p": float(u.pvalue),
        "med_dq": float(np.median(np.abs([qin[i]-qin[j] for i,j in matched_pairs]))),
        "med_ddaa": float(np.median(np.abs([daa[i]-daa[j] for i,j in matched_pairs])))}

# ---- R4 joint patch fit: Warsaw d_so + CODE d_of at scanned axis
# (each catalogue's own strongest leg, z-scored within catalogue)
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

wor_r = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
wosc_r = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tableb.dat"))
PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}
wor, wosc = {}, {}
for r in wor_r:
    k = r["desig"]
    if k not in wor or PREF.get(r["com"], 9) < PREF.get(wor[k]["com"], 9): wor[k] = r
for r in wosc_r:
    k = r["desig"]
    if k not in wosc or PREF_OSC.get(r["com"], 9) < PREF_OSC.get(wosc[k]["com"], 9): wosc[k] = r
wsub = []
for k, ro in wor.items():
    if k not in wosc or not (0 < ro["aa"] < 100): continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    ps = perih_dir(math.radians(wosc[k]["w"]), math.radians(wosc[k]["Om"]), math.radians(wosc[k]["i"]))
    wsub.append(dict(aph=-po, d=sep(-po, -ps)))

def scan(sub, dkey):
    best = (1.0, None, None)
    for lam in np.arange(0, 360, 10):
        for b in np.arange(-60, 61, 10):
            ax = lv(lam, b)
            t2 = np.array([sep(r["aph"], ax) for r in sub])
            v2 = np.array([r[dkey] for r in sub]); i2 = t2 < 60
            if i2.sum() < 4 or (~i2).sum() < 4: continue
            p2 = mannwhitneyu(v2[i2], v2[~i2], alternative="greater").pvalue
            if p2 < best[0]: best = (p2, lam, b)
    return best

bc = scan(mat, "d_of"); bw = scan(wsub, "d")
res["R4_joint_axes"] = {
    "code_best": {"p": float(bc[0]), "lam": bc[1], "beta": bc[2],
                  "sep_tno": float(sep(lv(bc[1], bc[2]), TNO))},
    "warsaw_best": {"p": float(bw[0]), "lam": bw[1], "beta": bw[2],
                    "sep_tno": float(sep(lv(bw[1], bw[2]), TNO))},
    "code_vs_warsaw_sep": float(sep(lv(bc[1], bc[2]), lv(bw[1], bw[2]))),
    "n_warsaw": len(wsub)}

out = str(RESULTS / "step_16_signed_channels.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)