#!/usr/bin/env python3
"""Step 079 -- cross-fitter element audit (lineage noise floor).

The deepest structural caveat is shared orbit-determination lineage:
Warsaw, CODE, and MW08 share fitting machinery, so the transit
replication is independent in sample rather than in method.  MPC's
CometEls is a fourth, separately maintained fit record.  For every
comet appearing in both CometEls and CODE osculating (or Warsaw
osculating), the periapsis-direction difference between the two
fitters is the direct, per-object lineage noise floor.  If the in-cap
orig->fut anomaly (~0.14 deg) sits far above that floor, and the floor
itself is cap-symmetric, a shared fitter convention cannot be the
explanation.

Inputs : data/raw/mpc/CometEls.txt,
         data/raw/code/code_osculating.html,
         data/raw/warsaw/warsaw_tableb.dat,
         results/step_b30_proper_time_slip.csv
Outputs: results/step_b44_cross_fitter.json
         results/step_b44_cross_fitter.csv
         results/figures/supplementary/step_b44_cross_fitter.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_079_cross_fitter")
tee_stdout(logger)
logger.header("Cross-fitter element audit -- lineage noise floor")

import json
import math
import re
import csv
import numpy as np
from scipy.stats import mannwhitneyu
from html.parser import HTMLParser

SEED = 20260918

def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

AXIS = lv(34.0, -13.0)   # transit-axis convention (comet cap, steps 030-086)
CAP  = 60.0

# ------------------------------------------------------------------
# Parsers
# ------------------------------------------------------------------

class TP(HTMLParser):
    def __init__(self):
        super().__init__(); self.rows=[]; self.cur=[]; self.buf=""; self.in_td=False
    def handle_starttag(self,t,a):
        if t=="tr": self.cur=[]
        elif t=="td": self.in_td=True; self.buf=""
    def handle_endtag(self,t):
        if t=="td": self.in_td=False; self.cur.append(self.buf.strip())
        elif t=="tr" and self.cur: self.rows.append(self.cur)
    def handle_data(self,d):
        if self.in_td: self.buf+=d

def parse_code_full(path):
    p=TP(); p.feed(open(path,encoding="utf-8",errors="replace").read())
    out={}
    for r in p.rows:
        if len(r)<14: continue
        try:
            out[r[0].strip()]=dict(desig=r[0].strip(),
                cls=re.sub(r"^\d","",r[3].strip()), T=r[7].strip(),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError,IndexError): continue
    return out

def parse_els(path):
    out={}
    for line in open(path):
        t=line.split()
        if len(t)<13: continue
        nm=re.sub(r"\s*\(.*","", " ".join(t[12:-1])).strip()
        try:
            out[nm]=dict(desig=nm, T=f"{t[1]}-{t[2]}-{t[3]}",
                q=float(t[4]), e=float(t[5]), w=float(t[6]),
                Om=float(t[7]), i=float(t[8]))
        except ValueError: continue
    return out

def parse_orbit_table(path):
    rows=[]
    for line in open(path):
        if len(line)<115: continue
        try:
            rows.append(dict(desig=line[5:17].strip(), com=line[3].strip(),
                q=float(line[42:56]), w=float(line[70:82]),
                Om=float(line[82:94]), i=float(line[94:106]),
                aa=float(line[106:115])))
        except ValueError: continue
    return rows

PREF_OSC={"a":0,"g":0,"d":1,"e":2,"f":2,"b":3,"c":3}
def dedup(rows, pref):
    out={}
    for r in rows:
        k=r["desig"]
        if k not in out or pref.get(r["com"],9)<pref.get(out[k]["com"],9):
            out[k]=r
    return out

code_osc  = parse_code_full(str(DATA_RAW/"code"/"code_osculating.html"))
code_orig = parse_code_full(str(DATA_RAW/"code"/"code_original.html"))
els       = parse_els(str(DATA_RAW/"mpc"/"CometEls.txt"))
war_osc   = dedup(parse_orbit_table(str(DATA_RAW/"warsaw"/"warsaw_tableb.dat")), PREF_OSC)
PREF      = {"a":0,"h":0,"e":1,"b":2}
war_orig  = dedup(parse_orbit_table(str(DATA_RAW/"warsaw"/"warsaw_tablec.dat")), PREF)

logger.info(f"fits: code_osc {len(code_osc)}, CometEls {len(els)}, warsaw_osc {len(war_osc)}")

# cap membership from step_b30 per-comet file (CODE class-1)
b30={}
with open(RESULTS/"step_b30_proper_time_slip.csv") as f:
    for r in csv.DictReader(f):
        b30[r["desig"]]=r

def cap_of(w,Om,i):
    d=-perih_dir(math.radians(w),math.radians(Om),math.radians(i))
    return float(np.degrees(np.arccos(np.clip(d@AXIS,-1,1))))

rows=[]
for k,co in code_osc.items():
    if k in els and k in b30:          # class-1 CODE sample
        rows.append(dict(desig=k, cat="code",
            w=co["w"],Om=co["Om"],i=co["i"],
            w2=els[k]["w"],Om2=els[k]["Om"],i2=els[k]["i"],
            theta=float(b30[k]["theta"])))
war_matched=0
for k,wo in war_osc.items():
    o=war_orig.get(k)
    if k in els and o is not None and 0<o["aa"]<100:
        rows.append(dict(desig=k, cat="warsaw",
            w=wo["w"],Om=wo["Om"],i=wo["i"],
            w2=els[k]["w"],Om2=els[k]["Om"],i2=els[k]["i"],
            theta=cap_of(o["w"],o["Om"],o["i"])))
        war_matched+=1
logger.info(f"matched: code {sum(r['cat']=='code' for r in rows)}, warsaw {war_matched}")

for r in rows:
    d1= perih_dir(math.radians(r["w"]),math.radians(r["Om"]),math.radians(r["i"]))
    d2= perih_dir(math.radians(r["w2"]),math.radians(r["Om2"]),math.radians(r["i2"]))
    r["dperi_deg"]=float(np.degrees(np.arccos(np.clip(d1@d2,-1,1))))
    r["in_cap"]=r["theta"]<CAP

# ------------------------------------------------------------------
# Statistics
# ------------------------------------------------------------------

res={"inputs":["data/raw/code/code_osculating.html",
               "data/raw/code/code_original.html",
               "data/raw/mpc/CometEls.txt",
               "data/raw/warsaw/warsaw_tableb.dat",
               "data/raw/warsaw/warsaw_tablec.dat",
               "results/step_b30_proper_time_slip.csv"],
     "method":"per-comet periapsis-direction difference between the "
       "MPC CometEls fit and the catalogue osculating fit (CODE or "
       "Warsaw) -- the direct inter-fitter noise floor; compared in "
       "and out of the 60-deg cap and against the measured anomaly.",
     "cap_deg":CAP,"axis":"34,-13","seed":SEED}

out={}
for cat in ("code","warsaw","pooled"):
    sub=[r for r in rows if cat=="pooled" or r["cat"]==cat]
    if not sub: continue
    inn=[r["dperi_deg"] for r in sub if r["in_cap"]]
    outt=[r["dperi_deg"] for r in sub if not r["in_cap"]]
    e={"n":len(sub),"n_in":len(inn),"n_out":len(outt),
       "median_dperi_deg":float(np.median([r["dperi_deg"] for r in sub])),
       "median_in":float(np.median(inn)) if inn else None,
       "median_out":float(np.median(outt)) if outt else None}
    if len(inn)>3 and len(outt)>3:
        u=mannwhitneyu(inn,outt)
        e["p_mwu_in_vs_out"]=float(u.pvalue)
        logger.info(f"{cat}: n={len(sub)} floor={e['median_dperi_deg']:.4f} deg "
                    f"in={e['median_in']:.4f} out={e['median_out']:.4f} "
                    f"p={u.pvalue:.4f}")
    out[cat]=e
res["groups"]=out

# does the per-comet floor track the per-comet anomaly?  A shared
# poorly-determined subset would correlate the two channels.
from scipy.stats import spearmanr
code_rows=[r for r in rows if r["cat"]=="code" and r["desig"] in b30]
if len(code_rows)>10:
    x=np.array([r["dperi_deg"] for r in code_rows])
    y=np.array([abs(float(b30[r["desig"]]["drot"])) for r in code_rows])
    for tag,m in (("in_cap",[r["in_cap"] for r in code_rows]),
                  ("out_cap",[not r["in_cap"] for r in code_rows]),
                  ("all",[True]*len(code_rows))):
        mm=np.array(m)
        rho,p=spearmanr(x[mm],y[mm])
        res.setdefault("floor_vs_anomaly",{})[tag]={
            "n":int(mm.sum()),"rho":float(rho),"p_2sided":float(p)}
        logger.info(f"floor vs |drot| {tag}: n={int(mm.sum())} rho={rho:+.2f} p={p:.3f}")

# direction of the inter-fitter displacement: a direction-dependent
# fitter bias would make the cross-fitter difference vectors point
# coherently inside the cap (and only there).
def disp_dir(r):
    d1=perih_dir(math.radians(r["w"]),math.radians(r["Om"]),math.radians(r["i"]))
    d2=perih_dir(math.radians(r["w2"]),math.radians(r["Om2"]),math.radians(r["i2"]))
    x=np.cross(d1,d2)
    n=np.linalg.norm(x)
    return x/n if n>1e-12 else None
dirs_in=[];dirs_out=[]
for r in rows:
    v=disp_dir(r)
    if v is None: continue
    (dirs_in if r["in_cap"] else dirs_out).append(v)
def rayleigh3(vs):
    from scripts.utils.statistics import rayleigh_3d
    return rayleigh_3d(vs)
if len(dirs_in)>3 and len(dirs_out)>3:
    Ri,pi=rayleigh3(dirs_in); Ro,po=rayleigh3(dirs_out)
    res["disp_dir"]={"n_in":len(dirs_in),"n_out":len(dirs_out),
        "R_in":Ri,"p_in":pi,"R_out":Ro,"p_out":po}
    logger.info(f"displacement-direction coherence: in R={Ri:.2f} p={pi:.3f}; "
                f"out R={Ro:.2f} p={po:.3f}")

with open(RESULTS/"step_b44_cross_fitter.json","w") as f:
    json.dump(res,f,indent=1)
with open(RESULTS/"step_b44_cross_fitter.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["desig","cat","theta","in_cap","dperi_deg",
                                   "w","Om","i","w2","Om2","i2"])
    w.writeheader()
    for r in rows: w.writerow({k:r[k] for k in w.fieldnames})
logger.data_save(RESULTS/'step_b44_cross_fitter.json')
logger.data_save(RESULTS/'step_b44_cross_fitter.csv')
# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig,axes=plt.subplots(1,2,figsize=(9,4.2))
ax=axes[0]
for cat,col in (("code","crimson"),("warsaw","darkorange")):
    vals=sorted(r["dperi_deg"] for r in rows if r["cat"]==cat)
    if vals:
        ax.plot(vals,np.linspace(0,1,len(vals)),color=col,lw=1.4,label=cat)
ax.axvline(0.14,color="k",ls=":",lw=1,label="in-cap anomaly")
ax.set_xscale("log"); ax.set_xlabel("|delta periapsis dir| between fitters (deg)")
ax.set_ylabel("CDF"); ax.legend(frameon=False,fontsize=8)
ax.set_title("inter-fitter noise floor vs anomaly",fontsize=10)
ax=axes[1]
dat=[ [r["dperi_deg"] for r in rows if r["cat"]=="code" and r["in_cap"]],
      [r["dperi_deg"] for r in rows if r["cat"]=="code" and not r["in_cap"]],
      [r["dperi_deg"] for r in rows if r["cat"]=="warsaw" and r["in_cap"]],
      [r["dperi_deg"] for r in rows if r["cat"]=="warsaw" and not r["in_cap"]] ]
ax.boxplot([np.log10(d) for d in dat],
           tick_labels=["code\nin","code\nout","war\nin","war\nout"],showfliers=False)
ax.set_ylabel("log10 |delta periapsis dir| (deg)")
ax.set_title("floor by cap membership",fontsize=10)
fig.tight_layout()
FIG=RESULTS/"figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG/"supplementary" / "step_b44_cross_fitter.png",dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b44_cross_fitter.png')