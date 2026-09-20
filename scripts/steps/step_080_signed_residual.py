#!/usr/bin/env python3
"""Step 080 -- signed residual coherence: precession vs anomaly.

Step 075 found the catalogue signed rotation directionally coherent
on both catalogues, but with a population-wide sign bias (70-84% one
sign).  The mundane source of a global bias is secular omega
precession by the giant planets across the ~600 yr leg separation --
a real, directed torque.  This step measures the SIGNED rotation the
N-body planets themselves produce on the same bidirectional legs
(shell 250 AU, step-063 convention), subtracts it from the catalogue
signed rotation, and asks whether the residual sign remains coherent
in-cap.  If the in-cap sign coherence survives subtraction it is a
directed phase structure beyond the planets; if it collapses to
50/50 the sign was planetary.

Inputs : data/raw/code/code_*.html, data/raw/warsaw/warsaw_tablec.dat,
         data/raw/spice/de440s.bsp, results/step_b40_signed_slip.csv
Outputs: results/step_b45_signed_residual.json
         results/step_b45_signed_residual.csv
         results/figures/supplementary/step_b45_signed_residual.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_080_signed_residual")
tee_stdout(logger)
logger.header("Signed residual coherence: secular precession vs anomaly")

import csv
import json
import math
import re
import numpy as np
import rebound
import spiceypy as sp
from html.parser import HTMLParser
from scipy.stats import binomtest, spearmanr

# ------------------------------------------------------------------
# CODE table parser (identical to step_063/071)
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

def parse_code(path):
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

SPK = DATA_RAW / "spice" / "de440s.bsp"
if not SPK.exists():
    raise FileNotFoundError(f"JPL ephemeris missing: {SPK}")
sp.furnsh(str(SPK))

AU_KM  = 149597870.7
DAY_YR = 365.25
EPS    = math.radians(23.4392911)
RX     = np.array([[1,0,0],
                   [0,math.cos(EPS),math.sin(EPS)],
                   [0,-math.sin(EPS),math.cos(EPS)]])

GM_SUN = 1.32712440018e11
GM = {"1": 2.2031868551e4, "2": 3.2485859200e5, "3": 4.0350323562e5,
      "4": 4.2828375814e4, "5": 1.2671276480e8, "6": 3.7940626000e7,
      "7": 5.7945490100e6, "8": 6.8365271006e6, "9": 1.0868657e3}
PLANET_IDS = list(GM.keys())

def jd_tt(y,m,d):
    if m<=2: y-=1; m+=12
    A=y//100; B=2-A+A//4
    return int(365.25*(y+4716))+int(30.6001*(m+1))+d+B-1524.5

def perihelion_et(T):
    m=re.match(r"(\d{4})\s+(\d{1,2})\s+(\d+\.?\d*)",T)
    if not m: return None
    return (jd_tt(int(m.group(1)),int(m.group(2)),float(m.group(3)))-2451545.0)*86400.0

def body_state(body,et):
    st,_=sp.spkezr(body,et,"J2000","NONE","0")
    return RX@np.array(st[:3])/AU_KM, RX@np.array(st[3:])/AU_KM*86400*DAY_YR

# ------------------------------------------------------------------
# Integration (step_071 machinery, single 250 AU shell)
# ------------------------------------------------------------------

MU     = 4*math.pi**2
SHELLS = [250.0]
T_MAX  = 40000.0
DT_OUT = 1.0

def init_sim(et):
    sim=rebound.Simulation()
    sim.units=("AU","yr","Msun")
    ps,vs=body_state("0",et)
    sim.add(x=ps[0],y=ps[1],z=ps[2],vx=vs[0],vy=vs[1],vz=vs[2],m=1.0)
    for pid in PLANET_IDS:
        p_,v_=body_state(pid,et)
        sim.add(x=p_[0],y=p_[1],z=p_[2],vx=v_[0],vy=v_[1],vz=v_[2],
                m=GM[pid]/GM_SUN)
    sim.move_to_com()
    return sim,ps,vs

def state_at_periapsis(ro):
    mu=4*math.pi**2
    q,e=ro["q"],ro["e"]
    Om,w,i=math.radians(ro["Om"]),math.radians(ro["w"]),math.radians(ro["i"])
    a=q/(1-e) if e<1 else -q/(e-1)
    n=math.sqrt(mu/abs(a)**3)
    vq=math.sqrt(mu*(1+e)/q)
    co,so,cO,sO,ci,si=np.cos(w),np.sin(w),np.cos(Om),np.sin(Om),np.cos(i),np.sin(i)
    P=np.array([cO*co-sO*so*ci, sO*co+cO*so*ci, so*si])
    Q=np.array([-cO*so-sO*co*ci, -sO*so+cO*co*ci, co*si])
    return q*P, vq*Q

def boundary_orbit(r_rel,v_rel,mu):
    r=np.linalg.norm(r_rel)
    h=np.cross(r_rel,v_rel)
    evec=(np.cross(v_rel,h)/mu)-r_rel/r
    en=np.linalg.norm(evec)
    phat=evec/en if en>1e-12 else r_rel/r
    return phat

def integrate_leg(ro,et,direction):
    sim,ps,vs=init_sim(et)
    rvec,vvec=state_at_periapsis(ro)
    sim.add(x=rvec[0]+ps[0],y=rvec[1]+ps[1],z=rvec[2]+ps[2],
            vx=vvec[0]+vs[0],vy=vvec[1]+vs[1],vz=vvec[2]+vs[2])
    nc=sim.N-1
    sim.integrator="ias15"
    t=direction*DT_OUT
    while abs(t)<T_MAX:
        sim.integrate(t,exact_finish_time=0)
        p=sim.particles
        r_rel=np.array([p[nc].x-p[0].x,p[nc].y-p[0].y,p[nc].z-p[0].z])
        rr=np.linalg.norm(r_rel)
        if rr>=SHELLS[0]:
            mtot=sum(pp.m for pp in sim.particles)
            rb=np.zeros(3); vb=np.zeros(3)
            for pp in sim.particles:
                rb+=pp.m*np.array([pp.x,pp.y,pp.z])
                vb+=pp.m*np.array([pp.vx,pp.vy,pp.vz])
            rb/=mtot; vb/=mtot
            r_b=np.array([p[nc].x,p[nc].y,p[nc].z])-rb
            v_b=np.array([p[nc].vx,p[nc].vy,p[nc].vz])-vb
            return boundary_orbit(r_b,v_b,MU*mtot)
        t+=direction*DT_OUT
    return None

def signed_rot(p_from,p_to,h):
    """signed angle p_from -> p_to about pole h (deg); + = along motion."""
    x=np.cross(p_from,p_to)
    s=np.dot(x,h)
    c=np.dot(p_from,p_to)
    return math.degrees(math.atan2(s,c))

def pole(Om,inc):
    return np.array([math.sin(inc)*math.sin(Om),
                     -math.sin(inc)*math.cos(Om),
                     math.cos(inc)])

def lv(l,b):
    l,b=math.radians(l),math.radians(b)
    return np.array([math.cos(b)*math.cos(l),math.cos(b)*math.sin(l),math.sin(b)])

def sep(a,b):
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b),-1,1))))

def perih_dir(om,Om,inc):
    co,so,cO,sO,ci,si=np.cos(om),np.sin(om),np.cos(Om),np.sin(Om),np.cos(inc),np.sin(inc)
    return np.array([cO*co-sO*so*ci, sO*co+cO*so*ci, so*si])

TNO=lv(34.0,-13.0)
CAP=60.0

# ------------------------------------------------------------------
# Sample (identical to step_063) + catalogue signed rotations (b40)
# ------------------------------------------------------------------

osc  = parse_code(str(DATA_RAW/"code"/"code_osculating.html"))
orig = parse_code(str(DATA_RAW/"code"/"code_original.html"))
fut  = parse_code(str(DATA_RAW/"code"/"code_future.html"))
warsaw={l[5:17].strip() for l in open(str(DATA_RAW/"warsaw"/"warsaw_tablec.dat")) if len(l)>115}

cat_srot={}
with open(RESULTS/"step_b40_signed_slip.csv") as f:
    for r in csv.DictReader(f):
        if r["cohort"]=="code":
            cat_srot[r["desig"]]=float(r["srot"])

sample=[]
for k,ro in orig.items():
    if k not in fut or k not in osc or k in warsaw or k not in cat_srot: continue
    if not (0<ro["aa"]<100 and ro["cls"] in ("1a","1a+","1b")): continue
    et=perihelion_et(osc[k]["T"])
    if et is None: continue
    sample.append((k,ro,osc[k],et))

logger.info(f"sample: {len(sample)} class-1 CODE comets")

def _run_comet(k,ro,oo,et,pb,pf_):
    if pb is None or pf_ is None:
        return None
    Om_,i_=math.radians(oo["Om"]),math.radians(oo["i"])
    h=pole(Om_,i_)
    srot_sim=signed_rot(pb,pf_,h)          # back -> fwd, same sense as orig->fut
    srot_cat=cat_srot[k]
    return dict(desig=k,
                theta=sep(-perih_dir(math.radians(ro["w"]),math.radians(ro["Om"]),i_),TNO),
                srot_cat=srot_cat,srot_sim=srot_sim,
                resid=srot_cat-srot_sim)

def _work_leg(job):
    # one leg per job: a few comets integrate far longer than the rest,
    # so splitting back/fwd legs keeps the pool balanced
    k,oo,et,direction=job
    try:
        return integrate_leg(oo,et,direction)
    except Exception as e:
        return ("__error__",k,str(e))

from scripts.utils.parallel import default_workers as _default_workers

if "--workers" in sys.argv:
    N_WORK=max(1,int(sys.argv[sys.argv.index("--workers")+1]))
else:
    N_WORK=_default_workers()

leg_jobs=[(k,oo,et,d) for (k,ro,oo,et) in sample for d in (-1,+1)]
if N_WORK>1 and len(leg_jobs)>1:
    import multiprocessing as mp
    ctx=mp.get_context("fork") if sys.platform!="win32" \
        else mp.get_context("spawn")

    def _init():
        # forked children inherit the parent's BSP fd; concurrent spkezr
        # reads through a shared descriptor corrupt each other -- reopen
        # the kernel so each worker holds its own file handle.
        sp.kclear()
        sp.furnsh(str(SPK))

    with ctx.Pool(min(N_WORK,len(leg_jobs)),initializer=_init) as pool:
        legs=pool.map(_work_leg,leg_jobs,chunksize=1)
else:
    legs=[_work_leg(j) for j in leg_jobs]

rows=[]; failed=[]
for idx,(k,ro,oo,et) in enumerate(sample):
    pb,pf_=legs[2*idx],legs[2*idx+1]
    for rec in (pb,pf_):
        if isinstance(rec,tuple) and rec[0]=="__error__":
            failed.append(rec[1]); logger.warning(f"{rec[1]}: {rec[2]}")
    if k in failed:
        continue
    if pb is None or pf_ is None:
        failed.append(k); logger.warning(f"{k}: shell not reached"); continue
    rows.append(_run_comet(k,ro,oo,et,pb,pf_))
    if (idx+1)%25==0: logger.info(f"  {idx+1}/{len(sample)}")

logger.info(f"integrated {len(rows)}, failed {len(failed)}")

# ------------------------------------------------------------------
# Statistics
# ------------------------------------------------------------------

res={"method":"signed simulated orig->fut rotation (250 AU shell, "
       "about the osculating pole) subtracted from the catalogue "
       "signed rotation of step_075; residual sign coherence tested "
       "in and out of the 60-deg cap.",
     "cap_deg":CAP,"axis":"34,-13"}

def signstats(sub,tag):
    x=np.array([r["resid"] for r in sub])
    npos=int((x>0).sum())
    bt=binomtest(npos,len(x),0.5)
    # declared direction: the boundary slip leaves a positive residual
    # (a lapse) once the planetary rotation is subtracted -> greater tail
    btg=binomtest(npos,len(x),0.5,alternative="greater")
    res[tag]={"n":len(x),"n_pos":npos,"frac_pos":npos/len(x),
              "median_resid":float(np.median(x)),
              "p_binom":float(btg.pvalue),
              "p_binom_2sided":float(bt.pvalue),
              "median_srot_sim":float(np.median([r["srot_sim"] for r in sub])),
              "median_srot_cat":float(np.median([r["srot_cat"] for r in sub]))}
    logger.info(f"{tag}: n={len(x)} resid pos={npos} ({npos/len(x):.2f}) "
                f"p={bt.pvalue:.4f} med_resid={np.median(x):+.4f} "
                f"(sim med={res[tag]['median_srot_sim']:+.4f}, "
                f"cat med={res[tag]['median_srot_cat']:+.4f})")

th=np.array([r["theta"] for r in rows])
inn=[r for r in rows if r["theta"]<CAP]
outt=[r for r in rows if r["theta"]>=CAP]
signstats(inn,"in_cap"); signstats(outt,"out_cap"); signstats(rows,"all")

# does the simulation predict the catalogue sign at all?
rho,p=spearmanr([r["srot_sim"] for r in rows],[r["srot_cat"] for r in rows])
res["sim_vs_cat_sign"]={"rho":float(rho),"p_2sided":float(p)}
logger.info(f"srot_sim vs srot_cat: rho={rho:+.3f} p={p:.4f}")

with open(RESULTS/"step_b45_signed_residual.json","w") as f:
    json.dump(res,f,indent=1)
with open(RESULTS/"step_b45_signed_residual.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["desig","theta","srot_cat","srot_sim","resid"])
    w.writeheader()
    for r in rows: w.writerow(r)
logger.data_save(RESULTS/'step_b45_signed_residual.json')
logger.data_save(RESULTS/'step_b45_signed_residual.csv')
# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig,axes=plt.subplots(1,2,figsize=(9,4.2))
ax=axes[0]
sc=[r["srot_cat"] for r in rows]; ss=[r["srot_sim"] for r in rows]
cc=["crimson" if r["theta"]<CAP else "0.6" for r in rows]
ax.scatter(ss,sc,s=9,c=cc,alpha=0.7)
ax.axhline(0,color="k",ls=":",lw=0.8); ax.axvline(0,color="k",ls=":",lw=0.8)
ax.plot([-1,1],[-1,1],color="k",ls="--",lw=0.8)
ax.set_xlabel("signed simulated rotation (deg)")
ax.set_ylabel("signed catalogue rotation (deg)")
ax.set_title("catalogue vs simulated sign",fontsize=10)
ax=axes[1]
for sub,col,lab in ((inn,"crimson","in-cap"),(outt,"0.5","out-of-cap")):
    vals=sorted(r["resid"] for r in sub)
    ax.plot(vals,np.linspace(0,1,len(vals)),color=col,lw=1.4,label=lab)
ax.axvline(0,color="k",ls=":",lw=0.8)
ax.set_xlabel("residual signed rotation (deg)")
ax.set_ylabel("CDF"); ax.legend(frameon=False,fontsize=8)
ax.set_title("residual sign coherence",fontsize=10)
fig.tight_layout()
FIG=RESULTS/"figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG/"supplementary" / "step_b45_signed_residual.png",dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b45_signed_residual.png')