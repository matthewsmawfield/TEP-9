"""step_10: the proper clock-channel decomposition.

The previous tests measured cap membership on the raw osc->orig
correction -- dominated by Jupiter's kick.  The proper TEP observable is
the NON-PLANETARY part of the correction, tested through channels where
planets and a proper-time boundary predict different signs:

  P1 rotation-per-kick: ddir/|D(1/a)| -- gravity rotates the orbit in
     proportion to the energy exchanged; a lapse boundary rotates the
     reconstruction without exchanging energy.  Ratio isolates the
     non-planetary rotation.
  P2 inbound vs outbound: D(orig->future) vs theta_orig vs theta_fut.
     An ingress boundary leaves the discrepancy tied to the ARRIVAL
     direction; planets leave it directionless on both legs.
  P3 CMB-frame motion: continuous correlation of the correction with
     the comet's velocity projection along the CMB dipole and the ISM
     inflow -- the field's rest frame, not a cap.
  P4 rotation-axis coherence: the unit vector of each comet's
     orbit-plane rotation (p_orig x p_osc).  Boundary crossing rotates
     about a shared axis; planetary kicks rotate about scattered axes.
  P5 NG residuals: A1 (radial), A2 (transverse), A3 (normal) recoil
     components vs axis -- out-of-plane A3 excess = a push normal to
     the orbital plane, the domain-wall signature.
  P6 epoch trend: correction vs perihelion year -- a moving boundary.

Data: Warsaw catalogue tables a1/b/b4/c/d (Krolikowska 2014), all real.
Pre-declared: axis (lam 34, beta -13), cap 60 deg, spike 0<1/a<100e-6.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_033_three_leg")
tee_stdout(logger)
logger.header("Three-leg orbit decomposition")
import json, math, os
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
from scipy.stats import spearmanr, mannwhitneyu, kendalltau

RNG = np.random.default_rng(11)


def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def lonlat_vec(l_deg, b_deg):
    l, b = math.radians(l_deg), math.radians(b_deg)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

def gal_vec(l_deg, b_deg):
    l, b = math.radians(l_deg), math.radians(b_deg)
    g = np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])
    return GAL2ECL @ g

TNO = lonlat_vec(34.0, -13.0)
CMB = gal_vec(264.02, 48.25)   # dipole apex in ecliptic
ISM = lonlat_vec(255.8, 5.16)

def parse_orbit_table(path):
    rows = []
    for line in open(path):
        if len(line) < 115: continue
        try:
            rows.append(dict(sample=line[0:2].strip(), com=line[3].strip(),
                desig=line[5:17].strip(),
                tyr=int(line[27:31]),
                q=float(line[42:56]), e=float(line[56:70]),
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
        if d:
            out.setdefault(d, []).append(dict(
                model=line[132:140].strip(), datat=line[122:132].strip(),
                qnew=line[156:159].strip()))
    return out

def parse_b4(path):
    out = {}
    for line in open(path):
        if len(line) < 90: continue
        d = line[3:15].strip()
        try:
            out[d] = dict(A1=float(line[17:29]), A2=float(line[41:52]),
                          A3=float(line[66:77]), eA1=float(line[29:41]),
                          eA2=float(line[52:64]), eA3=float(line[77:86]),
                          model=line[111:115].strip())
        except ValueError:
            continue
    return out

orig_r = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
osc_r  = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tableb.dat"))
fut_r  = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tabled.dat"))
meta   = parse_a1(str(DATA_RAW / "warsaw" / "warsaw_tablea1.dat"))
ng     = parse_b4(str(DATA_RAW / "warsaw" / "warsaw_tableb4.dat"))

PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
def dedup(rows, pref):
    out = {}
    for r in rows:
        k = r["desig"]
        if k not in out or pref.get(r["com"], 9) < pref.get(out[k]["com"], 9):
            out[k] = r
    return out
orig = dedup(orig_r, PREF)
PREF_FUT = {"i": 0, "l": 0, "j": 2, "k": 2}
fut  = dedup(fut_r, PREF_FUT)
PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}
osc  = dedup(osc_r, PREF_OSC)

comets = []
for k, ro in orig.items():
    if k not in osc: continue
    rs = osc[k]
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    ps = perih_dir(math.radians(rs["w"]), math.radians(rs["Om"]), math.radians(rs["i"]))
    rot = np.cross(-po, -ps)          # rotation axis of aphelion shift
    rn = np.linalg.norm(rot)
    c = dict(desig=k, sample=ro["sample"], tyr=ro["tyr"],
        aph=-po, p_orig=-po, p_osc=-ps,
        rotaxis=rot / rn if rn > 0 else None,
        aa_orig=ro["aa"], aa_osc=rs["aa"],
        ddir=sep(-po, -ps), daa=rs["aa"] - ro["aa"],
        i_orig=ro["i"], q=ro["q"],
        model=meta.get(k, [{}])[0].get("model", "?"),
        datat=meta.get(k, [{}])[0].get("datat", "?"),
        vdot_cmb=float(np.dot(po, CMB)),   # inbound velocity ~ +perih dir
        vdot_ism=float(np.dot(po, ISM)),
        theta=sep(-po, TNO))
    if k in fut:
        rf = fut[k]
        pf = perih_dir(math.radians(rf["w"]), math.radians(rf["Om"]), math.radians(rf["i"]))
        c["theta_fut"] = sep(-pf, TNO)
        c["ddaaf"] = rf["aa"] - ro["aa"]               # orig -> future energy
        c["ddirf"] = sep(-po, -pf)                     # orig -> future direction
    if k in ng:
        c.update({kk: ng[k][kk] for kk in ("A1", "A2", "A3", "eA3")})
    comets.append(c)

spike = [c for c in comets if 0 < c["aa_orig"] < 100]
res = {"n_spike": len(spike)}
th = np.array([c["theta"] for c in spike])
cap = th < 60

# P1: rotation per unit energy kick
r1 = {}
ratio = np.array([c["ddir"] / max(abs(c["daa"]), 1e-9) for c in spike])
r1["median_ratio_in"] = float(np.median(ratio[cap]))
r1["median_ratio_out"] = float(np.median(ratio[~cap]))
r1["p_greater"] = float(mannwhitneyu(ratio[cap], ratio[~cap], alternative="greater").pvalue)
# same on full-arc only
full_idx = np.array([c["datat"] == "full" for c in spike])
if full_idx.sum() > 10:
    u = mannwhitneyu(ratio[cap & full_idx], ratio[~cap & full_idx], alternative="greater")
    r1["fullarc_median_in"] = float(np.median(ratio[cap & full_idx]))
    r1["fullarc_median_out"] = float(np.median(ratio[~cap & full_idx]))
    r1["fullarc_p"] = float(u.pvalue)
res["P1_rotation_per_kick"] = r1

# P2: inbound vs outbound leg
p2 = {}
have_f = [c for c in spike if "ddaaf" in c]
if len(have_f) > 15:
    to = np.array([c["theta"] for c in have_f])
    tf = np.array([c["theta_fut"] for c in have_f])
    daaf = np.abs(np.array([c["ddaaf"] for c in have_f]))
    ddirf = np.array([c["ddirf"] for c in have_f])
    p2["n"] = len(have_f)
    p2["corr_absDaaf_vs_theta_orig"] = {"rho": float(spearmanr(to, daaf)[0]), "p": float(spearmanr(to, daaf)[1])}
    p2["corr_absDaaf_vs_theta_fut"]  = {"rho": float(spearmanr(tf, daaf)[0]), "p": float(spearmanr(tf, daaf)[1])}
    p2["corr_ddirf_vs_theta_orig"]   = {"rho": float(spearmanr(to, ddirf)[0]), "p": float(spearmanr(to, ddirf)[1])}
    p2["corr_ddirf_vs_theta_fut"]    = {"rho": float(spearmanr(tf, ddirf)[0]), "p": float(spearmanr(tf, ddirf)[1])}
    for lab, tt in [("orig", to), ("fut", tf)]:
        inc = tt < 60
        u = mannwhitneyu(ddirf[inc], ddirf[~inc], alternative="greater")
        p2[f"ddirf_cap60_{lab}"] = {"n_in": int(inc.sum()),
            "med_in": float(np.median(ddirf[inc])), "med_out": float(np.median(ddirf[~inc])),
            "p": float(u.pvalue)}
res["P2_leg_asymmetry"] = p2

# P3: CMB / ISM velocity projection (continuous)
p3 = {}
for key, lab in [("vdot_cmb", "cmb"), ("vdot_ism", "ism")]:
    v = np.array([c[key] for c in spike])
    d = np.array([c["ddir"] for c in spike])
    rho, p = spearmanr(v, d)
    p3[f"ddir_vs_vdot_{lab}"] = {"rho": float(rho), "p": float(p)}
res["P3_frame_motion"] = p3

# P4: rotation-axis coherence -- are correction axes non-random inside cap?
p4 = {}
rots = {lab: np.array([c["rotaxis"] for c in spike
        if c["rotaxis"] is not None and (c["theta"] < 60) == (lab == "in")])
        for lab in ("in", "out")}
for lab, rv in rots.items():
    if len(rv) >= 8:
        R = np.linalg.norm(np.sum(rv, axis=0)) / len(rv)
        p4[f"Rbar_{lab}"] = float(R)
        p4[f"n_{lab}"] = len(rv)
        # isotropy p for axial data: use |cos| mean against random
        cs = np.abs(np.sum(rv * np.mean(rv, axis=0) / np.linalg.norm(np.mean(rv, axis=0)), axis=1))
        mc = []
        for _ in range(10000):
            x = RNG.normal(size=(len(rv), 3)); x /= np.linalg.norm(x, axis=1)[:, None]
            mm = np.mean(x, axis=0); mm /= np.linalg.norm(mm)
            mc.append(np.abs(np.sum(x * mm, axis=1)).mean())
        p4[f"p_coherence_{lab}"] = float((np.sum(np.array(mc) >= cs.mean()) + 1) / 10001)
# also: do in-cap rotation axes align with the TNO axis or its tangent?
if "in" in rots and len(rots["in"]) >= 8:
    al = np.abs(rots["in"] @ TNO)
    mc2 = []
    for _ in range(10000):
        x = RNG.normal(size=(len(rots["in"]), 3)); x /= np.linalg.norm(x, axis=1)[:, None]
        mc2.append(np.abs(x @ TNO).mean())
    p4["abs_cos_roaxis_tno_in"] = float(al.mean())
    p4["p_vs_random"] = float((np.sum(np.array(mc2) >= al.mean()) + 1) / 10001)
res["P4_rotation_axis"] = p4

# P5: NG components vs axis
p5 = {}
ngc = [c for c in spike if "A3" in c]
p5["n_ng"] = len(ngc)
if len(ngc) >= 10:
    t5 = np.array([c["theta"] for c in ngc])
    for key in ("A1", "A2", "A3"):
        v = np.array([c[key] for c in ngc])
        rho, p = spearmanr(t5, np.abs(v))
        p5[f"|{key}|_vs_theta"] = {"rho": float(rho), "p": float(p)}
    a3 = np.array([c["A3"] for c in ngc])
    inc5 = t5 < 60
    if inc5.sum() >= 3 and (~inc5).sum() >= 3:
        p5["A3_in_cap"] = {"vals": [f"{c['desig']}:{c['A3']:.2f}@th{c['theta']:.0f}" for c in ngc if c["theta"] < 60]}
        p5["median_absA3_in"] = float(np.median(np.abs(a3[inc5])))
        p5["median_absA3_out"] = float(np.median(np.abs(a3[~inc5])))
        u = mannwhitneyu(np.abs(a3[inc5]), np.abs(a3[~inc5]), alternative="greater")
        p5["p_absA3_greater_in"] = float(u.pvalue)
    # significance-weighted: A3/eA3 outliers near axis
    p5["A3_over_err"] = [{"desig": c["desig"], "A3": c["A3"], "sig": c["A3"]/c["eA3"] if c["eA3"] else None,
                          "theta": c["theta"]} for c in ngc if c["eA3"] and abs(c["A3"]/c["eA3"]) > 2]
res["P5_ng_recoil"] = p5

# P6: epoch trend inside cap -- moving boundary?
p6 = {}
sub = [c for c in spike if c["theta"] < 60 and c["tyr"] > 1800]
if len(sub) >= 10:
    t6 = np.array([c["tyr"] for c in sub]); d6 = np.array([c["ddir"] for c in sub])
    rho, p = spearmanr(t6, d6)
    p6["n"] = len(sub); p6["rho_ddir_vs_year_incap"] = float(rho); p6["p"] = float(p)
    kt = kendalltau(t6, d6)
    p6["kendall"] = {"tau": float(kt[0]), "p": float(kt[1])}
sub_out = [c for c in spike if c["theta"] >= 60 and c["tyr"] > 1800]
if len(sub_out) >= 10:
    t6 = np.array([c["tyr"] for c in sub_out]); d6 = np.array([c["ddir"] for c in sub_out])
    p6["rho_ddir_vs_year_outcap"] = {"rho": float(spearmanr(t6, d6)[0]), "p": float(spearmanr(t6, d6)[1])}
res["P6_epoch_trend"] = p6

out = str(RESULTS / "step_10_proper_clock.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)