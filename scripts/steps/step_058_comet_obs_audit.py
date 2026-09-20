#!/usr/bin/env python3
"""
TEP-9 step 058 -- comet observational-leverage audit
=====================================================

The transit signature (steps 030-041) rests on the three-leg orbit
reconstruction discrepancy d_of being larger for comets whose original
aphelia fall inside the boundary cap.  The dangerous confound is
observational leverage: if in-cap comets happened to be worse
observed -- fewer observations, shorter arcs, or arcs concentrated on
one side of perihelion -- a standard orbit fit would reconstruct them
worse regardless of any boundary, manufacturing a directional signal
from nothing.

This step measures the observing circumstances directly, on both
catalogues:

  CODE (class-1 matched sample, N=54):
    - nobs, arc length, perihelion year, q, i in-cap vs out-cap
    - arc_frac = (T_peri - arc_start)/(arc_end - arc_start):
      the fraction of the observed arc lying BEFORE perihelion.
      frac ~ 0 or 1 means a one-sided arc (maximal leverage);
      ~0.5 means symmetric data around perihelion.
    - d_of retested on the balanced-arc subset (0.15<frac<0.85)
    - regression control: d_of residualized on (nobs, arc, arc_frac),
      then re-tested against the axis
    - matched-pairs control on (nobs, arc, q)

  Warsaw (spike sample, N=99):
    - tablea1 supplies datat (pre/post/full arc symmetry), Nobs,
      arc years, RMS, Nres, and the heliocentric distances of the
      first/last observations (dh1/dh2)
    - datat composition in-cap vs out-cap; the orbit-quality and
      residual channels audited the same way

Outputs: results/step_b23_comet_obs_audit.json,
         figures/supplementary/step_b23_comet_obs_audit.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_058_comet_obs_audit")
tee_stdout(logger)
logger.header("Comet observational-leverage audit")

import json
import math
import re
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact
from html.parser import HTMLParser

RNG = np.random.default_rng(20260918)



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
                     math.cos(b) * math.sin(l),
                     math.sin(b)])


# pre-declared comet-side axis (extreme-TNO cluster axis, a>250)
TNO = lv(34.0, -13.0)


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
            dates = re.findall(r"(\d{4})\s+(\d{2})\s+(\d{2})", r[5])
            arc_days = np.nan
            if len(dates) >= 2:
                (y1, m1, d1), (y2, m2, d2) = dates[0], dates[-1]
                arc_days = abs((int(y2) - int(y1)) * 365.25
                               + (int(m2) - int(m1)) * 30.44
                               + (int(d2) - int(d1)))
                arc_start = (int(y1), int(m1), int(d1))
            else:
                arc_start = None
            # perihelion date "YYYY MM DD.ddd"
            mt = re.match(r"\s*(\d{4})\s+(\d{2})\s+(\d+\.?\d*)", r[7])
            t_year = np.nan
            if mt:
                t_year = (int(mt.group(1)) + (int(mt.group(2)) - 1) / 12
                          + float(mt.group(3)) / 365.25)
            out[r[0].strip()] = dict(
                desig=r[0].strip(), model=r[1].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                nobs=float(r[4]) if r[4].strip() else np.nan,
                arc_days=arc_days,
                t_year=t_year,
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
            if arc_start is not None and np.isfinite(t_year):
                y, m, d = arc_start
                start_year = y + (m - 1) / 12 + d / 365.25
                arc_years = out[r[0].strip()]["arc_days"] / 365.25
                out[r[0].strip()]["arc_frac"] = (
                    (t_year - start_year) / arc_years
                    if arc_years > 0 else np.nan)
            else:
                out[r[0].strip()]["arc_frac"] = np.nan
        except (ValueError, IndexError):
            continue
    return out


def parse_warsaw_orbits(path):
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(sample=line[0:2].strip(), com=line[3].strip(),
                desig=line[5:17].strip(),
                q=float(line[42:56]), e=float(line[56:70]),
                w=float(line[70:82]), Om=float(line[82:94]),
                i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows


def parse_warsaw_a1(path):
    """tablea1: observational material per designation.

    The table lists several rows per comet -- the canonical full-arc
    solution (datat contains "full") plus PRE/POST/DIST subset rows.
    The full-arc row carries the observational material this audit
    needs; prefer it over any subset row seen earlier."""
    out = {}
    for line in open(path):
        if len(line) < 160:
            continue
        d = line[3:15].strip()
        if not d:
            continue
        try:
            rec = dict(
                nobs=int(line[88:93]), arcy=float(line[97:102]),
                dh1=float(line[106:111]), dh2=float(line[112:117]),
                datat=line[122:132].strip(), model=line[132:140].strip(),
                rms=float(line[163:167]), nres=int(line[167:173]),
                qnew=line[156:159].strip())
        except (ValueError, IndexError):
            continue
        if "full" in rec["datat"]:
            out[d] = rec
        elif d not in out:
            out[d] = rec
    return out


def mwu(vin, vout):
    if len(vin) < 4 or len(vout) < 4:
        return None
    u = mannwhitneyu(vin, vout, alternative="two-sided")
    return {"n_in": int(len(vin)), "n_out": int(len(vout)),
            "med_in": float(np.median(vin)),
            "med_out": float(np.median(vout)),
            "p_2sided": float(u.pvalue)}


def main():
    res = {}

    # ================= CODE matched sample =================
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
        pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]),
                       math.radians(fut[k]["i"]))
        rows.append(dict(desig=k, warsaw=k in warsaw,
            aa=ro["aa"], q=ro["q"], i=ro["i"], cls=ro["cls"],
            nobs=ro["nobs"], arc_days=ro["arc_days"],
            arc_frac=ro["arc_frac"], t_year=ro["t_year"],
            aph=-po, d_of=sep(-po, -pf)))

    sp = [r for r in rows
          if 0 < r["aa"] < 100 and r["q"] < 3.1
          and r["cls"] in ("1a", "1a+", "1b")]
    co = [r for r in sp if not r["warsaw"]]
    print(f"CODE-only matched: {len(co)}")
    res["n_code_only"] = len(co)

    th = np.array([sep(r["aph"], TNO) for r in co])
    inc = th < 60
    dof = np.array([r["d_of"] for r in co])
    res["cap_def"] = {"axis": "(34,-13) pre-declared extreme-TNO",
                      "cap_deg": 60, "n_in": int(inc.sum()),
                      "n_out": int((~inc).sum()),
                      "dof_p": float(mannwhitneyu(
                          dof[inc], dof[~inc],
                          alternative="greater").pvalue)}

    # ---- A1: circumstance audit, in vs out ----
    audit = {}
    for key in ("nobs", "arc_days", "arc_frac", "t_year", "q", "i", "aa"):
        vin = np.array([r[key] for r, m in zip(co, inc)
                        if m and np.isfinite(r[key])])
        vout = np.array([r[key] for r, m in zip(co, ~inc)
                         if m and np.isfinite(r[key])])
        t = mwu(vin, vout)
        if t:
            audit[key] = t
            print(f"  {key:9s} in={t['med_in']:.3g} out={t['med_out']:.3g} "
                  f"p={t['p_2sided']:.3f}")
    res["A1_circumstance_audit"] = {
        "rows": audit,
        "note": ("if the cap selected poorly-observed comets (fewer "
                 "nobs, shorter or one-sided arcs), the d_of excess "
                 "would be an observational artifact; all p_2sided "
                 "should be non-significant for the claim to stand")}

    # ---- A2: balanced-arc retest ----
    bal = np.array([np.isfinite(r["arc_frac"])
                    and 0.15 < r["arc_frac"] < 0.85 for r in co])
    res["A2_balanced_arc"] = {"n": int(bal.sum())}
    if inc[bal].sum() >= 4 and (~inc & bal).sum() >= 4:
        u = mannwhitneyu(dof[bal & inc], dof[bal & ~inc],
                         alternative="greater")
        res["A2_balanced_arc"].update(
            n_in=int((bal & inc).sum()), n_out=int((bal & ~inc).sum()),
            med_in=float(np.median(dof[bal & inc])),
            med_out=float(np.median(dof[bal & ~inc])),
            p=float(u.pvalue))
        print(f"A2 balanced-arc: {res['A2_balanced_arc']}")

    # ---- A3: regression on observing circumstances ----
    mask = np.array([np.isfinite(r["nobs"]) and np.isfinite(r["arc_days"])
                     and np.isfinite(r["arc_frac"]) for r in co])
    X = np.column_stack([
        np.array([r["nobs"] for r in co])[mask],
        np.array([r["arc_days"] for r in co])[mask],
        np.array([r["arc_frac"] for r in co])[mask]])
    X = np.column_stack([np.ones(len(X)), X])
    y = dof[mask]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    inc_m = inc[mask]
    u = mannwhitneyu(resid[inc_m], resid[~inc_m], alternative="greater")
    rho, pr = spearmanr(th[mask], resid)
    res["A3_regression_control"] = {
        "n": int(mask.sum()),
        "coefficients": {"nobs": float(beta[1]),
                         "arc_days": float(beta[2]),
                         "arc_frac": float(beta[3])},
        "residual_cap": {"med_in": float(np.median(resid[inc_m])),
                         "med_out": float(np.median(resid[~inc_m])),
                         "p": float(u.pvalue)},
        "residual_spearman": {"rho": float(rho), "p": float(pr)}}
    print(f"A3 regression residual cap p={u.pvalue:.4f} "
          f"rho={rho:.3f} p={pr:.4f}")

    # ---- A4: matched pairs on (nobs, arc_days, q) ----
    z = np.column_stack([
        (np.array([r["nobs"] for r in co]) - 200) / 400,
        np.array([r["arc_days"] for r in co]) / 3000,
        np.array([r["q"] for r in co]) / 2])
    pairs, used = [], set()
    for j in np.where(inc)[0]:
        best, bd = None, 1e9
        for k in np.where(~inc)[0]:
            if k in used:
                continue
            d = np.linalg.norm(z[j] - z[k])
            if d < bd:
                best, bd = k, d
        if best is not None:
            used.add(best)
            pairs.append((dof[j], dof[best]))
    pv = np.array(pairs)
    from scipy.stats import wilcoxon
    try:
        w = wilcoxon(pv[:, 0], pv[:, 1], alternative="greater")
        res["A4_matched_circumstance"] = {
            "n_pairs": len(pairs),
            "med_in": float(np.median(pv[:, 0])),
            "med_out": float(np.median(pv[:, 1])),
            "p_wilcoxon": float(w.pvalue)}
    except ValueError:
        res["A4_matched_circumstance"] = {"n_pairs": len(pairs)}
    print(f"A4 matched pairs: {res['A4_matched_circumstance']}")

    # ================= Warsaw spike sample =================
    a1 = parse_warsaw_a1(DATA_RAW / "warsaw" / "warsaw_tablea1.dat")
    orows = parse_warsaw_orbits(DATA_RAW / "warsaw" / "warsaw_tablec.dat")
    frows = parse_warsaw_orbits(DATA_RAW / "warsaw" / "warsaw_tabled.dat")
    PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
    PREF_FUT = {"i": 0, "l": 0, "j": 2, "k": 2}
    o_dedup, f_dedup = {}, {}
    for r in orows:
        k = r["desig"]
        if k not in o_dedup or PREF.get(r["com"], 9) < PREF.get(
                o_dedup[k]["com"], 9):
            o_dedup[k] = r
    for r in frows:
        k = r["desig"]
        if k not in f_dedup or PREF_FUT.get(r["com"], 9) < PREF_FUT.get(
                f_dedup[k]["com"], 9):
            f_dedup[k] = r

    orows_b = parse_warsaw_orbits(DATA_RAW / "warsaw" / "warsaw_tableb.dat")
    PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}
    s_dedup = {}
    for r in orows_b:
        k = r["desig"]
        if k not in s_dedup or PREF_OSC.get(r["com"], 9) < PREF_OSC.get(
                s_dedup[k]["com"], 9):
            s_dedup[k] = r

    wrows = []
    for k, ro in o_dedup.items():
        if k not in f_dedup or k not in a1:
            continue
        po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                       math.radians(ro["i"]))
        pf = perih_dir(math.radians(f_dedup[k]["w"]),
                       math.radians(f_dedup[k]["Om"]),
                       math.radians(f_dedup[k]["i"]))
        rec = dict(desig=k, aa=ro["aa"], q=ro["q"],
                   aph=-po, d_of=sep(-po, -pf), **a1[k])
        if k in s_dedup:
            ps = perih_dir(math.radians(s_dedup[k]["w"]),
                           math.radians(s_dedup[k]["Om"]),
                           math.radians(s_dedup[k]["i"]))
            rec["d_so"] = sep(-po, -ps)
        else:
            rec["d_so"] = np.nan
        wrows.append(rec)
    wsp = [r for r in wrows if 0 < r["aa"] < 100]
    wth = np.array([sep(r["aph"], TNO) for r in wsp])
    winc = wth < 60
    wdof = np.array([r["d_of"] for r in wsp])
    res["warsaw"] = {"n_spike": len(wsp)}

    waudit = {}
    for key in ("nobs", "arcy", "rms", "nres", "dh1", "dh2"):
        vin = np.array([r[key] for r, m in zip(wsp, winc) if m])
        vout = np.array([r[key] for r, m in zip(wsp, winc) if not m])
        t = mwu(vin, vout)
        if t:
            waudit[key] = t
            print(f"  W {key:6s} in={t['med_in']:.3g} "
                  f"out={t['med_out']:.3g} p={t['p_2sided']:.3f}")
    # datat composition: fraction 'full' (two-sided arc)
    dt = np.array([r["datat"] for r in wsp])
    full_in = int(((dt == "full") & winc).sum())
    nfull_in = int((winc & (dt != "full")).sum())
    full_out = int(((dt == "full") & ~winc).sum())
    nfull_out = int((~winc & (dt != "full")).sum())
    orr, pf_ = fisher_exact([[full_in, nfull_in],
                             [full_out, nfull_out]])
    res["warsaw"]["A5_circumstance_audit"] = waudit
    res["warsaw"]["A6_datat_composition"] = {
        "full_in_cap": [full_in, nfull_in],
        "full_out_cap": [full_out, nfull_out],
        "fisher_p_2sided": float(pf_),
        "note": ("'full' = data both sides of perihelion; if in-cap "
                 "comets were preferentially one-sided (pre/post), "
                 "the reconstruction gap would be leverage-driven")}
    print(f"  W datat full: {full_in}/{full_in+nfull_in} in vs "
          f"{full_out}/{full_out+nfull_out} out, fisher p={pf_:.3f}")

    # d_of and d_so within 'full'-arc subset only
    fmask = dt == "full"
    if (fmask & winc).sum() >= 4 and (fmask & ~winc).sum() >= 4:
        u = mannwhitneyu(wdof[fmask & winc], wdof[fmask & ~winc],
                         alternative="greater")
        res["warsaw"]["A7_full_arc_only"] = {
            "n_in": int((fmask & winc).sum()),
            "n_out": int((fmask & ~winc).sum()),
            "med_in": float(np.median(wdof[fmask & winc])),
            "med_out": float(np.median(wdof[fmask & ~winc])),
            "p_dof": float(u.pvalue)}
        wso = np.array([r["d_so"] for r in wsp])
        m2 = fmask & np.isfinite(wso)
        if (m2 & winc).sum() >= 4 and (m2 & ~winc).sum() >= 4:
            u2 = mannwhitneyu(wso[m2 & winc], wso[m2 & ~winc],
                              alternative="greater")
            res["warsaw"]["A7_full_arc_only"].update(
                d_so_n_in=int((m2 & winc).sum()),
                d_so_n_out=int((m2 & ~winc).sum()),
                d_so_med_in=float(np.median(wso[m2 & winc])),
                d_so_med_out=float(np.median(wso[m2 & ~winc])),
                p_dso=float(u2.pvalue))
        print(f"  W full-arc-only: {res['warsaw']['A7_full_arc_only']}")

    (RESULTS / "step_b23_comet_obs_audit.json").write_text(
        json.dumps(res, indent=1, default=float))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
    nb_in = np.array([r["nobs"] for r, m in zip(co, inc) if m])
    nb_out = np.array([r["nobs"] for r, m in zip(co, ~inc) if ~m])
    ax[0].hist(nb_in, bins=15, alpha=0.6, color="crimson",
               label=f"in-cap (n={len(nb_in)})")
    ax[0].hist(nb_out, bins=15, alpha=0.6, color="steelblue",
               label=f"out-cap (n={len(nb_out)})")
    ax[0].set(xlabel="n observations", ylabel="N",
              title="A1: data volume in vs out")
    ax[0].legend(fontsize=7)
    af_in = np.array([r["arc_frac"] for r, m in zip(co, inc)
                      if m and np.isfinite(r["arc_frac"])])
    af_out = np.array([r["arc_frac"] for r, m in zip(co, ~inc)
                       if ~m and np.isfinite(r["arc_frac"])])
    ax[1].hist(af_in, bins=np.linspace(0, 1, 16), alpha=0.6,
               color="crimson", label="in-cap")
    ax[1].hist(af_out, bins=np.linspace(0, 1, 16), alpha=0.6,
               color="steelblue", label="out-cap")
    ax[1].axvline(0.5, color="k", ls=":", lw=1)
    ax[1].set(xlabel="arc fraction before perihelion", ylabel="N",
              title="arc symmetry in vs out")
    ax[1].legend(fontsize=7)
    ax[2].scatter(th[~inc], dof[~inc], s=14, color="steelblue",
                  label="out-cap")
    ax[2].scatter(th[inc], dof[inc], s=14, color="crimson",
                  label="in-cap")
    ax[2].axvline(60, color="k", ls="--", lw=1)
    ax[2].set(xlabel="angle to axis (deg)",
              ylabel="$d_{of}$ (deg)",
              title="discrepancy vs axis distance")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(RESULTS / "figures" / "supplementary" / "step_b23_comet_obs_audit.png",
                dpi=300)
    logger.data_save(RESULTS / "step_b23_comet_obs_audit.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b23_comet_obs_audit.png")


if __name__ == "__main__":
    main()
