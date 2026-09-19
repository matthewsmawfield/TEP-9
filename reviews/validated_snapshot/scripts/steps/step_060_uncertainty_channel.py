#!/usr/bin/env python3
"""
TEP-9 step 060 -- the uncertainty channel and third-lineage energy check
=========================================================================

The one-apparition comet catalogue (Krolikowska 2014, A&A 571, A63;
VizieR J/A+A/571/A63) adds two pieces of information the Warsaw and
CODE samples lack:

  1. per-element formal uncertainties on the osculating, original and
     future orbit solutions (tables b1, d1, e1), and
  2. the original / osculating / future inverse semimajor axes as
     catalogued independently by Marsden & Williams (2008) -- a third
     orbit-determination lineage (tablec1, MW08 columns).

With uncertainties, the reconstruction discrepancy can be NORMALIZED:
z_rot = d_of / sigma_dir asks whether the angular mismatch exceeds
what the catalogues themselves admit as error.  With MW08, the
question becomes whether two independent pipelines disagree more for
comets crossing the boundary direction.

Tests (n=38 unique comets, 1901-1950 one-apparition cohort; GR rows
preferred where a comet has multiple solutions):

  U1  the quality confound: published uncertainties (e_arg, e_long,
      e_i, e_aa, rms, Nres) in-cap vs out-cap.  If the cap selected
      poorly-determined orbits, its discrepancies are catalogued-
      error artifacts.
  U2  normalized rotation |d_of|/sigma_dir in-cap vs out-cap and vs
      angular distance to the axis.
  U3  normalized energy z_E = (aaori - aafut)/sigma(1/a): direction
      selectivity of the ENERGY leg, in units of the published
      uncertainty.
  U4  MW08 third lineage: |MWori - MWfut| energy-history jump and
      cross-catalogue |aaori - MWori| disagreement vs the axis.

Axes: the same candidate directions as step_059 (both TNO axes plus
controls).

Outputs: results/step_b25_uncertainty_channel.json,
         figures/step_b25_uncertainty_channel.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_060_uncertainty_channel")
tee_stdout(logger)
logger.header("Uncertainty channel and third-lineage check")

import json
import math
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
import xml.etree.ElementTree as ET
from scipy.stats import mannwhitneyu, spearmanr



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
}


def load_vot(path):
    t = ET.parse(path)
    rows = []
    for el in t.getroot().iter():
        if el.tag.split("}")[-1] == "TR":
            rows.append([td.text for td in list(el)])
    fields = [el.get("name") for el in t.getroot().iter()
              if el.tag.split("}")[-1] == "FIELD"]
    out = []
    for r in rows:
        d = dict(zip(fields, r))
        out.append(d)
    return out


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def mwu(a, b):
    if len(a) < 4 or len(b) < 4:
        return None
    u = mannwhitneyu(a, b, alternative="two-sided")
    return {"n_in": int(len(a)), "n_out": int(len(b)),
            "med_in": float(np.median(a)),
            "med_out": float(np.median(b)),
            "p_2sided": float(u.pvalue)}


def main():
    lpc = DATA_RAW / "lpc"
    osc = load_vot(lpc / "lpc_osc_2006_2010.vot")
    org = load_vot(lpc / "lpc_orig_2006_2010.vot")
    fut = load_vot(lpc / "lpc_fut_2006_2010.vot")
    aat = load_vot(lpc / "lpc_aaori_2006_2010.vot")

    # index by comet; prefer GR rows for one-per-comet dedup
    def dedup(rows):
        d = {}
        for r in rows:
            k = r["Comet"].strip()
            if k not in d or (d[k]["Model"].strip() != "GR"
                              and r["Model"].strip() == "GR"):
                d[k] = r
        return d

    o_d, s_d, f_d, a_d = (dedup(org), dedup(osc), dedup(fut), dedup(aat))
    common = sorted(set(o_d) & set(s_d) & set(f_d) & set(a_d))
    print(f"one-apparition cohort: {len(o_d)} comets, "
          f"{len(common)} with all four legs")

    recs = []
    for k in common:
        ro, rs, rf, ra = o_d[k], s_d[k], f_d[k], a_d[k]
        po = perih_dir(math.radians(f(ro["arg"])),
                       math.radians(f(ro["long"])),
                       math.radians(f(ro["i"])))
        ps = perih_dir(math.radians(f(rs["arg"])),
                       math.radians(f(rs["long"])),
                       math.radians(f(rs["i"])))
        pf = perih_dir(math.radians(f(rf["arg"])),
                       math.radians(f(rf["long"])),
                       math.radians(f(rf["i"])))
        # direction-shift proxy sigma: quadrature of published angular
        # element errors (deg); documented proxy, not a full cov
        sig_o = math.hypot(f(ro["e_arg"]), f(ro["e_long"]), f(ro["e_i"]))
        sig_f = math.hypot(f(rf["e_arg"]), f(rf["e_long"]), f(rf["e_i"]))
        sig_dir = math.hypot(sig_o, sig_f)
        z_e = (f(ra["aaori"]) - f(ra["aafut"])) / math.hypot(
            f(ra["e_aaori"]), f(ra["e_aafut"]))
        recs.append(dict(
            desig=k, model=ro["Model"].strip(),
            aa_ori=f(ra["aaori"]), q=f(ro["q"]),
            aph=-po,
            d_of=sep(-po, -pf), d_so=sep(-po, -ps), d_sf=sep(-ps, -pf),
            sig_dir=sig_dir, z_rot=sep(-po, -pf) / max(sig_dir, 1e-6),
            z_e=z_e,
            e_arg=f(ro["e_arg"]), e_long=f(ro["e_long"]),
            e_i=f(ro["e_i"]), e_aa=f(ro["e_aa"]),
            rms=f(ra["rms"]), nres=f(ra["Nr"]),
            MWori=f(ra["MWori"]), MWosc=f(ra["MWosc"]),
            MWfut=f(ra["MWfut"]),
            xline=abs(f(ra["aaori"]) - f(ra["MWori"])),
            mw_jump=abs(f(ra["MWori"]) - f(ra["MWfut"]))))

    res = {"n_cohort": len(recs),
           "note": ("one-apparition cohort J/A+A/571/A63; axes as in "
                    "step_059; sig_dir is a quadrature proxy of the "
                    "published angular element errors")}

    for axname, ax in AXES.items():
        th = np.array([sep(r["aph"], ax) for r in recs])
        inc = th < 60
        block = {}

        # ---- U1 quality confound ----
        aud = {}
        for key in ("e_arg", "e_long", "e_i", "e_aa", "rms", "nres"):
            t = mwu(np.array([r[key] for r in recs])[inc],
                    np.array([r[key] for r in recs])[~inc])
            if t:
                aud[key] = t
        block["U1_uncertainty_audit"] = aud

        # ---- U2 normalized rotation ----
        zr = np.array([r["z_rot"] for r in recs])
        d_ = np.array([r["d_of"] for r in recs])
        t = mwu(zr[inc], zr[~inc])
        if t:
            t["p_greater"] = float(mannwhitneyu(
                zr[inc], zr[~inc], alternative="greater").pvalue)
            rho, p = spearmanr(th, zr)
            t["spearman_rho"], t["spearman_p"] = float(rho), float(p)
            block["U2_norm_rotation"] = t
        t = mwu(d_[inc], d_[~inc])
        if t:
            t["p_greater"] = float(mannwhitneyu(
                d_[inc], d_[~inc], alternative="greater").pvalue)
            block["U2b_raw_rotation"] = t

        # ---- U3 normalized energy ----
        ze = np.array([r["z_e"] for r in recs])
        t = mwu(np.abs(ze)[inc], np.abs(ze)[~inc])
        if t:
            rho, p = spearmanr(th, ze)
            t["signed_spearman_rho"] = float(rho)
            t["signed_spearman_p"] = float(p)
            block["U3_norm_energy"] = t

        # ---- U4 MW08 third lineage ----
        mj = np.array([r["mw_jump"] for r in recs])
        xl = np.array([r["xline"] for r in recs])
        u4 = {}
        t = mwu(mj[inc], mj[~inc])
        if t:
            t["p_greater"] = float(mannwhitneyu(
                mj[inc], mj[~inc], alternative="greater").pvalue)
            u4["mw_energy_jump"] = t
        t = mwu(xl[inc], xl[~inc])
        if t:
            t["p_greater"] = float(mannwhitneyu(
                xl[inc], xl[~inc], alternative="greater").pvalue)
            rho, p = spearmanr(th, xl)
            t["spearman_rho"], t["spearman_p"] = float(rho), float(p)
            u4["cross_lineage_disagreement"] = t
        block["U4_mw08_third_lineage"] = u4

        res[axname] = block
        print(f"  {axname}: U2 zrot p={block.get('U2_norm_rotation', {}).get('p_greater')}, "
              f"U3 ze med-in={block.get('U3_norm_energy', {}).get('med_in')}")

    (RESULTS / "step_b25_uncertainty_channel.json").write_text(
        json.dumps(res, indent=1, default=float))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax_ = AXES["tno_extreme_34_-13"]
    th = np.array([sep(r["aph"], ax_) for r in recs])
    inc = th < 60
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
    ax[0].scatter(th[~inc], np.array([r["z_rot"] for r in recs])[~inc],
                  s=16, color="steelblue", label="out-cap")
    ax[0].scatter(th[inc], np.array([r["z_rot"] for r in recs])[inc],
                  s=16, color="crimson", label="in-cap")
    ax[0].axvline(60, color="k", ls="--", lw=1)
    ax[0].set(xlabel="angle to axis (deg)", ylabel="$z_{rot}$",
              title="U2: normalized rotation discrepancy")
    ax[0].legend(fontsize=7)
    ax[1].scatter(th[~inc], np.array([r["z_e"] for r in recs])[~inc],
                  s=16, color="steelblue")
    ax[1].scatter(th[inc], np.array([r["z_e"] for r in recs])[inc],
                  s=16, color="crimson")
    ax[1].axhline(0, color="0.5", lw=0.7)
    ax[1].axvline(60, color="k", ls="--", lw=1)
    ax[1].set(xlabel="angle to axis (deg)", ylabel="$z_E$",
              title="U3: normalized energy jump")
    ax[2].scatter(th[~inc], np.array([r["xline"] for r in recs])[~inc],
                  s=16, color="steelblue")
    ax[2].scatter(th[inc], np.array([r["xline"] for r in recs])[inc],
                  s=16, color="crimson")
    ax[2].axvline(60, color="k", ls="--", lw=1)
    ax[2].set(xlabel="angle to axis (deg)",
              ylabel="|1/a$_{orig}$ Warsaw - MW08|",
              title="U4: cross-lineage disagreement")
    fig.tight_layout()
    fig.savefig(RESULTS / "figures" / "step_b25_uncertainty_channel.png",
                dpi=150)
    print("wrote results/step_b25_uncertainty_channel.json, "
          "figures/step_b25_uncertainty_channel.png")


if __name__ == "__main__":
    main()
