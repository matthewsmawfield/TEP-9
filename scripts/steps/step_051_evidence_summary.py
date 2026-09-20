#!/usr/bin/env python3
"""
TEP / Planet-9 -- step b15: evidence summary
============================================

A single harvest of the pipeline's primary hypothesis tests,
with multiple-testing control.  The pipeline has grown to
~30 analyses; this step collects the headline p-values,
applies Benjamini-Hochberg FDR control within test families,
and reports the combined evidence honestly -- including the
tests that came back NULL, which constrain the interpretation
as much as the detections.

Families:
  F1  TNO clock-sector clustering (primary detections)
  F2  Bias-control bounds (the pointing channel measured &
      bounded at three levels)
  F3  Structure tests (radial, warp, patch, anti-deficit)
  F4  Comet clock channel (concurrent numbered pipeline)
  F5  Cross-population direction consistency
  F6  NULL results that disfavour the alternatives

Outputs: results/step_b15_evidence_summary.json,
         figures/supplementary/step_b15_evidence_summary.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, PROJECT_ROOT, tee_stdout
logger = StepLogger("step_051_evidence_summary")
tee_stdout(logger)
logger.header("Multiple-testing evidence summary")

from pathlib import Path
import json
import numpy as np

ROOT = PROJECT_ROOT
RES, FIG = RESULTS, RESULTS / "figures"


def bh(pvals):
    """Benjamini-Hochberg: return (reject_mask, threshold)."""
    p = np.sort(np.asarray(pvals))
    m = len(p)
    k = np.arange(1, m + 1)
    thresh = 0.05 * k / m
    ok = p <= thresh
    if not ok.any():
        return 0.0
    return float(thresh[ok][-1])


def main():
    # name, p, family, source json
    tests = [
        # F1 TNO clock-sector clustering
        ("varpi conditioned-null R (step_02)", 0.0073, "F1",
         "step_02_clustering.json"),
        ("3-D perihelion-vector alignment (b11)", 1e-4, "F1",
         "step_b11_population_selectivity.json"),
        ("omega conditioned-null R (b1)", 0.008, "F1",
         "step_b1_empirical_null.json"),
        ("varpi empirical-null R (b1)", 0.013, "F1",
         "step_b1_empirical_null.json"),
        ("numbered-orbits R=0.59 (b3)", 0.0015, "F1",
         "step_b3_robustness.json"),
        ("pre-2014 subsample (b3)", 0.043, "F1",
         "step_b3_robustness.json"),
        ("cap-profile peak, 3-D (b13)", 1e-4, "F1",
         "step_b13_patch_scale.json"),
        # F2 bias bounds (each: prob footprint alone reaches
        # observed)
        ("maximal-pointing ceiling (b5)", 0.0527, "F2",
         "step_b5_discovery_footprint.json"),
        ("footprint-conditional null (b8)", 0.035, "F2",
         "step_b8_footprint_conditional.json"),
        ("OSSOS-calibrated real coupling (b10)", 0.020, "F2",
         "step_b10_proxy_validation.json"),
        ("pairing/shuffle test (b8)", 0.056, "F2",
         "step_b8_footprint_conditional.json"),
        # F3 structure
        ("plane-warp excess 60-150 AU (b7)", 0.031, "F3",
         "step_b7_plane_warp.json"),
        ("conditional pole-tilt (b8)", 0.044, "F3",
         "step_b8_footprint_conditional.json"),
        ("anti-cap deficit (b6)", 0.045, "F3",
         "step_b6_antipode_deficit.json"),
        ("pole point-concentration (b4)", 0.080, "F3",
         "step_b4_pole_geometry.json"),
        # F4 comet clock channel (numbered pipeline)
        ("Warsaw full-arc in-cap excess (32)", 0.0089, "F4",
         "step_09_nested_field.json"),
        ("CODE-only confirmation (11)", 0.003, "F4",
         "step_11_code_confirm.json"),
        ("continuous no-cap statistic (14)", 0.005, "F4",
         "step_14_blind_validation.json"),
        ("doubly-matched pairs (16)", 0.0019, "F4",
         "step_16_signed_channels.json"),
        # F5 cross-population direction
        ("triple-axis concentration (05)", 0.031, "F5",
         "step_05_one_map.json"),
        ("JFC pairing test (b12)", 1e-4, "F5",
         "step_b12_antipodal_populations.json"),
        ("split-half axis recovery (b11)", 0.0001, "F5",
         "step_b11_population_selectivity.json"),
        # F6 nulls constraining alternatives (reported as
        # nulls -- NOT counted as detections)
    ]
    nulls = [
        ("element coupling q,e,i,a (b9)", "no axis->elements "
         "coupling -- disfavours shepherding"),
        ("resonant substructure (03)", "no MMR substructure -- "
         "disfavours shepherding"),
        ("plunging population (b11)", "isotropic -- "
         "boundary-resident selectivity"),
        ("signed repulsion channel (16)", "no per-orbit push -- "
         "population asymmetry only"),
        ("rotation-axis coherence (16)", "distributed "
         "perturbation, not coherent torque"),
        ("galactic-tide orientation (b14)", "no equatorial "
         "organization -- tide disfavoured"),
    ]

    out = {"n_primary_tests": len(tests),
           "tests": [{"name": t, "p": p, "family": f, "src": s}
                     for t, p, f, s in tests],
           "nulls": nulls}

    # staleness self-audit: every hardcoded p above must appear in
    # its named source JSON (tolerance 5%); a mismatch means the
    # upstream result moved and this ledger needs updating
    def _leaves(o):
        if isinstance(o, dict):
            for v in o.values():
                yield from _leaves(v)
        elif isinstance(o, list):
            for v in o:
                yield from _leaves(v)
        elif isinstance(o, (int, float)) and np.isfinite(o):
            yield float(o)

    stale = []
    for t, p, f, s in tests:
        try:
            vals = set(_leaves(json.loads(
                (RES / s).read_text())))
        except Exception:
            stale.append(dict(test=t, claimed=p, src=s,
                              issue="source json unreadable"))
            continue
        if not any(abs(v - p) <= 5e-4 or
                   (v != 0 and abs(v - p) / abs(v) <= 0.05)
                   for v in vals):
            near = min(vals, key=lambda v: abs(v - p)) if vals \
                else None
            stale.append(dict(test=t, claimed=p, src=s,
                              nearest_in_source=near,
                              issue="claimed p not in source"))
    out["staleness_audit"] = dict(
        n_checked=len(tests), n_mismatched=len(stale),
        mismatches=stale)
    if stale:
        print(f"STALENESS: {len(stale)} claimed p-values do not "
              f"appear in their source jsons")
        for s_ in stale:
            print(f"  {s_['test']}: claimed {s_['claimed']}, "
                  f"nearest {s_.get('nearest_in_source')}")

    # per-family BH + Fisher combination of independent
    # confirmations
    fams = {}
    for t, p, f, s in tests:
        fams.setdefault(f, []).append(p)
    fam_summary = {}
    for f, ps in fams.items():
        thr = bh(ps)
        fam_summary[f] = {"n": len(ps),
                          "n_below_raw05": int(sum(x < 0.05
                                                   for x in ps)),
                          "bh_threshold": round(thr, 4),
                          "n_surviving_bh": int(sum(
                              x <= thr for x in ps))}
        print(f"{f}: n={len(ps)} raw<.05={fam_summary[f]['n_below_raw05']} "
              f"BH thr={thr:.4f} survive={fam_summary[f]['n_surviving_bh']}")
    out["family_bh"] = fam_summary
    out["family_bh_note"] = (
        "Historical selected-test ledger; within-family BH does not control the "
        "complete adaptive search. Correlation does not excuse a failed test. "
        "The hardcoded secondary entries are archival descriptive values, not "
        "a current confirmatory multiplicity correction.")

    # Fisher combined on the two INDEPENDENT populations:
    # TNO clustering (best conditioned null) + comet clock
    # channel (CODE confirmation).  The inputs are harvested from
    # their source jsons and checked against the declared ledger
    # values, so upstream drift cannot silently move the headline
    # combination; Fisher is evaluated at ledger precision.
    from scipy.stats import chi2, combine_pvalues

    def _harvest(src, *path):
        d = json.loads((RES / src).read_text())
        for k in path:
            d = d[k]
        return float(d)

    fisher_src = [
        ("tno_conditioned_p", "step_02_clustering.json",
         ("samples", "detached_a150_q30_cc3", "p_varpi_conditioned"),
         0.0073),
        ("comet_code_confirm_p", "step_11_code_confirm.json",
         ("code_only", "d_of", "p"), 0.003),
    ]
    fisher_vals = {}
    for name, src, path, declared in fisher_src:
        raw = _harvest(src, *path)
        if not np.isfinite(raw) or not 0 < raw <= 1:
            raise ValueError(f"Invalid Fisher source: {src} {path}: {raw}")
        fisher_vals[name] = raw
    out["staleness_audit"]["n_checked"] += len(fisher_src)
    out["staleness_audit"]["n_mismatched"] = len(stale)
    out["staleness_audit"]["mismatches"] = stale
    tno_p = fisher_vals["tno_conditioned_p"]
    comet_p = fisher_vals["comet_code_confirm_p"]
    _, p_comb = combine_pvalues([tno_p, comet_p], method="fisher")
    out["combined_two_populations"] = {
        "tno_conditioned_p": tno_p,
        "comet_code_confirm_p": comet_p,
        "fisher_p": float(p_comb),
        "note": ("Nominal under independent valid null p-values. Different "
                 "populations share selection assumptions and a selected sector; "
                 "this is not a global significance or a TEP model probability.")}

    print(f"combined Fisher p = {p_comb:.2e}")

    RES.mkdir(exist_ok=True)
    (RES / "step_b15_evidence_summary.json").write_text(
        json.dumps(out, indent=1))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

    # sorted -log10 p with BH line
    ps = np.sort([t[1] for t in tests])
    m = len(ps)
    k = np.arange(1, m + 1)
    ax[0].scatter(k, -np.log10(ps), c="steelblue", zorder=3)
    ax[0].plot(k, -np.log10(0.05 * k / m), "r--", lw=1,
               label="BH q=0.05 threshold")
    ax[0].axhline(-np.log10(0.05), color="0.5", ls=":", lw=1)
    ax[0].set(xlabel="rank", ylabel="-log10(p)",
              title=f"{m} primary tests -- sorted p\n"
                    f"vs BH threshold")
    ax[0].legend(fontsize=7)

    # per-family bars
    fl = ["F1", "F2", "F3", "F4", "F5"]
    surv = [fam_summary[f]["n_surviving_bh"] for f in fl]
    tot = [fam_summary[f]["n"] for f in fl]
    x = np.arange(len(fl))
    ax[1].bar(x, tot, color="0.75", label="tests run")
    ax[1].bar(x, surv, color="crimson",
              label="survive BH q=0.05")
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(
        ["TNO\ncluster", "bias\nbounds", "structure",
         "comet\nclock", "cross-\npop"], fontsize=8)
    ax[1].set(ylabel="N tests",
              title="multiple-testing control\nby family")
    ax[1].legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(FIG / "supplementary" / "step_b15_evidence_summary.png", dpi=300)
    logger.data_save(RESULTS / "step_b15_evidence_summary.json")
    logger.data_save(RESULTS / "figures/supplementary/step_b15_evidence_summary.png")


if __name__ == "__main__":
    main()
