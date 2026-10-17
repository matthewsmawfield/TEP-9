# Temporal Equivalence Principle: A Proper-Time Domain Boundary in the Outer Solar System

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22858190.svg)](https://doi.org/10.5281/zenodo.22858190)
![Temporal Equivalence Principle: A Proper-Time Domain Boundary in the Outer Solar System](./site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Addis Ababa)  
**First published:** 30 September 2026 · **Last updated:** 17 October 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.22858190](https://doi.org/10.5281/zenodo.22858190)  
**Website:** [https://mlsmawfield.com/tep/9/](https://mlsmawfield.com/tep/9/)  
**Paper Series:** TEP Series · Group G (VII. Domain-Boundary & Observational Anomaly Case Studies) · Paper 36 (`G36`)

## Abstract

<!-- Generated from site/components/1_abstract.html; update the manuscript source, then run tools/sync_current_corpus_readmes.py. -->

Adopting the Temporal Equivalence Principle (TEP) field framework established in the companion corpus, this paper tests an outer-solar-system realization: a localized proper-time domain boundary organizing detached trans-Neptunian residents and long-period comet transits about the same 60° sector of inertial sky. The two populations probe complementary secular and crossing responses.

The detached resident population clusters in perihelion longitude ($R_\varpi =$ 0.332, conditioned-null $p =$ 0.00735), and a Keplerian discovery-selection model retains a directional excess on the 52 class-1 objects within its stated visibility limits ($p =$ 0.0092). In the transit population, the primary coordinate-level observable — spatial aphelion anisotropy computed directly from osculating elements with no three-leg reconstruction — leans toward the shared sector across both historical and modern eras ($d_\parallel = +0.110$, $p < 10^{-4}$ pre-2018; $d_\parallel =$ +0.0949, $p =$ 0.00215 post-2017; Fisher-combined $p = 1.8\times10^{-6}$; full CometEls census $d_\parallel = +$0.10, $p =$ 2×10⁻³), remaining flat across all orbit-quality and arc-length strata. As a secondary, solver-sensitive observable, the CODE three-leg transit record carries a direction-dependent orbit-reconstruction rotation ($p =$ 0.00295), retained after planetary-encounter and observing-leverage controls. Independent two-leg fits of raw Minor Planet Center astrometry for 585 pre-2018 comets reproduce catalogue rotations ($\rho =$ +0.998) and retain the inbound-leg excess ($p =$ 0.0370; CODE-overlap $p =$ 0.00688), with no resolved in-cap energy excess. Under crossing-localized transport, the angular residual corresponds to an equivalent-time contrast of $+$4.9 yr; this is a model-conditioned conversion rather than a direct clock reading.

A ten-channel synthesis across six catalogue lineages gives $S_{10} =$ 186.2 at the resident axis, with conditional sky-area rank 0.0145 on the fixed observed catalogues. This measures directional convergence; it is not a global discovery probability. While the coordinate-level aphelion dipole persists unchanged in the post-2017 holdout, the estimator-dependent reconstruction channel reverses its declared-axis carrier and exhibits a displaced residual structure under short-arc solver absorption, restricting the simple stationary three-leg interpretation. Point-mass, tidal, drag and reduction explanations are tested through quantitative amplitude, localization and transfer-function controls. The canonical disformal benchmark specifies the required field excursion using actual conic entry and exit geometry; its source profile and full metric-to-observation response remain to be derived.

The results provide mutually reinforcing evidence for a shared direction-organized outer-system structure; its interpretation as an outer-system TEP domain boundary remains model-conditioned and testable through the prospective discoveries registered below. The registered prospective test predicts 59% of new detached objects in the cap, against a 42% characterized-footprint expectation, requiring approximately 123 discoveries for 95% power at a 1% test. This separates the present evidence from the next key validation.

Key words: Kuiper belt: general – comets: general – gravitation – celestial mechanics – astrometry

Current manuscript source: [1_abstract.html](site/components/1_abstract.html).

## Reproduce and validate

From the repository root:

```bash
python3 -m pytest -q tests
python3 scripts/audits/audit_tep9_observables.py
python3 scripts/audits/audit_publication.py
cd site
npm ci
npm run build:markdown
```

Generate the PDF with `python3 scripts/generate_site_pdf.py` from the root. The Python analysis requires the dependencies in `requirements.txt`; the PDF converter additionally uses Playwright and a Chromium installation. Active review versions are recorded by the publication audit and may differ from the declared requirements. Cross-channel steps read sibling TEP-GNSS-II, TEP-GNSS-MGEX and TEP-LLR result files, so a standalone checkout does not reproduce those external pipelines.

For the full registered pipeline:

```bash
python3 scripts/run_all.py
```

There are 152 registered steps. Acquisition steps can contact live services and change the raw-data snapshot. Preserve the current data and provenance when reproducing this version. Later steps may also fetch missing cached data. The full run is computationally expensive and is distinct from the focused review reruns listed in the review report. The exact publication audit runs last.

## Evidence and numerical safeguards

- `site/claims.json` binds headline numbers to specific JSON fields; missing, invalid or nonfinite values fail the publication build.
- Monte Carlo ranks retain the observed outcome: `(exceedances + 1)/(draws + 1)`; no rounding precedes a combination.
- `scripts/utils/coordinates.py` defines the verified J2000 frame transforms and stable small-angle geometry.
- The revised geometric comparison preserves the three directions of each ISO as one orbit block. It is explicitly retrospective.
- Cross-channel Fisher references use one GNSS-II direction term; catalogue overlap and selection are not claimed to be calibrated away.
- The ten-channel random-axis landscape describes the observed sky concentration. It is not a null-data maximum-statistic calibration.

## Data and outputs

Public catalogue observations are under `data/raw/`; source-specific provenance files record retrieval details and available hashes. The empirical record is separate from simulated nulls, bootstrap samples, orbit integrations and numerical validation controls. CODE exports select the published preferred solution and require matching model identities across the orbit tables; retrieval provenance records their hashes. MPC astrometry and SBDB seeds are cached for the refits.

Results are under `results/`, figures under `results/figures/`, and logs under `logs/`. `results/step_b78_claims_trace.json` records the current repository inventory, exact publication bindings and runtime versions. It validates traceability and structure, not every inference in every historical result file.
