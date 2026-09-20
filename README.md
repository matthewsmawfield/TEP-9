# TEP-9: outer-Solar-System domain-boundary tests

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22858191.svg)](https://doi.org/10.5281/zenodo.22858191)

**DOI:** [10.5281/zenodo.22858191](https://doi.org/10.5281/zenodo.22858191)

This repository develops the outer-Solar-System application of the Temporal Equivalence Principle (TEP), taking field existence and universal matter coupling from the companion programme. It investigates the domain geometry and linked responses of resident trans-Neptunian objects and transiting comets.

The resident angular concentration and the comet-aphelion dipole isolate the same 60-degree sector of inertial sky under measured discovery-footprint and tide-aware nulls. An independent two-leg refit of the raw Minor Planet Center astrometry — built without any Warsaw-lineage input — reproduces the per-comet boundary rotation at rho = +0.998, and the post-2017 prospective cohort's dominant structure replicates as a lapse-slip anomaly of the identical signature class at the opposite polarity of a bipolar axis anchored near the CMB rest frame. A ten-channel global synthesis spanning six catalogue lineages (SBDB, DES, MPCORB, CODE, Warsaw, CometEls) yields an omnibus statistic S_10 = 187.2 (p = 0.0118 against a 20,000-draw random-axis permutation null). These multi-lineage data provide mutually reinforcing evidence for a non-integrable dynamical time field in the outer solar system. See `reviews/CORPUS_SCOPE_CORRECTION.md` and `reviews/COMPANION_PREMISES.json` for the imported results and exact manuscript versions.

The manuscript source is `site/components/`. The full review and unresolved scientific requirements are in `reviews/FULL_REVIEW_20260919.md`. Pre-review manuscript components are preserved under `reviews/pre_review_components/` and are superseded.

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

There are 141 registered steps. Acquisition steps can contact live services and change the raw-data snapshot. Preserve the current data and provenance when reproducing this version. Later steps may also fetch missing cached data. The full run is computationally expensive and is distinct from the focused review reruns listed in the review report. The exact publication audit runs last.

## Evidence and numerical safeguards

- `site/claims.json` binds headline numbers to specific JSON fields; missing, invalid or nonfinite values fail the publication build.
- Monte Carlo ranks retain the observed outcome: `(exceedances + 1)/(draws + 1)`; no rounding precedes a combination.
- `scripts/utils/coordinates.py` defines the verified J2000 frame transforms and stable small-angle geometry.
- The revised geometric comparison preserves the three directions of each ISO as one orbit block. It is explicitly retrospective.
- Cross-channel Fisher references use one GNSS-II direction term; catalogue overlap and selection are not claimed to be calibrated away.
- The ten-channel random-axis landscape describes the observed sky concentration. It is not a null-data maximum-statistic calibration.

## Data and outputs

Public catalogue observations are under `data/raw/`; source-specific provenance files record retrieval details and available hashes. The empirical record is separate from simulated nulls, bootstrap samples, orbit integrations and numerical validation controls. CODE is a hash-pinned archival export. MPC astrometry and SBDB seeds are cached for the refits.

Results are under `results/`, figures under `results/figures/`, and logs under `logs/`. `results/step_b78_claims_trace.json` records the current repository inventory, exact publication bindings and runtime versions. It validates traceability and structure, not every inference in every historical result file.
