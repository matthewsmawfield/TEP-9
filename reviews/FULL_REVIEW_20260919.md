# TEP-9 full manuscript and analysis review

> Scope correction: `CORPUS_SCOPE_CORRECTION.md` supersedes this review's standalone treatment of TEP field existence. TEP-9 imports that premise from the companion programme. The numerical repairs and application-specific limitations below remain relevant.

Review date: 19 September 2026. Scope: the ten manuscript components, the registered pipeline, numerical utilities, stored results, publication build and selected numerical reruns. The accompanying validation manifest records the exact reviewed files. This is not a claim that all 128 steps were regenerated from a clean data acquisition, nor that every historical result has been scientifically validated.

## Outcome

The rewritten paper makes its strongest defensible case around a shared **spatial** structure: the resident angular concentration and the comet-aphelion dipole. It separates those observations from a physical identification with TEP. The available calculations do not prove TEP: the full conventional planetary propagation reproduces the catalogue rotations, the original temporal rotation prediction fails its holdout, and the conversion of angle to equivalent orbital time is not an observation of proper time.

These are substantive distinctions. Stronger causal wording cannot repair them. A convincing physical case requires a specified TEP response in the actual observations, a comparison with conventional alternatives under calibrated selection and uncertainty, and successful withheld-data prediction. Section 6 now supplies the mathematical observation-model structure and a concrete validation protocol.

Original manuscript components are preserved in `pre_review_components/`. They are superseded, not additional evidence. The reviewed paper is built from `site/components/`; headline numbers are bound to exact result fields.

## Evidence that survives, and its limits

| Observation | Reviewed result | Interpretation |
|---|---|---|
| Secure detached resident sample, 44 objects | Longitude resultant 0.3316; conditioned-null p = 0.00734963 | Concentration under the specified orbit-conditioned angular reference. A full survey simulator is still needed. |
| Deduplicated resident union, 76 objects | p = 0.0802688 against the adopted footprint baseline | The broader union is not significant at 0.05. Overlapping catalogue memberships cannot be counted as independent replications. |
| CODE-only deep plungers, 54 objects | In/out median rotation 0.212024/0.130745 degrees; p = 0.00295055 | Directional population contrast in a quantity reproduced by conventional dynamics. |
| Full conventional prediction versus catalogue | Direct residual cap p = 0.817444; reduced-regression proxy p = 0.00187695 | The proxy is not a residual after subtraction of the full dynamical prediction. |
| Original temporal rotation holdout | Medians 0.113920/0.171153 degrees; greater-than p = 0.955560 | Opposite ordering to the declared positive-excess prediction. The later displaced-axis search is a revised hypothesis. |
| Spatial aphelion dipole | Pre-2018 d = 0.1101; post-2017 d = 0.0949 | Persistent spatial target. The post-2017 Galactic and ecliptic references give p = 0.00214989 and 0.00224989. |
| Modern displaced residual by catalogue stratum | Uniform-modern p = 0.283223; mixed p = 6.45262e-7 | A material catalogue-composition confound remains. This contrast does not by itself prove a particular instrumental explanation. |
| Retrospective geometric comparison | Orbit-block joint p = 0.0775461; mirror p = 0.0733463 | Three directions from one ISO are treated as one orbit, not three independent discoveries. Geometry was not genuinely held out. |

These are component probabilities under their stated nulls, not posterior probabilities for TEP. The resident a > 150 au axis and the nested a > 250 au/comet-sector axis are now distinguished, as are longitude sectors and spherical caps. Selection assumptions, failures, catalogue overlap and choice of observable are stated next to the results they affect.

## Repaired numerical and inferential issues

1. **Cross-channel duplication and invalid combinations.** Step 135 previously used two geometric terms from the same GNSS-II direction and treated a product as a joint probability. It now uses one GNSS-II contribution and a proper Fisher reference. The external meridian product 0.00919981 corresponds to nominal Fisher p = 0.0523338. The four-channel Fisher reference is 8.97463e-9, explicitly uncalibrated for retrospective selection and channel dependence. It is not presented as global TEP significance.

2. **Silent substitution of the strongest possible simulated result.** The comet contribution in step 135 read a nonexistent JSON key, caught the exception and substituted 5e-5. It now reads the actual `T2_joint_opposite_polarity/p_J` result and fails when that input is missing. Its current numerical value happens to lie at the simulation floor; that coincidence did not make the fallback valid.

3. **Retrospective geometry described as prediction.** Step 139 now labels its geometry retrospective, excludes the resident direction used in constructing the sector from external validation, and excludes the MGEX direction flagged non-identifiable by its source. It rotates all directions belonging to each ISO together, preserving their internal geometry. The old joint value near 0.024 increases to 0.07755 under the corrected block comparison. A closest-of-two-cones tail now integrates the union of both eligible bands, including overlap.

4. **Selected components inside a purported global significance.** Step 126 selects objects in a cap before applying an unrestricted angular reference. Its random-axis landscape compares different directions on the same observed data; it does not generate null data and rerun the entire axis/model search. Its score and rank are retained as descriptive results, with that calibration limit made explicit. A genuine global test must simulate the selection and complete search procedure.

5. **Finite Monte Carlo and directional selection.** Steps 130–132 now use finite-simulation ranks `(k + 1)/(B + 1)`. Step 129's direction-selected one-sided combination is replaced with its two-sided reference. Step 051 reads source probabilities without a fabricated fallback and qualifies the historical multiple-testing family. Exact probabilities remain unrounded until display or combination.

6. **Orbit-fitting implementation.** Propagation sorted the requested epochs but did not restore output to the caller's order. It now returns states in the original requested order. A differential-correction iteration could retain a worse or nonfinite proposal after exhausting damping attempts; it now keeps the last accepted state. Residual-sign documentation was corrected. Nine stratified cached real-comet replays, spanning the pre-2018, post-2017 and historical-signature cohorts and observation-count quantiles, reproduce earlier RMS and rotation values within 1e-9. This is a regression check, not full-cohort regeneration or covariance validation.

7. **Machine verdicts stronger than their tests.** Steps 127, 128 and 134 now identify exploratory refit comparisons and expose measured comparison flags. Step 141 computes its sign statistic from the observed crossing counts rather than hard-coded counts. SCLK coefficients are described as correlation/calibration information, not as direct independent measurements of anomalous oscillator rates; dependent probes are not certified by a nominal sign combination.

8. **Physical identification.** The paper no longer identifies a reduced dynamical regression residual with a missing gravitational rotation or proper-time lapse. It explains why Keplerian phase displacement preserves periapsis direction, why a conservative deflection need not produce a net asymptotic energy change, and why an exact scalar differential has zero closed-loop integral. A linearized nuisance-fit projection shows which part of a specified TEP signal could survive orbit and calibration fitting. These corrections preserve testability while removing unsupported uniqueness claims.

The review also verifies earlier local fixes to the J2000 coordinate transform, small-angle leg-twist calculation, probability precision and log-unit descriptions. Regression tests compare the frame with SPICE and protect the small-angle geometry. Those corrections were already present when this review began and are not represented as newly discovered changes.

## Publication and reproducibility repairs

- Replaced permissive numeric-pool matching with 34 exact headline-claim bindings, used at 43 locations. Each binding names its result file and JSON path. Missing, nonfinite or invalid probability inputs fail the renderer; unrelated nearby numbers cannot satisfy a claim.
- Replaced step 114 with the exact publication audit and moved it to the end of the 128-step registry. The audit inventories source hashes, results, runtime versions and registration coverage. It does not infer provenance merely from a filename appearing in another script or from a modification time.
- Made HTML and Markdown use the same binding renderer; removed obsolete flyby-result injection. Updated the manuscript metadata, catalogue citation, data-availability account and README instructions.
- Fixed the PDF print stylesheet: an HTML `<style>` wrapper had been passed where raw CSS was required. Improved print scale, page flow, figure captions and bibliographic metadata; synchronized generated PDF download copies.
- Replaced the obsolete `live-server` dependency chain with a local Node HTTP preview and file watching. The final npm audit reports zero known vulnerabilities. Added the directly imported YAML dependency. The preview's content and traversal checks pass.
- Recorded active runtime versions separately from declared dependency pins. A reproducible environment must use a recorded lock/snapshot; the old requirements are not claimed to describe the review's actual interpreter.
- Preserved pre-existing worktree changes. No destructive reset, external publication, push or commit was performed.

## Verification scope

The publication audit checks 148 Python sources, all 128 registered entries, 122 stored step-result JSON files, and all 43 bound claim occurrences. Its reviewed run reports zero hard errors and zero nonfinite leaves in those step-result files. Passing this audit establishes structural consistency and exact claim sourcing, not scientific proof.

The focused reruns cover steps 010, 034, 051, 113, 124, 126, 129, 130, 131, 132, 135, 139, 140 and 141; these include the resident result, CODE contrast, lineage dipoles, corrected geometry and revised probability combinations. Refit summary runs for 127 and 128 used existing checkpoints. Step 134's final corrected summary is rerun separately. The observable audit independently checks the dynamical-residual distinction, coordinate sensitivity, simulation precision and phase interpretation. Fourteen regression tests pass. Nine fresh real-data refits agree with the earlier checkpointed values to 1e-9.

A concurrent Devin process began full `--redo` runs for steps 127 and 128 against this same directory during final validation. The review's duplicate checkpoint-writing process was stopped. The external runs were left running. Consequently the reviewed publication is identified by a frozen copy of its bound result inputs and their hashes, rather than by a claim that every subsequently changing file belongs to the reviewed snapshot. Those runs may require rebuilding downstream summaries and the paper when they finish. Refit checkpoints written during the overlap may contain repeated designations; the readers deduplicate by designation, but a clean fresh execution with exclusive ownership is preferable for the next provenance record.

The final `review_validation.json` records artifact hashes, runtime inventory, claim checks, regression output, observable-audit output and the frozen bound inputs. See `validated_snapshot/` for the reviewed inputs and publication source. This is a publication snapshot, not a full raw-data archive.

## Remaining requirements for a decisive physical claim

1. Specify a finite TEP field/response model that predicts the measured astrometry or clock observations, including how conventional nuisance fitting absorbs its signal.
2. Fit conventional and TEP alternatives through the same survey-selection and observation pipeline, with full covariance and plausible catalogue/station/epoch systematics.
3. Demonstrate integration and boundary-extraction convergence, injection recovery, uncertainty coverage and false-positive calibration for the proposed effect size.
4. Freeze the model, axis, cuts, sign, primary statistic and stopping rule before observing a new validation set. A revised bipolar/moving model needs its own successful validation; it cannot inherit the failed original prediction's status.
5. Calibrate the entire analysis search and cross-channel dependence before quoting a global significance. A small nominal Fisher reference does not perform this calibration.

These are unresolved scientific tasks, not editorial omissions that can be closed by stronger language. The reviewed paper now presents the spatial result clearly, carries its negative controls alongside it, and defines what further evidence would turn the proposed TEP explanation into a discriminating test.
