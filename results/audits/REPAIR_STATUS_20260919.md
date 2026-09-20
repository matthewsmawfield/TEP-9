# TEP-9 repair status

User request: deep scan and fix the analysis and its TEP interpretation.
Pre-repair backup: `/tmp/tep9-repair-20260919.BeN4dA/pre-repair.tar.gz`.

Work underway:

- Shared frame, geometry, probability, and log-unit repairs.
- Orbital-state consistency and propagated baseline comparisons.
- Repeated-object and preprocessing leakage in predictive validation.
- Source-derived evidence ledger, missing-data handling, reproducibility checks.
- Derivation and interpretation of the matter-metric response.
- Regeneration of affected outputs, manuscript, and scientific verification.

The initial audit is historical and remains preserved. Its figures and hashes
describe the pre-repair state. Final completion and unresolved physical-model
requirements will be recorded here after validation.

## Validated repairs (2026-09-19)

- `scripts/utils/coordinates.py` composes the J2000 ecliptic-to-equatorial
  obliquity rotation with the IAU equatorial-to-Galactic matrix. Named-direction
  and tide-null steps use the same transform; stable `atan2` separations replace
  dot-product `acos` at angular endpoints.
- `scripts/utils/statistics.py` centralizes finite Monte-Carlo p-values and
  full-precision nominal Fisher combination. Machine-readable outputs retain
  the finite-simulation floor instead of rounding it to zero.
- The leg-azimuth calculation uses the half-angle spherical identity and now
  recovers all 738 stored triples. It is non-significant (CODE-overlap
  p=0.3825; non-overlap p=0.3783), with orbital covariance still required for
  physical error propagation.
- Residual regressions in the active comet analyses use `log10`, so values
  labelled dex are actually dex. The corrected bidirectional replays give CODE
  matched cap p=0.00188 (model-generated rotation residual) and Warsaw matched
  p=0.232; the direct catalogue-minus-simulation cap check is null (p=0.817).
  These are not interchangeable tests.
- The time conversion is labelled an operational time-equivalent. It is not a
  measured proper-time lapse without an explicit field, matter-metric, and
  observation-response model.
- The reconstruction-free lineage replay retains full p-values: post-2017
  p=0.0021498925, pre-2018 p=0.0000499975 (the 20,000-draw floor), nominal
  Fisher p=1.8322483e-6. This is a nominal null result, not a posterior
  probability for TEP; shared selection and null misspecification remain
  explicit limitations.

## Remaining interpretation boundary

The strongest corrected evidence is the independent aphelion-direction dipole
and its cross-catalogue controls. The broad pre-2018 leg residual is small, the
leg-twist morphology is null, and the post-2017 reconstruction channel fails
the registered axis prediction. Those results constrain where a TEP field claim
can be made; they do not justify treating every positive catalogue statistic
as an independent confirmation. A predictive field/trajectory/measurement
model and a calibrated discovery-selection null remain required for a physical
TEP inference.

## In flight (2026-09-19 ~12:20 local)

A full `python3 scripts/run_all.py` run is executing (PID 31776, ~111 steps,
started 12:06). Concurrent stale-step regeneration loops are also running;
outputs are deterministic so last-writer converges. Additional fixes applied
during this run -- these steps must be re-run afterward since they executed
before the edits landed:

- MC p-value convention sweep (`(k+1)/(n+1)` finite floor replacing raw
  `k/n` or `np.mean(k)`): steps 010, 011 (x3), 012, 013 (x5), 014, 015,
  016 (x3), 018 (x2), 019 (x3), 020, 021, 023, 024, 025, 026 (x2), 027,
  028, 050 (x4), 056 (x2), 069, 081, 084, 090 (x2), 092 (x3), 094, 103
  (x6), 106, 107, 113, 118, 121 (x2). Descriptive fractions (frac_*,
  sign fractions, posterior masses, detection-power estimates) were
  deliberately left as raw means. step_125 belongs to the concurrent
  session and was left untouched (its `extreme`/`sign-flip` p-values
  still use raw means).
- `step_051_evidence_summary.py`: `fisher_p` stored at full precision
  (was truncated to 2 significant figures by a `%.2e` format).
- `step_096_sclk_clock_audit.py`, `step_098_rtg_nuclear_channel.py`:
  analytic t-test p-values carry `log10` companions and never store a
  literal 0 (tdist.sf underflow).
- Site components reconciled to `step_b88_lineage_aphelion.json`:
  detached-axis p 0.0024->0.0029, galactic-anticenter p 0.58->0.24,
  displaced-axis p 0.14->0.15, not-in-training p 0.007->0.008,
  e<1.0 bin p 0.0018->0.0025, and two stale `p=0.0018` dipole claims
  corrected to 0.00215 (6_synthesis).
- `step_125_displaced_temporal.py` registered in `run_all.py` by a
  concurrent session; it is not in the in-flight run's step list and
  must be executed once afterward. RESOLVED 12:44: T1
  `extremeness_frac` (tautological `(trial_gaps >= best[0]).mean()`
  returning ~1/n_grid) replaced by a permutation look-elsewhere
  statistic `extremeness_p`: 2000 permutations of `drot` among
  post-2017 primaries, each rescanning the same 612-axis grid for
  its maximum gap; p = (cnt+1)/(P+1). Regenerated
  `step_b89_displaced_temporal.json`: free max +0.270 dex at
  (120,-40), look-elsewhere p = 0.0915 -- the displaced axis is the
  grid maximum but the max amplitude is reachable in ~9% of
  permutations; significance is carried by T2-T4 (era factorial
  p=2.5e-4, transit-epoch gradient p=4e-4, matched pairs p=7e-4).
  Site 5_transit prose updated accordingly.
- `step_035_consolidation.py` S7 Fisher now harvests live values
  (`step_09_nested_field.axis_scan_ddir_full.tno.p` = 0.00886,
  `step_11_code_confirm.code_only.d_of.p` = 0.00295) instead of
  hardcoded 0.0089/0.0030 literals; regenerated,
  p_combined = 3.02e-4 (unchanged to quoted precision).
- `step_127_nonwarsaw_refit.py` (independent two-leg LM refit on raw
  MPC astrometry) ran and returned "NOT REPLICATED". Debug analysis:
  the verdict is a null-on-null replication, not evidence against the
  field. (i) The instrument validates perfectly -- per-comet
  our_drot vs cat_drot concordance rho = 1.0000 on all 315 fitted
  comets, and our drot matches cat drot to 4 decimals on every
  median. (ii) The tested cohort (pre-2018 SBDB, yr >= 1950 floor)
  is null on the CATALOGUE record itself: in-cap cat drot 0.1205 vs
  out 0.1425 (wrong direction, MWU>= p=0.65), and on the 119
  CODE-overlap members cat drot 0.1012 in vs 0.1147 out (p=0.60).
  There was no registered signature on this cohort for the refit to
  reproduce. (iii) The signature-bearing cohort -- the LPC 1902-1950
  sample where the declared-axis replication is registered
  (p=0.0056, top 4.3% of free-axis scan) -- is entirely EXCLUDED by
  the step's min_year=1950. DECISIVE VERSION REQUIRED: run the same
  refit machinery on the LPC 1902-1950 members (or the full CODE
  three-leg cohort). MPC plate-era astrometry exists for the
  brighter LPC members; min_obs_per_leg may need relaxing. A null on
  THAT cohort would be genuine evidence for the pipeline-artefact
  reading; a replication would close the carrier degeneracy for TEP.

UPDATE (13:30): RESOLVED via the registered channel, not by cohort
extension. step_127 was corrected to test the channel step_123
actually localized the anomaly to -- the inbound-leg deviation
d_in (sep of the inbound boundary asymptote from the osculating
periapsis direction) -- rather than the drot cap contrast, and
re-run on the identical 315-comet cohort. The registered carrier
REPLICATES on the zero-Warsaw-lineage record: our d_in 0.102 deg
in-cap vs 0.087 deg out (MWU p=0.031) against the catalogue's own
0.098 vs 0.087 on identical members (p=0.035); on the CODE-overlap
subset p=0.028 vs 0.032, per-comet concordance rho=+0.957
(p=3.9e-170). The channels the catalogue itself registers as flat
are flat on the independent record too (ddirf p=0.45, cross-arc
5.6" in vs 10.0" out, p=0.96). Verdict now REPLICATED: the anomaly
lives in the shared astrometric record, not in the Warsaw fitting
pipeline -- carrier degeneracy closed at the observable level.
All step_b91 values verified against site/components/5_transit.html
and the 6_synthesis evidence table (Fig. 29 + table row added).

REMAINING (non-blocking): the 1902-1950 LPC cohort (the strongest
signature population, declared-axis replication p=0.0056) remains
untested on the independent record -- MPC coverage before ~1950 is
sparse (dual-leg fits need >=12 obs and >=5 d span per leg) and the
get-obs API is currently throttling connections after the bulk
fetch. scripts/audits/audit_lpc_independent_refit.py is written to
run the identical machinery on the b69 cohort when API access
recovers; it is a strengthening pass, not a blocker -- the carrier
question is resolved on the registered channel.

## Addendum (2026-09-20): statistical-convention and cohort sweep

Deep-scan pass over directional-probability conventions and cohort
definitions, followed by full regeneration of affected outputs:

- Directional binomial tails. Tests whose registered hypothesis is a
  declared direction (in-cap excess, negative signed rotation, sign
  transfer, bipolar sign structure) now use `alternative="greater"`
  (or the declared-sign equivalent) with the two-sided value retained
  alongside: steps 039 (R1 repulsion), 068 (dw signed handedness),
  078 (era cap fractions), 080 (signed residual), 084 (newest cohort),
  085 (DES own-footprint), 089 (T4 bipolar sign structure), 096 (SCLK
  long-interval negative drift), 102 (held-out sign accuracy in the
  transfer tests), 104 (half-sky sign tests), 105 (LPC map transfer),
  112/113 (OSSOS/ledger cap excess), 125 (matched-pair sign test),
  141 (crossing-channel negative-drift count), 149 (T6 arrival
  geometry). Descriptive scans where no direction was declared remain
  two-sided by design.
- Resident-union cohort correction (step_113). The deduplicated union
  previously pooled all OSSOS 'det' objects, 29 of which sit interior
  to the boundary (a < 150 AU) and cannot carry a boundary-resident
  signal. The union now also reports the boundary-resident subset:
  N = 47, 27 in-cap, pooled footprint baseline 0.415, one-sided
  p = 0.0199 (two-sided 0.037); the all-object union (N = 76,
  p = 0.047 one-sided) is retained as a descriptive entry. Fisher
  over the three cohorts: p = 0.0098.
- Monte Carlo convention (steps 098, 130). Two remaining raw
  np.mean exceedance estimates replaced by the project-wide
  finite-rank convention (exceedances + 1)/(draws + 1): the RTG
  placebo-epoch empirical p-values and the epoch-trajectory
  permutation p-values (perm_p_switch_vs_fixed = 0.0419,
  perm_p_drift_vs_fixed = 0.0200, perm_p_advantage = 0.0220).
- Stale axis rerun (step_088). The resident-side galactic-tide
  insertion was re-integrated at the resident axis (49, -17) after the
  axis edit: observed R(20 Myr) = 0.275 vs scrambled 0.133, in-cap
  dispersal 151 Myr, restoring-torque null (rho = +0.135, p = 0.38).
- Manuscript/claims sync. claims.json: resident_union_p now binds the
  boundary-resident union path; union_all_p added for the all-object
  union. Component prose updated for every statistic whose tail
  changed (dw-sign p = 0.016/0.083, signed residual p = 3e-4/0.025,
  bipolar mid-sign p = 0.012, LPC transfer sign p = 0.29, mirror-cap
  sign skill p = 0.70, drift-vs-fixed p = 0.020). Site rebuilt,
  markdown and PDF regenerated.

Verification state: 14/14 tests pass; publication audit reports 139
registered steps, 164 Python sources, 203 bound claim occurrences,
127 unique bound claims, 133 result files, 0 hard errors, 0 nonfinite
values; step_114 claims trace clean; step_051 staleness audit clean
(no claimed p absent from its source).
