# TEP-9 observable and analysis audit — 19 September 2026

This focused audit confirms numerical errors and a material gap between the
computed orbital statistics and their interpretation as proper-time evidence.
One correction substantially strengthens a nominal spatial statistic. Other
corrections recover weak results or narrow what an existing result establishes.
They do not establish that every weak channel is a hidden TEP detection.

The repair pass changed the shared geometry/statistics utilities and the
affected source steps, then regenerated the directly affected result files.
This report is still a focused audit, not a complete review or rerun of every
step in the 111-step pipeline.

Reproduce from the repository root:

```bash
python3 scripts/audits/audit_tep9_observables.py
```

The companion `audit_tep9_observables.json` records input hashes,
full-precision results, coordinate consumers, and validation details.

## 1. A rounding bug discards the strongest era result

In `scripts/steps/step_124_lineage_aphelion.py:169`, Monte Carlo p-values
are rounded before subsequent calculations. The pre-2018 p-value is
`1/20001 = 0.00004999750012499375`, corresponding to zero exceedances in
20,000 draws. Rounding to four decimal places makes it zero.

At line 376, the Fisher combination removes values outside `0 < p < 1`,
silently dropping this era. The advertised two-era result therefore contains
only the post-2017 result. This is a demonstrable source of artificially weak
reported significance.

Replaying the original null function and random-call order reproduces the
stored primary p-values after rounding. Keeping the original precision yields:

| Calculation | Post-2017 p | Pre-2018 p | Nominal Fisher p |
|---|---:|---:|---:|
| Stored result | 0.0018 | 0, discarded | 0.0018 |
| Full precision, original frame | 0.00179991 | 0.0000499975 | 0.00000154997 |
| Full precision, corrected frame | 0.00214989 | 0.0000499975 | 0.00000183225 |

The repair retains full precision and exceedance counts in machine
results, round only display strings, keep all valid component tests, and
reject invalid p-values explicitly. A valid p-value of one must also remain
in the combination, with its zero contribution and its degrees of freedom.

These are **nominal** combined p-values under the existing nulls and Fisher's
independence assumptions. The older cohort reaches the Monte Carlo resolution
floor. Shared selection effects are not removed by using disjoint objects.
Neither number is a calibrated probability that TEP is true.

## 2. Equatorial coordinates are being used as ecliptic coordinates

The matrix called `ECL2GAL` in `scripts/utils/tep9_common.py:31` is the
equatorial-to-Galactic matrix. Applying it directly to the ecliptic orbital
vectors omits the obliquity rotation. The same literal occurs in 16 source
files, including the shared utility; the JSON enumerates them. Indirect
consumers and cached downstream results extend the affected scope.

Two concrete reference directions demonstrate the error:

| Direction | Current purported ecliptic longitude/latitude | Correct J2000 mean ecliptic longitude/latitude | Angular displacement |
|---|---|---|---:|
| Galactic centre | 266.405°, −28.936° | 266.840°, −5.536° | 23.403° |
| North Galactic pole | 192.859°, +27.128° | 180.023°, +29.811° | 11.591° |

The audit constructs the reference transform using Astropy's
`BarycentricMeanEcliptic(equinox="J2000")` and independently checks it against
the original matrix composed with the missing obliquity rotation. Agreement
is within `1e-6` per matrix entry; smaller differences reflect frame and
obliquity conventions. See the [Astropy coordinate-transform documentation](https://docs.astropy.org/en/stable/coordinates/transforming.html).

Correcting this changes the post-2017 primary dipole p-value from about
0.00180 to 0.00215 for the same random seed. The dipole amplitude itself,
`d_parallel = 0.09489578`, is unchanged: it is evaluated in ecliptic space.
The frame bug therefore does not erase this spatial result, but Galactic
control directions, tide orientations, and related claims require reanalysis.

The repaired calculation retains the original absolute-latitude resampling,
random hemisphere assignment, and uniform-longitude assumptions. Fixing the
rotation does not validate those assumptions as a survey-selection model.

## 3. The leg-twist diagnostic is numerically recoverable

In `scripts/steps/step_123_leg_coherence.py:181`, the cutoff
`sin(d_in)*sin(d_out) < 0.01` rejects every one of 738 comets. Their actual
denominators range from `3.36e-8` to `6.68e-4`. The claimed universal
numerical degeneracy is a consequence of this cutoff.

Using the equivalent half-angle expression avoids subtraction of cosines
close to one. Every stored angle triple is finite and geometrically
admissible. The recovered cosines agree with an independent long-double
replay to `1.51e-9`; the residual is set by decimal precision of the stored
angles in the near-singular regime.

The recovered one-sided twist comparison gives:

| Cohort | N | In-cap N | p, in-cap cosine smaller |
|---|---:|---:|---:|
| CODE overlap | 288 | 87 | 0.38249 |
| Non-overlap | 450 | 107 | 0.37831 |

This repairs a discarded test but does not recover a significant twist.
Numerical precision on stored angles is distinct from observational
uncertainty; covariance propagation is still needed to determine how well
the physical azimuth is measured.

## 4. A planetary-regression residual is being interpreted too strongly

In `scripts/steps/step_063_bidirectional_rotation.py:318`, the regression
outcome is `drot_sim`: rotation produced by the standard planetary
integration. The covariates summarize its energy kick, minimum sampled
encounter distance, perihelion distance, and inclination. A residual from
this reduced regression is not a residual from subtracting the full
planetary prediction from independent observations.

Recomputing both quantities on the same catalogue gives:

| CODE sample | N | Catalogue/simulation rank correlation | Median absolute catalogue−simulation rotation | Proxy-residual cap p | Direct catalogue−simulation cap p |
|---|---:|---:|---:|---:|---:|
| q < 3.1 AU | 54 | 0.999085 | 0.0006395° | 0.001877 | 0.81744 |
| All class 1 | 131 | 0.999183 | 0.0004897° | 0.012154 | 0.51271 |

The direct test uses the same 60° cap and the same greater-than alternative.
It shows that the significant proxy contrast is not a significant direct
catalogue-minus-integration contrast.

This does not independently exclude TEP: the catalogue boundary solutions
also come from orbital propagation. It establishes an interpretation limit.
The interesting remaining question is why this selected population samples
standard dynamical responses anisotropically. Answering it needs a population
and discovery-selection model. It cannot be settled by calling the residual
of the reduced regression a failure of standard dynamics.

## 5. Equivalent orbital time is not yet a derived TEP proper-time observable

In `scripts/steps/step_065_proper_time_slip.py:14`, the argument identifies a
periapsis-direction rotation with the same angular advance of the comet
along its orbit, then uses `delta_tau = delta_theta*r_b**2/h`.

The rate `h/r**2` is the true-anomaly sweep rate. Orbital phase and periapsis
orientation are different quantities; [REBOUND's orbital-element documentation](https://rebound.hanno-rein.de/particles/orbitalelements/)
distinguishes them. The audit demonstrates the distinction analytically:
a Kepler orbit advanced through 0.3 radians changes position by 17.19° while
its eccentricity-vector direction remains unchanged to numerical precision.

The conversion therefore defines an equivalent time scale. Interpreting it
as accumulated proper time requires a response model that connects the TEP
field to actual trajectories, light propagation, recorded measurements, and
the orbital fit. Multiplication by `r_b**2` also builds a quadratic dependence
on the chosen reference radius into the reported time amplitude.

The repository's canonical starting point is `core/definitions.md:9`:

`g_tilde = A(phi)^2 g + B(phi) dphi dphi`.

In the conformal limit `B=0`, along a specified timelike worldline,
`d_tau_tilde = A(phi) d_tau_g`. This relation does not itself equate a
fitted apsidal rotation to a clock offset. Matter trajectories, reference
clocks, and the observation model must be treated consistently.

The canonical screening definition is continuous and dependent on environment
and measurement response; it does not imply a strong detectable signal in
every channel. A common field can predict some weak or null channels. Those
responses must be specified independently of the observed outcomes if the
cross-channel comparison is to provide evidence.

## 6. Log-unit labels are now consistent

Steps 117–123 now use `np.log10` for residual models whose amplitudes are
reported in dex. The earlier natural-log labels were a unit error; converting
the old example `+0.3532556` gives `+0.1534170` dex, while rank-test p-values
are unchanged. The affected result files must be regenerated before numerical
amplitudes are quoted in the manuscript.

## What a consistent TEP analysis requires next

1. Centralize the verified coordinate transform, preserve p-value precision,
   replace the twist cutoff with numerical and covariance-aware handling,
   correct units, and regenerate all affected results and manuscript claims.
2. Specify an outer-Solar-System field and screening response within the
   canonical matter-metric framework. Derive each measured channel from that
   field, including sign, amplitude, and possible suppression.
3. Pass standard and TEP trajectories through the same measurement,
   discovery-selection, and orbit-fitting procedures. Demonstrate recovery
   of injected effects and a calibrated false-positive rate on the standard
   baseline. Synthetic validation should remain clearly separate from the
   real-data evidence ledger.
4. Fit the common physical parameters on a declared training sample and
   evaluate the fixed predictions on independent data. Keep the already
   failed prospective residual prediction in the ledger; a later spatial
   diagnostic answers a different question.

The spatial dipole remains an interesting anomaly under the tested nulls,
and its combined nominal significance was understated by a real bug. The
present calculations do not yet supply the common physical prediction needed
to claim strong, consistent evidence for TEP across the channels.
