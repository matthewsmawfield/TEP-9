# The Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19454862.svg)](https://doi.org/10.5281/zenodo.19454862)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![TEP-EFA: Earth Flyby Anomaly](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Yogyakarta)  
**First published:** 17 May 2026 · **Last updated:** 14 September 2026
**Status:** Preprint (Open for Collaboration)  
**DOI:** [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862)  
**Website:** [https://mlsmawfield.com/tep/efa/](https://mlsmawfield.com/tep/efa/)  
**Paper Series:** TEP Series: Paper 15 (Earth Flyby)

## Abstract


Twelve Earth gravity assist flybys spanning nine spacecraft are analyzed within the Temporal Equivalence Principle (TEP) framework. TEP posits that proper time is governed by a dynamical scalar field φ, whose conformal gradients and disformal velocity-dependent terms can produce screened trajectory-level residuals in Earth-vicinity motion. All non-gravitational matter couples universally to a causal matter metric through conformal coupling A(φ) = exp(β_A φ/M_Pl), producing a scalar force F = β_A,eff c² ∇φ/M_Pl on test masses, where β_A,eff = β_A × S_⊕(r) incorporates geometric screening via Temporal Topology. The screening factor S_⊕(r) encodes continuous suppression of Temporal Shear in density gradients, with a characteristic transition radius R_sol = R_T(M_⊕) = (3M_⊕/4πρ_T)^(1/3) ≈ 4146 km for ρ_T ≈ 20 g/cm³. Since the UCD calibration of ρ_T descends from the GNSS coherence-length measurement, this terrestrial length scale is a single shared input rather than an independent validation.

The scalar force manifests as a "Phantom Mass" artifact—velocity anomalies that mimic unmodeled gravitational mass distributions. The non-radial component, modulated by Earth's oblateness (J2, J3, J4), trajectory asymmetry, velocity-dependent disformal coupling, and perigee plasma environment, produces the observed flyby anomaly. Four published anomalies are retained in the catalog (NEAR, Galileo 1990, Rosetta 2005, Cassini), alongside five published nulls/bounds and three flybys without public anomaly reports. Under a full-catalog likelihood with an explicit geometry-spread systematic uncertainty term, information-criterion comparisons favour the restricted TEP model over the Null and improve on the Anderson empirical baseline. Because the sample size is small and the adopted systematic treatment is load-bearing, these comparisons are read as model-selection support rather than calibrated evidence.

A structural consistency check is provided by the gated n = 3 detection subset, in which the anomaly magnitudes are monotonically ordered with the trajectory-asymmetry factor (cosδ_in − cosδ_out). With n = 3 this is an ordering observation rather than a correlation measurement (six possible orderings give p = 1/6, a ~1σ ceiling), so it is reported as a qualitative consistency check rather than a primary discriminator. The inverse-variance weighted mean β_fit ≈ 1.73 × 10⁻³ (random-effects SE 7.85 × 10⁻⁴; between-flyby τ ≈ 1.49 × 10⁻³) quantifies the cross-flyby effective response scale; it is not the bare conformal coupling β_A. The per-flyby fits span 1.01 × 10⁻³ to 5.33 × 10⁻³ consistent with geometry-dependent modulation. The fitted amplitudes are trajectory-response coefficients, not the Cassini source charge; Cassini PPN compliance is evaluated through the solar source-charge projection α_eff = S_Σ α_0 (Section 4.6.1a). Cassini's predicted sign is negative where the published anomaly is positive; this sign mismatch is propagated from the Anderson et al. (2008) geometry factor and is retained as an open diagnostic stress test rather than a resolved confirmation. The geometry envelope carries seven nominal deterministic coefficients plus one fitted amplitude against four retained anomalies, so the model-selection comparisons are read as organisational support under the adopted convention rather than calibrated evidence.

This work shows that a restricted TEP Temporal-Shear model quantitatively organizes the published Earth flyby anomaly catalogue under the adopted Anderson trajectory-geometry convention, while remaining consistent with precision solar system constraints and identifying specific stress tests—principally independent reconstruction of the NEAR asymptotic-state geometry and raw DSN reanalysis—that are required to advance from model organisation to decisive confirmation.

The terrestrial screening factor S_⊕(r) constitutes the Earth-vicinity geometric realization of the abstract environmental operator S_Σ(E). Rather than a hard boundary, this radial profile traces the continuous gradient-suppression of the Temporal Shear within Earth's local potential well.
## Key Findings

Three gated Step 008 fits (NEAR, Galileo 1990, Rosetta 2005) with β<sub>fit</sub> spanning 1.01×10⁻³ to 5.33×10⁻³ (factor 5.3), inverse-variance weighted mean β<sub>fit</sub> = 1.73×10⁻³ ± 6.82×10⁻⁵, and random-effects β<sub>RE</sub> = 2.56×10⁻³ ± 7.85×10⁻⁴; Cassini stays in the literature catalog but is excluded from the sign-gated ensemble on sign mismatch at β<sub>ref</sub>. All gated screened couplings satisfy |γ − 1| ≈ 2β<sub>A,eff</sub>² < 2.3×10⁻⁵. Step 026 model comparison on the full n = 9 catalog strongly favors TEP restricted over Null (ΔBIC ≈ 714) and over Anderson empirical (ΔBIC ≈ 78). Step 039 extends the check to the null-result catalog: raw true positives on the three anomalies, raw true nulls on several published bounds, one deterministic fixed-amplitude warning case (Juno), and zero uncertainty-aware raw-tension cases after random-effects prediction scatter. Step 008 bootstrap (n = 10,000) and leave-one-out stability coefficient ≈ 0.148 on the gated trio.

---

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 1** | [TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS) | Global Time Echoes: Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) |
| **Paper 2** | [TEP-GNSS-II](https://github.com/matthewsmawfield/TEP-GNSS-II) | Global Time Echoes: 25-Year Temporal Evolution | [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) |
| **Paper 3** | [TEP-GNSS-RINEX](https://github.com/matthewsmawfield/TEP-GNSS-RINEX) | Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) |
| **Paper 4** | [TEP-GL](https://github.com/matthewsmawfield/TEP-GL) | Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations | [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) |
| **Paper 5** | [TEP-GTE](https://github.com/matthewsmawfield/TEP-GTE) | Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle | [10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) |
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064365](https://doi.org/10.5281/zenodo.18064365) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059250](https://doi.org/10.5281/zenodo.18059250) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064581](https://doi.org/10.5281/zenodo.18064581) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109761](https://doi.org/10.5281/zenodo.18109761) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | The Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | The Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102061](https://doi.org/10.5281/zenodo.19102061) |
| **Paper 15** | **TEP-EFA** (This repo) | Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862) |
| **Paper 16** | [TEP-J0437](https://github.com/matthewsmawfield/TEP-J0437) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029) |

## Directory Structure

```
TEP-EFA/
├── data/                          # Flyby trajectories and processed outputs
│   └── raw/                       # Original data sources
├── scripts/
│   ├── steps/                     # Sequential analysis pipeline
│   └── utils/                     # Shared utilities
├── results/                       # Analytical outputs and figures
├── site/
│   ├── components/                # Manuscript HTML sections
│   └── public/                    # Static assets
├── logs/                          # Pipeline execution logs
├── config/                        # Pipeline configuration
├── README.md                      # This file
└── requirements.txt               # Python dependencies
```

## Installation

```bash
# Clone repository
git clone https://github.com/matthewsmawfield/TEP-EFA.git
cd TEP-EFA

# Install dependencies
pip install -r requirements.txt
```

## Reproduction Steps

The flyby TEP analysis pipeline is fully automated. Run the complete pipeline:

```bash
python scripts/run_all.py
```

This will populate `results/` and `data/processed/` with fresh analysis outputs.

## Citation

```bibtex
@article{smawfield2026earthflyby,
  title={The Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.19454862},
  note={Preprint v0.2 (Yogyakarta)}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
