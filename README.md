# The Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19454863.svg)](https://doi.org/10.5281/zenodo.19454863)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![TEP-EFA: Earth Flyby Anomaly](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Yogyakarta)  
**First published:** 10 May 2026 · **Last updated:** 10 May 2026  
**Status:** Complete  
**DOI:** [10.5281/zenodo.19454863](https://doi.org/10.5281/zenodo.19454863)  
**Website:** [https://mlsmawfield.com/tep/efa/](https://mlsmawfield.com/tep/efa/)  
**Paper Series:** TEP Series: Paper 15 (Earth Flyby)

## Abstract

Twelve Earth gravity assist flybys spanning nine spacecraft (NEAR, Galileo 1990/1992, Cassini, Rosetta 2005/2007/2009, MESSENGER, Juno, Stardust, OSIRIS-REx, BepiColombo) are analyzed within the Temporal Equivalence Principle (TEP) framework. The TEP framework proposes that global simultaneity is inherently non-integrable, with the rate of time represented as a dynamical scalar field φ. All non-gravitational matter couples universally to a causal matter metric through conformal coupling A(φ) = exp(β φ/M_Pl), where β is a dimensionless coupling constant and M_Pl is the reduced Planck mass. This coupling produces a scalar force F = β_eff c² ∇φ/M_Pl on test masses, where β_eff = β × S_⊕(r) incorporates geometric screening via Temporal Topology. The screening factor S_⊕(r) encodes continuous suppression of Temporal Shear in density gradients, with a characteristic transition radius R_sol ≈ 4146 km derived from first-principles Universal Critical Density (UCD) saturation model and independently validated by GNSS atomic clock correlations (λ_TEP ≈ 4000 km).

The scalar force manifests as a "Phantom Mass" artifact—velocity anomalies that mimic unmodeled gravitational mass distributions through field-gradient couplings. The radial component of this force is indistinguishable from a small shift in GM and is absorbed by orbit determination programs. The non-radial component, modulated by Earth's oblateness (J2, J3, J4), trajectory asymmetry, and velocity-dependent disformal coupling, produces the observed flyby anomaly. Four primary detections are successfully fitted: NEAR (13.46 ± 0.01 mm/s), Galileo 1990 (3.92 ± 0.03 mm/s), Rosetta 2005 (1.82 ± 0.05 mm/s), and Cassini (0.11 ± 0.05 mm/s). The Cassini sign reversal is resolved through velocity-dependent disformal coupling that reverses the prediction sign for high-velocity (v > 16 km/s) anti-aligned trajectories.

Field values are computed self-consistently from the first-principles field equation φ(ρ) = Λ [n Λ^(n+4) M_Pl / (β ρ)]^(1/(n+1)), yielding φ_earth ≈ 2.4×10⁴ GeV and φ_space ≈ 2.0×10¹⁰ GeV for n = 3 and Λ = 10 MeV. Fitted β values span a factor of 24.0 (2.40×10⁻⁵ to 5.76×10⁻⁴), with the inverse-variance weighted mean β = 4.64×10⁻⁴ ± 2.32×10⁻⁵ (5% uncertainty). The heterogeneity in fitted values (I² ≈ 100%) reflects genuine geometry-dependent physical variation in effective coupling across flybys, consistent with the continuous Temporal Topology gradient structure. All fitted values satisfy solar system PPN constraints (|γ - 1| = 2β_eff² < 2.3×10⁻⁵) with safety margins exceeding 10³×. Bayesian model comparison strongly favors TEP over the null model (Bayes factor B₁₀ = 45.3, Akaike weight = 92.5%). The model achieves R² = 0.89 between predicted and observed anomalies, with residuals consistent with normal distribution (Shapiro-Wilk p = 0.45). Null results for high-altitude flybys (Galileo 1992, MESSENGER, Juno) are explained through trajectory geometry and orbit determination filtering, with numerical simulation demonstrating >50% TEP signal suppression by modern empirical acceleration methods.

This work bridges the gap between precision solar system tests and cosmological dynamics, demonstrating that the Temporal Equivalence Principle is a measurable physical effect with direct experimental validation and providing a new avenue for exploring the intersection of gravity, time, and matter, and ultimately shedding new light on the fundamental nature of spacetime.

## Key Findings

Four successful TEP fits to flyby anomalies (NEAR, Galileo 1990, Rosetta 2005, Cassini) with fitted β values spanning a factor of 24.0 (2.40×10⁻⁵ to 5.76×10⁻⁴). All fitted parameters satisfy solar system PPN constraints (|γ - 1| = 2β_eff²) with safety margins of 2-3 orders of magnitude. The model predicts null results for symmetric trajectories (Galileo 1992, MESSENGER) and explains most null results through trajectory geometry. Bootstrap (n=10,000) and leave-one-out cross-validation indicate stability.

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
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064366](https://doi.org/10.5281/zenodo.18064366) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059251](https://doi.org/10.5281/zenodo.18059251) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064582](https://doi.org/10.5281/zenodo.18064582) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109761](https://doi.org/10.5281/zenodo.18109761) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | The Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | The Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102062](https://doi.org/10.5281/zenodo.19102062) |
| **Paper 15** | **TEP-EFA** (This repo) | Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454863](https://doi.org/10.5281/zenodo.19454863) |

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
  doi={10.5281/zenodo.19454863},
  note={Preprint v0.1 (Yogyakarta)}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
