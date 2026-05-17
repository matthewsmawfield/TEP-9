# The Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19454863.svg)](https://doi.org/10.5281/zenodo.19454863)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![TEP-EFA: Earth Flyby Anomaly](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Yogyakarta)  
**First published:** 17 May 2026 · **Last updated:** 17 May 2026  
**Status:** Preprint (Open for Collaboration)  
**DOI:** [10.5281/zenodo.19454863](https://doi.org/10.5281/zenodo.19454863)  
**Website:** [https://mlsmawfield.com/tep/efa/](https://mlsmawfield.com/tep/efa/)  
**Paper Series:** TEP Series: Paper 15 (Earth Flyby)

## Abstract

Twelve Earth gravity assist flybys spanning nine spacecraft are analyzed within the Temporal Equivalence Principle (TEP) framework. TEP posits that global simultaneity is inherently non-integrable, with the rate of time represented as a dynamical scalar field φ. All non-gravitational matter couples universally to a causal matter metric through conformal coupling A(φ) = exp(β φ/M_Pl), producing a scalar force F = β_eff c² ∇φ/M_Pl on test masses, where β_eff = β × S_⊕(r) incorporates geometric screening via Temporal Topology. The screening factor S_⊕(r) encodes continuous suppression of Temporal Shear in density gradients, with a characteristic transition radius R_sol ≈ 4146 km derived from the UCD saturation model and independently validated by GNSS atomic clock correlations (λ_TEP ≈ 4000 km).

The scalar force manifests as a "Phantom Mass" artifact — velocity anomalies that mimic unmodeled gravitational mass distributions. The non-radial component, modulated by Earth's oblateness (J2, J3, J4), trajectory asymmetry, velocity-dependent disformal coupling, and perigee plasma environment, produces the observed flyby anomaly. Four published anomalies are retained in the catalog (NEAR, Galileo 1990, Rosetta 2005, Cassini), alongside five published nulls/bounds and three flybys without public anomaly reports.

Information-criterion comparison on the full catalog of nine flybys with published observations (n = 9) strongly favors the TEP restricted model over both the Null (ΔBIC ≈ 714) and the Anderson empirical baseline (ΔBIC ≈ 79), using a single fitted parameter β with λ_TEP, S_⊕, and v_trans pre-specified from independent data and first-principles derivations. The gated n = 3 subset yields weaker TEP-restricted vs Anderson separation (ΔBIC ≈ 19), confirming that the null flybys are essential for discriminating the physics-based model from the empirical baseline. The random-effects summary β_RE ≈ 2.56 × 10⁻³ ± 7.85 × 10⁻⁴ (SE; between-flyby τ ≈ 1.49 × 10⁻³) quantifies the honest cross-flyby amplitude scale, while the per-flyby fits span 1.01 × 10⁻³ to 5.33 × 10⁻³ consistent with geometry-dependent modulation. All fitted amplitudes satisfy PPN constraints via Temporal Topology screening (|γ − 1| ≈ 2β_eff² < 2.3 × 10⁻⁵). Cassini's small predicted anomaly in the mixed disformal regime, where conformal-gradient and disformal terms partially cancel at the reference coupling, is consistent with the velocity-dependent regime structure predicted by the TEP field equations (v_trans ≈ 16.8 km/s), though the literature geometry sign remains an open diagnostic stress test rather than a resolved confirmation.

This work demonstrates that the Temporal Equivalence Principle framework quantitatively explains the Earth flyby anomaly while remaining consistent with precision solar system constraints, establishing a new avenue for exploring the intersection of gravity, time, and matter.

## Key Findings

Three gated Step 008 fits (NEAR, Galileo 1990, Rosetta 2005) with β spanning 1.01×10⁻³ to 5.33×10⁻³ (factor 5.3), inverse-variance weighted mean β = 1.73×10⁻³ ± 6.82×10⁻⁵, and random-effects β_RE = 2.56×10⁻³ ± 7.85×10⁻⁴; Cassini stays in the literature catalog but is excluded from the sign-gated ensemble on sign mismatch at β_ref. All gated screened couplings satisfy |γ − 1| ≈ 2β_eff² < 2.3×10⁻⁵. Step 026 model comparison on the full n = 9 catalog strongly favors TEP restricted over Null (ΔBIC ≈ 714) and over Anderson empirical (ΔBIC ≈ 78). Step 039 extends the check to the null-result catalog: raw true positives on the three anomalies, raw true nulls on several published bounds, one deterministic fixed-amplitude warning case (Juno), and zero uncertainty-aware raw-tension cases after random-effects prediction scatter. Step 008 bootstrap (n = 10,000) and leave-one-out stability coefficient ≈ 0.148 on the gated trio.

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
