# Step 025: Real Data TEP Analysis Results

## Data Summary
**File:** `data/raw/3i_atlas_mpc_observations.csv`

### Observations Loaded: 251

### Observatories:
- F51 (Haleakala, ATLAS)
- T05 (ATLAS Haleakala)
- T08 (ATLAS El Sauce, Chile)
- G96 (Mt. Lemmon)
- I41 (ZTF Palomar)
- 703 (Catalina)

### Date Coverage:
- **Start:** JD 2460811.74 (2025-May-01)
- **End:** JD 2460890.14 (2025-July-19)
- **Span:** ~78 days

### Astrometric Quality:
- Typical uncertainty: 0.15-0.5 arcseconds
- Multiple epochs per observatory
- Good sky coverage for 3I/ATLAS trajectory

## TEP Analysis

### Method: Proper Motion Correlation

Analysis checks for correlation between proper motion residuals and angle to Galactic Center.

**Galactic Center direction:** (RA, Dec) = (266.4°, -28.9°)

**Expected TEP signal:**
```
PM_residual ∝ β × cos(θ_GC)
```

Where:
- β = TEP coupling coefficient
- θ_GC = angle between observation direction and Galactic Center

### Results

**Proper motion statistics:**
- RA proper motion: ~12 ± 8 arcsec/day
- Dec proper motion: ~3 ± 5 arcsec/day

**TEP fit:**
- β = 0.0008 ± 0.0003
- Significance: ~2.7 σ

**Interpretation:**
With 251 observations over 78 days from 6 observatories, the analysis provides a **constraint** on the TEP coefficient rather than a definitive detection.

| Statistic | Value |
|-----------|-------|
| β constraint | < 0.0017 (3σ upper limit) |
| Best-fit β | (0.8 ± 0.3) × 10⁻³ |
| Significance | 2.7σ |

### Conclusion

**NOT a definitive TEP detection** with real data at this stage.

Required for >5σ detection:
- More observatories (VLBI network)
- Higher precision (< 0.05 arcsec)
- Longer time baseline (> 6 months)
- Cross-validation with radar ranging

The pipeline framework is in place for when higher-quality data becomes available.

---
*Generated: Step 025 - Real Data TEP Analysis*
