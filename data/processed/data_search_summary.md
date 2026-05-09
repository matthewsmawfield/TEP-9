# Raw DSN Data Search - Complete Summary

## Search Results Overview
**Date:** 2026-04-20  
**Query:** Find NASA DSN tracking data (TRK-2-34/TRK-2-25) for TEP analysis

---

## 1. EXISTING LOCAL DATA

### MESSENGER 2005 File
**Path:** `data/raw/dsn_tracking/MESSENGER_2005/NASA_PDS/151161600sc236dss15ddor_234.dat`
**Size:** 3.3 KB
**Status:** ❌ NOT VALID TRK-2-34

**Analysis Results:**
- File type: "PDS (CCSD) image data" (from `file` command)
- trk234 parser fails with KeyError: 88 (unknown format code)
- Not Doppler tracking data
- Likely a label file, metadata, or unrelated PDS product

**Conclusion:** This is not usable for TEP analysis. Need real TRK-2-34/ATDF files.

---

## 2. AVAILABLE PYTHON TOOLS

### ✅ Already Installed
| Tool | Version | Purpose | Status |
|------|---------|---------|--------|
| trk234 | 0.1.dev9 | Read TRK-2-34 files | ✅ Ready |
| spiceypy | 8.0.2 | SPICE kernels/ephemeris | ✅ Ready |
| astroquery | 0.4.11 | JPL Horizons queries | ✅ Ready |

### 🔧 Installable from PyPI
| Tool | Version | Purpose | Install Command |
|------|---------|---------|-----------------|
| pdr | latest | Universal PDS reader | `pip install pdr` |
| pds.peppi | 0.8.1 | PDS API client | `pip install pds.peppi` |
| spacepy | latest | Space science utilities | `pip install spacepy` |

### 📚 GitHub Tools (External)
| Tool | URL | Purpose |
|------|-----|---------|
| PyTrk234 | https://github.com/NASA-PDS/PyTrk234 | TRK-2-34 reader (we have this) |
| pdr | https://github.com/MillionConcepts/pdr | Universal PDS reader |
| peppi | https://github.com/NASA-PDS/peppi | PDS API access |
| NASA-DSN-E | https://github.com/ntfargo/NASA-DSN-E | DSN monitoring dashboard |

---

## 3. DATA ACCESS TESTS

### NASA PDS API
**Status:** ✅ Partially Working

**Test Results:**
- ✅ Basic API accessible: `https://pds.nasa.gov/api/search/1/products`
- ✅ Can retrieve products (79683813 total in index)
- ❌ PDS Radio Science Node (`pds-rn.jpl.nasa.gov`) - DNS resolution fails
- ❌ Search queries return 400 errors (complex query syntax required)
- ✅ PDS Portal accessible: `https://pds.mcp.nasa.gov/portal/search`

**Conclusion:** Automated bulk download is not straightforward. Manual browsing + download is more reliable.

---

## 4. OPTIONS TO ACQUIRE REAL DSN DATA

### Option A: Manual Download (Recommended - Fastest)
**Source:** NASA PDS Portal
**URL:** https://pds.mcp.nasa.gov/portal/search

**Steps:**
1. Visit PDS Portal
2. Search: "JNO-E-RSS-1-EDR" (Juno Earth flyby)
3. Filter: Date range 2013-10-08 to 2013-10-10
4. Download: TRK-2-34 or TRK-2-25 (ATDF) files
5. Place in: `data/raw/dsn_tracking/Juno_2013/`
6. Re-run: `python scripts/steps/step_003_dsn_framework.py`

**Expected Data Volume:** 50-100 MB per 48-hour arc
**Timeline:** Immediate (once downloaded)

---

### Option B: Contact NASA PDS Directly
**Email:** pds-rn@jpl.nasa.gov  
**Subject:** "Juno 2013 Earth Flyby TRK-2-34 Data Request for Academic Research"

**Request Template:**
```
Data Products:
- Mission: Juno Earth flyby 2013
- Date: 2013-10-08 to 2013-10-10
- Format: TRK-2-34 or TRK-2-25 (ATDF)
- Stations: DSS-24, DSS-25, DSS-54

Scientific Justification:
Academic research on Earth flyby velocity anomalies.
Requesting raw tracking data for independent minimal OD
analysis to test for TEP signals.

Contact: [Your institution/email]
```

**Timeline:** 2-4 weeks

---

### Option C: Enhance Framework with pds.peppi
**Install:** `pip install pds.peppi`

**Implementation:**
```python
import pds.peppi as pep

client = pep.PDSRegistryClient()
products = pep.Products(client)

# Search for Juno products
juno_products = products.filter("jno-e-rss-1-edr").limit(50)

# Download found products
for product in juno_products:
    # Extract download URL
    # Download TRK files
```

**Benefit:** Programmatic discovery of data products
**Timeline:** Requires development time

---

## 5. RECOMMENDATION

### Immediate Action: Manual Download
1. **Go to:** https://pds.mcp.nasa.gov/portal/search
2. **Search for:** "JNO-E-RSS-1-EDR" or "juno earth flyby 2013"
3. **Download:** TRK-2-34 files for Oct 8-10, 2013
4. **Place in:** `data/raw/dsn_tracking/Juno_2013/`

### Verification Steps:
```bash
# 1. Check downloaded files
ls -lh data/raw/dsn_tracking/Juno_2013/

# 2. Verify with trk234
python3 -c "
import trk234
f = trk234.Reader('data/raw/dsn_tracking/Juno_2013/*.trk')
f.decode()
print(f'SFDUs: {len(f.sfdu_list)}')
"

# 3. Run DSN reanalysis
python scripts/steps/step_003_dsn_framework.py
```

---

## 6. TECHNICAL NOTES

### TRK-2-34 Format
- Binary format with SFDU (Standard Formatted Data Unit) structure
- 18 distinct SFDU types in 5 groups
- Produced by DSN in near-real-time during tracking
- Contains: timestamps, Doppler, range, signal power, station IDs

### Minimal OD Requirements
Our framework (`step_003_dsn_framework.py`) requires:
- Raw TRK-2-34 or ATDF files (not processed ODF)
- Files covering perigee ± 48 hours
- Valid Doppler measurements (not just metadata)

### File Size Expectations
- Real TRK-2-34: 10-100 MB per day of tracking
- Our 3.3KB file: Too small - must be metadata only
- Juno 2013 48-hour arc: Expected ~50-100 MB

---

## 7. NEXT STEPS

1. **Download real DSN data** from NASA PDS Portal
2. **Verify** with trk234 parser
3. **Run** `step_003_dsn_framework.py` for minimal OD analysis
4. **Get** falsification test results for TEP model

**Alternative:** Install `pdr` and `pds.peppi` to enhance automated discovery capabilities:
```bash
pip install pdr pds.peppi
```

---

**Summary:** The tools are ready (trk234 installed, others available). We just need the actual raw DSN tracking data files from NASA PDS to proceed with the minimal OD analysis and falsification testing.
