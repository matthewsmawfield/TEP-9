# NASA/PDS/DSN Data Access Tools - Comprehensive Search Results

## Summary of Findings
**Date:** 2026-04-20
**Status:** Multiple Python tools identified; trk234 already installed

---

## ✅ ALREADY INSTALLED (Ready to Use)

### 1. PyTrk234 - NASA PDS Official TRK-2-34 Reader
**Status:** ✅ INSTALLED (trk234 0.1.dev9)
**Source:** https://github.com/NASA-PDS/PyTrk234
**Use Case:** Reading raw DSN TRK-2-34 tracking files

**Usage:**
```python
import trk234

# Read TRK-2-34 file
f = trk234.Reader('filename.tnf')
f.decode()

# Access SFDUs
for s in f.sfdu_list:
    # Get timestamp
    timestamp = s.timestamp()
    
    # Get Doppler/downlink frequency
    if hasattr(s.trk_chdo, 'dl_freq'):
        frequency = s.trk_chdo.dl_freq
    
    # Get signal power
    if hasattr(s.trk_chdo, 'pcn0'):
        power = s.trk_chdo.pcn0
```

**Key Features:**
- Reads binary TRK-2-34/TNF format (SFDU structure)
- Extracts: timestamps, frequencies, signal power, station IDs
- Official NASA-PDS maintained tool
- Decodes all 18 SFDU types

---

### 2. SpiceyPy - SPICE Toolkit
**Status:** ✅ INSTALLED (spiceypy 8.0.2)
**Source:** https://github.com/AndrewAnnex/SpiceyPy
**Use Case:** Spacecraft ephemeris, geometry calculations

**Already used in pipeline:** `step_001a_download_spice.py`, `step_001b_spice_to_json.py`

---

### 3. AstroQuery
**Status:** ✅ INSTALLED (astroquery 0.4.11)
**Source:** https://github.com/astropy/astroquery
**Use Case:** JPL Horizons queries, online astronomical data

**Features:**
- JPL Horizons interface (already used in pipeline)
- Solar System Dynamics data
- Various NASA/archive queries

---

## 🔧 AVAILABLE (Installable from PyPI)

### 4. PDR - Planetary Data Reader
**Status:** 🔧 NOT INSTALLED (pip installable)
**Install:** `pip install pdr` or `conda install -c conda-forge pdr`
**Source:** https://github.com/MillionConcepts/pdr
**Docs:** https://pdr.readthedocs.io/

**Use Case:** Universal PDS3/PDS4 data reader

**Usage:**
```python
import pdr

# Read any PDS file
data = pdr.read('/path/to/file.LBL')

# Access data objects
print(data.keys())  # List available data objects
image = data['IMAGE']  # Get image as numpy array
table = data['TABLE']  # Get table as pandas DataFrame

# Access metadata
metadata = data.metadata
instrument = data.metaget('INSTRUMENT_HOST_NAME')
```

**Key Features:**
- Single function `pdr.read()` for all PDS data types
- Returns: NumPy arrays (images), pandas DataFrames (tables)
- Supports PDS3 and PDS4 formats
- NASA-funded project (Grant 80NSSC21K0885)

---

### 5. PDS.peppi - PDS API Client
**Status:** 🔧 NOT INSTALLED (pip installable)
**Install:** `pip install pds.peppi`
**Source:** https://github.com/NASA-PDS/peppi
**Docs:** https://nasa-pds.github.io/peppi/
**Version:** 0.8.1

**Use Case:** Programmatic access to NASA PDS API

**Usage:**
```python
import pds.peppi as pep

# Connect to PDS
client = pep.PDSRegistryClient()

# Search for products
products = pep.Products(client)
juno_products = products.filter("juno").limit(20)

# Access data
for product in juno_products:
    print(product.id, product.title)
```

**Key Features:**
- Pythonic interface to PDS Search API
- Query products by mission, target, instrument
- Download data programmatically
- Python 3.12+ required

---

## 📊 OTHER RELEVANT TOOLS FOUND

### 6. NASA-DSN-E - DSN Monitor/Downloader
**Source:** https://github.com/ntfargo/NASA-DSN-E
**Use Case:** Real-time DSN monitoring, historical data analysis

**Features:**
- Fetches real-time data from DSN Now
- Historical analysis capabilities
- Predictive analytics with ML
- Web dashboard interface
- Python 3.10+

**Note:** This is for real-time monitoring, not archival data retrieval

---

### 7. MonteCop - JPL Trajectory Tool Interface
**Source:** https://github.com/nasa-jpl/MonteCop
**Use Case:** Transfer spacecraft trajectories between NASA tools (Monte ↔ Copernicus)

**Features:**
- JPL Monte flight mechanics tool interface
- Trajectory solution conversion
- NASA official tool

---

### 8. pds3 - PDS3 Reader (Legacy)
**Source:** https://github.com/mkelley/pds3
**Use Case:** Read PDS3 format files

**Note:** Use `pdr` instead - it's more comprehensive and supports PDS4

---

## 🎯 RECOMMENDED ACTION PLAN

### Immediate (Use What's Installed)
1. **Use trk234** to parse existing MESSENGER .dat file
   - Test if the 3.3KB file contains valid tracking data
   - Extract any Doppler measurements available

2. **If trk234 works:** 
   - Download more TRK files from NASA PDS
   - Process with `step_003_dsn_framework.py`

### Short-term (Enhance Capabilities)
3. **Install pdr:**
   ```bash
   pip install pdr
   ```
   - Provides universal PDS file reading
   - May handle additional DSN data formats

4. **Install pds.peppi:**
   ```bash
   pip install pds.peppi
   ```
   - Enable programmatic PDS data discovery
   - Automate finding Juno 2013 TRK files

### Medium-term (Full Integration)
5. **Update requirements.txt:**
   ```
   pdr>=0.0.1
   pds.peppi>=0.8.1
   ```

6. **Enhance step_003_dsn_framework.py:**
   - Add pdr as fallback parser for unknown formats
   - Add pds.peppi for automated PDS product discovery
   - Add download capability via PDS API

---

## 🔍 TESTING trk234 ON EXISTING DATA

Let's test the installed trk234 on the MESSENGER file:

```python
import trk234
from pathlib import Path

# Test file
filepath = Path('/Users/matthewsmawfield/www/TEP-EFA/data/raw/dsn_tracking/MESSENGER_2005/NASA_PDS/151161600sc236dss15ddor_234.dat')

# Read with trk234
f = trk234.Reader(str(filepath))
f.decode()

print(f"Number of SFDUs: {len(f.sfdu_list)}")

# Try to extract data
for i, s in enumerate(f.sfdu_list[:5]):
    print(f"\nSFDU {i}:")
    if hasattr(s, 'timestamp'):
        print(f"  Time: {s.timestamp()}")
    if hasattr(s, 'trk_chdo') and s.trk_chdo:
        if hasattr(s.trk_chdo, 'dl_freq'):
            print(f"  DL Freq: {s.trk_chdo.dl_freq}")
        if hasattr(s.trk_chdo, 'doppler'):
            print(f"  Doppler: {s.trk_chdo.doppler}")
```

**Next Step:** Run this test to verify trk234 can parse the existing file.

---

## 📋 INSTALLATION COMMANDS SUMMARY

```bash
# Already installed (verified)
# - trk234
# - spiceypy
# - astroquery

# Recommended additions
pip install pdr                    # Universal PDS reader
pip install pds.peppi                # PDS API client

# Optional
pip install spacepy                  # Space science utilities
```

---

## 🌐 KEY URLS

| Resource | URL |
|----------|-----|
| PyTrk234 GitHub | https://github.com/NASA-PDS/PyTrk234 |
| PDR GitHub | https://github.com/MillionConcepts/pdr |
| PDR Docs | https://pdr.readthedocs.io/ |
| peppi GitHub | https://github.com/NASA-PDS/peppi |
| peppi Docs | https://nasa-pds.github.io/peppi/ |
| NASA PDS | https://pds.nasa.gov/ |
| PDS Portal | https://pds.mcp.nasa.gov/portal/search |
| DSN Now | https://eyes.nasa.gov/dsn/dsn.html |

---

**Recommendation:** 
1. Test trk234 on the existing MESSENGER file
2. If it works, download more TRK files from PDS
3. Install `pdr` and `pds.peppi` for enhanced capabilities
4. Update pipeline to use these tools
