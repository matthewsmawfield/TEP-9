"""
Step 030: Juno 2013 Raw DSN Reanalysis - Critical Validation Test

This step implements a raw-archive check by:
1. Downloading raw TRK-2-25 (ATDF) data from NASA PDS for Juno 2013 Earth flyby
2. Building a Doppler **pairwise-difference residual proxy** (not orbit determination)
3. Summarizing values in a ±2 h perigee window
4. Comparing to a TEP-predicted +2.25 mm/s scale and a ±0.08 mm/s decision threshold

On successful TRK processing, exports ``data/time_resolved_flyby_residuals/Juno.json``
for Step 042 when the pairwise proxy is Doppler-derived **mm/s**; archives full rows to
``results/step030_juno_pairwise_residual_series.json``. When only ``ramp_freq_hz`` is
available (typical NJPL TNF), archives Hz pairwise rows to
``results/step030_juno_pairwise_ramp_freq_series.json`` (not fed to Step 042).

**archive_presmoothing_sensitivity** (Step 030 JSON): real-archive only. Applies
uniform W-sample presmoothing to ``ramp_freq_hz`` per station (mimicking one leg of
operational count / Doppler averaging), then recomputes sequential pairwise deltas.
This answers an Occam-adjacent observable question—whether that stage materially moves
the pairwise proxy on the ingested product. It does **not** implement empirical
acceleration states or batch OD; MONTE / ODP-class adjudication remains out of scope.
Re-runs Step 042 when
``results/step038_trajectory_series.json`` exists so residual correlations pick up the
new sidecar.

A definitive minimal-OD pipeline (e.g. MONTE / JPL ODP-class fit with EGM-96,
no empirical accelerations) is **not** implemented here; ``intended_minimal_od_when_bound``
documents the target configuration for when that software is attached.

**Public ephemeris batch (automatic):** when ``data/raw/jpl_horizons/Juno_2013/Juno_2013_trajectory.json``
and ``results/step038_3d_state_vectors.json`` (Juno block) exist, Step 030 runs
:func:`scripts.utils.minimal_od_juno_horizons.run_horizons_public_ephemeris_batch`
and writes ``horizons_public_ephemeris_batch`` plus legacy ``horizons_minimal_od_batch``
(velocity-only summary) into ``results/step030_juno_dsn_reanalysis.json``. Set
``TEP_030_SKIP_HORIZONS_PUBLIC_OD=1`` to disable. This uses only public Horizons
``range_m`` / ``velocity_m_s`` — not TRK Doppler or MONTE-class OD.

Set ``TEP_030_PROBE_PUBLIC_URLS=1`` to append ``juno_2013_public_url_probe`` to the
Step 030 JSON (lightweight streamed GET per curated public URL). Requires network;
default is off so CI remains deterministic.

FALSIFICATION GATE (proxy-only):
--------------------------------
The perigee-window statistic is a **pairwise Doppler difference** proxy
(Δv ≈ c Δf/f between consecutive samples), **not** a post-fit OD residual.
Status values ``JUNO_PAIRWISE_PROXY_BELOW_0p08_MM_S_GATE`` /
``JUNO_PAIRWISE_PROXY_ABOVE_0p08_MM_S_GATE`` compare |mean proxy| only to
the fixed 0.08 mm/s operational gate used historically in this repository.
They are **not** commensurate with literature minimal-OD Δv at the +2.25 mm/s
TEP reference scale; OD-class adjudication requires MONTE / ODP-class residuals.

Data Source: NASA PDS Radio Science Node
- Primary: https://pds-rn.jpl.nasa.gov/data/jno-e-rss-1-edr/ (Earth flyby)
- Backup: https://pds-rn.jpl.nasa.gov/data/jno-j-rss-1-edr/ (Jupiter cruise)
- Verified HTTPS browse (outer cruise gravity; see ``data/raw/dsn_tracking/Juno_2013/DOWNLOAD_INSTRUCTIONS.txt`` for temporal caveats): https://pds-atmospheres.nmsu.edu/PDS/data/jnogrv_0001/DATA/TNF/

Author: TEP-EFA Pipeline
Date: 2026-04-19
"""

import os
import sys
import json
import subprocess
import numpy as np
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import requests
from scipy.ndimage import uniform_filter1d
from datetime import datetime, timedelta, timezone
import time
import re

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.dsn_pds_ingest import ingest_mission_tracking
from scripts.utils.dsn_tracking_discovery import (
    _label_path_for_data_file,
    discover_dsn_tracking_file,
    is_trk234_archive,
)
from scripts.utils.flyby_time_series_residuals import (
    build_juno_pairwise_042_sidecar_from_rows,
    save_juno_pairwise_step030_archive,
    save_juno_ramp_pairwise_step030_archive,
    write_juno_042_sidecar,
)
from scripts.utils.step_logger import StepLogger
from scripts.utils.trk234_extract import extract_trk234_measurements


from scripts.utils.pds3_lbl_time import parse_pds3_lbl_start_stop_utc  # noqa: E402 — re-export for tests


class JunoDSNReanalysis:
    """
    Complete raw DSN reanalysis pipeline for Juno 2013 Earth flyby.
    """
    
    # NASA PDS endpoints for Juno Radio Science
    PDS_RN_BASE = "https://pds-rn.jpl.nasa.gov/data/"
    JUNO_EARTH_FLYBY_COLLECTION = "jno-e-rss-1-edr/"
    # Verified HTTPS bulk TNF browse (outer-cruise OCRU; see DOWNLOAD_INSTRUCTIONS.txt).
    ATMOSPHERES_JNO_OCRU_TNF_INDEX = (
        "https://pds-atmospheres.nmsu.edu/PDS/data/jnogrv_0001/DATA/TNF/"
    )

    # Curated public HTTPS targets for ``TEP_030_PROBE_PUBLIC_URLS=1`` (id, url, purpose).
    JUNO_PUBLIC_URL_PROBE_TARGETS: Tuple[Tuple[str, str, str], ...] = (
        (
            "pds_rn_jno_e_rss_edr",
            "https://pds-rn.jpl.nasa.gov/data/jno-e-rss-1-edr/",
            "Juno Earth-phase RSS EDR browse (TRK target)",
        ),
        (
            "pds_rn_jno_j_rss_edr",
            "https://pds-rn.jpl.nasa.gov/data/jno-j-rss-1-edr/",
            "Juno cruise RSS EDR browse (alternate)",
        ),
        (
            "pds_atmospheres_jnogrv_tnf",
            "https://pds-atmospheres.nmsu.edu/PDS/data/jnogrv_0001/DATA/TNF/",
            "OCRU gravity TNF index (NJPL-class; post-flyby era)",
        ),
        (
            "pds_nasa_dataset_ocru",
            "https://pds.nasa.gov/ds-view/pds/viewDataset.jsp?dsid=JUNO-J-RSS-1-OCRU-V1.0",
            "PDS dataset landing page for OCRU",
        ),
        (
            "jpl_horizons_portal",
            "https://ssd.jpl.nasa.gov/horizons/",
            "JPL Horizons (ephemeris source for Step 004 JSON)",
        ),
        (
            "naif_juno_spk_browse",
            "https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/",
            "NAIF Juno SPICE SPK directory (navigation kernels)",
        ),
    )
    
    # Juno 2013 Earth flyby parameters
    FLYBY_DATE = datetime(2013, 10, 9, 19, 21, 0)  # Perigee UTC
    ANALYSIS_WINDOW_HOURS = 48  # ±24 hours around perigee
    
    # TEP predictions
    TEP_PREDICTED_SIGNAL = 2.25  # mm/s at global β̄
    FALSIFICATION_THRESHOLD = 0.08  # mm/s at 95% confidence
    
    def __init__(self):
        self.logger = StepLogger("step_030_juno_reanalysis", PROJECT_ROOT)
        self.data_dir = PROJECT_ROOT / 'data' / 'raw' / 'dsn_tracking' / 'Juno_2013'
        self.results_dir = PROJECT_ROOT / 'results'
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Track downloaded/processed data
        self.raw_files = []
        self.doppler_data = []
        self.analysis_results = {}

    def _project_root(self) -> Path:
        return self.data_dir.parent.parent.parent.parent

    def _horizons_public_ephemeris_optional(self) -> Optional[Dict[str, Any]]:
        """
        Best-effort fit to public JPL Horizons JSON + Step 038 Juno anchor.
        Skips when ``TEP_030_SKIP_HORIZONS_PUBLIC_OD=1`` or inputs missing.
        """
        if os.environ.get("TEP_030_SKIP_HORIZONS_PUBLIC_OD") == "1":
            return None
        root = self._project_root()
        traj = root / "data" / "raw" / "jpl_horizons" / "Juno_2013" / "Juno_2013_trajectory.json"
        s40 = root / "results" / "step038_3d_state_vectors.json"
        if not traj.is_file() or not s40.is_file():
            return None
        from scripts.utils.minimal_od_juno_horizons import run_horizons_public_ephemeris_batch

        return run_horizons_public_ephemeris_batch(root)

    def _trk_ingest_temporal_audit(self, file_paths: List[str]) -> Dict[str, Any]:
        """
        For each ingested NJPL archive with a sibling PDS3 ``.LBL``, test whether
        ``START_TIME``/``STOP_TIME`` overlaps the Step 030 flyby analysis window
        (± ``ANALYSIS_WINDOW_HOURS``/2 around ``FLYBY_DATE`` UTC).

        Uses :func:`parse_pds3_lbl_start_stop_utc` so day-of-year PDS3 times parse
        correctly (``dsn_pds_ingest._parse_pds_timestamp`` is ISO-first and can
        mis-handle ``YYYY-DOY`` tokens).
        """
        perigee = self.FLYBY_DATE.replace(tzinfo=timezone.utc)
        window_h = float(self.ANALYSIS_WINDOW_HOURS)
        win_lo = perigee - timedelta(hours=window_h / 2.0)
        win_hi = perigee + timedelta(hours=window_h / 2.0)
        rows: List[Dict[str, Any]] = []
        root = PROJECT_ROOT.resolve()
        for fp in file_paths:
            p = Path(fp).resolve()
            if not p.is_file() or not is_trk234_archive(p):
                continue
            try:
                rel = str(p.relative_to(root))
            except ValueError:
                rel = str(p)
            row: Dict[str, Any] = {"file": rel}
            lbl = _label_path_for_data_file(p)
            if lbl is None or not lbl.is_file():
                row["sibling_lbl"] = None
                row["overlaps_flyby_window_utc"] = None
                row["note"] = "No sibling .LBL next to archive; overlap not evaluated"
                rows.append(row)
                continue
            try:
                rel_lbl = str(lbl.resolve().relative_to(root))
            except ValueError:
                rel_lbl = str(lbl)
            row["sibling_lbl"] = rel_lbl
            try:
                lbl_text = lbl.read_text(encoding="utf-8", errors="ignore")
                span = parse_pds3_lbl_start_stop_utc(lbl_text)
                if span is None:
                    row["overlaps_flyby_window_utc"] = None
                    row["lbl_parse_error"] = "parse_pds3_lbl_start_stop_utc returned None"
                else:
                    t0, t1 = span
                    if t1 < t0:
                        t0, t1 = t1, t0
                    row["lbl_start_utc"] = t0.isoformat()
                    row["lbl_stop_utc"] = t1.isoformat()
                    row["overlaps_flyby_window_utc"] = bool(t0 <= win_hi and t1 >= win_lo)
            except (ValueError, OSError, RuntimeError) as exc:
                row["overlaps_flyby_window_utc"] = None
                row["lbl_parse_error"] = str(exc)
            rows.append(row)

        def _known_false(r: Dict[str, Any]) -> bool:
            return r.get("overlaps_flyby_window_utc") is False

        def _known_true(r: Dict[str, Any]) -> bool:
            return r.get("overlaps_flyby_window_utc") is True

        any_overlap = any(_known_true(r) for r in rows)
        all_evaluated_false = bool(rows) and all(_known_false(r) for r in rows)
        return {
            "perigee_utc": perigee.isoformat(),
            "window_total_hours": window_h,
            "window_utc": f"{win_lo.isoformat()} to {win_hi.isoformat()}",
            "files": rows,
            "any_file_overlaps_flyby_window": any_overlap if rows else None,
            "all_evaluated_files_miss_flyby_window": all_evaluated_false,
            "interpretation": (
                "Sibling PDS3 label interval vs Earth-flyby perigee window. "
                "OCRU outer-cruise TNFs typically return false; do not treat their "
                "pairwise ramp proxy as perigee Doppler."
            ),
        }

    def _build_evidence_tier_assessment(
        self,
        *,
        trk_temporal: Dict[str, Any],
        falsification_result: Dict[str, Any],
        proxy_results: Dict[str, Any],
        horizons_pub: Optional[Dict[str, Any]],
        perigee_matched: bool,
    ) -> Dict[str, Any]:
        """Machine-readable Tier I–III closure for Juno DSN / OD-filtering narrative."""
        proxy_kind = proxy_results.get("proxy_kind")
        trk_off_epoch = bool(trk_temporal.get("all_evaluated_files_miss_flyby_window"))
        trk_commensurate = bool(
            perigee_matched
            and proxy_kind == "doppler_pair_mm_s"
            and falsification_result.get("status")
            not in (None, "inconclusive", "skipped")
        )

        tier_iii_status = "incomplete"
        tier_iii_detail = (
            "No perigee-overlapping Doppler-populated TRK arc in the repository; "
            "pairwise proxy is not commensurate with published post-OD Δv."
        )
        if trk_commensurate:
            tier_iii_status = "partial"
            tier_iii_detail = (
                "Perigee-window Doppler pairwise proxy available; "
                "MONTE/ODP-class batch OD with TEP forces still required for navigation-scale closure."
            )
        elif trk_off_epoch:
            tier_iii_detail = (
                "Ingested PDS3 label intervals miss the 2013 Earth-flyby window "
                "(e.g. OCRU outer-cruise 2015 TNF); ramp-frequency Hz proxy only."
            )

        horizons_note: Optional[str] = None
        if horizons_pub is not None:
            vo = horizons_pub.get("velocity_only") or {}
            bls = vo.get("batch_least_squares") or {}
            rms = bls.get("post_fit_rms_rdot_residual_m_s")
            if rms is not None:
                horizons_note = (
                    f"Public Horizons geocentric deldot batch fit (pure Kepler, ±3 h): "
                    f"post-fit RMS ≈ {float(rms):.2f} m/s — ephemeris consistency check only, "
                    "not TRK Doppler or empirical-acceleration OD."
                )

        return {
            "tier_I_literature_geometry": {
                "status": "complete",
                "sources": [
                    "published Δv (literature)",
                    "Step 007 geometry envelope",
                    "Step 039 cross-catalog classification",
                ],
            },
            "tier_II_uncertainty_aware_null": {
                "status": "complete_for_juno",
                "step039_juno": (
                    "deterministic fixed-amplitude warning at pooled β; "
                    "uncertainty-aware raw layer: true_null when random-effects scatter propagated"
                ),
            },
            "tier_III_perigee_trk_minimal_od_mm_s": {
                "status": tier_iii_status,
                "trk_commensurate_with_published_delta_v": trk_commensurate,
                "trk_off_epoch_archive": trk_off_epoch,
                "falsification_test_status": falsification_result.get("status"),
                "proxy_kind": proxy_kind,
                "detail": tier_iii_detail,
            },
            "horizons_public_ephemeris_stress_test": horizons_note,
            "od_absorption_not_required_for_juno_null": (
                "Step 039 classifies Juno as uncertainty-aware compatible with zero at "
                "random-effects amplitude scatter; OD-absorption is an unvalidated conjecture "
                "(Step 021 withholds F_OD; Step 012 synthetic quarantined)."
            ),
        }

    def _probe_juno_public_urls_optional(self) -> Optional[Dict[str, Any]]:
        """
        Optional live check of public archive entry points (no bulk download).
        Enable with ``TEP_030_PROBE_PUBLIC_URLS=1``.
        """
        if os.environ.get("TEP_030_PROBE_PUBLIC_URLS") != "1":
            return None
        self.logger.subsection("Public URL probe (TEP_030_PROBE_PUBLIC_URLS=1)")
        session = requests.Session()
        session.headers.update(
            {
                "User-Agent": (
                    "TEP-EFA-Pipeline/1.0 (Academic research; Juno 2013 public data audit)"
                )
            }
        )
        probes: List[Dict[str, Any]] = []
        for probe_id, url, purpose in self.JUNO_PUBLIC_URL_PROBE_TARGETS:
            entry: Dict[str, Any] = {
                "id": probe_id,
                "url": url,
                "purpose": purpose,
            }
            try:
                with session.get(
                    url, timeout=25, allow_redirects=True, stream=True
                ) as resp:
                    entry["method"] = "GET_stream"
                    entry["status_code"] = int(resp.status_code)
                    entry["final_url"] = str(resp.url)
                    # Read at most one chunk so directory listings/HTML start to flow.
                    next(resp.iter_content(1024), b"")
                self.logger.info(
                    f"  {probe_id}: HTTP {entry['status_code']} {entry['final_url']}"
                )
            except requests.RequestException as exc:
                entry["status_code"] = None
                entry["error"] = f"{type(exc).__name__}: {exc}"
                self.logger.warning(f"  {probe_id}: {entry['error']}")
            probes.append(entry)
        return {
            "status": "completed",
            "note": (
                "Single streamed GET (≤1 KiB read) per URL; does not download TRK. "
                "Records which public entry points respond for operator audits."
            ),
            "probes": probes,
        }

    def query_pds_inventory(self) -> Dict:
        """
        Query NASA PDS Radio Science Node for Juno Earth flyby data inventory.
        
        Returns metadata about available TRK-2-25 files without downloading.
        """
        self.logger.subsection("Querying NASA PDS Radio Science Node")
        
        # PDS data products for Earth flyby (2013-10-09)
        # The data is organized by year/month/day
        year = self.FLYBY_DATE.year
        month = self.FLYBY_DATE.month
        day = self.FLYBY_DATE.day
        
        inventory = {
            'mission': 'Juno',
            'flyby_date': self.FLYBY_DATE.isoformat(),
            'pds_collection': self.PDS_RN_BASE + self.JUNO_EARTH_FLYBY_COLLECTION,
            'expected_data_products': [
                {
                    'type': 'TRK-2-25 (ATDF)',
                    'description': 'Archival Tracking Data File - Raw Doppler',
                    'url_pattern': f'{self.PDS_RN_BASE}{self.JUNO_EARTH_FLYBY_COLLECTION}{year}/{month:02d}/',
                    'expected_files': [
                        f'jnoe_{year}{month:02d}{day:02d}_*.trk',
                        f'jnoe_{year}{month:02d}{day-1:02d}_*.trk',
                        f'jnoe_{year}{month:02d}{day+1:02d}_*.trk'
                    ]
                },
                {
                    'type': 'ODF (Orbit Data File)',
                    'description': 'Processed orbit data',
                    'url_pattern': f'{self.PDS_RN_BASE}{self.JUNO_EARTH_FLYBY_COLLECTION}{year}/{month:02d}/',
                    'note': 'May contain pre-processed data - use TRK-2-25 for raw Doppler proxy'
                }
            ],
            'dss_stations': ['DSS-24', 'DSS-25', 'DSS-54', 'DSS-34', 'DSS-63'],
            'frequency_bands': ['X-band (8.4 GHz)', 'Ka-band (32 GHz)'],
            'access_method': 'HTTP download from PDS-RN',
            'data_volume_mb': 'Estimated 50-100 MB for 48-hour arc'
        }
        
        self.logger.info(f"PDS Collection: {inventory['pds_collection']}")
        self.logger.info(f"Expected data window: {self.ANALYSIS_WINDOW_HOURS} hours around perigee")
        self.logger.info(f"Target stations: {', '.join(inventory['dss_stations'])}")
        
        return inventory
    
    def download_trk225_data(self) -> Dict:
        """
        Download raw TRK-2-25 tracking data from NASA PDS.
        
        Tries multiple NASA data sources:
        1. PDS Radio Science Node (primary)
        2. PDS Geosciences Node
        3. JPL DSN Data Archive
        4. NASA CDDIS (if applicable)
        """
        self.logger.subsection("Acquiring Raw TRK-2-25 Data")

        ingest_manifest = ingest_mission_tracking(PROJECT_ROOT, "Juno_2013")
        self.logger.info(
            f"PDS ingest status: {ingest_manifest['status']} "
            f"({len(ingest_manifest.get('downloads', []))} files)"
        )

        existing_files = []
        for pattern in ('*.trk', '*.dat', '*.TNF', '*.tnf'):
            existing_files.extend(self.data_dir.rglob(pattern))
        existing_files = [
            path for path in existing_files
            if path.is_file() and is_trk234_archive(path)
        ]
        if existing_files:
            paths = [str(f) for f in existing_files]
            audit = self._trk_ingest_temporal_audit(paths)
            perigee_overlap = bool(audit.get("any_file_overlaps_flyby_window"))
            self.logger.info(f"Found {len(existing_files)} existing NJPL archive(s) on disk:")
            for f in existing_files:
                self.logger.info(f"  - {f.name}")
            self.raw_files = paths
            if perigee_overlap:
                self.logger.success(
                    "At least one local archive overlaps the 2013 Earth-flyby analysis window."
                )
                return {
                    "files_downloaded": len(existing_files),
                    "file_paths": paths,
                    "data_directory": str(self.data_dir),
                    "status": "success",
                    "source": "local_cache_perigee_window",
                    "perigee_window_matched": True,
                    "trk_ingest_temporal_audit": audit,
                }
            self.logger.warning(
                "Local NJPL archive(s) present but none overlap the 2013-10-09 perigee window "
                "(outer-cruise OCRU or missing .LBL). Proceeding as format-validation / proxy only."
            )
            return {
                "files_downloaded": len(existing_files),
                "file_paths": paths,
                "data_directory": str(self.data_dir),
                "status": "success_format_validation_only",
                "source": "local_cache_outside_perigee_window",
                "perigee_window_matched": False,
                "trk_ingest_temporal_audit": audit,
                "note": (
                    "TRK-2-34 parser and Hz pairwise proxy may run; do not treat as perigee Doppler OD."
                ),
            }
        
        # Calculate date range
        start_date = self.FLYBY_DATE - timedelta(hours=self.ANALYSIS_WINDOW_HOURS/2)
        end_date = self.FLYBY_DATE + timedelta(hours=self.ANALYSIS_WINDOW_HOURS/2)
        
        self.logger.info(f"Data window: {start_date.date()} to {end_date.date()}")
        self.logger.info(f"Perigee: {self.FLYBY_DATE} UTC")
        
        downloaded_files = []
        
        # Multiple PDS sources to try
        pds_sources = [
            ("PDS-RN Earth Flyby", f"{self.PDS_RN_BASE}{self.JUNO_EARTH_FLYBY_COLLECTION}"),
            ("PDS-RN Jupiter Cruise", f"{self.PDS_RN_BASE}jno-j-rss-1-edr/"),
        ]
        
        session = requests.Session()
        session.headers.update({
            'User-Agent': 'TEP-EFA-Pipeline/1.0 (Academic Research; Contact: tep-efa@example.edu)'
        })
        
        for source_name, base_collection in pds_sources:
            self.logger.info(f"Trying {source_name}...")
            
            current_date = start_date
            while current_date <= end_date:
                year = current_date.year
                month = current_date.month
                day = current_date.day
                
                # URL patterns
                urls_to_try = [
                    f"{base_collection}{year}/{month:02d}/{day:02d}/",
                    f"{base_collection}{year}/{month:02d}/",
                    f"{base_collection}data/{year}/{month:02d}/{day:02d}/",
                ]
                
                for url in urls_to_try:
                    try:
                        self.logger.info(f"  Checking: {url}")
                        response = session.get(url, timeout=10, allow_redirects=True)
                        
                        if response.status_code == 200:
                            # Look for tracking data files
                            file_patterns = [
                                r'href="([^"]*\.trk)"',
                                r'href="([^"]*\.dat)"',
                                r'href="([^"]*\.TNF)"',
                                r'href="([^"]*jnoe[^"]*\.[^"]*)"',  # Juno Earth flyby pattern
                            ]
                            
                            for pattern in file_patterns:
                                files = re.findall(pattern, response.text, re.IGNORECASE)
                                for filename in files:
                                    file_url = url + filename
                                    local_path = self.data_dir / filename
                                    
                                    if not local_path.exists():
                                        self.logger.info(f"    Downloading: {filename}")
                                        try:
                                            file_response = session.get(file_url, timeout=30)
                                            if file_response.status_code == 200:
                                                with open(local_path, 'wb') as f:
                                                    f.write(file_response.content)
                                                downloaded_files.append(str(local_path))
                                                self.logger.info(f"      Saved: {len(file_response.content)} bytes")
                                        except Exception as e:
                                            self.logger.warning(f"      Error: {e}")
                                    else:
                                        downloaded_files.append(str(local_path))
                                        
                    except requests.exceptions.ConnectionError as e:
                        self.logger.warning(f"  Connection failed (network restriction): {e}")
                    except requests.exceptions.Timeout as e:
                        self.logger.warning(f"  Timeout accessing {url}: {e}")
                    except (ValueError, KeyError, IndexError, AttributeError) as e:
                        self.logger.debug(f"  Error: {e}")
                
                current_date += timedelta(days=1)
        
        # Summary
        self.raw_files = downloaded_files
        
        if downloaded_files:
            audit = self._trk_ingest_temporal_audit(downloaded_files)
            return {
                'files_downloaded': len(downloaded_files),
                'file_paths': downloaded_files,
                'data_directory': str(self.data_dir),
                'status': 'success',
                'source': 'pds_download',
                'perigee_window_matched': bool(audit.get('any_file_overlaps_flyby_window')),
                'trk_ingest_temporal_audit': audit,
            }
        else:
            ocru_paths = self._download_atmospheres_ocru_tnfs(session)
            if ocru_paths:
                self.raw_files = [str(p) for p in ocru_paths]
            paths = [str(p) for p in ocru_paths]
            audit = self._trk_ingest_temporal_audit(paths)
            perigee_overlap = bool(audit.get("any_file_overlaps_flyby_window"))
            return {
                'files_downloaded': len(ocru_paths),
                'file_paths': paths,
                'data_directory': str(self.data_dir),
                'status': 'success' if perigee_overlap else 'success_format_validation_only',
                'source': 'atmospheres_jnogrv_0001_tnf',
                'perigee_window_matched': perigee_overlap,
                'trk_ingest_temporal_audit': audit,
                'note': (
                    'NMSU PDS Atmospheres JUNO-J-RSS-1-OCRU-V1.0 TNF subset; '
                    'see DOWNLOAD_INSTRUCTIONS.txt for perigee-window caveats.'
                ),
            }
            return {
                'files_downloaded': 0,
                'file_paths': [],
                'data_directory': str(self.data_dir),
                'status': 'no_data_accessible',
                'note': (
                    'PDS download requires network access or '
                    'TEP_030_ALLOW_ATMOSPHERES_OCRU_MIRROR=1 for the verified Atmospheres '
                    'TNF index. Manual download: DOWNLOAD_INSTRUCTIONS.txt'
                ),
            }

    def _filter_tnf_names_by_lbl_window(
        self,
        session: requests.Session,
        idx_url: str,
        tnf_names: List[str],
        max_n: int,
    ) -> List[str]:
        """
        Keep only ``.TNF`` basenames whose sibling ``.LBL`` START_TIME/STOP_TIME overlaps
        the Step 030 analysis window around ``FLYBY_DATE`` (``ANALYSIS_WINDOW_HOURS``).
        """
        flyby = self.FLYBY_DATE.replace(tzinfo=timezone.utc)
        win_lo = flyby - timedelta(hours=self.ANALYSIS_WINDOW_HOURS / 2)
        win_hi = flyby + timedelta(hours=self.ANALYSIS_WINDOW_HOURS / 2)
        scan_cap = int(os.environ.get("TEP_030_ATMOSPHERES_LBL_SCAN_CAP", "500"))
        overlap: List[str] = []
        scanned = 0
        for name in sorted(set(tnf_names)):
            if len(overlap) >= max_n:
                break
            if scanned >= scan_cap:
                self.logger.warning(
                    f"LBL time filter: hit TEP_030_ATMOSPHERES_LBL_SCAN_CAP={scan_cap} "
                    f"with {len(overlap)} overlap(s); raise cap or narrow index."
                )
                break
            scanned += 1
            stem = Path(name).stem
            lbl_url = idx_url + stem + ".LBL"
            try:
                lr = session.get(lbl_url, timeout=45)
            except requests.RequestException as exc:
                self.logger.debug(f"LBL GET {lbl_url}: {exc}")
                continue
            if lr.status_code != 200:
                continue
            try:
                span = parse_pds3_lbl_start_stop_utc(lr.text)
            except ValueError as exc:
                self.logger.debug(f"LBL parse {stem}.LBL: {exc}")
                continue
            if span is None:
                continue
            t0, t1 = span
            if t1 < win_lo or t0 > win_hi:
                continue
            overlap.append(name)
            self.logger.info(
                f"  LBL overlap: {name} ({t0.isoformat()}–{t1.isoformat()}) in flyby window"
            )
        return overlap

    def _download_atmospheres_ocru_tnfs(self, session: requests.Session) -> List[Path]:
        """
        Optional bulk ingest from the verified NMSU Atmospheres HTTPS mirror for
        JUNO-J-RSS-1-OCRU-V1.0 (NJPL-class TNF).

        Opt-in only: set environment variable ``TEP_030_ALLOW_ATMOSPHERES_OCRU_MIRROR=1``.
        Caps download count with ``TEP_030_ATMOSPHERES_MAX_TNF_FILES`` (default 30).
        Index URL: ``TEP_030_ATMOSPHERES_INDEX_URL`` (default class constant).
        When ``TEP_030_ATMOSPHERES_LBL_TIME_FILTER=1``, only TNFs whose sibling ``.LBL``
        interval overlaps the Juno 2013 flyby analysis window are kept; if none match,
        raises ``RuntimeError`` (strict; no fallback downloads).
        Otherwise prefers basenames containing ``2013`` when present; else lexicographic
        first names, with a volume ERRATA warning.
        """
        if os.environ.get("TEP_030_ALLOW_ATMOSPHERES_OCRU_MIRROR") != "1":
            return []
        max_n = int(os.environ.get("TEP_030_ATMOSPHERES_MAX_TNF_FILES", "30"))
        if max_n < 1:
            raise ValueError("TEP_030_ATMOSPHERES_MAX_TNF_FILES must be >= 1")

        self.logger.subsection(
            "Atmospheres OCRU TNF mirror (TEP_030_ALLOW_ATMOSPHERES_OCRU_MIRROR=1)"
        )
        dest = self.data_dir / "NASA_PDS"
        dest.mkdir(parents=True, exist_ok=True)

        idx_url = (
            os.environ.get(
                "TEP_030_ATMOSPHERES_INDEX_URL", self.ATMOSPHERES_JNO_OCRU_TNF_INDEX
            ).rstrip("/")
            + "/"
        )
        self.logger.info(f"GET index: {idx_url}")
        r = session.get(idx_url, timeout=60, allow_redirects=True)
        r.raise_for_status()
        hrefs = re.findall(
            r'href="([^"]+\.(?:TNF|tnf))"', r.text, flags=re.IGNORECASE
        )
        # De-duplicate, preserve order
        seen = set()
        tnf_names: List[str] = []
        for h in hrefs:
            name = h.split("/")[-1]
            if not name.lower().endswith(".tnf"):
                continue
            if name in seen:
                continue
            seen.add(name)
            tnf_names.append(name)

        if not tnf_names:
            self.logger.warning("Atmospheres index contained no .TNF href matches.")
            return []

        filter_lbl = os.environ.get("TEP_030_ATMOSPHERES_LBL_TIME_FILTER") == "1"
        if filter_lbl:
            overlap = self._filter_tnf_names_by_lbl_window(
                session, idx_url, tnf_names, max_n
            )
            if not overlap:
                flyby = self.FLYBY_DATE.replace(tzinfo=timezone.utc)
                win_lo = flyby - timedelta(hours=self.ANALYSIS_WINDOW_HOURS / 2)
                win_hi = flyby + timedelta(hours=self.ANALYSIS_WINDOW_HOURS / 2)
                raise RuntimeError(
                    "TEP_030_ATMOSPHERES_LBL_TIME_FILTER=1: no product .LBL START_TIME/"
                    "STOP_TIME overlaps the Juno 2013 Earth-flyby analysis window "
                    f"{win_lo.isoformat()}–{win_hi.isoformat()} (UTC). The jnogrv_0001 OCRU "
                    "volume may contain no perigee tracking; use JNO-E-RSS Earth-phase "
                    "holdings or another index (TEP_030_ATMOSPHERES_INDEX_URL)."
                )
            chosen = overlap
            self.logger.info(
                f"LBL time filter: {len(chosen)} TNF(s) overlap the flyby analysis window."
            )
        else:
            preferred = [n for n in tnf_names if "2013" in n]
            if preferred:
                chosen = preferred[:max_n]
                self.logger.info(
                    f"Selected {len(chosen)} TNF(s) with '2013' in filename (cap {max_n})."
                )
            else:
                chosen = sorted(tnf_names)[:max_n]
                self.logger.warning(
                    "No TNF basename contains '2013'; using lexicographically first "
                    f"{len(chosen)} file(s). Volume ERRATA: archive may start after 2013 "
                    "Earth flyby — pairwise statistics are not asserted as perigee Doppler."
                )

        out: List[Path] = []
        for name in chosen:
            file_url = idx_url + name
            local_path = dest / name
            if local_path.is_file() and is_trk234_archive(local_path):
                out.append(local_path)
                self.logger.info(f"  cache hit (NJPL): {name}")
                continue
            self.logger.info(f"  GET {file_url}")
            fr = session.get(file_url, timeout=120)
            fr.raise_for_status()
            local_path.write_bytes(fr.content)
            if not is_trk234_archive(local_path):
                local_path.unlink(missing_ok=True)
                raise RuntimeError(
                    f"Downloaded file failed NJPL TRK-2-34 magic check: {name}"
                )
            out.append(local_path)
            self.logger.info(f"  saved {len(fr.content)} bytes -> {local_path.name}")

        return out

    def parse_trk225_file(self, filepath: str) -> List[Dict]:
        """
        Parse TRK-2-34 / TRK-2-25-class NJPL SFDU archive for radiometric observables.

        Uses ``trk234.Reader.decode()`` so ``sfdu_list`` is populated. Many archival
        TNFs expose ``ramp_freq`` / phase CHDO fields rather than a populated
        ``doppler`` attribute; in that case we fall back to
        :func:`scripts.utils.trk234_extract.extract_trk234_measurements` and carry
        ``ramp_freq_hz`` for a pairwise **frequency-difference** proxy (not mm/s).
        """
        measurements = []
        
        try:
            import trk234
            
            reader = trk234.Reader(filepath)
            reader.decode()
            
            for sfdu in reader.sfdu_list:
                if hasattr(sfdu, 'trk_chdo') and sfdu.trk_chdo is not None:
                    trk = sfdu.trk_chdo
                    
                    meas = {
                        'source_file': Path(filepath).name,
                        'sfdu_type': getattr(sfdu.label, 'sfdutype', None),
                        'observable_basis': 'trk234_chdo',
                    }
                    
                    if hasattr(trk, 'doppler') and trk.doppler is not None:
                        meas['doppler_hz'] = float(trk.doppler)
                    
                    if hasattr(trk, 'recvtime') and trk.recvtime is not None:
                        meas['timestamp'] = str(trk.recvtime)
                    
                    if hasattr(trk, 'rxtonefreq') and trk.rxtonefreq is not None:
                        meas['frequency_hz'] = float(trk.rxtonefreq)
                    
                    if hasattr(trk, 'dss_id') and trk.dss_id is not None:
                        meas['station'] = f"DSS-{trk.dss_id}"
                    
                    if 'doppler_hz' in meas:
                        measurements.append(meas)

            if measurements:
                return measurements

            extracted = extract_trk234_measurements(Path(filepath))
            for row in extracted:
                su = row.get("station_uplink")
                station = f"DSS-{int(su)}" if su is not None else "UNKNOWN"
                rf = row.get("ramp_freq_hz")
                if rf is None:
                    continue
                measurements.append(
                    {
                        "source_file": row.get("source_file", Path(filepath).name),
                        "sfdu_index": row.get("sfdu_index"),
                        "timestamp": row.get("timestamp"),
                        "station": station,
                        "ramp_freq_hz": float(rf),
                        "observable_basis": "trk234_ramp_freq_hz",
                        "band_uplink": row.get("band_uplink"),
                        "band_downlink": row.get("band_downlink"),
                    }
                )
            if measurements:
                self.logger.info(
                    f"  {Path(filepath).name}: using TRK-2-34 ramp_freq_hz extract "
                    f"({len(measurements)} rows; no CHDO doppler populated)."
                )

        except ImportError:
            self.logger.error("trk234 library not available. TRK-2-25 parsing requires this library for scientific integrity.")
            self.logger.error("Install trk234: pip install trk234")
            raise RuntimeError("trk234 library required for TRK-2-25 parsing. No fallback available.")
            
        except (OSError, ValueError, AttributeError) as e:
            self.logger.error(f"Error parsing {filepath}: {e}")
            
        return measurements
    
    def _basic_trk_parse(self, filepath: str) -> List[Dict]:
        """Basic TRK-2-25 parsing without external library."""
        measurements = []
        
        try:
            with open(filepath, 'rb') as f:
                data = f.read()
                
            # Look for Doppler data patterns in binary
            # This is a simplified parser - full TRK-2-25 requires detailed format knowledge
            
            # Find SFDU markers (typically start with specific byte patterns)
            sfdu_markers = [m.start() for m in re.finditer(b'SFDU', data)]
            
            for i, marker in enumerate(sfdu_markers[:100]):  # Limit to first 100 SFDUs
                try:
                    # Extract timestamp and Doppler from SFDU
                    # This is format-dependent and would need refinement
                    chunk = data[marker:marker+256]
                    
                    # Look for time stamp (usually near beginning)
                    # and Doppler value (usually floating point)
                    
                    meas = {
                        'source_file': Path(filepath).name,
                        'sfdu_index': i,
                        'marker_position': marker,
                        'raw_bytes': chunk[:64].hex()
                    }
                    
                    measurements.append(meas)
                    
                except (ValueError, KeyError, AttributeError):
                    continue
                    
        except Exception as e:
            self.logger.error(f"Basic parsing failed for {filepath}: {e}")
            
        return measurements
    
    def extract_doppler_pair_residual_proxy(self, doppler_data: List[Dict]) -> Dict:
        """
        Build a per-station pairwise-difference proxy (not orbit determination).

        When CHDO ``doppler`` is populated, sequential ΔDoppler is converted with
        Δv ≈ c Δf/f to a **mm/s** proxy. When only ``ramp_freq_hz`` is available
        (typical NJPL TNF via ``extract_trk234_measurements``), the proxy is
        **per-station Δ ramp frequency in Hz** — not commensurate with the mm/s
        falsification gate (see :meth:`compute_falsification_test`).
        """
        self.logger.subsection("Pairwise radiometric proxy (not orbit determination)")
        
        if not doppler_data:
            return {
                'status': 'no_data',
                'message': 'No Doppler data available for analysis',
                'residual_series': [],
            }
        
        use_doppler = any('doppler_hz' in m for m in doppler_data)
        use_ramp = (not use_doppler) and any(
            m.get('ramp_freq_hz') is not None for m in doppler_data
        )
        if use_doppler:
            self.logger.info(f"Input: {len(doppler_data)} raw measurements (Doppler CHDO path)")
        elif use_ramp:
            self.logger.info(
                f"Input: {len(doppler_data)} raw measurements (ramp_freq_hz path; "
                "pairwise deltas are Hz, not mm/s)"
            )
        else:
            self.logger.info(f"Input: {len(doppler_data)} raw measurements (no usable fields)")
        
        stations = {}
        for meas in doppler_data:
            station = meas.get('station', 'UNKNOWN')
            stations.setdefault(station, []).append(meas)
        
        self.logger.info(f"Data from {len(stations)} stations: {list(stations.keys())}")
        
        residuals = []
        proxy_kind = 'doppler_pair_mm_s'
        analysis_method = 'doppler_pair_residual_proxy'
        analysis_note = (
            'Sequential per-station ΔDoppler converted with Δv ≈ c Δf/f; '
            'not post-fit OD residuals.'
        )
        
        for station_name, station_data in stations.items():
            station_data_sorted = sorted(
                station_data, key=lambda x: x.get('timestamp', '')
            )
            
            for i in range(1, len(station_data_sorted)):
                prev = station_data_sorted[i - 1]
                curr = station_data_sorted[i]
                
                if use_doppler:
                    if 'doppler_hz' not in prev or 'doppler_hz' not in curr:
                        continue
                    doppler_diff = curr['doppler_hz'] - prev['doppler_hz']
                    if 'frequency_hz' not in curr:
                        continue
                    freq = curr['frequency_hz']
                    c = 299792458  # m/s
                    velocity_mm_s = (doppler_diff / freq) * c * 1000  # mm/s
                    residuals.append({
                        'station': station_name,
                        'timestamp': curr.get('timestamp'),
                        'doppler_diff_hz': doppler_diff,
                        'velocity_mm_s': velocity_mm_s,
                        'proxy_kind': proxy_kind,
                    })
                elif use_ramp:
                    pf = prev.get('ramp_freq_hz')
                    cf = curr.get('ramp_freq_hz')
                    if pf is None or cf is None:
                        continue
                    delta_hz = float(cf) - float(pf)
                    residuals.append({
                        'station': station_name,
                        'timestamp': curr.get('timestamp'),
                        'ramp_freq_delta_hz': delta_hz,
                        'proxy_kind': 'ramp_freq_pair_delta_hz',
                    })
        
        self.logger.info(f"Computed {len(residuals)} pairwise-difference points")
        
        perigee_residuals = self._analyze_perigee_passage(residuals)
        
        if use_ramp:
            proxy_kind = 'ramp_freq_pair_delta_hz'
            analysis_method = 'ramp_freq_pair_delta_hz_proxy'
            analysis_note = (
                'Sequential per-station difference of archival ramp_freq_hz (Hz); '
                'not Doppler, not mm/s, not OD residuals; not gated against 0.08 mm/s.'
            )
        
        return {
            'status': 'success',
            'analysis_method': analysis_method,
            'analysis_method_note': analysis_note,
            'proxy_kind': proxy_kind,
            'commensurate_with_literature_delta_v_mm_s': False,
            'n_stations': len(stations),
            'n_residuals': len(residuals),
            'perigee_analysis': perigee_residuals,
            'residual_series': residuals,
            'intended_minimal_od_when_bound': {
                'gravity_field': 'EGM-96 (10×10)',
                'empirical_accelerations': 'DISABLED',
                'outlier_rejection': 'Disabled or 5σ (document in OD run)',
                'doppler_smoothing': 'Raw',
                'estimation_params': 'Initial state (6) + SRP (1)',
                'software_class': 'MONTE / JPL ODP-class batch least squares',
            },
        }
    
    def _analyze_perigee_passage(self, residuals: List[Dict]) -> Dict:
        """Analyze residuals around perigee passage for TEP signal."""
        
        perigee_window_start = self.FLYBY_DATE - timedelta(hours=2)
        perigee_window_end = self.FLYBY_DATE + timedelta(hours=2)
        
        perigee_residuals = []
        
        for res in residuals:
            if res.get('timestamp'):
                try:
                    ts = datetime.fromisoformat(res['timestamp'].replace('Z', '+00:00'))
                    
                    if perigee_window_start <= ts <= perigee_window_end:
                        perigee_residuals.append(res)
                        
                except (ValueError, KeyError, AttributeError):
                    continue
        
        if not perigee_residuals:
            return {
                'n_points': len(perigee_residuals),
                'message': 'No perigee residuals computed'
            }

        ramp_deltas = [
            r['ramp_freq_delta_hz'] for r in perigee_residuals
            if 'ramp_freq_delta_hz' in r
        ]
        velocities = [
            r['velocity_mm_s'] for r in perigee_residuals
            if 'velocity_mm_s' in r
        ]

        if ramp_deltas and not velocities:
            mean_v = float(np.mean(ramp_deltas))
            std_v = float(np.std(ramp_deltas))
            return {
                'n_points': len(perigee_residuals),
                'mean_ramp_freq_delta_hz': mean_v,
                'std_ramp_freq_delta_hz': std_v,
                'proxy_kind': 'ramp_freq_pair_delta_hz',
                'window_start': perigee_window_start.isoformat(),
                'window_end': perigee_window_end.isoformat(),
                'note': (
                    'Perigee-window mean is in Hz (ramp frequency pairwise differences); '
                    'not comparable to the 0.08 mm/s Doppler-proxy gate.'
                ),
            }

        if velocities:
            mean_v = float(np.mean(velocities))
            std_v = float(np.std(velocities))
            return {
                'n_points': len(perigee_residuals),
                'mean_velocity_mm_s': mean_v,
                'std_velocity_mm_s': std_v,
                'proxy_kind': 'doppler_pair_mm_s',
                'window_start': perigee_window_start.isoformat(),
                'window_end': perigee_window_end.isoformat(),
                'proxy_abs_mean_exceeds_0p08_mm_s_gate': abs(mean_v) > self.FALSIFICATION_THRESHOLD,
            }

        return {
            'n_points': len(perigee_residuals),
            'message': 'No perigee residuals with velocity_mm_s or ramp_freq_delta_hz',
        }
    
    def export_juno_residual_series_for_042(self, residuals: List[Dict]) -> Optional[Path]:
        """
        Export pairwise-Doppler proxy (mm/s) timestamps for Step 042 correlation.

        Uses only rows with valid timestamp and velocity_mm_s from real TRK processing.
        Also writes ``results/step030_juno_pairwise_residual_series.json`` (full rows)
        for reproducibility.
        """
        payload = build_juno_pairwise_042_sidecar_from_rows(residuals)
        if payload is None:
            self.logger.warning("Residual export: fewer than 3 valid timestamps; skipping 042 sidecar.")
            return None
        out_path = write_juno_042_sidecar(PROJECT_ROOT, payload)
        arch_path = save_juno_pairwise_step030_archive(PROJECT_ROOT, residuals)
        self.logger.info(
            f"Archived {len(residuals)} pairwise proxy rows for reproducibility: "
            f"{arch_path.relative_to(PROJECT_ROOT)}"
        )
        return out_path

    def _refresh_step042_after_juno_sidecar(self) -> None:
        """Re-run Step 042 so Juno time series correlations see the new sidecar."""
        traj = PROJECT_ROOT / "results" / "step038_trajectory_series.json"
        if not traj.is_file():
            self.logger.warning(
                "Skipping Step 042 refresh: missing results/step038_trajectory_series.json "
                "(run Step 004 then Step 038 first)."
            )
            return
        script = PROJECT_ROOT / "scripts" / "steps" / "step_042_time_resolved_cosmography.py"
        self.logger.subsection("Refreshing Step 042 after Juno sidecar export")
        proc = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(PROJECT_ROOT),
        )
        if proc.returncode != 0:
            self.logger.warning(
                f"Step 042 exited with code {proc.returncode}; "
                "see logs/step_042_time_resolved_cosmography.log"
            )
        else:
            self.logger.success("Step 042 refresh completed.")

    def compute_falsification_test(self, proxy_results: Dict) -> Dict:
        """
        Gate the perigee-window **pairwise-Doppler proxy** mean against a fixed
        0.08 mm/s repository threshold.

        This uses ``proxy_results['perigee_analysis']`` from
        :meth:`extract_doppler_pair_residual_proxy`; it does not substitute
        for a full minimal-OD residual time series and is **not** dimensionally
        aligned with published post-OD flyby Δv anomalies.
        """
        self.logger.subsection("PROXY GATE RESULT (not OD-scale TEP falsification)")
        
        perigee = proxy_results.get('perigee_analysis', {})
        basis = proxy_results.get('analysis_method', 'unknown')
        proxy_kind = proxy_results.get('proxy_kind') or perigee.get('proxy_kind')

        if 'mean_ramp_freq_delta_hz' in perigee:
            return {
                'status': 'NOT_COMMENSURATE_RAMP_FREQ_PROXY_NO_MM_S_GATE',
                'analysis_basis': basis,
                'reason': (
                    'Perigee-window statistic is mean pairwise Δ ramp_freq_hz (Hz), '
                    'not a Doppler-derived mm/s proxy; the 0.08 mm/s gate does not apply.'
                ),
                'perigee_hz_summary': {
                    'mean_ramp_freq_delta_hz': perigee.get('mean_ramp_freq_delta_hz'),
                    'std_ramp_freq_delta_hz': perigee.get('std_ramp_freq_delta_hz'),
                    'n_points': perigee.get('n_points'),
                },
                'recommendation': (
                    'Use Earth-flyby products with CHDO Doppler populated, or attach '
                    'MONTE-class OD for mm/s-commensurate residuals.'
                ),
                'claims_tep_od_scale_falsified': None,
            }

        if proxy_kind == 'ramp_freq_pair_delta_hz':
            return {
                'status': 'inconclusive',
                'analysis_basis': basis,
                'reason': (
                    'Ramp-frequency Hz pairwise proxy is defined, but no residuals fall '
                    'in the ±2 h UTC window around 2013-10-09 19:21 (Earth flyby), or Hz '
                    'summary is missing. Ingested arc may be outside perigee (e.g. OCRU '
                    'outer cruise) — see DOWNLOAD_INSTRUCTIONS.txt.'
                ),
                'recommendation': (
                    'Overlap perigee UTC with the archive, or use Doppler-populated '
                    'Earth-flyby TRK products for the mm/s gate.'
                ),
                'claims_tep_od_scale_falsified': None,
            }

        if 'mean_velocity_mm_s' not in perigee:
            reason = 'No perigee velocity measurement available'
            if basis == 'ramp_freq_pair_delta_hz_proxy':
                reason = (
                    'No perigee-window ramp-frequency pairwise residuals (UTC window '
                    '±2 h of 2013-10-09 19:21 UTC), or Hz summary not populated.'
                )
            return {
                'status': 'inconclusive',
                'analysis_basis': basis,
                'reason': reason,
                'recommendation': 'Check data quality and reprocess',
                'claims_tep_od_scale_falsified': None,
            }
        
        measured_v = perigee['mean_velocity_mm_s']
        measured_std = perigee.get('std_velocity_mm_s', 0.02)
        
        # Statistical significance of the proxy mean (internal consistency only)
        z_score = abs(measured_v) / measured_std if measured_std > 0 else 0
        
        below_gate = abs(measured_v) < self.FALSIFICATION_THRESHOLD
        if below_gate:
            status = 'JUNO_PAIRWISE_PROXY_BELOW_0p08_MM_S_GATE'
            result = {
                'status': status,
                'analysis_basis': basis,
                'conclusion': (
                    f'|mean proxy| = {abs(measured_v):.3f} mm/s < gate '
                    f'({self.FALSIFICATION_THRESHOLD:.3f} mm/s).'
                ),
                'interpretation': (
                    'The perigee-window pairwise-Doppler proxy mean lies below the '
                    'repository 0.08 mm/s gate. This is **not** equivalent to showing '
                    'that a minimal-OD Δv anomaly is absent at the literature scale; '
                    'the proxy is not OD-commensurate.'
                ),
                'measured_velocity_mm_s': float(measured_v),
                'falsification_threshold_mm_s': self.FALSIFICATION_THRESHOLD,
                'z_score': float(z_score),
                'reference_tep_scale_mm_s': self.TEP_PREDICTED_SIGNAL,
                'claims_tep_od_scale_falsified': False,
                'commensurate_with_literature_delta_v_mm_s': False,
            }
        else:
            status = 'JUNO_PAIRWISE_PROXY_ABOVE_0p08_MM_S_GATE'
            result = {
                'status': status,
                'analysis_basis': basis,
                'conclusion': (
                    f'|mean proxy| = {abs(measured_v):.3f} mm/s ≥ gate '
                    f'({self.FALSIFICATION_THRESHOLD:.3f} mm/s).'
                ),
                'interpretation': (
                    'The proxy mean exceeds the 0.08 mm/s gate. This still does not, '
                    'by itself, establish minimal-OD recovery of a physical Δv anomaly; '
                    'attach MONTE-class residuals for OD-scale inference.'
                ),
                'measured_velocity_mm_s': float(measured_v),
                'predicted_velocity_mm_s': self.TEP_PREDICTED_SIGNAL,
                'z_score': float(z_score),
                'agreement_with_od_scale_prediction': (
                    'Comparable magnitude'
                    if abs(abs(measured_v) - abs(self.TEP_PREDICTED_SIGNAL)) < 1.0
                    else 'Not comparable (proxy ≠ OD residual)'
                ),
                'claims_tep_od_scale_falsified': False,
                'commensurate_with_literature_delta_v_mm_s': False,
            }
        
        return result

    def compute_matched_filter_tep_template(
        self,
        residuals: List[Dict],
        perigee_altitude_km: float = 817.0,
        perigee_velocity_km_s: float = 14.79,
        lambda_tep_km: float = 4000.0,
        earth_radius_km: float = 6371.0,
    ) -> Dict:
        """
        Time-domain matched-filter test: correlate pairwise residuals with a
        predicted TEP impulse template centred at perigee.

        The scalar force peaks at closest approach where the field gradient is
        strongest.  For a hyperbolic trajectory the altitude near perigee scales
        as h(t) ≈ h_p + (v_p²/2r_p)(t−t_p)², so the gradient factor
        exp(−h/λ_TEP) is approximately Gaussian in time with width
        σ_t = sqrt(r_p · λ_TEP) / v_p.

        This test does **not** require OD calibration; it asks only whether the
        raw pairwise residual series contains a structured signal whose time
        profile matches the expected TEP impulse shape.
        """
        self.logger.subsection("Matched-filter TEP template correlation")

        # --- Select residuals with usable scalar values and timestamps ---
        values = []
        times_s = []
        for res in residuals:
            ts_str = res.get("timestamp")
            if not ts_str:
                continue
            try:
                ts = datetime.fromisoformat(ts_str.replace("Z", "+00:00"))
            except (ValueError, TypeError):
                continue
            # Use velocity_mm_s when available, otherwise ramp_freq_delta_hz
            val = res.get("velocity_mm_s")
            if val is None:
                val = res.get("ramp_freq_delta_hz")
            if val is None or not np.isfinite(float(val)):
                continue
            dt_s = (ts - self.FLYBY_DATE).total_seconds()
            values.append(float(val))
            times_s.append(dt_s)

        if len(values) < 5:
            return {
                "status": "insufficient_data",
                "n_points": len(values),
                "reason": "Fewer than 5 residuals with timestamps and scalar values",
            }

        values_arr = np.asarray(values, dtype=float)
        times_arr = np.asarray(times_s, dtype=float)

        # --- Restrict to ±30 min around perigee for the matched-filter window ---
        window_s = 30.0 * 60.0
        mask = np.abs(times_arr) <= window_s
        if mask.sum() < 5:
            return {
                "status": "insufficient_data",
                "n_points": int(mask.sum()),
                "reason": "Fewer than 5 residuals inside ±30 min perigee window",
            }

        t_win = times_arr[mask]
        v_win = values_arr[mask]

        # --- Build Gaussian TEP template ---
        r_p_m = (earth_radius_km + perigee_altitude_km) * 1000.0
        lambda_m = lambda_tep_km * 1000.0
        v_p_m_s = perigee_velocity_km_s * 1000.0
        sigma_t = float(np.sqrt(r_p_m * lambda_m)) / v_p_m_s  # seconds

        # Normalised template (zero mean, unit variance over the window)
        template = np.exp(-0.5 * (t_win / sigma_t) ** 2)
        template = template - np.mean(template)
        template_norm = np.linalg.norm(template)
        if template_norm == 0.0:
            return {
                "status": "degenerate_template",
                "sigma_t_seconds": sigma_t,
                "reason": "Template has zero norm",
            }
        template = template / template_norm

        # --- Correlation and simple SNR ---
        v_centred = v_win - np.mean(v_win)
        corr_coeff = float(np.dot(v_centred, template))
        # Pearson r via projection (template is already normalised)
        v_norm = np.linalg.norm(v_centred)
        pearson_r = corr_coeff / v_norm if v_norm > 0 else 0.0

        # Approximate p-value: t-statistic for Pearson r on n-2 d.o.f.
        n_eff = len(v_win)
        if n_eff > 2 and abs(pearson_r) < 1.0:
            t_stat = pearson_r * np.sqrt((n_eff - 2) / (1.0 - pearson_r ** 2))
            # Two-tailed p from t-distribution CDF
            from scipy import stats as sp_stats
            p_value = float(2.0 * (1.0 - sp_stats.t.cdf(abs(t_stat), n_eff - 2)))
        else:
            p_value = 1.0

        # Simple peak SNR: ratio of template-weighted sum to RMS of off-template
        # background (using a ±2σ exclusion band around perigee)
        bg_rms = None
        off_mask = np.abs(t_win) > 2.0 * sigma_t
        if off_mask.sum() >= 3:
            bg_rms = float(np.std(v_win[off_mask]))
        else:
            bg_rms = float(np.std(v_win))
        peak_snr = corr_coeff / bg_rms if bg_rms > 0 else 0.0

        self.logger.info(f"Matched-filter window: ±{window_s/60:.0f} min, n = {n_eff}")
        self.logger.info(f"Template σ_t = {sigma_t:.1f} s ({sigma_t/60:.1f} min)")
        self.logger.info(f"Pearson r = {pearson_r:+.3f}, p = {p_value:.3f} (n={n_eff})")
        self.logger.info(f"Peak SNR = {peak_snr:+.2f}")

        return {
            "status": "computed",
            "n_points": n_eff,
            "window_seconds": window_s,
            "sigma_t_seconds": sigma_t,
            "perigee_altitude_km": perigee_altitude_km,
            "perigee_velocity_km_s": perigee_velocity_km_s,
            "lambda_tep_km": lambda_tep_km,
            "pearson_r": pearson_r,
            "p_value": p_value,
            "peak_snr": peak_snr,
            "background_rms": bg_rms,
            "interpretation": (
                "Correlation between pairwise residuals and a Gaussian TEP impulse "
                "template centred at perigee.  |r| > 0.3 with p < 0.05 would suggest "
                "a structured perigee signal consistent with TEP; |r| ≈ 0 supports "
                "absence of such structure in the raw proxy data."
            ),
        }

    def compute_archive_presmoothing_sensitivity(
        self,
        measurements: List[Dict],
        window_samples: Tuple[int, ...] = (1, 5, 11, 31, 61),
    ) -> Dict:
        """
        Observable-domain sensitivity: presmooth archival ramp (or Doppler) per
        station, then re-form sequential pairwise differences.

        Interprets Occam vs OD-suppression only at the level of "does heavy averaging
        of the radiometric observable before differencing materially change this
        proxy on this archive?". Empirical accelerations and batch least squares are
        not modeled.
        """
        self.logger.subsection(
            "Archive presmoothing sensitivity (pairwise proxy, not OD)"
        )

        if not measurements:
            return {"status": "skipped", "reason": "no_measurements"}

        use_ramp = any(m.get("ramp_freq_hz") is not None for m in measurements)
        use_doppler = any("doppler_hz" in m for m in measurements) and not use_ramp

        if not use_ramp and not use_doppler:
            return {
                "status": "skipped",
                "reason": "no_ramp_freq_hz_or_doppler_hz",
            }

        basis = "ramp_freq_hz" if use_ramp else "doppler_hz"

        def _parse_ts(s: Optional[str]) -> Optional[datetime]:
            if not s:
                return None
            try:
                return datetime.fromisoformat(s.replace("Z", "+00:00"))
            except (ValueError, TypeError):
                return None

        stations: Dict[str, List[Dict]] = {}
        for m in measurements:
            st = m.get("station", "UNKNOWN")
            stations.setdefault(st, []).append(m)

        per_station: Dict[str, Dict] = {}
        for st_name, rows in stations.items():
            rows_sorted = sorted(rows, key=lambda x: str(x.get("timestamp", "")))
            ts_list: List[datetime] = []
            y_list: List[float] = []
            freq_hz: List[float] = []

            for row in rows_sorted:
                t = _parse_ts(row.get("timestamp"))
                if t is None:
                    continue
                if use_ramp:
                    rf = row.get("ramp_freq_hz")
                    if rf is None:
                        continue
                    y_list.append(float(rf))
                else:
                    if "doppler_hz" not in row:
                        continue
                    y_list.append(float(row["doppler_hz"]))
                    fh = row.get("frequency_hz")
                    freq_hz.append(float(fh) if fh is not None else float("nan"))

                ts_list.append(t)

            if len(y_list) < 10:
                continue

            t_sec = np.array(
                [(t - ts_list[0]).total_seconds() for t in ts_list], dtype=float
            )
            y = np.asarray(y_list, dtype=float)
            dt_med = float(np.median(np.diff(t_sec))) if len(t_sec) > 1 else 0.0

            win_stats: Dict[str, Dict[str, float]] = {}
            for W in window_samples:
                if W < 1:
                    continue
                if W == 1:
                    ys = y
                else:
                    ys = uniform_filter1d(y, size=int(W), mode="nearest")

                if use_ramp:
                    deltas = np.diff(ys)
                else:
                    # Sequential ΔDoppler → mm/s using downstream frequency (NaN-safe)
                    fh = np.asarray(freq_hz, dtype=float)
                    if fh.shape[0] != ys.shape[0]:
                        continue
                    dopp_d = np.diff(ys)
                    f_down = fh[1:]
                    c = 299792458.0
                    mask = np.isfinite(f_down) & (f_down != 0.0)
                    deltas = np.empty_like(dopp_d)
                    deltas[:] = np.nan
                    deltas[mask] = (dopp_d[mask] / f_down[mask]) * c * 1000.0

                fin = deltas[np.isfinite(deltas)]
                if fin.size < 3:
                    continue

                rms = float(np.sqrt(np.mean(fin**2)))
                m_abs = float(np.mean(np.abs(fin)))
                std = float(np.std(fin))
                win_stats[str(W)] = {
                    "n_pairwise": float(fin.size),
                    "rms": rms,
                    "mean_abs": m_abs,
                    "std": std,
                    "approx_median_dt_s": dt_med,
                    "approx_smoothing_span_s": float(W) * dt_med,
                }

            if not win_stats:
                continue

            raw_key = "1" if "1" in win_stats else str(window_samples[0])
            raw_mean_abs = win_stats.get(raw_key, {}).get("mean_abs")
            frac_change_vs_raw: Dict[str, Optional[float]] = {}
            for wk, stats in win_stats.items():
                if wk == raw_key or raw_mean_abs is None or raw_mean_abs == 0.0:
                    frac_change_vs_raw[wk] = None
                else:
                    frac_change_vs_raw[wk] = float(
                        (stats["mean_abs"] - raw_mean_abs) / raw_mean_abs
                    )

            # One-shot 3σ clip on W=1 pairwise deltas (mimics aggressive editing)
            clip_note: Optional[Dict[str, float]] = None
            if "1" in win_stats and use_ramp:
                ys1 = y
                d1 = np.diff(ys1)
                fin1 = d1[np.isfinite(d1)]
                if fin1.size >= 10:
                    mu = float(np.mean(fin1))
                    sig = float(np.std(fin1))
                    if sig > 0.0:
                        keep = np.abs(fin1 - mu) <= 3.0 * sig
                        n_drop = int(fin1.size - np.sum(keep))
                        fin1c = fin1[keep]
                        clip_note = {
                            "fraction_removed": float(n_drop) / float(fin1.size),
                            "mean_abs_before": float(np.mean(np.abs(fin1))),
                            "mean_abs_after_3sigma": float(np.mean(np.abs(fin1c)))
                            if fin1c.size
                            else float("nan"),
                        }

            per_station[st_name] = {
                "basis": basis,
                "n_points_arc": int(y.size),
                "median_sample_dt_s": dt_med,
                "by_uniform_presmooth_samples": win_stats,
                "fractional_change_in_mean_abs_vs_W1": frac_change_vs_raw,
                "three_sigma_clip_on_W1_deltas_hz": clip_note,
            }

        if not per_station:
            return {
                "status": "skipped",
                "reason": "insufficient_parsed_points_per_station",
                "basis": basis,
            }

        summary = (
            "Uniform presmoothing of the archival observable before pairwise differencing. "
            "Large |fractional_change_in_mean_abs_vs_W1| under W≈31–61 suggests the proxy "
            "is sensitive to averaging comparable to multi-second–minute counting; that "
            "supports observable-level plausibility of pipeline attenuation but does not "
            "substitute for MONTE-class OD with empirical states."
        )

        self.logger.info(
            f"Presmoothing sensitivity: {len(per_station)} station(s), "
            f"basis={basis}, windows={window_samples}"
        )

        return {
            "status": "computed",
            "basis": basis,
            "window_samples": [int(w) for w in window_samples],
            "per_station": per_station,
            "interpretation": summary,
        }

    def run_full_reanalysis(self) -> Dict:
        """Execute complete Juno 2013 raw DSN reanalysis with REAL DATA ONLY."""
        self.logger.header("STEP 030: JUNO 2013 RAW DSN REANALYSIS (CRITICAL TEST)")
        
        self.logger.info("="*70)
        self.logger.info("REQUIREMENT: REAL DSN DATA ONLY - NO SYNTHETIC SUBSTITUTES")
        self.logger.info("="*70)
        self.logger.info("Pairwise-Doppler proxy gate: ±0.08 mm/s on |mean proxy| (not OD-scale Δv).")
        self.logger.info("TEP reference scale +2.25 mm/s is not commensurate with this proxy without OD binding.")
        self.logger.info("="*70)
        
        # Step 1: Query PDS inventory
        inventory = self.query_pds_inventory()
        
        # Step 2: Download raw data
        download_result = self.download_trk225_data()
        
        ok_statuses = ("success", "success_format_validation_only")
        if download_result.get("status") not in ok_statuses or download_result["files_downloaded"] == 0:
            self.logger.warning("=" * 70)
            self.logger.warning("NO PERIGEE-WINDOW JUNO 2013 DSN ARCHIVE ON DISK")
            self.logger.warning("=" * 70)
            self.logger.warning(f"Download status: {download_result.get('status')}")
            self.logger.warning("Files found: 0")
            self.logger.info("")
            self.logger.info("MANUAL DOWNLOAD REQUIRED:")
            self.logger.info("1. Visit: https://pds.mcp.nasa.gov/portal/search")
            self.logger.info("2. Search: 'Juno Earth flyby 2013' or 'JNO-E-RSS-1-EDR'")
            self.logger.info("3. Download TRK-2-34/NJPL products whose .LBL overlaps 2013-10-08 to 2013-10-10 UTC")
            self.logger.info("4. Place in: data/raw/dsn_tracking/Juno_2013/")
            self.logger.info("5. Re-run: python scripts/steps/step_030_juno_reanalysis.py")
            self.logger.info("")
            self.logger.info("CONTACT: pds-rn@jpl.nasa.gov — 'Juno 2013 Earth Flyby TRK-2-25 Data Request'")
            self.logger.warning("=" * 70)

            horizons_pub = None
            try:
                horizons_pub = self._horizons_public_ephemeris_optional()
            except (FileNotFoundError, KeyError, ValueError, RuntimeError) as exc:
                self.logger.warning(f"Horizons public ephemeris batch skipped: {exc}")

            payload = {
                "status": "no_perigee_window_trk",
                "reason": "PERIGEE_WINDOW_TRK_REQUIRED",
                "message": "No NJPL archive on disk overlaps the 2013 Earth-flyby window",
                "resolution": "Manual perigee-window ingest from NASA PDS (see DOWNLOAD_INSTRUCTIONS.txt)",
                "instructions": str(self.data_dir / "DOWNLOAD_INSTRUCTIONS.txt"),
                "inventory": inventory,
                "download_attempt": download_result,
                "data_provenance": {"files": [], "n_measurements": 0},
            }
            if horizons_pub is not None:
                payload["horizons_public_ephemeris_batch"] = horizons_pub
            url_probe = self._probe_juno_public_urls_optional()
            if url_probe is not None:
                payload["juno_2013_public_url_probe"] = url_probe
            output_file = self.results_dir / "step030_juno_dsn_reanalysis.json"
            with open(output_file, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, default=str)
            return payload

        perigee_matched = bool(download_result.get("perigee_window_matched"))
        if perigee_matched:
            self.logger.success("=" * 70)
            self.logger.success(
                f"PERIGEE-WINDOW DSN ARCHIVE: {download_result['files_downloaded']} file(s)"
            )
            self.logger.success("=" * 70)
        else:
            self.logger.warning("=" * 70)
            self.logger.warning(
                f"FORMAT-VALIDATION ONLY: {download_result['files_downloaded']} NJPL file(s) "
                "outside 2013 flyby window"
            )
            self.logger.warning("=" * 70)
        for f in download_result["file_paths"]:
            self.logger.info(f"  - {Path(f).name}")
        
        # Step 3: Parse TRK-2-25 files
        self.logger.subsection("Parsing TRK-2-25 Data")
        all_measurements = []
        
        for filepath in self.raw_files:
            measurements = self.parse_trk225_file(filepath)
            all_measurements.extend(measurements)
            self.logger.info(f"  {Path(filepath).name}: {len(measurements)} measurements")
        
        self.doppler_data = all_measurements
        self.logger.info(f"Total measurements: {len(all_measurements)}")
        
        # STRICT CHECK: Must have measurements
        if len(all_measurements) == 0:
            self.logger.error("="*70)
            self.logger.error("DATA ERROR: Files present but no valid measurements extracted")
            self.logger.error("="*70)
            return {
                'status': 'FAILED',
                'reason': 'NO_VALID_MEASUREMENTS',
                'message': 'TRK-2-25 files present but parsing failed',
                'files': download_result['file_paths']
            }
        
        # Step 4: Doppler pair residual proxy (not orbit determination)
        proxy_results = self.extract_doppler_pair_residual_proxy(all_measurements)
        export_list = None
        if isinstance(proxy_results, dict):
            export_list = proxy_results.get("residual_series")
        if isinstance(proxy_results, dict) and "residual_series" in proxy_results:
            del proxy_results["residual_series"]

        if (
            isinstance(export_list, list)
            and len(export_list) >= 3
            and proxy_results.get("status") == "success"
            and proxy_results.get("proxy_kind") == "doppler_pair_mm_s"
        ):
            ep = self.export_juno_residual_series_for_042(export_list)
            if ep is not None:
                self.logger.success(
                    f"Time-series export for Step 042: {ep.relative_to(PROJECT_ROOT)} "
                    f"(archive: results/step030_juno_pairwise_residual_series.json)"
                )
                self._refresh_step042_after_juno_sidecar()
        elif (
            isinstance(export_list, list)
            and len(export_list) >= 3
            and proxy_results.get("status") == "success"
            and proxy_results.get("proxy_kind") == "ramp_freq_pair_delta_hz"
        ):
            ramp_path = save_juno_ramp_pairwise_step030_archive(PROJECT_ROOT, export_list)
            self.logger.success(
                f"Archived {len(export_list)} ramp-frequency pairwise rows (Hz, not 042 mm/s): "
                f"{ramp_path.relative_to(PROJECT_ROOT)}"
            )

        # Step 5: Gate statistic vs falsification threshold
        falsification_result = self.compute_falsification_test(proxy_results)

        # Step 5b: Matched-filter TEP template correlation (time-domain structure)
        matched_filter_result = {}
        if isinstance(export_list, list) and len(export_list) >= 5:
            matched_filter_result = self.compute_matched_filter_tep_template(export_list)
        else:
            matched_filter_result = {
                "status": "skipped",
                "reason": "Fewer than 5 residuals available for matched-filter test",
            }

        # Step 5c: Observable presmoothing sensitivity (real archive; not MONTE OD)
        archive_presmoothing = self.compute_archive_presmoothing_sensitivity(
            all_measurements
        )

        trk_temporal = self._trk_ingest_temporal_audit(download_result["file_paths"])
        if trk_temporal.get("all_evaluated_files_miss_flyby_window"):
            self.logger.warning("=" * 70)
            self.logger.warning(
                "TRK INGEST TEMPORAL AUDIT: every evaluated PDS3 label interval "
                "misses the Earth-flyby analysis window around "
                f"{self.FLYBY_DATE.isoformat()} UTC (±{self.ANALYSIS_WINDOW_HOURS/2:.0f} h). "
                "Pairwise ramp/Doppler proxy is not asserted as perigee-window TRK."
            )
            self.logger.warning("=" * 70)

        perigee_matched = bool(download_result.get("perigee_window_matched"))
        analysis_mode = (
            "perigee_window_trk"
            if perigee_matched
            else "format_validation_or_hz_proxy_only"
        )

        horizons_pub: Optional[Dict[str, Any]] = None
        try:
            horizons_pub = self._horizons_public_ephemeris_optional()
        except (FileNotFoundError, KeyError, ValueError, RuntimeError) as exc:
            self.logger.warning(f"Horizons public ephemeris batch skipped: {exc}")

        evidence_tier = self._build_evidence_tier_assessment(
            trk_temporal=trk_temporal,
            falsification_result=falsification_result,
            proxy_results=proxy_results,
            horizons_pub=horizons_pub,
            perigee_matched=perigee_matched,
        )

        # Compile final results
        final_results = {
            'step': '030_juno_reanalysis',
            'timestamp': datetime.now().isoformat(),
            'status': download_result.get('status', 'success'),
            'perigee_window_matched': perigee_matched,
            'analysis_mode': analysis_mode,
            'flyby_date': self.FLYBY_DATE.isoformat(),
            'data_source': 'NASA PDS Radio Science Node',
            'trk_ingest_temporal_audit': trk_temporal,
            'data_provenance': {
                'files': download_result['file_paths'],
                'n_measurements': len(all_measurements),
                'measurements': all_measurements[:100] if len(all_measurements) > 100 else all_measurements  # Limit output size
            },
            'inventory': inventory,
            'download': download_result,
            'evidence_tier_assessment': evidence_tier,
            'mission_dsn_minimal_od_status': (
                'not_performed_on_perigee_trk'
                if trk_temporal.get('all_evaluated_files_miss_flyby_window')
                else (
                    'perigee_proxy_only'
                    if perigee_matched
                    else 'format_validation_or_hz_proxy_only'
                )
            ),
            'data_processing': {
                'n_measurements': len(all_measurements),
                'trk_commensurate_with_published_delta_v': evidence_tier[
                    'tier_III_perigee_trk_minimal_od_mm_s'
                ]['trk_commensurate_with_published_delta_v'],
                'doppler_residual_proxy_analysis': proxy_results,
                'ramp_freq_pairwise_archive': (
                   str(self.results_dir / 'step030_juno_pairwise_ramp_freq_series.json')
                    if proxy_results.get('proxy_kind') == 'ramp_freq_pair_delta_hz'
                    and isinstance(export_list, list)
                    and len(export_list) >= 3
                    and proxy_results.get('status') == 'success'
                    else None
                ),
            },
            'falsification_test': falsification_result,
            'matched_filter_tep_template': matched_filter_result,
            'archive_presmoothing_sensitivity': archive_presmoothing,
        }

        if horizons_pub is not None:
            final_results["horizons_public_ephemeris_batch"] = horizons_pub
            vo = horizons_pub.get("velocity_only") or {}
            final_results["horizons_minimal_od_batch"] = {
                "status": vo.get("status"),
                "interpretation": (
                    "Legacy key: velocity-only subset of horizons_public_ephemeris_batch "
                    "(see velocity_only + range_and_velocity there)."
                ),
                "batch_least_squares": vo.get("batch_least_squares"),
                "scipy_least_squares": vo.get("scipy_least_squares"),
                "residuals_summary": vo.get("residuals_summary"),
                "arc": horizons_pub.get("arc"),
                "physics": horizons_pub.get("physics"),
            }
            self.logger.info(
                "Horizons public ephemeris batch: see horizons_public_ephemeris_batch "
                "(velocity_only + range_and_velocity joint fit when range_m present)."
            )

        url_probe = self._probe_juno_public_urls_optional()
        if url_probe is not None:
            final_results["juno_2013_public_url_probe"] = url_probe
        
        # Save results
        output_file = self.results_dir / 'step030_juno_dsn_reanalysis.json'
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(final_results, f, indent=2, default=str)
        
        self.logger.section("FINAL RESULTS SAVED")
        self.logger.info(f"Output: {output_file}")
        
        if 'status' in falsification_result:
            self.logger.info(f"Falsification test: {falsification_result['status']}")
            self.logger.info(f"Conclusion: {falsification_result.get('conclusion', 'N/A')}")
        
        return final_results


def main():
    """Execute Juno 2013 raw DSN reanalysis."""
    reanalysis = JunoDSNReanalysis()
    results = reanalysis.run_full_reanalysis()
    
    # Log summary
    logger = StepLogger("step_030_juno_reanalysis", PROJECT_ROOT)
    
    if results.get('falsification_test', {}).get('status') == 'JUNO_PAIRWISE_PROXY_BELOW_0p08_MM_S_GATE':
        logger.warning("="*70)
        logger.warning("JUNO: PAIRWISE-DOPPLER PROXY BELOW 0.08 mm/s GATE")
        logger.warning("="*70)
        logger.info("This is a proxy gate only; it does not assert OD-scale TEP falsification.")

    elif results.get('falsification_test', {}).get('status') == 'NOT_COMMENSURATE_RAMP_FREQ_PROXY_NO_MM_S_GATE':
        logger.info("="*70)
        logger.info("JUNO: RAMP-FREQUENCY PAIRWISE PROXY (Hz) — mm/s GATE NOT APPLIED")
        logger.info("="*70)
        logger.info(
            "Real TRK-2-34 data ingested; perigee statistic is Δ ramp_freq in Hz. "
            "See results/step030_juno_pairwise_ramp_freq_series.json when archived."
        )

    elif results.get('falsification_test', {}).get('status') == 'JUNO_PAIRWISE_PROXY_ABOVE_0p08_MM_S_GATE':
        logger.success("="*70)
        logger.success("JUNO: PAIRWISE-DOPPLER PROXY ABOVE 0.08 mm/s GATE")
        logger.success("="*70)
        logger.info("Interpret with MONTE-class minimal OD before claiming navigation-level recovery.")

    else:
        logger.info("="*70)
        logger.info("INCONCLUSIVE OR NON-GATED — see falsification_test.status")
        logger.info("="*70)
        fs = results.get('falsification_test', {}).get('status')
        basis_ft = results.get('falsification_test', {}).get('analysis_basis', '')
        if fs not in ('NOT_COMMENSURATE_RAMP_FREQ_PROXY_NO_MM_S_GATE',) and not (
            fs == 'inconclusive' and basis_ft == 'ramp_freq_pair_delta_hz_proxy'
        ):
            logger.info("Raw data download may be required.")
        logger.info("See data/raw/dsn_tracking/Juno_2013/DOWNLOAD_INSTRUCTIONS.txt")

    # Matched-filter summary
    mf = results.get('matched_filter_tep_template', {})
    if mf.get('status') == 'computed':
        logger.info("="*70)
        logger.info("MATCHED-FILTER TEP TEMPLATE")
        logger.info("="*70)
        logger.info(f"Pearson r = {mf.get('pearson_r', 0):+.3f}, p = {mf.get('p_value', 1):.3f}")
        logger.info(f"Peak SNR = {mf.get('peak_snr', 0):+.2f} (template σ_t = {mf.get('sigma_t_seconds', 0):.1f} s)")
        logger.info("|r| ≈ 0 with non-significant p supports absence of structured TEP signal in proxy.")
    
    logger.log_step_summary(0, "SUCCESS")
    status = results.get("status")
    if status == "FAILED":
        return 1
    if status == "no_perigee_window_trk":
        logger.warning(
            "No perigee-window TRK on disk; JSON written with Horizons/public probes only."
        )
    # Exit 0 whenever processing completed; falsification outcome is in JSON.
    return 0


if __name__ == '__main__':
    sys.exit(main())
