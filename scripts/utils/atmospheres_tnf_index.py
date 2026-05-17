"""NMSU PDS Atmospheres HTTPS TNF index browse (Juno OCRU and similar mirrors)."""

from __future__ import annotations

import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Optional

import requests

from scripts.utils.pds3_lbl_time import interval_overlaps_window, parse_pds3_lbl_start_stop_utc

DEFAULT_JUNO_OCRU_TNF_INDEX = (
    "https://pds-atmospheres.nmsu.edu/PDS/data/jnogrv_0001/DATA/TNF/"
)


def list_tnf_names_from_index_html(html: str) -> list[str]:
    hrefs = re.findall(r'href="([^"]+\.(?:TNF|tnf))"', html, flags=re.IGNORECASE)
    seen: set[str] = set()
    names: list[str] = []
    for href in hrefs:
        name = href.split("/")[-1]
        if not name.lower().endswith(".tnf") or name in seen:
            continue
        seen.add(name)
        names.append(name)
    return names


def fetch_tnf_index(
    session: requests.Session,
    index_url: str = DEFAULT_JUNO_OCRU_TNF_INDEX,
) -> tuple[str, list[str]]:
    idx = index_url.rstrip("/") + "/"
    response = session.get(idx, timeout=60, allow_redirects=True)
    response.raise_for_status()
    return idx, list_tnf_names_from_index_html(response.text)


def filter_tnf_names_by_lbl_window(
    session: requests.Session,
    index_url: str,
    tnf_names: list[str],
    window_start: datetime,
    window_end: datetime,
    *,
    max_matches: int = 30,
    scan_cap: int = 500,
) -> tuple[list[dict[str, Any]], int]:
    """
    Return metadata rows for TNFs whose sibling ``.LBL`` overlaps ``[window_start, window_end]``.
    """
    idx = index_url.rstrip("/") + "/"
    overlap_rows: list[dict[str, Any]] = []
    scanned = 0
    for name in sorted(set(tnf_names)):
        if len(overlap_rows) >= max_matches:
            break
        if scanned >= scan_cap:
            break
        scanned += 1
        stem = Path(name).stem
        lbl_url = idx + stem + ".LBL"
        try:
            lr = session.get(lbl_url, timeout=45)
        except requests.RequestException:
            continue
        if lr.status_code != 200:
            continue
        try:
            span = parse_pds3_lbl_start_stop_utc(lr.text)
        except ValueError:
            continue
        if span is None:
            continue
        t0, t1 = span
        if not interval_overlaps_window(t0, t1, window_start, window_end):
            continue
        overlap_rows.append(
            {
                "tnf_name": name,
                "lbl_url": lbl_url,
                "start_utc": t0.isoformat(),
                "stop_utc": t1.isoformat(),
            }
        )
    return overlap_rows, scanned


def probe_atmospheres_tnf_index(
    session: requests.Session,
    perigee: datetime,
    *,
    window_hours: float = 48.0,
    index_url: str = DEFAULT_JUNO_OCRU_TNF_INDEX,
    scan_cap: int = 500,
) -> dict[str, Any]:
    """Scan Atmospheres TNF index; report perigee-window overlaps (no download)."""
    if perigee.tzinfo is None:
        perigee = perigee.replace(tzinfo=timezone.utc)
    win_lo = perigee - timedelta(hours=window_hours / 2.0)
    win_hi = perigee + timedelta(hours=window_hours / 2.0)
    try:
        idx, names = fetch_tnf_index(session, index_url)
    except requests.RequestException as exc:
        return {
            "index_url": index_url,
            "status": "index_fetch_failed",
            "error": str(exc),
            "perigee_utc": perigee.isoformat(),
            "window_utc": f"{win_lo.isoformat()} to {win_hi.isoformat()}",
        }
    overlap, scanned = filter_tnf_names_by_lbl_window(
        session,
        idx,
        names,
        win_lo,
        win_hi,
        scan_cap=scan_cap,
    )
    return {
        "index_url": idx,
        "status": "complete",
        "tnf_listed": len(names),
        "lbl_scanned": scanned,
        "perigee_window_overlaps": len(overlap),
        "overlap_candidates": overlap,
        "perigee_utc": perigee.isoformat(),
        "window_utc": f"{win_lo.isoformat()} to {win_hi.isoformat()}",
        "interpretation": (
            "jnogrv_0001 OCRU volume is outer-cruise; perigee-window TNFs are rare. "
            "Zero overlaps is expected unless Earth-phase JNO-E-RSS products are ingested separately."
        ),
    }
