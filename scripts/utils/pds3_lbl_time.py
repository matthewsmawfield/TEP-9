"""Parse PDS3 product label START_TIME / STOP_TIME (including day-of-year UTC)."""

from __future__ import annotations

import re
from datetime import datetime, timezone
from typing import Optional, Tuple


def pds3_lbl_value_line(lbl_text: str, key: str) -> Optional[str]:
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(.+)$", lbl_text, re.MULTILINE | re.IGNORECASE)
    if not m:
        return None
    raw = m.group(1).strip()
    if "/*" in raw:
        raw = raw.split("/*", 1)[0].strip()
    return raw.strip().strip('"')


def parse_pds3_lbl_start_stop_utc(lbl_text: str) -> Optional[Tuple[datetime, datetime]]:
    """
    Parse PDS3 START_TIME / STOP_TIME from a product label (day-of-year form).

    Example: ``2015-041T03:10:30`` → 2015, day-of-year 41, UTC.
    """
    st = pds3_lbl_value_line(lbl_text, "START_TIME")
    et = pds3_lbl_value_line(lbl_text, "STOP_TIME")
    if not st or not et:
        return None

    def _parse_token(token: str) -> datetime:
        m = re.match(r"(\d{4})-(\d{3})T(\d{2}):(\d{2}):(\d{2})", token.strip())
        if not m:
            raise ValueError(f"Unrecognized PDS3 UTC time token: {token!r}")
        y, doy, hh, mm, ss = map(int, m.groups())
        return datetime.strptime(f"{y}{doy:03d}", "%Y%j").replace(
            hour=hh, minute=mm, second=ss, tzinfo=timezone.utc
        )

    t0 = _parse_token(st)
    t1 = _parse_token(et)
    if t1 < t0:
        t0, t1 = t1, t0
    return t0, t1


def interval_overlaps_window(
    start: datetime,
    stop: datetime,
    window_start: datetime,
    window_end: datetime,
) -> bool:
    if start.tzinfo is None:
        start = start.replace(tzinfo=timezone.utc)
    if stop.tzinfo is None:
        stop = stop.replace(tzinfo=timezone.utc)
    if window_start.tzinfo is None:
        window_start = window_start.replace(tzinfo=timezone.utc)
    if window_end.tzinfo is None:
        window_end = window_end.replace(tzinfo=timezone.utc)
    return start <= window_end and stop >= window_start
