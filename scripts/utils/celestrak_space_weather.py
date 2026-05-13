"""Celestrak SW-All daily space weather lookup (F10.7, Kp)."""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Dict

import requests

CELESTRAK_SW_ALL_URL = "https://celestrak.org/SpaceData/SW-All.csv"
WAYBACK_SW_ALL_URL = (
    "https://web.archive.org/web/20250101000000/https://celestrak.org/SpaceData/SW-All.csv"
)
CACHE_PATH = (
    Path(__file__).resolve().parent.parent.parent
    / "data"
    / "raw"
    / "space_weather"
    / "SW-All.csv"
)
_SW_ALL_BY_DATE: Dict[str, Dict[str, str]] | None = None
_SW_ALL_SOURCE: str | None = None


def _request_sw_all_text(url: str) -> str:
    response = requests.get(
        url,
        timeout=120,
        headers={
            "User-Agent": "TEP-EFA-Pipeline/1.0 (Academic Research)",
            "Accept": "text/csv,*/*",
        },
    )
    response.raise_for_status()
    if "DATE,BSRN" not in response.text[:200]:
        raise RuntimeError(f"Unexpected SW-All payload from {url}")
    return response.text


def _write_cache(text: str) -> None:
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    CACHE_PATH.write_text(text, encoding="utf-8")


def _load_sw_all_text() -> tuple[str, str]:
    if CACHE_PATH.is_file():
        return CACHE_PATH.read_text(encoding="utf-8"), str(CACHE_PATH)

    errors: list[str] = []
    for source_name, url in (
        ("Celestrak_SW-All.csv", CELESTRAK_SW_ALL_URL),
        ("Internet_Archive_SW-All.csv", WAYBACK_SW_ALL_URL),
    ):
        try:
            text = _request_sw_all_text(url)
            _write_cache(text)
            return text, source_name
        except requests.RequestException as exc:
            errors.append(f"{source_name}: {exc}")

    raise RuntimeError(
        "Unable to load Celestrak SW-All.csv from live Celestrak or Internet Archive; "
        + "; ".join(errors)
    )


def _load_sw_all_by_date() -> Dict[str, Dict[str, str]]:
    global _SW_ALL_BY_DATE, _SW_ALL_SOURCE
    if _SW_ALL_BY_DATE is not None:
        return _SW_ALL_BY_DATE

    text, source = _load_sw_all_text()
    by_date: Dict[str, Dict[str, str]] = {}
    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        date = row.get("DATE", "").strip()
        if date:
            by_date[date] = row

    if not by_date:
        raise RuntimeError("Celestrak SW-All.csv returned no rows")

    _SW_ALL_BY_DATE = by_date
    _SW_ALL_SOURCE = source
    return by_date


def lookup_space_weather(date_str: str) -> Dict[str, float | str]:
    """Return F10.7 and Kp for a calendar date (YYYY-MM-DD)."""
    from datetime import datetime

    datetime.strptime(date_str, "%Y-%m-%d")
    row = _load_sw_all_by_date().get(date_str)
    if row is None:
        raise RuntimeError(f"No Celestrak SW-All entry for {date_str}")

    f10_7_obs = float(row["F10.7_OBS"])
    if f10_7_obs < 0:
        f10_7 = float(row["F10.7_ADJ"])
        f10_7_field = "F10.7_ADJ"
    else:
        f10_7 = f10_7_obs
        f10_7_field = "F10.7_OBS"

    if f10_7 < 0:
        raise RuntimeError(f"No usable F10.7 value in Celestrak SW-All for {date_str}")

    kp_raw = row.get("KP_SUM", "").strip()
    if not kp_raw:
        raise RuntimeError(f"No usable KP_SUM value in Celestrak SW-All for {date_str}")
    kp_sum = float(kp_raw)
    if kp_sum < 0:
        raise RuntimeError(f"No usable KP_SUM value in Celestrak SW-All for {date_str}")

    source = _SW_ALL_SOURCE or "Celestrak_SW-All.csv"
    return {
        "f10_7": f10_7,
        "kp_sum": kp_sum,
        "kp": kp_sum / 10.0,
        "f10_7_field": f10_7_field,
        "data_source": source,
        "source_url": CELESTRAK_SW_ALL_URL,
    }
