"""Discover locally cached DSN tracking products under data/raw/dsn_tracking."""

from __future__ import annotations

import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Optional

from scripts.utils.dsn_pds_ingest import _parse_pds_timestamp


def _label_path_for_data_file(data_path: Path) -> Optional[Path]:
    for suffix in (".lbl", ".LBL"):
        candidate = data_path.with_suffix(suffix)
        if candidate.is_file():
            return candidate
    return None


def _read_label_field(label_path: Path, field_name: str) -> Optional[str]:
    pattern = re.compile(rf"^{re.escape(field_name)}\s*=\s*(.+)$", re.IGNORECASE)
    for line in label_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = pattern.match(line.strip())
        if match:
            return match.group(1).strip().strip('"')
    return None


def _label_time_window(label_path: Path) -> tuple[datetime, datetime]:
    start_value = _read_label_field(label_path, "START_TIME")
    stop_value = _read_label_field(label_path, "STOP_TIME")
    if not start_value or not stop_value:
        raise ValueError(f"{label_path.name} is missing START_TIME/STOP_TIME")
    return _parse_pds_timestamp(start_value), _parse_pds_timestamp(stop_value)


def label_covers_perigee(
    label_path: Path,
    perigee: datetime,
    window_hours: float = 48.0,
) -> bool:
    start, stop = _label_time_window(label_path)
    if stop < start:
        start, stop = stop, start
    window_start = perigee - timedelta(hours=window_hours / 2.0)
    window_end = perigee + timedelta(hours=window_hours / 2.0)
    if perigee.tzinfo is None:
        perigee = perigee.replace(tzinfo=timezone.utc)
    if window_start.tzinfo is None:
        window_start = window_start.replace(tzinfo=timezone.utc)
    if window_end.tzinfo is None:
        window_end = window_end.replace(tzinfo=timezone.utc)
    return start <= window_end and stop >= window_start


def is_trk234_archive(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(4) == b"NJPL"


def is_label_only_dat_file(data_path: Path) -> bool:
    if not data_path.is_file():
        return False
    if data_path.suffix.lower() != ".dat":
        return False
    header = data_path.read_bytes()[:64]
    return header.startswith(b"CCSD") or b"PDS_VERSION_ID" in header[:512]


def discover_dsn_tracking_file(
    mission_dir: Path,
    perigee: Optional[datetime] = None,
    window_hours: float = 48.0,
) -> Optional[Path]:
    if not mission_dir.is_dir():
        return None

    preferred = mission_dir / f"{mission_dir.name}_doppler.trk"
    if preferred.is_file():
        return preferred

    candidates: list[Path] = []
    for pattern in ("*.dat", "*.DAT", "*.trk", "*.TRK", "*.odf", "*.ODF", "*.RAW", "*.raw"):
        candidates.extend(mission_dir.rglob(pattern))

    valid: list[Path] = []
    for candidate in candidates:
        if is_label_only_dat_file(candidate):
            continue
        label_path = _label_path_for_data_file(candidate)
        if perigee is not None and label_path is not None:
            try:
                if not label_covers_perigee(label_path, perigee, window_hours=window_hours):
                    continue
            except ValueError:
                continue
        valid.append(candidate)

    if not valid:
        return None

    trk234_candidates = [path for path in valid if is_trk234_archive(path)]
    if trk234_candidates:
        return max(trk234_candidates, key=lambda path: path.stat().st_size)

    return max(valid, key=lambda path: path.stat().st_size)


def discover_trk234_file(
    project_root: Path,
    mission: str = "MESSENGER_2005",
    perigee: Optional[datetime] = None,
    window_hours: float = 48.0,
) -> Optional[Path]:
    mission_dir = project_root / "data" / "raw" / "dsn_tracking" / mission
    return discover_dsn_tracking_file(
        mission_dir,
        perigee=perigee,
        window_hours=window_hours,
    )
