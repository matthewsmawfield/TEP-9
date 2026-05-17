"""
Optional time series of flyby-related measurements for Step 042.

Real data only: JSON files under data/time_resolved_flyby_residuals/
See data/time_resolved_flyby_residuals/schema.json for the contract.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REQUIRED_TOP_LEVEL = {"mission", "quantity", "quantity_description", "reference", "utc_iso", "value_mm_s"}

JUNO_042_REFERENCE = "NASA PDS TRK ingest; Step 030 Juno 2013 reanalysis"


def build_juno_pairwise_042_sidecar_from_rows(
    residuals: List[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """
    Build the Step 042 sidecar payload from Step 030 pairwise proxy rows.

    Each row must carry ``timestamp`` and ``velocity_mm_s`` from real TRK processing.
    Returns None if fewer than three valid samples (caller should skip export).
    """
    utc_iso: List[str] = []
    value_mm_s: List[float] = []
    for row in residuals:
        ts_raw = row.get("timestamp")
        if ts_raw is None or "velocity_mm_s" not in row:
            continue
        try:
            ts = datetime.fromisoformat(str(ts_raw).replace("Z", "+00:00"))
            if ts.tzinfo is None:
                ts = ts.replace(tzinfo=timezone.utc)
            else:
                ts = ts.astimezone(timezone.utc)
        except (ValueError, TypeError):
            continue
        utc_iso.append(ts.strftime("%Y-%m-%dT%H:%M:%SZ"))
        value_mm_s.append(float(row["velocity_mm_s"]))
    if len(utc_iso) < 3:
        return None
    return {
        "mission": "Juno",
        "quantity": "pairwise_doppler_proxy_mm_s",
        "quantity_description": (
            "Sequential pairwise ΔDoppler per station converted with Δv ≈ c Δf/f at the "
            "receiver frequency recorded in each TRK measurement; not batch post-fit OD "
            "residual (see Step 030 documentation)."
        ),
        "reference": JUNO_042_REFERENCE,
        "utc_iso": utc_iso,
        "value_mm_s": value_mm_s,
        "sigma_mm_s": None,
    }


def write_juno_042_sidecar(project_root: Path, payload: Dict[str, Any]) -> Path:
    """Write validated Juno sidecar to data/time_resolved_flyby_residuals/Juno.json."""
    out_dir = residuals_dir(project_root)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "Juno.json"
    missing = REQUIRED_TOP_LEVEL - set(payload.keys())
    if missing:
        raise ValueError(f"Juno 042 payload missing keys: {sorted(missing)}")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    return out_path


def save_juno_pairwise_step030_archive(project_root: Path, rows: List[Dict[str, Any]]) -> Path:
    """
    Persist full pairwise proxy rows under results/ for reproducibility and for
    materializing ``data/time_resolved_flyby_residuals/Juno.json`` without re-parsing TRK.
    """
    out_path = project_root / "results" / "step030_juno_pairwise_residual_series.json"
    payload = {
        "schema_version": 1,
        "source_step": "030",
        "mission": "Juno",
        "proxy_kind": "doppler_pair_mm_s",
        "description": (
            "Pairwise Doppler proxy rows from real TRK-2-25 ingestion (Step 030)."
        ),
        "n_rows": len(rows),
        "rows": rows,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    return out_path


def save_juno_ramp_pairwise_step030_archive(project_root: Path, rows: List[Dict[str, Any]]) -> Path:
    """
    Persist per-station sequential Δ ``ramp_freq_hz`` rows (Hz, not mm/s).

    Step 042 ``Juno.json`` remains Doppler-mm/s only; this archive is for audit,
    replays, and any future Hz-aware consumer.
    """
    out_path = project_root / "results" / "step030_juno_pairwise_ramp_freq_series.json"
    payload = {
        "schema_version": 1,
        "source_step": "030",
        "mission": "Juno",
        "proxy_kind": "ramp_freq_pair_delta_hz",
        "unit": "Hz",
        "description": (
            "Pairwise differences of archival ramp_freq_hz per station from real "
            "TRK-2-34 decode (Step 030). Not Doppler-derived mm/s; not for the 042 "
            "sidecar schema (value_mm_s)."
        ),
        "n_rows": len(rows),
        "rows": rows,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    return out_path


def residuals_dir(project_root: Path) -> Path:
    return project_root / "data" / "time_resolved_flyby_residuals"


def load_validated_residual_file(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(str(path))
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"{path}: root must be an object")
    missing = REQUIRED_TOP_LEVEL - set(data.keys())
    if missing:
        raise ValueError(f"{path}: missing required keys: {sorted(missing)}")
    utc = data["utc_iso"]
    val = data["value_mm_s"]
    if not isinstance(utc, list) or not isinstance(val, list):
        raise ValueError(f"{path}: utc_iso and value_mm_s must be arrays")
    if len(utc) != len(val):
        raise ValueError(f"{path}: utc_iso length {len(utc)} != value_mm_s length {len(val)}")
    if len(utc) < 3:
        raise ValueError(f"{path}: need at least 3 samples, got {len(utc)}")
    sig = data.get("sigma_mm_s")
    if sig is not None:
        if not isinstance(sig, list) or len(sig) != len(val):
            raise ValueError(f"{path}: sigma_mm_s must be null or same length as value_mm_s")
    return data


def _parse_utc_iso(s: str) -> datetime:
    if not isinstance(s, str):
        raise TypeError("utc_iso entries must be strings")
    if s.endswith("Z"):
        return datetime.fromisoformat(s.replace("Z", "+00:00"))
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)


def build_residual_timeline(data: Dict[str, Any]) -> List[Tuple[datetime, float, Optional[float]]]:
    """Return list of (t_utc, value_mm_s, sigma_or_none)."""
    out: List[Tuple[datetime, float, Optional[float]]] = []
    sigs = data.get("sigma_mm_s")
    for i, ts in enumerate(data["utc_iso"]):
        t = _parse_utc_iso(str(ts))
        v = float(data["value_mm_s"][i])
        sig: Optional[float] = None
        if isinstance(sigs, list) and i < len(sigs) and sigs[i] is not None:
            sig = float(sigs[i])
        out.append((t, v, sig))
    out.sort(key=lambda x: x[0])
    return out


def align_geometry_epochs_to_residuals(
    epoch_utc_iso: List[str],
    residual_timeline: List[Tuple[datetime, float, Optional[float]]],
    max_delta_s: float,
) -> Tuple[List[int], List[float], List[float], int]:
    """
    For each geometry epoch (in order), take nearest residual by absolute time if within max_delta_s.

    Returns:
        epoch_indices: indices into epoch table that received a match
        y_residual_mm_s: matched residual values
        delta_t_s: signed (residual_time - geometry_time) for the match
        n_skipped_no_neighbor: count of epochs with no sample within window
    """
    if not residual_timeline:
        return [], [], [], len(epoch_utc_iso)

    res_times = [t for t, _, _ in residual_timeline]
    res_vals = [v for _, v, _ in residual_timeline]

    epoch_indices: List[int] = []
    y_residual: List[float] = []
    delta_t: List[float] = []
    skipped = 0

    for i, giso in enumerate(epoch_utc_iso):
        try:
            tg = _parse_utc_iso(str(giso))
        except (ValueError, TypeError):
            skipped += 1
            continue
        best_j = -1
        best_dt = float("inf")
        for j, tr in enumerate(res_times):
            d = abs((tr - tg).total_seconds())
            if d < best_dt:
                best_dt = d
                best_j = j
        if best_j < 0 or best_dt > max_delta_s:
            skipped += 1
            continue
        epoch_indices.append(i)
        y_residual.append(float(res_vals[best_j]))
        delta_t.append(float((res_times[best_j] - tg).total_seconds()))

    return epoch_indices, y_residual, delta_t, skipped
