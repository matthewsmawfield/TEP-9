"""Extract radiometric observables from TRK-2-34 SFDU archives."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List


def extract_trk234_measurements(filepath: Path) -> List[Dict[str, Any]]:
    """Decode a TRK-2-34 archive and return per-SFDU observables."""
    import trk234

    reader = trk234.Reader(str(filepath))
    if len(reader.index) < 2:
        return []

    reader.decode()
    measurements: List[Dict[str, Any]] = []

    for index, sfdu in enumerate(reader.sfdu_list):
        trk = getattr(sfdu, "trk_chdo", None)
        if trk is None:
            continue

        ramp_freq = getattr(trk, "ramp_freq", None)
        ul_lo_phs_cycles = getattr(trk, "ul_lo_phs_cycles", None)
        ul_hi_phs_cycles = getattr(trk, "ul_hi_phs_cycles", None)
        ul_frac_phs_cycles = getattr(trk, "ul_frac_phs_cycles", None)

        if ramp_freq is None and ul_lo_phs_cycles is None:
            continue

        uplink_dss, downlink_dss = sfdu.dss_id()
        uplink_band, downlink_band = sfdu.radio_band()

        measurements.append(
            {
                "sfdu_index": index,
                "timestamp": sfdu.timestamp().isoformat(),
                "station_uplink": uplink_dss,
                "station_downlink": downlink_dss,
                "band_uplink": uplink_band,
                "band_downlink": downlink_band,
                "tracking_mode": sfdu.tracking_mode(),
                "ramp_freq_hz": float(ramp_freq) if ramp_freq is not None else None,
                "ul_lo_phs_cycles": int(ul_lo_phs_cycles) if ul_lo_phs_cycles is not None else None,
                "ul_hi_phs_cycles": int(ul_hi_phs_cycles) if ul_hi_phs_cycles is not None else None,
                "ul_frac_phs_cycles": int(ul_frac_phs_cycles) if ul_frac_phs_cycles is not None else None,
                "data_source": "TRK-2-34",
                "source_file": filepath.name,
            }
        )

    return measurements
