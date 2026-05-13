"""Extract TRK-2-18 ODF observables from PDS3 label-guided binary tables."""

from __future__ import annotations

import re
import struct
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

ODF_EPOCH = datetime(1950, 1, 1, tzinfo=timezone.utc)


@dataclass(frozen=True)
class _Column:
    name: str
    start_byte: int
    nbytes: int


@dataclass(frozen=True)
class _Table:
    name: str
    start_record: int
    rows: int
    row_bytes: int
    columns: tuple[_Column, ...]


def _label_path_for_data_file(data_path: Path) -> Path:
    for suffix in (".lbl", ".LBL"):
        candidate = data_path.with_suffix(suffix)
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"No PDS label found for {data_path}")


def _read_label_field(label_text: str, field_name: str) -> str:
    pattern = re.compile(rf"^{re.escape(field_name)}\s*=\s*(.+)$", re.IGNORECASE | re.MULTILINE)
    match = pattern.search(label_text)
    if not match:
        raise ValueError(f"Label field {field_name} not found")
    return match.group(1).strip().strip('"')


def _parse_table_locations(label_text: str) -> dict[str, int]:
    locations: dict[str, int] = {}
    for match in re.finditer(
        r"^\^(\w+)\s*=\s*\(\"[^\"]+\",\s*(\d+)\)",
        label_text,
        flags=re.MULTILINE,
    ):
        locations[match.group(1)] = int(match.group(2))
    return locations


def _parse_table_definition(label_text: str, table_name: str) -> _Table:
    block_match = re.search(
        rf"OBJECT\s*=\s*{re.escape(table_name)}\s*(.*?)\nEND_OBJECT\s*=\s*{re.escape(table_name)}",
        label_text,
        flags=re.DOTALL,
    )
    if not block_match:
        raise ValueError(f"Table {table_name} not defined in label")

    block = block_match.group(1)
    rows = int(re.search(r"ROWS\s*=\s*(\d+)", block).group(1))
    row_bytes = int(re.search(r"ROW_BYTES\s*=\s*(\d+)", block).group(1))
    columns: list[_Column] = []
    for column_match in re.finditer(
        r"OBJECT\s*=\s*COLUMN\s*(.*?)\n\s*END_OBJECT\s*=\s*COLUMN",
        block,
        flags=re.DOTALL,
    ):
        column_block = column_match.group(1)
        name_match = re.search(r'NAME\s*=\s*"([^"]+)"', column_block)
        start_match = re.search(r"START_BYTE\s*=\s*(\d+)", column_block)
        bytes_match = re.search(r"BYTES\s*=\s*(\d+)", column_block)
        if not name_match or not start_match or not bytes_match:
            continue
        columns.append(
            _Column(
                name=name_match.group(1),
                start_byte=int(start_match.group(1)),
                nbytes=int(bytes_match.group(1)),
            )
        )

    if not columns:
        raise ValueError(f"Table {table_name} has no readable columns")

    locations = _parse_table_locations(label_text)
    if table_name not in locations:
        raise ValueError(f"Missing ^{table_name} location in label")

    return _Table(
        name=table_name,
        start_record=locations[table_name],
        rows=rows,
        row_bytes=row_bytes,
        columns=tuple(columns),
    )


def _odf_timestamp(seconds: int, nanoseconds: int) -> datetime:
    return ODF_EPOCH + timedelta(seconds=seconds, microseconds=nanoseconds / 1000.0)


def _read_column(row: bytes, column: _Column) -> int:
    start = column.start_byte - 1
    end = start + column.nbytes
    if column.nbytes == 4:
        return struct.unpack(">I", row[start:end])[0]
    if column.nbytes == 8:
        return struct.unpack(">Q", row[start:end])[0]
    raise ValueError(f"Unsupported column width for {column.name}: {column.nbytes}")


def _read_table_rows(
    data_path: Path,
    record_bytes: int,
    table: _Table,
) -> list[bytes]:
    offset = (table.start_record - 1) * record_bytes
    row_span = table.rows * table.row_bytes
    with data_path.open("rb") as handle:
        handle.seek(offset)
        payload = handle.read(row_span)
    if len(payload) < row_span:
        raise ValueError(
            f"Table {table.name} requires {row_span} bytes at offset {offset}, found {len(payload)}"
        )
    return [
        payload[index : index + table.row_bytes]
        for index in range(0, row_span, table.row_bytes)
    ]


def _column_value(row: bytes, columns: tuple[_Column, ...], target_name: str) -> int:
    for column in columns:
        if column.name == target_name:
            return _read_column(row, column)
    raise KeyError(target_name)


def extract_trk218_measurements(filepath: Path) -> list[dict[str, Any]]:
    """Decode TRK-2-18 ODF summary records from a PDS3 label-guided binary file."""
    label_path = _label_path_for_data_file(filepath)
    label_text = label_path.read_text(encoding="utf-8", errors="ignore")
    record_bytes = int(_read_label_field(label_text, "RECORD_BYTES"))
    summary_table = _parse_table_definition(label_text, "ODF7B_TABLE")
    rows = _read_table_rows(filepath, record_bytes, summary_table)

    measurements: list[dict[str, Any]] = []
    for row_index, row in enumerate(rows):
        first_seconds = _column_value(row, summary_table.columns, "FIRST SAMPLE TIME - INTEGER PART")
        first_fraction = _column_value(
            row,
            summary_table.columns,
            "FIRST SAMPLE TIME - FRACTIONAL PART",
        )
        last_seconds = _column_value(row, summary_table.columns, "LAST SAMPLE TIME - INTEGER PART")
        last_fraction = _column_value(
            row,
            summary_table.columns,
            "LAST SAMPLE TIME - FRACTIONAL PART",
        )
        station_id = _column_value(row, summary_table.columns, "STATION ID")
        sample_count = _column_value(row, summary_table.columns, "NUMBER OF SAMPLES")
        doppler_id = _column_value(row, summary_table.columns, "NETWORK OR DOPPLER ID")
        band_id = _column_value(row, summary_table.columns, "BAND ID")
        data_type_id = _column_value(row, summary_table.columns, "DATA TYPE ID")

        first_timestamp = _odf_timestamp(first_seconds, first_fraction)
        last_timestamp = _odf_timestamp(last_seconds, last_fraction)
        measurements.append(
            {
                "record_index": row_index,
                "timestamp": first_timestamp.isoformat(),
                "last_timestamp": last_timestamp.isoformat(),
                "station_id": station_id,
                "doppler_channel_id": doppler_id,
                "band_id": band_id,
                "data_type_id": data_type_id,
                "sample_count": sample_count,
                "data_source": "TRK-2-18",
                "source_file": filepath.name,
            }
        )

    return measurements
