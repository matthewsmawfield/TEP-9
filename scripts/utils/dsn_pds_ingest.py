"""Attempt perigee-matched NASA PDS ingest for DSN tracking products."""

from __future__ import annotations

import json
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Optional
from urllib.parse import quote, urlparse

import requests

PDS_SEARCH_API = "https://pds.mcp.nasa.gov/api/search/1/products"
TRACKING_EXTENSIONS = (".dat", ".DAT", ".trk", ".TRK", ".tnf", ".TNF", ".odf", ".ODF")
PREDICTION_MISSION_ALIASES = {
    "NEAR_1998": "NEAR",
    "Juno_2013": "Juno",
    "Cassini_1999": "Cassini",
    "Stardust_2001": "Stardust",
}
MISSION_DIR_ALIASES = {
    "NEAR": "NEAR_1998",
    "Juno": "Juno_2013",
    "Cassini": "Cassini_1999",
    "Stardust": "Stardust_2001",
}
TRACKING_KEYWORDS = (
    "trk",
    "tnf",
    "doppler",
    "tracking",
    "radio science",
    "rss",
    "ddor",
    "odf",
    "atdf",
)
REFERENCE_PRODUCT_LIDS = {
    "DART_TRK234_REFERENCE": [
        "urn:nasa:pds:dart:data_trk234:dart_hga_tnf_20211124t032621_v01::1.0",
    ],
}
# PDS MCP search rejects the `start` offset parameter; results cap at `limit` (default 100)
# with no `next` link when total hits exceed the page size. Subdivide `lid like` queries
# until each leaf returns a complete page (hits == len(data) < limit) or hits == 0.
PDS_MCP_RESULT_CAP = 100
# Only run full LID-tree enumeration when the root query is small enough to stay bounded.
PDS_EXHAUSTIVE_LID_MAX_ROOT_HITS = 400
MISSION_COLLECTION_LIDS: dict[str, list[str]] = {
    "NEAR_1998": [
        "urn:nasa:pds:near_rss_raw:data_odf::1.0",
    ],
    "MESSENGER_2005": [
        "urn:nasa:pds:mess-rs-raw:data-tnf::1.0",
    ],
    "Galileo_1990": [
        "urn:nasa:pds:go-j-rss-1-edr::1.0",
        "urn:nasa:pds:go-j-rss-2-edr::1.0",
    ],
    "Galileo_1992": [
        "urn:nasa:pds:go-j-rss-1-edr::1.0",
        "urn:nasa:pds:go-j-rss-2-edr::1.0",
    ],
    "Juno_2013": [
        "urn:nasa:pds:jno-e-rss-1-edr::1.0",
        "urn:nasa:pds:jno-j-rss-1-edr::1.0",
    ],
    "Cassini_1999": [
        "urn:nasa:pds:co-rss-1-edr::1.0",
        "urn:nasa:pds:co-rss-2-edr::1.0",
    ],
}


def _parse_pds_timestamp(value: str) -> datetime:
    cleaned = value.strip().strip('"')
    if cleaned.upper() == "N/A":
        raise ValueError("missing timestamp")
    # PDS3 day-of-year form: 2015-041T03:10:30 (must precede ISO ``fromisoformat``)
    match = re.match(r"^(\d{4})-(\d{3})T", cleaned)
    if match:
        year = int(match.group(1))
        doy = int(match.group(2))
        remainder = cleaned.split("T", 1)[1]
        base = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1)
        hour, minute, second = remainder.split(":")
        return base.replace(
            hour=int(hour),
            minute=int(minute),
            second=int(float(second)),
        )
    if cleaned.endswith("Z"):
        return datetime.fromisoformat(cleaned.replace("Z", "+00:00"))
    if "T" in cleaned:
        return datetime.fromisoformat(cleaned).replace(tzinfo=timezone.utc)
    return datetime.fromisoformat(cleaned).replace(tzinfo=timezone.utc)


def resolve_prediction_key(mission: str) -> str:
    return PREDICTION_MISSION_ALIASES.get(mission, mission)


def resolve_tracking_mission_dir(mission_name: str) -> str:
    return MISSION_DIR_ALIASES.get(mission_name, mission_name)


def flyby_tracking_missions(project_root: Path) -> list[str]:
    catalog_file = project_root / "results" / "step003_archival_flyby_catalog.json"
    if not catalog_file.exists():
        raise FileNotFoundError(
            "Step 003 flyby catalog is required for DSN tracking mission selection"
        )
    with open(catalog_file, encoding="utf-8") as handle:
        catalog = json.load(handle)

    missions: list[str] = []
    seen: set[str] = set()
    for flyby in catalog.get("flybys", []):
        if not flyby.get("dsn_tracking_available"):
            continue
        mission_dir = resolve_tracking_mission_dir(flyby["mission_name"])
        mission_path = project_root / "data" / "raw" / "dsn_tracking" / mission_dir
        if not mission_path.is_dir():
            continue
        if mission_dir in seen:
            continue
        seen.add(mission_dir)
        missions.append(mission_dir)
    return missions


def _relative_path(project_root: Path, path: Path) -> str:
    try:
        return str(path.relative_to(project_root))
    except ValueError:
        return str(path)


def _perigee_title_tokens(
    perigee: datetime,
    window_start: datetime,
    window_end: datetime,
) -> list[str]:
    tokens = {
        window_start.strftime("%Y-%m-%d"),
        window_end.strftime("%Y-%m-%d"),
        perigee.strftime("%Y-%m-%d"),
        perigee.strftime("%Y-%j"),
        window_start.strftime("%Y-%j"),
        window_end.strftime("%Y-%j"),
    }
    return sorted(tokens)


def load_perigee_datetime(project_root: Path, mission: str) -> datetime:
    trajectory_file = project_root / "data" / "raw" / "flyby_trajectories" / f"{mission}_trajectory.json"
    if trajectory_file.exists():
        with open(trajectory_file, encoding="utf-8") as handle:
            perigee = json.load(handle).get("perigee", {})
        datetime_str = str(perigee.get("datetime", "")).replace("A.D. ", "").strip()
        if not datetime_str:
            raise KeyError(f"{_relative_path(project_root, trajectory_file)} missing perigee.datetime")
        for fmt in ("%Y-%b-%d %H:%M:%S.%f", "%Y-%b-%d %H:%M:%S", "%Y-%m-%d"):
            try:
                parsed = datetime.strptime(datetime_str, fmt)
                return parsed.replace(tzinfo=timezone.utc)
            except ValueError:
                continue
        raise ValueError(f"Could not parse perigee datetime for {mission}: {datetime_str}")

    horizons_trajectory = (
        project_root
        / "data"
        / "raw"
        / "jpl_horizons"
        / mission
        / f"{mission}_trajectory.json"
    )
    if horizons_trajectory.exists():
        with open(horizons_trajectory, encoding="utf-8") as handle:
            data = json.load(handle)
        timestamps = data.get("timestamp")
        ranges = data.get("range_m")
        if not isinstance(timestamps, list) or not isinstance(ranges, list) or len(timestamps) != len(ranges):
            raise ValueError(
                f"{_relative_path(project_root, horizons_trajectory)} missing timestamp/range_m arrays"
            )
        try:
            min_index, _ = min(enumerate(ranges), key=lambda item: float(item[1]))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{_relative_path(project_root, horizons_trajectory)} contains non-numeric range_m entries"
            ) from exc
        ts = timestamps[min_index]
        if not isinstance(ts, str) or not ts.strip():
            raise ValueError(
                f"{_relative_path(project_root, horizons_trajectory)} contains invalid timestamp at perigee index"
            )
        for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
            try:
                parsed = datetime.strptime(ts.strip(), fmt)
                return parsed.replace(tzinfo=timezone.utc)
            except ValueError:
                continue
        raise ValueError(f"Could not parse Horizons perigee timestamp for {mission}: {ts}")

    predictions_file = project_root / "results" / "step007_tep_predictions.json"
    if not predictions_file.exists():
        raise FileNotFoundError(
            "Perigee-matched DSN ingest requires perigee times from "
            f"{_relative_path(project_root, trajectory_file)} (Step 004 output) "
            "or results/step007_tep_predictions.json."
        )
    with open(predictions_file, encoding="utf-8") as handle:
        predictions = json.load(handle).get("predictions", {})
    prediction_key = resolve_prediction_key(mission)
    if prediction_key not in predictions:
        raise KeyError(f"Mission {mission} not found in Step 007 predictions")

    datetime_str = predictions[prediction_key]["perigee"]["datetime"].replace("A.D. ", "").strip()
    for fmt in ("%Y-%b-%d %H:%M:%S.%f", "%Y-%b-%d %H:%M:%S", "%Y-%m-%d"):
        try:
            parsed = datetime.strptime(datetime_str, fmt)
            return parsed.replace(tzinfo=timezone.utc)
        except ValueError:
            continue
    raise ValueError(f"Could not parse perigee datetime for {mission}: {datetime_str}")


def load_mission_metadata(project_root: Path, mission: str) -> dict[str, Any]:
    metadata_file = project_root / "data" / "raw" / "dsn_tracking" / mission / "metadata.json"
    if not metadata_file.exists():
        return {}
    with open(metadata_file, encoding="utf-8") as handle:
        return json.load(handle)


def _normalize_collection_lid(value: str) -> str:
    cleaned = value.strip().rstrip("/")
    if cleaned.startswith("urn:"):
        return cleaned if "::" in cleaned else f"{cleaned}::1.0"
    return cleaned


def mission_collection_lids(mission: str, metadata: dict[str, Any]) -> list[str]:
    lids: list[str] = []
    seen: set[str] = set()
    for lid in MISSION_COLLECTION_LIDS.get(mission, []):
        normalized = _normalize_collection_lid(lid)
        if normalized not in seen:
            seen.add(normalized)
            lids.append(normalized)
    archive_info = metadata.get("archive_info", {})
    for entry in archive_info.get("pds_collections", []):
        if not isinstance(entry, str):
            continue
        if entry.startswith("urn:"):
            normalized = _normalize_collection_lid(entry)
            if normalized not in seen:
                seen.add(normalized)
                lids.append(normalized)
    return lids


def metadata_archive_search_tokens(metadata: dict[str, Any]) -> list[str]:
    tokens: list[str] = []
    seen: set[str] = set()
    archive_info = metadata.get("archive_info", {})
    for entry in archive_info.get("pds_collections", []):
        if not isinstance(entry, str) or entry.startswith("urn:"):
            continue
        slug = entry.rstrip("/").split("/")[-1].lower()
        if not slug or slug in seen:
            continue
        seen.add(slug)
        tokens.append(slug.replace(".", "-"))
        if "-edr-" in slug or "-rss-" in slug:
            tokens.append(slug.split("-edr", 1)[0])
    return tokens


def _product_time_bounds(product: dict[str, Any]) -> tuple[datetime, datetime]:
    properties = product.get("properties", {})
    start_value = (
        product.get("start_date_time")
        or properties.get("pds:Time_Coordinates.pds:start_date_time", [None])[0]
    )
    stop_value = (
        product.get("stop_date_time")
        or properties.get("pds:Time_Coordinates.pds:stop_date_time", [None])[0]
    )
    if not isinstance(start_value, str) or not isinstance(stop_value, str):
        raise ValueError("product missing PDS time coordinates")
    start = _parse_pds_timestamp(start_value)
    stop = _parse_pds_timestamp(stop_value)
    if start.tzinfo is None:
        start = start.replace(tzinfo=timezone.utc)
    if stop.tzinfo is None:
        stop = stop.replace(tzinfo=timezone.utc)
    return start, stop


def _product_overlaps_window(
    product: dict[str, Any],
    window_start: datetime,
    window_end: datetime,
) -> bool:
    try:
        start, stop = _product_time_bounds(product)
    except ValueError:
        return False
    return start <= window_end and stop >= window_start


def _mission_search_tokens(mission: str) -> list[str]:
    base = mission.split("_", 1)[0]
    tokens = {base, mission.replace("_", " ")}
    if base.lower() == "near":
        tokens.add("NEAR")
    if base.lower() == "galileo":
        tokens.add("Galileo")
    if base.lower() == "cassini":
        tokens.add("Cassini")
    if base.lower() == "messenger":
        tokens.add("MESSENGER")
    if base.lower() == "juno":
        tokens.add("Juno")
    if base.lower() == "stardust":
        tokens.add("Stardust")
    return sorted(tokens)


def _build_time_window_query(
    start: datetime,
    end: datetime,
    title_token: Optional[str] = None,
) -> str:
    start_iso = start.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")
    end_iso = end.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")
    clauses = [
        f'(pds:Time_Coordinates.pds:start_date_time le "{end_iso}")',
        f'(pds:Time_Coordinates.pds:stop_date_time ge "{start_iso}")',
    ]
    if title_token:
        clauses.append(
            f'(pds:Identification_Area.pds:title like "*{title_token}*")'
        )
    return "(" + " and ".join(clauses) + ")"


def _product_download_urls(product: dict[str, Any]) -> list[str]:
    properties = product.get("properties", {})
    refs = properties.get("ops:Label_File_Info.ops:file_ref", [])
    names = properties.get("ops:Label_File_Info.ops:file_name", [])
    file_names = properties.get("pds:File.pds:file_name", [])
    urls: list[str] = []
    for ref in refs:
        if isinstance(ref, str) and ref.startswith("http"):
            urls.append(ref)
    label_url = product.get("metadata", {}).get("label_url")
    for name in list(names) + list(file_names):
        if not isinstance(name, str):
            continue
        if name.lower().endswith(TRACKING_EXTENSIONS):
            if label_url:
                parsed = urlparse(label_url)
                urls.append(f"{parsed.scheme}://{parsed.netloc}{Path(parsed.path).with_name(name)}")
    for url in re.findall(r"https?://[^\"\\s<>]+", json.dumps(product)):
        if any(url.lower().endswith(ext.lower()) for ext in TRACKING_EXTENSIONS):
            urls.append(url)
    return sorted(set(urls))


def _is_collection_product_lid(lid: str) -> bool:
    base = lid.split("::", 1)[0]
    return base.endswith(":data_odf") or base.endswith(":data_tnf")


def _is_collection_product(product: dict[str, Any]) -> bool:
    return _is_collection_product_lid(str(product.get("id", "")))


def _looks_like_tracking_product(product: dict[str, Any]) -> bool:
    title = str(product.get("title", "")).lower()
    properties = product.get("properties", {})
    names = properties.get("ops:Label_File_Info.ops:file_name", [])
    name_blob = " ".join(str(name).lower() for name in names)
    blob = f"{title} {name_blob}"
    return any(keyword in blob for keyword in TRACKING_KEYWORDS)


def _rewrite_download_url(url: str) -> str:
    if url.startswith("https://pds.nasa.gov/"):
        return url.replace("https://pds.nasa.gov/", "https://pds.mcp.nasa.gov/", 1)
    if url.startswith("/"):
        return f"https://pds.mcp.nasa.gov{url}"
    return url


def _download_file(
    session: requests.Session,
    url: str,
    destination: Path,
    project_root: Optional[Path] = None,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    response = session.get(url, timeout=120, stream=True)
    response.raise_for_status()
    with open(destination, "wb") as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)
    path_value = (
        _relative_path(project_root, destination)
        if project_root is not None
        else str(destination)
    )
    return {
        "url": url,
        "path": path_value,
        "bytes": destination.stat().st_size,
    }


def fetch_product(session: requests.Session, product_lid: str) -> dict[str, Any]:
    response = session.get(
        f"{PDS_SEARCH_API}/{quote(product_lid, safe='')}",
        timeout=60,
    )
    response.raise_for_status()
    return response.json()


def download_product_files(
    session: requests.Session,
    product: dict[str, Any],
    mission_dir: Path,
    project_root: Optional[Path] = None,
) -> list[dict[str, Any]]:
    downloads: list[dict[str, Any]] = []
    for url in _product_download_urls(product):
        if not any(url.lower().endswith(ext.lower()) for ext in TRACKING_EXTENSIONS):
            continue
        destination = mission_dir / "NASA_PDS" / Path(urlparse(url).path).name
        if destination.exists() and destination.stat().st_size > 0:
            path_value = (
                _relative_path(project_root, destination)
                if project_root is not None
                else str(destination)
            )
            downloads.append(
                {
                    "url": url,
                    "path": path_value,
                    "bytes": destination.stat().st_size,
                    "status": "cached",
                }
            )
            continue
        downloads.append(
            _download_file(
                session,
                _rewrite_download_url(url),
                destination,
                project_root=project_root,
            )
        )
    return downloads


def _collection_is_indexed(session: requests.Session, collection_lid: str) -> bool:
    response = session.get(
        f"{PDS_SEARCH_API}/{quote(collection_lid, safe='')}",
        timeout=60,
    )
    if response.status_code == 404:
        return False
    response.raise_for_status()
    return True


def search_collection_members_by_title(
    session: requests.Session,
    collection_lid: str,
    title_tokens: list[str],
    limit: int = 50,
) -> tuple[list[dict[str, Any]], int]:
    products: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    pages_scanned = 0
    if not _collection_is_indexed(session, collection_lid):
        return products, pages_scanned
    for token in title_tokens:
        response = session.get(
            f"{PDS_SEARCH_API}/{quote(collection_lid, safe='')}/members",
            params={"q": f'(pds:Identification_Area.pds:title like "*{token}*")', "limit": limit},
            timeout=120,
        )
        if response.status_code == 404:
            continue
        response.raise_for_status()
        payload = response.json()
        pages_scanned += 1
        url = payload.get("links", {}).get("next")
        while True:
            for member in payload.get("data", []):
                if not _looks_like_tracking_product(member):
                    continue
                member_id = str(member.get("id", ""))
                if member_id in seen_ids:
                    continue
                seen_ids.add(member_id)
                products.append(fetch_product(session, member_id))
            if not url:
                break
            response = session.get(url, timeout=120)
            response.raise_for_status()
            payload = response.json()
            pages_scanned += 1
            url = payload.get("links", {}).get("next")
    return products, pages_scanned


def search_collection_members(
    session: requests.Session,
    collection_lid: str,
    window_start: datetime,
    window_end: datetime,
    limit: int = 200,
) -> tuple[list[dict[str, Any]], int]:
    products: list[dict[str, Any]] = []
    pages_scanned = 0
    if not _collection_is_indexed(session, collection_lid):
        return products, pages_scanned
    response = session.get(
        f"{PDS_SEARCH_API}/{quote(collection_lid, safe='')}/members",
        params={"limit": limit},
        timeout=120,
    )
    if response.status_code == 404:
        return products, pages_scanned
    response.raise_for_status()
    payload = response.json()
    url = None
    while True:
        pages_scanned += 1
        for member in payload.get("data", []):
            if not _looks_like_tracking_product(member):
                continue
            if not _product_overlaps_window(member, window_start, window_end):
                continue
            member_id = str(member.get("id", ""))
            if not member_id:
                continue
            products.append(fetch_product(session, member_id))
        url = payload.get("links", {}).get("next")
        if not url:
            break
        response = session.get(url, timeout=120)
        response.raise_for_status()
        payload = response.json()
    return products, pages_scanned


def search_lid_family_products(
    session: requests.Session,
    lid_prefix: str,
    window_start: datetime,
    window_end: datetime,
    limit: int = 200,
) -> list[dict[str, Any]]:
    start_iso = window_start.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")
    end_iso = window_end.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")
    query = (
        f'((pds:Time_Coordinates.pds:start_date_time le "{end_iso}") and '
        f'(pds:Time_Coordinates.pds:stop_date_time ge "{start_iso}") and '
        f'(lid like "{lid_prefix}*"))'
    )
    response = session.get(
        PDS_SEARCH_API,
        params={"q": query, "limit": limit},
        timeout=120,
    )
    response.raise_for_status()
    products: list[dict[str, Any]] = []
    for product in response.json().get("data", []):
        if _is_collection_product(product):
            continue
        if not _looks_like_tracking_product(product):
            continue
        products.append(product)
    return products


def search_lid_root_exhaustive(
    session: requests.Session,
    lid_root: str,
    charset: str,
    *,
    max_api_calls: int = 500,
    limit: int = PDS_MCP_RESULT_CAP,
    max_root_hits: Optional[int] = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """
    Enumerate every product under lid_root despite the MCP 100-row cap.

    lid_root must end with ':' so each query is (lid like "{lid_root}{suffix}*").
    The MCP gateway rejects the `start` offset parameter, so when summary.hits
    exceeds the returned rows the only recovery is to subdivide the LID pattern.
    """
    if not lid_root.endswith(":"):
        raise ValueError("lid_root must end with ':' for exhaustive LID queries")
    if not charset:
        raise ValueError("charset must be non-empty")

    from collections import deque

    diagnostics: dict[str, Any] = {
        "exhaustive_lid_root": lid_root,
        "exhaustive_charset_size": len(charset),
        "exhaustive_api_calls": 0,
        "exhaustive_unique_products": 0,
        "exhaustive_cap_hit": False,
    }

    queue: deque[str] = deque([""])
    collected: dict[str, dict[str, Any]] = {}

    while queue and diagnostics["exhaustive_api_calls"] < max_api_calls:
        suffix = queue.popleft()
        query = f'(lid like "{lid_root}{suffix}*")'
        response = session.get(
            PDS_SEARCH_API,
            params={"q": query, "limit": limit},
            timeout=120,
        )
        response.raise_for_status()
        diagnostics["exhaustive_api_calls"] += 1
        payload = response.json()
        hits = int(payload.get("summary", {}).get("hits") or 0)
        data = list(payload.get("data", []))

        if suffix == "":
            diagnostics["exhaustive_root_hits"] = hits

        if suffix == "" and max_root_hits is not None and hits > max_root_hits:
            diagnostics["exhaustive_skipped_due_to_hit_count"] = hits
            return [], diagnostics

        if hits == 0:
            continue

        truncated = hits > len(data) or (len(data) == limit and hits > limit)
        if truncated:
            if diagnostics["exhaustive_api_calls"] + len(charset) > max_api_calls:
                diagnostics["exhaustive_cap_hit"] = True
                break
            for ch in charset:
                queue.append(suffix + ch)
            continue

        for product in data:
            product_id = str(product.get("id", ""))
            if product_id:
                collected[product_id] = product

    diagnostics["exhaustive_unique_products"] = len(collected)
    return list(collected.values()), diagnostics


def _exhaustive_lid_charset_for_mission(mission: str) -> Optional[str]:
    if mission == "NEAR_1998":
        return "0123456789abcdef"
    if mission == "MESSENGER_2005":
        return "0123456789abcdefghijklmnopqrstuvwxyz_"
    return None


def search_collection_tracking_products(
    session: requests.Session,
    mission: str,
    perigee: datetime,
    metadata: dict[str, Any],
    window_hours: float = 48.0,
    limit: int = 200,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    window_start = perigee - timedelta(hours=window_hours / 2.0)
    window_end = perigee + timedelta(hours=window_hours / 2.0)
    products: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    title_pages = 0
    member_pages = 0
    unindexed_collections: list[str] = []
    exhaustive_reports: list[dict[str, Any]] = []
    exhaustive_roots: set[str] = set()

    for collection_lid in mission_collection_lids(mission, metadata):
        if not _collection_is_indexed(session, collection_lid):
            unindexed_collections.append(collection_lid)
            continue
        title_tokens = _perigee_title_tokens(perigee, window_start, window_end)
        title_products, title_page_count = search_collection_members_by_title(
            session,
            collection_lid,
            title_tokens,
            limit=limit,
        )
        title_pages += title_page_count
        for product in title_products:
            product_id = str(product.get("id", ""))
            if product_id in seen_ids:
                continue
            seen_ids.add(product_id)
            products.append(product)

        member_products, member_page_count = search_collection_members(
            session,
            collection_lid,
            window_start,
            window_end,
            limit=limit,
        )
        member_pages += member_page_count
        for product in member_products:
            product_id = str(product.get("id", ""))
            if product_id in seen_ids:
                continue
            seen_ids.add(product_id)
            products.append(product)

        lid_prefix = collection_lid.split("::", 1)[0]
        if not lid_prefix.endswith(":"):
            lid_prefix = f"{lid_prefix}:"
        for product in search_lid_family_products(
            session,
            lid_prefix,
            window_start,
            window_end,
            limit=min(limit, PDS_MCP_RESULT_CAP),
        ):
            product_id = str(product.get("id", ""))
            if product_id in seen_ids:
                continue
            seen_ids.add(product_id)
            products.append(product)

        charset = _exhaustive_lid_charset_for_mission(mission)
        if charset and lid_prefix not in exhaustive_roots:
            exhaustive_roots.add(lid_prefix)
            exhaustive, exdiag = search_lid_root_exhaustive(
                session,
                lid_prefix,
                charset,
                max_api_calls=500,
                max_root_hits=PDS_EXHAUSTIVE_LID_MAX_ROOT_HITS,
            )
            report = {"lid_root": lid_prefix}
            report.update(exdiag)
            if "exhaustive_skipped_due_to_hit_count" not in exdiag:
                for product in exhaustive:
                    if _is_collection_product(product):
                        continue
                    if not _looks_like_tracking_product(product):
                        continue
                    if not _product_overlaps_window(product, window_start, window_end):
                        continue
                    product_id = str(product.get("id", ""))
                    if product_id in seen_ids:
                        continue
                    seen_ids.add(product_id)
                    products.append(product)
            exhaustive_reports.append(report)

    return products, {
        "title_search_pages": title_pages,
        "member_search_pages": member_pages,
        "unindexed_collection_lids": unindexed_collections,
        "exhaustive_lid_searches": exhaustive_reports,
    }


def _infer_product_lid(collection_lid: str, label_stem: str) -> str:
    base = collection_lid.split("::", 1)[0]
    return f"{base}:{label_stem.lower()}::1.0"


def download_label_companion_files(
    session: requests.Session,
    mission_dir: Path,
    metadata: dict[str, Any],
    mission: str,
    perigee: Optional[datetime] = None,
    window_hours: float = 48.0,
    project_root: Optional[Path] = None,
) -> list[dict[str, Any]]:
    from scripts.utils.dsn_tracking_discovery import label_covers_perigee

    downloads: list[dict[str, Any]] = []
    collection_lids = mission_collection_lids(mission, metadata)
    if not collection_lids:
        return downloads

    for label_path in sorted(mission_dir.rglob("*.lbl")):
        if perigee is not None:
            try:
                if not label_covers_perigee(label_path, perigee, window_hours=window_hours):
                    continue
            except ValueError:
                continue
        data_path = label_path.with_suffix(".dat")
        if data_path.exists() and data_path.stat().st_size > 0:
            continue
        product = fetch_product(session, _infer_product_lid(collection_lids[0], label_path.stem))
        downloads.extend(
            download_product_files(
                session,
                product,
                mission_dir,
                project_root=project_root,
            )
        )
    return downloads


def search_tracking_products(
    session: requests.Session,
    mission: str,
    perigee: datetime,
    metadata: dict[str, Any] | None = None,
    window_hours: float = 48.0,
    limit: int = 100,
) -> list[dict[str, Any]]:
    start = perigee - timedelta(hours=window_hours / 2.0)
    end = perigee + timedelta(hours=window_hours / 2.0)
    products: list[dict[str, Any]] = []
    seen_ids: set[str] = set()

    title_tokens = _perigee_title_tokens(perigee, start, end)
    search_tokens = _mission_search_tokens(mission)
    if metadata:
        search_tokens.extend(metadata_archive_search_tokens(metadata))

    for token in sorted(set(search_tokens + title_tokens)):
        query = _build_time_window_query(start, end, token)
        response = session.get(
            PDS_SEARCH_API,
            params={"q": query, "limit": limit},
            timeout=60,
        )
        response.raise_for_status()
        payload = response.json()
        url = payload.get("links", {}).get("next")
        for _page in range(50):
            for product in payload.get("data", []):
                product_id = str(product.get("id", ""))
                if product_id in seen_ids:
                    continue
                if not _looks_like_tracking_product(product):
                    continue
                if not _product_overlaps_window(product, start, end):
                    continue
                seen_ids.add(product_id)
                products.append(product)
            if not url:
                break
            response = session.get(url, timeout=60)
            response.raise_for_status()
            payload = response.json()
            url = payload.get("links", {}).get("next")
    return products


def ingest_mission_tracking(
    project_root: Path,
    mission: str,
    window_hours: float = 48.0,
) -> dict[str, Any]:
    mission_dir = project_root / "data" / "raw" / "dsn_tracking" / mission
    mission_dir.mkdir(parents=True, exist_ok=True)
    metadata = load_mission_metadata(project_root, mission)
    if mission in REFERENCE_PRODUCT_LIDS:
        perigee = None
    else:
        perigee = load_perigee_datetime(project_root, mission)

    session = requests.Session()
    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": "TEP-EFA-Pipeline/1.0 (Academic Research)",
        }
    )

    manifest: dict[str, Any] = {
        "mission": mission,
        "perigee_utc": perigee.isoformat() if perigee is not None else None,
        "window_hours": window_hours,
        "pds_search_api": PDS_SEARCH_API,
        "metadata_source": (
            _relative_path(project_root, mission_dir / "metadata.json")
            if metadata
            else None
        ),
        "downloads": [],
        "candidate_products": [],
        "search_diagnostics": {
            "collection_lids": mission_collection_lids(mission, metadata),
            "title_search_pages": 0,
            "member_search_pages": 0,
        },
        "status": "no_indexed_products",
    }

    try:
        if perigee is not None:
            products = search_tracking_products(
                session,
                mission,
                perigee,
                metadata=metadata,
                window_hours=window_hours,
            )
            if not products:
                products, collection_diag = search_collection_tracking_products(
                    session,
                    mission,
                    perigee,
                    metadata,
                    window_hours=window_hours,
                )
                manifest["search_diagnostics"].update(collection_diag)
        else:
            products = []
    except requests.RequestException as exc:
        manifest["status"] = "pds_search_failed"
        manifest["error"] = str(exc)
        local_njpl = _scan_local_njpl_trk234_files(mission_dir, project_root)
        if local_njpl:
            manifest["local_njpl_trk234_artifacts"] = local_njpl
            manifest["status"] = "pds_search_failed_local_njpl_verified"
        _write_manifest(mission_dir, manifest)
        return manifest

    if not products and mission in REFERENCE_PRODUCT_LIDS:
        for product_lid in REFERENCE_PRODUCT_LIDS[mission]:
            products.append(fetch_product(session, product_lid))

    manifest["candidate_products"] = [
        {
            "id": product.get("id"),
            "title": product.get("title"),
            "start_date_time": product.get("start_date_time"),
            "stop_date_time": product.get("stop_date_time"),
        }
        for product in products
    ]

    for product in products:
        if _is_collection_product(product):
            continue
        try:
            manifest["downloads"].extend(
                download_product_files(
                    session,
                    product,
                    mission_dir,
                    project_root=project_root,
                )
            )
        except requests.RequestException as exc:
            manifest.setdefault("download_errors", []).append(
                {"product_id": product.get("id"), "error": str(exc)}
            )

    if manifest["downloads"]:
        manifest["status"] = "downloaded"
    elif manifest["candidate_products"]:
        manifest["status"] = "candidates_without_downloadable_files"
    try:
        manifest["downloads"].extend(
            download_label_companion_files(
                session,
                mission_dir,
                metadata,
                mission,
                perigee=perigee,
                window_hours=window_hours,
                project_root=project_root,
            )
        )
        if manifest["downloads"] and manifest["status"] != "downloaded":
            manifest["status"] = "downloaded"
    except requests.RequestException as exc:
        manifest.setdefault("label_download_errors", []).append(str(exc))

    local_njpl = _scan_local_njpl_trk234_files(
        mission_dir, project_root, perigee=perigee, window_hours=window_hours
    )
    if local_njpl:
        manifest["local_njpl_trk234_artifacts"] = local_njpl
        if manifest["status"] == "no_indexed_products":
            manifest["status"] = "no_indexed_products_local_njpl_verified"
        elif manifest["status"] == "pds_search_failed":
            manifest["status"] = "pds_search_failed_local_njpl_verified"
        if any(row.get("overlaps_perigee_window") for row in local_njpl):
            manifest["perigee_window_local_njpl"] = True

    if mission == "Juno_2013" and perigee is not None:
        try:
            from scripts.utils.atmospheres_tnf_index import probe_atmospheres_tnf_index

            manifest["atmospheres_index_probe"] = probe_atmospheres_tnf_index(
                session, perigee, window_hours=window_hours
            )
        except requests.RequestException as exc:
            manifest["atmospheres_index_probe"] = {
                "status": "probe_failed",
                "error": str(exc),
            }

    _write_manifest(mission_dir, manifest)
    return manifest


def ingest_catalog_missions(
    project_root: Path,
    missions: list[str],
    window_hours: float = 48.0,
) -> dict[str, Any]:
    summary = {"missions": {}, "status": "complete"}
    for mission in missions:
        summary["missions"][mission] = ingest_mission_tracking(
            project_root,
            mission,
            window_hours=window_hours,
        )
    return summary


def _write_manifest(mission_dir: Path, manifest: dict[str, Any]) -> None:
    manifest_path = mission_dir / "ingest_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


def _scan_local_njpl_trk234_files(
    mission_dir: Path,
    project_root: Path,
    *,
    perigee: Optional[datetime] = None,
    window_hours: float = 48.0,
) -> list[dict[str, Any]]:
    """
    Discover TRK-2-34 / NJPL-class tracking files already on disk under ``mission_dir``.

    These are not synthetic; they are listed for provenance when the PDS search API
    returns no indexed products (common for some collection LIDs).
    """
    from scripts.utils.dsn_tracking_discovery import (
        _label_path_for_data_file,
        is_label_only_dat_file,
        is_trk234_archive,
    )
    from scripts.utils.pds3_lbl_time import interval_overlaps_window, parse_pds3_lbl_start_stop_utc

    win_lo = win_hi = None
    if perigee is not None:
        if perigee.tzinfo is None:
            perigee = perigee.replace(tzinfo=timezone.utc)
        win_lo = perigee - timedelta(hours=window_hours / 2.0)
        win_hi = perigee + timedelta(hours=window_hours / 2.0)

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for pattern in ("*.dat", "*.DAT", "*.trk", "*.TRK", "*.tnf", "*.TNF"):
        for path in mission_dir.rglob(pattern):
            if not path.is_file():
                continue
            if is_label_only_dat_file(path):
                continue
            if not is_trk234_archive(path):
                continue
            rel = path.resolve().relative_to(project_root.resolve())
            key = str(rel)
            if key in seen:
                continue
            seen.add(key)
            st = path.stat()
            row: dict[str, Any] = {
                "path": key.replace("\\", "/"),
                "bytes": int(st.st_size),
                "format_magic": "NJPL",
                "note": "Verified on disk; TRK-2-34 SFDU archive (not from PDS API ingest row).",
            }
            if win_lo is not None and win_hi is not None:
                lbl = _label_path_for_data_file(path)
                if lbl is not None and lbl.is_file():
                    try:
                        span = parse_pds3_lbl_start_stop_utc(
                            lbl.read_text(encoding="utf-8", errors="ignore")
                        )
                        if span is not None:
                            t0, t1 = span
                            row["lbl_start_utc"] = t0.isoformat()
                            row["lbl_stop_utc"] = t1.isoformat()
                            row["overlaps_perigee_window"] = interval_overlaps_window(
                                t0, t1, win_lo, win_hi
                            )
                        else:
                            row["overlaps_perigee_window"] = None
                    except (ValueError, OSError):
                        row["overlaps_perigee_window"] = None
                else:
                    row["overlaps_perigee_window"] = None
            rows.append(row)
    rows.sort(key=lambda r: r["path"])
    return rows
