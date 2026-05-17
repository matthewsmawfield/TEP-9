"""
Step 031: PDS API Search for Juno 2013 Earth Flyby Data

Locates Juno 2013 Earth-flyby tracking products via:
- Shared MCP ingest queries (Steps 005/028)
- NMSU Atmospheres OCRU TNF index LBL time-window probe
- Optional Peppi library (iterator slice; no deprecated ``.limit()``)
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.atmospheres_tnf_index import probe_atmospheres_tnf_index
from scripts.utils.dsn_pds_ingest import PDS_SEARCH_API, search_tracking_products
from scripts.utils.step_logger import StepLogger

JUNO_PERIGEE_UTC = datetime(2013, 10, 9, 19, 21, tzinfo=timezone.utc)
FLYBY_WINDOW_HOURS = 48.0


def search_with_mcp_api(logger: StepLogger) -> dict[str, Any]:
    """Search PDS MCP using the same ingest queries as Steps 005/028."""
    import requests

    logger.subsection("PDS MCP Search API (shared ingest)")
    logger.info(f"Endpoint: {PDS_SEARCH_API}")

    session = requests.Session()
    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": "TEP-EFA-Pipeline/1.0",
        }
    )

    metadata_path = PROJECT_ROOT / "data/raw/dsn_tracking/Juno_2013/metadata.json"
    metadata: dict[str, Any] = {}
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    products = search_tracking_products(
        session,
        "Juno_2013",
        JUNO_PERIGEE_UTC,
        metadata=metadata,
        window_hours=FLYBY_WINDOW_HOURS,
        limit=100,
    )
    logger.info(f"search_tracking_products: {len(products)} perigee-window candidate(s)")

    mcp_rows = [
        {
            "id": product.get("id", "N/A"),
            "title": product.get("title", "N/A"),
            "start_date": product.get("start_date_time"),
            "stop_date": product.get("stop_date_time"),
            "strategy": "search_tracking_products",
        }
        for product in products
    ]

    atm_probe = probe_atmospheres_tnf_index(
        session, JUNO_PERIGEE_UTC, window_hours=FLYBY_WINDOW_HOURS
    )
    logger.info(
        f"Atmospheres index: {atm_probe.get('tnf_listed', 0)} TNF listed, "
        f"{atm_probe.get('perigee_window_overlaps', 0)} LBL overlap(s) in flyby window"
    )

    return {
        "api": PDS_SEARCH_API,
        "perigee_utc": JUNO_PERIGEE_UTC.isoformat(),
        "mcp_perigee_candidates": len(mcp_rows),
        "mcp_results": mcp_rows,
        "atmospheres_index_probe": atm_probe,
    }


def search_with_peppi(logger: StepLogger) -> dict[str, Any]:
    """Optional Peppi search (slice iterable results; no ``.limit()``)."""
    logger.subsection("Peppi library (optional)")
    try:
        import pds.peppi as pep
    except ImportError:
        logger.info("Peppi not installed (pip install pds.peppi)")
        return {"status": "peppi_not_installed", "search_results": [], "n_total": 0}

    search_results: list[dict[str, Any]] = []
    try:
        client = pep.PDSRegistryClient()
        products = pep.Products(client)
        queries = [
            products.has_target("Earth"),
            products.filter('title like "*juno*"'),
            products.filter('title like "*radio*science*"'),
        ]
        for query in queries:
            try:
                for product in islice(query, 50):
                    search_results.append(
                        {
                            "id": getattr(product, "id", "N/A"),
                            "title": getattr(product, "title", "N/A"),
                            "type": getattr(product, "product_type", "N/A"),
                            "strategy": "peppi_islice",
                        }
                    )
            except Exception as exc:
                logger.warning(f"Peppi query failed: {exc}")
    except Exception as exc:
        logger.warning(f"Peppi search failed: {exc}")
        return {"status": "error", "error": str(exc), "search_results": [], "n_total": 0}

    logger.info(f"Peppi products collected: {len(search_results)}")
    return {
        "step": "031_pds_api_search_juno",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "tool": "peppi",
        "search_results": search_results,
        "n_total": len(search_results),
    }


def main() -> int:
    logger = StepLogger("step_031_pds_api_search_juno", PROJECT_ROOT)
    logger.header("STEP 031: PDS API SEARCH FOR JUNO 2013 DATA")

    peppi_results = search_with_peppi(logger)
    mcp_results = search_with_mcp_api(logger)

    output_file = PROJECT_ROOT / "results" / "step031_pds_api_search.json"
    payload = {
        **peppi_results,
        "mcp_primary": mcp_results,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    with open(output_file, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    combined_file = PROJECT_ROOT / "results" / "step031_pds_api_search_combined.json"
    with open(combined_file, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "step": "031_pds_api_search_combined",
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "peppi": peppi_results,
                "mcp_api": mcp_results,
                "status": "complete",
            },
            handle,
            indent=2,
        )

    n_mcp = mcp_results.get("mcp_perigee_candidates", 0)
    n_atm = mcp_results.get("atmospheres_index_probe", {}).get("perigee_window_overlaps", 0)
    logger.section("SEARCH SUMMARY")
    logger.info(f"MCP perigee-window candidates: {n_mcp}")
    logger.info(f"Atmospheres LBL overlaps in flyby window: {n_atm}")
    logger.info(f"Peppi products: {peppi_results.get('n_total', 0)}")
    logger.info(f"Wrote {output_file}")

    if n_mcp == 0 and n_atm == 0:
        logger.warning(
            "No automated perigee-window products found — manual PDS-RN ingest required "
            "(see data/raw/dsn_tracking/Juno_2013/DOWNLOAD_INSTRUCTIONS.txt)"
        )

    logger.log_step_summary(0, "SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
