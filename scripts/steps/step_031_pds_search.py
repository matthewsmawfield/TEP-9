"""
Step 031: PDS API Search for Juno 2013 Earth Flyby Data

This step uses the NASA PDS Search API to locate Juno 2013 Earth flyby
tracking data (TRK-2-25/TRK-2-34 format) from the PDS archives.

Uses:
- PDS API: https://nasa-pds.github.io/pds-api/
- Peppi Python library: pip install pds.peppi
- Direct API calls to https://pds.nasa.gov/api/search/1/

Data Products Sought:
- Collection: JNO-E-RSS-1-EDR (Juno Earth Radio Science)
- Format: TRK-2-25, TRK-2-34, or TNF
- Date: 2013-10-08 to 2013-10-10 (perigee 2013-10-09 19:21 UTC)

Author: TEP-EFA Pipeline
Date: 2026-04-19
"""

import sys
import json
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


def search_with_peppi():
    """Search for Juno data using the Peppi Python library."""
    logger = StepLogger("step_031_pds_api_search_juno", PROJECT_ROOT)
    logger.header("STEP 031: PDS API SEARCH FOR JUNO 2013 DATA")
    
    logger.info("Attempting to search using PDS Peppi library...")
    logger.info("API Reference: https://nasa-pds.github.io/pds-api/")
    logger.info("Peppi Docs: https://nasa-pds.github.io/peppi/")
    
    try:
        import pds.peppi as pep
        
        logger.success("Peppi library available")
        logger.info("Connecting to PDS Registry...")
        
        # Connect to PDS
        client = pep.PDSRegistryClient()
        
        # Search for Juno context
        logger.subsection("Searching for Juno Mission Context")
        context = pep.Context()
        
        try:
            juno_results = context.INSTRUMENT_HOSTS.search("juno")
            logger.info(f"Found {len(juno_results)} Juno instrument host entries")
            
            if juno_results:
                juno = juno_results[0]
                logger.info(f"Juno LID: {juno.lid}")
        except Exception as e:
            logger.warning(f"Could not search instrument hosts: {e}")
            juno = None
        
        # Search for products
        logger.subsection("Searching for Juno Earth Flyby Products")
        products = pep.Products(client)
        
        # Try different search strategies
        search_results = []
        
        # Strategy 1: Search by target (Earth) and mission
        logger.info("Strategy 1: Search by target=Earth, mission=Juno...")
        try:
            earth_products = products.has_target("Earth")
            if juno:
                earth_products = earth_products.has_instrument_host(juno.lid)
            
            # Filter by date range if possible
            results = list(earth_products.limit(50))
            logger.info(f"Found {len(results)} Earth-target products")
            
            for product in results:
                search_results.append({
                    'id': getattr(product, 'id', 'N/A'),
                    'title': getattr(product, 'title', 'N/A'),
                    'type': getattr(product, 'product_type', 'N/A'),
                    'strategy': 'target_earth'
                })
        except Exception as e:
            logger.warning(f"Strategy 1 failed: {e}")
        
        # Strategy 2: Search by collection type (radio science)
        logger.info("Strategy 2: Search for radio science data...")
        try:
            # Search for RSS or radio science
            rss_products = products.filter("radio science")
            results = list(rss_products.limit(20))
            logger.info(f"Found {len(results)} radio science products")
            
            for product in results:
                search_results.append({
                    'id': getattr(product, 'id', 'N/A'),
                    'title': getattr(product, 'title', 'N/A'),
                    'type': getattr(product, 'product_type', 'N/A'),
                    'strategy': 'radio_science'
                })
        except Exception as e:
            logger.warning(f"Strategy 2 failed: {e}")
        
        # Strategy 3: Direct LID search for known collection
        logger.info("Strategy 3: Search for JNO-E-RSS-1-EDR collection...")
        known_collections = [
            "urn:nasa:pds:jno-e-rss-1-edr",
            "urn:nasa:pds:jno-j-rss-1-edr",
            "urn:nasa:pds:juno_radio_science"
        ]
        
        for collection in known_collections:
            try:
                coll_products = products.filter(collection)
                results = list(coll_products.limit(10))
                if results:
                    logger.success(f"Found {len(results)} products in {collection}")
                    for product in results:
                        search_results.append({
                            'id': getattr(product, 'id', 'N/A'),
                            'title': getattr(product, 'title', 'N/A'),
                            'collection': collection,
                            'strategy': 'known_collection'
                        })
            except Exception as e:
                logger.debug(f"Collection {collection} not found: {e}")
        
        # Save results
        output = {
            'step': '031_pds_api_search_juno',
            'timestamp': datetime.now().isoformat(),
            'tool': 'peppi',
            'search_results': search_results,
            'n_total': len(search_results)
        }
        
        output_file = PROJECT_ROOT / 'results' / 'step031_pds_api_search.json'
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=2)
        
        logger.section("SEARCH RESULTS")
        logger.info(f"Total products found: {len(search_results)}")
        logger.info(f"Output: {output_file}")
        
        return output
        
    except ImportError:
        logger.error("Peppi library not installed")
        logger.info("Install with: pip install pds.peppi")
        return {'status': 'peppi_not_installed', 'error': 'pip install pds.peppi'}
        
    except Exception as e:
        logger.error(f"Search failed: {e}")
        return {'status': 'error', 'error': str(e)}


def search_with_direct_api():
    """Search using direct PDS API REST calls."""
    logger = StepLogger("step_031_pds_api_search_juno", PROJECT_ROOT)
    logger.header("STEP 031b: DIRECT PDS API SEARCH")
    
    import requests
    
    # PDS Search API base URL - CORRECT ENDPOINT IS /products
    PDS_API_BASE = "https://pds.nasa.gov/api/search/1/products"
    
    logger.info(f"API Endpoint: {PDS_API_BASE}")
    logger.info("Reference: https://nasa-pds.github.io/pds-api/guides/search/quickstart.html")
    
    try:
        session = requests.Session()
        session.headers.update({
            'Accept': 'application/json',
            'User-Agent': 'TEP-EFA-Pipeline/1.0'
        })
        
        # Test API connection
        logger.subsection("Testing API Connection")
        response = session.get(f"{PDS_API_BASE}?limit=1", timeout=10)
        logger.info(f"API Status: {response.status_code}")
        
        if response.status_code == 200:
            logger.success("PDS API accessible")
            
            # Search strategies
            all_results = []
            
            # Strategy 1: Search for Juno
            logger.subsection("Strategy 1: Search for Juno products")
            search_configs = [
                {'q': 'juno', 'limit': 50},
                {'q': 'juno earth', 'limit': 50},
                {'q': 'juno 2013', 'limit': 50},
                {'q': 'JNO-E-RSS', 'limit': 50},
                {'q': 'radio science', 'limit': 50},
            ]
            
            for config in search_configs:
                try:
                    logger.info(f"  Query: {config['q']}")
                    search_response = session.get(PDS_API_BASE, params=config, timeout=30)
                    
                    if search_response.status_code == 200:
                        data = search_response.json()
                        products = data.get('data', [])
                        logger.info(f"    Found: {len(products)} products")
                        
                        for product in products:
                            all_results.append({
                                'id': product.get('id', 'N/A'),
                                'title': product.get('title', 'N/A'),
                                'type': product.get('type', 'N/A'),
                                'query': config['q'],
                                'start_date': product.get('start_date_time'),
                                'stop_date': product.get('stop_date_time')
                            })
                    else:
                        logger.debug(f"    Query failed: {search_response.status_code}")
                        
                except Exception as e:
                    logger.debug(f"    Query error: {e}")
            
            # Filter for Earth flyby date
            logger.subsection("Filtering for 2013 Earth Flyby")
            from datetime import datetime
            
            flyby_products = []
            for result in all_results:
                start = result.get('start_date', '')
                if start and '2013-10' in start:
                    flyby_products.append(result)
                    logger.info(f"  Found 2013-10 product: {result['title'][:60]}...")
            
            logger.info(f"Total products found: {len(all_results)}")
            logger.info(f"2013-10 products: {len(flyby_products)}")
            
            return {
                'total_products': len(all_results),
                'flyby_products': len(flyby_products),
                'all_results': all_results,
                'flyby_candidates': flyby_products
            }
        else:
            logger.warning(f"API not accessible: {response.status_code}")
            logger.warning(f"Response: {response.text[:200]}")
            
    except Exception as e:
        logger.error(f"API connection failed: {e}")
        import traceback
        logger.debug(traceback.format_exc())
    
    return None


def main():
    """Execute PDS API search for Juno data."""
    logger = StepLogger("step_031_pds_api_search_juno", PROJECT_ROOT)
    
    logger.header("STEP 031: PDS API SEARCH FOR JUNO 2013 DATA")
    logger.info("="*70)
    logger.info("Searching NASA PDS for Juno 2013 Earth Flyby TRK-2-25 data")
    logger.info("="*70)
    
    all_results = {}
    
    # Try Peppi first
    logger.section("Method 1: Peppi Python Library")
    peppi_results = search_with_peppi()
    all_results['peppi'] = peppi_results
    
    # Try direct API
    logger.section("Method 2: Direct PDS REST API")
    api_results = search_with_direct_api()
    all_results['direct_api'] = api_results
    
    # Summary
    logger.section("SEARCH SUMMARY")
    
    if peppi_results.get('n_total', 0) > 0:
        logger.success(f"✓ Peppi found {peppi_results['n_total']} products")
    else:
        logger.warning("✗ Peppi search returned no results")
    
    if api_results:
        logger.success("✓ Direct API connection successful")
    else:
        logger.warning("✗ Direct API not accessible or returned no data")
    
    logger.info("")
    logger.info("NEXT STEPS:")
    logger.info("1. Review search results in results/step031_pds_api_search.json")
    logger.info("2. If products found, download URLs will be provided")
    logger.info("3. If no results, try manual search at https://pds.mcp.nasa.gov/portal/search")
    logger.info("")
    logger.info("Alternative direct contact:")
    logger.info("  pds-rn@jpl.nasa.gov")
    logger.info("  Subject: 'Juno 2013 Earth Flyby TRK-2-25 Data Request'")
    
    # Save combined results
    output_file = PROJECT_ROOT / 'results' / 'step031_pds_api_search_combined.json'
    with open(output_file, 'w') as f:
        json.dump({
            'step': '031_pds_api_search_combined',
            'timestamp': datetime.now().isoformat(),
            'peppi': peppi_results,
            'direct_api': api_results,
            'status': 'complete'
        }, f, indent=2)
    
    logger.info(f"")
    logger.info(f"Combined results: {output_file}")
    
    logger.log_step_summary(0, "SUCCESS")
    return 0


if __name__ == '__main__':
    sys.exit(main())
