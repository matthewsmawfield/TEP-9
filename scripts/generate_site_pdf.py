#!/usr/bin/env python3
"""
Generate PDF from TEP Site
==========================

Generates a high-quality PDF from the built TEP site static HTML.
Uses the html_to_pdf.py converter with optimized settings for academic manuscripts.

Usage:
    python scripts/generate_site_pdf.py
    python scripts/generate_site_pdf.py --quality high --wait-time 5
"""

import asyncio
import argparse
import subprocess
import sys
import re
from pathlib import Path

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent / 'utils'))

try:
    import yaml
except ImportError:
    yaml = None

from html_to_pdf import HTMLToPDFConverter, create_preset_configs


def load_citation_metadata():
    """Load version and codename from CITATION.cff."""
    base_dir = Path(__file__).parent.parent
    citation_file = base_dir / 'CITATION.cff'
    
    if not citation_file.exists():
        print("⚠️  CITATION.cff not found, using defaults")
        return {'version': 'v0.1', 'codename': 'Addis Ababa', 'title': 'TEP-9'}
    
    try:
        if yaml:
            with open(citation_file, 'r') as f:
                data = yaml.safe_load(f)
            version_str = data.get('version', 'v0.1')
        else:
            # Parse manually if yaml not available
            with open(citation_file, 'r') as f:
                content = f.read()
            version_match = re.search(r'version:\s*"?([^"\n]+)"?', content)
            version_str = version_match.group(1).strip() if version_match else 'v0.1'
        
        # Parse version string like 'v0.1 (Addis Ababa)'
        pattern = r'^(v?[\d.]+)(?:\s*\(([^)]+)\))?$'
        match = re.match(pattern, version_str.strip())
        
        if match:
            version = match.group(1).lstrip('v')
            codename = match.group(2) or 'Addis Ababa'
        else:
            version = version_str.lstrip('v')
            codename = 'Addis Ababa'
        
        # Paper number: VERSION.json "paper" field (e.g. "36"), fallback to
        # manifest paper_series ("TEP Series: Paper 36 ...") then no prefix.
        paper = ''
        version_file = Path(__file__).parent.parent / 'VERSION.json'
        try:
            import json as _json
            paper = str(_json.load(open(version_file)).get('paper', '')).strip()
        except Exception:
            pass
        if not paper:
            try:
                import json as _json
                mfile = Path(__file__).parent.parent / 'site' / 'manifest.json'
                m = re.search(r'Paper\s+(\d+)', _json.load(open(mfile)).get('paper_series', ''))
                paper = m.group(1) if m else ''
            except Exception:
                pass

        return {'version': version, 'codename': codename, 'paper': paper}
        
    except Exception as e:
        print(f"⚠️  Error parsing CITATION.cff: {e}, using defaults")
        return {'version': 'v0.1', 'codename': 'Addis Ababa', 'title': 'TEP-9'}


def build_static_site():
    """Build the static site first."""
    print("🔨 Building static site...")
    site_dir = Path(__file__).parent.parent / 'site'
    build_script = site_dir / 'build.js'
    
    if not build_script.exists():
        print("❌ Build script not found:", build_script)
        return False
    
    try:
        result = subprocess.run(
            ['node', str(build_script)],
            cwd=site_dir,
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print("❌ Site build failed:")
        print(e.stderr)
        return False


def copy_pdf_to_docs(source_pdf: Path, docs_dir: Path):
    """Copy PDF to the docs directory with proper naming."""
    docs_dir.mkdir(parents=True, exist_ok=True)
    
    # Load metadata from CITATION.cff
    metadata = load_citation_metadata()
    version_str = f"v{metadata['version']}-{metadata['codename'].replace(' ', '-')}"
    
    # Primary PDF name (corpus convention: {paper}-TEP-9-{version}-{codename}.pdf)
    prefix = f"{metadata['paper']}-" if metadata.get('paper') else ''
    target_name = f"{prefix}TEP-9-{version_str}.pdf"
    target_path = docs_dir / target_name
    
    # Copy the file
    import shutil
    shutil.copy2(source_pdf, target_path)
    
    print(f"📄 Copied to: {target_path}")
    print(f"   Size: {target_path.stat().st_size / (1024*1024):.2f} MB")
    
    return target_path


def copy_pdf_to_root(source_pdf: Path, base_dir: Path):
    """Copy PDF to the project root directory."""
    # Load metadata from CITATION.cff
    metadata = load_citation_metadata()
    version_str = f"v{metadata['version']}-{metadata['codename'].replace(' ', '-')}"
    
    # Primary PDF name (corpus convention: {paper}-TEP-9-{version}-{codename}.pdf)
    prefix = f"{metadata['paper']}-" if metadata.get('paper') else ''
    target_name = f"{prefix}TEP-9-{version_str}.pdf"
    target_path = base_dir / target_name
    
    # Copy the file
    import shutil
    shutil.copy2(source_pdf, target_path)
    
    print(f"📄 Copied to root: {target_path}")
    print(f"   Size: {target_path.stat().st_size / (1024*1024):.2f} MB")
    
    return target_path


def process_pdf_with_metadata(pdf_path: Path):
    """Run the PDF processing script to add metadata and compress."""
    import json
    from pypdf import PdfReader, PdfWriter
    base = Path(__file__).resolve().parent.parent
    manifest = json.loads((base / 'site/manifest.json').read_text())
    reader = PdfReader(pdf_path)
    writer = PdfWriter()
    writer.clone_document_from_reader(reader)
    writer.add_metadata({
        '/Title': manifest['title'],
        '/Author': manifest['author'],
        '/Subject': 'Outer-Solar-System catalogue tests and TEP domain-boundary hypothesis',
        '/Keywords': 'TEP, trans-Neptunian objects, comet orbits, spatial selection, numerical audit',
    })
    temporary = pdf_path.with_suffix('.metadata.tmp.pdf')
    with temporary.open('wb') as stream:
        writer.write(stream)
    temporary.replace(pdf_path)
    print('PDF bibliographic metadata embedded')


async def generate_pdf(quality: str = 'high', wait_time: float = 5.0, skip_build: bool = False):
    """Generate PDF from the TEP site."""
    
    # Build site first unless skipped
    if not skip_build:
        if not build_static_site():
            print("❌ Cannot generate PDF - site build failed")
            return False
    
    # Paths
    base_dir = Path(__file__).parent.parent
    dist_dir = base_dir / 'site' / 'dist'
    docs_dir = base_dir / 'site' / 'public' / 'docs'
    
    if not dist_dir.exists():
        print(f"❌ Dist directory not found: {dist_dir}")
        print("   Run site build first: node site/build.js")
        return False
    
    html_file = dist_dir / 'index.html'
    if not html_file.exists():
        print(f"❌ HTML file not found: {html_file}")
        return False
    
    # Output path (temp, will be copied to docs)
    output_pdf = dist_dir / 'manuscript.pdf'
    
    # Select preset based on quality
    # Note: scale controls content zoom (affects page count)
    # device_scale_factor controls pixel density (affects image/text sharpness)
    presets = create_preset_configs()
    if quality == 'maximum':
        # Highest resolution for archival/print quality
        options = presets['high_quality'].copy()
        options['scale'] = 0.90  # Larger content, target <20 pages
        options['device_scale_factor'] = 3.0  # High pixel density for sharpness
        options['viewport'] = {'width': 1920, 'height': 1080}
        options['prefer_css_page_size'] = True
    elif quality == 'high':
        options = presets['high_quality'].copy()
        options['scale'] = 0.90
        options['device_scale_factor'] = 2.5
        options['viewport'] = {'width': 1920, 'height': 1080}
        options['prefer_css_page_size'] = True
    elif quality == 'print':
        options = presets['print_ready'].copy()
        options['scale'] = 0.90
        options['device_scale_factor'] = 2.0
    else:
        options = presets['web_optimized'].copy()
    
    options['custom_css'] = '''
    @media print {
        .header-image, .header-image img { max-height: 135px !important; object-fit: cover; }
        .abstract { break-inside: avoid !important; }
        body, main, article { padding-bottom: 0 !important; margin-bottom: 0 !important; }
        .manuscript-section:last-child, .manuscript-section:last-child > p:last-child { margin-bottom: 0 !important; padding-bottom: 0 !important; }
        p { orphans: 3; widows: 3; }
    }
    '''
    options['wait_time'] = wait_time
    options['format'] = 'A4'
    options['margin_top'] = '1.2cm'
    options['margin_bottom'] = '1.5cm'  # Increased for footer
    options['margin_left'] = '1cm'
    options['margin_right'] = '1cm'
    
    # Add page numbers footer
    options['display_header_footer'] = True
    options['header_template'] = '<div></div>'  # Empty header
    options['footer_template'] = '''
        <div style="font-size:9px; text-align:center; width:100%; color:#555555; font-family:system-ui,-apple-system,sans-serif; padding-bottom:5mm;">
            Page <span class="pageNumber"></span> of <span class="totalPages"></span>
        </div>
    '''
    
    print(f"\n📄 Generating PDF from: {html_file}")
    print(f"   Quality: {quality}")
    print(f"   Wait time: {wait_time}s (for MathJax rendering)")
    
    async with HTMLToPDFConverter() as converter:
        success = await converter.convert_file(
            str(html_file),
            str(output_pdf),
            options
        )
        
        if not success:
            print("❌ PDF generation failed")
            return False
        
        print(f"✅ PDF generated: {output_pdf}")
        print(f"   Size: {output_pdf.stat().st_size / (1024*1024):.2f} MB")
        
        # Load metadata for display
        metadata = load_citation_metadata()
        version_str = f"v{metadata['version']}_{metadata['codename'].replace(' ', '_')}"
        print(f"   Version: {version_str}")
        
        # Copy to docs directory
        final_pdf = copy_pdf_to_docs(output_pdf, docs_dir)
        
        # Copy to root directory (from uncompressed source)
        copy_pdf_to_root(output_pdf, base_dir)
        
        # Process with metadata (compresses the docs version)
        process_pdf_with_metadata(final_pdf)
        
        # Synchronize every downloadable copy after metadata processing
        if final_pdf.exists():
            import shutil
            root_pdf = base_dir / final_pdf.name
            shutil.copy2(final_pdf, root_pdf)
            for folder in (dist_dir / 'docs', dist_dir / 'public' / 'docs'):
                folder.mkdir(parents=True, exist_ok=True)
                shutil.copy2(final_pdf, folder / final_pdf.name)
            shutil.copy2(final_pdf, output_pdf)
            print(f"PDF copies synchronized: {root_pdf}")
            print(f"   Size: {root_pdf.stat().st_size / (1024*1024):.2f} MB")
        
        print(f"\n✅ Complete! PDF available at:")
        print(f"   {final_pdf}")
        print(f"   {base_dir / final_pdf.name}")
        
        return True


def main():
    parser = argparse.ArgumentParser(
        description='Generate PDF from TEP site',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate with default settings
  python scripts/generate_site_pdf.py
  
  # Maximum quality (highest resolution)
  python scripts/generate_site_pdf.py --quality maximum --wait-time 10
  
  # High quality with longer wait for MathJax
  python scripts/generate_site_pdf.py --quality high --wait-time 10
  
  # Skip site build (use existing dist/)
  python scripts/generate_site_pdf.py --skip-build
        """
    )
    
    parser.add_argument(
        '--quality',
        choices=['maximum', 'high', 'print', 'web'],
        default='high',
        help='PDF quality preset (default: high, maximum for highest resolution)'
    )
    parser.add_argument(
        '--wait-time',
        type=float,
        default=5.0,
        help='Seconds to wait for dynamic content (MathJax) to render'
    )
    parser.add_argument(
        '--skip-build',
        action='store_true',
        help='Skip site build and use existing dist/ directory'
    )
    
    args = parser.parse_args()
    
    try:
        success = asyncio.run(generate_pdf(
            quality=args.quality,
            wait_time=args.wait_time,
            skip_build=args.skip_build
        ))
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n⚠️  Cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
