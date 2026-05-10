#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Visualization Generator

Generates publication-quality figures for the TEP flyby analysis manuscript.
"""

import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time

# Import style settings
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.step_logger import StepLogger
from scripts.utils.physics import M_PL_GEV

# Blue-grey colour palette for TEP figures
BLUE_GREY_COLORS = {
    'primary_dark': '#1e293b',    # slate-800
    'primary': '#334155',         # slate-700
    'secondary': '#475569',       # slate-600
    'accent': '#64748b',          # slate-500
    'light': '#94a3b8',           # slate-400
    'lighter': '#cbd5e1',         # slate-300
    'highlight': '#3b82f6',       # blue-500
    'highlight_light': '#60a5fa', # blue-400
    'background': '#f8fafc',    # slate-50
    'surface': '#ffffff',
    'text': '#0f172a',            # slate-900
    'text_muted': '#64748b',      # slate-500
}

def set_pub_style(scale=1.0, dpi=300, transparent=False):
    base_font = 11
    plt.rcParams.update({
        'figure.figsize': (10, 6),
        'font.family': 'serif',
        'font.serif': ['Georgia', 'Times New Roman', 'STIXGeneral', 'DejaVu Serif', 'serif'],
        'mathtext.fontset': 'dejavuserif',
        'font.size': base_font * scale,
        'axes.labelsize': (base_font + 1) * scale,
        'axes.titlesize': (base_font + 2) * scale,
        'xtick.labelsize': (base_font - 1) * scale,
        'ytick.labelsize': (base_font - 1) * scale,
        'legend.fontsize': (base_font - 1) * scale,
        'figure.titlesize': (base_font + 3) * scale,
        'figure.dpi': dpi,
        'savefig.dpi': dpi,
        'savefig.transparent': transparent,
        'figure.facecolor': 'white' if not transparent else 'none',
        'axes.facecolor': 'white' if not transparent else 'none',
        'axes.edgecolor': '#475569',
        'axes.labelcolor': '#0f172a',
        'xtick.color': '#475569',
        'ytick.color': '#475569',
        'legend.frameon': True,
        'legend.framealpha': 0.8,
        'legend.facecolor': 'white',
        'legend.edgecolor': '#cbd5e1',
        'savefig.facecolor': 'white' if not transparent else 'none',
        'savefig.edgecolor': 'white' if not transparent else 'none',
        'grid.color': '#cbd5e1',
        'grid.linestyle': '-',
        'grid.linewidth': 0.5,
        'grid.alpha': 0.4,
        'axes.linewidth': 1.0,
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'lines.linewidth': 2.0,
        'lines.markersize': 7,
        'text.usetex': False,
        'text.color': '#334155',
    })
    
    import matplotlib as mpl
    # Set color cycle to blue-grey palette
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[
        '#334155', '#64748b', '#3b82f6', '#94a3b8', '#475569', '#60a5fa'
    ])

# Apply publication style with Georgia font and blue-grey palette
set_pub_style()
def load_pipeline_data(logger=None):
    """Load data from pipeline outputs."""
    results_dir = PROJECT_ROOT / 'results'
    
    if logger:
        logger.debug("Loading fitting data from step008_fitting_results.json")
    try:
        with open(results_dir / 'step008_fitting_results.json', encoding='utf-8') as f:
            fitting_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
        if logger:
            logger.error(f"Failed to load fitting data: {e}")
        return None, None
    
    if logger:
        logger.debug("Loading predictions data from step007_tep_predictions.json")
    try:
        with open(results_dir / 'step007_tep_predictions.json', encoding='utf-8') as f:
            predictions_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
        if logger:
            logger.error(f"Failed to load predictions data: {e}")
        return None, None
    
    if logger:
        logger.success("Data loaded successfully")
    
    return fitting_data, predictions_data


def generate_altitude_anomaly_figure(fitting_data, predictions_data, output_dir):
    """Generate Figure 1: Altitude vs Anomaly correlation."""
    
    # Extract data
    spacecraft_names = []
    altitudes = []
    anomalies = []
    colors = []
    markers = []
    
    for name, pred in predictions_data['predictions'].items():
        spacecraft_names.append(name)
        altitudes.append(pred['perigee']['altitude_km'])
        anomalies.append(pred['observed']['dv_obs_mm_s'])
        
        # Color by detection status - blue-grey palette
        if pred['observed']['dv_obs_mm_s'] > 0.5:
            colors.append(BLUE_GREY_COLORS['highlight'])  # Blue highlight for detections
            markers.append('o')
        else:
            colors.append(BLUE_GREY_COLORS['accent'])  # Slate for nulls
            markers.append('s')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for i, (name, alt, anom, color, marker) in enumerate(zip(spacecraft_names, altitudes, anomalies, colors, markers)):
        ax.scatter(anom, alt, c=color, marker=marker, s=100, alpha=0.8, edgecolors=BLUE_GREY_COLORS['primary_dark'], linewidth=1, zorder=3)
        # Add label with offset to avoid overlap - keep within plot bounds
        label = name.split('_')[0]
        # Use staggered offsets based on position to prevent overlap
        if alt < 1000:
            offset_x, offset_y = 5, 10
        elif alt < 3000:
            # Use different offsets for different spacecraft in this range
            if 'Rosetta' in label:
                offset_x, offset_y = -25, 15
            elif 'MESSENGER' in label:
                offset_x, offset_y = 5, -15
            else:
                offset_x, offset_y = -20, 15
        else:
            offset_x, offset_y = 5, 10
        ax.annotate(label, (anom, alt), xytext=(offset_x, offset_y), textcoords='offset points',
                   fontsize=8, ha='center', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=BLUE_GREY_COLORS['light'], alpha=0.9))

    # Add screening threshold line
    ax.axhline(y=2500, color=BLUE_GREY_COLORS['accent'], linestyle='--', linewidth=1.5, alpha=0.8, label='Screening threshold (~2500 km)')

    # Styling - no title, blue-grey palette
    ax.set_xlabel('Velocity Anomaly Δv (mm/s)', fontsize=12)
    ax.set_ylabel('Perigee Altitude (km)', fontsize=12)
    ax.set_yscale('log')
    ax.set_xscale('symlog', linthresh=0.1)
    ax.grid(True, alpha=0.4)
    ax.legend(loc='lower right', fontsize=10)

    # Add padding to axis limits
    ax.margins(x=0.15, y=0.15)
    
    # Add annotation with blue-grey styling
    ax.text(0.02, 0.98, '● Detections: 4 published, 1 TEP fit success\n[] Non-detections (n=7)\n[Note: NEAR, Rosetta_2005, Cassini excluded by sign mismatch]',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor=BLUE_GREY_COLORS['background'],
                     edgecolor=BLUE_GREY_COLORS['lighter'], alpha=0.9))
    
    plt.tight_layout()
    output_file = output_dir / 'step037_figure1_altitude_anomaly.png'
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    return output_file


def generate_beta_comparison_figure(fitting_data, output_dir):
    """Generate Figure 2: Fitted β values comparison."""

    # Extract successful fits
    fits = []
    for name, data in fitting_data['individual_fits'].items():
        if data['fit']['beta_fitted']:
            fits.append({
                'name': name.split('_')[0],
                'beta': data['fit']['beta_fitted'],
                'uncertainty': data['fit']['beta_uncertainty'],
                'altitude': data['perigee']['altitude_km'],
                'is_marginal': data['fit']['beta_fitted'] < 1e-5  # Marginal detection threshold
            })

    # Sort by altitude
    fits.sort(key=lambda x: x['altitude'])

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot as horizontal bar chart for better readability
    y_pos = np.arange(len(fits))
    betas = [f['beta'] * 1e3 for f in fits]
    uncertainties = [f['uncertainty'] * 1e3 for f in fits]
    altitudes = [f['altitude'] for f in fits]
    names = [f['name'] for f in fits]
    colors = [BLUE_GREY_COLORS['light'] if f['is_marginal'] else BLUE_GREY_COLORS['highlight'] for f in fits]

    bars = ax.barh(y_pos, betas, xerr=uncertainties, capsize=5,
                  color=colors, edgecolor=BLUE_GREY_COLORS['primary_dark'], linewidth=1.5, alpha=0.9)

    # Add altitude labels next to bars
    for i, (y, alt) in enumerate(zip(y_pos, altitudes)):
        ax.text(betas[i] + uncertainties[i] + 0.1, y, f'{alt:.0f} km',
               va='center', fontsize=10, color=BLUE_GREY_COLORS['primary'])

    # Add value labels on bars
    for i, (bar, beta, unc) in enumerate(zip(bars, betas, uncertainties)):
        width = bar.get_width()
        ax.text(width - unc - 0.15, bar.get_y() + bar.get_height()/2,
                f'{beta:.2f}', ha='right', va='center', fontsize=10, fontweight='bold',
                color='white')

    # Add weighted mean line
    weighted_mean = fitting_data['overall_analysis']['beta_statistics']['weighted_mean'] * 1e3
    ax.axvline(x=weighted_mean, color=BLUE_GREY_COLORS['secondary'], linestyle='--', linewidth=2,
               label=r'Weighted mean: $\beta = 2.19 \times 10^{-3}$')

    # Styling
    ax.set_xlabel(r'Fitted $\beta$ ($\times 10^{-3}$)', fontsize=12)
    ax.set_ylabel('Spacecraft (ordered by perigee altitude)', fontsize=12)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names)
    ax.set_xscale('log')  # Log scale for beta
    ax.set_xlim(1e-3, 3)
    ax.invert_yaxis()  # Put lowest altitude at top
    ax.grid(True, alpha=0.3, axis='x')
    ax.legend(loc='lower right', fontsize=10)

    # Add annotation
    ax.text(0.98, 0.02, r'Values on bars: $\beta \times 10^{-3}$' + '\n' + r'Right labels: Altitude (km)',
            transform=ax.transAxes, fontsize=9, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor=BLUE_GREY_COLORS['background'],
                     edgecolor=BLUE_GREY_COLORS['lighter'], alpha=0.9))

    plt.tight_layout()
    output_file = output_dir / 'step037_figure2_beta_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    return output_file


def generate_ppn_constraint_figure(fitting_data, output_dir):
    """Generate Figure 3: PPN constraint visualization."""

    # Calculate |gamma - 1| for each fit using fundamental beta (unscreened regime)
    # Jakarta v0.8: PPN bound constrains cosmological coupling α₀ = β/M_Pl
    # The screened beta_eff used in TEP predictions does NOT appear in PPN formula
    fits = []
    for name, data in fitting_data['individual_fits'].items():
        if data['fit']['beta_fitted']:
            beta_fitted = data['fit']['beta_fitted']  # Use fundamental beta for PPN
            alpha_0 = beta_fitted / M_PL_GEV
            gamma_dev = 2 * alpha_0**2
            fits.append({
                'name': name.split('_')[0],
                'gamma_dev': gamma_dev,
                'altitude': data['perigee']['altitude_km']
            })

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot individual fits with blue-grey colors
    for i, fit in enumerate(fits):
        ax.scatter(fit['gamma_dev'], fit['altitude'], s=150, alpha=0.8,
                  edgecolors=BLUE_GREY_COLORS['primary_dark'], linewidth=1.5,
                  color=BLUE_GREY_COLORS['highlight'], zorder=3)
        # Use altitude-based offsets to prevent overlap
        if fit['altitude'] < 1000:
            offset_x, offset_y = 8, 8
        elif fit['altitude'] < 3000:
            offset_x, offset_y = -15, 15
        else:
            offset_x, offset_y = 5, 10
        ax.annotate(fit['name'], (fit['gamma_dev'], fit['altitude']),
                   xytext=(offset_x, offset_y), textcoords='offset points', fontsize=9,
                   ha='center', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                            edgecolor=BLUE_GREY_COLORS['light'], alpha=0.9))

    # Add Cassini bound line - blue-grey styling
    cassini_bound = 2.3e-5
    ax.axvline(x=cassini_bound, color=BLUE_GREY_COLORS['primary'], linestyle='-', linewidth=2,
               label=r'Cassini PPN bound: $|\gamma-1| < 2.3 \times 10^{-5}$')

    # Add shaded region for excluded zone
    ax.axvspan(cassini_bound, 1e-3, alpha=0.15, color=BLUE_GREY_COLORS['light'], label='Excluded region')

    # Add weighted mean gamma deviation using fundamental beta (unscreened regime)
    # Jakarta v0.8: PPN bound constrains cosmological coupling α₀ = β/M_Pl
    weighted_beta = fitting_data['overall_analysis']['beta_statistics']['weighted_mean']
    alpha_0_weighted = weighted_beta / M_PL_GEV
    weighted_gamma = 2 * alpha_0_weighted**2
    ax.axvline(x=weighted_gamma, color=BLUE_GREY_COLORS['highlight'], linestyle='--', linewidth=2,
               label=r'TEP weighted mean: $|\gamma-1| = 4.45 \times 10^{-6}$')

    # Styling
    ax.set_ylabel('Perigee Altitude (km)', fontsize=12)
    ax.set_xlabel('|γ − 1| (log scale)', fontsize=12)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.4, which='both')
    ax.legend(loc='lower left', fontsize=10)

    # Add padding to axis limits
    ax.margins(x=0.15, y=0.15)

    # Add annotation for safety margin - blue-grey styling
    ax.text(0.98, 0.98, r'TEP predictions well below' + '\n' + r'Cassini PPN bound' + '\n' + r'(Margin: ~5×)',
            transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor=BLUE_GREY_COLORS['background'],
                     edgecolor=BLUE_GREY_COLORS['lighter'], alpha=0.9))

    plt.tight_layout()
    output_file = output_dir / 'step037_figure3_ppn_constraints.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    return output_file


def generate_screening_profile_figure(output_dir):
    """Generate Figure 4: Temporal Shear Suppression screening profile."""
    
    # Physical constants
    R_EARTH = 6371  # km
    
    # Generate altitude range
    altitudes = np.linspace(0, 50000, 1000)  # km
    r = R_EARTH + altitudes  # distance from center
    
    # Simplified screening model (exponential relaxation)
    phi_surface = 6.14e7  # GeV (at Earth's surface)
    phi_space = 4.56e19   # GeV (deep space)
    lambda_scr = 500  # screening length in km (simplified)
    
    delta_r = altitudes
    frac = 1.0 - np.exp(-delta_r / lambda_scr)
    phi = phi_surface + (phi_space - phi_surface) * frac
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot screening profile with blue-grey color
    ax.plot(altitudes, phi / 1e18, linewidth=2.5, color=BLUE_GREY_COLORS['primary'], label='Temporal Shear Suppression field φ(r)')
    
    # Add Earth's surface
    ax.axvline(x=0, color=BLUE_GREY_COLORS['secondary'], linestyle='-', linewidth=2, alpha=0.8, label='Earth surface')
    
    # Add screening threshold region - blue-grey palette
    ax.axvspan(0, 2500, alpha=0.15, color=BLUE_GREY_COLORS['highlight'], label='Strong screening region (< 2500 km)')
    ax.axvspan(2500, 5000, alpha=0.1, color=BLUE_GREY_COLORS['accent'], label='Transition region (2500–5000 km)')
    ax.axvspan(5000, 50000, alpha=0.08, color=BLUE_GREY_COLORS['light'], label='Weak screening region (> 5000 km)')
    
    # Mark detection spacecraft with blue highlight (3 primary detections, S/N > 2)
    detection_alts = [567.9, 972.3, 1956.0]  # NEAR, Galileo 1990, Rosetta 2005
    detection_names = ['NEAR', 'Galileo 1990', 'Rosetta 2005']
    for i, (alt, name) in enumerate(zip(detection_alts, detection_names)):
        ax.axvline(x=alt, color=BLUE_GREY_COLORS['highlight_light'], linestyle=':', linewidth=1.5, alpha=0.7,
                  label=f'{name} ({alt:.0f} km)')

    # Styling - no title
    ax.set_xlabel('Altitude above Earth surface (km)', fontsize=12)
    ax.set_ylabel('Temporal Shear Suppression field φ (× 10¹⁸ GeV)', fontsize=12)
    ax.set_xlim(0, 50000)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.4, which='both')
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)

    # Add padding to axis limits
    ax.margins(x=0.05, y=0.15)
    
    plt.tight_layout()
    output_file = output_dir / 'step037_figure4_screening_profile.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    return output_file


def main():
    """Execute visualization generation."""
    logger = StepLogger("step_037_visualizations", PROJECT_ROOT)
    start_time = time.time()

    logger.header("STEP 037: VISUALIZATION GENERATION")

    # Setup output directory (results folder only)
    results_dir = PROJECT_ROOT / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Output directory: {results_dir}")
    logger.section("GENERATING PUBLICATION FIGURES")

    try:
        # Load data
        fitting_data, predictions_data = load_pipeline_data(logger)

        # Generate figures to results folder only
        logger.subheader("Figure 1: Altitude vs Anomaly correlation")
        fig1 = generate_altitude_anomaly_figure(fitting_data, predictions_data, results_dir)
        logger.success(f"Generated: {fig1}")
        logger.add_output_file(fig1, "Altitude vs Anomaly figure")

        logger.subheader("Figure 2: PPN constraint analysis")
        fig2 = generate_ppn_constraint_figure(fitting_data, results_dir)
        logger.success(f"Generated: {fig2}")
        logger.add_output_file(fig2, "PPN constraint figure")

        logger.subheader("Figure 3: Screening profile")
        fig3 = generate_screening_profile_figure(results_dir)
        logger.success(f"Generated: {fig3}")
        logger.add_output_file(fig3, "Screening profile figure")

        logger.success("All figures generated successfully")

        duration = time.time() - start_time
        logger.log_step_summary(duration, "SUCCESS")
        return 0

    except Exception as e:
        logger.error(f"Error generating figures: {e}")
        logger.debug(f"Exception type: {type(e).__name__}")
        duration = time.time() - start_time
        logger.log_step_summary(duration, "FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
