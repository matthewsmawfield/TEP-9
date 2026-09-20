#!/usr/bin/env python3
"""TEP-9 shared publication figure style.

Single source of truth for colours, fonts, sizes and resolution across all
step figures, matched to the site palette (styles.css) and the conventions
used by the other TEP-* repositories (TEP-DR4, TEP-HC).

Usage: imported and applied automatically by tep9_common.py; step scripts
may also call apply_style() after creating subplots if needed, or import
COLORS / FIG_SIZE for per-figure choices.
"""

TEP_COLORS = {
    'primary': '#1C2E4A',       # deep navy (headers)
    'primary_light': '#2C3E50',
    'secondary': '#5D6D7E',     # slate
    'accent': '#84a3aa',        # muted teal accent
    'highlight': '#1A5276',     # link/highlight blue
    'signal': '#b43b4e',        # reserved signal red (use sparingly)
    'text': '#17202A',
    'gray': '#566573',
    'grid': '#EBEDEF',
    'background': '#FFFFFF',
}

# Canonical sizes (inches) matching the web/PDF templates used across TEP-*.
FIG_SIZE = {
    'single_column': (3.5, 2.6),
    'double_column': (7.2, 4.5),
    'full_width': (10.0, 6.0),
    'web_standard': (7.5, 4.6),
    'web_tall': (7.5, 5.4),
    'web_two_panel': (7.5, 4.0),
    'web_quad': (7.5, 5.6),
    'web_dense': (7.5, 6.6),
    'web_grid_3x3': (7.5, 7.2),
}

DPI_PUB = 300

# Ordered colour cycle: navy, slate, accent teal, deep blue, signal red.
PROP_CYCLE = ['#1C2E4A', '#5D6D7E', '#84a3aa', '#1A5276', '#b43b4e']


def apply_style(scale=1.0, dpi=DPI_PUB):
    """Apply the TEP-9 publication style to matplotlib rcParams.

    Idempotent; safe to call at import time or per-figure.  Manuscript
    (hero) figures carry no rendered titles — that text lives in the
    manuscript captions — while supplementary diagnostics keep their
    titles for standalone readability.
    """
    import matplotlib as mpl
    base = 13
    mpl.rcParams.update({
        'figure.figsize': FIG_SIZE['full_width'],
        'figure.dpi': dpi,
        'savefig.dpi': dpi,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
        'savefig.edgecolor': 'white',
        'savefig.transparent': False,
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif',
                       'serif'],
        'mathtext.fontset': 'dejavuserif',
        'font.size': base * scale,
        'axes.labelsize': (base + 1) * scale,
        'axes.titlesize': (base + 2) * scale,
        'figure.titlesize': (base + 3) * scale,
        'xtick.labelsize': (base - 1) * scale,
        'ytick.labelsize': (base - 1) * scale,
        'legend.fontsize': (base - 1) * scale,
        'text.color': '#17202A',
        'axes.labelcolor': '#17202A',
        'axes.edgecolor': '#475569',
        'xtick.color': '#17202A',
        'ytick.color': '#17202A',
        'axes.linewidth': 1.2,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'lines.linewidth': 2.0,
        'lines.markersize': 7,
        'axes.grid': True,
        'grid.color': '#475569',
        'grid.linestyle': '-',
        'grid.linewidth': 0.5,
        'grid.alpha': 0.2,
        'legend.frameon': True,
        'legend.framealpha': 0.8,
        'legend.facecolor': 'white',
        'legend.edgecolor': '#475569',
        'text.usetex': False,
    })
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=PROP_CYCLE)

    return TEP_COLORS
