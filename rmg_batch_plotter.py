"""
rmg_batch_plotter.py
--------------------
Python translation of RmgBatchPlotter.m by Robert E. Kopp (2008).

This is the core plotting engine that the matRockmag GUI's "Run" button
called.  It handles all combinations of:
  - per-sample vs. multi-sample overlay (Multisample flag)
  - individual figures vs. subplot grid (Subplots flag)
  - autosave to PNG and/or EPS
  - a custom file prefix
  - selective routine execution (any subset of the 9 routine types)

Routine names (matching the GUI's listRoutine):
    'IRM'          – SIRM acquisition + AF demagnetisation
    'dIRM'         – dIRM/dB derivative
    'ARM'          – ARM acquisition
    'LowrieFuller' – Lowrie-Fuller (ARM + IRM)
    'AF'           – Lowrie-Fuller multi-AF variant
    'Fuller'       – Fuller plot
    'RRM'          – RRM
    'Backfield'    – Backfield IRM
    'StatBox'      – Statistics text box

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import os
import math
import matplotlib
import matplotlib.pyplot as plt

from rmg_plots import (
    rmg_sirm_plot,
    rmg_sirm_derivative_plot,
    rmg_arm_plot,
    rmg_lowrie_fuller_plot,
    rmg_fuller_plot,
    rmg_rrm_plot,
    rmg_backfield_plot,
    rmg_stat_box,
)

# All routine names exactly as they appear in the MATLAB GUI
ALL_ROUTINES = ['IRM', 'dIRM', 'ARM', 'LowrieFuller', 'AF',
                'Fuller', 'RRM', 'Backfield', 'StatBox']


# ─────────────────────────────────────────────────────────────────────────────
#  Internal: dispatch one routine → one axes
# ─────────────────────────────────────────────────────────────────────────────

def _make_plot(data_list, routine, af_level=None, ax=None):
    """Render a single routine into ax (or current axes if None)."""
    if ax is None:
        ax = plt.gca()

    r = routine.upper() if routine not in ('LowrieFuller', 'StatBox') else routine

    if r in ('IRM', 'IRMZ', 'IRMX'):
        rmg_sirm_plot(data_list, ax=ax)
    elif r == 'DIRM':
        rmg_sirm_derivative_plot(data_list, ax=ax)
    elif r == 'ARM':
        rmg_arm_plot(data_list, ax=ax)
    elif r == 'LowrieFuller':
        rmg_lowrie_fuller_plot(data_list, ax=ax)
    elif r == 'AF':
        # 'MultiAF' variant: ARM-only Lowrie-Fuller
        rmg_lowrie_fuller_plot(data_list, single_mode=True, ax=ax)
    elif r == 'FULLER':
        rmg_fuller_plot(data_list, ax=ax)
    elif r == 'RRM':
        rmg_rrm_plot(data_list, ax=ax)
    elif r == 'BACKFIELD':
        rmg_backfield_plot(data_list, ax=ax)
    elif r == 'StatBox':
        rmg_stat_box(data_list, ax=ax)


def _make_multiplot(data_list, routines, af_level=None):
    """Create a subplot grid with one panel per routine."""
    n   = len(routines)
    rows = math.ceil(math.sqrt(n))
    cols = math.ceil(n / rows)
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
    axes_flat = [axes] if n == 1 else list(
        axes.flat if hasattr(axes, 'flat') else [axes])
    for i, routine in enumerate(routines):
        _make_plot(data_list, routine, af_level, ax=axes_flat[i])
    # hide any unused panels
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis('off')
    plt.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
#  Internal: autosave helper
# ─────────────────────────────────────────────────────────────────────────────

def _autosave(fig, stem, formats):
    """Save fig to stem.png and/or stem.eps depending on formats list."""
    for fmt in formats:
        if fmt in ('png', 'PNG'):
            fig.savefig(f'{stem}.png', dpi=150, bbox_inches='tight')
            print(f'  saved → {stem}.png')
        elif fmt in ('eps', 'EPS'):
            fig.savefig(f'{stem}.eps', format='eps', bbox_inches='tight')
            print(f'  saved → {stem}.eps')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_batch_plotter  (RmgBatchPlotter.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_batch_plotter(data_list, routines=None, *,
                      multisample=False,
                      subplots=False,
                      autosave_formats=None,
                      af_level=None,
                      file_prefix=''):
    """
    Core plotting engine – equivalent to the MATLAB GUI's "Run" button.

    Parameters
    ----------
    data_list        : list of RmgData  (or single RmgData)
    routines         : list of str, default ALL_ROUTINES
                       Any subset of: 'IRM','dIRM','ARM','LowrieFuller',
                       'AF','Fuller','RRM','Backfield','StatBox'
    multisample      : bool  – overlay all samples on each plot (default False)
    subplots         : bool  – put all routines in one subplot grid (default False)
    autosave_formats : list of str, e.g. ['png'], ['eps'], ['png','eps']
    af_level         : float or list of floats (T) – AF level(s); auto-detected if None
    file_prefix      : str  – prefix for saved filenames

    Returns
    -------
    list of matplotlib Figure objects
    """
    if not isinstance(data_list, list):
        data_list = [data_list]
    if routines is None:
        routines = ALL_ROUTINES
    if autosave_formats is None:
        autosave_formats = []

    figs = []

    if multisample:
        # All samples overlaid on shared plots
        if subplots:
            fig = _make_multiplot(data_list, routines, af_level)
            figs.append(fig)
            _autosave(fig, f'{file_prefix}plots', autosave_formats)
        else:
            for routine in routines:
                fig, ax = plt.subplots(figsize=(7, 5))
                _make_plot(data_list, routine, af_level, ax=ax)
                figs.append(fig)
                _autosave(fig, f'{file_prefix}{routine}', autosave_formats)

    else:
        # One set of plots per sample
        for i, data in enumerate(data_list):
            # Per-sample AF level
            if af_level is not None and hasattr(af_level, '__len__') and len(af_level) > 1:
                w_af = af_level[i]
            else:
                w_af = af_level

            if subplots:
                fig = _make_multiplot([data], routines, w_af)
                figs.append(fig)
                stem = f'{file_prefix}{data.samplename}-plots'
                _autosave(fig, stem, autosave_formats)
            else:
                for routine in routines:
                    fig, ax = plt.subplots(figsize=(7, 5))
                    _make_plot([data], routine, w_af, ax=ax)
                    figs.append(fig)
                    stem = f'{file_prefix}{data.samplename}-{routine}'
                    _autosave(fig, stem, autosave_formats)

    return figs
