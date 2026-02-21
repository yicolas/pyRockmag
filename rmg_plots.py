"""
rmg_plots.py
------------
Python translations of all plot functions from matRockmag:
  - rmg_sirm_plot()
  - rmg_sirm_derivative_plot()
  - rmg_arm_plot()
  - rmg_lowrie_fuller_plot()
  - rmg_fuller_plot()
  - rmg_rrm_plot()
  - rmg_backfield_plot()
  - rmg_stat_box()
  - rmg_data_full_analysis()   (the 3×3 dashboard)

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from rmg_curves import (
    rmg_sirm_curve,
    rmg_arm_curve,
    rmg_irm_backfield_curve,
    rmg_rrm_curve,
    rmg_lowrie_fuller_curves,
    rmg_fuller_curves,
    rmg_af_level_find,
    _safe_interp1,
)
from rmg_stats import rmg_stats


LINE_COLORS = ['k', 'g', 'b', 'r', 'c', 'm', 'y']
LINE_SYMS   = ['.', 'x', 'o', '+', '*', 's', 'd', 'v', '^', '>', '<']


def _style(i, colors=None, syms=None):
    colors = colors or LINE_COLORS
    syms   = syms   or LINE_SYMS
    return colors[i % len(colors)] + syms[i % len(syms)] + '-'


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_sirm_plot
# ─────────────────────────────────────────────────────────────────────────────

def rmg_sirm_plot(data_list, ax: Axes = None):
    """IRM acquisition and AF demagnetisation curves (semi-log)."""
    if not isinstance(data_list, list):
        data_list = [data_list]
    if ax is None:
        ax = plt.gca()

    irm_list = []
    for d in data_list:
        curves = rmg_sirm_curve(d)
        best = curves[-1]
        for c in curves:
            if c.doesExist and c.AF.doesExist:
                best = c
                break
        irm_list.append(best)

    for i, irm in enumerate(irm_list):
        if not irm.doesExist:
            continue
        dc   = irm.IRM.treatmentDCFields
        frac = irm.IRM.fracmags
        # skip zero-field baseline steps before semilogx
        pos = dc > 0
        if np.sum(pos) > 1:
            ax.plot(dc[pos] * 1000, frac[pos],
                        _style(i), label=data_list[i].samplename)
        if irm.AF.doesExist and hasattr(irm.AF, 'fracmags'):
            af_f = irm.AF.treatmentAFFields
            af_frac = irm.AF.fracmags
            pos_af = af_f > 0
            if np.sum(pos_af) > 1:
                ax.plot(af_f[pos_af] * 1000, af_frac[pos_af], _style(i))

    if any(irm.doesExist for irm in irm_list):
        ax.set_xlabel('B (mT)')
        ax.set_ylabel('f$_{SIRM}$')
        title_suffix = '' if len(data_list) > 1 else f': {data_list[0].samplename}'
        ax.set_title(f'IRM{title_suffix}')
        if len(data_list) > 1:
            ax.legend(loc='upper left', fontsize='small')
    else:
        ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_sirm_derivative_plot
# ─────────────────────────────────────────────────────────────────────────────

def rmg_sirm_derivative_plot(data_list, ax: Axes = None,
                              mode='both', show_points=True):
    """
    dIRM/dB log-derivative plot.
    mode: 'both', 'acq', 'af'
    """
    if not isinstance(data_list, list):
        data_list = [data_list]
    if ax is None:
        ax = plt.gca()

    irm_list = []
    for d in data_list:
        curves = rmg_sirm_curve(d)
        best = curves[-1]
        for c in curves:
            if c.doesExist and c.AF.doesExist:
                best = c
                break
        irm_list.append(best)

    IRM_COLORS = ['g', 'b', 'r']
    AF_COLORS  = ['y', 'c', 'm']
    SINGLE_COLORS = ['g', 'b', 'r', 'y', 'c', 'm']

    title_suffix = '' if len(data_list) > 1 else f': {data_list[0].samplename}'

    for i, irm in enumerate(irm_list):
        if not irm.doesExist:
            continue
        ci = i % 3
        symi = LINE_SYMS[i % len(LINE_SYMS)]

        has_af = (irm.AF.doesExist and
                  hasattr(irm.AF, 'logDerivFields') and
                  len(irm.AF.logDerivFields) > 0)

        if (mode in ('both', 'acq') and hasattr(irm.IRM, 'logDerivFields')
                and len(irm.IRM.logDerivFields) > 0):
            x_acq  = 10**(irm.IRM.logDerivFields + 3)  # Convert log to linear mT
            ld     = irm.IRM.logderiv
            lds    = irm.IRM.logderivSmooth
            finite = np.isfinite(ld)
            if show_points and np.any(finite):
                ax.plot(x_acq[finite], ld[finite], IRM_COLORS[ci] + symi)
            finite_s = np.isfinite(lds)
            if np.any(finite_s):
                ax.plot(x_acq[finite_s], lds[finite_s],
                        IRM_COLORS[ci] + '-',
                        label=f'{data_list[i].samplename} Acq')

        if (mode in ('both', 'af') and has_af
                and len(irm.AF.logDerivFields) > 0):
            x_af   = 10**(irm.AF.logDerivFields + 3)  # Convert log to linear mT
            ld     = irm.AF.logderiv
            lds    = irm.AF.logderivSmooth
            finite = np.isfinite(ld)
            if show_points and np.any(finite):
                ax.plot(x_af[finite], ld[finite], AF_COLORS[ci] + symi)
            finite_s = np.isfinite(lds)
            if np.any(finite_s):
                ax.plot(x_af[finite_s], lds[finite_s],
                        AF_COLORS[ci] + '-',
                        label=f'{data_list[i].samplename} AF')

    if any(irm.doesExist for irm in irm_list):
        ax.set_xlabel('B (mT)')  # Linear scale now
        ax.set_ylabel('df$_{IRM}$/dB')
        ax.set_title(f'dIRM/dB{title_suffix}')
        ax.set_ylim(bottom=0)
        ax.legend(loc='upper left', fontsize='small')  # Always show legend
    else:
        ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_arm_plot
# ─────────────────────────────────────────────────────────────────────────────

def rmg_arm_plot(data_list, af_levels=None, ax: Axes = None):
    """ARM acquisition curves. Plots only the last replicate if multiple exist."""
    if not isinstance(data_list, list):
        data_list = [data_list]
    if ax is None:
        ax = plt.gca()
    
    # Reference curves from RmgARMPlot.m
    # Chiton: biogenic single-domain magnetite (Kirschvink & Lowenstam 1979)
    # MS1: synthetic PSD magnetite standard
    _emuConvert   = 1e-3
    _GaussConvert = 1e-4   # Oe → T; ×1000 → mT
    
    _chiton_DC_mT  = np.arange(0, 21) * _GaussConvert * 1000
    _chiton_mz     = _emuConvert * np.array([
        1.02e-06, 1.745e-05, 3.845e-05, 5.815e-05, 7.521e-05,
        9.431e-05, 1.107e-04, 1.297e-04, 1.466e-04, 1.596e-04,
        1.751e-04, 1.916e-04, 2.064e-04, 2.214e-04, 2.379e-04,
        2.520e-04, 2.652e-04, 2.832e-04, 2.941e-04, 3.066e-04, 3.190e-04
    ])
    _chiton_IRM100 = 0.001936 * _emuConvert
    
    _MS1_DC_mT  = np.arange(0, 21) * _GaussConvert * 1000
    _MS1_mz     = _emuConvert * np.array([
        4.002e-04, 1.099e-02, 2.553e-02, 3.462e-02, 4.089e-02,
        4.560e-02, 4.913e-02, 5.138e-02, 5.335e-02, 5.458e-02,
        5.582e-02, 5.672e-02, 5.742e-02, 5.808e-02, 5.862e-02,
        5.909e-02, 5.947e-02, 5.977e-02, 6.013e-02, 6.040e-02, 6.072e-02
    ])
    _MS1_IRM100 = 0.06725 * _emuConvert
    
    # Clip to 0-1 mT to match plot x-axis
    _mask_c = _chiton_DC_mT <= 1.0
    _mask_m = _MS1_DC_mT    <= 1.0
    
    # Plot reference curves
    ax.plot(_chiton_DC_mT[_mask_c], _chiton_mz[_mask_c] / _chiton_IRM100,
            color='grey', linewidth=1.5, linestyle='--',
            zorder=1, label='Chiton (SD ref.)')
    ax.plot(_MS1_DC_mT[_mask_m], _MS1_mz[_mask_m] / _MS1_IRM100,
            color='lightgrey', linewidth=1.5, linestyle='--',
            zorder=1, label='MS1 (PSD ref.)')

    for i, d in enumerate(data_list):
        arms = rmg_arm_curve(d)
        # Filter to only existing ARM curves with fracmags
        valid_arms = [arm for arm in arms if arm.doesExist and hasattr(arm, 'fracmags')]
        
        # Plot only the LAST replicate (most recent measurement)
        if valid_arms:
            arm = valid_arms[-1]
            af_lev = arm.treatmentAFFields[0] if len(arm.treatmentAFFields) > 0 else np.nan
            lbl = f'{arm.samplename} ({round(af_lev * 1000)} mT)' if np.isfinite(af_lev) else arm.samplename
            ax.plot(arm.treatmentDCFields * 1000, arm.fracmags, _style(i), label=lbl)

    if ax.lines:
        ax.set_xlabel('B$_{DC}$ (mT)')
        ax.set_xlim(0, 1)  # X-axis 0-1 mT
        ax.set_ylim(0, 1)  # Y-axis 0-1
        title_suffix = '' if len(data_list) > 1 else f': {data_list[0].samplename}'
        ax.set_ylabel('fIRM100')  # Fractional IRM at 100 mT
        ax.set_title(f'ARM{title_suffix}')
        ax.legend(loc='lower right', fontsize='small')
    else:
        ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_lowrie_fuller_plot
# ─────────────────────────────────────────────────────────────────────────────

def rmg_lowrie_fuller_plot(data_list, af_levels=None,
                            single_mode=False, ax: Axes = None):
    """
    Lowrie-Fuller plot: normalised AF demagnetisation of ARM and IRM.
    single_mode=True plots only the ARM curve (MATLAB '1' flag).
    Plots only the last replicate if multiple measurements exist.
    """
    if not isinstance(data_list, list):
        data_list = [data_list]
    if ax is None:
        ax = plt.gca()

    NRM_COLORS = ['k', 'b', 'r']
    ARM_COLORS = ['g', 'c', 'm']

    count_arm = 0
    count_irm = 0

    for i, d in enumerate(data_list):
        curves = rmg_lowrie_fuller_curves(d)
        # Filter to existing curves
        valid_curves = [lf for lf in curves if lf.doesExist]
        
        # Plot only the LAST replicate
        if valid_curves:
            lf = valid_curves[-1]
            arm_af = lf.ARMAF
            irm_af = lf.IRMAF
            if arm_af.doesExist:
                ax.plot(arm_af.treatmentAFFields * 1000,
                        arm_af.mmagNormalized,
                        ARM_COLORS[count_arm % 3] + LINE_SYMS[count_arm % len(LINE_SYMS)] + '-',
                        label=f'{arm_af.samplename} ARM')
                count_arm += 1
            if not single_mode and irm_af.doesExist:
                ax.plot(irm_af.treatmentAFFields * 1000,
                        irm_af.mmagNormalized,
                        NRM_COLORS[count_irm % 3] + LINE_SYMS[count_irm % len(LINE_SYMS)] + '-',
                        label=f'{irm_af.samplename} IRM')
                count_irm += 1

    if ax.lines:
        ax.set_xlabel('B (mT)')
        ax.set_ylabel('Fractional Mag.')  # Changed from 'Normalised moment'
        ax.set_xlim(right=150)  # X-axis ends at 150 mT
        title_suffix = '' if len(data_list) > 1 else f': {data_list[0].samplename}'
        ax.set_title(f'Lowrie-Fuller{title_suffix}')
        ax.legend(loc='upper right', fontsize='small')
    else:
        ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_fuller_plot
# ─────────────────────────────────────────────────────────────────────────────

def rmg_fuller_plot(data_list, af_levels=None, ax: Axes = None):
    """Fuller plot: NRM/ARM vs IRM (log-log). Plots only the last ARM replicate."""
    if not isinstance(data_list, list):
        data_list = [data_list]
    if ax is None:
        ax = plt.gca()

    NRM_COLORS = ['k', 'b', 'r']
    ARM_COLORS = ['g', 'c', 'm']
    SYMS = ['.', 'o', 's', 'x']

    count_nrm = 0
    count_arm = 0

    for i, d in enumerate(data_list):
        fc = rmg_fuller_curves(d)
        if not fc.doesExist:
            continue

        if fc.NRM.doesExist and len(fc.calcNRM) > 0:
            ax.loglog(fc.calcIRM, fc.calcNRM,
                      NRM_COLORS[count_nrm % 3] + SYMS[(count_nrm * 2 + 1) % 4] + '-',
                      label=f'{fc.NRM.samplename} NRM')
            count_nrm += 1

        # Plot only the LAST ARM replicate
        valid_arms = [(arm, calc) for arm, calc in zip(fc.ARM, fc.calcARM) if arm.doesExist]
        if valid_arms:
            arm, calc_arm = valid_arms[-1]
            ax.loglog(fc.calcIRM, calc_arm,
                      ARM_COLORS[count_arm % 3] + SYMS[(count_arm * 2) % 4] + '-',
                      label=f'{arm.samplename} ARM')
            count_arm += 1

    if ax.lines:
        # Reference lines
        std = np.array([1e-15, 1])
        for frac, lbl in [(1e-4, '1:10⁴'), (1e-3, '1:1000'),
                          (1e-2, '1:100'),  (1e-1, '1:10'), (1, '1:1')]:
            ax.loglog(std, std * frac, 'k:', linewidth=0.5)

        title_suffix = '' if len(data_list) > 1 else f': {data_list[0].samplename}'
        ax.set_xlabel('IRM (Am²)')
        ax.set_ylabel('m (Am²)')
        ax.set_title(f'Fuller{title_suffix}')
        ax.legend(loc='best', fontsize='small')
    else:
        ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_rrm_plot
# ─────────────────────────────────────────────────────────────────────────────

def rmg_rrm_plot(data_list, af_levels=None, ax: Axes = None):
    """RRM plot: B_rrm vs spin rate."""
    if not isinstance(data_list, list):
        data_list = [data_list]
    if ax is None:
        ax = plt.gca()
    
    # Reference curves from RmgARMPlot.m
    # Chiton: biogenic single-domain magnetite (Kirschvink & Lowenstam 1979)
    # MS1: synthetic PSD magnetite standard
    _emuConvert   = 1e-3
    _GaussConvert = 1e-4   # Oe → T; ×1000 → mT
    
    _chiton_DC_mT  = np.arange(0, 21) * _GaussConvert * 1000
    _chiton_mz     = _emuConvert * np.array([
        1.02e-06, 1.745e-05, 3.845e-05, 5.815e-05, 7.521e-05,
        9.431e-05, 1.107e-04, 1.297e-04, 1.466e-04, 1.596e-04,
        1.751e-04, 1.916e-04, 2.064e-04, 2.214e-04, 2.379e-04,
        2.520e-04, 2.652e-04, 2.832e-04, 2.941e-04, 3.066e-04, 3.190e-04
    ])
    _chiton_IRM100 = 0.001936 * _emuConvert
    
    _MS1_DC_mT  = np.arange(0, 21) * _GaussConvert * 1000
    _MS1_mz     = _emuConvert * np.array([
        4.002e-04, 1.099e-02, 2.553e-02, 3.462e-02, 4.089e-02,
        4.560e-02, 4.913e-02, 5.138e-02, 5.335e-02, 5.458e-02,
        5.582e-02, 5.672e-02, 5.742e-02, 5.808e-02, 5.862e-02,
        5.909e-02, 5.947e-02, 5.977e-02, 6.013e-02, 6.040e-02, 6.072e-02
    ])
    _MS1_IRM100 = 0.06725 * _emuConvert
    
    # Clip to 0-1 mT to match plot x-axis
    _mask_c = _chiton_DC_mT <= 1.0
    _mask_m = _MS1_DC_mT    <= 1.0
    
    # Plot reference curves
    ax.plot(_chiton_DC_mT[_mask_c], _chiton_mz[_mask_c] / _chiton_IRM100,
            color='grey', linewidth=1.5, linestyle='--',
            zorder=1, label='Chiton (SD ref.)')
    ax.plot(_MS1_DC_mT[_mask_m], _MS1_mz[_mask_m] / _MS1_IRM100,
            color='lightgrey', linewidth=1.5, linestyle='--',
            zorder=1, label='MS1 (PSD ref.)')

    for i, d in enumerate(data_list):
        curves = rmg_rrm_curve(d)
        for rrm in curves:
            if rrm.doesExist:
                ax.plot(rrm.spins, rrm.Brrm * 1e6, _style(i),
                        label=rrm.samplename)

    if ax.lines:
        ax.set_xlabel('Spin (Hz)')
        ax.set_ylabel('B$_{RRM}$ (μT)')
        title_suffix = '' if len(data_list) > 1 else f': {data_list[0].samplename}'
        ax.set_title(f'RRM{title_suffix}')
        if len(data_list) > 1:
            ax.legend(loc='upper left', fontsize='small')
    else:
        ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_backfield_plot
# ─────────────────────────────────────────────────────────────────────────────

def rmg_backfield_plot(data_list, ax: Axes = None):
    """Backfield IRM demagnetisation curve."""
    if not isinstance(data_list, list):
        data_list = [data_list]
    if ax is None:
        ax = plt.gca()

    irm_list = []
    for d in data_list:
        curves = rmg_irm_backfield_curve(d)
        if curves and curves[-1].doesExist:
            irm_list.append((d, curves[-1]))

    for i, (d, irm) in enumerate(irm_list):
        mz = irm.Demag.mvector[2, :]
        if mz[0] < 0:
            mz = -mz
        ax.plot(irm.Demag.treatmentDCFields * 1000, mz / irm.Demag.mass,
                _style(i), label=d.samplename)

    if irm_list:
        # zero line
        dc = irm_list[-1][1].Demag.treatmentDCFields
        ax.plot(dc * 1000, np.zeros_like(dc), '-k', linewidth=0.8)
        title_suffix = '' if len(data_list) > 1 else f': {data_list[0].samplename}'
        ax.set_xlabel('-B (mT)')
        ax.set_ylabel('Am²/kg')
        ax.set_title(f'Backfield IRM{title_suffix}')
        if len(data_list) > 1:
            ax.legend(loc='best', fontsize='small')
    else:
        ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_stat_box
# ─────────────────────────────────────────────────────────────────────────────

def rmg_stat_box(data_list, af_levels=None, ax: Axes = None):
    """Text box showing key rock-mag statistics."""
    if not isinstance(data_list, list):
        data_list = [data_list]
    if ax is None:
        ax = plt.gca()

    stats_list = rmg_stats(data_list)

    lines = []
    for s in stats_list:
        lines.append(f"Sample: {s['sample']}")
        lines.append('')

        def _fmt(key, label, unit=''):
            v = s.get(key, np.nan)
            if np.isfinite(v):
                lines.append(f'{label} = {v:.4g} {unit}'.strip())

        _fmt('susceptibility',       'χ',              s['units']['susceptibility'])
        _fmt('sIRM',                 'sIRM',           s['units']['sIRM'])
        _fmt('sIRMperkg',            'sIRM/kg',        s['units']['sIRMperkg'])
        _fmt('Hcr',                  'H_cr',           s['units']['Hcr'])
        _fmt('CisowskiR',            'R')
        _fmt('dfIRMdB',              f"(df_IRM/dB)_{s.get('dfIRMdBatField',''):.4g}")
        _fmt('MDFofIRM',             'MDF_IRM',        s['units']['MDF'])
        _fmt('ARMsusceptibilityToIRM', 'k_ARM/IRM',    s['units']['ARMsusceptibilityToIRM'])
        _fmt('ARMtoIRMat100uT',      '(ARM/IRM)_0.1mT')
        _fmt('MDFofARM',             'MDF_ARM',        s['units']['MDF'])
        lines.append('')

    ax.text(0.05, 0.95, '\n'.join(lines),
            transform=ax.transAxes,
            verticalalignment='top',
            fontsize=7,
            family='monospace')
    ax.axis('off')


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_data_full_analysis  (3×3 dashboard)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_data_full_analysis(data_list, af_levels=None):
    """
    9-Panel Rock-Mag Dashboard: IRM, ARM, Backfield, Lowrie-Fuller, etc.

    Parameters
    ----------
    data_list : list of RmgData (or single RmgData)
    af_levels : list of floats, optional  (auto-detected if None)

    Returns
    -------
    matplotlib Figure
    """
    if not isinstance(data_list, list):
        data_list = [data_list]

    if af_levels is None:
        af_levels = rmg_af_level_find(data_list)
    # ensure same length
    if len(af_levels) < len(data_list):
        af_levels = (af_levels * len(data_list))[:len(data_list)]

    sample_title = data_list[0].samplename if len(data_list) == 1 else 'Multiple samples'
    fig, axes = plt.subplots(3, 3, figsize=(15, 11))
    fig.suptitle(f'8-Panel Rock-Mag Analysis – {sample_title}', fontsize=14, y=0.995)

    rmg_sirm_plot(data_list,            ax=axes[0, 0])
    rmg_sirm_derivative_plot(data_list, ax=axes[0, 1])
    rmg_arm_plot(data_list,             ax=axes[0, 2])
    rmg_lowrie_fuller_plot(data_list,   ax=axes[1, 0])
    rmg_lowrie_fuller_plot(data_list, single_mode=True, ax=axes[1, 1])
    rmg_fuller_plot(data_list,          ax=axes[1, 2])
    rmg_backfield_plot(data_list,       ax=axes[2, 0])
    rmg_stat_box(data_list,             ax=axes[2, 1])
    axes[2, 2].axis('off')  # Empty plot where RRM was

    plt.tight_layout(h_pad=3.0, w_pad=2.0)  # 0.5x more spacing
    return fig
