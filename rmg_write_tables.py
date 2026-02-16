"""
rmg_write_tables.py
-------------------
Python translations of the table-export functions:

  - rmg_stats_write_table()    ← RmgStatsWriteTable.m
      Writes a tab-delimited .asc summary of all rock-mag statistics.

  - rmg_irm_fits_write_table() ← RmgIRMFitsWriteTable.m
      Writes a tab-delimited .asc table of Gaussian/SGG component fit parameters.

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import os
import numpy as np
from rmg_stats import rmg_stats


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_stats_write_table  (RmgStatsWriteTable.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_stats_write_table(data_list, filename):
    """
    Compute rock-mag statistics and write a tab-delimited .asc file.

    Parameters
    ----------
    data_list : list of RmgData (or single RmgData)
    filename  : str  – output filename WITHOUT extension (.asc added automatically)

    Output file columns
    -------------------
    Sample, Susceptibility, sIRM, sIRM/kg, Hcr, R,
    IRM30/IRM100, IRM100/IRM300,
    Integrated frac. IRMAcq-IRMAF,
    dfIRM/dB (field), dfIRM/dB,
    MAF of IRM, DP of IRM Acq, Skewness of IRM Acq,
    MDF of IRM, DP of IRM AF, Skewness of IRM AF,
    ARM susceptibility, ARM susceptibility/IRM,
    ARM/IRM @ 100uT, ARM/IRM @ 500uT, ARM/IRM @ 500uT - Dev from Tanh,
    MDF of ARM, DP of ARM, Skewness of ARM,
    MDF of IRM at ARM, DP of IRM at ARM, Skewness of IRM at ARM,
    NRM, NRM-directional, MDF of NRM
    """
    if not isinstance(data_list, list):
        data_list = [data_list]

    stats_list = rmg_stats(data_list)
    u = stats_list[0]['units']

    header_parts = [
        'Sample',
        f'Susceptibility ({u["susceptibility"]})',
        f'sIRM ({u["sIRM"]})',
        f'sIRM/kg ({u["sIRMperkg"]})',
        f'Hcr ({u["Hcr"]})',
        'R',
        'IRM30/IRM100',
        'IRM100/IRM300',
        f'Integrated frac. IRMAcq-IRMAF ({u["Hcr"]})',
        'dfIRM/dB (field)',
        'dfIRM/dB',
        f'MAF of IRM ({u["MDF"]})',
        'DP of IRM Acq',
        'Skewness of IRM Acq',
        f'MDF of IRM ({u["MDF"]})',
        'DP of IRM AF',
        'Skewness of IRM AF',
        f'ARM susceptibility ({u["ARMsusceptibility"]})',
        f'ARM susceptibility/IRM ({u["ARMsusceptibilityToIRM"]})',
        'ARM/IRM @ 100uT',
        'ARM/IRM @ 500uT',
        'ARM/IRM @ 500uT - Deviation from Tanh',
        f'MDF of ARM ({u["MDF"]})',
        'DP of ARM',
        'Skewness of ARM',
        f'MDF of IRM at ARM ({u["MDF"]})',
        'DP of IRM at ARM',
        'Skewness of IRM at ARM',
        f'NRM ({u["sIRM"]})',
        f'NRM-directional ({u["sIRM"]})',
        f'MDF of NRM ({u["MDF"]})',
    ]

    # Keys in the stats dict that correspond to the columns (in order)
    stat_keys = [
        'susceptibility', 'sIRM', 'sIRMperkg', 'Hcr', 'CisowskiR',
        'IRM30toIRM100', 'IRM100toIRM300', 'IntegratedIRMMinusAF',
        'dfIRMdBatField', 'dfIRMdB',
        'MAFofIRM', 'DPofIRMAcq', 'SkewnessofIRMAcq',
        'MDFofIRM', 'DPofIRMAF', 'SkewnessofIRMAF',
        'ARMsusceptibility', 'ARMsusceptibilityToIRM',
        'ARMtoIRMat100uT', 'ARMtoIRMat500uT', 'ARMtoIRMat500uTDeviationfromTanh',
        'MDFofARM', 'DPofARM', 'SkewnessofARM',
        'MDFofIRMatARM', 'DPofIRMatARM', 'SkewnessofIRMatARM',
        'NRM', 'NRMdirectional', 'MDFofNRM',
    ]

    out_path = filename if filename.endswith('.asc') else filename + '.asc'
    with open(out_path, 'w') as fid:
        fid.write('\t'.join(header_parts) + '\n')

        for s in stats_list:
            row = [s['sample']]
            for k in stat_keys:
                v = s.get(k, np.nan)
                row.append(_fmt(v))
            fid.write('\t'.join(row) + '\n')

    print(f'Stats table written → {out_path}')
    return out_path


def _fmt(v):
    """Format a value for a tab-delimited file (empty string if NaN)."""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ''
    return f'{v:.6g}'


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_irm_fits_write_table  (RmgIRMFitsWriteTable.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_irm_fits_write_table(datasets, fitsets, filename, pair_tags=None):
    """
    Write a tab-delimited .asc table of IRM component fit parameters.

    Parameters
    ----------
    datasets   : list of lists of RmgData (or single-item lists)
    fitsets    : list of lists of fit dicts (from rmg_sirm_derivative_curve_fits or similar)
                 Each fit dict must have:
                   .goodness['sse'], .goodness['adjrsquare']
                   .popt  (length 3*n_components for Gaussian: [a,b,c,...])
    filename   : str – output WITHOUT extension (.asc added automatically)
    pair_tags  : list of str, optional  – one tag per dataset group
    """
    if pair_tags is None:
        pair_tags = []

    out_path = filename if filename.endswith('.asc') else filename + '.asc'

    header = '\t'.join([
        'Tag', 'Sample', 'Component',
        'a', 'a error', 'b', 'b error', 'c', 'c error',
        'SSE', 'r²'
    ])

    with open(out_path, 'w') as fid:
        fid.write(header + '\n')

        for i, (dataset, fits) in enumerate(zip(datasets, fitsets)):
            tag = pair_tags[i] if i < len(pair_tags) else ''
            for j, (sample, fit) in enumerate(zip(dataset, fits)):
                # support both Gaussian popt (length 3k) and generic popt
                popt = fit.get('popt', np.array([]))
                good = fit.get('goodness', {})
                sse  = good.get('sse', np.nan)
                r2   = good.get('adjrsquare', np.nan)

                # Try to get 95% CI from covariance; fallback to zeros
                pcov = fit.get('pcov', None)
                if pcov is not None:
                    perr = 1.96 * np.sqrt(np.abs(np.diag(pcov)))
                else:
                    perr = np.zeros_like(popt)

                n_comps = len(popt) // 3
                for k in range(n_comps):
                    a = popt[k * 3];    a_e = perr[k * 3]
                    b = popt[k * 3 + 1]; b_e = perr[k * 3 + 1]
                    c = popt[k * 3 + 2]; c_e = perr[k * 3 + 2]
                    sname = getattr(sample, 'samplename', str(sample))
                    row = [
                        tag, sname, str(k + 1),
                        _fmt(a), _fmt(a_e),
                        _fmt(b), _fmt(b_e),
                        _fmt(c), _fmt(c_e),
                        _fmt(sse), _fmt(r2),
                    ]
                    fid.write('\t'.join(row) + '\n')

    print(f'IRM fits table written → {out_path}')
    return out_path
