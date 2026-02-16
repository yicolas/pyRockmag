"""
rmg_stats.py
------------
Python translation of RmgStats.m by Robert E. Kopp (2008).
Computes rock-magnetic statistics from an RmgData object.

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import numpy as np
from rmg_curves import (
    rmg_fuller_curves,
    rmg_lowrie_fuller_curves,
    rmg_sirm_curve,
    rmg_rrm_curve,
    rmg_arm_curve,
    rmg_irm_backfield_curve,
)


def rmg_stats(data_list):
    """
    Compute rock-magnetic statistics for a list of RmgData objects.

    Parameters
    ----------
    data_list : list of RmgData  (or single RmgData)

    Returns
    -------
    list of dicts, one per sample
    """
    if not isinstance(data_list, list):
        data_list = [data_list]

    units = {
        'susceptibility':           'Am²/T',
        'sIRM':                     'Am²',
        'sIRMperkg':                'Am²/kg',
        'Hcr':                      'mT',
        'MDF':                      'mT',
        'ARMsusceptibility':        'Am²/T',
        'ARMsusceptibilityToIRM':   '1/T',
    }

    results = []
    for data in data_list:
        s = {'sample': data.samplename, 'units': units}

        # Susceptibility (last positive value)
        if hasattr(data, 'suscep') and len(data.suscep) > 0:
            pos = np.where(data.suscep > 0)[0]
            s['susceptibility'] = float(data.suscep[pos[-1]]) if len(pos) > 0 else np.nan
        else:
            s['susceptibility'] = np.nan

        # ── Lowrie-Fuller curves ─────────────────────────────────────────────
        lf_curves = rmg_lowrie_fuller_curves(data)
        lf = lf_curves[-1] if any(c.doesExist for c in lf_curves) else None

        if lf is not None and lf.doesExist:
            arm_af = lf.ARMAF
            irm_af = lf.IRMAF
            s['MDFofARM']        = getattr(arm_af, 'MDF', np.nan) * 1e3
            s['MDFofIRMatARM']   = getattr(irm_af, 'MDF', np.nan) * 1e3
            s['DPofARM']         = getattr(arm_af, 'dispersion', np.nan)
            s['DPofIRMatARM']    = getattr(irm_af, 'dispersion', np.nan)
            s['SkewnessofARM']   = getattr(arm_af, 'skewness', np.nan)
            s['SkewnessofIRMatARM'] = getattr(irm_af, 'skewness', np.nan)
        else:
            for k in ('MDFofARM', 'MDFofIRMatARM', 'DPofARM', 'DPofIRMatARM',
                      'SkewnessofARM', 'SkewnessofIRMatARM'):
                s[k] = np.nan

        # ── SIRM curve ────────────────────────────────────────────────────────
        sirm_curves = rmg_sirm_curve(data)
        # prefer one that has an AF demagnetisation
        best_sirm = sirm_curves[-1] if sirm_curves else None
        for c in sirm_curves:
            if c.doesExist and c.AF.doesExist:
                best_sirm = c
                break

        if best_sirm is not None and best_sirm.doesExist:
            sc = best_sirm
            s['sIRM']       = sc.sIRM
            s['sIRMperkg']  = sc.sIRMperkg
            s['dfIRMdBatField'] = getattr(sc, 'dfIRMdBatField', np.nan)
            s['dfIRMdB']    = getattr(sc, 'dfIRMdB', np.nan)
            s['IRM30toIRM100']  = getattr(sc, 'IRM30toIRM100', np.nan)
            s['IRM100toIRM300'] = getattr(sc, 'IRM100toIRM300', np.nan)
            s['MAFofIRM']   = getattr(sc, 'MAF', np.nan) * 1e3
            s['DPofIRMAcq'] = getattr(sc.IRM, 'dispersion', np.nan)

            if sc.AF.doesExist:
                s['Hcr']                = getattr(sc, 'Hcr', np.nan) * 1e3
                s['CisowskiR']          = getattr(sc, 'R', np.nan)
                s['SkewnessofIRMAcq']   = getattr(sc.IRM, 'skewness', np.nan)
                s['IntegratedIRMMinusAF'] = getattr(sc.diff, 'integrateddelta', np.nan)
                s['MDFofIRM']           = getattr(sc, 'MDF', np.nan) * 1e3
                s['SkewnessofIRMAF']    = getattr(sc.AF, 'skewness', np.nan)
                s['DPofIRMAF']          = getattr(sc.AF, 'dispersion', np.nan)
            else:
                for k in ('Hcr', 'CisowskiR', 'SkewnessofIRMAcq',
                          'IntegratedIRMMinusAF', 'MDFofIRM',
                          'SkewnessofIRMAF', 'DPofIRMAF'):
                    s[k] = np.nan
        else:
            for k in ('sIRM', 'sIRMperkg', 'Hcr', 'CisowskiR',
                      'IRM30toIRM100', 'IRM100toIRM300', 'IntegratedIRMMinusAF',
                      'MAFofIRM', 'MDFofIRM', 'dfIRMdBatField', 'dfIRMdB',
                      'DPofIRMAcq', 'DPofIRMAF', 'SkewnessofIRMAcq',
                      'SkewnessofIRMAF'):
                s[k] = np.nan

        # ── ARM curve ─────────────────────────────────────────────────────────
        arm_curves = rmg_arm_curve(data)
        arm = next((c for c in arm_curves if c.doesExist), None)

        if arm is not None:
            s['ARMsusceptibility']      = getattr(arm, 'ARMsusceptibility', np.nan)
            s['ARMsusceptibilityToIRM'] = getattr(arm, 'ARMsusceptibilityToIRM', np.nan)
            s['ARMtoIRMat100uT']        = getattr(arm, 'ARMtoIRMat100uT', np.nan)
            s['ARMtoIRMat500uT']        = getattr(arm, 'ARMtoIRMat500uT', np.nan)
            s['ARMtoIRMat500uTDeviationfromTanh'] = getattr(
                arm, 'ARMtoIRMat500uTDeviationfromTanh', np.nan)
        else:
            for k in ('ARMsusceptibility', 'ARMsusceptibilityToIRM',
                      'ARMtoIRMat100uT', 'ARMtoIRMat500uT',
                      'ARMtoIRMat500uTDeviationfromTanh'):
                s[k] = np.nan

        # NRM
        s['NRM'] = np.nan
        s['NRMdirectional'] = np.nan
        s['MDFofNRM'] = np.nan

        results.append(s)

    return results
