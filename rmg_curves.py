"""
rmg_curves.py
-------------
Python translations of the curve-extraction MATLAB functions:
  - moving()
  - rmg_af_level_find()
  - rmg_extract_af_of_step()
  - rmg_sirm_curve()
  - rmg_arm_curve()
  - rmg_irm_backfield_curve()
  - rmg_rrm_curve()
  - rmg_lowrie_fuller_curves()
  - rmg_fuller_curves()

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import numpy as np
from scipy.interpolate import interp1d


MU0 = 4 * np.pi * 1e-7


# ─────────────────────────────────────────────────────────────────────────────
#  Utility helpers
# ─────────────────────────────────────────────────────────────────────────────

def moving(y, span):
    """
    Moving average matching MATLAB's moving() / smooth() behaviour.
    Handles edge effects with triangular weighting at ends.
    """
    y = np.asarray(y, dtype=float).ravel()
    span = int(np.floor(span))
    n = len(y)
    span = min(span, n)
    width = span - 1 + (span % 2)   # force odd
    if width <= 1:
        return y.copy()

    # main filter
    kernel = np.ones(width) / width
    c = np.convolve(y, kernel, mode='full')[:n]

    # beginning edge
    cbegin_vals = []
    for i in range(1, width - 1, 2):
        cbegin_vals.append(np.sum(y[:i + 1]) / (i + 1))
    cbegin = np.array(cbegin_vals)

    # end edge
    cend_vals = []
    for i in range(1, width - 1, 2):
        cend_vals.append(np.sum(y[n - i - 1:]) / (i + 1))
    cend = np.array(cend_vals[::-1])

    # combine
    mid = c[width - 1:]
    result = np.concatenate([cbegin, mid, cend])
    # safety: match length
    if len(result) != n:
        result = np.interp(np.arange(n), np.linspace(0, n - 1, len(result)), result)
    return result


def _safe_interp1(x, y, xi, kind='linear', fill_value='extrapolate'):
    """Wrapper around scipy interp1d that handles edge cases."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    valid = np.isfinite(x) & np.isfinite(y)
    if np.sum(valid) < 2:
        return np.full_like(xi, np.nan, dtype=float)
    try:
        f = interp1d(x[valid], y[valid], kind=kind,
                     bounds_error=False, fill_value='extrapolate')
        return f(xi)
    except Exception:
        return np.full_like(xi, np.nan, dtype=float)


# ─────────────────────────────────────────────────────────────────────────────
#  _compute_af_stats  (shared logic for AF and IRM acquisition stats)
# ─────────────────────────────────────────────────────────────────────────────

def _compute_af_stats(fields, normalized):
    """
    Given treatment fields and a normalized magnetisation curve, compute
    MDF, dispersion, skewness, log-derivative etc.
    Returns a dict of stats (added as attributes elsewhere).
    """
    stats = {}
    fields = np.asarray(fields, dtype=float)
    normalized = np.asarray(normalized, dtype=float)

    if len(fields) < 2:
        return stats

    # filter out zero / negative fields before taking log10
    pos_mask = fields > 0
    if np.sum(pos_mask) < 2:
        return stats
    fields     = fields[pos_mask]
    normalized = normalized[pos_mask]

    # log derivatives
    dfields     = np.diff(fields)
    dlogfields  = np.diff(np.log10(fields))
    # guard against any remaining zero dlogfields (duplicate field values)
    valid = dlogfields != 0
    if np.sum(valid) < 1:
        return stats
    logDerivFields = np.log10(fields[:-1]) + 0.5 * dlogfields
    with np.errstate(invalid='ignore', divide='ignore'):
        logderiv = np.where(valid, -np.diff(normalized) / dlogfields, np.nan)

    smoothSpan = max(1, int(np.floor(len(logDerivFields) * 0.03)) * 2 + 1)
    logderivSmooth = moving(logderiv, smoothSpan)

    stats['logDerivFields']  = logDerivFields
    stats['logderiv']        = logderiv
    stats['logderivSmooth']  = logderivSmooth

    # MDF
    try:
        mdf = _safe_interp1(normalized[::-1], fields[::-1], 0.5)[0]
    except Exception:
        mdf = np.nan
    stats['MDF'] = mdf

    # mean, dispersion, skewness (log-normal moments)
    w = logderivSmooth.copy()
    w[w < 0] = 0
    total = np.trapezoid(w, logDerivFields)
    if total > 0 and np.isfinite(total):
        mean_log = np.trapezoid(w * logDerivFields, logDerivFields) / total
        disp = np.sqrt(np.trapezoid(w * (logDerivFields - mean_log) ** 2, logDerivFields) / total)
        skew = np.trapezoid(w * (logDerivFields - mean_log) ** 3, logDerivFields) / (disp ** 3 * total) if disp > 0 else np.nan
    else:
        mean_log = disp = skew = np.nan

    stats['meanlogfield'] = mean_log
    stats['dispersion']   = disp
    stats['skewness']     = skew
    return stats


# ─────────────────────────────────────────────────────────────────────────────
#  RmgAFLevelFind
# ─────────────────────────────────────────────────────────────────────────────

def rmg_af_level_find(data_list, mode=None):
    """
    Find the AF treatment level for each RmgData in data_list.
    Returns a list of AF levels (in Tesla).
    """
    if not isinstance(data_list, list):
        data_list = [data_list]

    af_levels = []
    for data in data_list:
        af_level = 1000e-4   # default 100 mT in Tesla
        if mode == 'RRM' and len(data.stepsRRM) > 0:
            levels = data.treatmentAFFields[data.stepsRRM]
            af_level = levels[-1]
        if len(data.stepsARM) > 0:
            levels = data.treatmentAFFields[data.stepsARM]
            af_level = levels[-1]
        af_levels.append(af_level)
    return af_levels


# ─────────────────────────────────────────────────────────────────────────────
#  RmgExtractAFOfStep
# ─────────────────────────────────────────────────────────────────────────────

def rmg_extract_af_of_step(data, stepnum, subtract_afmax=True):
    """
    Extract the AF demagnetisation sequence following step index stepnum (0-based).
    Returns an object with doesExist, mvector, treatmentAFFields, etc.
    """
    class AF:
        pass
    y = AF()
    y.experiment = 'AF'
    y.samplename = data.samplename
    y.mass       = data.mass
    y.doesExist  = False

    ts = y.targetStep = type('TS', (), {})()
    ts.stepnum          = stepnum
    ts.steptype         = data.steptypes[stepnum]
    ts.level            = data.levels[stepnum]
    ts.bias             = data.bias[stepnum]
    ts.spin             = data.spin[stepnum]
    ts.mvector          = data.mvector[:, stepnum]
    ts.treatmentDCField = data.treatmentDCFields[stepnum]
    ts.treatmentAFField = data.treatmentAFFields[stepnum]
    ts.treatmentTemp    = data.treatmentTemps[stepnum]

    n = len(data.levels)
    if stepnum + 1 >= n:
        return y

    currentBlock = data.stepBlock[stepnum]   # 1-based

    # Next block must be 'AF'
    next_block_idx = int(currentBlock)       # block IDs are 1-based, list is 0-based
    if next_block_idx >= len(data.BlockType):
        return y
    if data.BlockType[next_block_idx].upper() != 'AF':
        return y

    # The next step must be in the next block
    if data.stepBlock[stepnum + 1] == currentBlock:
        return y

    # Find baseline (AFMax block adjacent to current block)
    prev_block_idx = currentBlock - 2   # 0-based index for block (currentBlock-1)
    next_next_idx  = currentBlock + 1   # 0-based for block (currentBlock+2)

    magBaselineStep = -1
    if prev_block_idx >= 0 and prev_block_idx < len(data.BlockType):
        if data.BlockType[prev_block_idx].upper() == 'AFMAX':
            steps_in_prev = np.where(data.stepBlock == (currentBlock - 1))[0]
            if len(steps_in_prev) > 0:
                magBaselineStep = steps_in_prev[-1]
    if magBaselineStep < 0 and next_next_idx < len(data.BlockType):
        if data.BlockType[next_next_idx].upper() == 'AFMAX':
            steps_in_nn = np.where(data.stepBlock == (currentBlock + 2))[0]
            if len(steps_in_nn) > 0:
                magBaselineStep = steps_in_nn[0]

    if magBaselineStep >= 0:
        y.baselineVector = data.mvector[:, magBaselineStep]
    else:
        y.baselineVector = np.zeros(3)

    y.doesExist = True
    af_block_id = currentBlock + 1
    y.AFSteps = np.where(data.stepBlock == af_block_id)[0]

    y.mvectorUnsub = data.mvector[:, y.AFSteps]

    if subtract_afmax:
        y.mvector = y.mvectorUnsub - y.baselineVector[:, np.newaxis]
    else:
        y.mvector = y.mvectorUnsub.copy()

    y.mmag  = np.sqrt(np.sum(y.mvector ** 2, axis=0))
    dmv = np.diff(y.mvector, axis=1)
    y.dmmag = np.sqrt(np.sum(dmv ** 2, axis=0))
    c = np.cumsum(y.dmmag[::-1])
    y.mdirectional = np.append(c[::-1], 0) + y.mmag[-1]

    y.mz      = y.mvector[2, :]
    y.mzNormalized          = y.mz / y.mz[0] if y.mz[0] != 0 else y.mz * np.nan
    y.mmagNormalized        = y.mmag / y.mmag[0] if y.mmag[0] != 0 else y.mmag * np.nan
    y.mdirectionalNormalized = y.mdirectional / y.mdirectional[0] if y.mdirectional[0] != 0 else y.mdirectional * np.nan

    y.treatmentAFFields  = data.treatmentAFFields[y.AFSteps]
    y.IRMDemagStepCount  = len(y.treatmentAFFields)

    if y.IRMDemagStepCount > 1:
        stats = _compute_af_stats(y.treatmentAFFields, y.mzNormalized)
        for k, v in stats.items():
            setattr(y, k, v)

    return y


# ─────────────────────────────────────────────────────────────────────────────
#  RmgSIRMCurve
# ─────────────────────────────────────────────────────────────────────────────

def rmg_sirm_curve(data):
    """
    Returns a list of SIRM curve objects for all IRM acquisition blocks.
    """
    class SIRMCurve:
        experiment = 'SIRM'

    def _empty():
        y = SIRMCurve()
        y.doesExist = False
        y.AF = type('AF', (), {'doesExist': False})()
        return y

    # candidate blocks: IRM/IRMz/IRMx with BlockSize > 1
    irm_types = {'irm', 'irmz', 'irmx'}
    candidate_blocks = [
        bid for bid, bt in enumerate(data.BlockType, start=1)
        if bt.lower() in irm_types and data.BlockSize[bid - 1] > 1
    ]

    if not candidate_blocks:
        return [_empty()]

    results = []
    for bid in candidate_blocks:
        blocksteps = np.where(data.stepBlock == bid)[0]
        if len(blocksteps) == 0:
            continue

        y = SIRMCurve()
        y.doesExist  = True
        y.samplename = data.samplename
        y.mass       = data.mass

        y.AF = rmg_extract_af_of_step(data, blocksteps[-1])
        y.baselineVector = y.AF.baselineVector if y.AF.doesExist else np.zeros(3)

        y.IRM = type('IRM', (), {})()
        y.IRM.IRMSteps         = blocksteps
        y.IRM.treatmentDCFields = data.treatmentDCFields[blocksteps]
        y.IRM.mvectorUnsub     = data.mvector[:, blocksteps]
        y.IRM.mvector          = y.IRM.mvectorUnsub - y.baselineVector[:, np.newaxis]
        y.IRM.IRMperkg         = data.mzperkg[blocksteps]

        sIRM_end = y.IRM.mvector[2, -1]
        y.IRM.fracmags = np.abs(y.IRM.mvector[2, :] / sIRM_end) if sIRM_end != 0 else np.zeros(len(blocksteps))
        y.sIRM       = np.abs(sIRM_end)
        y.sIRMperkg  = y.sIRM / data.mass if data.mass > 0 else np.nan

        if len(blocksteps) > 1:
            dc = y.IRM.treatmentDCFields
            # filter zero / negative DC fields before log10
            pos = dc > 0
            dc_pos = dc[pos]
            frac_pos = y.IRM.fracmags[pos]
            dfields    = np.diff(dc)
            y.IRM.dfields     = dfields
            y.IRM.derivFields = dc[:-1] + 0.5 * dfields
            if len(dc_pos) > 1:
                dlogfields = np.diff(np.log10(dc_pos))
                valid = dlogfields != 0
                y.IRM.dlogfields     = dlogfields
                y.IRM.logDerivFields = np.log10(dc_pos[:-1]) + 0.5 * dlogfields
                with np.errstate(invalid='ignore', divide='ignore'):
                    y.IRM.derivFracmags  = np.diff(y.IRM.fracmags[pos]) / np.diff(dc_pos)
                    y.IRM.logderiv       = np.where(valid,
                                                    np.diff(frac_pos) / dlogfields,
                                                    np.nan)
                smoothSpan = max(1, int(np.floor(len(y.IRM.logDerivFields) * 0.03)) * 2 + 1)
                y.IRM.logderivSmooth = moving(y.IRM.logderiv, smoothSpan)
            else:
                y.IRM.dlogfields     = np.array([])
                y.IRM.logDerivFields = np.array([])
                y.IRM.derivFracmags  = np.array([])
                y.IRM.logderiv       = np.array([])
                y.IRM.logderivSmooth = np.array([])

            # MDF, MAF, dispersion, skewness
            try:
                y.MAF = _safe_interp1(y.IRM.fracmags, dc, 0.5)[0]
            except Exception:
                y.MAF = np.nan

            stats = _compute_af_stats(dc, y.IRM.fracmags)
            y.IRM.dispersion = stats.get('dispersion', np.nan)
            y.IRM.skewness   = stats.get('skewness', np.nan)
            y.IRM.MDF        = stats.get('MDF', np.nan)

        # Hcr, Cisowski R from backfield overlap
        if y.AF.doesExist:
            irm_fields  = y.IRM.treatmentDCFields
            irm_fracs   = y.IRM.fracmags
            af_fields   = y.AF.treatmentAFFields
            af_fracs    = y.AF.fracmags if hasattr(y.AF, 'fracmags') else np.abs(y.AF.mzNormalized)

            # Build af fracmag relative to sIRM acquisition
            af_mz = y.AF.mvector[2, :]
            af_fracmags = np.abs(af_mz / sIRM_end) if sIRM_end != 0 else af_mz * np.nan
            y.AF.fracmags = af_fracmags

            try:
                y.Hcr = _safe_interp1(af_fracmags, af_fields, 0.5)[0]
            except Exception:
                y.Hcr = np.nan

            # Cisowski R: fraction of IRM fields where IRM_acq > IRM_AF (normalized to same max)
            try:
                common_fields = np.sort(np.unique(np.concatenate([irm_fields, af_fields])))
                irm_interp = _safe_interp1(irm_fields, irm_fracs, common_fields)
                af_interp  = _safe_interp1(af_fields,  af_fracmags, common_fields)
                y.R = float(np.sum(irm_interp < af_interp)) / len(common_fields)
            except Exception:
                y.R = np.nan

            # MDF of AF demagnetisation
            if y.AF.IRMDemagStepCount > 1:
                af_stats = _compute_af_stats(af_fields, af_fracmags)
                y.MDF        = af_stats.get('MDF', np.nan)
                y.AF.dispersion = af_stats.get('dispersion', np.nan)
                y.AF.skewness   = af_stats.get('skewness', np.nan)
                y.AF.logDerivFields = af_stats.get('logDerivFields', np.array([]))
                y.AF.logderiv       = af_stats.get('logderiv', np.array([]))
                y.AF.logderivSmooth = af_stats.get('logderivSmooth', np.array([]))
            else:
                y.MDF = np.nan

            # Integrated fractional IRM acq – AF
            try:
                common_f = np.sort(np.unique(np.concatenate([irm_fields, af_fields])))
                irm_i = _safe_interp1(irm_fields, irm_fracs, common_f)
                af_i  = _safe_interp1(af_fields,  af_fracmags, common_f)
                y.diff = type('D', (), {})()
                y.diff.integrateddelta = np.trapezoid(irm_i - af_i, common_f)
            except Exception:
                y.diff = type('D', (), {'integrateddelta': np.nan})()

        # IRM30/IRM100, IRM100/IRM300
        dc = y.IRM.treatmentDCFields
        fracs = y.IRM.fracmags
        try:
            irm30  = _safe_interp1(dc, fracs, 30e-3)[0]
            irm100 = _safe_interp1(dc, fracs, 100e-3)[0]
            irm300 = _safe_interp1(dc, fracs, 300e-3)[0]
            y.IRM30toIRM100  = irm30 / irm100 if irm100 != 0 else np.nan
            y.IRM100toIRM300 = irm100 / irm300 if irm300 != 0 else np.nan
            y.IRM.logderiv  # ensure exists
            y.dfIRMdBatField = dc[-1] * 1000
            y.dfIRMdB        = y.IRM.logderiv[-1]
        except Exception:
            y.IRM30toIRM100  = np.nan
            y.IRM100toIRM300 = np.nan

        results.append(y)

    return results if results else [_empty()]


# ─────────────────────────────────────────────────────────────────────────────
#  RmgARMCurve
# ─────────────────────────────────────────────────────────────────────────────

def _corresponding_irm(data, irm_level):
    """Find the IRM mvector at a given DC field level."""
    irm_types = {'irm', 'irmz', 'irmx'}
    all_irm_steps = [i for i, s in enumerate(data.steptypes) if s.lower() in irm_types]
    if not all_irm_steps:
        return np.array([np.nan, np.nan, np.nan])

    all_irm_steps = np.array(all_irm_steps)
    exact = all_irm_steps[np.where(data.treatmentDCFields[all_irm_steps] == irm_level)[0]]
    if len(exact) > 0:
        return data.mvector[:, exact[-1]]

    # Interpolate
    fd = data.treatmentDCFields[all_irm_steps] - irm_level
    sign_changes = np.where(fd[:-1] * fd[1:] < 0)[0]
    if len(sign_changes) == 0:
        return np.array([np.nan, np.nan, np.nan])

    idx = sign_changes[0]
    s0, s1 = all_irm_steps[idx], all_irm_steps[idx + 1]
    x0, x1 = np.log10(data.treatmentDCFields[s0]), np.log10(data.treatmentDCFields[s1])
    xi = np.log10(irm_level)
    result = np.zeros(3)
    for row in range(3):
        result[row] = np.interp(xi, [x0, x1], [data.mvector[row, s0], data.mvector[row, s1]])
    return result


def rmg_arm_curve(data):
    """Returns a list of ARM curve objects."""
    class ARMCurve:
        experiment = 'ARM'

    def _empty():
        y = ARMCurve()
        y.doesExist = False
        return y

    # ARM blocks followed immediately by AF block
    candidate_blocks = []
    for bid in range(1, len(data.BlockType)):
        if (data.BlockType[bid - 1].upper() == 'ARM' and
                bid < len(data.BlockType) and
                data.BlockType[bid].upper() == 'AF'):
            candidate_blocks.append(bid)   # 1-based block id

    if not candidate_blocks:
        return [_empty()]

    results = []
    for bid in candidate_blocks:
        blocksteps = np.where(data.stepBlock == bid)[0]
        if len(blocksteps) == 0:
            continue

        y = ARMCurve()
        y.doesExist   = True
        y.samplename  = data.samplename
        y.mass        = data.mass

        y.ARMsteps          = blocksteps
        y.bias              = data.bias[blocksteps]
        y.treatmentDCFields = data.treatmentDCFields[blocksteps]
        y.treatmentAFFields = data.treatmentAFFields[blocksteps]
        y.mvectorUnsub      = data.mvector[:, blocksteps]
        y.mvector           = y.mvectorUnsub - y.mvectorUnsub[:, 0:1]  # subtract first step

        denom = y.mvector[2, :]
        dc100 = interp1d(y.treatmentDCFields, np.abs(y.mvector[2, :]),
                         bounds_error=False, fill_value='extrapolate')(1e-4)
        y.ARMsusceptibility = dc100 / (1e-4 / MU0)

        unique_af = np.unique(y.treatmentAFFields)
        af_level = unique_af[0] if len(unique_af) == 1 else np.nan

        if np.isfinite(af_level):
            irm_mv = _corresponding_irm(data, af_level)
            y.mvectorIRM = irm_mv - y.mvectorUnsub[:, 0]
            irm_mz = y.mvectorIRM[2]
        else:
            irm_mz = np.nan

        if np.isfinite(irm_mz) and irm_mz != 0:
            y.fracmags = np.abs(y.mvector[2, :] / irm_mz)
            y.ARMsusceptibilityToIRM = y.ARMsusceptibility / np.abs(irm_mz)
            y.ARMtoIRMat100uT  = y.ARMsusceptibilityToIRM * (1e-4 / MU0)
            y.ARMtoIRMat500uT  = float(interp1d(
                y.treatmentDCFields, np.abs(y.mvector[2, :]),
                bounds_error=False, fill_value='extrapolate')(5e-4)) / np.abs(irm_mz)
            y.ARMtoIRMat500uTDeviationfromTanh = (
                y.ARMtoIRMat500uT - np.tanh((5e-4 / MU0) * y.ARMsusceptibilityToIRM))
        else:
            # no matching IRM — plot raw ARM Mz normalised to its own maximum
            mz_abs = np.abs(y.mvector[2, :])
            mz_max = mz_abs[-1] if mz_abs[-1] != 0 else (np.max(mz_abs) if np.max(mz_abs) > 0 else 1.0)
            y.fracmags = mz_abs / mz_max
            y.mvectorIRM = np.array([np.nan, np.nan, np.nan])
            y.ARMsusceptibilityToIRM = np.nan
            y.ARMtoIRMat100uT = np.nan
            y.ARMtoIRMat500uT = np.nan
            y.ARMtoIRMat500uTDeviationfromTanh = np.nan

        dfields = np.diff(y.treatmentDCFields)
        y.derivfields   = y.treatmentDCFields[:-1] + 0.5 * dfields
        with np.errstate(invalid='ignore', divide='ignore'):
            y.fracmagsderiv = np.diff(y.fracmags) / dfields

        results.append(y)

    return results if results else [_empty()]


# ─────────────────────────────────────────────────────────────────────────────
#  RmgIRMBackfieldCurve
# ─────────────────────────────────────────────────────────────────────────────

def rmg_irm_backfield_curve(data):
    """Returns a list of IRM backfield curve objects."""
    class BackfieldCurve:
        experiment = 'IRMBackfield'

    def _empty():
        y = BackfieldCurve()
        y.doesExist = False
        y.samplename = data.samplename
        return y

    irm_sequences = rmg_sirm_curve(data)
    if not any(s.doesExist for s in irm_sequences):
        return [_empty()]

    # Find sign changes in DC field within each sequence
    candidates = {}
    for i, seq in enumerate(irm_sequences):
        if not seq.doesExist:
            continue
        dc = seq.IRM.treatmentDCFields
        dsign = np.diff(np.sign(dc))
        cands = np.where(np.abs(dsign) == 2)[0] + 1  # step AFTER sign change
        if len(cands) > 0:
            candidates[i] = cands

    if not candidates:
        return [_empty()]

    results = []
    for i, cand_idxs in candidates.items():
        seq = irm_sequences[i]
        for ci in cand_idxs:
            # These are positions within the IRM block; map to global step
            global_step = seq.IRM.IRMSteps[ci]
            demag = _extract_dc_demag(data, global_step)
            if demag.doesExist:
                obj = BackfieldCurve()
                obj.doesExist  = True
                obj.samplename = data.samplename
                obj.Demag      = demag
                results.append(obj)

    return results if results else [_empty()]


def _extract_dc_demag(data, stepnum):
    """Extract DC demagnetisation steps following a given step."""
    class Demag:
        pass
    y = Demag()
    y.samplename = data.samplename
    y.mass       = data.mass
    y.doesExist  = False

    n = len(data.levels)
    if stepnum + 1 >= n:
        return y

    ts = type('TS', (), {})()
    ts.stepnum          = stepnum
    ts.steptype         = data.steptypes[stepnum]
    ts.level            = data.levels[stepnum]
    ts.bias             = data.bias[stepnum]
    ts.spin             = data.spin[stepnum]
    ts.mvector          = data.mvector[:, stepnum]
    ts.treatmentDCField = data.treatmentDCFields[stepnum]
    ts.treatmentTemp    = data.treatmentTemps[stepnum]
    y.targetStep = ts

    irm_types = {'irm', 'irmz', 'irmx'}
    next_type = data.steptypes[stepnum + 1].lower()
    if next_type not in irm_types:
        return y

    currentBlock = data.stepBlock[stepnum]
    demagBlock   = data.stepBlock[stepnum + 1]

    # Find baseline from adjacent AFMax block
    magBaselineStep = -1
    prev_bi = currentBlock - 2
    nnxt_bi = demagBlock
    if prev_bi >= 0 and data.BlockType[prev_bi].upper() == 'AFMAX':
        steps_prev = np.where(data.stepBlock == (currentBlock - 1))[0]
        if len(steps_prev) > 0:
            magBaselineStep = steps_prev[-1]
    elif nnxt_bi < len(data.BlockType) and data.BlockType[nnxt_bi].upper() == 'AFMAX':
        steps_nn = np.where(data.stepBlock == (demagBlock + 1))[0]
        if len(steps_nn) > 0:
            magBaselineStep = steps_nn[0]

    if magBaselineStep >= 0:
        y.baselineVector = data.mvector[:, magBaselineStep]
    else:
        y.baselineVector = np.zeros(3)

    y.doesExist  = True
    demag_steps  = np.intersect1d(
        np.where(data.stepBlock == demagBlock)[0],
        np.arange(stepnum, n)
    )
    y.DemagSteps  = demag_steps
    y.mvectorUnsub = data.mvector[:, demag_steps]
    y.mvector      = y.mvectorUnsub - y.baselineVector[:, np.newaxis]
    y.treatmentDCFields = data.treatmentDCFields[demag_steps]
    y.mass = data.mass
    return y


# ─────────────────────────────────────────────────────────────────────────────
#  RmgRRMCurve
# ─────────────────────────────────────────────────────────────────────────────

def rmg_rrm_curve(data):
    """Returns a list of RRM curve objects."""
    class RRMCurve:
        experiment = 'RRM'

    def _empty():
        y = RRMCurve()
        y.doesExist = False
        return y

    candidate_blocks = [bid for bid, bt in enumerate(data.BlockType, start=1)
                        if bt.upper() == 'RRM']
    if not candidate_blocks:
        return [_empty()]

    # ARM info for susceptibility
    arm_curves = rmg_arm_curve(data)
    arm_fields = []
    for arm in arm_curves:
        if arm.doesExist:
            uq = np.unique(arm.treatmentAFFields)
            arm_fields.append(uq[0] if len(uq) == 1 else np.nan)
        else:
            arm_fields.append(np.nan)

    results = []
    for bid in candidate_blocks:
        blocksteps = np.where(data.stepBlock == bid)[0]
        if len(blocksteps) == 0:
            continue

        y = RRMCurve()
        y.doesExist  = True
        y.samplename = data.samplename
        y.mass       = data.mass
        y.RRMsteps   = blocksteps
        y.spins      = data.spin[blocksteps]

        # baseline at spin=0 step
        zero_idx = np.where(y.spins == 0)[0]
        if len(zero_idx) == 0:
            y.doesExist = False
            continue
        y.baselineVector = data.mvector[:, blocksteps[zero_idx[0]]]
        y.deltas = data.mvector[:, blocksteps] - y.baselineVector[:, np.newaxis]

        # Find matching ARM
        target_af = data.treatmentAFFields[blocksteps[0]]
        arm_susc = np.nan
        for j, af in enumerate(arm_fields):
            if np.isclose(af, target_af, rtol=1e-6):
                arm_susc = arm_curves[j].ARMsusceptibility
                break

        y.ARMsusceptibility = arm_susc
        if arm_susc != 0 and np.isfinite(arm_susc):
            y.Brrm = MU0 * y.deltas[2, :] / arm_susc
        else:
            y.Brrm = np.full(len(blocksteps), np.nan)

        results.append(y)

    return results if results else [_empty()]


# ─────────────────────────────────────────────────────────────────────────────
#  RmgLowrieFullerCurves
# ─────────────────────────────────────────────────────────────────────────────

def rmg_lowrie_fuller_curves(data):
    """Returns a list of Lowrie-Fuller curve objects."""
    class LFCurve:
        experiment = 'Lowrie-Fuller'

    def _empty():
        y = LFCurve()
        y.doesExist = False
        return y

    irm_types = {'irm', 'irmz', 'irmx'}

    # Candidate ARM blocks: ARM followed by AF
    cand_arm = [bid for bid in range(1, len(data.BlockType))
                if data.BlockType[bid - 1].upper() == 'ARM'
                and data.BlockType[bid].upper() == 'AF']

    # Candidate IRM blocks: IRM followed by AF
    cand_irm = [bid for bid in range(1, len(data.BlockType))
                if data.BlockType[bid - 1].lower() in irm_types
                and data.BlockType[bid].upper() == 'AF']

    if not cand_arm or not cand_irm:
        return [_empty()]

    # IRM levels
    irm_steps_list = {}
    irm_levels     = {}
    for bid in cand_irm:
        steps = np.where(data.stepBlock == bid)[0]
        if len(steps) > 0:
            irm_steps_list[bid] = steps
            irm_levels[bid]     = abs(data.treatmentDCFields[steps[-1]])

    # ARM levels & match to IRM
    results = []
    for arm_bid in cand_arm:
        arm_steps = np.where(data.stepBlock == arm_bid)[0]
        if len(arm_steps) == 0:
            continue
        arm_af_level = data.treatmentAFFields[arm_steps[-1]]

        # Find matching IRM block (IRM level == ARM AF level)
        match_irm = None
        for irm_bid, irm_lev in irm_levels.items():
            if np.isclose(irm_lev, arm_af_level, rtol=1e-6):
                match_irm = irm_bid
                break
        if match_irm is None:
            continue

        y = LFCurve()
        y.doesExist = True
        y.ARMAF  = rmg_extract_af_of_step(data, arm_steps[-1])
        y.IRMAF  = rmg_extract_af_of_step(data, irm_steps_list[match_irm][-1])

        if y.ARMAF.doesExist and y.IRMAF.doesExist:
            arm_mdf = getattr(y.ARMAF, 'MDF', np.nan)
            irm_mdf = getattr(y.IRMAF, 'MDF', np.nan)
            y.deltaMDF = arm_mdf - irm_mdf
        else:
            y.deltaMDF = np.nan

        results.append(y)

    return results if results else [_empty()]


# ─────────────────────────────────────────────────────────────────────────────
#  RmgFullerCurves
# ─────────────────────────────────────────────────────────────────────────────

def rmg_fuller_curves(data):
    """Returns a Fuller Curves object."""
    class FullerCurve:
        experiment = 'Fuller'

    y = FullerCurve()
    y.doesExist = False

    irm_types = {'irm', 'irmz', 'irmx'}

    # NRM blocks followed by AF
    cand_nrm = [bid for bid in range(1, len(data.BlockType))
                if data.BlockType[bid - 1].upper() == 'NRM'
                and data.BlockType[bid].upper() == 'AF']
    if not cand_nrm:
        return y

    nrm_steps = np.where(data.stepBlock == cand_nrm[0])[0]
    y.NRM = rmg_extract_af_of_step(data, nrm_steps[-1])

    # IRM blocks followed by AF
    cand_irm = [bid for bid in range(1, len(data.BlockType))
                if data.BlockType[bid - 1].lower() in irm_types
                and data.BlockType[bid].upper() == 'AF']
    if not cand_irm:
        return y

    irm_steps = np.where(data.stepBlock == cand_irm[-1])[0]
    y.IRM = rmg_extract_af_of_step(data, irm_steps[-1])

    if not y.IRM.doesExist:
        return y

    y.doesExist   = True
    y.trialFields = y.IRM.treatmentAFFields
    y.calcIRM     = y.IRM.mmag

    if y.NRM.doesExist:
        y.calcNRM = np.abs(_safe_interp1(y.NRM.treatmentAFFields, y.NRM.mmag, y.trialFields))
        y.NRMtoIRM = y.calcNRM[0] / y.calcIRM[0] if y.calcIRM[0] != 0 else np.nan
    else:
        y.calcNRM  = np.zeros_like(y.trialFields)
        y.NRMtoIRM = np.nan

    # ARM blocks followed by AF
    cand_arm = [bid for bid in range(1, len(data.BlockType))
                if data.BlockType[bid - 1].upper() == 'ARM'
                and data.BlockType[bid].upper() == 'AF']

    y.ARM     = []
    y.calcARM = []
    for arm_bid in cand_arm:
        arm_steps = np.where(data.stepBlock == arm_bid)[0]
        if len(arm_steps) > 0:
            arm_af = rmg_extract_af_of_step(data, arm_steps[-1])
            y.ARM.append(arm_af)
            if arm_af.doesExist:
                interped = np.abs(_safe_interp1(arm_af.treatmentAFFields, arm_af.mmag, y.trialFields))
                y.calcARM.append(interped)
            else:
                y.calcARM.append(np.zeros_like(y.trialFields))

    if y.calcARM:
        y.calcARM = np.array(y.calcARM)
        y.ARMtoIRM = (y.calcARM[:, 0] / y.calcIRM[0]
                      if y.calcIRM[0] != 0 else np.full(len(y.ARM), np.nan))
    else:
        y.calcARM  = np.zeros((0, len(y.trialFields)))
        y.ARMtoIRM = np.array([])

    return y
