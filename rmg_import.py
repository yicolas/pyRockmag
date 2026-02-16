"""
rmg_import.py
-------------
Python translation of RMGImport.m by Robert E. Kopp (2008).
Parses a 2G .rmg file and returns an RmgData object.

Observed file format (Shoemaker 2G Magnetometer):
  Row 0  : title  e.g. "    Basalt1.1,,Vol: 1, "
  Row 1  : header e.g. " ,Level,Bias Field (G),Spin Speed (rps),Hold Time (s),Mz (emu),..."
  Row 2+ : data rows — step type is the FIRST comma-delimited token; columns follow.
           Rows whose step type starts with "Instrument" or "Time:" are metadata and skipped.

Column layout (0-based AFTER the step-type token):
  0  Level            – applied field magnitude (Oe/Gauss)
  1  Bias Field (G)   – DC bias / ARM bias (Gauss)
  2  Spin Speed (rps)
  3  Hold Time (s)
  4  Mz (emu)
  5  Std Dev Z
  6  Mz/Vol           – moment per volume
  7  Moment Susceptibility (emu/Oe)
  8  Mx (emu)
  9  Std Dev X
  10 My (emu)
  ...

Unit conversions applied here:
  Fields  : Oe / Gauss → Tesla   (× 1e-4)
  Moments : emu        → Am²     (× 1e-3)

Step-type vocabulary in this file:
  NRM          natural remanent magnetisation (baseline)
  AF           AF demagnetisation of NRM
  AFmax        maximum-field baseline (before/after IRM block)
  IRMz         IRM acquisition step (z-axis, field in Level column)
  AFz          AF demagnetisation of IRM/ARM (z-axis, field in Level column)
  ARM          ARM acquisition step (AF in Level, DC bias in Bias column)

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import re
import numpy as np

MU0          = 4 * np.pi * 1e-7
OE_TO_T      = 1e-4      # Oersted / Gauss → Tesla
EMU_TO_AM2   = 1e-3      # emu → Am²


class RmgData:
    """Container matching the MATLAB Rmg structure fields."""
    pass


def _f(s):
    """Parse a string to float, return nan on failure."""
    try:
        return float(s.strip())
    except (ValueError, TypeError, AttributeError):
        return np.nan


def rmg_import(filename: str) -> RmgData:
    """
    Parse a .rmg file and return an RmgData object.

    Parameters
    ----------
    filename : str  – path to the .rmg file

    Returns
    -------
    RmgData
    """
    with open(filename, 'r', errors='replace') as fh:
        lines = [ln.rstrip('\r\n') for ln in fh.readlines()]

    # drop trailing blank lines
    while lines and not lines[-1].strip():
        lines.pop()

    y = RmgData()

    # ── Row 0: title ──────────────────────────────────────────────────────────
    y.titlerow   = lines[0]
    title_toks   = [t.strip() for t in lines[0].split(',')]
    y.samplename = title_toks[0].strip() if title_toks else 'unknown'
    y.sampledesc = title_toks[1].strip() if len(title_toks) > 1 else ''

    # Extract mass/volume from title (e.g. "Vol: 1" or "Mass: 0.152")
    y.mass = np.nan
    for tok in title_toks:
        m = re.search(r'(?:Vol|Mass|vol|mass)[:\s]+([\d.eE+\-]+)', tok)
        if m:
            val = float(m.group(1))
            if val != 1.0:          # "Vol: 1" is a placeholder, not a real mass
                y.mass = val
            break

    # ── Row 1: header ─────────────────────────────────────────────────────────
    y.headerrow   = lines[1]
    hdr           = [t.strip() for t in lines[1].split(',')]
    # hdr[0] is blank (step-type label); hdr[1] = "Level" → data index 0, etc.

    # Column index within the data portion (i.e. hdr[k] → data index k-1)
    col = {}
    for k, h in enumerate(hdr):
        if k == 0:
            continue
        di  = k - 1
        hl  = h.lower()
        if 'level' in hl:
            col['level'] = di
        elif 'bias' in hl:
            col['bias'] = di
        elif 'spin' in hl:
            col['spin'] = di
        elif hl == 'mz (emu)' or hl == 'mz':
            col.setdefault('mz', di)
        elif hl == 'mx (emu)' or hl == 'mx':
            col.setdefault('mx', di)
        elif hl == 'my (emu)' or hl == 'my':
            col.setdefault('my', di)
        elif 'mz/vol' in hl or 'mz/mass' in hl:
            col['mzperkg'] = di
        elif 'suscept' in hl:
            col['suscep'] = di
        elif 'date' in hl or 'time' in hl:
            col['datetime'] = di

    ncols = len(hdr) - 1     # number of data columns (excluding step-type)

    # ── Rows 2+: data ─────────────────────────────────────────────────────────
    # Skip metadata rows (Instrument:, Time:)
    SKIP_PREFIXES = ('instrument', 'time:')

    steptypes = []
    raw_data  = []    # list of lists: raw_data[row][col_idx]

    for line in lines[2:]:
        if not line.strip():
            continue
        parts   = [p.strip() for p in line.split(',')]
        stype   = parts[0]
        if any(stype.lower().startswith(p) for p in SKIP_PREFIXES):
            continue
        steptypes.append(stype)
        vals = parts[1:]          # everything after the step-type token
        # pad or trim to ncols
        row = [vals[i] if i < len(vals) else '' for i in range(ncols)]
        raw_data.append(row)

    n = len(steptypes)
    y.steptypes = steptypes

    def get_col(name):
        if name not in col:
            return np.full(n, np.nan)
        di = col[name]
        return np.array([_f(raw_data[r][di]) for r in range(n)], dtype=float)

    # Raw CGS values
    levels_raw  = get_col('level')     # Oe / Gauss
    bias_raw    = get_col('bias')      # Gauss
    mz_raw      = get_col('mz')        # emu
    mx_raw      = get_col('mx')        # emu
    my_raw      = get_col('my')        # emu
    mzperkg_raw = get_col('mzperkg')   # emu/cm³ or emu/g
    suscep_raw  = get_col('suscep')
    y.spin      = get_col('spin')

    # Convert to SI
    levels_T = levels_raw * OE_TO_T    # Tesla
    bias_T   = bias_raw   * OE_TO_T    # Tesla
    mz       = mz_raw     * EMU_TO_AM2 # Am²
    mx       = mx_raw     * EMU_TO_AM2
    my       = my_raw     * EMU_TO_AM2
    mzperkg  = mzperkg_raw * EMU_TO_AM2

    y.mz      = mz
    y.mx      = mx
    y.my      = my
    y.mzperkg = mzperkg
    y.suscep  = suscep_raw

    # ── Step-type index arrays ────────────────────────────────────────────────
    def match(pattern, exact=False):
        pat = pattern.lower()
        out = []
        for i, s in enumerate(steptypes):
            sl = s.lower()
            if exact:
                if sl == pat:
                    out.append(i)
            else:
                if sl.startswith(pat):
                    out.append(i)
        return np.array(out, dtype=int)

    y.stepsNRM     = match('nrm')
    y.stepsAFmax   = match('afmax')
    y.stepsAFz     = match('afz')
    y.stepsAF_nrm  = match('af', exact=True)  # pure 'AF' = NRM demagnetisation
    y.stepsAF      = np.sort(np.union1d(y.stepsAFz, y.stepsAF_nrm))
    y.stepsIRM     = match('irm')   # IRMz, IRMx, IRM
    y.stepsARM     = match('arm')
    y.stepsRRM     = match('rrm')
    y.stepsThermal = np.sort(np.union1d(match('tt'), match('trm')))

    # ── Treatment field arrays (SI) ───────────────────────────────────────────
    y.treatmentAFFields = np.zeros(n)
    y.treatmentDCFields = np.zeros(n)
    y.bias              = np.zeros(n)

    # AF demagnetisation steps: Level = the AF field applied
    for idx in np.concatenate([y.stepsAF, y.stepsAFmax]):
        y.treatmentAFFields[idx] = levels_T[idx]

    # ARM: Level = AF field, Bias = DC bias
    for idx in y.stepsARM:
        y.treatmentAFFields[idx] = levels_T[idx]
        y.treatmentDCFields[idx] = bias_T[idx]
        y.bias[idx]              = bias_T[idx]

    # IRM acquisition: Level = the DC field applied
    for idx in y.stepsIRM:
        y.treatmentDCFields[idx] = levels_T[idx]

    # RRM
    for idx in y.stepsRRM:
        y.treatmentAFFields[idx] = levels_T[idx]

    # ── mvector (3 × n) in Am² ────────────────────────────────────────────────
    y.mvector = np.vstack([mx, my, mz])    # shape (3, n)

    if np.isfinite(y.mass) and y.mass > 0:
        y.mvectorperkg = y.mvector / y.mass
    else:
        y.mvectorperkg = y.mvector * np.nan

    # ── Block detection ───────────────────────────────────────────────────────
    # Normalise variant step types into canonical block types so that
    # e.g. IRMz/IRMx are treated as the same block as IRM, and AFz as AF.
    def _canon(s):
        sl = s.lower()
        if sl in ('afz',):               return 'AF'
        if sl in ('irmz', 'irmx'):       return 'IRM'
        if sl in ('afmax',):             return 'AFMAX'
        if sl in ('nrm',):               return 'NRM'
        if sl in ('arm',):               return 'ARM'
        if sl in ('af',):                return 'AF'
        return s.upper()

    norm = [_canon(s) for s in steptypes]

    stepBlock = np.zeros(n, dtype=int)
    blockType = []
    blockSize = []

    if n > 0:
        stepBlock[0] = 1
        blockType.append(norm[0])
        blockSize.append(1)
        for i in range(1, n):
            if norm[i] != norm[i - 1]:
                new_id = stepBlock[i - 1] + 1
                stepBlock[i] = new_id
                blockType.append(norm[i])
                blockSize.append(1)
            else:
                stepBlock[i] = stepBlock[i - 1]
                blockSize[-1] += 1

    y.stepBlock = stepBlock       # 1-based block IDs (matches MATLAB convention)
    y.BlockType = blockType       # list indexed by (block_id - 1)
    y.BlockSize = np.array(blockSize, dtype=int)

    # ── Compatibility attributes expected by rmg_curves.py ──────────────────
    y.levels         = levels_raw          # raw Oe/Gauss field labels
    # y.bias already set above as per-step DC bias in Tesla
    # treatmentTemps: room temperature for all steps (no thermal data in this file)
    y.treatmentTemps = np.full(n, 298.0)
    for idx in y.stepsThermal:
        y.treatmentTemps[idx] = levels_T[idx]   # thermal steps use Level as temp
    return y

