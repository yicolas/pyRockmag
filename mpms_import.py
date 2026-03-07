"""
mpms_import.py
--------------
Python port of MPMSImport.m (Robert E. Kopp, 2008)

Imports a Quantum Design MPMS low-temperature magnetometry file (.dat/.csv)
and computes FC/ZFC/LTC curves plus key derived parameters.

The MPMS .dat format has 31 header lines, then comma-delimited data where:
  col 3 (index 2) = applied field (Oe)   → converted to T  (/1e4)
  col 4 (index 3) = temperature (K)
  col 5 (index 4) = moment (emu)         → converted to A m²  (/1e3)

Usage
-----
    from mpms_import import mpms_import
    d = mpms_import('sample.dat')
    print(d.sample, d.dFC, d.memory)
"""

import os
import numpy as np


class MPMSData:
    """Container for MPMS import results. All attributes mirror the
    MATLAB struct fields produced by MPMSImport.m."""
    pass


def _interpolate(x, y, xi):
    """1-D linear interpolation at a single point xi, ignoring NaNs."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return np.nan
    return float(np.interp(xi, x[mask], y[mask]))


def mpms_import(filename: str) -> MPMSData:
    """
    Parse an MPMS .dat file and return an MPMSData object.

    Parameters
    ----------
    filename : str
        Path to the MPMS .dat (or .csv) file.

    Returns
    -------
    MPMSData
        Object with attributes:
          sample            – sample name (filename stem)
          data              – raw 2-D array (all columns)
          field             – applied field in T
          temperature       – temperature in K
          moment            – moment in A m²
          RTPoints          – indices where T ≈ 300 K
          FCTemp / FCMoment – field-cooled curve
          FCdMoment         – dM/dT for FC
          ZFCTemp / ZFCMoment – zero-field-cooled curve
          ZFCdMoment        – dM/dT for ZFC
          LTCTemp / LTCMoment – low-temperature cycling curve
          FCMomentNorm      – FC moment normalised to value at 80 K
          ZFCMomentNorm     – ZFC moment normalised to FC value at 80 K
          LTCMomentNorm     – LTC moment normalised to starting value
          LTCReversalPoint  – index of temperature minimum in LTC
          LTCCoolingdMoment – dM/dT on LTC cooling leg
          LTCWarmingdMoment – dM/dT on LTC warming leg
          dFC               – (FC80-FC150)/FC80
          dZFC              – (ZFC80-ZFC150)/ZFC80
          dFCtodZFC         – dFC / dZFC
          memory            – LTC end moment / LTC start moment
    """
    d = MPMSData()

    # ── sample name: everything before the last '.' ───────────────────────
    base = os.path.basename(filename)
    dot  = base.rfind('.')
    d.sample = base[:dot] if dot > 0 else base

    # ── read data: 31 header lines, comma-delimited ───────────────────────
    rows = []
    with open(filename, 'r') as fh:
        for _ in range(31):
            fh.readline()
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            try:
                row = [float(p) if p.strip() != '' else np.nan
                       for p in parts]
                rows.append(row)
            except ValueError:
                continue

    d.data = np.array(rows)

    # ── extract primary columns ───────────────────────────────────────────
    d.field       = d.data[:, 2] / 1e4          # Oe  → T
    d.temperature = d.data[:, 3]                 # K
    d.moment      = d.data[:, 4] / 1e3          # emu → A m²

    # ── room-temperature waypoints (T ≈ 300 K) ───────────────────────────
    d.RTPoints = np.where(np.abs(d.temperature - 300) < 0.5)[0]

    if len(d.RTPoints) < 2:
        raise ValueError(
            f"mpms_import: fewer than 2 room-temperature points found in "
            f"'{filename}'. Check that the file contains FC + ZFC + LTC segments."
        )

    # ── FC segment: index 0 → first RT point ─────────────────────────────
    d.FCTemp   = d.temperature[: d.RTPoints[0] + 1]
    d.FCMoment = d.moment     [: d.RTPoints[0] + 1]
    d.FCdMoment = np.diff(d.FCMoment) / np.diff(d.FCTemp)

    # ── ZFC segment: first RT point+1 → second RT point ──────────────────
    d.ZFCTemp   = d.temperature[d.RTPoints[0] + 1 : d.RTPoints[1] + 1]
    d.ZFCMoment = d.moment     [d.RTPoints[0] + 1 : d.RTPoints[1] + 1]
    d.ZFCdMoment = np.diff(d.ZFCMoment) / np.diff(d.ZFCTemp)

    # ── LTC segment: second-to-last RT point → last RT point ─────────────
    d.LTCTemp   = d.temperature[d.RTPoints[-2] : d.RTPoints[-1] + 1]
    d.LTCMoment = d.moment     [d.RTPoints[-2] : d.RTPoints[-1] + 1]

    # ── interpolated reference values ────────────────────────────────────
    FC80  = _interpolate(d.FCTemp,  d.FCMoment,  80.0)
    FC150 = _interpolate(d.FCTemp,  d.FCMoment,  150.0)
    ZFC80  = _interpolate(d.ZFCTemp, d.ZFCMoment, 80.0)
    ZFC150 = _interpolate(d.ZFCTemp, d.ZFCMoment, 150.0)

    # ── normalised curves ─────────────────────────────────────────────────
    d.FCMomentNorm  = d.FCMoment  / FC80  if np.isfinite(FC80)  and FC80  != 0 else d.FCMoment * np.nan
    d.ZFCMomentNorm = d.ZFCMoment / FC80  if np.isfinite(FC80)  and FC80  != 0 else d.ZFCMoment * np.nan
    d.LTCMomentNorm = (d.LTCMoment / d.LTCMoment[0]
                       if d.LTCMoment[0] != 0 else d.LTCMoment * np.nan)

    # ── LTC reversal (minimum temperature) ───────────────────────────────
    rev_idx = int(np.argmin(d.LTCTemp))
    d.LTCReversalPoint = rev_idx

    d.LTCCoolingdMoment = (np.diff(d.LTCMoment[: rev_idx + 1]) /
                            np.diff(d.LTCTemp  [: rev_idx + 1]))
    d.LTCWarmingdMoment = (np.diff(d.LTCMoment[rev_idx:]) /
                            np.diff(d.LTCTemp  [rev_idx:]))

    # ── scalar derived parameters ─────────────────────────────────────────
    d.dFC     = (FC80 - FC150) / FC80    if np.isfinite(FC80)  and FC80  != 0 else np.nan
    d.dZFC    = (ZFC80 - ZFC150) / ZFC80 if np.isfinite(ZFC80) and ZFC80 != 0 else np.nan
    d.dFCtodZFC = d.dFC / d.dZFC         if d.dZFC != 0 else np.nan
    d.memory  = (d.LTCMoment[-1] / d.LTCMoment[0]
                 if d.LTCMoment[0] != 0 else np.nan)

    return d
