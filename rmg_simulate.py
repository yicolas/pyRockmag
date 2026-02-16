"""
rmg_simulate.py
---------------
Python translations of simulation and specialised analysis functions:

  - rmg_irm_simulate()                      ← RmgIRMSimulate.m
  - rmg_simulate_switching_field_dist()     ← RmgSimulateSwitchingFieldDistribution.m
  - rmg_fractional_arm()                    ← RmgFractionalARM.m
  - rmg_arm_curve_subtract()               ← RmgARMCurveSubtract.m
  - codica_error_smoothing()               ← CODICAErrorSmoothing.m

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import butter, filtfilt

from rmg_curves import rmg_arm_curve, rmg_af_level_find


MU0 = 4 * np.pi * 1e-7


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_irm_simulate  (RmgIRMSimulate.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_irm_simulate(B_values, Ban_values, linewidth_values):
    """
    Simulate IRM acquisition using a distribution of anisotropy fields
    convolved with a Gaussian linewidth.

    Parameters
    ----------
    B_values        : 1-D array  – applied fields (T)
    Ban_values      : 1-D array  – anisotropy field values (T)
    linewidth_values: 1-D array  – Gaussian linewidth per Ban (T)

    Returns
    -------
    absorption : ndarray, shape (len(B_values), len(linewidth_values), len(Ban_values))
    """
    B_values         = np.asarray(B_values,         dtype=float).ravel()
    Ban_values       = np.asarray(Ban_values,        dtype=float).ravel()
    linewidth_values = np.asarray(linewidth_values,  dtype=float).ravel()

    dtheta = 0.001 * np.pi
    theta  = np.arange(0, 0.5 * np.pi, dtheta)

    # switching field for each (Ban, theta) pair
    Ban_g, theta_g = np.meshgrid(Ban_values, theta, indexing='ij')  # (nBan, nTheta)
    sw_field = np.abs(Ban_g / np.cos(theta_g))                       # (nBan, nTheta)
    weightings = 2 * np.pi * np.sin(theta_g) * dtheta * np.cos(theta_g)  # (nBan, nTheta)

    # For each (B, linewidth, Ban, theta):
    # contrib = weightings * N(B; sw_field, linewidth^2)
    # sum over theta → absorption
    nB  = len(B_values)
    nLW = len(linewidth_values)
    nBan = len(Ban_values)

    absorption = np.zeros((nB, nLW, nBan))

    for ib, B in enumerate(B_values):
        for ilw, lw in enumerate(linewidth_values):
            if lw <= 0:
                continue
            contrib = weightings * (1 / (np.sqrt(2 * np.pi) * lw)) * \
                      np.exp(-0.5 * ((B - sw_field) / lw) ** 2)
            absorption[ib, ilw, :] = np.sum(contrib, axis=1)

    return absorption


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_simulate_switching_field_dist  (RmgSimulateSwitchingFieldDistribution.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_simulate_switching_field_dist(B_values, Ban_values, linewidth_values):
    """
    Simulate the switching field distribution (SFD).

    Same interface as rmg_irm_simulate but uses
    sw_field = 0.5 * Ban / cos(theta - pi/2)  per the MATLAB source.

    Returns
    -------
    absorption : ndarray, shape (len(B_values), len(linewidth_values), len(Ban_values))
    """
    B_values         = np.asarray(B_values,         dtype=float).ravel()
    Ban_values       = np.asarray(Ban_values,        dtype=float).ravel()
    linewidth_values = np.asarray(linewidth_values,  dtype=float).ravel()

    dtheta = 0.01 * np.pi
    theta  = np.arange(0, 0.5 * np.pi, dtheta)

    Ban_g, theta_g = np.meshgrid(Ban_values, theta, indexing='ij')
    # Note: cos(theta - pi/2) = sin(theta)
    sw_field   = np.abs(0.5 * Ban_g / np.sin(theta_g + 1e-30))
    weightings = 2 * np.pi * np.sin(theta_g) * dtheta * np.cos(theta_g)

    nB   = len(B_values)
    nLW  = len(linewidth_values)
    nBan = len(Ban_values)
    absorption = np.zeros((nB, nLW, nBan))

    for ib, B in enumerate(B_values):
        for ilw, lw in enumerate(linewidth_values):
            if lw <= 0:
                continue
            contrib = weightings * (1 / (np.sqrt(2 * np.pi) * lw)) * \
                      np.exp(-0.5 * ((B - sw_field) / lw) ** 2)
            absorption[ib, ilw, :] = np.sum(contrib, axis=1)

    return absorption


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_fractional_arm  (RmgFractionalARM.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_fractional_arm(data, af_level=None):
    """
    Extract the fractional ARM acquisition curve at a given AF level.

    Parameters
    ----------
    data      : RmgData object
    af_level  : float (T), optional – defaults to rmg_af_level_find()

    Returns
    -------
    dict with doesExist, bias, treatmentDCField, treatmentACField,
              mvector, fracmags
    """
    if af_level is None:
        af_level = rmg_af_level_find([data])[0]

    # Find ARM steps at this AF level with zero DC bias
    arm_steps = data.stepsARM
    cand = arm_steps[
        np.isclose(data.treatmentAFFields[arm_steps], af_level, rtol=1e-6) &
        np.isclose(data.treatmentDCFields[arm_steps], 0.0,       atol=1e-12)
    ]

    if len(cand) == 0:
        return {'doesExist': False}

    stepnum = cand[-1]

    # Find subsequent contiguous ARM steps
    subseq = np.arange(stepnum + 1, len(data.levels))
    subseq_arm    = np.intersect1d(subseq, data.stepsARM)
    subseq_nonarm = np.setdiff1d(subseq, data.stepsARM)
    if len(subseq_nonarm) == 0:
        last = len(data.levels)
    else:
        last = subseq_nonarm[0]
    subseq_arm = subseq_arm[subseq_arm < last]

    # Find matching IRM step at same field
    irm_types = {'irm', 'irmz', 'irmx'}
    irm_steps = np.array([i for i, s in enumerate(data.steptypes)
                           if s.lower() in irm_types])
    cand_irm = irm_steps[
        np.isclose(data.treatmentDCFields[irm_steps], af_level, rtol=1e-6)
    ]
    if len(cand_irm) == 0:
        return {'doesExist': False}
    stepnum_irm = cand_irm[-1]

    mv_arm   = data.mvector[:, subseq_arm]
    mv_unsub = mv_arm - mv_arm[:, 0:1]
    mv_irm   = data.mvector[:, stepnum_irm] - mv_arm[:, 0]

    irm_mz = mv_irm[2]
    fracmags = np.abs(mv_unsub[2, :] / irm_mz) if irm_mz != 0 else mv_unsub[2, :] * np.nan

    return {
        'doesExist':        True,
        'bias':             data.bias[subseq_arm],
        'treatmentDCField': data.treatmentDCFields[subseq_arm],
        'treatmentACField': af_level,
        'mvectorUnsub':     data.mvector[:, subseq_arm],
        'mvector':          mv_unsub,
        'mvectorIRM':       mv_irm,
        'fracmags':         fracmags,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_arm_curve_subtract  (RmgARMCurveSubtract.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_arm_curve_subtract(curve1, curve2, weights=None):
    """
    Subtract one fractional ARM curve from another (weighted).

    Parameters
    ----------
    curve1, curve2 : dicts from rmg_fractional_arm()
    weights        : tuple (w1, w2), default (1, 1)

    Returns
    -------
    dict with doesExist, treatmentDCFields, Mz, fracmags,
              ARMsusceptibility, ARMsusceptibilityToIRM, ARMtoIRMat100uT
    """
    w1, w2 = (1.0, 1.0) if weights is None else weights

    dc1 = curve1['treatmentDCField']
    dc2 = curve2['treatmentDCField']
    common_dc = np.intersect1d(dc1, dc2)

    def _interp(dc, mz, targets):
        f = interp1d(dc, mz, bounds_error=False, fill_value='extrapolate')
        return f(targets)

    Mz1 = _interp(dc1, curve1['mvector'][2, :], common_dc) * w1
    Mz2 = _interp(dc2, curve2['mvector'][2, :], common_dc) * w2
    Mz  = Mz1 - Mz2

    mvIRM1 = curve1['mvectorIRM'] * w1
    mvIRM2 = curve2['mvectorIRM'] * w2
    mvIRM  = mvIRM1 - mvIRM2

    irm_mz = mvIRM[2]
    fracmags = np.abs(Mz / irm_mz) if irm_mz != 0 else Mz * np.nan

    dfields = np.diff(common_dc)
    deriv_fields   = common_dc[:-1] + 0.5 * dfields
    fracmags_deriv = np.diff(fracmags) / dfields

    arm_susc = float(interp1d(common_dc, np.abs(Mz),
                              bounds_error=False, fill_value='extrapolate')(1e-4)) / (1e-4 / MU0)
    arm_susc_to_irm = arm_susc / abs(irm_mz) if irm_mz != 0 else np.nan
    arm_to_irm_100  = arm_susc_to_irm * (1e-4 / MU0)

    return {
        'doesExist':              True,
        'treatmentDCFields':      common_dc,
        'treatmentACField1':      curve1.get('treatmentACField'),
        'treatmentACField2':      curve2.get('treatmentACField'),
        'mvectorIRM1':            mvIRM1,
        'mvectorIRM2':            mvIRM2,
        'Mz1':                    Mz1,
        'Mz2':                    Mz2,
        'Mz':                     Mz,
        'mvectorIRM':             mvIRM,
        'fracmags':               fracmags,
        'derivfields':            deriv_fields,
        'fracmagsderiv':          fracmags_deriv,
        'ARMsusceptibility':      arm_susc,
        'ARMsusceptibilityToIRM': arm_susc_to_irm,
        'ARMtoIRMat100uT':        arm_to_irm_100,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  codica_error_smoothing  (CODICAErrorSmoothing.m)
# ─────────────────────────────────────────────────────────────────────────────

def codica_error_smoothing(x, y, smooth_span=0, smooth_order=3):
    """
    CODICA-style error smoothing following Egli (2003).

    Fits a tanh model to the main trend, removes it to isolate residuals,
    then applies a Butterworth filter to smooth those residuals.

    Parameters
    ----------
    x, y          : 1-D array-like (IRM acquisition curve)
    smooth_span   : int, optional – Butterworth filter cutoff (default: 6% of len)
    smooth_order  : int, optional – filter order (default: 3)

    Returns
    -------
    dict with:
      residual   – y minus tanh fit
      smoothed   – Butterworth-filtered residual
      rescaled   – tanh fit + smoothed residual
      xscaling   – dict with tanh fit parameters
      errscaling – dict with sine fit to residuals
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()

    # ── tanh scaling fit ─────────────────────────────────────────────────────
    xscaling = _xscale_tanh(x, y)
    y_tanh   = xscaling['model'](x)
    residual = y - y_tanh

    # ── sine fit to residuals ────────────────────────────────────────────────
    errscaling = _xscale_sin(x, residual)

    # ── Butterworth filter ───────────────────────────────────────────────────
    if smooth_span == 0:
        smooth_span = int(np.floor(len(x) * 0.06)) * 2 + 1
    smooth_span = max(3, smooth_span)

    residual_smooth = _butter_filter(residual, smooth_span, smooth_order)
    rescaled = y_tanh + residual_smooth

    return {
        'residual':   residual,
        'smoothed':   residual_smooth,
        'rescaled':   rescaled,
        'xscaling':   xscaling,
        'errscaling': errscaling,
    }


def _xscale_tanh(x, y):
    """Fit yres + 0.5*ymax*(1 - tanh(a*(x^p - xmedian^p))) to (x,y)."""
    maxy_idx = np.argmax(y)
    miny_idx = np.argmin(y)
    xmax = x[maxy_idx]; xmin = x[miny_idx]
    miny = y[miny_idx]; maxy = y[maxy_idx]

    def tanh_model(xv, yres, ymax, a, p, xmedian):
        return yres + 0.5 * ymax * (1 - np.tanh(a * (xv ** p - xmedian ** p)))

    try:
        p0 = [miny, maxy, 5 / max((xmax - xmin) ** 0.05, 1e-10), 0.05, (xmax + xmin) / 2]
        lo = [0, 0, -np.inf, 0, -np.inf]
        hi = [np.inf, np.inf, np.inf, 1.0, np.inf]
        popt, _ = curve_fit(tanh_model, x, y, p0=p0, bounds=(lo, hi), maxfev=10000)
    except Exception:
        popt = [miny, maxy, 1.0, 0.05, np.mean(x)]

    return {
        'yres': popt[0], 'ymax': popt[1],
        'a': popt[2], 'p': popt[3], 'xmedian': popt[4],
        'model': lambda xv, p=popt: tanh_model(xv, *p),
        'popt': popt,
    }


def _xscale_sin(x, y):
    """Fit a1*sin(b1*x^p + c1) to (x,y) as error model."""
    def sin_model(xv, a1, b1, c1, p):
        return a1 * np.sin(b1 * xv ** p + c1)

    try:
        p0 = [max(y) * 0.5, 1.0, 0.0, 0.05]
        lo = [0, -np.inf, -np.inf, 0]
        hi = [np.inf, np.inf, np.inf, 1.0]
        popt, _ = curve_fit(sin_model, x, y, p0=p0, bounds=(lo, hi), maxfev=10000)
    except Exception:
        popt = [0.0, 1.0, 0.0, 0.05]

    return {
        'a1': popt[0], 'b1': popt[1], 'c1': popt[2], 'p': popt[3],
        'model': lambda xv, p=popt: sin_model(xv, *p),
        'popt': popt,
    }


def _butter_filter(data, cutoff, order):
    """
    Low-pass Butterworth filter in the frequency domain, matching
    the MATLAB FFT-based implementation in CODICAErrorSmoothing.m.
    """
    n = len(data)
    ft = np.fft.fft(data)
    v  = np.arange(1, n // 2 + 1, dtype=float)
    bfilter_half = 1 - 1.0 / ((1 + (cutoff / v) ** order) ** (0.5 * order))

    if n % 2 == 1:
        bfilter = np.concatenate([[1.0], bfilter_half, bfilter_half[::-1]])
    else:
        bfilter = np.concatenate([[1.0], bfilter_half, bfilter_half[-2::-1]])

    # pad/trim to match ft length
    bfilter = bfilter[:n]
    if len(bfilter) < n:
        bfilter = np.pad(bfilter, (0, n - len(bfilter)))

    filtered = np.fft.ifft(ft * bfilter).real
    return filtered
