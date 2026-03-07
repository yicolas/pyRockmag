"""
rmg_curve_fits.py
-----------------
Python translations of the IRM/AF derivative curve-fitting functions:

  - rmg_sirm_derivative_curve_fits()      ← RmgSIRMDerivativeCurveFits.m
      Fits 1–5 Gaussian components to IRM acquisition and AF demagnetisation
      log-derivative curves using F-ratio model selection.

  - rmg_sirm_derivative_curve_fits_sgg()  ← RmgSIRMDerivativeCurveFitsSGG.m
      Same but uses SGG (Skewed Generalized Gaussian) components via fit_sgg_multi.

  - rmg_sirm_derivative_curve_fit_comps() ← RmgSIRMDerivativeCurveFitComps.m
      Fits only amplitudes with user-supplied means and stds (Gaussian).

  - rmg_af_derivative_curve_fits()        ← RmgAFDerivativeCurveFits.m
      Fits 1–5 Gaussians to an AF demagnetisation log-derivative curve.

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import f as f_dist
from scipy.stats import norm

from fit_sgg import fit_sgg_multi


# ─────────────────────────────────────────────────────────────────────────────
#  Internal helpers: multi-Gaussian fitting with F-ratio model selection
# ─────────────────────────────────────────────────────────────────────────────

def _gauss_n(n):
    """Return a sum-of-n-Gaussians model: params = [a1,b1,c1, a2,b2,c2, ...]."""
    def model(x, *p):
        out = np.zeros_like(x, dtype=float)
        for i in range(n):
            a, b, c = p[i * 3], p[i * 3 + 1], p[i * 3 + 2]
            out += a * np.exp(-0.5 * ((x - b) / c) ** 2)
        return out
    return model


def _fit_gauss_n(x, y, n, p0=None):
    """Fit sum of n Gaussians to (x,y). Returns (popt, goodness_dict)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    model = _gauss_n(n)

    if p0 is None:
        a0 = max(y) / n
        b0 = np.mean(x)
        c0 = np.std(x) * 0.5 or 0.1
        p0 = np.tile([a0, b0, c0], n)

    lo = np.tile([0.0, -np.inf, 1e-9], n)
    hi = np.tile([np.inf, np.inf, np.inf], n)

    try:
        popt, pcov = curve_fit(model, x, y, p0=p0,
                               bounds=(lo, hi),
                               maxfev=10000)
    except Exception:
        popt = np.array(p0, dtype=float)
        pcov = np.diag(np.full(len(p0), np.nan))

    y_pred = model(x, *popt)
    sse = float(np.sum((y - y_pred) ** 2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    dfe = max(1, len(x) - len(popt))
    r2  = 1 - sse / sst if sst > 0 else np.nan

    good = {'sse': sse, 'dfe': dfe, 'rsquare': r2,
            'adjrsquare': 1 - (1 - r2) * (len(x) - 1) / dfe if r2 is not np.nan else np.nan}
    return popt, good, model


def _auto_gauss_fit(x, y, max_n=5):
    """
    Try Gaussian fits from 1 to max_n components, stop when F-ratio p > 0.01.
    Returns a list of fit dicts and the index of the best one.
    """
    fits = []
    best_n = 1

    popt1, good1, mdl1 = _fit_gauss_n(x, y, 1)
    fits.append({'popt': popt1, 'goodness': good1, 'model': mdl1, 'n': 1,
                 'fit': lambda xv, p=popt1, m=mdl1: m(xv, *p)})

    for n in range(2, max_n + 1):
        # warm-start: carry over previous params, add a new component
        prev_p = fits[-1]['popt']
        new_p  = np.concatenate([prev_p,
                                 [max(y) / n, np.mean(x), np.std(x) * 0.5]])
        popt_n, good_n, mdl_n = _fit_gauss_n(x, y, n, p0=new_p)
        fits.append({'popt': popt_n, 'goodness': good_n, 'model': mdl_n, 'n': n,
                     'fit': lambda xv, p=popt_n, m=mdl_n: m(xv, *p)})

        prev_g = fits[-2]['goodness']
        curr_g = fits[-1]['goodness']
        sse_p, dfe_p = prev_g['sse'], prev_g['dfe']
        sse_c, dfe_c = curr_g['sse'], curr_g['dfe']

        if sse_c > 0 and dfe_c > 0 and dfe_p > dfe_c:
            Fr = ((sse_p - sse_c) / sse_c) / ((dfe_p - dfe_c) / dfe_c)
            pv = 1 - f_dist.cdf(Fr, dfe_c, dfe_p)
        else:
            Fr = np.nan
            pv = 1.0

        fits[-1]['Fratio'] = Fr
        fits[-1]['pvalue'] = pv

        if pv > 0.01:
            best_n = n - 1
            break
        else:
            best_n = n

    return fits, best_n


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_sirm_derivative_curve_fits  (RmgSIRMDerivativeCurveFits.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_sirm_derivative_curve_fits(sirm_curve):
    """
    Fit multi-Gaussian models to the log-derivative of an SIRM curve.

    Parameters
    ----------
    sirm_curve : object returned by rmg_sirm_curve() – must have
                 .IRM.logDerivFields, .IRM.logderivSmooth,
                 .AF.logDerivFields,  .AF.logderivSmooth

    Returns
    -------
    The input object with added attributes:
      .IRM.logderivfit   – list of fit dicts (index 0 = 1-component, etc.)
      .IRM.logderivbestfit, .IRM.logderivbestfitgaussians
      .AF.logderivfit
      .AF.logderivbestfit, .AF.logderivbestfitgaussians
    """
    import copy
    y = copy.copy(sirm_curve)

    # ── IRM acquisition ───────────────────────────────────────────────────────
    irm_x = np.asarray(sirm_curve.IRM.logDerivFields, dtype=float)
    irm_y = np.asarray(sirm_curve.IRM.logderivSmooth, dtype=float).ravel()
    irm_fits, irm_best_n = _auto_gauss_fit(irm_x, irm_y)
    y.IRM.logderivfit           = irm_fits
    y.IRM.logderivbestfit       = irm_fits[irm_best_n - 1]
    y.IRM.logderivbestfitgaussians = irm_best_n

    # ── AF demagnetisation ────────────────────────────────────────────────────
    if (sirm_curve.AF.doesExist and
            hasattr(sirm_curve.AF, 'logDerivFields') and
            len(sirm_curve.AF.logDerivFields) > 1):
        af_x = np.asarray(sirm_curve.AF.logDerivFields, dtype=float)
        af_y = np.asarray(sirm_curve.AF.logderivSmooth, dtype=float).ravel()
        af_fits, af_best_n = _auto_gauss_fit(af_x, af_y)
        y.AF.logderivfit           = af_fits
        y.AF.logderivbestfit       = af_fits[af_best_n - 1]
        y.AF.logderivbestfitgaussians = af_best_n

    return y


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_sirm_derivative_curve_fits_sgg  (RmgSIRMDerivativeCurveFitsSGG.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_sirm_derivative_curve_fits_sgg(sirm_curve):
    """
    Fit multi-SGG models (auto-selected) to IRM and AF log-derivative curves.

    Returns a dict with keys 'IRMfit' and 'AFfit', each being the output of
    fit_sgg_multi().
    """
    irm_x = np.asarray(sirm_curve.IRM.logDerivFields, dtype=float)
    irm_y = np.asarray(sirm_curve.IRM.logderivSmooth, dtype=float).ravel()
    irm_fit = fit_sgg_multi(irm_x, irm_y)

    af_fit = None
    if (sirm_curve.AF.doesExist and
            hasattr(sirm_curve.AF, 'logDerivFields') and
            len(sirm_curve.AF.logDerivFields) > 1):
        af_x = np.asarray(sirm_curve.AF.logDerivFields, dtype=float)
        af_y = np.asarray(sirm_curve.AF.logderivSmooth, dtype=float).ravel()
        af_fit = fit_sgg_multi(af_x, af_y)

    return {'IRMfit': irm_fit, 'AFfit': af_fit}


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_sirm_derivative_curve_fit_comps  (RmgSIRMDerivativeCurveFitComps.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_sirm_derivative_curve_fit_comps(sirm_curve, means, stds):
    """
    Fit amplitudes only to IRM and AF log-derivative curves, with normal
    (Gaussian) components whose means and stds are fixed.

    Parameters
    ----------
    sirm_curve : SIRM curve object
    means, stds : array-like of length n  (in log-field units, i.e. log10(mT))

    Returns
    -------
    dict with keys 'AF' and 'IRM', each containing 'fit', 'goodness', 'popt'
    """
    means = np.asarray(means, dtype=float)
    stds  = np.asarray(stds,  dtype=float)
    n     = len(means)

    def _norm_model_factory(mu, si):
        def model(x, *amps):
            out = np.zeros_like(x, dtype=float)
            for i in range(n):
                out += amps[i] * norm.pdf(x, mu[i], si[i])
            return out
        return model

    model = _norm_model_factory(means, stds)
    p0  = np.full(n, 1.0)
    lo  = np.zeros(n)
    hi  = np.full(n, np.inf)

    result = {'means': means, 'stds': stds}

    for tag, fields_attr, smooth_attr, obj in [
            ('IRM', 'logDerivFields', 'logderivSmooth', sirm_curve.IRM),
            ('AF',  'logDerivFields', 'logderivSmooth', sirm_curve.AF),
    ]:
        if tag == 'AF' and (not sirm_curve.AF.doesExist or
                             not hasattr(sirm_curve.AF, 'logDerivFields')):
            result['AF'] = {'popt': np.full(n, np.nan), 'goodness': {}}
            continue

        xv = np.asarray(getattr(obj, fields_attr), dtype=float)
        yv = np.asarray(getattr(obj, smooth_attr),  dtype=float).ravel()

        try:
            popt, _ = curve_fit(model, xv, yv, p0=p0,
                                bounds=(lo, hi), maxfev=10000)
        except Exception:
            popt = p0.copy()

        y_pred = model(xv, *popt)
        sse = float(np.sum((yv - y_pred) ** 2))
        sst = float(np.sum((yv - np.mean(yv)) ** 2))
        dfe = max(1, len(xv) - n)
        r2  = 1 - sse / sst if sst > 0 else np.nan

        result[tag] = {
            'popt':      popt,
            'fit':       lambda xv2, p=popt: model(xv2, *p),
            'goodness':  {'sse': sse, 'dfe': dfe, 'adjrsquare': r2},
        }

    return result


# ─────────────────────────────────────────────────────────────────────────────
#  rmg_af_derivative_curve_fits  (RmgAFDerivativeCurveFits.m)
# ─────────────────────────────────────────────────────────────────────────────

def rmg_af_derivative_curve_fits(af_curve):
    """
    Fit 1–5 Gaussians to the mz log-derivative of an AF demagnetisation curve.

    Parameters
    ----------
    af_curve : object from rmg_extract_af_of_step() – must have
               .log10treatmentAFderivFields, .mzlogderivSmooth

    Returns
    -------
    The input object (copy) with added attributes:
      .mzlogderivfit           – list of fit dicts
      .mzlogderivbestfit
      .mzlogderivbestfitgaussians
    """
    import copy
    y = copy.copy(af_curve)

    if not hasattr(af_curve, 'log10treatmentAFderivFields'):
        return y

    xv = np.asarray(af_curve.log10treatmentAFderivFields, dtype=float)
    yv = np.asarray(af_curve.mzlogderivSmooth, dtype=float).ravel()

    fits, best_n = _auto_gauss_fit(xv, yv)
    y.mzlogderivfit               = fits
    y.mzlogderivbestfit           = fits[best_n - 1]
    y.mzlogderivbestfitgaussians  = best_n
    return y


# =============================================================================
#  NEW: fit_sgg_comps                         (fitSGGComps.m)
#  NEW: rmg_sirm_derivative_curve_fit_comps_sgg  (component-constrained SGG)
# =============================================================================

def fit_sgg_comps(x, y, means, stds, qs, ps=None):
    """
    Fit a weighted sum of SGG components where shapes are fixed and only
    amplitudes (a_i) are free parameters.

    Port of fitSGGComps.m (Kopp 2008).

    Parameters
    ----------
    x, y   : array-like — data to fit
    means  : array-like — fixed component means
    stds   : array-like — fixed component standard deviations
    qs     : array-like — fixed SGG skewness parameters
    ps     : array-like, optional — SGG shape exponents (default 2 = Gaussian)

    Returns
    -------
    dict with: a, m, s, q, p  (a fitted, rest fixed)
               a_error95, m_error95, s_error95, q_error95, p_error95
               x, expected
               y  (n_points × n_components array of component curves)
               totalArea, mean, dispersion, skewness  (per component)
               goodness  (sse, rsquare, dfe)
    """
    _trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz

    x     = np.asarray(x,     dtype=float).ravel()
    y     = np.asarray(y,     dtype=float).ravel()
    means = np.asarray(means, dtype=float)
    stds  = np.asarray(stds,  dtype=float)
    qs    = np.asarray(qs,    dtype=float)
    n     = len(means)
    ps    = np.full(n, 2.0) if ps is None else np.asarray(ps, dtype=float)

    try:
        from sgg import sgg as _sgg
    except ImportError:
        def _sgg(xd, m, s, q, p=2):
            z = (xd - m) / max(abs(s), 1e-30)
            return np.exp(-0.5 * np.abs(z) ** p)

    def _model(xdata, *amplitudes):
        out = np.zeros_like(xdata, dtype=float)
        for i, a in enumerate(amplitudes):
            out += a * _sgg(xdata, means[i], stds[i], qs[i], ps[i])
        return out

    start = np.full(n, max(float(_trapz(y, x)) / n, 1e-30))

    try:
        popt, pcov = curve_fit(_model, x, y, p0=start,
                               bounds=(np.zeros(n), np.full(n, np.inf)),
                               maxfev=10_000)
    except Exception:
        popt = start.copy()
        pcov = np.full((n, n), np.nan)

    diag = np.diag(pcov)
    perr = (np.sqrt(np.where(diag >= 0, diag, np.nan))
            if not np.any(np.isnan(pcov)) else np.full(n, np.nan))

    expected = _model(x, *popt)

    y_comps  = np.zeros((len(x), n))
    totalArea  = np.zeros(n)
    comp_mean  = np.zeros(n)
    comp_disp  = np.zeros(n)
    comp_skew  = np.zeros(n)

    for i in range(n):
        comp = popt[i] * _sgg(x, means[i], stds[i], qs[i], ps[i])
        y_comps[:, i] = comp
        area = float(_trapz(comp, x))
        totalArea[i] = area
        if area > 0:
            comp_mean[i] = float(_trapz(comp * x, x)) / area
            comp_disp[i] = float(np.sqrt(_trapz(comp * (x - comp_mean[i])**2, x) / area))
            if comp_disp[i] > 0:
                comp_skew[i] = float(_trapz(comp * (x - comp_mean[i])**3, x) /
                                     (area * comp_disp[i]**3))

    residuals = y - expected
    sse = float(np.sum(residuals**2))
    sst = float(np.sum((y - np.mean(y))**2))
    dfe = max(len(x) - n, 1)

    return {
        'a':          popt,
        'm':          means.copy(),
        's':          stds.copy(),
        'q':          qs.copy(),
        'p':          ps.copy(),
        'a_error95':  1.96 * perr,
        'm_error95':  np.zeros(n),
        's_error95':  np.zeros(n),
        'q_error95':  np.zeros(n),
        'p_error95':  np.zeros(n),
        'x':          x,
        'expected':   expected,
        'y':          y_comps,
        'totalArea':  totalArea,
        'mean':       comp_mean,
        'dispersion': comp_disp,
        'skewness':   comp_skew,
        'goodness':   {
            'sse':     sse,
            'rsquare': 1.0 - sse / sst if sst > 0 else np.nan,
            'dfe':     dfe,
        },
    }


def rmg_sirm_derivative_curve_fit_comps_sgg(sirm_curve, means, stds, qs=None):
    """
    Fit IRM and AF derivative curves with component shapes fixed externally
    using SGG components (amplitudes only are free).

    Port of RmgSIRMDerivativeCurveFitComps.m — extended to support SGG shapes.
    Pass qs=None for symmetric Gaussian behaviour identical to MATLAB original.

    Parameters
    ----------
    sirm_curve : SIRM curve object from rmg_sirm_curve()
    means      : array-like — component means in log10(mT)
    stds       : array-like — component standard deviations
    qs         : array-like, optional — SGG skewness (default: zeros)

    Returns
    -------
    dict with keys 'IRM' and 'AF', each containing a fit_sgg_comps result dict,
    plus 'means', 'stds', 'qs'
    """
    import copy

    means = np.asarray(means, dtype=float)
    stds  = np.asarray(stds,  dtype=float)
    qs    = np.zeros(len(means)) if qs is None else np.asarray(qs, dtype=float)

    result = {'means': means, 'stds': stds, 'qs': qs}

    if hasattr(sirm_curve, 'AF') and getattr(sirm_curve.AF, 'doesExist', False):
        af_x = np.asarray(sirm_curve.AF.logDerivFields,  dtype=float).ravel()
        af_y = np.asarray(sirm_curve.AF.logderivSmooth,  dtype=float).ravel()
        result['AF'] = fit_sgg_comps(af_x, af_y, means, stds, qs)
    else:
        result['AF'] = None

    if hasattr(sirm_curve, 'IRM') and len(getattr(sirm_curve.IRM, 'logDerivFields', [])) > 0:
        irm_x = np.asarray(sirm_curve.IRM.logDerivFields, dtype=float).ravel()
        irm_y = np.asarray(sirm_curve.IRM.logderivSmooth, dtype=float).ravel()
        result['IRM'] = fit_sgg_comps(irm_x, irm_y, means, stds, qs)
    else:
        result['IRM'] = None

    return result
