"""
fit_sgg.py
----------
Python translations of the SGG (Skewed Generalized Gaussian) fitting tools:
  - fit_sgg()        ← fitSGG.m          (free-parameter SGG fit, 1–5 components)
  - fit_sgg_comps()  ← fitSGGComps.m     (fixed-shape fit, only amplitudes free)
  - fit_sgg_multi()  ← fitSGGMulti.m     (auto-select best number of components)

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.stats import f as f_dist
from scipy.special import gamma

from sgg import sgg


# ─────────────────────────────────────────────────────────────────────────────
#  Helper: build a multi-component SGG model function
# ─────────────────────────────────────────────────────────────────────────────

def _make_sgg_model(n):
    """Return a callable f(x, *params) summing n SGG components.
    params layout: [a1, m1, s1, q1,  a2, m2, s2, q2, ...]  (4 params per comp).
    """
    def model(x, *params):
        out = np.zeros_like(x, dtype=float)
        for i in range(n):
            a, m, s, q = params[i * 4: i * 4 + 4]
            out += a * sgg(x, m, s, q, p=2)
        return out
    return model


def _make_sgg_comps_model(n, means, stds, qs):
    """Return f(x, *amplitudes) with shapes fixed, only amplitudes free."""
    def model(x, *amps):
        out = np.zeros_like(x, dtype=float)
        for i in range(n):
            out += amps[i] * sgg(x, means[i], stds[i], qs[i], p=2)
        return out
    return model


def _component_stats(x, y_comp):
    """Compute mean, dispersion (std), skewness from a component curve."""
    total = np.trapezoid(y_comp, x)
    if total <= 0 or not np.isfinite(total):
        return np.nan, np.nan, np.nan
    mu  = np.trapezoid(y_comp * x, x) / total
    var = np.trapezoid(y_comp * (x - mu) ** 2, x) / total
    sig = np.sqrt(var) if var >= 0 else np.nan
    if sig > 0:
        skw = np.trapezoid(y_comp * (x - mu) ** 3, x) / (sig ** 3 * total)
    else:
        skw = np.nan
    return mu, sig, skw


def _goodness(y_obs, y_pred):
    """Compute SSE, dfe, adj-R² matching MATLAB goodness struct."""
    n  = len(y_obs)
    sse = float(np.sum((y_obs - y_pred) ** 2))
    sst = float(np.sum((y_obs - np.mean(y_obs)) ** 2))
    r2  = 1 - sse / sst if sst > 0 else np.nan
    dfe = n - 1   # residual degrees of freedom (1 free param baseline)
    return {'sse': sse, 'dfe': dfe, 'rsquare': r2}


# ─────────────────────────────────────────────────────────────────────────────
#  fit_sgg  (fitSGG.m)
# ─────────────────────────────────────────────────────────────────────────────

def fit_sgg(x, y, n=1, start_point=None):
    """
    Fit a sum of n SGG components to (x, y).

    Parameters
    ----------
    x, y       : 1-D array-like
    n          : int (1–5)   number of components
    start_point: array-like, length 4*n, layout [a,m,s,q] per comp.
                 Auto-detected from a single Gaussian fit if None.

    Returns
    -------
    dict with keys:
      a, m, s, q, p           – component parameters (arrays of length n)
      a_error95, m_error95,   – 95% confidence half-widths
        s_error95, q_error95
      x, expected             – original x, model prediction
      y                       – per-component curves (n × len(x))
      totalArea, mean,
        dispersion, skewness  – per-component moments
      goodness                – dict with sse, dfe, rsquare
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    n = int(n)

    # ── initial guess ────────────────────────────────────────────────────────
    if start_point is None:
        # simple Gaussian estimate for location and width
        total = np.trapezoid(y, x)
        a0 = max(1e-10, total / n)
        m0 = np.average(x, weights=np.clip(y, 0, None)) if np.sum(y > 0) > 0 else np.mean(x)
        s0 = max(0.01, np.std(x) * 0.5)
        start_point = np.tile([a0, m0, s0, 1.0], n)
    else:
        start_point = np.asarray(start_point, dtype=float)

    p0 = start_point.copy()

    # ── bounds: a≥0.01, s≥0, q∈[0.5,1.5] ───────────────────────────────────
    lo = np.tile([0.01, -10.0, 1e-6, 0.5], n)
    hi = np.tile([np.inf, 10.0, np.inf, 1.5], n)

    model = _make_sgg_model(n)

    try:
        popt, pcov = curve_fit(model, x, y, p0=p0,
                               bounds=(lo, hi),
                               maxfev=10000,
                               ftol=1e-8, xtol=1e-8)
        perr_95 = 1.96 * np.sqrt(np.diag(np.abs(pcov)))
    except Exception:
        popt = p0.copy()
        perr_95 = np.full_like(p0, np.nan)

    y_pred = model(x, *popt)
    good   = _goodness(y, y_pred)
    good['dfe'] = len(x) - len(popt)   # proper dfe

    # ── unpack parameters ────────────────────────────────────────────────────
    a_arr = popt[0::4];  m_arr = popt[1::4]
    s_arr = popt[2::4];  q_arr = popt[3::4]

    a_err = perr_95[0::4];  m_err = perr_95[1::4]
    s_err = perr_95[2::4];  q_err = perr_95[3::4]

    # ── per-component curves & moments ───────────────────────────────────────
    y_comps    = np.zeros((n, len(x)))
    totalArea  = np.zeros(n)
    comp_mean  = np.zeros(n)
    comp_disp  = np.zeros(n)
    comp_skew  = np.zeros(n)

    for i in range(n):
        y_comps[i] = a_arr[i] * sgg(x, m_arr[i], s_arr[i], q_arr[i], p=2)
        totalArea[i] = np.trapezoid(y_comps[i], x)
        comp_mean[i], comp_disp[i], comp_skew[i] = _component_stats(x, y_comps[i])

    f = {
        'a': a_arr, 'm': m_arr, 's': s_arr, 'q': q_arr,
        'p': np.full(n, 2.0),
        'a_error95': a_err, 'm_error95': m_err,
        's_error95': s_err, 'q_error95': q_err,
        'p_error95': np.zeros(n),
        'x': x, 'expected': y_pred,
        'y': y_comps,
        'totalArea': totalArea,
        'mean': comp_mean, 'dispersion': comp_disp, 'skewness': comp_skew,
        'goodness': good,
        'model': model, 'popt': popt,
    }
    return f


# ─────────────────────────────────────────────────────────────────────────────
#  fit_sgg_comps  (fitSGGComps.m)
# ─────────────────────────────────────────────────────────────────────────────

def fit_sgg_comps(x, y, means, stds, qs, ps=None):
    """
    Fit amplitudes only, with shapes (means, stds, qs) fixed.

    Parameters
    ----------
    x, y         : 1-D array-like
    means, stds, qs, ps : shape parameters (length n each)

    Returns
    -------
    dict  (same structure as fit_sgg)
    """
    x  = np.asarray(x,     dtype=float).ravel()
    y  = np.asarray(y,     dtype=float).ravel()
    means = np.asarray(means, dtype=float)
    stds  = np.asarray(stds,  dtype=float)
    qs    = np.asarray(qs,    dtype=float)
    n  = len(means)
    if ps is None:
        ps = np.full(n, 2.0)

    model = _make_sgg_comps_model(n, means, stds, qs)
    total = max(1e-12, np.trapezoid(y, x))
    p0  = np.full(n, total / n)
    lo  = np.zeros(n)
    hi  = np.full(n, np.inf)

    try:
        popt, pcov = curve_fit(model, x, y, p0=p0,
                               bounds=(lo, hi),
                               maxfev=10000)
        perr_95 = 1.96 * np.sqrt(np.diag(np.abs(pcov)))
    except Exception:
        popt = p0.copy()
        perr_95 = np.full(n, np.nan)

    y_pred = model(x, *popt)
    good   = _goodness(y, y_pred)
    good['dfe'] = len(x) - n

    y_comps   = np.zeros((n, len(x)))
    totalArea = np.zeros(n)
    comp_mean = np.zeros(n)
    comp_disp = np.zeros(n)
    comp_skew = np.zeros(n)
    for i in range(n):
        y_comps[i] = popt[i] * sgg(x, means[i], stds[i], qs[i], p=2)
        totalArea[i] = np.trapezoid(y_comps[i], x)
        comp_mean[i], comp_disp[i], comp_skew[i] = _component_stats(x, y_comps[i])

    return {
        'a': popt, 'm': means, 's': stds, 'q': qs, 'p': ps,
        'a_error95': perr_95,
        'm_error95': np.zeros(n), 's_error95': np.zeros(n),
        'q_error95': np.zeros(n), 'p_error95': np.zeros(n),
        'x': x, 'expected': y_pred,
        'y': y_comps,
        'totalArea': totalArea,
        'mean': comp_mean, 'dispersion': comp_disp, 'skewness': comp_skew,
        'goodness': good,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  fit_sgg_multi  (fitSGGMulti.m)
# ─────────────────────────────────────────────────────────────────────────────

def fit_sgg_multi(x, y, max_n=None):
    """
    Auto-select the best number of SGG components using an F-ratio test (p>0.01).

    Parameters
    ----------
    x, y   : 1-D array-like
    max_n  : int, optional  maximum number of components to try (default: min((len-1)//5, 5))

    Returns
    -------
    dict with all fit_sgg fields plus:
      fitSGG_list  – list of fit dicts for each n tried
      Fratio, pvalue – arrays of F-ratio and p-value for each additional component
      bestfitSGGn  – chosen number of components
      bestfit      – the chosen fit dict
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()

    if max_n is None:
        max_n = min((len(x) - 1) // 5, 5)
    max_n = max(1, int(max_n))

    fits   = []
    Fratios = []
    pvalues = []

    # fit n=1
    fits.append(fit_sgg(x, y, n=1))
    best_n = 1

    for n in range(2, max_n + 1):
        # warm-start from previous fit
        prev = fits[-1]
        # interleave prev params and add a new component at the centre
        sp_prev = np.zeros((n - 1) * 4)
        sp_prev[0::4] = prev['a'] * (n - 1) / n
        sp_prev[1::4] = prev['m']
        sp_prev[2::4] = prev['s'] * (n - 1) / n
        sp_prev[3::4] = np.ones(n - 1)
        new_comp = [np.mean(prev['a']), np.mean(prev['m']),
                    np.mean(prev['s']), 1.0]
        sp = np.concatenate([sp_prev, new_comp])

        fits.append(fit_sgg(x, y, n=n, start_point=sp))

        prev_g = fits[-2]['goodness']
        curr_g = fits[-1]['goodness']

        sse_prev, dfe_prev = prev_g['sse'], prev_g['dfe']
        sse_curr, dfe_curr = curr_g['sse'], curr_g['dfe']

        if sse_curr > 0 and dfe_curr > 0 and dfe_prev > dfe_curr:
            Fr = ((sse_prev - sse_curr) / sse_curr) / \
                 ((dfe_prev - dfe_curr) / dfe_curr)
            pv = 1 - f_dist.cdf(Fr, dfe_curr, dfe_prev)
        else:
            Fr = np.nan
            pv = 1.0   # not significant → stop

        Fratios.append(Fr)
        pvalues.append(pv)

        if pv > 0.01:
            best_n = n - 1
            break
        else:
            best_n = n

    best_fit = fits[best_n - 1]

    result = dict(best_fit)   # copy all best-fit fields to top level
    result.update({
        'fitSGG_list': fits,
        'Fratio': np.array(Fratios),
        'pvalue': np.array(pvalues),
        'bestfitSGGn': best_n,
        'bestfit': best_fit,
    })
    return result
