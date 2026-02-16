"""
sgg.py
------
Python translation of SGG.m by Robert E. Kopp (2008).

The SGG (Skewed Generalized Gaussian) probability density function,
used for fitting IRM component distributions.

Copyright (C) 2008 Robert E. Kopp (original MATLAB code), GNU GPL v.3
Python port 2025.
"""

import numpy as np


def sgg(x, m, s, q, p=2):
    """
    Skewed Generalized Gaussian (SGG) probability density function.

    Parameters
    ----------
    x : array-like   – evaluation points
    m : float        – location (mean in log-field space)
    s : float        – scale (width)
    q : float        – skewness parameter  (q' in [0,2], mapped to q in [-1,1])
    p : float        – shape parameter (p=2 → generalized Gaussian)

    Returns
    -------
    y : ndarray, shape (len(x),)

    Notes
    -----
    qprime is used such that (1,2) maps to (-1,0):
        q = mod(qprime + 1, 2) - 1
    """
    x = np.asarray(x, dtype=float).ravel()

    # map qprime → q
    qv = (q + 1) % 2 - 1

    # avoid exact zeros (numerical stability)
    EPS = 1e-17
    sv = s   + (s   == 0) * EPS
    qv = qv  + (qv  == 0) * EPS
    pv = p   + (p   == 0) * EPS

    sv = np.where(sv == 0, EPS, sv)
    pv = np.where(pv == 0, EPS, pv)
    qv = np.where(qv == 0, EPS, qv)

    z = x - m
    term1 = qv * np.exp(qv * z / sv) + np.exp(z / (qv * sv)) / qv
    arg   = np.log((np.exp(qv * z / sv) + np.exp(z / (qv * sv))) / 2)
    term2 = np.exp(-np.abs(arg) ** pv / 2)

    denom = (2 ** (1 + 1.0 / pv)) * sv * _gamma_safe(1 + 1.0 / pv) * \
            (np.exp(qv * z / sv) + np.exp(z / (qv * sv)))

    y = np.abs(term1) * term2 / denom
    return y


def _gamma_safe(x):
    from scipy.special import gamma
    return gamma(x)
