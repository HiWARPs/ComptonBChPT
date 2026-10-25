"""bhigh.py -- two-point tensor coefficients of any rank (LoopTools/Denner basis, mu = 1, Delta_UV = 0) from the
Feynman-parameter form  B_{(00)^j 1^r}(p2, m0, m1) = int_0^1 dx (-x)^r D^j/(2^j j!) (Delta_UV - ln D + H_j),
D(x) = x^2 p2 - x (p2 + m0 - m1) + m0 - i eps.  Used for the rank-4 bubbles that LoopTools does not provide."""
import numpy as np
from math import factorial
_xg, _wg = np.polynomial.legendre.leggauss(48)

def _nodes(p2, m0, m1):
    """quadrature nodes on [0,1], split at the real zeros of D (log singularities), with x = a + (b-a) u^2-type
    clustering towards the split points"""
    a, b, c = p2, -(p2 + m0 - m1), m0
    cuts = [0.0, 1.0]
    if abs(a) > 1e-300:
        disc = b * b - 4 * a * c
        if disc > 0:
            for r in ((-b - np.sqrt(disc)) / (2 * a), (-b + np.sqrt(disc)) / (2 * a)):
                if 0 < r < 1: cuts.append(r)
    elif abs(b) > 0:
        r = -c / b
        if 0 < r < 1: cuts.append(r)
    cuts = sorted(cuts)
    X, W = [], []
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        mid = 0.5 * (lo + hi)
        for s, e in ((lo, mid), (mid, hi)):           # cluster towards both ends: x = s + (e-s) (1-cos(pi u))/2 style
            u = 0.5 * (_xg + 1); w = 0.5 * _wg
            x = s + (e - s) * u ** 3 if s == lo else e - (e - s) * (1 - u) ** 3
            jac = (e - s) * 3 * u ** 2 if s == lo else (e - s) * 3 * (1 - u) ** 2
            X.append(x); W.append(w * jac)
    return np.concatenate(X), np.concatenate(W)

def bcoef(j, r, p2, m0, m1):
    """(finite part, UV coefficient) of B_{(00)^j 1^r}"""
    x, w = _nodes(p2, m0, m1)
    D = x * x * p2 - x * (p2 + m0 - m1) + m0
    lnD = np.log(np.abs(D)) - 1j * np.pi * (D < 0)
    H = sum(1.0 / k for k in range(1, j + 1))
    pre = (-x) ** r * D ** j / (2 ** j * factorial(j))
    return np.sum(w * pre * (H - lnD)), np.sum(w * pre)

def names_rank4():
    return {'bb0000': (2, 0), 'bb0011': (1, 2), 'bb1111': (0, 4), 'bb00001': (2, 1), 'bb00111': (1, 3), 'bb11111': (0, 5)}
