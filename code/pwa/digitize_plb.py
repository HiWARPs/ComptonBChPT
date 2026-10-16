"""Recover the data points used in Krupina-Lensky-Pascalutsa PLB 782 (2018) 34 from the vector
graphics of its Fig. 5 (UnpolXS_new_v2.pdf, Mathematica export): marker centroids and error-bar
end points, mapped to (theta_lab, dsigma/dOmega_lab) with the panel tick marks.
Output: results/plb_fig5_points.csv  (bin, experiment(colour), theta_lab, value, err, x/y diagnostics)"""
import numpy as np
import pandas as pd
import pdfvec

PDF = '<release>/code/refs/KLP2018/src/UnpolXS_new_v2.pdf'
COLOURS = {(1.0, 0.0, 0.6): 'OlmosdeLeon', (0.0, 0.0, 1.0): 'Federspiel', (0.0, 0.0, 0.0): 'Hallin',
           (1.0, 0.0, 0.0): 'MacGibbon', (0.0, 1.0, 0.0): 'Baranov', (0.0, 0.502, 0.0): 'Goldansky',
           (1.0, 0.6, 0.0): 'Oxley', (0.8, 1.0, 0.0): 'Hyman', (0.4, 0.502, 0.0): 'Pugh',
           (0.502, 0.502, 1.0): 'Bernardini', (0.6, 0.0, 0.537): 'Zieger'}
# panel geometry (page points): x0 of left/right column, y0 of each row (value 0), scales
X0 = {0: 50.866, 1: 469.135}
XS = 69.0212 / 30.0            # pt per degree
Y0 = [1075.39, 871.32, 667.26, 463.19, 259.118, 55.049]
YS = 37.0495 / 5.0             # pt per nb/sr
PANELS = {(0, 0): 59, (1, 0): 69, (1, 1): 79, (2, 0): 89, (2, 1): 99, (3, 0): 109, (3, 1): 117,
          (4, 0): 127, (4, 1): 135, (5, 0): 143, (5, 1): 149}
W, H = 414.13, 200.07


def panel_of(x, y):
    col = 0 if x < 467 else 1
    for r, y0 in enumerate(Y0):
        if y0 - 15 <= y <= y0 + H + 15:
            return r, col
    return None


def main():
    st = pdfvec.streams(PDF)
    items, texts = pdfvec.interpret(st[0])
    markers, bars = [], []
    for it in items:
        if it['op'].startswith('f') and it['fill'] in COLOURS:
            for s in it['sub']:
                if 3 <= len(s) <= 6:
                    markers.append((COLOURS[it['fill']], s.mean(axis=0)[0], s.mean(axis=0)[1], len(s),
                                    np.ptp(s[:, 0]), np.ptp(s[:, 1])))
        if it['op'] == 'S' and it['stroke'] in COLOURS:
            for s in it['sub']:
                if len(s) == 2 and abs(s[0][0] - s[1][0]) < 1e-6 and abs(s[0][1] - s[1][1]) > 0:
                    bars.append((COLOURS[it['stroke']], s[0][0], min(s[0][1], s[1][1]), max(s[0][1], s[1][1])))
    M = pd.DataFrame(markers, columns=['exp', 'x', 'y', 'nv', 'wx', 'wy'])
    # remove duplicate markers (each marker drawn as fill + outline fill)
    M['xr'] = M.x.round(2); M['yr'] = M.y.round(2)
    M = M.drop_duplicates(['exp', 'xr', 'yr'])
    # drop black squares that are legend/other (none expected) and markers outside frames
    B = pd.DataFrame(bars, columns=['exp', 'x', 'ylo', 'yhi'])
    rows = []
    for _, m in M.iterrows():
        p = panel_of(m.x, m.y)
        if p is None or p not in PANELS:
            continue
        r, c = p
        th = (m.x - X0[c]) / XS
        val_marker = (m.y - Y0[r]) / YS
        # error bar: vertical segments of same colour at the marker x (bars may be clipped by the frame)
        bb = B[(B.exp == m.exp) & (abs(B.x - m.x) < 0.6) & (B.ylo <= m.y + 0.6) & (B.yhi >= m.y - 0.6)]
        if len(bb) == 0:
            # bars drawn as two halves (marker->top, bottom->marker)
            bb = B[(B.exp == m.exp) & (abs(B.x - m.x) < 0.6) & ((abs(B.ylo - m.y) < 0.6) | (abs(B.yhi - m.y) < 0.6))]
        # if several markers share this x (same exp), keep only segments touching this marker
        others = M[(M.exp == m.exp) & (abs(M.x - m.x) < 0.6) & (abs(M.y - m.y) > 0.6)]
        if len(others) and len(bb) > 1:
            touch = bb[(abs(bb.ylo - m.y) < 0.6) | (abs(bb.yhi - m.y) < 0.6)]
            if len(touch):
                bb = touch
        if len(bb):
            ylo, yhi = bb.ylo.min(), bb.yhi.max()
            clipped = (ylo < Y0[r] - 0.5) or (yhi > Y0[r] + H + 0.5)
            if len(bb) == 1 and not (abs(bb.ylo.iloc[0] - m.y) < 0.6 or abs(bb.yhi.iloc[0] - m.y) < 0.6):
                val = (0.5 * (ylo + yhi) - Y0[r]) / YS        # full bar through the marker
                err = 0.5 * (yhi - ylo) / YS
            else:
                # half-bars: error = larger half-length (symmetric errors)
                err = max(yhi - m.y, m.y - ylo) / YS
                val = (m.y - Y0[r]) / YS
            xbar = bb.x.mean()
        else:
            val, err, xbar, clipped = val_marker, np.nan, m.x, False
        rows.append(dict(bin=PANELS[p], exp=m.exp, theta_lab=round((xbar - X0[c]) / XS, 3), value=round(val, 4),
                         err=round(err, 4), value_marker=round(val_marker, 4), theta_marker=round(th, 3),
                         nvert=m.nv, nbar=len(bb), clipped=clipped))
    D = pd.DataFrame(rows).sort_values(['bin', 'exp', 'theta_lab'])
    return D


if __name__ == '__main__':
    D = main()
    D.to_csv('results/plb_fig5_points.csv', index=False)
    pd.set_option('display.width', 200)
    print(D.groupby('exp').size())
    print(len(D))
    print(D.to_string())


CURVES = {(0.0, 0.667, 0.667): 'Born', (0.0, 0.0, 1.0): 'fit1', (0.667, 0.0, 0.667): 'fit1p',
          (1.0, 0.0, 0.0): 'fit1pp', (0.0, 0.0, 0.0): 'BChPT'}


def curves():
    """theory curves of PLB Fig. 5: dict (bin, name) -> (theta_lab array, dsigma_lab array)."""
    st = pdfvec.streams(PDF)
    items, texts = pdfvec.interpret(st[0])
    out = {}
    for it in items:
        if it['op'] != 'S' or it['stroke'] not in CURVES:
            continue
        for s in it['sub']:
            if len(s) < 15:
                continue
            p = panel_of(s[:, 0].mean(), s[:, 1].mean())
            if p is None or p not in PANELS:
                continue
            r, c = p
            th = (s[:, 0] - X0[c]) / XS
            v = (s[:, 1] - Y0[r]) / YS
            key = (PANELS[p], CURVES[it['stroke']])
            out.setdefault(key, []).append((th, v))
    res = {}
    for k, lst in out.items():
        if len(lst) == 1:
            res[k] = lst[0]
        else:   # fit 1: central curve + two band edges -> central = median at each angle
            th0 = lst[0][0]
            vv = np.array([np.interp(th0, t, v) for t, v in lst])
            res[k] = (th0, np.median(vv, axis=0))
            res[(k[0], k[1] + '_lo')] = (th0, vv.min(axis=0))
            res[(k[0], k[1] + '_hi')] = (th0, vv.max(axis=0))
    return res
