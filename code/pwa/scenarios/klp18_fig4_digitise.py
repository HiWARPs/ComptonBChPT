"""Recover the published dynamical polarisabilities of Krupina, Lensky, Pascalutsa, PLB 782 (2018) 34
(KLP18), Fig. 4 (arXiv source file polarizabilities_NEW_v2.pdf, Mathematica vector export), for
Fit 1 (blue) and Fit 1'' (red): every point is drawn as a vertical error bar with two caps; the value is
the mid-point of the bar and the error its half-length (symmetric errors).  Axes are calibrated with
the major tick marks of each panel (exact vector coordinates, no pixel digitisation).
Output: results/klp18_fig4_points.csv  (fit, key, E_lab [MeV], value, err)."""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PWA = os.path.dirname(HERE)
sys.path.insert(0, PWA)
import pdfvec  # noqa: E402

PDF = '<release>/code/refs/KLP2018/src/polarizabilities_NEW_v2.pdf'
OUT = os.path.join(PWA, 'results', 'klp18_fig4_points.csv')
COL = {(0.0, 0.0, 1.0): "Fit 1", (1.0, 0.0, 0.0): "Fit 1''"}
# panel frames (page points): column x0 (E = 0), row y-range; y calibration from two major ticks
XE0 = {0: 75.37, 1: 547.30}
XSC = (458.68 - 75.37) / 150.0                       # pt per MeV (ticks at 0, 50, 100, 150 MeV)
ROWS = {0: (546.19, 785.38), 1: (302.22, 541.41), 2: (58.25, 297.44)}
# (row, col) -> (key, (y_tick_a, val_a), (y_tick_b, val_b))
PANELS = {(0, 0): ('alphaE1', (546.19, 0.0), (738.31, 20.0)),
          (0, 1): ('betaM1', (546.19, 0.0), (730.18, 10.0)),
          (1, 0): ('gE1E1', (350.06, -15.0), (529.45, 0.0)),
          (1, 1): ('gM1M1', (302.22, 0.0), (501.88, 5.0)),
          (2, 0): ('gE1M2', (112.15, 0.0), (280.59, 1.0)),
          (2, 1): ('gM1E2', (84.89, 0.0), (244.70, 3.0))}


def digitise():
    items, _ = pdfvec.interpret(pdfvec.streams(PDF)[0])
    bars = set()
    for it in items:
        if it['op'] != 'S' or it['stroke'] not in COL:
            continue
        for s in it['sub']:
            if len(s) == 2 and abs(s[0][0] - s[1][0]) < 1e-6 and abs(s[0][1] - s[1][1]) > 1e-6:
                bars.add((COL[it['stroke']], round(s[0][0], 3), round(min(s[0][1], s[1][1]), 3),
                          round(max(s[0][1], s[1][1]), 3)))
    rows = []
    for fit, x, ylo, yhi in sorted(bars):
        col = 0 if x < 527 else 1
        row = [r for r, (a, b) in ROWS.items() if a - 1 <= 0.5 * (ylo + yhi) <= b + 1]
        if len(row) != 1:
            continue
        key, (ya, va), (yb, vb) = PANELS[(row[0], col)]
        sc = (vb - va) / (yb - ya)
        E = (x - XE0[col]) / XSC
        rows.append(dict(fit=fit, key=key, E_lab=round(E, 2), value=va + (0.5 * (ylo + yhi) - ya) * sc,
                         err=0.5 * (yhi - ylo) * sc, clipped=bool(ylo <= ROWS[row[0]][0] + 0.01 or yhi >= ROWS[row[0]][1] - 0.01)))
    d = pd.DataFrame(rows).sort_values(['fit', 'key', 'E_lab']).reset_index(drop=True)
    d['E_lab'] = d.E_lab.round(0)
    return d


if __name__ == '__main__':
    d = digitise()
    d.to_csv(OUT, index=False)
    print(d.groupby(['fit', 'key']).size())
    print(d.pivot_table(index='E_lab', columns=['fit', 'key'], values='value').round(2).to_string())
