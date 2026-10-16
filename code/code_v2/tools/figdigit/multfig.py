"""Digitise LMP2015 draft Fig. (Im of l=1 Compton multipoles): red solid (with gamma N Delta vertex
corrections), red dashed (without), black dotted (MAID). Axis calibration from the major ticks
(x: 25 MeV spacing, tick at 175 MeV = first interior tick; y: from labelled ticks)."""
import numpy as np
from pdfpaths import paths
FN = '<release>/LP_ms3/fig_pion_production_multipoles_l1_LS_v2.pdf'
# panels: (row, col) -> name, and y tick calibration (path y of two ticks and their values)
PANELS = {(0, 0): 'EE1+', (0, 1): 'EE1-', (1, 0): 'MM1+', (1, 1): 'MM1-', (2, 0): 'EM1+', (2, 1): 'ME1+'}

def load():
    P = paths(FN)
    ticks = []
    curves = []
    for p in P:
        if p['op'] != 'S':
            continue
        pts = np.array([pt for s in p['pts'] for pt in s])
        if p['scol'] == (0.0, 0.0, 0.0) and p['dash'] == '[ ] 0':
            for s in p['pts']:
                if len(s) == 2:
                    ticks.append((s[0][0], s[0][1], s[1][0], s[1][1]))
        elif len(pts) > 20 and p['scol'] in [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0)]:
            kind = {('[ ] 0', (1.0, 0.0, 0.0)): 'vc', ('[ 4 4 ] 0', (1.0, 0.0, 0.0)): 'novc',
                    ('[ 1 4 ] 0', (0.0, 0.0, 0.0)): 'maid'}.get((p['dash'], p['scol']))
            if kind:
                curves.append((kind, pts))
    return np.array(ticks), curves

def calib():
    T, curves = load()
    vert = T[np.abs(T[:, 0] - T[:, 2]) < 1e-6]
    hor = T[np.abs(T[:, 1] - T[:, 3]) < 1e-6]
    Lv = np.abs(vert[:, 1] - vert[:, 3]); Lh = np.abs(hor[:, 0] - hor[:, 2])
    xt = sorted(set(np.round(vert[(Lv > 7) & (Lv < 20), 0], 2)))
    out = {}
    # x: columns. interior ticks every 25 MeV; the first interior tick of each column = 175 MeV
    cols = [[x for x in xt if x < 580], [x for x in xt if x > 580]]
    xcal = []
    for c in cols:
        inner = c[1:-1]          # drop frame edges
        dx = np.mean(np.diff(inner))
        xcal.append((inner[0], dx))
    return xcal, curves, hor[(Lh > 7) & (Lh < 20)]

def ycal_from_labels():
    # label values (from pdftotext) mapped to path coordinates: y_path = y_label_center + 10.1
    # determined from the tick/label pairs ('1' <-> 439.1, '0' <-> 564.6 in panel MM1+)
    return None

if __name__ == '__main__':
    xcal, curves, ht = calib()
    print('x calibration (tick at 175 MeV, px per 25 MeV):', xcal)
    for kind, pts in curves:
        col = 0 if pts[:, 0].mean() < 580 else 1
        x0, dx = xcal[col]
        om = 175 + (pts[:, 0] - x0) / dx * 25
        i = np.argmin(pts[:, 1])   # y increases downward -> max value at min y
        print(kind, 'col', col, 'y-range %.1f..%.1f' % (pts[:, 1].min(), pts[:, 1].max()), ' omega range %.1f..%.1f' % (om.min(), om.max()),
              ' omega at extreme (min y): %.1f' % om[i])

# y calibration: (path y of two labelled major ticks, values) per panel, in units 1e-3/m_pi
YCAL = {'EE1+': ((87.8, 0.05), (284.5, 0.0)), 'EE1-': ((52.3, 0.5), (287.7, 0.0)),
        'MM1+': ((439.1, 1.0), (564.6, 0.0)), 'MM1-': ((439.1, 0.05), (564.6, 0.0)),
        'EM1+': ((640.4, 0.01), (851.4, 0.0)), 'ME1+': ((709.3, 0.05), (848.0, 0.0))}

def digitise():
    """returns {panel: {kind: (omega_cm[MeV], value[1e-3/mpi])}}"""
    T, curves = load()
    out = {}
    for kind, pts in curves:
        col = 0 if pts[:, 0].mean() < 580 else 1
        ym = pts[:, 1].mean()
        row = 0 if ym < 300 else (1 if ym < 590 else 2)
        name = PANELS[(row, col)]
        x0 = 187.4 if col == 0 else 659.2      # tick of omega = 150 MeV (label x + 57)
        dx = 63.27
        om = 150 + (pts[:, 0] - x0) / dx * 25
        (ya, va), (yb, vb) = YCAL[name]
        val = va + (pts[:, 1] - ya) * (vb - va) / (yb - ya)
        o = np.argsort(om)
        out.setdefault(name, {})[kind] = (om[o], val[o])
    return out
