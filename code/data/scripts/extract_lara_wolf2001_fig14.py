#!/usr/bin/env python3
"""Recover the LARA (Bonn) angular distributions of Wolf et al., EPJA 12 (2001) 231 [Wolf:2001ha] from the vector
graphics of Fig. 14 of the arXiv source (nucl-ex/0109013, figures/fig14.eps, written by HIGZ): eight panels of
dsigma/dOmega_cm versus theta_cm at the tagger energies 385.8 ... 454.9 MeV, "one-dimensional" (missing-time)
analysis, about 23 angular points per panel between 60 and 160 deg.  The data were never tabulated.

Method: the PostScript tokens are parsed in order.  Clip boxes ('w h x y C') delimit the panels; the filled circles
('x y m20') are the data points; the error bars are drawn as vertical segments with horizontal caps of 22-23 units
('x-11 y m 23 X s'), whose outermost y positions give the bar ends; the numeric tick labels (x labels centred,
'xs 2 div neg 0 t'; y labels right-justified, 'xs neg 0 t') are matched to the nearest major tick marks and fitted
linearly, which calibrates the axes of every panel (all panels have the same frame and the same ranges).
The energy of a panel is the '(xxx.x MeV)' label placed inside its clip box.
The error bars are the individual (random) errors as plotted: counting statistics (+) point-to-point systematics
(detection efficiency, geometry, background subtraction); the common scale error of 3 % (tagging 2 %, target 2 %)
goes to norm_err_frac, as for the Fig. 4 data of Galler et al. (db_lara_digitized.csv).

Output: CSV with the columns of work/data/db_old.csv -- PLOTTING ONLY (work/delta_plots, dsigma/dOmega vs theta),
never used in fits or chi^2 tables -- and a check figure (<out>.png) to compare by eye with the rendered original.

Run (from anywhere; the EPS is an untrusted download, so it is only read):
  python3 -I extract_lara_wolf2001_fig14.py <path/to/Wolf2001_arXiv/src/figures/fig14.eps> <out_csv>
"""
import re, sys, csv, math
import numpy as np

MP = 938.272


def th_cm2lab(E, thc):
    b = E / (E + MP)
    c = math.cos(math.radians(thc))
    return math.degrees(math.acos(max(-1.0, min(1.0, (c + b) / (1 + b * c)))))


fn, out_csv = sys.argv[1], sys.argv[2]
txt = open(fn, encoding='latin-1').read()
txt = txt[txt.index('%%EndProlog'):]
# tokens: PostScript strings '(...)' kept as one token
tok = re.findall(r'\([^()]*\)|[^\s]+', txt)
num = lambda t: re.fullmatch(r'-?\d+(\.\d+)?', t) is not None

boxes, markers, hsegs, vsegs, labels = [], [], [], [], []
box = None
gs_pos, align = None, 'L'
i = 0
while i < len(tok):
    t = tok[i]
    if t == 'C' and i >= 4 and all(num(x) for x in tok[i - 4:i]):
        w, h, x, y = (float(v) for v in tok[i - 4:i])
        box = (x, y, w, h); boxes.append(box)
    elif t == 'm20' and i >= 2 and num(tok[i - 1]) and num(tok[i - 2]):
        markers.append((float(tok[i - 2]), float(tok[i - 1]), box))
    elif t == 's' and i >= 4 and tok[i - 1] in ('X', 'Y') and tok[i - 3] == 'm' and num(tok[i - 2]) \
            and num(tok[i - 4]) and num(tok[i - 5]):
        x, y, d = float(tok[i - 5]), float(tok[i - 4]), float(tok[i - 2])
        (hsegs if tok[i - 1] == 'X' else vsegs).append((x, y, d))
    elif t == 'gsave' and i + 2 < len(tok) and num(tok[i + 1]) and num(tok[i + 2]):
        gs_pos = (float(tok[i + 1]), float(tok[i + 2])); align = 'L'
    elif t == 'xs' and tok[i + 1:i + 5] == ['2', 'div', 'neg', '0']:
        align = 'C'
    elif t == 'xs' and tok[i + 1:i + 3] == ['neg', '0']:
        align = 'R'
    elif t.startswith('(') and i + 1 < len(tok) and tok[i + 1] == 'show' and gs_pos is not None:
        labels.append((t[1:-1], gs_pos, align))
    i += 1

boxes = sorted(set(boxes))
# ---- calibration from the tick labels (nearest major tick mark at the frame edge), panel by panel
xlab = [(float(s), p) for s, p, a in labels if a == 'C' and num(s)]
ylab = [(float(s), p) for s, p, a in labels if a == 'R' and num(s)]
xt = [(x, y) for x, y, d in vsegs if abs(d) in (26, 27)]              # vertical tick marks (x, y of the frame edge)
yt = [(x, y) for x, y, d in hsegs if abs(d) in (26, 27)]              # horizontal tick marks
inbox = lambda b, X, Y: b[2] < 1000 and b[0] - 2 <= X <= b[0] + b[2] + 2 and b[1] - 2 <= Y <= b[1] + b[3] + 2
per = {}
for v, (X, Y) in xlab:
    k = min(xt, key=lambda u: abs(u[0] - X) + 0.01 * abs(u[1] - Y)); assert abs(k[0] - X) < 15, (v, X, k)
    b = [b for b in boxes if inbox(b, *k)][0]
    per.setdefault(('x', b), []).append((k[0], v))
for v, (X, Y) in ylab:
    k = min(yt, key=lambda u: abs(u[1] - (Y + 16)) + 0.01 * abs(u[0] - X)); assert abs(k[1] - Y - 16) < 15, (v, Y, k)
    b = [b for b in boxes if inbox(b, *k)][0]
    per.setdefault(('y', b), []).append((k[1], v))
rngs = {'x': [], 'y': []}
for (ax_name, b), pts in per.items():
    p, v = zip(*pts)
    sl, ic = np.polyfit(p, v, 1)
    res = np.max(np.abs(np.polyval([sl, ic], p) - v))
    lo, hi = (b[0], b[0] + b[2]) if ax_name == 'x' else (b[1], b[1] + b[3])
    rngs[ax_name].append(np.polyval([sl, ic], [lo, hi]))
    print('%s calibration of panel at (%g, %g): %d labels, max residual %.3f, frame range %.3f ... %.3f' %
          (ax_name, b[0], b[1], len(p), res, *rngs[ax_name][-1]))
for k in rngs:
    r = np.array(rngs[k]); assert np.ptp(r, axis=0).max() < (0.5 if k == 'x' else 1.0), r
X0, X1 = np.mean(rngs['x'], axis=0); Y0, Y1 = np.mean(rngs['y'], axis=0)
print('frame ranges: theta %.2f ... %.2f deg, dsigma %.2f ... %.2f nb/sr' % (X0, X1, Y0, Y1))
assert abs(X1 - X0 - 150) < 1.0 and abs(Y1 - Y0 - 350) < 1.5, (X0, X1, Y0, Y1)


def cal(b, x, y):
    return X0 + (x - b[0]) / b[2] * (X1 - X0), Y0 + (y - b[1]) / b[3] * (Y1 - Y0)


# ---- panel energies
E_of_box = {}
for s, (X, Y), a in labels:
    m = re.fullmatch(r'(\d+\.\d) MeV', s)
    if m:
        bb = [b for b in boxes if inbox(b, X, Y)]
        assert len(bb) == 1, (s, X, Y, bb)
        E_of_box[bb[0]] = float(m.group(1))
print('panels:', sorted(E_of_box.values()))

# ---- points with error bars (caps: horizontal segments of 22-23 units centred on the marker x)
caps = [(x + d / 2, y) for x, y, d in hsegs if abs(d) in (22, 23)]
rows = []
for x, y, b in markers:
    assert b in E_of_box, (x, y, b)
    cc = [cy for cx, cy in caps if abs(cx - x) <= 1.5 and abs(cy - y) <= 200]
    lo = max([c for c in cc if c < y], default=None); hi = min([c for c in cc if c > y], default=None)
    th, v = cal(b, x, y)
    if lo is None or hi is None:
        # no caps (bar hidden by the frame?): the ends of the vertical segments below and above the marker
        vv = [(y1, y1 + d) for x1, y1, d in vsegs if abs(x1 - x) <= 1.5 and abs(y1 - y) <= 200]
        lo = max([min(a) for a in vv if min(a) < y], default=None); hi = min([max(a) for a in vv if max(a) > y], default=None)
        print('  no caps for point at E=%.1f theta=%.1f: bar from the vertical segments %s' % (E_of_box[b], th, (lo, hi)))
    if lo is None or hi is None:
        err = float('nan')
    else:
        err = 0.5 * (hi - lo) / b[3] * (Y1 - Y0)
        if abs((hi - y) - (y - lo)) > 2:
            print('  asymmetric bar at E=%.1f theta=%.1f: -%d +%d units' % (E_of_box[b], th, y - lo, hi - y))
    rows.append((E_of_box[b], th, v, err))
rows.sort()
print('points: %d (%s per panel)' % (len(rows), {E: sum(1 for r in rows if r[0] == E) for E in sorted(E_of_box.values())}))

src = 'Wolf et al. EPJA 12,231 (2001) Fig. 14 (arXiv nucl-ex/0109013 figures/fig14.eps, HIGZ vector graphics)'
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['dataset', 'texkey', 'obs', 'E_lab_MeV', 'dE_MeV', 'theta_lab_deg', 'theta_cm_deg',
                                      'frame_of_cs', 'value', 'stat_err', 'sys_err_pt', 'norm_err_frac', 'units', 'source',
                                      'notes'])
    w.writeheader()
    for E, th, v, e in rows:
        w.writerow(dict(dataset='Wolf2001_LARA_angdist', texkey='Wolf:2001ha;Galler:2001ht', obs='CS', E_lab_MeV=E, dE_MeV='',
                        theta_lab_deg=round(th_cm2lab(E, th), 3), theta_cm_deg=round(th, 2), frame_of_cs='cm',
                        value=round(v, 1), stat_err=round(e, 1), sys_err_pt='', norm_err_frac=0.03, units='nb/sr',
                        source=src,
                        notes='DIGITIZED from the EPS vector graphics (extract_lara_wolf2001_fig14.py); one-dimensional '
                              '(missing-time) analysis; orig_angle=cm; error = individual (stat (+) point-to-point syst) as '
                              'plotted; scale 3% (tagging 2%, target 2%); PLOTTING ONLY, not used in fits'))

# ---- check figure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
Es = sorted(E_of_box.values())
fig, axs = plt.subplots(3, 3, figsize=(9, 9), sharex=True, sharey=True)
for ax, E in zip(axs.ravel(), Es):
    r = np.array([(th, v, e) for EE, th, v, e in rows if EE == E])
    ax.errorbar(r[:, 0], r[:, 1], r[:, 2], fmt='o', ms=3, color='k')
    ax.set_title('%.1f MeV' % E); ax.set_xlim(30, 180); ax.set_ylim(0, 350)
    ax.set_xticks([60, 90, 120, 150]); ax.grid(ls=':')
axs.ravel()[-1].axis('off')
fig.savefig(out_csv.replace('.csv', '_check.png'), dpi=90)
print('written', out_csv)
