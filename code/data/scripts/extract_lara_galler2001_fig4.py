#!/usr/bin/env python3
"""Recover the LARA (Bonn) proton Compton cross sections of Galler et al., PLB 503 (2001) 245 [Galler:2001ht]
(same data in Wolf et al., EPJA 12 (2001) 231 [Wolf:2001ha]) from the vector graphics of Fig. 4 of the arXiv source
(nucl-ex/0102003, figures/fig4.eps, written by xmgr): three stacked panels, dsigma/dOmega_cm at theta_cm = 75, 90,
125 deg (angular bins of 15 deg) versus the photon lab energy, 265-780 MeV.  The data were never tabulated.

Method: the PostScript coordinates of the filled-circle symbols ('r x y fa', the LARA points) and of the vertical
error-bar segments ('x y1 m x y2 l') are read from the EPS; the axes are calibrated with the tick labels (x labels
centred 'CS', y labels right-justified 'RJ').  Calibration residuals <= 0.13 MeV and <= 0.15 nb/sr; the recovered
energies agree with the tagger energies quoted in Wolf et al. (264.9, 275.6, ..., 375.0 MeV) to <= 0.2 MeV.
The error bars are statistics (+) point-to-point systematics, as plotted (Galler et al., text to Figs. 4-6); the
common scale error of +-3 % (tagging 2 %, target 2 %) goes to norm_err_frac.

Output: CSV with the columns of work/data/db_old.csv -- PLOTTING ONLY (work/delta_plots, dsigma/dOmega vs E), never
used in fits or chi^2 tables -- and a check figure (<out>.png) to compare by eye with the rendered original.

Run (from anywhere; the EPS is an untrusted download, so it is only read):
  python3 -I extract_lara_galler2001_fig4.py <path/to/Galler2001_arXiv/src/figures/fig4.eps> <out_csv>
"""
import re, sys, csv, math
import numpy as np

MP = 938.272
PANEL_THETA = {2: 75.0, 1: 90.0, 0: 125.0}       # panel index counted from the bottom of the figure


def th_cm2lab(E, thc):
    b = E / (E + MP)
    c = math.cos(math.radians(thc))
    return math.degrees(math.acos(max(-1.0, min(1.0, (c + b) / (1 + b * c)))))


fn, out_csv = sys.argv[1], sys.argv[2]
lines = open(fn, encoding='latin-1').read().splitlines()
labels, syms, segs = [], [], []
last_translate, prev = None, None
for i, ln in enumerate(lines):
    m = re.match(r'^(-?\d+) (-?\d+) translate$', ln)
    if m:
        last_translate = (int(m.group(1)), int(m.group(2)))
    m = re.match(r'^\((.*)\) (CS|RJ)$', ln)
    if m and last_translate:
        labels.append((last_translate[0], last_translate[1], m.group(1), m.group(2)))
    m = re.match(r'^(\d+) (-?\d+) (-?\d+) (fa|oa|osq|ot1|ot2|ot3|ot4|odi)$', ln)
    if m:
        syms.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4), i))
    m = re.match(r'^(-?\d+) (-?\d+) l$', ln)
    if m and prev:
        pm = re.match(r'^(-?\d+) (-?\d+) m$', prev)
        if pm:
            segs.append((int(pm.group(1)), int(pm.group(2)), int(m.group(1)), int(m.group(2)), i))
    prev = ln

# x calibration: the centred numeric labels of the (shared) bottom axis
xl = [(x, float(t)) for x, y, t, j in labels if j == 'CS' and re.match(r'^\d+$', t)]
xs = np.array([a for a, b in xl]); xv = np.array([b for a, b in xl])
px = np.polyfit(xs, xv, 1)
print('x labels', [int(v) for v in xv], '; E = %.5f*x + %.3f ; resid' % tuple(px), np.round(np.polyval(px, xs) - xv, 3))
# y calibration per panel: right-justified numeric labels, panels separated where the label value drops
yl = sorted([(y, float(t), x) for x, y, t, j in labels if j == 'RJ' and re.match(r'^\d+$', t)], key=lambda q: q[0])
panels, cur = [], [yl[0]]
for q in yl[1:]:
    if q[1] <= cur[-1][1]:
        panels.append(cur); cur = [q]
    else:
        cur.append(q)
panels.append(cur)
cal = []
for k, p in enumerate(panels):
    ys = np.array([q[0] for q in p]); vs = np.array([q[1] for q in p])
    py = np.polyfit(ys, vs, 1)
    print('panel %d (bottom=0): y labels %s ; v = %.5f*y + %.3f ; resid %s' % (k, [int(q[1]) for q in p], py[0], py[1], np.round(np.polyval(py, ys) - vs, 3)))
    cal.append((py, ys.min(), ys.max()))
assert len(cal) == len(PANEL_THETA), 'expected %d panels' % len(PANEL_THETA)
bounds = []
for k in range(len(cal)):
    lo = -1e9 if k == 0 else 0.5 * (cal[k - 1][2] + cal[k][1])
    hi = 1e9 if k == len(cal) - 1 else 0.5 * (cal[k][2] + cal[k + 1][1])
    bounds.append((lo, hi))
vsegs = [s for s in segs if s[0] == s[2] and s[1] != s[3]]
out = []
for r, x, y, kind, li in syms:
    if kind != 'fa':                      # filled circles = LARA; the open symbols are the earlier experiments
        continue
    k = [j for j, (lo, hi) in enumerate(bounds) if lo <= y < hi][0]
    py = cal[k][0]
    E = float(np.polyval(px, x)); v = float(np.polyval(py, y))
    cand = [s for s in vsegs if abs(s[0] - x) <= 1 and min(s[1], s[3]) - 2 <= y <= max(s[1], s[3]) + 2]
    note = ''
    if cand:
        cand.sort(key=lambda s: abs(0.5 * (s[1] + s[3]) - y))     # the bar centred on the symbol
        s = cand[0]
        mid = 0.5 * (s[1] + s[3]); err = 0.5 * abs(s[3] - s[1]) * py[0]
        if abs(mid - y) > 2:
            note = 'error-bar midpoint off by %d units' % int(mid - y)
    else:
        err = float('nan'); note = 'no error bar found'
    th = PANEL_THETA[k]
    out.append(dict(dataset='Galler2001_LARA', texkey='Galler:2001ht;Wolf:2001ha', obs='CS', E_lab_MeV=round(E, 1), dE_MeV='',
                    theta_lab_deg=round(th_cm2lab(E, th), 3), theta_cm_deg=th, frame_of_cs='cm', value=round(v, 1),
                    stat_err=round(err, 1), sys_err_pt='', norm_err_frac=0.03, units='nb/sr',
                    source='Galler et al. PLB 503,245 (2001) Fig. 4 (arXiv nucl-ex/0102003 figures/fig4.eps, xmgr vector graphics)',
                    notes='DIGITIZED from the EPS vector graphics (extract_lara_galler2001_fig4.py; eps line %d); '
                          'orig_angle=cm, 15-deg c.m. angular bin; error = stat (+) point-to-point syst as plotted; '
                          'scale 3%% (tagging 2%%, target 2%%); PLOTTING ONLY, not used in fits%s' % (li + 1, ('; ' + note) if note else '')))
out.sort(key=lambda o: (-o['theta_cm_deg'], o['E_lab_MeV']))
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
print('LARA points:', len(out), ' per angle:', {t: sum(1 for o in out if o['theta_cm_deg'] == t) for t in PANEL_THETA.values()})
print('notes:', sum(1 for o in out if 'off by' in o['notes'] or 'no error' in o['notes']))

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, axs = plt.subplots(3, 1, figsize=(6.5, 8), sharex=True)
for ax, th in zip(axs, (75.0, 90.0, 125.0)):
    o = [q for q in out if q['theta_cm_deg'] == th]
    ax.errorbar([q['E_lab_MeV'] for q in o], [q['value'] for q in o], yerr=[q['stat_err'] for q in o], fmt='o', ms=4, color='k', capsize=0)
    ax.set_xlim(175, 870); ax.set_ylim(0, 300); ax.grid(True, lw=0.3)
    ax.text(0.5, 0.85, r'$\theta_{cm}=%g^\circ$ (%d points)' % (th, len(o)), transform=ax.transAxes, ha='center')
axs[-1].set_xlabel(r'$E_\gamma$ [MeV]')
fig.savefig(out_csv.replace('.csv', '_check.png'), dpi=110)
print('check figure:', out_csv.replace('.csv', '_check.png'))
