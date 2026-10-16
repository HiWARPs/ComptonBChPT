#!/usr/bin/env python3
"""
Unpolarised dsigma/dOmega_lab in the Delta region vs LMP2015 Fig. 14 (published EPJC PDF, p. 16):
the red solid curves (full calculation incl. the gamma N Delta vertex corrections) are digitised
from a 600-dpi rendering (pdftoppm; axis calibration from the major tick marks: x 303.5 px/100 MeV,
y 178.5 px/100 nb/sr (upper row) and 198.5 px/50 nb/sr (lower row)) and compared with v2 with and
without the vertex corrections at E = 150..400 MeV.  Digitisation uncertainty: ~0.5 nb/sr in y,
~0.5 MeV in x (matters on the steep flanks of the resonance).
usage: validate_fig14.py [path/to/s10052-015-3791-0.pdf]
"""
import os, sys, subprocess, tempfile
import numpy as np
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
PDF = sys.argv[1] if len(sys.argv) > 1 else os.path.abspath(os.path.join(HERE, '..', '..', '..', 's10052-015-3791-0.pdf'))
DRV = os.path.join(HERE, '..', 'driver', 'run_points.sh')


def digitise():
    d = tempfile.mkdtemp(prefix='f14_')
    subprocess.run(['pdftoppm', '-f', '16', '-l', '16', '-r', '600', '-png', PDF, os.path.join(d, 'p')], check=True)
    fn = [f for f in os.listdir(d) if f.endswith('.png')][0]
    im = np.array(Image.open(os.path.join(d, fn)).convert('RGB')).astype(int)
    y0, x0 = (540 + 40) * 6, 100 * 6
    sub = im[y0:y0 + 2200, x0:x0 + 4000]
    red = (sub[:, :, 0] > 200) & (sub[:, :, 1] < 70) & (sub[:, :, 2] < 70)
    cols = [(190.5, 1420), (1433.5, 2663), (2676.5, 3906)]
    rows = [(51, 922, 905.5, 178.5 / 100), (1012, 1883, 1883.5, 198.5 / 50)]
    panels = {45: (0, 0), 60: (0, 1), 90: (0, 2), 120: (1, 0), 135: (1, 1), 180: (1, 2)}
    out = {}
    for th, (r, c) in panels.items():
        cx0, cx1 = cols[c]; ry0, ry1, yzero, ppu = rows[r]
        E = []; V = []
        for xp in range(int(cx0) + 3, int(cx1) - 3):
            rr = np.where(red[ry0 + 3:ry1 - 3, xp])[0]
            if len(rr) == 0:
                continue
            g = max(np.split(rr, np.where(np.diff(rr) > 2)[0] + 1), key=len)
            if len(g) < 3:
                continue
            E.append((xp - cx0) / 303.5 * 100); V.append((yzero - (g.mean() + ry0 + 3)) / ppu)
        out[th] = (np.array(E), np.array(V))
    return out


def main():
    D = digitise()
    Es = np.arange(150, 401, 10)
    d = tempfile.mkdtemp(prefix='f14r_')
    p = os.path.join(d, 'p.txt')
    with open(p, 'w') as f:
        for th in D:
            for E in Es:
                f.write('cs_%d_%d %d %d CS\n' % (th, E, E, th))
    res = {}
    for tag, key in [('noVC', 'dvertex=f'), ('VC', 'dvertex=t')]:
        o = os.path.join(d, tag + '.out')
        subprocess.run([DRV, p, o, key], check=True, capture_output=True)
        res[tag] = {l.split()[0]: float(l.split()[5]) for l in open(o) if not l.startswith('#')}
    print('dsigma/dOmega_lab [nb/sr]: LMP2015 Fig. 14 red curve (digitised) vs v2 with / without vertex corrections')
    allr = {'VC': [], 'noVC': []}
    for th in D:
        E, V = D[th]
        line = 'theta_lab=%3d:' % th
        for e in [200, 250, 300, 330, 350, 400]:
            m = np.abs(E - e) < 1.0
            line += '  %d: %6.1f|%6.1f|%6.1f' % (e, np.median(V[m]) if m.any() else np.nan,
                                                 np.interp(e, Es, [res['VC']['cs_%d_%d' % (th, x)] for x in Es]),
                                                 np.interp(e, Es, [res['noVC']['cs_%d_%d' % (th, x)] for x in Es]))
        print(line)
        for tag in ('VC', 'noVC'):
            for e in Es:
                m = np.abs(E - e) < 1.0
                if m.any():
                    fig = np.median(V[m]); v = res[tag]['cs_%d_%d' % (th, e)]
                    allr[tag].append(abs(v / fig - 1))
    for tag in ('VC', 'noVC'):
        a = np.array(allr[tag])
        print('%-5s: |v2/LMP - 1| over %d points (150-400 MeV, 6 angles): median %.4f, 90%% quantile %.4f, max %.4f' %
              (tag, len(a), np.median(a), np.quantile(a, 0.9), a.max()))


if __name__ == '__main__':
    main()
