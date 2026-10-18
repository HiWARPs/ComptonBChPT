#!/usr/bin/env python3
"""Build db_forward_obs_above.csv: the empirical forward Compton observables ABOVE the pion-production threshold, from the
Supplemental Material of Gryniuk, Hagelstein, Pascalutsa, PRD 94 (2016) 034043 [Gryniuk:2016gnm] (ten text files, format
"X  Y  dY" with X = nu in GeV; downloaded with the user's permission on 2026-10-17 into work/data/raw/GHP2016_supplement/;
all ten files provided by the user from the APS supplement pages of papers I and II; the paper-I files have no error column):
  ds_1.txt, ds_2.txt           unpolarised forward cross section dsigma/dOmega_lab(0) [nb/sr] for the two photoabsorption fits
  Sigma2z_ds_1.txt, _2.txt     polarised forward cross section -2 Re(f g*) [nb/sr] for the two fits (-> Sigma_2z(0) = ratio to ds_n)
  f_1.txt, f_2.txt, g.txt, cs_abs_1/2.txt, heldiffcs_abs.txt   amplitudes and photoabsorption cross sections (not needed here)
Like the below-threshold forward points (build_db_forward.py), the value is the MEAN of fits I and II with half their
difference added in quadrature to the mean quoted error; with only one fit available the file is marked 'fitI only'.
Energies kept: the panel energies of the Delta-region figures (work/delta_plots FIG2, and the Sigma_2z bins), linearly
interpolated on the 5-MeV grid of the supplement, plus a coarse 10-MeV grid 155-400 MeV for the notes.
PLOTTING ONLY (open stars at theta = 0 in Figs. deltacstheta and deltaS2); never used in a fit.

SINCE 2026-10-18 (default input): the 2026 update of the empirical forward amplitudes (public repository
github.com/HiWARPs/ForwardRCS, files copied to work/data/raw/ForwardRCS2026/).  When <raw_dir> contains obs.txt, the rows are
built from obs.txt (dsigma/dOmega_lab(0) [nb/sr], Sigma_2z(0), exp+fit errors on the 5-MeV grid), obs_selection.txt (their
data-selection systematic, from f) and the PWA-model systematic of g below 0.2 GeV (pol_g_pwa_model.txt, propagated to
dsigma(0) and Sigma_2z(0) with the tabulated f, g of unpol_f.txt / pol_g.txt; signed modes, so linear propagation); the three
are added in quadrature; values and errors interpolated linearly to the panel energies as before.  The GHP-based file
(<raw_dir> = work/data/raw/GHP2016_supplement) is kept as db_forward_obs_above_ghp2016.csv.

usage: python3 -I build_db_forward_above.py <raw_dir> <out_csv>
"""
import os, sys, csv, math
RAW, OUT = sys.argv[1], sys.argv[2]
PANELS = [160.0, 182.0, 213.1, 224.3, 234.7, 244.8, 254.1, 265.1, 275.6, 286.5, 298.2, 310.1, 323.8, 333.8, 350.5, 370.5, 389.0, 444.4, 454.9]
S2Z_BINS = [275.0, 295.0]                      # centres of the A2 Sigma_2z bins 265-285 and 285-305 MeV
COARSE = []                                    # (a 10-MeV grid was added until 2026-10-17; it produced a second star in the panels)


def read(name):
    """'nu  Y  [dY]' tables (the 2016 supplement has errors; the 2015 supplement of PRD 92, 074031 -- same ds_1/ds_2, f_1/f_2
    -- has none: then dY = None and the relative error of the other fit is used by combine())"""
    fn = os.path.join(RAW, name)
    if not os.path.exists(fn):
        return None
    d = {}
    for line in open(fn):
        p = line.split()
        if not p or p[0].startswith('#'):
            continue
        try:
            x, y = float(p[0]), float(p[1])
            e = float(p[2]) if len(p) > 2 else None
        except ValueError:
            continue
        d[round(1000 * x, 3)] = (y, e)         # key: nu in MeV
    return d


def interp(d, E):
    """linear interpolation of (value, error) of the 5-MeV grid to E (the panel energies of the figures); None outside"""
    ks = sorted(d)
    if E < ks[0] or E > ks[-1]:
        return None
    hi = min(k for k in ks if k >= E - 1e-9); lo = max(k for k in ks if k <= E + 1e-9)
    if hi == lo:
        return d[lo]
    t = (E - lo) / (hi - lo)
    e = None if (d[lo][1] is None or d[hi][1] is None) else d[lo][1] + t * (d[hi][1] - d[lo][1])
    return (d[lo][0] + t * (d[hi][0] - d[lo][0]), e)


def combine(d1, d2, E):
    """mean of fits I and II at the energy E (interpolated on the 5-MeV grid), error = mean quoted error (+) half-difference
    in quadrature; None if E is outside the grid"""
    v1 = interp(d1, E) if d1 else None
    v2 = interp(d2, E) if d2 else None
    if v1 and v2:
        (y1, e1), (y2, e2) = v1, v2
        note = 'mean of fits I and II; err = mean quoted (+) half-difference; interpolated on the 5-MeV grid'
        if e1 is None and e2 is None:
            raise ValueError('no errors for either fit at E = %g' % E)
        if e2 is None:
            e2 = e1 * abs(y2 / y1); note += '; fit-II error taken as the relative error of fit I'
        if e1 is None:
            e1 = e2 * abs(y1 / y2); note += '; fit-I error taken as the relative error of fit II'
        return E, 0.5 * (y1 + y2), math.sqrt((0.5 * (e1 + e2)) ** 2 + (0.5 * (y1 - y2)) ** 2), note
    v = v1 or v2
    if v is None or v[1] is None:
        return None
    return E, v[0], v[1], 'fit %s only; interpolated on the 5-MeV grid' % ('I' if v1 else 'II')


def build_2026():
    """rows from the ForwardRCS 2026 tables (see the module docstring)"""
    def rd(fn):
        return {round(1000 * float(l.split()[0]), 3): [float(x) for x in l.split()[1:]] for l in open(os.path.join(RAW, fn))
                if l.strip() and l.split()[0][0].isdigit()}
    o, osel, f, g, gp = rd('obs.txt'), rd('obs_selection.txt'), rd('unpol_f.txt'), rd('pol_g.txt'), rd('pol_g_pwa_model.txt')
    K = 1000 / 389.379                          # nb/sr per (mub GeV)^2
    cs, s2 = {}, {}
    for E in sorted(o):
        if E not in osel or E not in gp:
            continue
        c, ec, n, en, s, es = o[E]
        rf, imf, rg, img = f[E][0], f[E][2], g[E][0], g[E][2]
        dRg, dIg = gp[E][0], gp[E][1]                       # signed PWA-model mode of g
        dc = 2 * K * (rg * dRg + img * dIg)                  # d|g|^2 -> dsigma(0)
        dn = -2 * K * (rf * dRg + imf * dIg)                 # d(-2 Re f g*) [nb/sr]
        ds = dn / c - s * dc / c
        cs[E] = (c, math.sqrt(ec ** 2 + osel[E][0] ** 2 + dc ** 2))
        s2[E] = (s, math.sqrt(es ** 2 + osel[E][2] ** 2 + ds ** 2))
    src = 'ForwardRCS 2026 update (github.com/HiWARPs/ForwardRCS), obs.txt (+) obs_selection.txt (+) PWA model of g (pol_g_pwa_model.txt)'
    out = []
    for obs, tab, Es, nd in (('CS', cs, PANELS, 4), ('S2Z', s2, S2Z_BINS, 5)):
        for E in Es:
            v = interp(tab, E)
            if v is None:
                continue
            out.append(dict(obs=obs, E_lab_MeV=E, E_panel=E, theta_lab_deg=0.0, value=round(v[0], nd), stat_err=round(v[1], nd),
                            units='nb/sr' if obs == 'CS' else '1', frame='lab', source=src,
                            notes=('dsigma/dOmega_lab(0) = |f|^2 + |g|^2' if obs == 'CS' else 'Sigma_2z(0) = -2 Re(f g*)/(|f|^2+|g|^2), sigma_3/2 - sigma_1/2')
                            + ' above threshold; error = exp+fit (+) data selection (+) PWA model of g; interpolated linearly on the 5-MeV grid; PLOTTING ONLY'))
    return out


if os.path.exists(os.path.join(RAW, 'obs.txt')):
    rows = build_2026()
    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['obs', 'E_lab_MeV', 'E_panel', 'theta_lab_deg', 'value', 'stat_err', 'units', 'frame', 'source', 'notes'])
        w.writeheader(); w.writerows(rows)
    print('2026 tables; rows:', len(rows), ' CS:', sum(r['obs'] == 'CS' for r in rows), ' S2Z:', sum(r['obs'] == 'S2Z' for r in rows))
    sys.exit(0)

ds1, ds2 = read('ds_1.txt'), read('ds_2.txt')
s1, s2 = read('Sigma2z_ds_1.txt'), read('Sigma2z_ds_2.txt')
rows = []
src = 'Gryniuk:2016gnm Supplemental Material (PRD 94, 034043), files ds_1/2.txt, Sigma2z_ds_1/2.txt'
if ds1 or ds2:
    for E in sorted(set(PANELS + COARSE)):
        c = combine(ds1, ds2, E)
        if c is None:
            continue
        k, y, e, note = c
        rows.append(dict(obs='CS', E_lab_MeV=k, E_panel=E, theta_lab_deg=0.0, value=round(y, 4), stat_err=round(e, 4), units='nb/sr', frame='lab',
                         source=src, notes='dsigma/dOmega_lab(0) = |f|^2 + |g|^2 above threshold; ' + note + '; PLOTTING ONLY'))
if s1 or s2:
    # the 2016 files Sigma2z_ds_1/2.txt tabulate the polarised forward cross section -2 Re(f g*) [nb/sr] (header checked);
    # Sigma_2z(0) = (-2 Re f g*)/(|f|^2 + |g|^2) is formed PER FIT and then the two fits are combined as for the cross section.
    def ratio(sd, dd):
        if not (sd and dd):
            return None
        out = {}
        for k, (y, e) in sd.items():
            if k in dd and dd[k][0] > 0:
                yd, ed = dd[k]
                r = y / yd
                out[k] = (r, abs(r) * math.sqrt((e / y) ** 2 + (ed / yd) ** 2) if y != 0 else e / yd)
        return out
    r1, r2 = ratio(s1, ds1), ratio(s2, ds2)
    for E in sorted(set(S2Z_BINS + COARSE)):
        c = combine(r1, r2, E)
        if c is None:
            continue
        k, y, e, note = c
        rows.append(dict(obs='S2Z', E_lab_MeV=k, E_panel=E, theta_lab_deg=0.0, value=round(y, 5), stat_err=round(e, 5), units='1', frame='lab',
                         source=src, notes='Sigma_2z(0) = -2 Re(f g*)/(|f|^2+|g|^2) from Sigma2z_ds_n.txt / ds_n.txt per fit, sigma_3/2 - sigma_1/2 convention as below threshold; ' + note + '; PLOTTING ONLY'))
with open(OUT, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['obs', 'E_lab_MeV', 'E_panel', 'theta_lab_deg', 'value', 'stat_err', 'units', 'frame', 'source', 'notes'])
    w.writeheader(); w.writerows(rows)
print('files found:', [n for n, d in (('ds_1', ds1), ('ds_2', ds2), ('Sigma2z_ds_1', s1), ('Sigma2z_ds_2', s2)) if d])
print('rows:', len(rows), ' CS:', sum(r['obs'] == 'CS' for r in rows), ' S2Z:', sum(r['obs'] == 'S2Z' for r in rows))
