"""Forward (theta=0) 'data points' Re f, Re g and the derived observables below the pion threshold.

usage: python3 -I build_db_forward.py [out_dir] [source]      source = '2026' (default) | 'ghp2016'

source '2026' (default since 2026-10-18): the 2026 update of the empirical forward amplitudes, public repository
  github.com/HiWARPs/ForwardRCS (O. Gryniuk, V. Pascalutsa), file forward_points.json copied to
  work/data/raw/ForwardRCS2026/: Re f, Re g [mub GeV] (unrounded) at the 12 energies, exp+fit covariances cov_f, cov_g
  (f and g independent) and the separate systematics as signed modes sel_f (data selection of f) and pwa_g (PWA
  model of g below 0.2 GeV).  Error model per energy i:
     total_i   = sqrt(cov_ii + sys_i^2)                              (sys = sel_f resp. pwa_g)
     coherent_i = sqrt(v_i^2 + sys_i^2),  v = leading eigenmode of cov (sqrt(lambda_0) x eigenvector)
     err_unc_i = sqrt(total_i^2 - coherent_i^2) = sqrt(cov_ii - v_i^2)  (uncorrelated remainder)
  -> db_forward_amp.csv, db_forward_obs.csv (columns as before + err_unc; stat_err = total).  The coherent part is
  carried by ONE nuisance per amplitude in the fits (work/fits/lib/forward.py, work/pwa/sumrules.py); this is exact
  when the leading mode and the systematic have the same energy shape (checked below: 99.9 % / 97.5 % of the exp+fit
  variance in the leading mode of f / g; |sel_f|/mode constant to 3 %, |pwa_g|/mode rising from 1.03 to 1.40, so the
  one-nuisance model is exact for f and approximate for g).
source 'ghp2016': Gryniuk, Hagelstein, Pascalutsa, PRD 92 (2015) 074031 [Gryniuk:2015eza] (Re f) and PRD 94 (2016)
  034043 [Gryniuk:2016gnm] (Re g); values transcribed from Table 'tab:values_forward_amps' of arXiv:1604.00789 (units
  mub*GeV); Re f = average of fit I and fit II, with half their difference added in quadrature to the quoted error
  (the 2-decimal rounding is the uncorrelated part, 0.005/sqrt(3), in forward.py).  -> db_forward_amp_ghp2016.csv,
  db_forward_obs_ghp2016.csv (identical to the db_forward_{amp,obs}.csv used until 2026-10-18).
Below pion threshold (nu < 150.7 MeV) the amplitudes are real, so
  dsigma/dOmega_lab(0) = |f|^2 + |g|^2 ,  Sigma_2z(0) = -2 Re(f g*)/(|f|^2+|g|^2)
[Eqs. (dsi), (Si2z) of arXiv:1604.00789].
"""
import csv, math, sys, os, json
MUB_GEV2 = 389.379  # 1 GeV^-2 = 389.379 mub ; |f|^2 [mub] = f[mub GeV]^2 / 389.379
Mp = 938.272
rows = [  # nu MeV, Re f (fit I), err, Re f (fit II), err, Re g, err
 (50, -2.84,.01, -2.85,.01, -0.26,.01),
 (60, -2.76,.01, -2.76,.01, -0.32,.01),
 (70, -2.65,.01, -2.66,.01, -0.37,.01),
 (79, -2.54,.01, -2.56,.01, -0.42,.01),
 (90, -2.38,.01, -2.40,.01, -0.48,.01),
 (100,-2.21,.01, -2.23,.01, -0.54,.01),
 (110,-2.00,.01, -2.03,.01, -0.59,.01),
 (120,-1.75,.01, -1.78,.01, -0.64,.01),
 (130,-1.45,.02, -1.48,.02, -0.68,.01),
 (135,-1.27,.02, -1.30,.02, -0.69,.01),
 (145,-0.79,.03, -0.83,.03, -0.65,.02),
 (150,-0.39,.03, -0.42,.03, -0.52,.02),
]
out = sys.argv[1] if len(sys.argv) > 1 else '.'
SOURCE = sys.argv[2] if len(sys.argv) > 2 else '2026'
assert SOURCE in ('2026', 'ghp2016'), SOURCE
RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'raw', 'ForwardRCS2026')
hdr = ['dataset','texkey','obs','E_lab_MeV','dE_MeV','theta_lab_deg','theta_cm_deg','frame_of_cs',
       'value','stat_err','sys_err_pt','norm_err_frac','units','source','notes']


def derived(f, g, ef, eg):
    """CS(0) [nb/sr, lab], Sigma_2z(0) and their errors for independent errors ef, eg of f, g"""
    cs = 1000*(f*f+g*g)/MUB_GEV2
    ecs = 1000*math.sqrt((2*f*ef)**2 + (2*g*eg)**2)/MUB_GEV2
    D = f*f+g*g; s2z = -2*f*g/D
    ds_df = (-2*g*D + 2*f*g*2*f)/D**2; ds_dg = (-2*f*D + 2*f*g*2*g)/D**2
    return cs, ecs, s2z, math.hypot(ds_df*ef, ds_dg*eg)


amp, obs = [], []
if SOURCE == 'ghp2016':
    src = 'arXiv:1604.00789 Table (tab:values_forward_amps)'
    for nu, f1, e1, f2, e2, g, eg in rows:
        f = 0.5*(f1+f2); ef = math.sqrt((0.5*(e1+e2))**2 + (0.5*(f1-f2))**2)
        amp.append(['GHP_forward','Gryniuk:2015eza','FWD_REF',nu,'',0,0,'na',f'{f:.4f}',f'{ef:.4f}',0,0,'mub*GeV',src,
                    f'Re f, mean of fit I ({f1}) and fit II ({f2}); err = quoted (+) half-difference'])
        amp.append(['GHP_forward','Gryniuk:2016gnm','FWD_REG',nu,'',0,0,'na',f'{g:.4f}',f'{eg:.4f}',0,0,'mub*GeV',src,'Re g'])
        cs, ecs, s2z, es2z = derived(f, g, ef, eg)
        obs.append(['GHP_forward','Gryniuk:2016gnm','CS',nu,'',0,0,'lab',f'{cs:.4f}',f'{ecs:.4f}',0,0,'nb/sr',src,
                    'derived: |f|^2+|g|^2 ; errors of different energies strongly correlated (same sigma_abs fits)'])
        obs.append(['GHP_forward','Gryniuk:2016gnm','S2Z',nu,'',0,0,'na',f'{s2z:.5f}',f'{es2z:.5f}',0,0,'1',src,
                    'derived: -2 f g/(f^2+g^2); Sigma_2z sign convention of arXiv:1604.00789 (3/2 - 1/2)'])
    names = ('db_forward_amp_ghp2016.csv', 'db_forward_obs_ghp2016.csv')
else:
    import numpy as np
    J = json.load(open(os.path.join(RAW, 'forward_points.json')))
    E = np.array(J['E'], float); F = np.array(J['re_f']); G = np.array(J['re_g'])
    assert list(E) == [r[0] for r in rows], E

    def split(cov, sys_, frac):
        w, v = np.linalg.eigh(np.array(cov))
        assert abs(w[-1] / w.sum() - frac) < 1e-6, (w[-1] / w.sum(), frac)      # = modes_*[1] of the json
        lead = v[:, -1] * np.sqrt(w[-1])
        lead *= np.sign(lead.sum())                                               # one sign at all energies
        assert np.all(lead > 0), 'leading mode changes sign'
        tot = np.sqrt(np.diag(cov) + sys_ ** 2)
        coh = np.sqrt(lead ** 2 + sys_ ** 2)
        unc = np.sqrt(np.clip(tot ** 2 - coh ** 2, 0, None))
        # same energy shape of mode and systematic (justifies one nuisance per amplitude): |sys|/lead nearly constant
        r = np.abs(sys_) / lead
        return tot, coh, unc, r
    tf, cf, uf, rf = split(J['cov_f'], np.array(J['sel_f']), J['modes_f'][1])
    tg, cg, ug, rg = split(J['cov_g'], np.array(J['pwa_g']), J['modes_g'][1])
    assert np.allclose(np.sqrt(np.diag(J['cov_f'])), J['err_f']) and np.allclose(np.sqrt(np.diag(J['cov_g'])), J['err_g'])
    src = 'ForwardRCS 2026 update (github.com/HiWARPs/ForwardRCS; raw/ForwardRCS2026/forward_points.json)'
    ntf = ('Re f (2026 global fit of sigma_tot with normalization nuisances); stat_err = total = exp+fit (+) data-selection '
           'systematic; err_unc = uncorrelated remainder after the leading exp+fit eigenmode (+) selection mode')
    ntg = ('Re g (2026 polarized fit, A2 2026 data, Bass et al. 2018 Regge); stat_err = total = exp+fit (+) PWA-model '
           'systematic (< 0.2 GeV); err_unc = uncorrelated remainder after the leading exp+fit eigenmode (+) PWA mode')
    for i, nu in enumerate(J['E']):
        f, g = float(F[i]), float(G[i])
        amp.append(['ForwardRCS2026', '', 'FWD_REF', nu, '', 0, 0, 'na', repr(f), repr(float(tf[i])), 0, 0, 'mub*GeV', src, ntf,
                    repr(float(uf[i]))])
        amp.append(['ForwardRCS2026', '', 'FWD_REG', nu, '', 0, 0, 'na', repr(g), repr(float(tg[i])), 0, 0, 'mub*GeV', src, ntg,
                    repr(float(ug[i]))])
        cs, ecs, s2z, es2z = derived(f, g, tf[i], tg[i])
        _, ucs, _, us2z = derived(f, g, uf[i], ug[i])
        obs.append(['ForwardRCS2026', '', 'CS', nu, '', 0, 0, 'lab', repr(cs), repr(float(ecs)), 0, 0, 'nb/sr', src,
                    'derived: |f|^2+|g|^2 from Re f, Re g; error = total of f (+) total of g; err_unc = uncorrelated part; '
                    'errors of different energies strongly correlated (one nuisance per amplitude in the fits)', repr(float(ucs))])
        obs.append(['ForwardRCS2026', '', 'S2Z', nu, '', 0, 0, 'na', repr(s2z), repr(float(es2z)), 0, 0, '1', src,
                    'derived: -2 f g/(f^2+g^2); Sigma_2z convention sigma_3/2 - sigma_1/2 (as GHP16); err_unc = uncorrelated part',
                    repr(float(us2z))])
    # self-checks against the tabulated observables of the repository (obs.txt, 5-MeV grid; exp+fit errors only);
    # above the pi0 p threshold (144.7 MeV) obs.txt includes Im f, Im g, neglected here: 4e-5 (relative) at 150 MeV
    ot = {round(1000 * float(l.split()[0])): [float(x) for x in l.split()[1:]] for l in open(os.path.join(RAW, 'obs.txt'))
          if l.strip() and l.split()[0][0].isdigit()}
    for i, nu in enumerate(J['E']):
        if nu in ot:
            cs, _, s2z, _ = derived(F[i], G[i], 0, 0)
            tol = 1e-4 if nu > 144.69 else 2e-6
            assert abs(cs / ot[nu][0] - 1) < tol and abs(s2z - ot[nu][4]) < tol, (nu, cs, ot[nu][0], s2z, ot[nu][4])
    print('E     tot_f   coh_f   unc_f  |sel|/mode   tot_g   coh_g   unc_g  |pwa|/mode   [mub GeV]')
    for i, nu in enumerate(J['E']):
        print('%4d  %.5f %.5f %.5f  %.3f   %.5f %.5f %.5f  %.3f' % (nu, tf[i], cf[i], uf[i], rf[i], tg[i], cg[i], ug[i], rg[i]))
    hdr = hdr + ['err_unc']
    names = ('db_forward_amp.csv', 'db_forward_obs.csv')
for name, data in zip(names, (amp, obs)):
    with open(os.path.join(out, name), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(hdr); w.writerows(data)
    print('wrote', name, len(data))
