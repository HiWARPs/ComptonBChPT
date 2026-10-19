"""Table: the parameter-free NNLO prediction (P = Pref, no free polarisability) versus each data set below threshold
(E_lab <= 150 MeV, minus MAIN_DROP; the HIgS 2022 sets, which are NOT part of the main database since 2026-10-17, are
listed as well, as a test of the prediction) and the forward points (E <= FWD_EMAX, nuisance model).
Columns per set: N; chi2 with n = 1; chi2 with the normalisation floated within its uncertainty (data + penalty);
the same with the truncation covariance (Appendix C, default settings); fitted n in the last case.
Also: global values for the MAIN database (pre-2015 + A2 2022 minus MAIN_DROP, + forward points; keys GLOBAL,
GLOBAL_angular) and for the main database + HIgS (keys GLOBAL_withHIGS, GLOBAL_angular_withHIGS) with all
normalisations/nuisances floated, with and without the truncation covariance; GLOBAL_allpts (main database with the
3-sigma outlier), GLOBAL_pre2015 (data prior to 2015), GLOBAL_A2only and GLOBAL_HIGSonly (the subset fits), each with the
forward points, for the fits table.
Writes scenarios/nnlo_prediction_table.{txt,json,tex}"""
import os, sys, json, numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..', 'lib')); os.chdir(os.path.join(here, '..'))
from common import *  # noqa

FWD_EMAX = 145.0 if PROFILE == 'v2' else 135.0
df = load_db(); T = outputs(df, 'NNLO')
angular_sets = [d for d in df.dataset.unique() if d not in FWD]
sLH = select(df, angular_sets, Emax=150.0)                 # all angular sets incl. HIgS (minus MAIN_DROP)
sLH = sLH[~sLH.pid.isin(MAIN_DROP)]
sL = sLH[~sLH.dataset.isin(HIGS)]                          # MAIN database (angular part): no HIgS
assert set(sL.dataset) == set(MAIN_L_SETS), set(sL.dataset) ^ set(MAIN_L_SETS)
sF = select(df, FWD, Emax=FWD_EMAX)


def evaluate(s, thcov=False):
    """prediction vs rows s: returns dict with chi2 at n=1 (shifts 0), chi2 with nuisances floated, fitted n."""
    C_th = truncation_for(df, s, tag='db')[0] if thcov else None
    F = build_fit(s, T, PREF_NNLO, [], thcov=C_th, fwd_mode='nuisance')
    nq = len(F.norm_sets) + len(F.shifts)
    R0 = F.summarize(np.zeros(nq), np.zeros((nq, nq)))          # n = 1, shifts = 0
    if nq:
        R = F.run(decompose=False)
    else:
        R = R0
    return dict(N=int(F.N), chi2_fixed=float(R0.chi2_data), chi2_float=float(R.chi2), chi2_float_data=float(R.chi2_data),
                n={k: (float(v['n']), float(v['dn'])) for k, v in R.perset.items()})


out = {}
for ds, g in sLH.groupby('dataset'):                       # per set: main-database sets AND the HIgS test sets
    a = evaluate(g, thcov=False); b = evaluate(g, thcov=True)
    out[ds] = dict(N=a['N'], fixed=a['chi2_fixed'], float=a['chi2_float'], th=b['chi2_float'],
                   n_float=a['n'][ds][0], n_th=b['n'][ds][0], dn_th=b['n'][ds][1], in_main=ds not in HIGS)
# forward points together (nuisance shifts common to both forward sets)
a = evaluate(sF, thcov=False); b = evaluate(sF, thcov=True)
out['GHP_forward'] = dict(N=a['N'], fixed=a['chi2_fixed'], float=a['chi2_float'], th=b['chi2_float'])
# global: main database (no HIgS) and main database + HIgS
for sfx, sLx in (('', sL), ('_withHIGS', sLH)):
    sAll = pd.concat([sLx, sF])
    A = evaluate(sAll, thcov=False); B = evaluate(sAll, thcov=True)
    out['GLOBAL' + sfx] = dict(N=A['N'], fixed=A['chi2_fixed'], float=A['chi2_float'], th=B['chi2_float'], th_data=B['chi2_float_data'])
    AL = evaluate(sLx, thcov=False); BL = evaluate(sLx, thcov=True)
    out['GLOBAL_angular' + sfx] = dict(N=AL['N'], fixed=AL['chi2_fixed'], float=AL['chi2_float'], th=BL['chi2_float'])
# the database variants of the fits table: all points incl. the 3-sigma outlier, and the data prior to 2015 only (each with the forward points)
sAP = select(df, list(MAIN_L_SETS), Emax=150.0)
sOL = select(df, list(OLD_L), Emax=150.0); sOL = sOL[~sOL.pid.isin(MAIN_DROP)]
sNL = select(df, [k for k in MAIN_L_SETS if k not in ('Pugh1957', 'Oxley1958', 'Hyman1959', 'Goldansky1960', 'Bernardini1960')], Emax=150.0)
sNL = sNL[~sNL.pid.isin(MAIN_DROP)]                        # main database without the bremsstrahlung data of 1957-1960
sA2 = select(df, A2, Emax=150.0); sHI = select(df, HIGS, Emax=150.0)   # the subset fits (A2 alone, HIgS alone)
for key, sx in (('GLOBAL_allpts', sAP), ('GLOBAL_pre2015', sOL), ('GLOBAL_noLegacy', sNL), ('GLOBAL_A2only', sA2),
                ('GLOBAL_HIGSonly', sHI)):
    A = evaluate(pd.concat([sx, sF]), thcov=False); B = evaluate(pd.concat([sx, sF]), thcov=True)
    out[key] = dict(N=A['N'], fixed=A['chi2_fixed'], float=A['chi2_float'], th=B['chi2_float'], th_data=B['chi2_float_data'])

lines = ['Pref = ' + ' '.join('%.3f' % x for x in PREF_NNLO),
         'main database = pre-2015 + A2 2022 minus %s (HIgS 2022 excluded since 2026-10-17, listed as a test); forward points E <= %.0f MeV'
         % (list(MAIN_DROP), FWD_EMAX),
         '%-26s %4s %9s %9s %9s %8s' % ('set', 'N', 'chi2(n=1)', 'n float', 'n+theory', 'n_th')]
for k, v in out.items():
    lines.append('%-26s %4d %9.1f %9.1f %9.1f %8s' % (k, v['N'], v['fixed'], v['float'], v['th'],
                                                     ('%.3f' % v['n_th']) if 'n_th' in v else ''))
open('scenarios/nnlo_prediction_table.txt', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
json.dump(out, open('scenarios/nnlo_prediction_table.json', 'w'), indent=1)

# LaTeX rows for the paper (grouped)
GROUPS = [('A2 (2022), $\\dd\\si/\\dd\\W$~\\cite{A2:2021wmm}', ['A2_2022_CS']),
          ('A2 (2022), $\\Si_3$~\\cite{A2:2021wmm}', ['A2_2022_S3']),
          ('HI$\\ga$S (2022), $\\dd\\si/\\dd\\W$~\\cite{Li:2022vnz}', ['HIGS_2022_CIRC']),
          ('HI$\\ga$S (2022), $\\dd\\si_\\parallel$, $\\dd\\si_\\perp$~\\cite{Li:2022vnz}', ['HIGS_2022_LIN']),
          ('TAPS (2001)~\\cite{deLeon:2001dnx}', ['OlmosdeLeon2001']),
          ('Hallin \\etal~\\cite{Hallin:1993ft}', ['Hallin1993_TabII']),
          ('MacGibbon \\etal~\\cite{MacGibbon:1995in}', ['MacGibbon1995_tagged', 'MacGibbon1995_untagged']),
          ('Federspiel \\etal~\\cite{Federspiel:1991yd}', ['Federspiel1991']),
          ('Zieger \\etal~\\cite{Zieger:1992jq}', ['Zieger1992']),
          ('Baranov \\etal~\\cite{Baranov:1974ec,Baranov:1975ju}', ['Baranov1974_75_th90', 'Baranov1974_75_th150']),
          ('Bremsstrahlung data, 1957--1960~\\cite{Pugh:1957zz,Oxley:1958zz,Hyman:1959zz,Goldansky:1960zz,Bernardini:1960wya}',
           ['Pugh1957', 'Oxley1958', 'Hyman1959', 'Goldansky1960', 'Bernardini1960']),
          ('Forward points, \\Eqref{forward}', ['GHP_forward'])]
tex = ['% main database = pre-2015 + A2 2022 (minus MAIN_DROP) + forward points; the HIgS rows are a test (not in the main database since 2026-10-17)']
for lab, keys in GROUPS:
    N = sum(out[k]['N'] for k in keys); f = sum(out[k]['fixed'] for k in keys)
    fl = sum(out[k]['float'] for k in keys); th = sum(out[k]['th'] for k in keys)
    tex.append('%s & %d & %.1f & %.1f & %.1f \\\\' % (lab, N, f, fl, th))
tex.append('\\hline')
for lab, key in (('All angular data (main database)', 'GLOBAL_angular'), ('All data (main database)', 'GLOBAL'),
                 ('All angular data incl.\\ HI$\\ga$S', 'GLOBAL_angular_withHIGS'), ('All data incl.\\ HI$\\ga$S', 'GLOBAL_withHIGS')):
    g = out[key]
    tex.append('%s & %d & %.1f & %.1f & %.1f \\\\' % (lab, g['N'], g['fixed'], g['float'], g['th']))
open('scenarios/nnlo_prediction_table.tex', 'w').write('\n'.join(tex) + '\n'); print('\n'.join(tex))
