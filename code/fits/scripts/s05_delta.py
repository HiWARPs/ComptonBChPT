"""Step 5: Delta region D (150 < E_lab <= 350 MeV).
(i)  chi^2 per data set of the parameter-free NNLO prediction and of the F1/F2 best fits (alpha, beta from region L).
(ii) exploratory fit of the four spin polarisabilities to Sigma_2x (Martel) + Sigma_2z (Paudyal) (+ LEGS Sigma_3
     below the 2pi threshold), with gamma_0 = -0.929(105) [GHP], alpha+beta = 14.0(2) and alpha-beta from F1;
     with and without floating gM (interpolation over a grid of complete NNLO runs).
     the 2015 code lacks the gamma-N-Delta vertex corrections (LMP2015 Sect. 2.2): Sigma_2x off by 0.03-0.06 at
     286-334 MeV; a 'vertex' theory error (0.03 Sigma_2x, 0.02 Sigma_2z, 0.01 Sigma_3) is added in a variant.
     v2 (profile 'v2'): vertex corrections included (dvertex=t); the same fits are repeated WITHOUT them ('_novtx',
     dvertex=f) to isolate their effect.
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.stats import chi2 as chi2dist
from common import *  # noqa
import bchpt
from s02_run_driver import GM_GRID

df = load_db()
T = outputs(df, 'NNLO')
HAVE_NV = (PROFILE == 'v2')
TX = {True: T, False: outputs(df, 'novtx') if HAVE_NV else None}
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
MAIN = {'F1': FL['F1_d_all_thB'], 'F2': FL['F2_d_all_thB'], 'F1_noth': FL['F1_d_all'], 'F2_noth': FL['F2_d_all']}
D_ALL = ['A2_Martel2015_S2x', 'A2_Paudyal2020_S2z', 'LEGS2001_S3', 'LEGS2001_CS', 'Hunger1997_75', 'Hunger1997_90',
         'Wissmann1999', 'Camen2002', 'Hallin1993_TabII', 'Hallin1993_TabII_D', 'OlmosdeLeon2001', 'Gray1967',
         'Genzel1972', 'Baranov1966', 'Kolb2000', 'Nagashima1964']
TH_VERTEX = {'S2X': 0.03, 'S2Z': 0.02, 'S3': 0.01}
GMT = {True: {g: outputs(df, 'gM%.2f' % g) for g in GM_GRID},
       False: {g: outputs(df, 'gM%.2fnv' % g) for g in GM_GRID} if HAVE_NV else None}
OUT = {}


def chi2_table(P, label, Tx=None):
    Tx = T if Tx is None else Tx
    s = df[df.dataset.isin(D_ALL) & (df.E > 150) & (df.E <= 350)].copy()
    s.loc[s.dataset == 'Hallin1993_TabII', 'dataset'] = 'Hallin1993_TabII_D'
    idx = Tx.subset(s.pid.values)
    th = Tx.predict_array(P, idx)
    rows = []
    for ds in dict.fromkeys(s.dataset):
        m = (s.dataset == ds).values
        y, e, t = s.value.values[m], s.err_pt.values[m], th[m]
        c0 = float(np.sum(((y - t) / e) ** 2))
        dn = float(s[m].norm.iloc[0])
        if dn > 0:      # optimal normalisation with penalty (analytic: linear in n)
            a = np.sum(t * t / e ** 2) + 1 / dn ** 2
            b = np.sum(t * y / e ** 2) + 1 / dn ** 2
            n = b / a
            c1 = float(np.sum(((y - n * t) / e) ** 2) + ((n - 1) / dn) ** 2)
        else:
            n, c1 = 1.0, c0
        rows.append(dict(set=ds, obs='/'.join(sorted(set(s[m].dobs))), N=int(m.sum()), Erange='%.0f-%.0f' % (s[m].E.min(), s[m].E.max()),
                         chi2_n1=c0, chi2_nfit=c1, n=n, chi2N=c1 / m.sum()))
        if ds == 'LEGS2001_S3':
            mb = m & (s.E.values < E2PI)
            rows.append(dict(set='LEGS2001_S3 (E<2pi thr.)', obs='S3', N=int(mb.sum()), Erange='<%.0f' % E2PI,
                             chi2_n1=float(np.sum(((s.value.values[mb] - th[mb]) / s.err_pt.values[mb]) ** 2)),
                             chi2_nfit=float(np.sum(((s.value.values[mb] - th[mb]) / s.err_pt.values[mb]) ** 2)), n=1.0,
                             chi2N=float(np.sum(((s.value.values[mb] - th[mb]) / s.err_pt.values[mb]) ** 2)) / mb.sum()))
    tab = pd.DataFrame(rows)
    OUT['chi2_' + label] = tab.to_dict('records')
    return tab


E2PI = 320.7


def gm_theory(idx, vtx=True):
    """observable at total polarisabilities P and coupling gM by cubic-spline interpolation over the gM grid"""
    grid = np.array(GM_GRID)
    G = GMT[vtx]

    def f(P, ex):
        O = np.array([G[g].predict_array(P, idx) for g in GM_GRID])
        return CubicSpline(grid, O, axis=0)(ex['gM'])
    return f


def spin_fit(name, sets, ab_from='F2', float_gM=False, vertex_err=False, gpi=None, legs_Emax=E2PI, free=None,
             gM_prior=(2.97, 0.10), vtx=True):
    m = df.dataset.isin(sets) & (df.E > 150) & (df.E <= 350)
    if 'LEGS2001_S3' in sets:
        m &= ~((df.dataset == 'LEGS2001_S3') & (df.E >= legs_Emax))
    s = df[m].copy()
    if vertex_err:
        s['err_pt'] = np.sqrt(s.err_pt ** 2 + np.array([TH_VERTEX.get(o, 0) for o in s.dobs]) ** 2)
    F1 = MAIN[ab_from]
    P0 = np.array(F1['P']).copy()
    amb, damb = P0[0] - P0[1], np.sqrt(F1['covP'][0][0] + F1['covP'][1][1] - 2 * F1['covP'][0][1])
    cons = [(C_BALDIN, BALDIN[0], BALDIN[1], 'Baldin'), (np.array([1, -1, 0, 0, 0, 0.]), amb, damb, 'a-b(F1)'),
            (C_GAMMA0, GHP_G0[0], GHP_G0[1], 'gamma0')]
    if gpi is not None:
        cons.append((np.array([0, 0, -1, 1, -1, 1.]), gpi[0], gpi[1], 'gammapi'))
    free = free or ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']
    Tx = TX[vtx]
    idx = Tx.subset(s.pid.values)
    tf = gm_theory(idx, vtx) if float_gM else None
    F = build_fit(s, Tx, P0, free, constraints=cons, theory_fn=tf,
                  extra=dict(gM=2.97) if float_gM else None,
                  extra_prior=dict(gM=gM_prior) if (float_gM and gM_prior is not None) else None)
    R = F.run(decompose=False)
    out = dict(name=name, sets=sets, N=len(s), P=R.P, dP=R.dP, covP=R.covP, chi2=R.chi2, dof=R.dof, p=R.pvalue,
               g0=gamma0(R.P), gpi=gammapi(R.P), perset=R.perset, extra=R.extra, pids=list(s.pid), model=R.model,
               chi2_data=R.chi2_data)
    J = np.array([0, 0, -1, 1, -1, 1.]); out['dgpi'] = float(np.sqrt(J @ R.covP @ J))
    J0 = np.array([0, 0, -1, -1, -1, -1.]); out['dg0'] = float(np.sqrt(J0 @ R.covP @ J0))
    OUT[name] = out
    gm = (' gM=%.3f(%.3f)' % R.extra['gM']) if float_gM else ''
    print('%-30s N=%3d gEE=%5.2f(%.2f) gMM=%5.2f(%.2f) gEM=%5.2f(%.2f) gME=%5.2f(%.2f) g0=%5.2f gpi=%5.2f(%.2f) a=%.2f b=%.2f%s chi2=%.1f/%d'
          % (name, len(s), R.P[2], R.dP[2], R.P[3], R.dP[3], R.P[4], R.dP[4], R.P[5], R.dP[5], gamma0(R.P), gammapi(R.P),
             out['dgpi'], R.P[0], R.P[1], gm, R.chi2, R.dof))
    return out


if __name__ == '__main__':
    pd.set_option('display.width', 200)
    tabs = [('NNLO', PREF_NNLO, None), ('F2', np.array(MAIN['F2']['P']), None), ('F1', np.array(MAIN['F1']['P']), None),
            ('F2_noth', np.array(MAIN['F2_noth']['P']), None)]
    if HAVE_NV:
        tabs += [('NNLO_novtx', PREF_NNLO, TX[False]), ('F2_novtx', np.array(MAIN['F2']['P']), TX[False])]
    for lab, P, Tx in tabs:
        tab = chi2_table(P, lab, Tx)
        print('=== chi2 per data set in the Delta region (150-350 MeV): %s ===' % lab)
        print(tab.round(3).to_string(index=False))
    # gM dependence of the chi^2 per data set (NNLO static pols of each gM, i.e. parameter-free predictions)
    rows = []
    for vtx in ([True, False] if HAVE_NV else [True]):
      for g in GM_GRID:
        Tg = GMT[vtx][g]
        s = df[df.dataset.isin(D_ALL) & (df.E > 150) & (df.E <= 350)]
        idx = Tg.subset(s.pid.values)
        th = Tg.predict_array(np.r_[np.array(MAIN['F2']['P'])[:2], Tg.Pref[2:]], idx)
        for ds in ['A2_Martel2015_S2x', 'A2_Paudyal2020_S2z', 'LEGS2001_S3', 'LEGS2001_CS', 'Hunger1997_75', 'Hunger1997_90',
                   'Wissmann1999', 'Camen2002']:
            mm = (s.dataset == ds).values
            rows.append(dict(gM=g, vtx=vtx, set=ds, chi2=float(np.sum(((s.value.values[mm] - th[mm]) / s.err_pt.values[mm]) ** 2)), N=int(mm.sum())))
    OUT['gM_scan'] = rows
    R_ = pd.DataFrame(rows)
    for vtx in R_.vtx.unique():
        print('--- gM scan, vertex corrections =', vtx)
        print(R_[R_.vtx == vtx].pivot(index='gM', columns='set', values='chi2').round(1).to_string())

    print('=== exploratory spin-polarisability fits (alpha-beta from main fit F2, alpha+beta = 14.0(2); gamma_0 = -0.929(105)) ===')
    A = ['A2_Martel2015_S2x', 'A2_Paudyal2020_S2z']
    AL = A + ['LEGS2001_S3']
    two = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1']
    GPI = (8.0, 1.8)      # MAMI dispersive gamma_pi (Camen 2002), as used by Martel 2015 / Paudyal 2020
    spin_fit('S_A2_2free', A, free=two)
    spin_fit('S_A2', A)                                   # 4 free: under-determined (shown for completeness)
    spin_fit('S_A2_gpi', A, gpi=GPI)
    spin_fit('S_A2_LEGS_2free', AL, free=two)
    spin_fit('S_A2_LEGS', AL)
    spin_fit('S_A2_LEGS_gpi', AL, gpi=GPI)
    for nm, st, kw in [('S_A2_2free', A, dict(free=two)), ('S_A2_gpi', A, dict(gpi=GPI)), ('S_A2_LEGS', AL, {}),
                       ('S_A2_LEGS_gpi', AL, dict(gpi=GPI))]:
        spin_fit(nm + '_vtx', st, vertex_err=True, **kw)
        spin_fit(nm + '_gM', st, float_gM=True, **kw)
        spin_fit(nm + '_gM_vtx', st, float_gM=True, vertex_err=True, **kw)
        spin_fit(nm + '_gMwide', st, float_gM=True, gM_prior=(2.97, 0.30), **kw)
    spin_fit('S_A2_LEGS_gpi_F1ab', AL, gpi=GPI, ab_from='F1')
    spin_fit('S_A2_LEGS_gpi_noth_ab', AL, gpi=GPI, ab_from='F2_noth')
    spin_fit('S_A2_LEGS_all_gpi', A + ['LEGS2001_S3'], gpi=GPI, legs_Emax=400.0)
    # gM alone (spin pols at NNLO; alpha, beta constrained as above), wide prior 2.97(0.5)
    W = (2.97, 0.5)
    spin_fit('gM_only_A2asym', A, float_gM=True, free=['alphaE1', 'betaM1'], gM_prior=W)
    spin_fit('gM_only_LEGS_S3', ['LEGS2001_S3'], float_gM=True, free=['alphaE1', 'betaM1'], gM_prior=W)
    spin_fit('gM_only_MAMIcs', ['Hunger1997_75', 'Hunger1997_90', 'Wissmann1999', 'Camen2002'], float_gM=True,
             free=['alphaE1', 'betaM1'], gM_prior=W)
    spin_fit('gM_only_LEGScs', ['LEGS2001_CS'], float_gM=True, free=['alphaE1', 'betaM1'], gM_prior=W)
    spin_fit('gM_only_Hallin', ['Hallin1993_TabII_D'], float_gM=True, free=['alphaE1', 'betaM1'], gM_prior=W)
    if HAVE_NV:
        print('=== same without gamma-N-Delta vertex corrections (dvertex=f) ===')
        spin_fit('S_A2_2free_novtx', A, free=two, vtx=False)
        spin_fit('S_A2_gpi_novtx', A, gpi=GPI, vtx=False)
        spin_fit('S_A2_LEGS_novtx', AL, vtx=False)
        spin_fit('S_A2_LEGS_gpi_novtx', AL, gpi=GPI, vtx=False)
        spin_fit('S_A2_LEGS_gpi_gMwide_novtx', AL, gpi=GPI, float_gM=True, gM_prior=(2.97, 0.30), vtx=False)
        for nm, st in [('gM_only_A2asym', A), ('gM_only_LEGS_S3', ['LEGS2001_S3']),
                       ('gM_only_MAMIcs', ['Hunger1997_75', 'Hunger1997_90', 'Wissmann1999', 'Camen2002']),
                       ('gM_only_LEGScs', ['LEGS2001_CS']), ('gM_only_Hallin', ['Hallin1993_TabII_D'])]:
            spin_fit(nm + '_novtx', st, float_gM=True, free=['alphaE1', 'betaM1'], gM_prior=W, vtx=False)
    save_json(OUT, 'delta_region.json')
