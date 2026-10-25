"""t_dpole.py -- Delta-pole tree graphs from the Feynman rules of delta.py vs s_Delta2 of the code (8 structures).
Fixes the sign SEPS of eps^{mu nu rho lam} = SEPS i gamma^[mu nu rho lam] gamma5 through the gM gE interference."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'py'))
import formtools as ft, kin, devalamp as dv, delta
from dgraphs import DTok, DGraph
from graphs import P_, Q_, QP_
from sdelta2_poly import sdelta2_poly

S8 = [ft.compile_proj(ft.parse_fortran(open(os.path.join(os.path.dirname(__file__), '..', 'out', f'gram{k+1}.log')).read()))
      for k in range(8)]

def vtx(c, p_delta, idx, k, e):
    if c == 'gM': return DTok('EPS4', p_delta, idx, k, e, g5=True), 1.5 * delta.SEPS * 1j
    return DTok('SCI', p_delta, k, e, idx, g5=True), 1.5j

def pole_graphs(c1, c2, rs='RS'):
    """c1 at the incoming-photon vertex A, c2 at the outgoing B; s- and u-channel"""
    gs = []
    # s-channel: N(p) -A-> Delta(p+q) -B-> N(pp)
    pd = P_ + Q_
    vB, fB = vtx(c2, pd, 'a1', -QP_, 'ep'); vA, fA = vtx(c1, pd, 'a2', Q_, 'ee')
    g = DGraph(f'dpole_s_{c1}{c2}_{rs}'); g.toks = [vB, DTok(rs, pd, 'a1', 'a2'), vA]; g.tree = ['ixD']
    g.pref = fA * fB * 1j * (-1j); gs.append(g)
    # u-channel: N(p) -B-> Delta(p-qp) -A-> N(pp)
    pd = P_ - QP_
    vA, fA = vtx(c1, pd, 'a1', Q_, 'ee'); vB, fB = vtx(c2, pd, 'a2', -QP_, 'ep')
    g = DGraph(f'dpole_u_{c1}{c2}_{rs}'); g.toks = [vA, DTok(rs, pd, 'a1', 'a2'), vB]; g.tree = ['iyD']
    g.pref = fA * fB * 1j * (-1j); gs.append(g)
    for g in gs:
        g.coup = {c1: 1, c2: 1, 'iMp': 2} if c1 != c2 else {c1: 2, 'iMp': 2}
        g.iso = {'1': (delta.T[2] @ delta.Td[2])[0, 0]}
    return gs

def code_amp(x, t, cs, GM, GE):
    y = -x - t / 2
    Md = cs['MD']
    def one(xx):
        s = 1 + 2 * xx
        G1 = 1.5 / (1 + Md) * GM; G2 = 1.5 / (1 + Md) * GE
        return -(2 / 3) / (s - Md ** 2) * np.array(sdelta2_poly(xx, t, G1, G2, Md))
    As, Au = one(x), one(y)
    A = As.copy(); A[[0, 1, 7]] += Au[[0, 1, 7]]; A[2:7] -= Au[2:7]
    return A

if __name__ == '__main__':
    cs = dv.dconsts()
    for seps in (1, -1):
        delta.SEPS = seps
        allg = []
        for c1 in ('gM', 'gE'):
            for c2 in ('gM', 'gE'):
                allg += pole_graphs(c1, c2)
        allg += pole_graphs('gM', 'gM', 'RSR') + pole_graphs('gE', 'gM', 'RSR')
        for side in 'LR':
            for gg in pole_graphs('gM', 'gE', 'RSH'):
                gg.toks[1] = DTok('RSH', *gg.toks[1].args, side); gg.name += side; allg.append(gg)
        for g in allg: g.name += f'_s{seps}'
        comp = dv.compile_graphs(allg, nproc=8)
        print(f'--- SEPS = {seps}')
        for (x, t) in [(0.05, -0.01), (0.12, -0.08), (0.2, -0.2)]:
            en = dv.env(x, t, cs)
            en8 = dict(en, T={(0, 0, 0, 0): 1.0}, U={})
            G8 = np.array([[eval(S8[k][j], en8) for k in range(8)] for j in range(6)])
            for (gm, ge) in [(cs['gM'], 0.), (0., cs['gE']), (cs['gM'], cs['gE'])]:
                cc = dict(cs, gM=gm, gE=ge)
                V = sum(c.raw(en) * c.coup(cc) * c.isof(cc) for n, c in comp.items() if n.split('_')[3] == 'RS')
                Vc = G8 @ code_amp(x, t, cs, gm, ge)
                print(f'x={x} t={t} gM={gm:.2f} gE={ge:.2f}: rel dev {np.max(abs(V - Vc)) / np.max(abs(Vc)):.2e}')
            V1 = sum(c.raw(en) * c.coup(cs) * c.isof(cs) for n, c in comp.items() if 'RSR' in n)
            V2 = sum(c.raw(en) * c.coup(cs) * c.isof(cs) for n, c in comp.items()
                     if n.split('_')[3] == 'RS' and ('gMgM' in n or 'gEgM' in n))
            print(f'   RSR vs RS (gMgM + gEgM): {np.max(abs(V1 - V2)) / np.max(abs(V2)):.1e}')
            V3 = sum(c.raw(en) * c.coup(cs) * c.isof(cs) for n, c in comp.items() if n.split('_')[3] == 'RS' and 'gMgE' in n)
            for side in 'LR':
                V4 = sum(c.raw(en) * c.coup(cs) * c.isof(cs) for n, c in comp.items() if n.split('_')[3] == 'RSH' + side)
                print(f'   RSH{side} vs RS (gMgE): {np.max(abs(V3 - V4)) / np.max(abs(V3)):.1e}')
