#!/usr/bin/env python3
"""
Absorptive parts of the one-loop (piN) gamma N Delta vertex corrections,
Pascalutsa & Vanderhaeghen, PRD 73 (2006) 034003, Appendix B ("Absorptive parts of the
vertex corrections"), Eqs. (gipi) [photon on the pion, graph a] and (giN) [photon on the
nucleon, graph b], transcribed literally and evaluated for real photons by taking Q^2 -> 0
numerically in high precision (mpmath; the individual coefficients contain 1/Q^2 terms that
cancel).  Mapping to the Lagrangian couplings (PV2006 Eq. gmecgi, Q^2 = 0):
   g_M = -(2/3) M_N (M_D + M_N)/M_D * g_1 ,   g_E = (2/3) M_N (M_D + M_N)/M_D * (g_1 - M_D g_2)
(optionally with M_D -> W off shell).  Used only as an independent cross-check of the
unitarity (Watson) construction implemented in As.f90 (SUBROUTINE vertex_im).
"""
import mpmath as mp

mp.mp.dps = 60


def im_g12(W, MN=0.93827, mpi=0.13957, MD=1.232, gA=1.27, hA=2.85, fpi=0.09221, Q2t=mp.mpf('1e-30')):
    W = mp.mpf(W)
    rt = MN / W; mut = mpi / W; rD = MD / W
    lam = mp.sqrt(((1 + rt)**2 - mut**2) * ((1 - rt)**2 - mut**2)) / 2
    bet = (1 - rt**2 + mut**2) / 2
    alp = 1 - bet
    kap = mp.sqrt(((1 + rt)**2 + Q2t) * ((1 - rt)**2 + Q2t)) / 2
    bg = (1 - rt**2 - Q2t) / 2
    wpi = 2 * kap * lam / (Q2t + 2 * bet * bg)
    wN = 2 * kap * lam / (Q2t + 2 * alp * bg)
    Lp = mp.log((1 + wpi) / (1 - wpi))
    LN = mp.log((1 + wN) / (1 - wN))
    f = 1 / (4 * mp.pi)
    Q = Q2t
    # ---- pion-photon graph (a)
    c1p_pi = f / (4 * kap**2) * (lam * ((1 - 3 * Q / kap**2) * (2 * bet * bg + Q) + 2 * Q)
                                 + 1 / (2 * kap) * Lp * (lam**2 * (2 * bg**2 * (1 - 1 / wpi**2) + 4 * Q / wpi**2)
                                                         - Q * (2 * bet * bg + Q)))
    d1p_pi = f * bg / (8 * kap**4) * (lam * (2 * lam**2 * (mp.mpf(2) / 3 * (bg**2 - Q) + 5 * Q / wpi**2)
                                             - Q * (4 * bet**2 + 2 * Q + 2 * bet * bg))
                                      + Q / (2 * kap) * Lp * (lam**2 * (bg * (2 * bet - bg) * (3 - 5 / wpi**2)
                                                                        + kap**2 * (1 - 1 / wpi**2))
                                                              + bet * (2 * bet - bg) * (2 * bet * bg + Q)))
    d2_pi = f * lam**3 / (2 * kap**2) * (mp.mpf(2) / 3 - 1 / wpi**2 + 1 / (2 * wpi) * Lp * (-1 + 1 / wpi**2))
    d2p_pi = f / (8 * kap**4) * (lam * (2 * lam**2 * (bg**2 * (-mp.mpf(4) / 3 + 2 / wpi**2) - 3 * Q / wpi**2)
                                        + mp.mpf(3) / 2 * Q * (2 * bet * bg + Q))
                                 + lam**2 / (2 * kap) * Lp * (4 * kap**2 * bet * bg * (1 - 1 / wpi**2)
                                                              - bg * Q * (2 * bet - bg) * (3 - 5 / wpi**2)))
    # ---- nucleon-photon graph (b)
    a1N = f * lam * alp
    a2N = f * (-lam**3 / 3)
    a3N = f * lam * (rt**2 + 4 * lam**2 / 3)
    b1N = f / (2 * kap) * LN
    b2N = f / (4 * kap**2) * (4 * lam * bg + Q * (1 + rt**2 + Q + 2 * (rt**2 - mut**2)) / (2 * kap) * LN)
    b3N = f / (2 * kap**2) * (-2 * lam + (2 * alp * bg + Q) / (2 * kap) * LN)
    c1N = 1 / (4 * kap**2) * (-2 * bg * a1N + 2 * (bg**2 * rt**2 - Q * lam**2) * b1N + Q * bg * b2N
                              + Q * (kap**2 + bg**2) * b3N)
    c3N = 1 / kap**2 * (3 * c1N + lam**2 * b1N)
    c1pN = -(1 - 3 * Q / kap**2) * c1N + Q * lam**2 / kap**2 * b1N - Q / 2 * b3N
    c2N = -2 * c1N - Q * c3N + rt**2 * b1N + Q * b3N
    d1N = 1 / (2 * kap**2) * (-bg * (2 * a2N + a3N) + (rt**2 * (kap**2 - Q) - Q * lam**2) * b2N
                              + Q / (2 * bg) * (2 * kap**2 - Q) * c1pN + Q / 2 * bg * (c2N + 2 * c1N))
    d1pN = (-mp.mpf(5) / 2 * d1N + rt**2 / 2 * b2N - (1 - rt**2) / (4 * kap**2) * a3N
            + Q / (2 * kap**2) * (5 * d1N + a3N / 2 - rt**2 * b2N + bg / 2 * (c1pN + c2N)))
    d2N = bg * d1N / Q + c1N / 2 - a2N / Q
    d2pN = bg * (d1N + 2 * d1pN) / (2 * Q) + c1pN / 4
    K = gA * hA * MN / (4 * fpi**2)
    g1a = K * d2_pi
    g2a = -K / W * ((rD + rt) / bg * (2 * d2p_pi - c1p_pi) + 2 * rD / bg * d1p_pi)
    g1b = -K * (d2N + rt * (rD - rt) / 2 * b3N + (1 - rt * rD) / (2 * bg) * c1pN - Q / 2 * c3N)
    g2b = K / W * (-rt * b3N - (rD - rt) / bg * c1pN + 2 * rD / bg * d1pN + 2 * (rD + rt) / bg * d2pN)
    return dict(g1a=g1a, g2a=g2a, g1b=g1b, g2b=g2b, g1=g1a + g1b, g2=g2a + g2b)


def im_gME(W, MN=0.93827, MD=1.232, offshell=False, **kw):
    r = im_g12(W, MN=MN, MD=MD, **kw)
    Md = mp.mpf(W) if offshell else mp.mpf(MD)
    c = mp.mpf(2) / 3 * MN * (Md + MN) / Md
    return float(-c * r['g1']), float(c * (r['g1'] - Md * r['g2'])), r


if __name__ == '__main__':
    for Q2 in ['1e-10', '1e-20', '1e-30']:
        gm, ge, _ = im_gME(1.232, Q2t=mp.mpf(Q2))
        print('Q2t=%s  Im gM=%.10f  Im gE=%.10f' % (Q2, gm, ge))
    print('  W     Im gM    Im gE   (on-shell mapping)   Im gM  Im gE (M_D->W)')
    for W in [1.10, 1.15, 1.20, 1.232, 1.26, 1.30, 1.35]:
        a = im_gME(W); b = im_gME(W, offshell=True)
        print('%6.3f %8.4f %8.4f      %8.4f %8.4f' % (W, a[0], a[1], b[0], b[1]))
