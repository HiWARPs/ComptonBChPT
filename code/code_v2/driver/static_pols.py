#!/usr/bin/env python3
"""Extract static dipole polarisabilities from static_As.out (ampl.exe &fitdriver staticmode=t).

Method: invariant amplitudes A_1..A_8 of the code (basis of Lensky-Pascalutsa 2010 / LMP 2015 Eq.(7))
-> L'vov amplitudes A^L_i via the matrix C^{A->A^L} of LMP 2015 App. A (Eq. A2),
-> a_i = A^L_i(nu->0, t=0) by a polynomial fit in nu^2,
-> polarisabilities via LMP 2015 Eq. (B3) [first six lines].
Normalisation: with the code's A_i (nucleon-mass units, alpha_em stripped) the L'vov relations must be
multiplied by 1/2; this is fixed independently by the forward (Thomson-ratio) identity
f(w) = -alpha_em/M + (alpha+beta) w^2, and checked on the code's own contact terms (component 6).
"""
import sys, numpy as np
smm, hbarc, aem = 0.93827, 0.19733, 1.0/137.036
inv = hbarc/smm
p2l = 1.0/(aem*inv**3*1e4)          # polza2LECs of ampl.f90
to_fm3 = 1.0/p2l                    # internal (M^-3/alpha_em) -> 1e-4 fm^3
to_fm4 = inv/p2l                    # internal (M^-4/alpha_em) -> 1e-4 fm^4
NUMIN = 0.019                       # lowest nu (M units) used: LoopTools loses precision below ~0.0125
NUMAX = 0.121                       # highest nu (pi threshold at nu=0.149)

def CAL(nu, t, M=1.0):
    eta = t - 4*(nu - t/4)*(-nu - t/4)            # (M^4 - s u)/M^2 with M=1
    R = [[-nu*(4*M**2-t)/(4*M**3), 0, -nu/M, 0, 0, 0],
         [-eta*nu/(2*M), 0, -nu*t/(2*M), 0, 0, 0],
         [nu**2/M**2, 0, 1, 1, (t-4*M**2)/(4*M**2), (t-4*M**2)/(4*M**2)],
         [0, -nu**2*(4*M**2-t)/(2*M**2), 0.5*(4*nu**2+t), t/2, 0.5*(4*nu**2-eta), -eta/2],
         [nu**2*t/(2*M**2), -nu**2*(4*M**2-t)/(2*M**2), 2*nu**2+t, t, 2*nu**2-eta, 2*nu**2-eta],
         [-nu**2*t/(2*M**2), nu**2*(4*M**2-t)/(2*M**2), 2*nu**2, 0, -2*nu**2, -2*nu**2],
         [0, eta*nu**2, -0.25*t*(4*nu**2+t), -t**2/4, eta*t/4, eta*t/4],
         [0, -2*nu**3/M, 0, 0, 0, 0]]
    return np.array(R)/nu**3          # shape (8,6): AL_i = sum_j R[j,i] A_j

def pols_from_a(a):
    a1,a2,a3,a4,a5,a6 = a
    p = np.array([-a1-a3-a6, a1-a3-a6, (a2-a4+2*a5+a6)/2, (-a2-a4-2*a5+a6)/2,
                  (a2-a4-a6)/2, (-a2-a4-a6)/2]) * 0.5
    p[:2] *= to_fm3; p[2:] *= to_fm4
    return p   # alphaE1, betaM1, gE1E1, gM1M1, gE1M2, gM1E2

def main(fn):
    D = np.loadtxt(fn)
    nu, t = D[:,0], D[:,1]
    names = ['Born','anomaly','piN loops','Delta pole','piDelta loops','contact(LEC)']
    fr = np.round(-t/nu**2, 6)
    flist = [fv for fv in np.unique(fr) if fv > 0]     # t=0 (forward) is singular for the loop integrals
    res = {}
    for c in range(6):
        Are = D[:, 2+16*c: 2+16*c+8]; Aim = D[:, 2+16*c+8: 2+16*c+16]
        AL = np.array([Are[k] @ CAL(nu[k], t[k]) for k in range(len(nu))])
        a = []; spread = []
        for i in range(6):
            # for every line t=-f nu^2: A^L_i = a_i + O(nu^2); fits of increasing order on the
            # low-nu points; the spread over fits and over f is reported as numerical uncertainty
            ests = []
            for fv in flist:
                # LoopTools loses precision for nu < ~0.0125 (E_gamma < ~12 MeV): use nu >= NUMIN
                sel = (fr == fv) & (nu >= NUMIN) & (nu <= NUMAX)
                x2 = nu[sel]**2; y = AL[sel, i]
                for deg in (4, 5, 6):
                    cf = np.polyfit(x2, y, deg)
                    ests.append(cf[-1])
            ests = np.array(ests)
            a.append(np.median(ests)); spread.append(np.max(ests)-np.min(ests))
        p = pols_from_a(np.array(a))
        dp = np.abs(pols_from_a(np.array(spread)))
        res[names[c]] = (p, dp, np.max(np.abs(Aim)))
    hdr = '%-15s %9s %9s %9s %9s %9s %9s %9s %9s' % ('component','alphaE1','betaM1','gE1E1','gM1M1','gE1M2','gM1E2','gamma0','gammapi')
    print(hdr)
    tot = np.zeros(6)
    for nm in names[2:5]+['contact(LEC)']:
        p, dp, im = res[nm]
        g0 = -(p[2]+p[3]+p[4]+p[5]); gp = p[3]+p[5]-p[2]-p[4]
        print('%-15s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f   (max fit spread %.1e)' %
              (nm, *p, g0, gp, dp.max()))
        tot += p
    g0 = -(tot[2]+tot[3]+tot[4]+tot[5]); gp = tot[3]+tot[5]-tot[2]-tot[4]
    print('%-15s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f' % ('TOTAL', *tot, g0, gp))
    try:
        from analytic_pols import pols
        consts = open(fn).readline() and [l for l in open(fn) if l.startswith('# constants')]
        ana = pols()
        print('\nanalytic (LMP2015 App. C, default constants of in_driver.template):')
        for nm in ['piN loops', 'Delta pole', 'piDelta loops', 'TOTAL']:
            v = ana[nm]; g0 = -(v[2]+v[3]+v[4]+v[5]); gp = v[3]+v[5]-v[2]-v[4]
            print('%-15s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f' % (nm, *v, g0, gp))
        print('NOTE: the 2015 code: alphaE1(piDelta) from the amplitude is biased by ~+0.03 by the Thomson-limit')
        print('      violation of the piDelta loops (one-sided p^2=M^2+1e-4 subtraction); fixed here.')
    except Exception as e:
        print('analytic comparison skipped:', e)
    p,_,_ = res['anomaly']
    print('%-15s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f   (pi0 pole; excluded from polarisabilities by convention)' % ('anomaly', *p))
    return res, tot

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'static_As.out')
