#!/usr/bin/env python3
"""Static dipole polarisabilities of the proton at NNLO (O(p^3) + O(p^4/Delta)) from the analytic
expressions of Lensky, McGovern, Pascalutsa, EPJC 75 (2015) 604, Appendix C (draft LP_ms3 tex),
evaluated with a given parameter set.  Units: alpha,beta [1e-4 fm^3], gammas [1e-4 fm^4].
Ragusa -> multipole: gE1E1=-g1-g3, gM1M1=g4, gE1M2=g3, gM1E2=g2+g4."""
import numpy as np
from scipy.integrate import quad

def pols(gA=1.27, F=0.09221, MN=0.93827, mpi=0.13957, MD=1.232, hA=2.85, gM=2.97, gE=-1.0,
         aem=1/137.036, hbarc=0.19733, Xi=0.0):
    e2 = 4*np.pi*aem
    mu = mpi/MN
    s = np.sqrt(4-mu**2); ac = np.arccos(mu**2/2-1); L = np.log(mu**2)
    cN = e2*gA**2/(192*np.pi**3*F**2*MN)
    aN = cN/(mu*s**5)*(mu*(-18*mu**6+157*mu**4-407*mu**2+304)*s + mu*(9*mu**4-20*mu**2+9)*s**5*L
                       -(9*mu**10-110*mu**8+479*mu**6-870*mu**4+590*mu**2-80)*ac)
    bN = -cN/(mu*s**3)*(mu*(-54*mu**4+235*mu**2-127)*s - mu*(27*mu**4-50*mu**2+9)*s**3*L
                        -(27*mu**8-212*mu**6+471*mu**4-246*mu**2+2)*ac)
    cg = e2*gA**2/(192*np.pi**3*F**2*MN**2)
    g1N = cg/(mu**2*s**5)*(s*(18*mu**8-173*mu**6+524*mu**4-452*mu**2-(4-mu**2)**2*(9*mu**4-28*mu**2+13)*mu**2*L+32)
                          + mu*(9*mu**10-118*mu**8+563*mu**6-1150*mu**4+848*mu**2-104)*ac)
    g2N = 2*cg/(mu**2*s**5)*(2*s*(18*mu**8-149*mu**6+348*mu**4-179*mu**2-(4-mu**2)**2*(9*mu**4-16*mu**2+4)*mu**2*L+4)
                          + mu*(18*mu**10-212*mu**8+868*mu**6-1400*mu**4+701*mu**2-44)*ac)
    g3N = -cg/2/(mu**2*s**5)*(s*(18*mu**8-141*mu**6+284*mu**4-4*mu**2-3*(mu**2-4)**2*(3*mu**4-4*mu**2-1)*mu**2*L-16)
                          + mu*(9*mu**10-102*mu**8+387*mu**6-510*mu**4+46*mu**2+56)*ac)
    g4N = cg/2/(mu**2*s**5)*(s*(-126*mu**8+1047*mu**6-2462*mu**4+1308*mu**2+(mu**2-4)**2*(63*mu**4-114*mu**2+29)*mu**2*L-16)
                          + (-63*mu**11+744*mu**9-3059*mu**7+4970*mu**5-2534*mu**3+152*mu)*ac)
    # piDelta loops
    d = (MD-MN)/MN
    D = lambda x: x**2+(1+d)**2-x*(2+2*d+d**2-mu**2)
    cD = e2*hA**2/(1728*np.pi**3*F**2*MD**2)
    I = lambda f: quad(f, 0, 1, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
    aD = cD*MN*(8*d+35/3+I(lambda x: 18*(x-1)**2*x**4*(x+d+1)/D(x)**2
                 - x**2*(120*x**3+4*(27*d-16)*x**2-(166*d+115)*x+57*(d+1))/D(x)
                 + 6*x*(x*(22*x+15*d-1)-12*(d+1))*(Xi-np.log(D(x)))))
    bD = cD*MN*(-8*d-65/6+I(lambda x: (24*x**2-32*x+9)*x**2*(x+d+1)/D(x)
                 - 6*(5*x-4)*x*(2*x+3*d+3)*(Xi-np.log(D(x)))))
    g1D = cD*(2*d+7/4+I(lambda x: 8*(x-1)**3*x**5*(x+d+1)/D(x)**3
                 - 2*(x-1)**2*x**3*(22*x**2+8*d*x+3*d+3)/D(x)**2
                 + (x-1)*x**2*(105*x**2-2*(5*d+38)*x+33*(d+1))/D(x)
                 + 2*x*(x*(-16*x+9*d+21)-12*(d+1))*(Xi-np.log(D(x)))))
    g2D = cD*(5/12+I(lambda x: 2*(4*x**2-7*x+3)*x**3*(x+d+1)/D(x)**2
                 - x**2*(21*x**3-4*(3*d+5)*x**2+(37*d+26)*x-24*(d+1))/D(x)
                 + 12*(3*x-2)*x**2*(Xi-np.log(D(x)))))
    g3D = cD*(-d-2/3+I(lambda x: (x-1)*x**5*(x+d+1)/D(x)**2
                 - x**2*(24*x**3+(2*d-11)*x**2+(8*d-1)*x-9*(d+1))/(2*D(x))
                 + x*(x*(34*x-9*d-33)+12*(d+1))*(Xi-np.log(D(x)))))
    g4D = cD*(d+2/3+I(lambda x: -(x**3+5*x**2-12*x+6)*x**3*(x+d+1)/D(x)**2
                 + x**2*(48*x**3-(22*d+59)*x**2+(80*d+71)*x-57*(d+1))/(2*D(x))
                 + x*(x*(-34*x+9*d+33)-12*(d+1))*(Xi-np.log(D(x)))))
    # Delta pole
    Mp = MN+MD; Dl = MD-MN
    aP = -aem*2*gE**2/Mp**3
    bP = aem*2*gM**2/(Mp**2*Dl)
    k = aem/(Mp**2*MN)
    g1P = k*(gM**2/Dl - 2*gE*gM/Mp + gE**2*MD/Mp**2)
    g2P = -k*(gM**2*MD/Dl**2 + gE*gM/Dl)
    g3P = -k*(gM**2/(2*Dl) - 3*gE*gM/(2*Mp))
    g4P = k*(gM**2*MD/Dl**2 - gE*gM/(2*Dl) + gE**2/(2*Mp))
    f3 = hbarc**3*1e4; f4 = hbarc**4*1e4     # GeV^-3 -> 1e-4 fm^3 etc.
    out = {}
    for nm,(a,b,g1,g2,g3,g4) in {'piN loops':(aN,bN,g1N,g2N,g3N,g4N),
                                  'Delta pole':(aP,bP,g1P,g2P,g3P,g4P),
                                  'piDelta loops':(aD,bD,g1D,g2D,g3D,g4D)}.items():
        gEE, gMM, gEM, gME = -g1-g3, g4, g3, g2+g4
        v = np.array([a*f3, b*f3, gEE*f4, gMM*f4, gEM*f4, gME*f4])
        out[nm] = v
    out['TOTAL'] = sum(out.values())
    return out

def table(out, title=''):
    print(title)
    print('%-15s %9s %9s %9s %9s %9s %9s %9s %9s'%('component','alphaE1','betaM1','gE1E1','gM1M1','gE1M2','gM1E2','gamma0','gammapi'))
    for nm,v in out.items():
        g0 = -(v[2]+v[3]+v[4]+v[5]); gp = v[3]+v[5]-v[2]-v[4]
        print('%-15s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f'%(nm,*v,g0,gp))

if __name__ == '__main__':
    table(pols(), 'Current in.dat constants (gA=1.27, fpi=92.21, mpi=139.57, M=938.27, MD=1232, hA=2.85, alpha=1/137.036):')
    table(pols(gA=1.267, F=0.0924, mpi=0.139, MD=1.31*0.93827),
          '\n2012 constants (gA=1.267, fpi=92.4, mpi=139.0, MD=1.31 M as hard-coded in old CorrIso loop):')
    table(pols(gA=1.267, F=0.0924, mpi=0.139, MD=1.31*0.93827, aem=1/137.04),
          '\n2012 constants, alpha_em=1/137.04:')
