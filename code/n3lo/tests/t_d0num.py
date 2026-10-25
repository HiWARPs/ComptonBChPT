import sys, os, numpy as np
from scipy import integrate
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lt'))
import ltpy
mp = (139.57/938.27)**2
def D0num(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4):
    # k0=0,k1=p1,k2=p1+p2,k3=p1+p2+p3 ; (ki-kj)^2
    K = {(0,1):p1,(1,2):p2,(2,3):p3,(0,3):p4,(0,2):s12,(1,3):s23}
    m = [m1,m2,m3,m4]
    def f(z, y, x):
        xs = [x, y, z, 1-x-y-z]
        if xs[3] < 0: return 0.
        D = sum(xs[i]*m[i] for i in range(4)) - sum(xs[i]*xs[j]*K[(i,j)] for (i,j) in K)
        return 1.0/D**2
    v, e = integrate.nquad(f, [lambda y, x: [0, 1-x-y], lambda x: [0, 1-x], [0, 1]], opts={'epsabs':1e-10,'epsrel':1e-9})
    return v
for nu, f in [(0.051, 1.0), (0.092, 0.5), (0.133, 2.0)]:
    t = -f*nu**2; x = nu - t/4; s = 1+2*x; u = 2-s-t
    args = (0., 1., 0., 1., u, s, mp, mp, 1., 1.)
    ltpy.setversion(0); ltpy.clear(); a, _ = ltpy.D(*args)
    ltpy.setversion(64); ltpy.clear(); b, _ = ltpy.D(*args)
    print(f'nu={nu} f={f}: FF {a["dd0"].real:.8g}  Denner {b["dd0"].real:.8g}  numeric {D0num(*args):.8g}')
