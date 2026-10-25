import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import photoprod as pp, evalamp as ev
c0 = ev.consts()
zs, wz = np.polynomial.legendre.leggauss(48); zs = 0.5*(zs+1); wz = 0.5*wz
gA, F, kv, ks = c0['gA'], c0['F'], c0['kp']-c0['kn'], c0['kp']+c0['kn']
for mu in (0.1487, 0.05, 0.02, 0.01, 0.005):
    c = dict(c0, mpi2=mu**2); nuth = mu + mu**2/2
    g0lo = g0k = 0.
    for z, w in zip(zs, wz):
        sc = 2*mu; nu = nuth + sc*z/(1-z)**2; jac = sc*(1+z)/(1-z)**3
        r = pp.helicity_xs(np.sqrt(1+2*nu), c, nth=20, nph=8)
        g0lo += -w*jac*(r[('lo','h32')]-r[('lo','h12')])/nu**3/(4*np.pi**2)
        g0k  += -w*jac*(r[('int','h32')]-r[('int','h12')])/nu**3/(4*np.pi**2)
    print(f'mu={mu}: mu^2 g0(LO^2) = {mu**2*g0lo:.5f} [HB {gA**2/(24*np.pi**2*F**2*4*np.pi) * 4*np.pi:.5f}]   mu g0(kappa) = {mu*g0k:.5f} [HB {-(gA**2/(192*np.pi*F**2))*(3*kv+ks):.5f}]')
