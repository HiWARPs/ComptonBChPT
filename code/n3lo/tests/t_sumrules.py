import sys, os, numpy as np, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import photoprod as pp, evalamp as ev
from scipy import integrate
c = ev.consts()
m = np.sqrt(c['mpi2']); nuth = m + m**2/2
t0 = time.time(); r = pp.helicity_xs(1.3, c, nth=24, nph=8); print('one W', time.time()-t0, r)
# integrate over nu with variable change nu = nuth/(1-z)^... use z in (0,1): nu = nuth + w*z/(1-z)
def xs(nu):
    W = np.sqrt(1 + 2*nu); return pp.helicity_xs(W, c, nth=24, nph=8)
zs, wz = np.polynomial.legendre.leggauss(40); zs = 0.5*(zs+1); wz = 0.5*wz
I = {k: 0. for k in ('ab_lo','g0_lo','ab_int','g0_int')}
for z, w in zip(zs, wz):
    sc = 0.3
    nu = nuth + sc*z/(1-z)**2; jac = sc*(1+z)/(1-z)**3
    r = xs(nu)
    for part in ('lo', 'int'):
        stot = 0.5*(r[(part,'h12')] + r[(part,'h32')]); dsig = r[(part,'h32')] - r[(part,'h12')]
        I['ab_'+part] += w*jac*stot/nu**2/(2*np.pi**2)
        I['g0_'+part] += -w*jac*dsig/nu**3/(4*np.pi**2)
conv3 = (0.19733/0.93827)**3*1e4; conv4 = conv3*(0.19733/0.93827)
print('alpha+beta  LO^2 :', I['ab_lo']*conv3, '  (p^3 loops: 5.0975)')
print('gamma0      LO^2 :', I['g0_lo']*conv4, '  (p^3 loops: 2.0135)')
print('alpha+beta  kappa:', I['ab_int']*conv3, '  (p^4 kappa loops: 0.6543)')
print('gamma0      kappa:', I['g0_int']*conv4, '  (p^4 kappa loops: -0.7823)')
