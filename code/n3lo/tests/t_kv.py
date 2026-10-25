"""kappa_v and kappa_s parts of the O(p^4) kappa-insertion alpha, beta (and gammas) vs m_pi/M."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev, statpol, photoprod as pp
base = ev.consts()
zs, wz = np.polynomial.legendre.leggauss(40); zs = 0.5*(zs+1); wz = 0.5*wz
conv3 = (0.19733/0.93827)**3*1e4
for mu in (0.1487, 0.10, 0.07):
    out = {}
    for lab, kp, kn in (('kv', 0.5, -0.5), ('ks', 0.5, 0.5)):     # unit kappa_v or kappa_s
        cs = dict(base, mpi2=mu**2, kp=kp, kn=kn)
        m = assemble.Model(4, cs)
        nus = np.linspace(0.13, 0.8, 9)*mu
        p, dp = statpol.extract(lambda x, t: m.A(x, t), nus=nus)
        # sum rule for alpha+beta of the same piece
        nuth = mu + mu**2/2; ab = 0.
        for z, w in zip(zs, wz):
            sc = 2*mu; nu = nuth + sc*z/(1-z)**2; jac = sc*(1+z)/(1-z)**3
            r = pp.helicity_xs(np.sqrt(1+2*nu), cs, nth=20, nph=8)
            ab += w*jac*0.5*(r[('int','h12')]+r[('int','h32')])/nu**2/(2*np.pi**2)
        out[lab] = (p, ab*conv3, dp.max())
    for lab in ('kv', 'ks'):
        p, ab, sp_ = out[lab]
        print(f'mu={mu}: per unit {lab}: alpha={p[0]:8.4f} beta={p[1]:8.4f}  a+b={p[0]+p[1]:8.4f} (SR {ab:8.4f})  a-b={p[0]-p[1]:8.4f}  g0={statpol.g0pi(p)[0]:8.4f} gpi={statpol.g0pi(p)[1]:8.4f}  spread {sp_:.1e}')
