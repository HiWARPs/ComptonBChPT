import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import photoprod as pp, evalamp as ev
base = ev.consts()
zs, wz = np.polynomial.legendre.leggauss(56); zs = 0.5*(zs+1); wz = 0.5*wz
conv3 = (0.19733/0.93827)**3*1e4
C_over = (1/137.036)/(96*np.pi**2*base['F']**2)*base['gA']**2/np.pi * conv3      # C g^2/(pi m), M=1
print('unit of HB O(q^4) terms (C gA^2/(pi M)) =', C_over)
for mu in (0.1487, 0.07, 0.03, 0.01, 0.003, 0.001):
    r_ = {}
    for lab, kp, kn in (('kv', 0.5, -0.5), ('ks', 0.5, 0.5)):
        cs = dict(base, mpi2=mu**2, kp=kp, kn=kn); nuth = mu + mu**2/2; ab = 0.
        for z, w in zip(zs, wz):
            sc = 2*mu; nu = nuth + sc*z/(1-z)**2; jac = sc*(1+z)/(1-z)**3
            r = pp.helicity_xs(np.sqrt(1+2*nu), cs, nth=20, nph=8)
            ab += w*jac*0.5*(r[('int','h12')]+r[('int','h32')])/nu**2/(2*np.pi**2)
        r_[lab] = ab*conv3
    print(f'mu={mu}: (a+b) per unit kv = {r_["kv"]:.5f}   per unit ks = {r_["ks"]:.5f}   [HB ks: 3 C/pi g^2/M (ln mu + 1/2) = {C_over*3*(np.log(mu)+0.5):.5f} + const]')
