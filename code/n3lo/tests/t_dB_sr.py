"""t_dB_sr.py -- forward sum rules for the gamma N Delta loop graphs (Bloop): their s-channel pi N cut is the
interference of the LO Born photoproduction amplitude with u-channel Delta exchange (photoprod_delta.py, explicit
eps tensor).  alpha+beta = 1/(2pi^2) int sigma/nu^2,  gamma0 = -1/(4pi^2) int (sigma_3/2 - sigma_1/2)/nu^3, pi N cut
only (the pi Delta cut opens at W = M + MD + ...; listed separately up to nu_max to show convergence)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import photoprod_delta as pd, devalamp as dv
c = dv.dconsts()
m = np.sqrt(c['mpi2']); nuth = m + m ** 2 / 2
conv3 = (0.19733 / 0.93827) ** 3 * 1e4; conv4 = conv3 * (0.19733 / 0.93827)
zs, wz = np.polynomial.legendre.leggauss(32); zs = 0.5 * (zs + 1); wz = 0.5 * wz
for numax in (1.0, 3.0, 10.0, 30.0):
    ab = g0 = 0.0
    # nu = nuth + (numax-nuth) * z^2  (dense near threshold)
    for z, w in zip(zs, wz):
        nu = nuth + (numax - nuth) * z ** 2; jac = 2 * (numax - nuth) * z
        r = pd.helicity_xs_int(np.sqrt(1 + 2 * nu), c, ('u',), nth=20, nph=8)
        stot = 0.5 * (r['h12'] + r['h32']); dsig = r['h32'] - r['h12']
        ab += w * jac * stot / nu ** 2 / (2 * np.pi ** 2)
        g0 += -w * jac * dsig / nu ** 3 / (4 * np.pi ** 2)
    print(f'nu_max = {numax:5.1f} M: alpha+beta = {ab * conv3:8.4f}   gamma0 = {g0 * conv4:8.4f}', flush=True)
