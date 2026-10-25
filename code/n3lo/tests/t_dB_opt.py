"""t_dB_opt.py -- optical theorem for the gamma N Delta loop graphs (Bloop) between the pi N and pi Delta thresholds:
Im f(nu) ~ nu sigma_tot,  Im g(nu) ~ nu (sigma_1/2 - sigma_3/2), with sigma the interference of the LO Born and the
u-channel Delta photoproduction amplitudes (photoprod_delta.py, independent of FORM).  The normalisation of the
forward L'vov combinations is calibrated with the O(p^4) kappa loops (validated by sum rules) and their kappa
interference (photoprod.py)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import statpol, ddelta, devalamp as dv, photoprod as pp, photoprod_delta as pd, n3lo

def fwd(Afun, nu):
    """forward L'vov combinations -(A3+A6)/nu and A4 at t -> 0 (linear extrapolation from t = -1e-4, -2e-4)"""
    out = []
    for t in (-1e-4, -2e-4):
        AL = statpol.lvov(Afun(nu - t / 4, t), nu, t)
        out.append(np.array([-(AL[2] + AL[5]) / nu, AL[3]]))
    return 2 * out[0] - out[1]

cs = dv.dconsts()
P4 = n3lo.P4()
D = ddelta.Delta45(cs, which=('B',))
for nu in (0.2, 0.3, 0.4, 0.5):
    W = np.sqrt(1 + 2 * nu)
    rk = pp.helicity_xs(W, cs, nth=24, nph=8)
    sk = np.array([0.5 * (rk[('int', 'h12')] + rk[('int', 'h32')]), rk[('int', 'h12')] - rk[('int', 'h32')]])
    rd = pd.helicity_xs_int(W, cs, ('u',), nth=24, nph=8)
    sd = np.array([0.5 * (rd['h12'] + rd['h32']), rd['h12'] - rd['h32']])
    Fk = fwd(lambda x, t: P4.parts(x, t)['kappa'], nu)
    Fd = fwd(lambda x, t: D.m['Bloop'].A(x, t), nu)
    ck = Fk.imag / (nu * sk); cd = Fd.imag / (nu * sd)
    print(f'nu={nu}: Im/(nu sigma)  kappa: f {ck[0]:.5f} g {ck[1]:.5f}   gamma N Delta: f {cd[0]:.5f} g {cd[1]:.5f}')
