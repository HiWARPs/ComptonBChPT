"""HB limit of the kappa-insertion spin polarisabilities: mu = m_pi/M -> 0 at fixed F (M units)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev, statpol
base = ev.consts()
to4 = statpol.pols_from_a.__globals__['to_fm4']
res = []
for mu in (0.12, 0.08, 0.05, 0.03, 0.02):
    cs = dict(base, mpi2=mu**2)
    m = assemble.Model(4, cs)
    nus = np.linspace(0.12, 0.8, 9) * mu
    p, dp = statpol.extract(lambda x, t: m.A(x, t), nus=nus)
    g = np.array(list(p[2:]) + list(statpol.g0pi(p))) / to4      # internal units alpha/M^4
    res.append((mu, g, p[:2] / statpol.pols_from_a.__globals__['to_fm3']))
    print(f'mu={mu}: mu*[gEE gMM gEM gME g0 gpi] = {np.round(mu*g, 5)}   [alpha, beta] = {np.round(res[-1][2], 4)}  spread {dp.max():.1e}')
gA, F, kv, ks = base['gA'], base['F'], base['kp'] - base['kn'], base['kp'] + base['kn']
print('HB prediction mu*gamma0 -> -(gA^2/(192 pi F^2))(3 kv + ks) =', -(gA**2 / (192 * np.pi * F**2)) * (3 * kv + ks))
