import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import statpol, n3lo
m = n3lo.P4()
cache = {}
def parts(x, t):
    if (x, t) not in cache: cache[(x, t)] = m.parts(x, t)
    return cache[(x, t)]
print(f'{"":10s} {"alphaE1":>9s} {"betaM1":>9s} {"gE1E1":>9s} {"gM1M1":>9s} {"gE1M2":>9s} {"gM1E2":>9s} {"gamma0":>9s} {"gammapi":>9s}')
for key in ('kappa', 'c1', 'c2', 'c3', 'total'):
    p, dp = statpol.extract(lambda x, t: parts(x, t)[key])
    print(f'{key:10s}', ' '.join(f'{v:9.4f}' for v in list(p) + list(statpol.g0pi(p))), f'  (spread {dp.max():.1e})')
