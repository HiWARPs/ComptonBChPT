import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import n3lo
t0 = time.time(); m = n3lo.P4(); print('setup', round(time.time()-t0,1), 's')
t0 = time.time()
for i in range(5): r = m.parts(0.1 + 0.02*i, -0.05)
print('per point', round((time.time()-t0)/5, 3), 's')
for k, v in r.items(): print(k, np.round(v.real, 6), np.round(v.imag, 6))
