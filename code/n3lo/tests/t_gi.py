import sys, os, numpy as np, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev
cs = ev.consts()
for order in (3, 4):
    t0 = time.time()
    m = assemble.Model(order, cs); m.gi_compile()
    print(f'order {order}: built in {time.time()-t0:.0f} s')
    for x, t in [(0.1, -0.05), (0.2, -0.1), (0.3, -0.3)]:
        W, big = m.GI(x, t)
        print(f'   x={x} t={t}: max|W_total| / max|W_graph| = {np.max(np.abs(W))/big:.2e}')
