import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev
cs = ev.consts()
def lec(x, t, da, db):
    y = -x - t/2
    return np.array([0.5*(da+db)*(x**2+y**2)+0.5*t*db, -db-(da+db)*t/8, 0, 0, 0, 0])
pts = [(0.1, -0.05), (0.2, -0.1), (0.3, -0.3), (0.05, -0.01)]
for order in (3, 4):
    m = assemble.Model(order, cs)
    rows, rhs = [], []
    for x, t in pts:
        A = m.A(x, t); Ap = m.A(x, t, pole=True)
        print(f'order {order} x={x} t={t}: |A_pole|/|A| = {np.max(abs(Ap))/np.max(abs(A)):.2e}   A_pole = {np.round(Ap.real, 6)}')
        rows += list(np.vstack([lec(x, t, 1, 0), lec(x, t, 0, 1)]).T); rhs += list(Ap.real)
    c, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    res = np.array(rhs) - np.array(rows) @ c
    print(f'   fit to the delta-alpha, delta-beta contact form: (da, db) = {c},  max residual {np.max(abs(res)):.2e}')
