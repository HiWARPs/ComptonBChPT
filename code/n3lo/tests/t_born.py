import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'py'))
import formtools as ft, kin
KS = '#procedure KSUB\n#endprocedure\nSymbols Q,kap;\n'
born = ft.run_form('born', KS + """
L GAM = (-ixs)*(Q*g_(1,ep)+kap/2/M*(g_(1,ep,qp)-ep.qp))*(g_(1,p)+g_(1,q)+M*gi_(1))*(Q*g_(1,ee)-kap/2/M*(g_(1,ee,q)-ee.q))
       +(-iys)*(Q*g_(1,ee)-kap/2/M*(g_(1,ee,q)-ee.q))*(g_(1,p)-g_(1,qp)+M*gi_(1))*(Q*g_(1,ep)+kap/2/M*(g_(1,ep,qp)-ep.qp));
.sort
#call Proj
#call Finish
""")
S = ['(-ep.ee)*gi_(1)', '(q.ep*qp.ee)*gi_(1)', '(-(g_(1,ep,ee)-ep.ee))', '(ep.ee*(g_(1,qp,q)-qp.q))',
     '(q.ep*(g_(1,qp,ee)-qp.ee)-qp.ee*(g_(1,q,ep)-q.ep))', '(q.ep*(g_(1,q,ee)-q.ee)-qp.ee*(g_(1,qp,ep)-qp.ep))']
gram = [ft.compile_proj(ft.run_form(f'gram{k+1}', KS + f'L GAM = {S[k]};\n.sort\n#call Proj\n#call Finish\n')) for k in range(6)]
cb = ft.compile_proj(born)
def snucl(x, t, Q, kap):
    Qk2 = (Q + kap) ** 2; k2 = kap ** 2
    a1 = x*t*k2/8 + x*Q**2 + x**2*k2/2 + t*Qk2/4
    a2 = x*k2/4 + t*k2/16 + Q**2/2 - Qk2/2
    a3 = a1; a4 = x*k2/4 + Qk2/2; a5 = -Qk2/2; a6 = Q**2/4 + Qk2/4 - k2/4; a7 = -k2/8; a8 = -a4
    return -np.array([a1, a2, a3, a4, a5, a6, a7, a8]) / (2 * x)
def born_code(x, t, Q, kap):
    y = -x - t / 2
    As = snucl(x, t, Q, kap); Au = snucl(y, t, Q, kap)
    A = As.copy(); A[[0, 1, 7]] += Au[[0, 1, 7]]; A[2:7] -= Au[2:7]
    # 8 -> 6 (identity in computeAs)
    A[2] += -t*(x**2+t/2)*A[6] - t*(x/2)*A[7]; A[4] += (x**2+t/2)*A[6] + (x/2)*A[7]; A[5] += -x**2*A[6] - (x/2)*A[7]
    return A[:6]
T = {(0, 0, 0, 0): 1.0}
for (x, t) in [(0.05, -0.03), (0.1, -0.15), (0.2, -0.05)]:
    for Q, kap in [(1., 0.), (1., 1.793), (0., -1.913)]:
        en = kin.env(x, t, 0.02); en.update(T=T, U={}, Q=Q, kap=kap)
        G = np.array([[eval(gram[k][j], en) for k in range(6)] for j in range(6)])
        V = np.array([eval(cb[j], en) for j in range(6)])
        A = np.linalg.solve(G, V)
        Ac = born_code(x, t, Q, kap)
        print(f'x={x} t={t} Q={Q} k={kap}: max|A-Acode| = {np.max(np.abs(A-Ac)):.2e}   A={np.round(A,5)}')

print('--- physical comparison: projections of the code Born (8 structures) vs FORM Born')
exec(open(os.path.join(os.path.dirname(__file__), 't_basis.py')).read().split("gram = ")[0].split("KS = ")[1].join(["KS = ",""]) if False else '')
import importlib
sys.argv=['x']
S8 = None
