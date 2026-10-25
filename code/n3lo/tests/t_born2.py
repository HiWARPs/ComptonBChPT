import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'py'))
import formtools as ft, kin
KS = '#procedure KSUB\n#endprocedure\nSymbols Q,kap;\n'
born = ft.compile_proj(ft.run_form('born', KS + """
L GAM = (-ixs)*(Q*g_(1,ep)+kap/2/M*(g_(1,ep,qp)-ep.qp))*(g_(1,p)+g_(1,q)+M*gi_(1))*(Q*g_(1,ee)-kap/2/M*(g_(1,ee,q)-ee.q))
       +(-iys)*(Q*g_(1,ee)-kap/2/M*(g_(1,ee,q)-ee.q))*(g_(1,p)-g_(1,qp)+M*gi_(1))*(Q*g_(1,ep)+kap/2/M*(g_(1,ep,qp)-ep.qp));
.sort
#call Proj
#call Finish
"""))
gram = [ft.compile_proj(ft.run_form(f'gram{k+1}', '')) if False else None for k in range(8)]
import glob
gram=[]
for k in range(8):
    txt=open(f'out/gram{k+1}.log').read(); gram.append(ft.compile_proj(ft.parse_fortran(txt)))
def snucl(x, t, Q, kap):
    Qk2 = (Q + kap) ** 2; k2 = kap ** 2
    a1 = x*t*k2/8 + x*Q**2 + x**2*k2/2 + t*Qk2/4
    a2 = x*k2/4 + t*k2/16 + Q**2/2 - Qk2/2
    a3 = a1; a4 = x*k2/4 + Qk2/2; a5 = -Qk2/2; a6 = Q**2/4 + Qk2/4 - k2/4; a7 = -k2/8; a8 = -a4
    return -np.array([a1, a2, a3, a4, a5, a6, a7, a8]) / (2 * x)
for (x, t) in [(0.05, -0.008), (0.1, -0.03), (0.2, -0.05)]:
  for Q, kap in [(1., 0.), (1., 1.793), (0., -1.913)]:
    y = -x - t/2
    As = snucl(x, t, Q, kap); Au = snucl(y, t, Q, kap)
    A = As.copy(); A[[0, 1, 7]] += Au[[0, 1, 7]]; A[2:7] -= Au[2:7]
    en = kin.env(x, t, 0.02); en.update(T={(0,0,0,0):1.0}, U={}, Q=Q, kap=kap)
    G = np.array([[eval(gram[k][j], en) for k in range(8)] for j in range(6)])
    Vc = G @ A
    V = np.array([eval(born[j], en) for j in range(6)])
    print(f'x={x} t={t} Q={Q} k={kap}: rel dev = {np.max(np.abs(V-Vc))/np.max(np.abs(V)):.2e}')
