import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'py'))
import formtools as ft, kin
KS = '#procedure KSUB\n#endprocedure\n'
S = ['(-ep.ee)*gi_(1)', '(q.ep*qp.ee)*gi_(1)', '(-(g_(1,ep,ee)-ep.ee))', '(ep.ee*(g_(1,qp,q)-qp.q))',
     '(q.ep*(g_(1,qp,ee)-qp.ee)-qp.ee*(g_(1,q,ep)-q.ep))', '(q.ep*(g_(1,q,ee)-q.ee)-qp.ee*(g_(1,qp,ep)-qp.ep))',
     '(q.ep*qp.ee*(g_(1,qp,q)-qp.q))']
# S8 = ag(q,qp,ep,ee) antisymmetrised product of 4 (FORM T8 of s_nucl.frm), via explicit antisymmetrisation
import itertools
def perm_sign(p):
    s=1; p=list(p)
    for i in range(len(p)):
        for j in range(i+1,len(p)):
            if p[i]>p[j]: s=-s
    return s
vs=['q','qp','ep','ee']; terms=[]
for p in itertools.permutations(range(4)):
    terms.append(('+' if perm_sign(p)>0 else '-')+'g_(1,'+','.join(vs[i] for i in p)+')')
S.append('(1/24)*('+''.join(terms)+')')
gram = [ft.compile_proj(ft.run_form(f'gram{k+1}', KS + f'L GAM = {S[k]};\n.sort\n#call Proj\n#call Finish\n')) for k in range(8)]
for (x, t) in [(0.05, -0.008), (0.1, -0.03), (0.2, -0.05)]:
    en = kin.env(x, t, 0.02); en.update(T={(0,0,0,0):1.0}, U={})
    G = np.array([[eval(gram[k][j], en) for k in range(8)] for j in range(6)])
    c7 = np.linalg.solve(G[:, :6], G[:, 6]); c8 = np.linalg.solve(G[:, :6], G[:, 7])
    print('S7 =', np.round(c7, 6), ' identity:', [0, 0, -t*(x**2+t/2), 0, x**2+t/2, -x**2])
    print('S8 =', np.round(c8, 6), ' identity:', [0, 0, -t*x/2, 0, x/2, -x/2])
