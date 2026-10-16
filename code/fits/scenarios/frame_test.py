"""Test the frame assignment of Delta-region cross-section sets whose frame is inferred/assumed:
for each point compute the v2 theory (main-fit alpha,beta) under (A) the adopted frame and (B) the other frame
(angle and cross section both reinterpreted), and (C) adopted angle frame but cross section in the other frame."""
import os, sys, numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..', 'lib')); os.chdir(os.path.join(here, '..'))
import bchpt
d = pd.read_csv('../delta_plots/results/data_with_theory.csv')
sets = ['Wissmann1999', 'Camen2002', 'Kolb2000', 'Hallin1993_TabII_D', 'Hallin1993_TabII', 'Genzel1972', 'Gray1967',
        'Nagashima1964', 'Baranov1966', 'LEGS2001_CS', 'Hunger1997_75', 'Hunger1997_90', 'OlmosdeLeon2001']
d = d[d.dataset.isin(sets) & d.is_cs & (d.E > 150)].reset_index(drop=True)
pts = []
for i, r in d.iterrows():
    adopted_lab = r.frame == 'LAB'
    th = r.theta
    # A: adopted
    pts.append((f'A{i}', r.E, th, 'CS' if adopted_lab else 'CSCM', 'LAB' if adopted_lab else 'CM'))
    # B: both reinterpreted (angle in the other frame, cross section in the other frame)
    pts.append((f'B{i}', r.E, th, 'CSCM' if adopted_lab else 'CS', 'CM' if adopted_lab else 'LAB'))
    # C: angle as adopted, cross section in the other frame
    pts.append((f'C{i}', r.E, th, 'CSCM' if adopted_lab else 'CS', 'LAB' if adopted_lab else 'CM'))
out = bchpt.run(pts, 'NNLO', name='frame_test', verbose=False)
P = np.array(list(out.Pref)); P[0], P[1] = 11.672, 2.876
pred = out.predict_array(P)
val = {pid: pred[out.index[pid]] for pid in out.ids}
rows = []
for i, r in d.iterrows():
    rows.append(dict(dataset=r.dataset, frame=r.frame, E=r.E, theta=r.theta, value=r.value, err=r.err,
                     A=val[f'A{i}'], B=val[f'B{i}'], C=val[f'C{i}']))
t = pd.DataFrame(rows)
print(f"{'dataset':22s} {'adopted':7s} {'N':>3s} | median data/theory:  A(adopted)  B(other frame)  C(xs only other) | chi2/N A  B  C")
for ds, g in t.groupby('dataset', sort=False):
    s = lambda k: (np.median(g.value / g[k]), np.mean(((g.value - g[k]) / g.err) ** 2))
    a, b, c = s('A'), s('B'), s('C')
    print(f"{ds:22s} {g.frame.iloc[0]:7s} {len(g):3d} | {a[0]:10.3f} {b[0]:12.3f} {c[0]:14.3f}      | {a[1]:6.1f} {b[1]:6.1f} {c[1]:6.1f}")
t.to_csv('scenarios/frame_test.csv', index=False)
