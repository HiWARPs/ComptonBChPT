#!/usr/bin/env python3
"""Compare the numerical output files (x*.output) of two main-program runs of ampl.exe.

usage: cmp_outputs.py DIR_A DIR_B [files...]
Prints, per file, the max relative and absolute difference of all numeric columns
(relative to max(|a|,|b|,floor)); exits 0.  Used by run_regression.sh.
"""
import sys, os, glob
import numpy as np


def load(fn):
    rows = []
    for line in open(fn):
        line = line.strip().replace('{', ' ').replace('}', ' ').replace(',', ' ')
        if not line:
            continue
        try:
            rows.append([float(v) for v in line.split()])
        except ValueError:
            continue
    if not rows:
        return np.zeros((0, 0))
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


def main():
    da, db = sys.argv[1], sys.argv[2]
    files = sys.argv[3:] or sorted(os.path.basename(f) for f in glob.glob(os.path.join(da, 'x*.output')))
    worst = 0.0
    for f in files:
        A, B = load(os.path.join(da, f)), load(os.path.join(db, f))
        if A.shape != B.shape:
            print('%-36s SHAPE MISMATCH %s %s' % (f, A.shape, B.shape))
            continue
        if A.size == 0:
            print('%-36s (empty)' % f)
            continue
        floor = 1e-6 * np.max(np.abs(A), axis=0, keepdims=True) + 1e-300
        rel = np.abs(A - B) / np.maximum(np.maximum(np.abs(A), np.abs(B)), floor)
        k = np.unravel_index(np.argmax(rel), rel.shape)
        worst = max(worst, rel.max())
        print('%-36s max rel diff %9.2e  (row %d col %d: %.8g vs %.8g)   max abs diff %9.2e' %
              (f, rel.max(), k[0], k[1], A[k], B[k], np.max(np.abs(A - B))))
    print('WORST max rel diff: %.3e' % worst)


if __name__ == '__main__':
    main()
