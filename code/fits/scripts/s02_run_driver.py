"""Step 2: compute the quadratic-form coefficients for all database points and all configurations
(order-by-order, parameter variations, gM grid).  Executable from $AMPL_EXE (default work/code/execs/ampl.exe)."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import bchpt
from common import load_db, outputs

GM_GRID = [2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 2.97, 3.00, 3.05, 3.10, 3.15, 3.20, 3.30, 3.40, 3.50]
for g in GM_GRID:
    bchpt.CONFIGS['gM%.2f' % g] = dict(gM=g)
    bchpt.CONFIGS['gM%.2fnv' % g] = dict(gM=g, dvertex=False)        # v2 without vertex corrections
V2CFG = ['PDG24', 'hA2.84', 'hA2.91', 'singlemass', 'novtx'] if bchpt.PROFILE == 'v2' else ['gA1.2754', 'fpi92.07', 'PDGconst']

if __name__ == '__main__':
    t0 = time.time()
    print('executable:', bchpt.exe_path())
    df = load_db()
    grid = ['gM%.2f' % g for g in GM_GRID] + (['gM%.2fnv' % g for g in GM_GRID] if bchpt.PROFILE == 'v2' else [])
    for cfg in ['NNLO', 'LO', 'NLO', 'BornPi0', 'NNLO_noPiD', 'NNLO_noDpole', 'gM2.87', 'gM3.07'] + V2CFG + grid:
        T = outputs(df, cfg)
        print('%-14s %d points, %d extrapolated' % (cfg, len(T.ids), T.flag.sum()))
    dfp = load_db(per_beamtime=True)
    for cfg in ['NNLO', 'LO', 'NLO', 'BornPi0']:
        outputs(dfp, cfg, tag='dbpb')
    print('done in %.1f s' % (time.time() - t0))
