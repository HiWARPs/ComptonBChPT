"""Role of the asymmetry data in the BChPT fits: main (F2+trunc) and F3 (spin pols free), with/without A2 Sigma3,
and with Sigma3 + HIGS polarized only.  Uses s03_fits_L.do_fit (cache from `make`)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'scripts')); sys.path.insert(0, os.path.join(here, '..', 'lib'))
os.chdir(os.path.join(here, '..'))
import s03_fits_L as S
ALL = S.OLD_L + S.A2 + S.HIGS
print('A2 sets:', S.A2, ' HIGS sets:', S.HIGS)
noS3 = [k for k in ALL if 'S3' not in k]
for kind in ('F2', 'F3'):
    for th in (True, False):
        sfx = '_thB' if th else ''
        S.do_fit(f'{kind}_all{sfx}', ALL, kind, thcov=th, decompose=False, store=False)
        S.do_fit(f'{kind}_noA2S3{sfx}', noS3, kind, thcov=th, decompose=False, store=False)
