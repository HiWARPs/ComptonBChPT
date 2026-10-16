"""chi^2 of the modern experiments (A2 2022 dsigma + Sigma3, HIGS 2022, TAPS 2001, Saskatoon = Hallin 1993 +
MacGibbon 1995) in the PWA fit of the paper vs the BChPT Baldin fit without truncation covariance (fits_L.json
'F1_d_all'), both on the MAIN DATABASE and with the same error model and floating normalisations (data chi^2 only,
normalisation penalties excluded).  This is the comparison of Sect. IV F of the paper ("Which of the two frameworks
describes the data better?").  Writes results/modern_chi2_comparison.{json,md}.  Run from work/pwa/."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PWA = os.path.dirname(HERE)
os.chdir(PWA)
FITS_L = os.path.join(os.path.dirname(os.path.dirname(PWA)), 'work', 'fits', 'results', 'fits_L.json')
BCHPT_KEY = sys.argv[1] if len(sys.argv) > 1 else 'F1_d_all'

# (label, PWA set names, BChPT set names)
GROUPS = [('A2 2022 dsigma', ['A2_2022_CS'], ['A2_2022_CS']),
          ('A2 2022 Sigma3', ['A2_2022_S3'], ['A2_2022_S3']),
          ('HIGS 2022', ['HIGS_2022_circ', 'HIGS_2022_lin'], ['HIGS_2022_CIRC', 'HIGS_2022_LIN']),
          ('TAPS 2001', ['TAPS2001'], ['OlmosdeLeon2001']),
          ('Saskatoon (Hallin 93 + MacGibbon 95)', ['Hallin1993', 'MacGibbon1995_tagged', 'MacGibbon1995_untagged'],
           ['Hallin1993_TabII', 'MacGibbon1995_tagged', 'MacGibbon1995_untagged'])]

pwa = json.load(open('results/main_fits_legacynorm.json'))['MAIN']['per_dataset']
bch = json.load(open(FITS_L))[BCHPT_KEY]['perset']
rows, out = [], {}
NP = NB = CP = CB = 0
for lab, sp, sb in GROUPS:
    np_, cp = sum(pwa[s]['N'] for s in sp), sum(pwa[s]['chi2'] for s in sp)
    nb, cb = sum(bch[s]['N'] for s in sb), sum(bch[s]['chi2'] for s in sb)
    out[lab] = dict(N_pwa=np_, chi2_pwa=cp, N_bchpt=nb, chi2_bchpt=cb)
    rows.append(f'| {lab} | {np_} | {cp:.1f} | {nb} | {cb:.1f} |')
    NP += np_; NB += nb; CP += cp; CB += cb
out['total'] = dict(N_pwa=NP, chi2_pwa=CP, N_bchpt=NB, chi2_bchpt=CB, bchpt_key=BCHPT_KEY)
md = [f'### chi^2 of the modern experiments: PWA fit (main database) vs BChPT `{BCHPT_KEY}` (Baldin sum rule, no truncation covariance; main database)\n',
      'Data chi^2 only (normalisation penalties excluded); both fits float the normalisations and use the same error model.\n',
      '| data set | N (PWA) | chi2 (PWA) | N (BChPT) | chi2 (BChPT) |', '|---|---|---|---|---|'] + rows + \
     [f'| **total** | **{NP}** | **{CP:.1f}** | **{NB}** | **{CB:.1f}** |']
assert NP == NB, (NP, NB)
open('results/modern_chi2_comparison.md', 'w').write('\n'.join(md) + '\n')
json.dump(out, open('results/modern_chi2_comparison.json', 'w'), indent=1)
print('\n'.join(md))
