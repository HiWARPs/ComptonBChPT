"""Convert the 2018 PWA database (Krupina, Lensky, Pascalutsa, PLB 782 (2018) 34),
exported from the author's notebook
  <home>/MySci/Math/ComptonScattering/PWA/Final_NB/data_gather.nb  (list `datafull`)
to the standard db CSV format.  138 points = PLB Table I.  Notes:
 - lab angles, lab dsigma/dOmega (nb/sr), single (total) error as used in the PLB fits;
 - E_lab is the value stored in `datafull`; for some sets (e.g. Baranov, Hallin) this is
   the PLB energy-bin label rather than the measured energy (e.g. Baranov 80.9 -> 79);
 - Hallin and Bernardini were converted c.m. -> lab in the notebook.
Usage: python3 -I build_db_klp2018.py <datafull.csv> <outdir>
"""
import csv, sys, os, math
Mp = 938.272
def th_cm(E, thl):
    E, c, M = E/1000.0, math.cos(math.radians(thl)), Mp/1000.0
    Ep = E/(1 + E/M*(1 - c)); t = -2*E*Ep*(1 - c); s = M*M + 2*M*E
    w = (s - M*M)/(2*math.sqrt(s)); return math.degrees(math.acos(max(-1, min(1, 1 + t/(2*w*w)))))
tex = {'OlmosDeLeon':'deLeon:2001dnx','Oxley':'Oxley:1958zz','Hyman':'Hyman:1959zz','Goldansky':'Goldansky:1960zz',
       'Pugh':'Pugh:1957zz','Federspiel':'Federspiel:1991yd','McGibbon':'MacGibbon:1995in','Baranov':'Baranov:1974ec',
       'Zieger':'Zieger:1992jq','Bernardini':'Bernardini:1960wya','Hallin':'Hallin:1993ft'}
hdr = ['dataset','texkey','obs','E_lab_MeV','dE_MeV','theta_lab_deg','theta_cm_deg','frame_of_cs',
       'value','stat_err','sys_err_pt','norm_err_frac','units','source','notes']
rows = []
for r in csv.DictReader(open(sys.argv[1])):
    E, th = float(r['E_lab_MeV']), float(r['theta_lab_deg'])
    rows.append([f"KLP2018_{r['name']}", tex[r['name']], 'CS', E, '', th, f'{th_cm(E, th):.3f}', 'lab',
                 r['dsigma_nb_sr'], r['err_nb_sr'], 0, '', 'nb/sr',
                 'author notebook PWA/Final_NB/data_gather.nb (datafull) = PLB 782 (2018) 34 Table I',
                 'total error as used in KLP2018 fits; E may be PLB bin label'])
with open(os.path.join(sys.argv[2], 'db_klp2018.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(hdr); w.writerows(rows)
print(len(rows))
