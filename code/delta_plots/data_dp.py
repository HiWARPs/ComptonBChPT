"""
data_dp.py -- display/chi^2 data table for the Delta-region figures.

Sources (read-only):
  * work/fits/db/fitdb.csv  -- assembled database with the frame / observable conventions and error models
                               used in the fits (fit_report Sect. 2); all sets except the GHP forward points.
  * work/data/db_old.csv    -- sets NOT in fitdb:
        Hallin1993_TabI    (SAL 93 Table I: angular distributions at 149/182/230/286 MeV, 15-MeV bins;
                            NOT independent of Table II -> used only as a chi^2 cross-check, not plotted)
        Baranov1966_refined (two forward points 235.0/55.4 deg cm and 239.5/69.5 deg cm; the 4 backward
                            points duplicate Baranov1966 T4 and are dropped, cf. db_old_README)
        Barbiellini1968    (Frascati Sigma at 90 deg cm, 300-335 MeV; source sign convention opposite to
                            LEGS/Babusci -> SIGN FLIPPED here for display, flagged)
  * Martel PhD thesis (UMass 2013) Table 7.5, 'Summed' column: Sigma_2x at 315.25-345.94 MeV, theta_lab =
    70..150 deg -- UNPUBLISHED (work/data/extracted/Martel2013_thesis.txt l. 7725ff., db_new_README Sect. 1).
Columns of the output table (results/data_display.csv):
  pid dataset label obs frame E Elo Ehi theta theta_lab theta_cm value err norm in_fitdb published note
  obs = driver observable (CS: lab dsigma/dOmega_lab, CSCM: c.m., S3, S2X, S2Z); frame = frame of theta.
"""
import os
import numpy as np
import pandas as pd
import config as C

MP = C.MP


def th_lab2cm(E, thl):
    b = E / (E + MP)
    c = np.cos(np.radians(thl))
    return np.degrees(np.arccos(np.clip((c - b) / (1 - b * c), -1, 1)))


def th_cm2lab(E, thc):
    b = E / (E + MP)
    c = np.cos(np.radians(thc))
    return np.degrees(np.arccos(np.clip((c + b) / (1 + b * c), -1, 1)))


def jac_cm_to_lab(E, theta_lab):
    """dsigma/dOmega_lab = J * dsigma/dOmega_cm  (massless outgoing photon), J = dOmega_cm/dOmega_lab"""
    b = E / (E + MP)
    c = np.cos(np.radians(theta_lab))
    return (1 - b * b) / (1 - b * c) ** 2


# legend label (experiment) per dataset; one label = one marker style
LABEL = {
    'LEGS2001_CS': 'LEGS 01', 'LEGS2001_S3': 'LEGS 01',
    'Hunger1997_75': 'MAMI 97', 'Hunger1997_90': 'MAMI 97', 'Wissmann1999': 'MAMI 99', 'Camen2002': 'MAMI 02',
    'Zieger1992': 'MAMI 92', 'OlmosdeLeon2001': 'TAPS 01', 'A2_2022_CS': 'A2 22', 'A2_2022_S3': 'A2 22',
    'A2_Martel2015_S2x': 'A2 15', 'Martel_thesis_S2x': 'Martel thesis (unpubl.)', 'A2_Paudyal2020_S2z': 'A2 20',
    'Hallin1993_TabII': 'SAL 93', 'Hallin1993_TabII_D': 'SAL 93', 'Hallin1993_TabI': 'SAL 93 (Tab. I)',
    'MacGibbon1995_tagged': 'SAL 95', 'MacGibbon1995_untagged': 'SAL 95', 'Kolb2000': 'SAL 00',
    'HIGS_2022_CIRC': 'HI$\\gamma$S 22', 'HIGS_2022_LIN': 'HI$\\gamma$S 22',
    'Genzel1972': 'Bonn 72', 'Gray1967': 'Illinois 67', 'Baranov1966': 'Moscow 66', 'Baranov1966_refined': 'Moscow 66',
    'Nagashima1964': 'Tokyo 64', 'Barbiellini1968': 'Frascati 68',
    'Federspiel1991': 'Illinois 91', 'Bernardini1960': 'Illinois 60', 'Baranov1974_75_th90': 'Moscow 74',
    'Baranov1974_75_th150': 'Moscow 74', 'Goldansky1960': 'Moscow 60', 'Pugh1957': 'MIT 57', 'Hyman1959': 'MIT 59',
    'Oxley1958': 'Chicago 58',
    'Galler2001_LARA': 'Bonn 01',
}
# references (for the notes / captions)
REF = {
    'LEGS 01': 'Blanpied et al., PRC 64 (2001) 025203 [Blanpied:2001ae]',
    'MAMI 97': 'Hünger et al., NPA 620 (1997) 385 [Hunger:1997pt] (75° and 90° c.m.)',
    'MAMI 99': 'Wissmann et al., NPA 660 (1999) 232 [Wissmann:1999vi] (130.7° lab)',
    'MAMI 02': 'Camen et al., PRC 65 (2002) 032202 / Kossert et al., EPJA 16 (2003) 259 [Camen:2001st, Kossert:2002ws] (136.2° lab)',
    'MAMI 92': 'Zieger et al., PLB 278 (1992) 34 [Zieger:1992jq]',
    'TAPS 01': 'Olmos de León et al., EPJA 10 (2001) 207 [deLeon:2001dnx]',
    'A2 22': 'Mornacchi et al., PRL 128 (2022) 132503 [A2:2021wmm]',
    'A2 15': 'Martel et al., PRL 114 (2015) 112501 [A2:2014iky]',
    'Martel thesis (unpubl.)': 'P. P. Martel, PhD thesis, UMass Amherst (2013), Table 7.5 (unpublished)',
    'A2 20': 'Paudyal et al., PRC 102 (2020) 035205 [A2:2019bqm]',
    'SAL 93': 'Hallin et al., PRC 48 (1993) 1497 [Hallin:1993ft] (Table II)',
    'SAL 93 (Tab. I)': 'Hallin et al., PRC 48 (1993) 1497 [Hallin:1993ft] (Table I)',
    'SAL 95': 'MacGibbon et al., PRC 52 (1995) 2097 [MacGibbon:1995in]',
    'SAL 00': 'Kolb et al., PRL 85 (2000) 1388 [Kolb:2000ix]',
    'HI$\\gamma$S 22': 'Li et al., PRL 128 (2022) 132502 [Li:2022vnz]',
    'Bonn 72': 'Genzel et al., Lett. Nuovo Cim. 4 (1972) 695 [Genzel:1972ba]',
    'Illinois 67': 'Gray and Hanson, PR 160 (1967) 1212 [Gray:1967zzd]',
    'Moscow 66': 'Baranov et al., Sov. J. Nucl. Phys. 3 (1966) 791; JETP 23 (1966) 242 [Baranov:1966aaa, Baranov:1966oqa]',
    'Tokyo 64': 'Nagashima, INS report 81 (1964) [Nagashima:1964zz]',
    'Frascati 68': 'Barbiellini et al., PR 174 (1968) 1665 [Barbiellini:1968rlx]',
    'Illinois 91': 'Federspiel et al., PRL 67 (1991) 1511 [Federspiel:1991yd]',
    'Illinois 60': 'Bernardini et al., Nuovo Cim. 18 (1960) 1203 [Bernardini:1960wya]',
    'Moscow 74': 'Baranov et al., PLB 52 (1974) 122; Yad. Fiz. 21 (1975) 689 [Baranov:1974ec, Baranov:1975ju]',
    'Moscow 60': 'Goldansky et al., Nucl. Phys. 18 (1960) 473',
    'MIT 57': 'Pugh et al., PR 105 (1957) 982 [Pugh:1957zz]',
    'MIT 59': 'Hyman et al., PRL 3 (1959) 93 [Hyman:1959zz]',
    'Chicago 58': 'Oxley, PR 110 (1958) 733 [Oxley:1958zz]',
    'Bonn 01': 'Galler et al., PLB 503 (2001) 245 [Galler:2001ht]; Wolf et al., EPJA 12 (2001) 231 [Wolf:2001ha] (LARA; '
               '75, 90, 125 deg c.m., digitized from the vector graphics of Fig. 4 of Galler et al.; plotting only)',
}
# marker style per label: (marker, color, filled)
TEAL, ORANGE, GREY, MAG, GREEN, PURPLE, BROWN = '#138d90', '#d9800b', '#3d3d3d', '#c2185b', '#2e7d32', '#6a3d9a', '#8c564b'
STYLE = {
    'LEGS 01': ('D', ORANGE, False), 'MAMI 97': ('o', TEAL, True), 'MAMI 99': ('^', TEAL, True),
    'MAMI 02': ('v', TEAL, True), 'MAMI 92': ('*', TEAL, True), 'TAPS 01': ('h', TEAL, True),
    'A2 22': ('s', TEAL, True), 'A2 15': ('o', TEAL, True), 'Martel thesis (unpubl.)': ('o', TEAL, False),
    'A2 20': ('s', TEAL, True),
    'SAL 93': ('s', GREY, False), 'SAL 93 (Tab. I)': ('s', GREY, False), 'SAL 95': ('o', GREY, False),
    'SAL 00': ('P', GREY, True),
    'HI$\\gamma$S 22': ('p', MAG, True),
    'Bonn 72': ('h', GREEN, False), 'Illinois 67': ('<', GREEN, True), 'Moscow 66': ('>', GREEN, True),
    'Tokyo 64': ('p', GREEN, False), 'Frascati 68': ('X', GREEN, True),
    'Illinois 91': ('^', PURPLE, False), 'Illinois 60': ('X', PURPLE, True), 'MIT 57': ('d', PURPLE, False),
    'MIT 59': ('d', PURPLE, True),
    'Moscow 74': ('v', BROWN, False), 'Moscow 60': ('<', BROWN, False), 'Chicago 58': ('>', BROWN, False),
    'Bonn 01': ('o', BROWN, True),
}
# legend order (Delta-region sets first)
ORDER = ['LEGS 01', 'MAMI 97', 'MAMI 99', 'MAMI 02', 'Bonn 01', 'SAL 93', 'SAL 00', 'Bonn 72', 'Illinois 67', 'Moscow 66',
         'Tokyo 64', 'Frascati 68', 'TAPS 01', 'A2 22', 'HI$\\gamma$S 22', 'SAL 95', 'Illinois 91', 'Moscow 74',
         'MIT 57', 'MIT 59', 'Moscow 60', 'Chicago 58', 'Illinois 60', 'MAMI 92', 'A2 15', 'Martel thesis (unpubl.)',
         'A2 20', 'SAL 93 (Tab. I)']

# energy bins of the double-polarisation sets (for horizontal bars in the vs-E figure)
EBINS = {'A2_Martel2015_S2x': (272.73, 303.32), 'Martel_thesis_S2x': (315.25, 345.94)}


def _f(x):
    try:
        v = float(x)
        return 0.0 if np.isnan(v) else v
    except Exception:
        return 0.0


def build():
    rows = []
    fdb = pd.read_csv(C.FITDB, keep_default_na=False)
    fdb = fdb[~fdb.dataset.str.startswith('GHP')]
    for _, r in fdb.iterrows():
        E = float(r.E)
        Elo, Ehi = EBINS.get(r.dataset, (np.nan, np.nan))
        if r.dataset == 'A2_Paudyal2020_S2z':
            Elo, Ehi = (265.0, 285.0) if E < 285 else (285.0, 305.0)
        rows.append(dict(pid=r.pid, dataset=r.dataset, label=LABEL[r.dataset], obs=r.dobs, frame=r.frame, E=E,
                         Elo=Elo, Ehi=Ehi, theta=float(r.theta), theta_lab=float(r.theta_lab),
                         theta_cm=float(r.theta_cm), value=float(r.value), err=float(r.err_pt), norm=float(r.norm),
                         in_fitdb=True, published=True, note=r.prov))
    old = pd.read_csv(os.path.join(C.DATA, 'db_old.csv'))
    # Hallin Table I (chi^2 cross-check only)
    for i, r in old[old.dataset == 'Hallin1993_TabI'].reset_index().iterrows():
        v, st, sy = float(r.value), _f(r.stat_err), _f(r.sys_err_pt)
        err = np.sqrt(st ** 2 + max(sy ** 2 - (0.036 * v) ** 2, 0.0))
        rows.append(dict(pid='Hallin1993_TabI_%03d' % (i + 1), dataset='Hallin1993_TabI', label=LABEL['Hallin1993_TabI'],
                         obs='CSCM', frame='CM', E=float(r.E_lab_MeV), Elo=np.nan, Ehi=np.nan, theta=float(r.theta_cm_deg),
                         theta_lab=float(r.theta_lab_deg), theta_cm=float(r.theta_cm_deg), value=v, err=err, norm=0.036,
                         in_fitdb=False, published=True,
                         note='SAL 93 Table I (15-MeV bins), not independent of Table II; c.m. angle and c.m. dsigma; '
                              'err = stat (+) [sys (-) 3.6% scale]; chi^2 cross-check only, not plotted'))
    # Baranov 1966 refined: the two forward points not contained in Baranov1966 T4
    br = old[old.dataset == 'Baranov1966_refined']
    br = br[br.theta_cm_deg < 80]
    for i, r in br.reset_index().iterrows():
        rows.append(dict(pid='Baranov1966_refined_%03d' % (i + 1), dataset='Baranov1966_refined',
                         label=LABEL['Baranov1966_refined'], obs='CSCM', frame='CM', E=float(r.E_lab_MeV), Elo=np.nan,
                         Ehi=np.nan, theta=float(r.theta_cm_deg), theta_lab=float(r.theta_lab_deg),
                         theta_cm=float(r.theta_cm_deg), value=float(r.value), err=_f(r.stat_err), norm=0.20,
                         in_fitdb=False, published=True,
                         note='Baranov 1966 refinement (JETP 23, 242), forward points only (the backward ones duplicate '
                              'Baranov1966 T4); c.m. frame assumed; norm 20% as Baranov1966'))
    # Barbiellini 1968, sign flipped to the LEGS/Babusci convention Sigma3 = (s_par - s_perp)/(s_par + s_perp)
    for i, r in old[old.dataset == 'Barbiellini1968'].reset_index().iterrows():
        rows.append(dict(pid='Barbiellini1968_%03d' % (i + 1), dataset='Barbiellini1968', label=LABEL['Barbiellini1968'],
                         obs='S3', frame='CM', E=float(r.E_lab_MeV), Elo=300.0, Ehi=335.0, theta=float(r.theta_cm_deg),
                         theta_lab=float(r.theta_lab_deg), theta_cm=float(r.theta_cm_deg), value=-float(r.value),
                         err=_f(r.stat_err), norm=0.0, in_fitdb=False, published=True,
                         note='HEPData ins55887: Sigma = %.2f(%.2f) in the source convention (opposite to LEGS); '
                              'sign flipped for display and chi^2' % (float(r.value), _f(r.stat_err))))
    # Martel thesis Table 7.5 (315.25-345.94 MeV), phi-summing method (as the published 273-303 MeV points)
    thesis = [(70, -0.160, 0.081), (90, -0.168, 0.038), (110, -0.202, 0.052), (130, -0.308, 0.114), (150, -0.074, 0.108)]
    Eb = (315.25, 345.94); E = 0.5 * (Eb[0] + Eb[1])
    for i, (th, v, e) in enumerate(thesis):
        rows.append(dict(pid='Martel_thesis_S2x_%03d' % (i + 1), dataset='Martel_thesis_S2x', label=LABEL['Martel_thesis_S2x'],
                         obs='S2X', frame='LAB', E=E, Elo=Eb[0], Ehi=Eb[1], theta=float(th), theta_lab=float(th),
                         theta_cm=float(th_lab2cm(E, th)), value=v, err=e, norm=0.0, in_fitdb=False, published=False,
                         note='Martel PhD thesis (UMass 2013) Table 7.5, phi-summing column; E bin 315.25-345.94 MeV '
                              '(above the 2pi0 threshold); error = stat (+) polarisation syst as for the published bin; '
                              'UNPUBLISHED'))
    df = pd.DataFrame(rows)
    # excluded: frames unknown (user decision 2026-10-14), see config.EXCLUDE_SETS
    df = df[~df.dataset.isin(C.EXCLUDE_SETS)].reset_index(drop=True)
    df['is_cs'] = df.obs.isin(['CS', 'CSCM'])
    return df


# plotting-only extension 400 < E <= 450 MeV (work/data/db_old_ext450.csv, built by
# work/data/scripts/build_db_old_ext450.py with the transcription code of db_old.csv): used ONLY in the figure of
# dsigma/dOmega vs E_gamma (Fig. 1 of the notes); never in the chi^2 tables.  Conventions as in work/fits fitdb
# (lib/database.py): c.m./lab frame from the source, error = stat (+) point syst, norm = quoted scale.
EXT450 = os.path.join(C.DATA, 'db_old_ext450.csv')
EXT_NAME = {'Hunger1997_th75cm': 'Hunger1997_75', 'Hunger1997_th90cm': 'Hunger1997_90', 'Wissmann1999': 'Wissmann1999',
            'Camen2002_Kossert2003': 'Camen2002'}


def build_ext450():
    old = pd.read_csv(EXT450)
    rows = []
    for i, r in old.iterrows():
        name = EXT_NAME[r.dataset]            # KeyError for any set not foreseen here (deliberate)
        cmframe = r.frame_of_cs == 'cm'
        orig_cm = str(r.notes).startswith('orig_angle=cm')
        fr = 'CM' if (cmframe or orig_cm) else 'LAB'
        thl, thc = float(r.theta_lab_deg), float(r.theta_cm_deg)
        rows.append(dict(pid='%s_ext450_%03d' % (name, i + 1), dataset=name, label=LABEL[name],
                         obs='CSCM' if cmframe else 'CS', frame=fr, E=float(r.E_lab_MeV), Elo=np.nan, Ehi=np.nan,
                         theta=thc if fr == 'CM' else thl, theta_lab=thl, theta_cm=thc, value=float(r.value),
                         err=float(np.hypot(_f(r.stat_err), _f(r.sys_err_pt))), norm=_f(r.norm_err_frac),
                         in_fitdb=False, published=True, note='db_old_ext450 (plotting only); ' + str(r.source)))
    df = pd.DataFrame(rows)
    df = df[~df.dataset.isin(C.EXCLUDE_SETS)].reset_index(drop=True)
    df['is_cs'] = df.obs.isin(['CS', 'CSCM'])
    return df


# plotting-only: the LARA (Bonn) cross sections of Galler et al. 2001 / Wolf et al. 2001 at theta_cm = 75, 90 and 125 deg
# (15-deg bins) from 265 to 780 MeV, never tabulated, recovered from the vector graphics of Fig. 4 of Galler et al.
# (work/data/db_lara_digitized.csv, built by work/data/scripts/extract_lara_galler2001_fig4.py).  Used ONLY in the
# figure of dsigma/dOmega vs E_gamma (Fig. 1 of the notes, to 800 MeV since 2026-10-17); never in the chi^2 tables.
LARA = os.path.join(C.DATA, 'db_lara_digitized.csv')


def build_lara():
    old = pd.read_csv(LARA)
    rows = []
    for i, r in old.iterrows():
        rows.append(dict(pid='Galler2001_LARA_%03d' % (i + 1), dataset='Galler2001_LARA', label=LABEL['Galler2001_LARA'],
                         obs='CSCM', frame='CM', E=float(r.E_lab_MeV), Elo=np.nan, Ehi=np.nan, theta=float(r.theta_cm_deg),
                         theta_lab=float(r.theta_lab_deg), theta_cm=float(r.theta_cm_deg), value=float(r.value),
                         err=float(np.hypot(_f(r.stat_err), _f(r.sys_err_pt))), norm=_f(r.norm_err_frac),
                         in_fitdb=False, published=True, note='db_lara_digitized (plotting only, digitized); ' + str(r.source)))
    df = pd.DataFrame(rows)
    df['is_cs'] = True
    return df


def driver_points(df):
    return [(r.pid, r.E, r.theta, r.obs, r.frame) for r in df.itertuples()]


def as_frame(df, frame):
    """unpolarised cross sections expressed in `frame` ('LAB' or 'CM'): returns (theta, value, err, converted)"""
    th = np.where(frame == 'LAB', df.theta_lab.values, df.theta_cm.values)
    J = jac_cm_to_lab(df.E.values, df.theta_lab.values)
    v = df.value.values.copy(); e = df.err.values.copy()
    conv = np.zeros(len(df), bool)
    if frame == 'LAB':
        m = (df.obs == 'CSCM').values
        v[m] *= J[m]; e[m] *= J[m]; conv |= m
    else:
        m = (df.obs == 'CS').values
        v[m] /= J[m]; e[m] /= J[m]; conv |= m
    return th, v, e, conv


if __name__ == '__main__':
    df = build()
    os.makedirs(C.RESULTS, exist_ok=True)
    df.to_csv(os.path.join(C.RESULTS, 'data_display.csv'), index=False)
    print(df.groupby('dataset').agg(N=('E', 'size'), Emin=('E', 'min'), Emax=('E', 'max'), obs=('obs', 'first'),
                                    frame=('frame', 'first'), in_fitdb=("in_fitdb", "first")).to_string())
