"""Forward (theta=0) 'data points' from the empirical sum-rule evaluation of
Gryniuk, Hagelstein, Pascalutsa, PRD 92 (2015) 074031 [Gryniuk:2015eza] (Re f)
and PRD 94 (2016) 034043 [Gryniuk:2016gnm] (Re g); values transcribed from
Table 'tab:values_forward_amps' of arXiv:1604.00789 (units mub*GeV).
Below pion threshold (nu < 150.7 MeV) the amplitudes are real, so
  dsigma/dOmega_lab(0) = |f|^2 + |g|^2 ,  Sigma_2z(0) = -2 Re(f g*)/(|f|^2+|g|^2)
[Eqs. (dsi), (Si2z) of arXiv:1604.00789].  Re f = average of fit I and fit II,
with half their difference added in quadrature to the quoted error.
"""
import csv, math, sys, os
MUB_GEV2 = 389.379  # 1 GeV^-2 = 389.379 mub ; |f|^2 [mub] = f[mub GeV]^2 / 389.379
Mp = 938.272
rows = [  # nu MeV, Re f (fit I), err, Re f (fit II), err, Re g, err
 (50, -2.84,.01, -2.85,.01, -0.26,.01),
 (60, -2.76,.01, -2.76,.01, -0.32,.01),
 (70, -2.65,.01, -2.66,.01, -0.37,.01),
 (79, -2.54,.01, -2.56,.01, -0.42,.01),
 (90, -2.38,.01, -2.40,.01, -0.48,.01),
 (100,-2.21,.01, -2.23,.01, -0.54,.01),
 (110,-2.00,.01, -2.03,.01, -0.59,.01),
 (120,-1.75,.01, -1.78,.01, -0.64,.01),
 (130,-1.45,.02, -1.48,.02, -0.68,.01),
 (135,-1.27,.02, -1.30,.02, -0.69,.01),
 (145,-0.79,.03, -0.83,.03, -0.65,.02),
 (150,-0.39,.03, -0.42,.03, -0.52,.02),
]
out = sys.argv[1] if len(sys.argv) > 1 else '.'
hdr = ['dataset','texkey','obs','E_lab_MeV','dE_MeV','theta_lab_deg','theta_cm_deg','frame_of_cs',
       'value','stat_err','sys_err_pt','norm_err_frac','units','source','notes']
src = 'arXiv:1604.00789 Table (tab:values_forward_amps)'
amp, obs = [], []
for nu, f1, e1, f2, e2, g, eg in rows:
    f = 0.5*(f1+f2); ef = math.sqrt((0.5*(e1+e2))**2 + (0.5*(f1-f2))**2)
    amp.append(['GHP_forward','Gryniuk:2015eza','FWD_REF',nu,'',0,0,'na',f'{f:.4f}',f'{ef:.4f}',0,0,'mub*GeV',src,
                f'Re f, mean of fit I ({f1}) and fit II ({f2}); err = quoted (+) half-difference'])
    amp.append(['GHP_forward','Gryniuk:2016gnm','FWD_REG',nu,'',0,0,'na',f'{g:.4f}',f'{eg:.4f}',0,0,'mub*GeV',src,'Re g'])
    cs = 1000*(f*f+g*g)/MUB_GEV2                                   # nb/sr, lab
    ecs = 1000*math.sqrt((2*f*ef)**2 + (2*g*eg)**2)/MUB_GEV2
    D = f*f+g*g; s2z = -2*f*g/D
    ds_df = (-2*g*D + 2*f*g*2*f)/D**2; ds_dg = (-2*f*D + 2*f*g*2*g)/D**2
    es2z = math.hypot(ds_df*ef, ds_dg*eg)
    obs.append(['GHP_forward','Gryniuk:2016gnm','CS',nu,'',0,0,'lab',f'{cs:.4f}',f'{ecs:.4f}',0,0,'nb/sr',src,
                'derived: |f|^2+|g|^2 ; errors of different energies strongly correlated (same sigma_abs fits)'])
    obs.append(['GHP_forward','Gryniuk:2016gnm','S2Z',nu,'',0,0,'na',f'{s2z:.5f}',f'{es2z:.5f}',0,0,'1',src,
                'derived: -2 f g/(f^2+g^2); Sigma_2z sign convention of arXiv:1604.00789 (3/2 - 1/2)'])
for name, data in (('db_forward_amp.csv', amp), ('db_forward_obs.csv', obs)):
    with open(os.path.join(out, name), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(hdr); w.writerows(data)
    print('wrote', name, len(data))
