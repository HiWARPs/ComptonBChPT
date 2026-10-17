# Running Delta parameters: NNLO fit (F2_d_all_thB type)

Produced by `running2_fit.py`.  alpha, beta [1e-4 fm^3]; gamma_M1M1 held at P0[3]; chi^2/dof.
Suffixes: `_own` = P0 = own Pref of the running amplitude, covariance recomputed; `_P0` = P0 = default Pref (contact terms
restore the NNLO gammas), covariance recomputed; `_own_covdef` = default-configuration covariance; `F2_d_all_*` = no theory covariance.

| fit | N | alpha_E1 | beta_M1 | alpha+beta | alpha-beta | gamma_M1M1 (fixed) | chi^2/dof | d alpha | d beta | d chi^2 |
|---|---|---|---|---|---|---|---|---|---|---|
| F2_d_all_thB | 283 | 11.546(0.453) | 2.978(0.418) | 14.524(0.248) | 8.568(0.836) | 2.939 | 220.73/263 | +0.000 | +0.000 | +0.00 |
| F2_d_all_thB_running2_own | 283 | 11.455(0.439) | 3.216(0.406) | 14.671(0.234) | 8.239(0.813) | 2.564 | 227.04/263 | -0.091 | +0.238 | +6.32 |
| F2_d_all_thB_running2_P0 | 283 | 11.478(0.441) | 3.164(0.408) | 14.642(0.234) | 8.314(0.817) | 2.939 | 221.89/263 | -0.068 | +0.186 | +1.17 |
| F2_d_all_thB_running2_own_covdef | 283 | 11.492(0.451) | 3.202(0.416) | 14.694(0.247) | 8.290(0.831) | 2.564 | 225.27/263 | -0.054 | +0.224 | +4.54 |
| F2_d_all_thB_running3_own | 283 | 11.532(0.445) | 3.117(0.412) | 14.650(0.232) | 8.415(0.825) | 2.939 | 221.73/263 | -0.014 | +0.139 | +1.00 |
| F2_d_all_thB_running3_P0 | 283 | 11.532(0.445) | 3.117(0.412) | 14.650(0.232) | 8.415(0.825) | 2.939 | 221.73/263 | -0.014 | +0.139 | +1.00 |
| F2_d_all_thB_running3_own_covdef | 283 | 11.566(0.452) | 3.107(0.418) | 14.673(0.248) | 8.459(0.835) | 2.939 | 220.69/263 | +0.020 | +0.129 | -0.04 |
| F2_d_all | 283 | 11.393(0.245) | 3.536(0.236) | 14.929(0.072) | 7.857(0.476) | 2.939 | 813.64/263 | +0.000 | +0.000 | +0.00 |
| F2_d_all_running2_own | 283 | 11.742(0.244) | 3.743(0.235) | 15.486(0.071) | 7.999(0.474) | 2.564 | 1206.10/263 | +0.349 | +0.207 | +392.46 |
| F2_d_all_running2_P0 | 283 | 11.413(0.245) | 3.705(0.236) | 15.118(0.073) | 7.707(0.476) | 2.939 | 760.37/263 | +0.019 | +0.169 | -53.27 |
| F2_d_all_running3_own | 283 | 11.489(0.245) | 3.652(0.236) | 15.140(0.073) | 7.837(0.476) | 2.939 | 762.95/263 | +0.095 | +0.116 | -50.69 |
| F2_d_all_running3_P0 | 283 | 11.489(0.245) | 3.652(0.236) | 15.140(0.073) | 7.837(0.476) | 2.939 | 762.95/263 | +0.095 | +0.116 | -50.69 |

running2: truncation sigma_th recomputed: max relative change 1.9e-01 (mean -7.8e-02).

running3: truncation sigma_th recomputed: max relative change 2.1e-01 (mean -5.7e-02).
