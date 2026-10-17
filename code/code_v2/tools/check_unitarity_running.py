#!/usr/bin/env python3
"""P33 unitarity (Watson's theorem) of the complete amplitude with the running Delta parameters (&v2 irun = 1, 2, 3).

Without running (regression test c): Im f_MM^{1+}(code) - Im f_MM^{1+}(piN loops) = k (|M_B + M_Delta|^2 - |M_B|^2) to 1e-9,
with M_B the I=3/2 Born M1+ multipole (CG-normalised) and M_Delta = -sqrt(a M_D Gamma/k) G/D the photoproduction Delta-pole
multipole built from the code's own Delta pole, G = gM - kappa gE + i Im G (Im G = M_B sqrt(k M_D Gamma/a) by the Watson
construction of the vertex Im parts), D = s - M_D^2 + i M_D Gamma, a(W) = alpha w^2 (W+M)^2 (W+M_D)/(8 (M+M_D)^2).

With running real parts the code's s-channel Delta pole is f_Delta = -a'(W) Z(s) G(s)^2/D(s) with G(s) = Re gM(s) - kappa Re gE(s)
+ i Im G, D(s) = s - M_D(s)^2 + i M_D Gamma(s), Z = 1 + Re tau_R and a' = a (W + M_D(s))/(W + M_D) (running numerator mass);
the imaginary part changes only through these running quantities.  Watson's theorem then holds exactly for a photoproduction
amplitude M_Delta = -sqrt(a' M_D Gamma Z/k) G(s)/D(s) PROVIDED Im G = M_B sqrt(k M_D Gamma/(a' Z)); the code keeps Im G at its
Z = 1, a' = a value (the vertex Im parts are NLO, the running is beyond the order), so the Born-Delta interference term is off
by the factor sqrt(a/(a' Z)): residual = 2 a' Im G (1 - sqrt(a'Z/a)) sqrt(a'Z/a) [Re G Re D + Im G Im D]/|D|^2 ... This script quantifies:
  (1) the identity f_Delta(code, s-channel Delta pole only) = -a' Z G^2/D with the running quantities of As.f90 (tools/running_table);
  (2) the unitarity residual of the complete amplitude relative to the generalised Watson prediction (code's Im G), the algebra
      check that the residual vanishes with Im G rescaled by sqrt(a/(a'Z)), and the residual relative to the naive prediction
      with fixed couplings (as in test c);
for irun = 0 (reference), 1, 2, 3 at E_lab = 200, 250, 300, 330, 360 MeV.
"""
import os, sys, subprocess, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
from pwa_run import run_pwa
from delta_pole_pw import Dfun, edem, MD, SMM
from born_multipoles import calibrated

AEM = 1 / 137.036
MPI = 0.13957
ES = [200.0, 250.0, 300.0, 330.0, 360.0]


def running_table(Ws_GeV, vxmode=2):
    """columns per (W, irun): ReG1 ReG2 Mrun ztau ImgM ImgE  (nucleon-mass units) from As.f90 delta_running"""
    d = tempfile.mkdtemp(prefix='rt_', dir=os.environ.get('V2SCRATCH'))
    exe = os.path.join(d, 'running_table')
    env = dict(os.environ, SDKROOT=subprocess.run(['xcrun', '--sdk', 'macosx', '--show-sdk-path'], capture_output=True, text=True).stdout.strip())
    subprocess.run(['/opt/homebrew/bin/gfortran-16', '-O2', '-ffree-line-length-none', '-std=legacy', os.path.join(HERE, 'running_table.f90'),
                    os.path.join(ROOT, 'build-gfortran-native', 'As.o'), '-o', exe], check=True, env=env)
    out = subprocess.run([exe, str(vxmode)] + ['%.10f' % w for w in Ws_GeV], capture_output=True, text=True).stdout
    R = {}
    for row in np.loadtxt(out.splitlines()):
        R[(round(row[0] * SMM, 6), int(round(row[1])))] = row[2:]
    return R


def main(verbose=True):
    B = calibrated(kp=0.0, kn=0.0)
    Ws = [np.sqrt(SMM ** 2 + 2 * SMM * E / 1000) for E in ES]
    RT = running_table(Ws)
    piN = run_pwa(ES, oneD=False, anomaly=False, delta=(False, False))
    worst = {}
    for ir in (0, 1, 2, 3):
        v2 = {'irun': ir}
        full = run_pwa(ES, v2=v2)
        dpole = run_pwa(ES, nucl=(False, False), anomaly=False, oneN=False, oneD=False, delta=(True, False), v2=v2)   # s-channel pole only
        if verbose:
            print('irun=%d' % ir)
        dev = dict(identity=0.0, watson_Z=0.0, watson_exact=0.0, naive=0.0)
        for E, w in zip(ES, Ws):
            Wm = w / SMM
            Eq = (w * w + MPI ** 2 - SMM ** 2) / (2 * w); k = np.sqrt(Eq * Eq - MPI ** 2)
            m32, _, _ = B.mult32(w); MB = np.sqrt(2 / 3) * m32['M1+'].real
            om = (Wm * Wm - 1) / (2 * Wm)
            a = AEM * om ** 2 * (Wm + 1) ** 2 * (Wm + MD) / (8 * (1 + MD) ** 2) / SMM
            kap = (Wm - 1) / (2 * (Wm + 1)); MG = MD * edem(Wm)                 # M_D Gamma(s) [M units]
            if ir == 0:
                ReG1, ReG2, Mrun, Z = 2.97, -1.0, MD, 1.0
            else:
                ReG1, ReG2, Mrun, Z, _, _ = RT[(round(w, 6), ir)]
            ap = a * (Wm + Mrun) / (Wm + MD)
            GR = ReG1 - kap * ReG2
            GI = MB * np.sqrt(k * MG / a)                                      # code's Watson Im part (Z = 1, a with M_D)
            D = Wm * Wm - Mrun ** 2 + 1j * MG
            G = GR + 1j * GI
            fD = -ap * Z * G * G / D                                           # expected Delta-pole multipole
            got = dpole[E]['MM1+']
            dev['identity'] = max(dev['identity'], abs(got / fD - 1))
            sc = np.sqrt(a / (ap * Z))                                        # Im G -> M_B sqrt(k M_D Gamma/(a' Z)): exact Watson
            ImfD = full[E]['MM1+'].imag - piN[E]['MM1+'].imag
            MDel = -np.sqrt(ap * MG * Z / k) * G / D
            pw = k * (abs(MB + MDel) ** 2 - abs(MB) ** 2)
            Gx = GR + 1j * GI * sc
            MDelx = -np.sqrt(ap * MG * Z / k) * Gx / D
            fDx = -ap * Z * Gx * Gx / D
            px = k * (abs(MB + MDelx) ** 2 - abs(MB) ** 2)
            MDel0 = -np.sqrt(a * MG / k) * (2.97 + kap * 1.0 + 1j * MB * np.sqrt(k * MG / a)) / Dfun(Wm)
            p0 = k * (abs(MB + MDel0) ** 2 - abs(MB) ** 2)
            dev['watson_Z'] = max(dev['watson_Z'], abs(ImfD / pw - 1))
            dev['watson_exact'] = max(dev['watson_exact'], abs(fDx.imag / px - 1))
            dev['naive'] = max(dev['naive'], abs(ImfD / p0 - 1))
            if verbose:
                print('  E=%3.0f W=%.4f  ReG=%.4f ImG=%.4f M_D(s)=%.4f Z=%.4f | |f_Delta code/formula-1| = %.1e | Im f_Delta / Watson(running) '
                      '- 1 = %+.2e [Im G mismatch sqrt(a/(a\'Z)) = %.4f; algebra check %+.1e] | / naive(fixed couplings) - 1 = %+.2e'
                      % (E, w, GR, GI, Mrun, Z, abs(got / fD - 1), ImfD / pw - 1, sc, fDx.imag / px - 1, ImfD / p0 - 1))
        worst[ir] = dev
        print('irun=%d: max |f_Delta(code)/(-a\' Z G^2/D) - 1| = %.1e; unitarity residual vs generalised Watson = %.2e '
              '(algebra check with Im G rescaled by sqrt(a/(a'Z)): %.1e); vs naive fixed-coupling prediction = %.2e'
              % (ir, dev['identity'], dev['watson_Z'], dev['watson_exact'], dev['naive']))
    return worst


if __name__ == '__main__':
    main()
