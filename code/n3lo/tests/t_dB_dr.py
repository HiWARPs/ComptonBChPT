"""t_dB_dr.py -- analyticity check of the gamma N Delta loop graphs (Bloop, CovDeltaGND of fortran/test_dl.exe), of
the pi Delta loops with the complete minimal coupling (A_full, CovDeltaFull) and of the O(p^4) piN loops (test_p4.exe): once-subtracted fixed-t dispersion relations
of the crossing-even L'vov amplitudes below the pion threshold,
    Re A^L_i(nu,t) - A^L_i(0,t) = (2 nu^2/pi) int Im A^L_i(nu',t) / (nu' (nu'^2 - nu^2)) dnu',
with both s-channel cuts (pi N from nu_piN; pi Delta from nu_piDelta).  The unsubtracted relation of gamma_0 (= a_4) is
not used: its UV pole is not exactly zero, Im A^L_4 tends to a constant and the integral grows logarithmically.
t = -2e-3 M^2 (t = 0 is Gram-singular); nu' up to NUMAX (LoopTools at larger nu' loses precision), tail ~ 1/NUMAX^2."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
from statpol import lvov
EXE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'fortran', 'test_dl.exe')
EXE4 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'fortran', 'test_p4.exe')   # part 3: O(p^4) piN loops
m, MD, NUMAX = 139.57 / 938.27, 1232.0 / 938.27, 2000.0
t = float(os.environ.get('TDR', '-2e-3'))                       # fixed t [M^2]
part = int(sys.argv[1]) if len(sys.argv) > 1 else 2          # 1: A_full, 2: Bloop, 3: O(p^4) piN loops (CovOneLoopP4)


def amps(nus):
    inp = '\n'.join('%.15e %.15e %d' % (nu - t / 4, t, part if part < 3 else 0) for nu in nus)          # x = nu - t/4
    out = subprocess.run([EXE if part < 3 else EXE4], input=inp, capture_output=True, text=True).stdout
    rows = [np.array([float(z) for z in l.split()]) for l in out.split('\n') if len(l.split()) == 12]
    return np.array([lvov(r[0::2] + 1j * r[1::2], nu, t) for r, nu in zip(rows, nus)])


nu_pin = ((1 + m) ** 2 - 1) / 2 + t / 4; nu_pid = ((MD + m) ** 2 - 1) / 2 + t / 4
z, w = np.polynomial.legendre.leggauss(80); z = 0.5 * (z + 1); w = 0.5 * w
L = np.log(NUMAX / nu_pid)
# pi N cut: nu = nu_pin + (nu_pid - nu_pin) z^2 ; pi Delta cut: nu = nu_pid exp(L z^2) (sqrt onsets)
nus = np.r_[nu_pin + (nu_pid - nu_pin) * z ** 2, nu_pid * np.exp(L * z ** 2)]
wts = np.r_[w * 2 * (nu_pid - nu_pin) * z, w * nu_pid * np.exp(L * z ** 2) * 2 * L * z]
cut = np.r_[np.zeros(80), np.ones(80)].astype(int)
ImA = amps(nus).imag
nl = np.linspace(0.02, 0.12, 11); low = amps(nl).real                      # A^L_i(0,t) from nu -> 0
fit = np.array([np.polyfit(nl ** 2, low[:, i], 4) for i in range(6)]); a0, a2 = fit[:, -1], fit[:, -2]
print('part %d (%s), t = %g M^2; |Im A^L_i| / nu^2 at nu = %.0f M: %s' % (part, ['A_full', 'Bloop', 'O(p^4) piN'][part - 1], t, nus[-1],
      ' '.join('%.1e' % v for v in np.abs(ImA[-1]) / nus[-1] ** 2)))
worst = 0.0
for nu0 in (0.06, 0.10, 0.14):
    re = amps([nu0]).real[0] - a0
    k = 2 * nu0 ** 2 / np.pi * wts / (nus * (nus ** 2 - nu0 ** 2))
    dr = np.array([(k * ImA[:, i] * (cut == c)).sum() for c in (0, 1) for i in range(6)]).reshape(2, 6)
    scale = np.maximum(np.abs(re), 1e-3 * np.max(np.abs(re)))
    dev = np.abs(dr.sum(0) - re) / scale
    worst = max(worst, dev.max())
    print('nu = %.2f M:  Re A^L(nu) - A^L(0) = %s' % (nu0, ' '.join('%10.3e' % v for v in re)))
    print('              DR (pi N + pi Delta) = %s' % ' '.join('%10.3e' % v for v in dr.sum(0)))
    print('              of it pi Delta cut   = %s   rel. dev. max %.1e' % (' '.join('%10.3e' % v for v in dr[1]), dev.max()))
print('worst relative deviation %.1e (once subtracted)' % worst)
# extrapolation-free form: divided differences in nu^2 of Re A^L at physical points below threshold, n-th order
# F[x_1..x_{n+1}] = (2/pi) int nu' Im A^L(nu') / prod_k (nu'^2 - x_k) dnu'  (a polynomial of degree < n in nu^2 drops out)
def divdiff(xs, fs):
    fs = list(fs)
    for lev in range(1, len(xs)):
        fs = [(fs[j + 1] - fs[j]) / (xs[j + lev] - xs[j]) for j in range(len(fs) - 1)]
    return fs[0]


worst2 = {}
for pts in ((0.05, 0.10, 0.145), (0.04, 0.08, 0.12, 0.145)):
    x = np.array(pts) ** 2
    re = divdiff(x, amps(list(pts)).real)
    k = 2 / np.pi * wts * nus / np.prod([nus ** 2 - xk for xk in x], axis=0)
    dr = np.array([(k * ImA[:, i]).sum() for i in range(6)])
    dev = np.abs(dr - re) / np.abs(re)
    worst2[len(pts) - 1] = dev
    print('order-%d divided difference at nu = %s:' % (len(pts) - 1, pts))
    print('   Re = %s' % ' '.join('%11.4e' % v for v in re))
    print('   DR = %s   rel. dev. %s' % (' '.join('%11.4e' % v for v in dr), ' '.join('%.1e' % v for v in dev)))
