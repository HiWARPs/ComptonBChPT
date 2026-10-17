# LoopTools-dq: LoopTools 2.16 with a quad-precision core

A drop-in replacement for the double-precision LoopTools library used by the Compton code. The calls and
`looptools.h` stay the same; all one-loop integrals and tensor coefficients are evaluated in quad precision
(REAL*16) and returned in double precision.

## Why

LoopTools 2.16 (6 Nov 2024 tarball, MD5 `9503320a00e7c16b4ac1fd20749d3feb`) is the newest release and is the
version the code was built with. Its Passarino–Veltman reduction of the box tensor coefficients (up to rank 5
here) divides by Gram determinants. These vanish at θ = 0° and 180° and scale like ν⁴ at low photon energy, so
double precision fails there: exactly at 0° and 180° the result is NaN or ~1e55, within a few degrees of them the
cross sections are corrupted, and below about 8 MeV the loops are useless. code_v2 works around both problems
(`&v2 loopsafe`: t-interpolation near 0°/180°, static-polarisability form below 15 MeV). LoopTools-dq removes the
cause.

Raw loop evaluation (`loopsafe=f`), lab dσ in nb/sr and Σ2z:

| point | double | LoopTools-dq | safeguarded (default) |
|---|---|---|---|
| 50 MeV, 179° | 3312 | 19.957 | 19.957 |
| 100 MeV, 179° | 21.37 | 18.7165 | 18.7165 |
| 20 MeV, 2° | 50.47 | 23.0909 | 23.0909 |
| 0.5 MeV, 150° | 1.7×10⁶ | 20.5688 | 20.5688 |
| 1 MeV, 150° | 362.5 | 20.5283 | 20.5283 |
| Σ2z, 60 MeV, 179° | 0.220 | 0.631197 | 0.631197 |

With LoopTools-dq, the raw and the safeguarded evaluations agree to ≤ 2×10⁻⁹ near 0°/180° (20–150 MeV) and to
5×10⁻⁷ at 0.5–5 MeV. The second number is a genuine low-energy residual of the code's amplitude, independent of
`epsos`, Ξ and μ, and is not a LoopTools effect.

The default output (`loopsafe=t`) changes by ≤ 4×10⁻⁶ relative. All regression tests pass, and several
precision checks improve:

| check | double | dq |
|---|---|---|
| πN unitarity (8 multipoles, 160–350 MeV) | 1.4×10⁻⁶ | 4.2×10⁻¹¹ |
| small multipoles near threshold, per channel | 2×10⁻² | 1.7×10⁻⁹ |
| driver quadratic-form selftest | 2.8×10⁻⁸ | 1.4×10⁻¹⁵ |
| static πN polarisabilities, numeric vs analytic | 0.0032 | 0.0006 |

The cost is speed: the default run takes 1.8 s instead of 0.5 s (×3.5).

Exactly 0° and 180° remain singular in any precision: the Gram determinant is zero there. code_v2's interpolation
(flag 1) handles those points as before. Close to them, quad precision holds the rank-5 coefficients to ≤ 10⁻⁷ from
about 1° away (30 MeV). In double precision this needs ≳ 10°.

## Build and use

```
./build_looptools_dq.sh $HOME/.local/looptools-dq [LoopTools-2.16.tar.gz]   # ~1.5 min; FC=... selects the compiler
cd .. && make -f makefile.gfortran-native                                     # links $HOME/.local/looptools-dq
```

`makefile.gfortran-native` links LoopTools-dq by default (since 2026-10-18). For the plain double-precision
library, pass `LT=<its prefix>` (e.g. `$HOME/.local/looptools-gfortran-native`).

The script checks the MD5 of the official tarball and downloads it if no path is given. It builds LoopTools twice
with `./configure --quad`. The second build prefixes every public name with `q`, via cpp defines appended to
`externals.h` and `cexternals.h`; `D0` and `E0` are remapped in `defs.h` because a global define would also hit
literals such as `1D0`. The script then checks that no unprefixed public symbol is left, adds the double front end
`ltdq.f90` and runs `test_dq.f90`. That self-check takes the 4th differences of all 80 box coefficients 1–2° from
the forward and backward directions: plain double LoopTools gives 0.8 and 1.1, LoopTools-dq gives 2×10⁻¹¹ and 3×10⁻⁸.

## Limits

- Only the routines the Compton code calls are wrapped: `A0 B0 B1 DB1 B0i C0 C0i D0 D0i`, `ltini ltexi
  clearcache`, `set/getdelta mudim lambda uvdiv minmass versionkey`. Anything else (`Xget`/`Xval` cache access,
  complex-mass `..C` routines, five-point functions) gives a link error, never a silent precision mix.
- Requires a compiler with REAL*16; tested with gfortran 16.1 on macOS arm64. On Linux, add `-lquadmath` to the link
  line if the linker asks for it (untested).
- LoopTools (T. Hahn, with the FF library of G. J. van Oldenborgh) is licensed under the LGPL v3 (`COPYING` in the tarball). This directory ships no
  LoopTools source, only the build script and the front end.
