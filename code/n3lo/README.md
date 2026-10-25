# Generator of the N3LO and N4LO graphs of proton Compton scattering in covariant BχPT

This directory produces the Fortran sources `code_v2/p4/` (the O(p^4) πN loops, N3LO) and `code_v2/dl/` (the Δ graphs
beyond N2LO: the photon on the Δ charge in the πΔ loops, N3LO; the γNΔ vertex inside the πN loops and the Δ pole
with the one-loop γNΔ vertex, N4LO), and it contains the independent check path (Python + LoopTools) and the tests
of these graphs.  The code in `code_v2` does not need anything from here to run; the generated sources are shipped.
The physics, the checks and the results are described in V. Pascalutsa, "Chiral EFT of nucleon Compton scattering
from N²LO to N⁴LO" (in preparation).

## Using the orders in the Fortran code

`code_v2`, namelist `&v2` of `in.dat` or the key `order=` of `code_v2/driver/run_points.sh`:

| order | content |
|---|---|
| `3.5` = `N2LO` (default) | Born, π⁰ pole, O(p³) πN loops, Δ pole, πΔ loops (Lensky, Pascalutsa 2010; Lensky, McGovern, Pascalutsa 2015) |
| `4` = `N3LO` | + O(p⁴) πN loops (Pauli insertions; c₁, c₂, c₃ graphs) + photon on the Δ charge in the πΔ loops |
| `4.5` = `N4LO` | + γNΔ vertex inside the πN loops + Δ pole with the one-loop γNΔ vertex (running g_M, g_E) |

Parts can be switched separately in `&v2`: `dloop=0|1`, `dgnd=0|1`, `irun=0|1|2`, `c1p4`, `c2p4`, `c3p4`, `p4part`
(see `code_v2/driver/README.md`).  Orders 4 and 4.5 are implemented for the proton.  A kinematic point costs about
1 s at orders 4 and 4.5 (5 ms without the Δ graphs, `dloop=0`); the reference static polarisabilities of the fit
driver are extracted numerically at start-up.

## Layout

| step | where |
|---|---|
| FORM headers: kinematics, polarisation sums, projections on A1..A6, tensor reduction in the LoopTools basis, D-expansion with the rational terms, on-shell vertex and self-energy projections; for the Δ graphs: Rarita–Schwinger tokens, contraction of the Δ indices, rank reduction by cancelling loop-momentum products against propagators per sub-topology | `form/common.h`, `form/dcommon.h` |
| graph generators: πN graphs of orders 3 and 4 and their subtraction terms; c_i graphs; Δ graphs (Feynman rules, enumeration of the photon placements, isospin) | `py/graphs.py`, `py/piN.py`, `py/tchan.py`, `py/delta.py`, `py/dgraphs.py` |
| FORM runs (cached in `out/`), evaluation with LoopTools-dq (`lt/`), renormalised amplitudes | `py/formtools.py`, `py/runform.py`, `py/evalamp.py`, `py/assemble.py`, `py/n3lo.py`, `py/devalamp.py`, `py/dassemble.py`, `py/ddelta.py`, `py/bhigh.py` |
| static polarisabilities (L'vov route of `code_v2/driver/static_pols.py`) | `py/statpol.py` |
| independent tree-level photoproduction (sum rules, optical theorem) | `py/photoprod.py`, `py/photoprod_delta.py` |
| Fortran generation (FORM O2 optimisation) | `py/build_fortran.py` → `fortran/p4/`, `py/build_fortran_delta.py` → `fortran/dl/` |

## Requirements and regeneration

FORM (version 5.0.2 was used; `form`, `tform`), Python 3 with numpy, sympy, mpmath, and LoopTools-dq (`code_v2/looptools`, installed
under `$HOME/.local/looptools-dq` or `$LT`).  Build the Python interface to LoopTools once, for example:

```
cd lt && gfortran -O2 -shared -fPIC -o libltwrap.dylib ltwrap.f90 -I$LT/include $LT/lib/libooptools.a
```

Then

```
python3 py/build_fortran.py              # -> fortran/p4/   (minutes with an empty FORM cache)
python3 py/build_fortran_delta.py af bl  # -> fortran/dl/   (hours with an empty FORM cache; FORM_BIN='tform -w8')
cp fortran/p4/*.F90 ../code_v2/p4/ ; cp fortran/dl/*.F90 ../code_v2/dl/
cd ../code_v2 && make -j 8 -f makefile.gfortran-native
```

The generators rewrite a source only when its content changes, so that `make` recompiles only what changed.

## Checks (tests/)

| test | result |
|---|---|
| `t_born*.py` | Born graphs = those of the code, 1e-15 |
| `t_p3.py` | O(p³) loops (pseudovector form) = the code's loops (pseudoscalar form) up to the code's on-shell offset, 1e-6 |
| `t_gi.py`, `t_tchan.py`, `t_dgi.py` | gauge invariance of every set, 1e-15 |
| `t_uv.py`, `t_c2pole.py`, `t_dA.py`, `t_dB.py` | UV structure: Pauli insertions, c₁, c₃ finite; c₂ and Δ graphs polynomial poles |
| `t_let.py`, `t_statpol.py`, `t_p4pols.py`, `t_dstatpol.py` | low-energy theorems; static polarisabilities (O(p³): analytic formulas of LMP2015 to 5e-4) |
| `t_sumrules.py`, `t_born_sr.py` | Baldin and γ₀ sum rules of the Pauli insertions with two tree-level codes |
| `t_hb*.py`, `t_kv*.py` | heavy-baryon limits (spin polarisabilities; c_i terms) |
| `t_dpole.py`, `t_lpdelta.py` | Δ pole and the N2LO πΔ loops of the code reproduced (1e-15, 1e-6; harness `build_harness.sh`) |
| `t_dB_opt.py` | optical theorem of the γNΔ-in-loop graphs vs `photoprod_delta.py`, 1e-6 |
| `t_dB_dr.py 1|2|3` | fixed-t dispersion relations (divided differences) of the three new sets, both cuts: 1e-5 |
| `t_fortran.py`, `t_fortran_dl.py`, `t_fortran_graphs.py` | generated Fortran = Python path, 1e-10 to 1e-13 (harnesses in `fortran/`) |
| `code_v2/tests/regression.py i j` | static polarisabilities of each set, additivity, order switch, low-energy limit, driver exactness |
