!* test_dq.f90 -- stability self-check: all 80 tensor coefficients of the pi-N box of Compton scattering
!* (graph 1 of CovOneLoop: D(0,M^2,M^2,0; s,t; mpi^2,mpi^2,M^2,mpi^2), M = 1) at E_lab = 30 MeV on five
!* equally spaced t points 1.0-2.3 deg from the forward and from the backward direction, where the
!* Gram determinant almost vanishes. The coefficients are smooth there, so their 4th finite differences must
!* be tiny; plain double-precision LoopTools fails this by many orders of magnitude, LoopTools-dq passes.
program test_dq
  implicit none
#include "looptools.h"
  double precision, parameter :: mpi2 = (139.57039d0/938.27209d0)**2, dt = 3.d-7
  double precision s, t0(2), worst(2)
  double complex v(5)
  integer id, k, side
  call ltini
  call setdelta(0.d0)
  s = 1.d0 + 2.d0*30.d0/938.27209d0
  t0 = (/ 0.d0, -(s-1.d0)**2/s /)          ! forward t = 0, backward s*u = M^4
  worst = 0.d0
  do side = 1, 2
     do id = dd0, Ndd-2, 3
        do k = 1, 5
           v(k) = D0i(id, 0.d0, 1.d0, 1.d0, 0.d0, s, t0(side) + (3-2*side)*(-k*dt), mpi2, mpi2, 1.d0, mpi2)
        enddo
        worst(side) = max(worst(side), abs(v(1)-4*v(2)+6*v(3)-4*v(4)+v(5))/maxval(abs(v)))
     enddo
  enddo
  call ltexi
  write(*,'(a,2es10.2)') ' test_dq: max |4th difference|/|D| forward, backward:', worst
  if(maxval(worst) .gt. 1.d-6 .or. any(worst .ne. worst)) then
     write(*,*) 'test_dq: FAILED (unstable tensor coefficients)'
     stop 1
  endif
  write(*,*) 'test_dq: passed'
end program
