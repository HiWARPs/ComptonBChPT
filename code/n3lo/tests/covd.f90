! covd.f90 -- NNLO pi Delta loops of the 2015/v2 code (CorrIsoCovOneDeltaLoop, s + crossed u, MS-bar)
! reads lines "x t mpi2 gafac", writes A(1..8) (Re, Im)
program covd
  implicit none
  real*8 x, t, mpi2, gafac, y, dmloop
  integer ios
  complex*16 A(8), A0(8)
  common /v2loopmass/ dmloop
  call ltini
  dmloop = 0.d0
  do
     read(*,*,iostat=ios) x, t, mpi2, gafac
     if (ios /= 0) exit
     y = -x - t/2
     call setdelta(0.d0)
     call CorrIsoCovOneDeltaLoop(1.d0+2.d0*x, t, mpi2, gafac, A, 0.d0)
     call CorrIsoCovOneDeltaLoop(1.d0+2.d0*y, t, mpi2, gafac, A0, 0.d0)
     A(1) = A(1) + A0(1); A(2) = A(2) + A0(2); A(3:7) = A(3:7) - A0(3:7); A(8) = A(8) + A0(8)
     write(*,'(16es24.15)') (dble(A(ios)), dimag(A(ios)), ios=1,8)
  end do
end program
