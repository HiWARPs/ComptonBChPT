! covp3.f90 -- total O(p^3) covariant piN loop amplitude of the 2015/v2 code (CovOneLoop, s + crossed u)
! reads lines "x t mpi2 gafac", writes A(1..8) (Re, Im)
program covp3
  implicit none
  real*8 x, t, mpi2, gafac, y, mpin2_c, mnuc2_c
  integer ipimode, ios
  complex*16 A(8), A0(8)
  common /v2pimass/ mpin2_c, ipimode, mnuc2_c
  call ltini
  do
     read(*,*,iostat=ios) x, t, mpi2, gafac
     if (ios /= 0) exit
     mpin2_c = mpi2; ipimode = 0; mnuc2_c = 1.d0
     y = -x - t/2
     call setdelta(0.d0)
     call CovOneLoop(1.d0+2.d0*x, t, mpi2, gafac, A)
     call CovOneLoop(1.d0+2.d0*y, t, mpi2, gafac, A0)
     A(1) = A(1) + A0(1); A(2) = A(2) + A0(2); A(3:7) = A(3:7) - A0(3:7); A(8) = A(8) + A0(8)
     write(*,'(16es24.15)') (dble(A(ios)), dimag(A(ios)), ios=1,8)
  end do
end program
