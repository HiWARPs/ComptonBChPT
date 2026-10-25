! test_dl.f90 -- generated Delta-loop routines (dl/CovDeltaLoops.F90) vs the Python check path (tests/t_fortran_dl.py)
! reads "x t part" (part 1: CovDeltaFull, 2: CovDeltaGND), writes A(1..6) (Re, Im)
program test_dl
  implicit none
  real*8 x, t, gA, F, mpi2, M, MD
  complex*16 A(8)
  integer ios, part
  M=938.27d0
  gA=1.27d0; F=92.21d0/M; mpi2=(139.57d0/M)**2; MD=1232.d0/M
  call ltini
  call dl_setup(gA,F,2.85d0,MD,2.97d0,-1.d0,mpi2)
  do
     read(*,*,iostat=ios) x, t, part
     if (ios /= 0) exit
     if (part == 1) call CovDeltaFull(x,t,A)
     if (part == 2) call CovDeltaGND(x,t,A)
     write(*,'(12es24.15)') (dble(A(ios)), dimag(A(ios)), ios=1,6)
  end do
end program
