! covdsec.f90 -- sections of CorrIsoCovOneDeltaLoop (Box, KR, KP, VR, SE, TP), s + crossed u
! reads "x t mpi2 gafac isec" (isec=0: all), writes A(1..8)
program covdsec
  implicit none
  real*8 x, t, mpi2, gafac, y, dmloop, secf(6)
  integer ios, isec
  complex*16 A(8), A0(8)
  common /v2loopmass/ dmloop
  common /secfac/ secf
  call ltini
  dmloop = 0.d0
  do
     read(*,*,iostat=ios) x, t, mpi2, gafac, isec
     if (ios /= 0) exit
     secf = 0.d0
     if (isec == 0) then
        secf = 1.d0
     else
        secf(isec) = 1.d0
     end if
     y = -x - t/2
     call setdelta(0.d0)
     call SecCorrIso(1.d0+2.d0*x, t, mpi2, gafac, A, 0.d0)
     call SecCorrIso(1.d0+2.d0*y, t, mpi2, gafac, A0, 0.d0)
     A(1) = A(1) + A0(1); A(2) = A(2) + A0(2); A(3:7) = A(3:7) - A0(3:7); A(8) = A(8) + A0(8)
     write(*,'(16es24.15)') (dble(A(ios)), dimag(A(ios)), ios=1,8)
  end do
end program
