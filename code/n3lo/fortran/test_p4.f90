program test_p4
  use p4mod
  implicit none
  real*8 x, t, gA, F, mpi2, M
  complex*16 A(8)
  integer ios, part
  M=938.27d0
  gA=1.27d0; F=92.21d0/M; mpi2=(139.57d0/M)**2
  call ltini
  call p4_setup(gA,F,mpi2,1.793d0,-1.913043d0,-0.99d0*M*1.d-3,1.38d0*M*1.d-3,-2.33d0*M*1.d-3)
  write(*,'(a,8es14.6)') '# vc1:', vc1(1:4)
  write(*,'(a,8es14.6)') '# vc2:', vc2(1:4)
  write(*,'(a,4es14.6)') '# sL0,sL1:', sL0, sL1
  do
     read(*,*,iostat=ios) x, t, part
     if (ios /= 0) exit
     p4_part=part
     call CovOneLoopP4(x,t,A)
     write(*,'(12es24.15)') (dble(A(ios)), dimag(A(ios)), ios=1,6)
  end do
end program
