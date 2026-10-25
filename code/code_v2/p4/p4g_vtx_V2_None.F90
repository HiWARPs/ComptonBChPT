SUBROUTINE p4g_vtx_V2_None(x,t,dl,ss,mpi2,F1s,F2s,L0s,L1s,W)
  IMPLICIT NONE
  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2
  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s
  COMPLEX*16, INTENT(out) :: W(3)
  REAL*8 ixs,iys,iPq,PP,p1,p2,p3,p4,p12,p23,ma(4)
  COMPLEX*16 Z(2)
  COMPLEX*16 Tbb00,Tbb11,Ubb00
  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i
  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t
  p1 = 1
  ma(1) = mpi2
  ma(2) = 1.d0
  Tbb00 = B0i(7,p1,ma(1),ma(2))
  Tbb11 = B0i(10,p1,ma(1),ma(2))
  Ubb00 = B0i(7+1,p1,ma(1),ma(2))

      Z(1)= - 4*Tbb00 + 2*Ubb00 - Tbb11
      Z(2)=dl - 2

      W(1) = 4*Z(2)*Z(1)

      Z(1)=4*Tbb00 - 2*Ubb00 + Tbb11
      Z(2)=dl + 4

      W(2) = 4*Z(2)*Z(1)


      W(3) =  0

END SUBROUTINE p4g_vtx_V2_None
