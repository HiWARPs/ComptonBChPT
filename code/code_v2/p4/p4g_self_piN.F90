SUBROUTINE p4g_self_piN(x,t,dl,ss,mpi2,F1s,F2s,L0s,L1s,W)
  IMPLICIT NONE
  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2
  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s
  COMPLEX*16, INTENT(out) :: W(2)
  REAL*8 ixs,iys,iPq,PP,p1,p2,p3,p4,p12,p23,ma(4)
  COMPLEX*16 Z(2)
  COMPLEX*16 Tbb00,Tbb001,Tbb11,Tbb111,Ubb00,Ubb001
  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i
  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t
  p1 = ss
  ma(1) = mpi2
  ma(2) = 1.d0
  Tbb00 = B0i(7,p1,ma(1),ma(2))
  Tbb001 = B0i(13,p1,ma(1),ma(2))
  Tbb11 = B0i(10,p1,ma(1),ma(2))
  Tbb111 = B0i(16,p1,ma(1),ma(2))
  Ubb00 = B0i(7+1,p1,ma(1),ma(2))
  Ubb001 = B0i(13+1,p1,ma(1),ma(2))

      Z(1)=Ubb00 - 2*Tbb00
      Z(2)= - ss*Tbb11

      W(1) = 2*Z(1) + Z(2)

      Z(1)= - Ubb001 + Ubb00 - Tbb00 + 3*Tbb001
      Z(2)=Tbb11 + Tbb111
      Z(2)=ss*Z(2)
      Z(1)=2*Z(1) + Z(2)

      W(2) = ss*Z(1)

END SUBROUTINE p4g_self_piN
