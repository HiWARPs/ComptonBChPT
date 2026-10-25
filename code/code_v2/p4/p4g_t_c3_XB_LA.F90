SUBROUTINE p4g_t_c3_XB_LA(x,t,dl,ss,mpi2,F1s,F2s,L0s,L1s,W)
  IMPLICIT NONE
  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2
  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s
  COMPLEX*16, INTENT(out) :: W(6)
  REAL*8 ixs,iys,iPq,PP,p1,p2,p3,p4,p12,p23,ma(4)
  COMPLEX*16 Z(4)
  COMPLEX*16 Tbb00
  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i
  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t
  p1 = 0
  ma(1) = mpi2
  ma(2) = mpi2
  Tbb00 = B0i(7,p1,ma(1),ma(2))

      Z(1)=PP*t**2
      Z(2)=PP*t
      Z(3)= - iPq*Z(2)**2
      Z(1)= - 2*Z(1) + Z(3)
      Z(3)= - t + 4
      Z(3)=Tbb00*Z(3)
      Z(1)=Z(3)*Z(1)
      Z(3)=8*Z(3)
      Z(2)= - x*Z(2)*Z(3)
      Z(1)=Z(2) + Z(1)
      Z(1)=Z(1)*iPq**3
      Z(1)= - Z(3) + Z(1)

      W(1) = 2*Z(1)

      Z(1)= - t + 4
      Z(1)=Tbb00*Z(1)
      Z(2)=Z(1)*iPq**3
      Z(2)=3*Z(2)
      Z(3)= - PP*Z(1)*iPq**4
      Z(3)= - Z(2) + Z(3)
      Z(3)=PP*Z(3)
      Z(1)=Z(1)*iPq**2
      Z(3)= - 2*Z(1) + Z(3)
      Z(3)=Z(3)*t**3
      Z(2)= - PP*Z(2)
      Z(2)= - 4*Z(1) + Z(2)
      Z(2)=Z(2)*t**2
      Z(1)=x*t*Z(1)
      Z(1)= - 8*Z(1) + Z(2)
      Z(1)=x*Z(1)

      W(2) = 4*Z(1) + Z(3)


      W(3) =  0

      Z(1)=PP*iPq**4
      Z(2)=iPq**3
      Z(3)= - 4*Z(2) - Z(1)
      Z(4)=4*PP
      Z(3)=Z(4)*x*Z(3)
      Z(1)= - 2*Z(2) - Z(1)
      Z(1)=t*PP*Z(1)
      Z(1)=Z(3) + Z(1)
      Z(1)=t*Z(1)
      Z(2)= - x**2*Z(2)*Z(4)
      Z(2)= - 1 + Z(2)
      Z(1)=8*Z(2) + Z(1)
      Z(1)=t*Z(1)
      Z(1)= - 32*x + Z(1)

      W(4) = t*Tbb00*Z(1)

      Z(1)=x**2
      Z(2)= - 12*x - t
      Z(2)=t*Z(2)
      Z(2)= - 48*Z(1) + Z(2)
      Z(2)=t*Z(2)
      Z(3)=x**3
      Z(2)= - 64*Z(3) + Z(2)
      Z(2)=t*Z(2)
      Z(3)= - 8*x - t
      Z(3)=t*Z(3)
      Z(1)= - 16*Z(1) + Z(3)
      Z(1)=t**2*PP*Z(1)
      Z(3)= - 4*x - t
      Z(3)=iPq*t**3*Z(3)*PP**2
      Z(1)=3*Z(1) + Z(3)
      Z(1)=iPq*Z(1)
      Z(1)=2*Z(2) + Z(1)

      W(5) = iPq**2*Tbb00*Z(1)

      Z(1)=x**2
      Z(2)= - 12*x - t
      Z(2)=t*Z(2)
      Z(2)= - 48*Z(1) + Z(2)
      Z(2)=t*Z(2)
      Z(3)=x**3
      Z(2)= - 64*Z(3) + Z(2)
      Z(2)=t*Z(2)
      Z(3)= - 8*x - t
      Z(3)=t*Z(3)
      Z(1)= - 16*Z(1) + Z(3)
      Z(1)=t**2*PP*Z(1)
      Z(3)= - 4*x - t
      Z(3)=iPq*t**3*Z(3)*PP**2
      Z(1)=3*Z(1) + Z(3)
      Z(1)=iPq*Z(1)
      Z(1)=2*Z(2) + Z(1)

      W(6) = iPq**2*Tbb00*Z(1)

END SUBROUTINE p4g_t_c3_XB_LA
