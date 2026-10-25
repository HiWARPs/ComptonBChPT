SUBROUTINE p4g_t_c2_XAB_L(x,t,dl,ss,mpi2,F1s,F2s,L0s,L1s,W)
  IMPLICIT NONE
  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2
  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s
  COMPLEX*16, INTENT(out) :: W(6)
  REAL*8 ixs,iys,iPq,PP,p1,p2,p3,p4,p12,p23,ma(4)
  COMPLEX*16 Z(6)
  COMPLEX*16 Taa0
  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i
  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t
  ma(1) = mpi2
  Taa0 = ma(1)*B0i(1,0.d0,0.d0,ma(1))

      Z(1)= - 1 + 1./4.*t
      Z(2)=2*x
      Z(3)=Z(2) + Z(1)
      Z(3)=t*Z(3)
      Z(4)= - 2 + x
      Z(4)=x*Z(4)
      Z(3)=4*Z(4) + Z(3)
      Z(3)=t*Z(3)
      Z(4)=x**2
      Z(3)= - 16*Z(4) + Z(3)
      Z(3)=t*Z(3)
      Z(4)=x + Z(1)
      Z(4)=t*Z(4)
      Z(2)= - Z(2) + 1./2.*Z(4)
      Z(2)=t**2*PP*Z(2)
      Z(1)=iPq*t**3*Z(1)*PP**2
      Z(1)=3*Z(2) + 1./2.*Z(1)
      Z(1)=iPq*Z(1)
      Z(1)=Z(3) + Z(1)

      W(1) = iPq**2*Taa0*Z(1)

      Z(1)=t*Taa0
      Z(2)= - Taa0 + 1./4.*Z(1)
      Z(3)=iPq**2
      Z(4)=Z(3)*Z(2)
      Z(5)=iPq**3
      Z(6)=PP*iPq**4
      Z(6)=Z(5) + 1./4.*Z(6)
      Z(2)=PP*Z(2)*Z(6)
      Z(2)=Z(4) + Z(2)
      Z(2)=Z(2)*t**4
      Z(1)= - Z(1) + 4*Taa0
      Z(3)=Z(1)*Z(3)
      Z(1)= - PP*Z(5)*Z(1)
      Z(1)= - 2*Z(3) + Z(1)
      Z(1)=Z(1)*t**3
      Z(3)=x*Z(3)*t**2
      Z(1)= - 4*Z(3) + Z(1)
      Z(1)=x*Z(1)

      W(2) = Z(1) + Z(2)


      W(3) =  0

      Z(1)=1./2.*t
      Z(2)=1./4.*t
      Z(3)= - 3*x - Z(2)
      Z(3)=Z(3)*Z(1)
      Z(4)=x**2
      Z(3)= - 6*Z(4) + Z(3)
      Z(3)=t*Z(3)
      Z(5)=x**3
      Z(3)= - 8*Z(5) + Z(3)
      Z(3)=Z(3)*t**2
      Z(5)= - x - 1./8.*t
      Z(1)=Z(1)*Z(5)
      Z(1)= - Z(4) + Z(1)
      Z(1)=t**3*PP*Z(1)
      Z(2)= - x - Z(2)
      Z(2)=iPq*t**4*Z(2)*PP**2
      Z(1)=3*Z(1) + 1./4.*Z(2)
      Z(1)=iPq*Z(1)
      Z(1)=Z(3) + Z(1)

      W(4) = iPq**2*Taa0*Z(1)

      Z(1)=x**2
      Z(2)=1./4.*t
      Z(3)= - 3*x - Z(2)
      Z(3)=t*Z(3)
      Z(3)= - 12*Z(1) + Z(3)
      Z(3)=t*Z(3)
      Z(4)=x**3
      Z(3)= - 16*Z(4) + Z(3)
      Z(3)=Z(3)*t**2
      Z(4)= - 2*x - Z(2)
      Z(4)=t*Z(4)
      Z(1)= - 4*Z(1) + Z(4)
      Z(1)=t**3*PP*Z(1)
      Z(2)= - x - Z(2)
      Z(2)=iPq*t**4*Z(2)*PP**2
      Z(1)=Z(1) + 1./4.*Z(2)
      Z(1)=iPq*Z(1)
      Z(1)=Z(3) + Z(1)

      W(5) = iPq**2*Taa0*Z(1)

      Z(1)=x**2
      Z(2)=1./4.*t
      Z(3)= - 3*x - Z(2)
      Z(3)=t*Z(3)
      Z(3)= - 12*Z(1) + Z(3)
      Z(3)=t*Z(3)
      Z(4)=x**3
      Z(3)= - 16*Z(4) + Z(3)
      Z(3)=Z(3)*t**2
      Z(4)= - 2*x - Z(2)
      Z(4)=t*Z(4)
      Z(1)= - 4*Z(1) + Z(4)
      Z(1)=t**3*PP*Z(1)
      Z(2)= - x - Z(2)
      Z(2)=iPq*t**4*Z(2)*PP**2
      Z(1)=Z(1) + 1./4.*Z(2)
      Z(1)=iPq*Z(1)
      Z(1)=Z(3) + Z(1)

      W(6) = iPq**2*Taa0*Z(1)

END SUBROUTINE p4g_t_c2_XAB_L
