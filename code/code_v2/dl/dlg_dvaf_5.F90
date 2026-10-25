SUBROUTINE dlg_dvaf_5(x,t,dl,ss,mpi2,MD,F1s,F2s,L0s,L1s,W)
  IMPLICIT NONE
  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2,MD
  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s
  COMPLEX*16, INTENT(out) :: W(3)
  REAL*8 ixs,iys,iPq,PP,ixD,iyD,iMD,y,p1,p2,p3,p4,p12,p23,ma(4)
  COMPLEX*16 Z(8)
  COMPLEX*16 T3bb00,T3bb001,U3bb00,U3bb001
  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i
  COMPLEX*16 a0
  ! topology 3 (B)
  p1 = 1
  ma(1) = mpi2
  ma(2) = MD**2
  T3bb00 = B0i(7,p1,ma(1),ma(2))
  U3bb00 = B0i(7+1,p1,ma(1),ma(2))
  T3bb001 = B0i(13,p1,ma(1),ma(2))
  U3bb001 = B0i(13+1,p1,ma(1),ma(2))

      Z(1)=MD**(-1)
      Z(2)=2*Z(1) + 1 - MD
      Z(2)=U3bb00*Z(2)
      Z(2)=Z(2) - U3bb001
      Z(3)=T3bb00*MD
      Z(4)= - T3bb00 - T3bb001 + 4.D0/3.D0*U3bb001
      Z(4)=Z(1)*Z(4)
      Z(2)=Z(4) + Z(3) + T3bb001 + 2.D0/3.D0*Z(2)
      Z(2)=dl*Z(2)
      Z(3)= - U3bb00 + T3bb00
      Z(4)=MD + 1
      Z(3)=Z(4)*Z(3)
      Z(3)=T3bb001 - U3bb001 + Z(3)
      Z(2)=2*Z(3) + Z(2)

      W(1) = 8*Z(2)

      Z(1)=MD**(-1)
      Z(2)=38 + 5*dl
      Z(3)=1.D0/9.D0*dl
      Z(2)=Z(2)*Z(3)
      Z(4)=Z(2) + 8
      Z(4)= - U3bb00*Z(4)
      Z(5)=dl + 7
      Z(6)=1.D0/3.D0*dl
      Z(5)=Z(5)*Z(6)
      Z(5)=Z(5) + 4
      Z(5)=2*Z(5)
      Z(7)=T3bb00*Z(5)
      Z(4)=Z(4) + Z(7)
      Z(4)=MD*Z(4)
      Z(7)=dl + 4
      Z(7)=Z(7)*Z(1)
      Z(8)= - Z(6)*Z(7)
      Z(5)=Z(5) + Z(8)
      Z(5)=T3bb001*Z(5)
      Z(8)= - 22 - dl
      Z(3)=Z(8)*Z(3)
      Z(8)=dl*Z(7)
      Z(8)= - 8 + 4.D0/9.D0*Z(8)
      Z(3)=Z(3) + Z(8)
      Z(3)=U3bb00*Z(3)
      Z(7)= - Z(7) + 10 + dl
      Z(6)=Z(6)*Z(7)
      Z(6)=8 + Z(6)
      Z(6)=T3bb00*Z(6)
      Z(2)= - Z(2) + Z(8)
      Z(2)=U3bb001*Z(2)
      Z(2)=Z(5) + Z(2) + Z(4) + Z(3) + Z(6)

      W(2) = 4*Z(2)

      Z(1)=MD**(-1)
      Z(2)=1 + 5*MD
      Z(3)=1.D0/3.D0*U3bb00
      Z(2)=Z(2)*Z(3)
      Z(3)=U3bb00 + U3bb001
      Z(3)=T3bb00 + T3bb001 - 4.D0/3.D0*Z(3)
      Z(3)=Z(3)*Z(1)
      Z(4)=1 + 2*MD
      Z(4)=Z(4)*T3bb00
      Z(2)=5.D0/3.D0*U3bb001 - 2*T3bb001 + Z(3) + Z(2) - Z(4)
      Z(3)= - 4 - dl

      W(3) = 4.D0/3.D0*dl*Z(3)*Z(2)

END SUBROUTINE dlg_dvaf_5
