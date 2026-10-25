SUBROUTINE dlg_dvaf_4(x,t,dl,ss,mpi2,MD,F1s,F2s,L0s,L1s,W)
  IMPLICIT NONE
  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2,MD
  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s
  COMPLEX*16, INTENT(out) :: W(3)
  REAL*8 ixs,iys,iPq,PP,ixD,iyD,iMD,y,p1,p2,p3,p4,p12,p23,ma(4)
  COMPLEX*16 Z(4)
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

      Z(1)=MD + 1
      Z(2)=Z(1)*U3bb00
      Z(2)=Z(2) + U3bb001
      Z(1)=Z(1)*T3bb00
      Z(1)=Z(1) + T3bb001
      Z(3)=Z(1) - Z(2)
      Z(1)=Z(1) - 2.D0/3.D0*Z(2)
      Z(1)=dl*Z(1)
      Z(1)=2*Z(3) + Z(1)

      W(1) = 8*Z(1)

      Z(1)=MD + 1
      Z(2)=Z(1)*U3bb00
      Z(2)=Z(2) + U3bb001
      Z(1)=Z(1)*T3bb00
      Z(1)=Z(1) + T3bb001
      Z(3)=7*Z(1) - 19.D0/3.D0*Z(2)
      Z(4)=2*Z(1) - 5.D0/3.D0*Z(2)
      Z(4)=dl*Z(4)
      Z(3)=2*Z(3) + Z(4)
      Z(3)=dl*Z(3)
      Z(1)=Z(1) - Z(2)
      Z(1)=8*Z(1) + 1.D0/3.D0*Z(3)

      W(2) = 4*Z(1)

      Z(1)= - 5.D0/3.D0*U3bb00 + 2*T3bb00
      Z(2)=MD + 1
      Z(1)=Z(1)*Z(2)
      Z(1)=Z(1) - 5.D0/3.D0*U3bb001 + 2*T3bb001
      Z(2)=4 + dl

      W(3) = 4.D0/3.D0*dl*Z(2)*Z(1)

END SUBROUTINE dlg_dvaf_4
