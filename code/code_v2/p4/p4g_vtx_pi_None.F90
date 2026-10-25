SUBROUTINE p4g_vtx_pi_None(x,t,dl,ss,mpi2,F1s,F2s,L0s,L1s,W)
  IMPLICIT NONE
  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2
  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s
  COMPLEX*16, INTENT(out) :: W(3)
  REAL*8 ixs,iys,iPq,PP,p1,p2,p3,p4,p12,p23,ma(4)
  COMPLEX*16 Z(6)
  COMPLEX*16 Tcc00,Tcc0000,Tcc001,Tcc0011,Tcc0012,Tcc002,Tcc0022,Tcc1,Tcc11,Tcc111,Tcc1112,Tcc112,Tcc1122,Tcc12,Tcc122,Tcc1222,Tcc222,Tcc2222,Ucc00,Ucc0000,Ucc001,Ucc0012,Ucc002,Ucc0022
  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i
  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t
  p1 = -dl
  p2 = 1
  p12 = 1
  ma(1) = mpi2
  ma(2) = mpi2
  ma(3) = 1.d0
  Tcc00 = C0i(10,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc0000 = C0i(40,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc001 = C0i(22,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc0011 = C0i(43,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc0012 = C0i(46,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc002 = C0i(25,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc0022 = C0i(49,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc1 = C0i(4,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc11 = C0i(13,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc111 = C0i(28,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc1112 = C0i(55,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc112 = C0i(31,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc1122 = C0i(58,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc12 = C0i(16,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc122 = C0i(34,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc1222 = C0i(61,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc222 = C0i(37,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc2222 = C0i(64,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc00 = C0i(10+1,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc0000 = C0i(40+1,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc001 = C0i(22+1,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc0012 = C0i(46+1,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc002 = C0i(25+1,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc0022 = C0i(49+1,p1,p2,p12,ma(1),ma(2),ma(3))

      Z(1)= - 2 - dl
      Z(1)=Tcc002*Z(1)
      Z(2)=dl - 2
      Z(3)=Tcc0000*Z(2)
      Z(4)=Tcc00*dl
      Z(1)= - Z(4) + Z(1) - Z(3)
      Z(3)=Tcc001 + Tcc0011 + Tcc0012
      Z(4)=2*dl
      Z(3)=Z(3)*Z(2)*Z(4)
      Z(5)=Tcc112 + Tcc12
      Z(4)=Z(4)*Z(5)
      Z(5)= - Tcc1122 + Tcc122 - Tcc1222
      Z(5)=dl*Z(5)
      Z(4)=Tcc2222 + Z(5) + Z(4)
      Z(5)=dl + 4
      Z(4)=Z(5)*Z(4)
      Z(6)= - Ucc0022 + 2*Ucc002
      Z(5)=Z(5)*Z(6)
      Z(2)=Ucc0000*Z(2)
      Z(2)=2*Z(2) + Z(5)
      Z(5)=6 + dl
      Z(5)=Tcc0022*Z(5)
      Z(1)=6*Z(5) + 2*Z(2) + Z(4) + Z(3) + 12*Z(1)

      W(1) = 4*Z(1)

      Z(1)=Tcc112 + Tcc12
      Z(2)=2*dl
      Z(1)=Z(1)*Z(2)
      Z(3)= - Ucc0022 + 2*Ucc002
      Z(4)= - Tcc1122 + Tcc122 - Tcc1222
      Z(4)=dl*Z(4)
      Z(1)=2*Z(3) + Tcc2222 + Z(4) + Z(1)
      Z(3)=dl + 8
      Z(3)=Z(3)*dl
      Z(3)=Z(3) + 16
      Z(1)=Z(3)*Z(1)
      Z(3)= - Tcc00 - Tcc001 - Tcc0011 - Tcc0012
      Z(4)=4*dl
      Z(3)=Z(4)*Z(3)
      Z(3)=24*Tcc0000 - 8*Ucc0000 + Z(3)
      Z(4)=dl + 4
      Z(3)=Z(4)*Z(3)
      Z(2)=17 + Z(2)
      Z(2)=dl*Z(2)
      Z(2)=36 + Z(2)
      Z(2)=Tcc0022*Z(2)
      Z(4)= - 7 - dl
      Z(4)=dl*Z(4)
      Z(4)= - 12 + Z(4)
      Z(4)=Tcc002*Z(4)
      Z(1)=8*Z(4) + 4*Z(2) + Z(3) + Z(1)

      W(2) = 2*Z(1)

      Z(1)=Ucc00 - Ucc0012
      Z(2)=4*dl
      Z(1)=Z(2)*Z(1)
      Z(2)=Tcc0022 + Ucc001
      Z(3)=8*dl
      Z(2)=Z(3)*Z(2)
      Z(3)= - Tcc002 + Ucc002 - Ucc0022
      Z(4)=2*dl
      Z(3)=Z(4)*Z(3)
      Z(4)=Tcc12 + Tcc112 - Tcc1122
      Z(5)=Tcc1 - Tcc1112
      Z(4)=6*Tcc11 + 4*Tcc111 + 2*Z(5) + 3*Z(4)
      Z(4)=Z(4)*dl**2
      Z(5)= - 10*Tcc00 - 20*Tcc001 + 16*Tcc0012 + Tcc222 + Tcc2222
      Z(5)=dl*Z(5)
      Z(1)=Z(4) + Z(5) + Z(3) + Z(2) + Z(1)
      Z(2)=dl + 4
      Z(1)=Z(2)*Z(1)
      Z(2)= - 2 - dl
      Z(2)=dl*Z(2)
      Z(2)=8 + Z(2)
      Z(2)=Tcc1222*dl*Z(2)
      Z(1)=Z(2) + Z(1)

      W(3) = 2*Z(1)

END SUBROUTINE p4g_vtx_pi_None
