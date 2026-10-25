SUBROUTINE p4g_vtx_int_mag(x,t,dl,ss,mpi2,F1s,F2s,L0s,L1s,W)
  IMPLICIT NONE
  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2
  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s
  COMPLEX*16, INTENT(out) :: W(3)
  REAL*8 ixs,iys,iPq,PP,p1,p2,p3,p4,p12,p23,ma(4)
  COMPLEX*16 Z(9)
  COMPLEX*16 Tcc00,Tcc0000,Tcc001,Tcc0011,Tcc0012,Tcc002,Tcc0022,Tcc111,Tcc1111,Tcc1112,Tcc112,Tcc1122,Tcc12,Tcc122,Tcc1222,Tcc222,Tcc2222,Ucc00,Ucc0000,Ucc001,Ucc0011,Ucc0012,Ucc002,Ucc0022
  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i
  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t
  p1 = 1
  p2 = -dl
  p12 = 1
  ma(1) = mpi2
  ma(2) = 1.d0
  ma(3) = 1.d0
  Tcc00 = C0i(10,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc0000 = C0i(40,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc001 = C0i(22,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc0011 = C0i(43,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc0012 = C0i(46,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc002 = C0i(25,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc0022 = C0i(49,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc111 = C0i(28,p1,p2,p12,ma(1),ma(2),ma(3))
  Tcc1111 = C0i(52,p1,p2,p12,ma(1),ma(2),ma(3))
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
  Ucc0011 = C0i(43+1,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc0012 = C0i(46+1,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc002 = C0i(25+1,p1,p2,p12,ma(1),ma(2),ma(3))
  Ucc0022 = C0i(49+1,p1,p2,p12,ma(1),ma(2),ma(3))

      Z(1)=Tcc001 + Tcc002
      Z(2)=6*dl
      Z(1)=Z(2)*Z(1)
      Z(3)= - Ucc001 - Ucc002
      Z(4)=2*dl
      Z(3)=Z(4)*Z(3)
      Z(5)=Tcc111 + Tcc222
      Z(5)=dl*Z(5)
      Z(6)= - Tcc12*dl**2
      Z(1)=Z(6) + Z(3) + Z(1) + Z(5)
      Z(3)=dl + 4
      Z(1)=Z(3)*Z(1)
      Z(5)=Tcc1112 + Tcc1222
      Z(2)=Z(2)*Z(5)
      Z(5)=12*dl
      Z(6)= - Ucc0012*Z(5)
      Z(7)=36*dl
      Z(8)=Tcc0012*Z(7)
      Z(2)=Z(8) + Z(6) + Z(2)
      Z(6)=dl + 2
      Z(2)=Z(6)*Z(2)
      Z(6)= - Ucc0011 - Ucc0022
      Z(5)=Z(5)*Z(6)
      Z(6)=Tcc0011 + Tcc0022
      Z(6)=Z(7)*Z(6)
      Z(3)=dl*Z(3)
      Z(3)=6 + Z(3)
      Z(3)=Tcc1122*Z(3)
      Z(3)=Tcc1111 + Tcc2222 + Z(3)
      Z(7)=3*dl
      Z(3)=Z(7)*Z(3)
      Z(7)=Tcc112 + Tcc122
      Z(8)=dl + 7
      Z(8)=Z(8)*dl
      Z(8)=Z(8) + 12
      Z(7)=Z(8)*Z(7)
      Z(7)=72*Tcc0000 - 24*Ucc00 - 60*Ucc0000 + Z(7)
      Z(7)=dl*Z(7)
      Z(8)=2 - dl
      Z(4)=Tcc00*Z(8)*Z(4)
      Z(1)=Z(4) + Z(1) + Z(7) + Z(2) + Z(3) + Z(6) + Z(5)

      W(1) = 4*Z(1)

      Z(1)=Tcc112 + Tcc122
      Z(2)=dl + 11
      Z(2)=Z(2)*dl
      Z(2)=Z(2) + 40
      Z(2)=Z(2)*dl
      Z(2)=Z(2) + 48
      Z(1)=Z(2)*Z(1)
      Z(2)=Tcc111 + Tcc222
      Z(3)=dl + 8
      Z(3)=Z(3)*dl
      Z(4)=Z(3) + 16
      Z(2)=Z(4)*Z(2)
      Z(5)=24*Tcc0000 + Tcc1111 + Tcc2222 - 8*Ucc00 - 20*Ucc0000
      Z(6)=dl + 4
      Z(5)=Z(6)*Z(5)
      Z(3)=22 + Z(3)
      Z(3)=dl*Z(3)
      Z(3)=24 + Z(3)
      Z(3)=Tcc1122*Z(3)
      Z(1)=Z(3) + Z(5) + Z(2) + Z(1)
      Z(1)=dl*Z(1)
      Z(2)= - Ucc0011 - Ucc0022
      Z(3)=4*dl
      Z(2)=Z(3)*Z(2)
      Z(5)=Tcc0011 + Tcc0022
      Z(7)=12*dl
      Z(5)=Z(7)*Z(5)
      Z(2)=Z(5) + Z(2)
      Z(2)=Z(6)*Z(2)
      Z(5)= - Ucc001 - Ucc002
      Z(5)=Z(4)*Z(5)
      Z(6)=Tcc1112 + Tcc1222
      Z(8)=dl + 6
      Z(8)=Z(8)*dl
      Z(8)=Z(8) + 8
      Z(6)=Z(8)*Z(6)
      Z(9)=10 + dl
      Z(9)=dl*Z(9)
      Z(9)=24 + Z(9)
      Z(9)=Tcc00*Z(9)
      Z(5)=Z(9) + Z(6) + Z(5)
      Z(6)=2*dl
      Z(5)=Z(6)*Z(5)
      Z(6)=Tcc001 + Tcc002
      Z(9)=6*dl
      Z(6)=Z(9)*Z(6)
      Z(9)=Tcc12*dl**2
      Z(6)=Z(9) + Z(6)
      Z(4)=Z(4)*Z(6)
      Z(3)= - Ucc0012*Z(3)
      Z(6)=Tcc0012*Z(7)
      Z(3)=Z(6) + Z(3)
      Z(3)=Z(8)*Z(3)
      Z(1)=Z(1) + Z(5) + Z(2) + Z(3) + Z(4)

      W(2) = 2*Z(1)


      W(3) =  0

END SUBROUTINE p4g_vtx_int_mag
