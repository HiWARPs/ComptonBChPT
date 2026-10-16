!!*** NAMING NOTE (added 2026-09-07) ***
!!* looptools.h is a C-PREPROCESSOR header, not a Fortran one: the 237 tensor-
!!* coefficient indices are cpp object-like macros (the lower-case b-b-0-0 spelling
!!* expands to 7, b-b-1-1 to 10, and so on).  cpp is CASE-SENSITIVE; Fortran is
!!* CASE-INSENSITIVE.  This file used to store the B-function results in variables
!!* whose names differed from those macros only in case, which worked solely
!!* because cpp left the upper-case spellings untouched.  They are now prefixed
!!* Bf... (Bf00, Bf11, Bf001, Bf111, DBf1o, DEBf00, ...) so that no Fortran
!!* variable here shares a name with a LoopTools macro in any case.
!!* Verified semantics-preserving: preprocessed source and compiled .text section
!!* are byte-identical to the previous version.
!!*
!!* Two things to keep in mind when editing:
!!*   - do not name a local after a LoopTools index, in any case;
!!*   - do not compile this file without -cpp, and do not rename it to .f90 --
!!*     the lower-case index arguments would silently become undeclared
!!*     implicit COMPLEX*16 variables instead of integer constants.
SUBROUTINE CovOneLoop(s,t,mpi2,gafac,A)
  USE trigd
!!   IMPLICIT COMPLEX*16 B, C, D
  implicit COMPLEX*16 (B-D)
#include "looptools.h"
   REAL*8, INTENT(in) :: s,t,mpi2,gafac
   COMPLEX*16,INTENT(out) :: A(8)
   COMPLEX*16, PARAMETER :: CI=(0.d0,1.d0)
   REAL*8 x,nuclm
   REAL*8 m1s,m2s,m3s,m4s
   REAL*8 p1s,p2s,p2so,p3s,p4s,p1p2s,p1p2so,p2p3s,ps,pso
   INTEGER i,j
   COMPLEX*16 Sum1, Sum2, Sum3
   !* [v2] physical pion masses per isospin channel (proton target; &v2 pimass='phys'):
   !*   photons on the (charged) pion  -> pi+ n intermediate states: mass^2 = mpi2 (argument)
   !*   photons on the nucleon         -> pi0 p intermediate states: mass^2 = mpin2
   !*   nucleon self-energy (graph 4)  -> both: isospin weights 2 (pi+ n) : 1 (pi0 p)
   !*   two-photon/isoscalar pi-pi graphs -> charged pions only (mpi2)
   !* ipimode=0 (default): mpin2 = mpi2, i.e. the v1 single-mass calculation (unchanged code path).
   !*   with physical masses the pi+ n loops also use the neutron mass for the internal nucleon
   !*   (mN2c = (M_n/M_p)^2), so that the pi+ n cusp sits at the physical threshold (151.4 MeV)
   REAL*8 mpin2_c, mpin2, mnuc2_c, mN2c
   INTEGER ipimode
   COMMON /v2pimass/ mpin2_c, ipimode, mnuc2_c
   COMPLEX*16 Bc0,Bc1,Bc0o,Bc1o,DBc1o
   !*   ipimode=2 / 3 (diagnostics): only the pi0 p / only the pi+ n channel (graph weights wn, wc)
   REAL*8 wc, wn
   COMPLEX*16 Asv(8)

   A=0.d0
   mpin2=mpi2
   mN2c=1.d0                 ! = nuclm**2 (M_p^2 in nucleon-mass units)
   IF(ipimode.GE.1) THEN
      mpin2=mpin2_c
      mN2c=mnuc2_c
   ENDIF
   wc=1.d0; wn=1.d0
   IF(ipimode.EQ.2) wc=0.d0
   IF(ipimode.EQ.3) wn=0.d0

   Sum1=0.d0
   Sum2=0.d0
   Sum3=0.d0

   nuclm=1.d0   !* it's assumed that all is given in nucleon-mass units
   x=(s-1.d0)/2.d0

   p1s=0.0d0
   p2s=nuclm**2
   p3s=nuclm**2
   p4s=0.0d0
   p2p3s=t
   p1p2s=s

   !* Graph 1:
   !*  Box graph with photons attached to the pions
   m1s=mpi2
   m4s=m1s
   m2s=m1s
   m3s=mN2c                  ! [v2] internal neutron (pi+ n channel); = M^2 by default

   D00=D0i(dd00,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12=D0i(dd12,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D13=D0i(dd13,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22=D0i(dd22,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D123=D0i(dd123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D112=D0i(dd112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D122=D0i(dd122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D113=D0i(dd113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D222=D0i(dd222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D001=D0i(dd001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D002=D0i(dd002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D003=D0i(dd003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
      

   A(1) = + D002 * (-8.D0) + t*D00 * (-2.D0) + t*D002 * (-2.D0) + t*D001 * (-4.D0)&
        + x*D00 * (-8.D0) + x*D002 * (-8.D0) + x*D001 * (-16.D0)
   
   A(2) = + D222 * ( 2.D0 ) + D123 * ( 8.D0 ) + D122 * ( 8.D0 ) + t*D22 * ( 1.D0/2.D0 )&
        + t*D13 * ( 2.D0 ) + t*D12 * ( 2.D0 ) + t*D222 * ( 1.D0/2.D0 ) + t*D123 * ( 4.D0 )&
        + t*D122 * ( 3.D0 ) + t*D113 * ( 4.D0 ) + t*D112 * ( 2.D0 ) + x*D22 * ( 2.D0 )&
        + x*D13 * ( 8.D0 ) + x*D12 * ( 8.D0 ) + x*D222 * ( 2.D0 ) + x*D123 * ( 16.D0 )&
        + x*D122 * ( 12.D0 ) + x*D113 * ( 16.D0 ) + x*D112 * ( 8.D0 )

   A(3) = 0.d0

   A(4) = D00 * (  - 4.D0 ) + D002 * (  - 4.D0 ) + D001 * (  - 8.D0 )

   A(5) = D002 * (  - 2.D0 ) + D001 * (  - 4.D0 )

   A(6) = D002 * ( 2.D0 ) + D001 * ( 4.D0 )

   A(7) = D22 * (  - 1.D0 ) + D13 * (  - 4.D0 )+ D12 * (  - 4.D0 )+ D222 * (  - 1.D0 )&
        + D123 * (  - 8.D0 ) + D122 * (  - 6.D0 ) + D113 * (  - 8.D0 ) + D112 * (  - 4.D0 )

   A(8) = 0.d0


!   write(*,*) 'Managed to calculate graph 1'
!   RETURN

   !* Graph 2:
   !*  Sum of reducible graphs with photons attached to the pions, renormalized
   p1s=0.d0
   p2s=s
   p2so=nuclm**2+1.D-6
   p1p2s=nuclm**2

   m1s=mpi2
   m2s=m1s
   m3s=mN2c                  ! [v2] pi+ n channel

   C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C22=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)

   C00o=C0i(cc00,p1s,p2so,p1p2s,m1s,m2s,m3s)
   C12o=C0i(cc12,p1s,p2so,p1p2s,m1s,m2s,m3s)
   C22o=C0i(cc22,p1s,p2so,p1p2s,m1s,m2s,m3s)

!   Sum1=4.D0*C22o+4.D0*C00o
!   write(*,*) 4.D0*C22o, 4.D0*C00o
!   Sum3=-4.D0*C12o

   A(1) = A(1) &
        +x**(-1)*t*C00o*(1.D0)+x**(-1)*t*C00*(-1.D0)+C22o*(4.D0)+C00o*(4.D0)+C00*(-4.D0)

   A(2) = A(2) &
        +x**(-1)*C22o*(2.D0)+x**(-1)*C22*(-2.D0)

   A(3) = A(3) &
        +x**(-1)*t*C00o*(1.D0)+x**(-1)*t*C00*(-1.D0)+C22o*(4.D0)+C00o*(4.D0)+C00*(-4.D0)

   A(4) = A(4) &
        +x**(-1)*C00o*(2.D0)+x**(-1)*C00*(-2.D0)

   A(5) = A(5) &
        +x**(-1)*C00o*(-2.D0)+x**(-1)*C00*(2.D0)+C12*(-1.D0)

   A(6) = A(6) &
        +x**(-1)*C22o*(1.D0)+x**(-1)*C00o*(2.D0)+x**(-1)*C22*(-1.D0)+x**(-1)*C00*(-2.D0)+C12*(1.D0)

   A(7) = A(7)

   A(8) = A(8) &
        +x**(-1)*C00o*(-2.D0)+x**(-1)*C00*(2.D0)

!   write(*,*) 'Managed to calculate graph 2'
   IF(ipimode.GE.2) A=wc*A          ! [v2] graphs 1+2: pi+ n channel
   Asv=A

   !* Graph 3:
   !*  Sum of reducible graphs with photons attached to the nucleons inside the loops, renormalized
   p1s=nuclm**2
   p2s=0.d0
   p1p2s=s
   p1p2so=nuclm**2+1.D-6
   m1s=mpin2                 ! [v2] pi0 p channel
   m2s=nuclm**2
   m3s=nuclm**2

   C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1=C0i(cc1,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C2=C0i(cc2,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C11=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C22=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)

   C00o=C0i(cc00,p1s,p2s,p1p2so,m1s,m2s,m3s)
   C1o=C0i(cc1,p1s,p2s,p1p2so,m1s,m2s,m3s)
   C2o=C0i(cc2,p1s,p2s,p1p2so,m1s,m2s,m3s)
   C11o=C0i(cc11,p1s,p2s,p1p2so,m1s,m2s,m3s)
   C12o=C0i(cc12,p1s,p2s,p1p2so,m1s,m2s,m3s)
   C22o=C0i(cc22,p1s,p2s,p1p2so,m1s,m2s,m3s)

!   Sum1=Sum1-(C22o+2.D0*C12o+C11o)+2.D0*C00o
!   write(*,*) C22o+2.D0*C12o+C11o, 2.D0*C00o

!   Sum3=Sum3-2.D0*C2o-2.D0*C22o-C12o

   A(1) = A(1) &
        +x**(-1)*t*C22o*(1.D0/4.D0)+x**(-1)*t*C12o*(1.D0/2.D0)+x**(-1)*t*C11o*(1.D0/4.D0)+x**(-1)*t*C00o*(1.D0/2.D0)&
        +x**(-1)*t*C22*(-1.D0/4.D0)+x**(-1)*t*C12*(-1.D0/2.D0)+x**(-1)*t*C11*(-1.D0/4.D0)+x**(-1)*t*C00*(-1.D0/2.D0)&
        +C22o*(-1.D0)+C12o*(-2.D0)+C11o*(-1.D0)+C00o*(2.D0)+C2*(-2.D0)+C1*(-2.D0)+C22*(-1.D0)+C12*(-2.D0)+C11*(-1.D0)&
        +C00*(-2.D0)+t*C2*(-1.D0/2.D0)+t*C22*(-1.D0/2.D0)+t*C12*(-1.D0/2.D0)+x*C2*(-2.D0)+x*C22*(-2.D0)+x*C12*(-2.D0)

   A(2) = A(2) &
        +x**(-1)*C22o*(-1.D0)+x**(-1)*C12o*(-2.D0)+x**(-1)*C11o*(-1.D0)+x**(-1)*C22*(1.D0)+x**(-1)*C12*(2.D0)+x**(-1)*C11*(1.D0)

   A(3) = A(3) &
        +x**(-1)*t*C22o*(1.D0/4.D0)+x**(-1)*t*C12o*(1.D0/2.D0)+x**(-1)*t*C11o*(1.D0/4.D0)+x**(-1)*t*C00o*(1.D0/2.D0)&
        +x**(-1)*t*C22*(-1.D0/4.D0)+x**(-1)*t*C12*(-1.D0/2.D0)+x**(-1)*t*C11*(-1.D0/4.D0)+x**(-1)*t*C00*(-1.D0/2.D0)&
        +C22o*(-1.D0)+C12o*(-2.D0)+C11o*(-1.D0)+C00o*(2.D0)+C2*(-2.D0)+C1*(-2.D0)+C22*(-1.D0)+C12*(-2.D0)+C11*(-1.D0)&
        +C00*(-2.D0)+t*C2*(-1.D0/2.D0)+t*C22*(-1.D0/2.D0)+t*C12*(-1.D0/2.D0)+x*C2*(-2.D0)+x*C22*(-2.D0)+x*C12*(-2.D0)

   A(4) = A(4) &
        +x**(-1)*C22o*(1.D0/2.D0)+x**(-1)*C12o*(1.D0)+x**(-1)*C11o*(1.D0/2.D0)+x**(-1)*C00o*(1.D0)+x**(-1)*C22*(-1.D0/2.D0)&
        +x**(-1)*C12*(-1.D0)+x**(-1)*C11*(-1.D0/2.D0)+x**(-1)*C00*(-1.D0)+C2*(-1.D0)+C22*(-1.D0)+C12*(-1.D0)

   A(5) = A(5) &
        +x**(-1)*C22o*(-1.D0/2.D0)+x**(-1)*C12o*(-1.D0)+x**(-1)*C11o*(-1.D0/2.D0)+x**(-1)*C00o*(-1.D0)+x**(-1)*C22*(1.D0/2.D0)&
        +x**(-1)*C12*(1.D0)+x**(-1)*C11*(1.D0/2.D0)+x**(-1)*C00*(1.D0)+C2*(1.D0/2.D0)+C22*(1.D0/2.D0)+C12*(1.D0/2.D0)

   A(6) = A(6) &
        +x**(-1)*C00o*(1.D0)+x**(-1)*C00*(-1.D0)+C2*(-1.D0/2.D0)+C22*(-1.D0/2.D0)+C12*(-1.D0/2.D0)

   A(7) = A(7)

   A(8) = A(8) &
        +x**(-1)*C22o*(-1.D0/2.D0)+x**(-1)*C12o*(-1.D0)+x**(-1)*C11o*(-1.D0/2.D0)+x**(-1)*C00o*(-1.D0)+x**(-1)*C22*(1.D0/2.D0)&
        +x**(-1)*C12*(1.D0)+x**(-1)*C11*(1.D0/2.D0)+x**(-1)*C00*(1.D0)+C2*(1.D0)+C22*(1.D0)+C12*(1.D0)

!   write(*,*) 'Managed to calculate graph 3'
   IF(ipimode.GE.2) A=Asv+wn*(A-Asv) ! [v2] graph 3: pi0 p channel
    !* Graph 4:
    !*  reducible graph with nucleon self energy, renormalized
    ps=s
    pso=nuclm**2
    m1s=mpi2
    m2s=nuclm**2

    Bf0=B0(ps,m1s,m2s)
    Bf1=B1(ps,m1s,m2s)

    Bf0o=B0(pso,m1s,m2s)
    Bf1o=B1(pso,m1s,m2s)

    DBf1o=DB1(pso,m1s,m2s)

    IF(ipimode.GE.1) THEN     ! [v2] isospin weights: (2*(pi+ n) + 1*(pi0 p))/3, coefficients below carry the 3
       m2s=mN2c                ! pi+ n part with the internal neutron mass
       Bc0=B0(ps,m1s,m2s); Bc1=B1(ps,m1s,m2s)
       Bc0o=B0(pso,m1s,m2s); Bc1o=B1(pso,m1s,m2s); DBc1o=DB1(pso,m1s,m2s)
       m2s=nuclm**2
       m1s=mpin2
       Bf0=(2.D0*wc*Bc0+wn*B0(ps,m1s,m2s))/3.D0
       Bf1=(2.D0*wc*Bc1+wn*B1(ps,m1s,m2s))/3.D0
       Bf0o=(2.D0*wc*Bc0o+wn*B0(pso,m1s,m2s))/3.D0
       Bf1o=(2.D0*wc*Bc1o+wn*B1(pso,m1s,m2s))/3.D0
       DBf1o=(2.D0*wc*DBc1o+wn*DB1(pso,m1s,m2s))/3.D0
    ENDIF

!    write(*,*) 3.D0*Bf0o, 3.D0*Bf1o, 6.D0* DBf1o

!    Sum2=3.D0*(Bf0o+Bf1o+2.D0*DBf1o)

!    write(*,*) Sum1,Sum2,Sum2/Sum1
!    stop

    A(1) = A(1) &
         +x**(-2)*t*Bf1o*(-3.D0/8.D0)+x**(-2)*t*Bf1*(3.D0/8.D0)+x**(-1)*t*DBf1o*(-3.D0/4.D0)+x**(-1)*t*Bf1o*(-3.D0/8.D0)&
         +x**(-1)*t*Bf1*(3.D0/8.D0)+x**(-1)*t*Bf0o*(-3.D0/8.D0)+x**(-1)*t*Bf0*(3.D0/8.D0)+DBf1o*(-3.D0)+Bf1o*(-3.D0/2.D0)&
         +Bf1*(3.D0/2.D0)+Bf0o*(-3.D0/2.D0)+Bf0*(3.D0/2.D0)

    A(2) = A(2)

    A(3) = A(3) &
         +x**(-2)*t*Bf1o*(-3.D0/8.D0)+x**(-2)*t*Bf1*(3.D0/8.D0)+x**(-1)*t*DBf1o*(-3.D0/4.D0)+x**(-1)*t*Bf1o*(-3.D0/8.D0)&
         +x**(-1)*t*Bf1*(3.D0/8.D0)+x**(-1)*t*Bf0o*(-3.D0/8.D0)+x**(-1)*t*Bf0*(3.D0/8.D0)+DBf1o*(-3.D0)+Bf1o*(-3.D0/2.D0)&
         +Bf1*(3.D0/2.D0)+Bf0o*(-3.D0/2.D0)+Bf0*(3.D0/2.D0)

    A(4) = A(4) &
         +x**(-2)*Bf1o*(-3.D0/4.D0)+x**(-2)*Bf1*(3.D0/4.D0)+x**(-1)*DBf1o*(-3.D0/2.D0)+x**(-1)*Bf1o*(-3.D0/4.D0)&
         +x**(-1)*Bf1*(3.D0/4.D0)+x**(-1)*Bf0o*(-3.D0/4.D0)+x**(-1)*Bf0*(3.D0/4.D0)

    A(5) = A(5) &
         +x**(-2)*Bf1o*(3.D0/4.D0)+x**(-2)*Bf1*(-3.D0/4.D0)+x**(-1)*DBf1o*(3.D0/2.D0)+x**(-1)*Bf1o*(3.D0/4.D0)&
         +x**(-1)*Bf1*(-3.D0/4.D0)+x**(-1)*Bf0o*(3.D0/4.D0)+x**(-1)*Bf0*(-3.D0/4.D0)

    A(6) = A(6) &
         +x**(-2)*Bf1o*(-3.D0/4.D0)+x**(-2)*Bf1*(3.D0/4.D0)+x**(-1)*DBf1o*(-3.D0/2.D0)+x**(-1)*Bf1o*(-3.D0/4.D0)&
         +x**(-1)*Bf1*(3.D0/4.D0)+x**(-1)*Bf0o*(-3.D0/4.D0)+x**(-1)*Bf0*(3.D0/4.D0)

    A(7) = A(7)

    A(8) = A(8) &
         +x**(-2)*Bf1o*(3.D0/4.D0)+x**(-2)*Bf1*(-3.D0/4.D0)+x**(-1)*DBf1o*(3.D0/2.D0)+x**(-1)*Bf1o*(3.D0/4.D0)&
         +x**(-1)*Bf1*(-3.D0/4.D0)+x**(-1)*Bf0o*(3.D0/4.D0)+x**(-1)*Bf0*(-3.D0/4.D0)

!    write(*,*) 'Managed to calculate graph 4'
   Asv=A
   !* Graph 1:
   !*  Box graph with photons attached to the nucleons
   p1s=nuclm**2
   p2s=0.d0
   p3s=0.d0
   p4s=nuclm**2
   p2p3s=t
   p1p2s=s

   m1s=mpin2                 ! [v2] pi0 p channel (photons on the nucleons)
   m2s=nuclm**2
   m3s=m2s
   m4s=m2s

   D00=D0i(dd00,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11=D0i(dd11,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12=D0i(dd12,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D13=D0i(dd13,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22=D0i(dd22,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D111=D0i(dd111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D123=D0i(dd123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D112=D0i(dd112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D122=D0i(dd122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D113=D0i(dd113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D222=D0i(dd222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D001=D0i(dd001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D002=D0i(dd002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D003=D0i(dd003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)

   A(1) = A(1) &
        +D22*(2.D0)+D13*(4.D0)+D12*(8.D0)+D11*(4.D0)+D00*(2.D0)+D222*(1.D0)+D002*(2.D0)+D123*(6.D0)+D122*(6.D0)&
        +D113*(6.D0)+D112*(6.D0)+D111*(2.D0)+D001*(4.D0)+t*D22*(1.D0/4.D0)+t*D13*(-1.D0/2.D0)+t*D11*(-1.D0/2.D0)&
        +t*D00*(-1.D0/2.D0)+t*D222*(1.D0/4.D0)+t*D002*(1.D0/2.D0)+t*D123*(-1.D0/2.D0)+t*D122*(1.D0)+t*D113*(-2.D0)&
        +t*D112*(1.D0/2.D0)+t**2*D13*(1.D0/4.D0)+t**2*D123*(-1.D0/4.D0)+x*D22*(5.D0)+x*D13*(2.D0)+x*D12*(12.D0)&
        +x*D11*(2.D0)+x*D00*(-2.D0)+x*D222*(3.D0)+x*D002*(2.D0)+x*D123*(6.D0)+x*D122*(12.D0)+x*D112*(6.D0)+x*t*D22*(1.D0/2.D0)&
        +x*t*D13*(1.D0)+x*t*D222*(1.D0/2.D0)+x*t*D123*(-1.D0)+x*t*D122*(1.D0)+x**2*D22*(2.D0)+x**2*D222*(2.D0)+x**2*D122*(4.D0)

   A(2) = A(2) &
        +D22*(2.D0)+D12*(4.D0)+D222*(1.D0)+D123*(2.D0)+D122*(2.D0)+D113*(2.D0)+D112*(-2.D0)+D111*(-2.D0)+t*D22*(1.D0/4.D0)&
        +t*D12*(-1.D0/2.D0)+t*D222*(1.D0/4.D0)+t*D123*(1.D0/2.D0)+t*D112*(-1.D0/2.D0)+x*D22*(1.D0)+x*D12*(-2.D0)+x*D222*(1.D0)&
        +x*D123*(2.D0)+x*D112*(-2.D0)

   A(3) = A(3) &
        +D22*(2.D0)+D13*(4.D0)+D12*(8.D0)+D11*(4.D0)+D00*(2.D0)+D222*(1.D0)+D002*(6.D0)+D123*(6.D0)+D122*(6.D0)+D113*(6.D0)&
        +D112*(6.D0)+D111*(2.D0)+D001*(12.D0)+t*D22*(1.D0/4.D0)+t*D13*(-1.D0/2.D0)+t*D11*(-1.D0/2.D0)+t*D00*(1.D0/2.D0)&
        +t*D222*(1.D0/4.D0)+t*D002*(3.D0/2.D0)+t*D123*(-1.D0/2.D0)+t*D122*(1.D0)+t*D113*(-2.D0)+t*D112*(1.D0/2.D0)+t**2*D13*(1.D0/4.D0)&
        +t**2*D123*(-1.D0/4.D0)+x*D22*(5.D0)+x*D13*(2.D0)+x*D12*(12.D0)+x*D11*(2.D0)+x*D00*(2.D0)+x*D222*(3.D0)+x*D002*(6.D0)&
        +x*D123*(6.D0)+x*D122*(12.D0)+x*D112*(6.D0)+x*t*D22*(1.D0/2.D0)+x*t*D13*(1.D0)+x*t*D222*(1.D0/2.D0)+x*t*D123*(-1.D0)&
        +x*t*D122*(1.D0)+x**2*D22*(2.D0)+x**2*D222*(2.D0)+x**2*D122*(4.D0)

   A(4) = A(4) &
        +D22*(1.D0/2.D0)+D13*(-1.D0)+D11*(-1.D0)+D00*(-1.D0)+D222*(1.D0/2.D0)+D002*(1.D0)+D123*(1.D0)+D122*(2.D0)+D112*(1.D0)&
        +t*D13*(1.D0/2.D0)+t*D123*(-1.D0/2.D0)+x*D22*(1.D0)+x*D222*(1.D0)+x*D122*(2.D0)

   A(5) = A(5) &
        +D22*(1.D0/2.D0)+D12*(1.D0)+D002*(-1.D0)+D123*(-1.D0)+D122*(-1.D0/2.D0)+D113*(-3.D0/2.D0)+D112*(-1.D0)+D111*(-1.D0/2.D0)&
        +D001*(-3.D0)+t*D13*(-1.D0/2.D0)+t*D113*(1.D0/2.D0)+x*D13*(-2.D0)+x*D12*(-1.D0)+x*D123*(-1.D0)+x*D122*(-1.D0)+x*D112*(-1.D0)

   A(6) = A(6) &
        +D22*(1.D0/2.D0)+D12*(1.D0)+D002*(1.D0)+D123*(1.D0)+D122*(1.D0/2.D0)+D113*(3.D0/2.D0)+D112*(1.D0)+D111*(1.D0/2.D0)&
        +D001*(3.D0)+t*D13*(1.D0/2.D0)+t*D113*(-1.D0/2.D0)+x*D13*(2.D0)+x*D12*(1.D0)+x*D123*(1.D0)+x*D122*(1.D0)+x*D112*(1.D0)

   A(7) = A(7) &
        +D22*(-1.D0/2.D0)+D12*(1.D0)+D222*(-1.D0/2.D0)+D123*(-1.D0)+D112*(1.D0)

   A(8) = A(8) &
        +D22*(-1.D0/2.D0)+D13*(1.D0)+D11*(1.D0)+D00*(-1.D0)+D222*(-1.D0/2.D0)+D002*(-3.D0)+D123*(-1.D0)+D122*(-2.D0)+D112*(-1.D0)&
        +t*D13*(-1.D0/2.D0)+t*D123*(1.D0/2.D0)+x*D22*(-1.D0)+x*D222*(-1.D0)+x*D122*(-2.D0)


   IF(ipimode.GE.2) A=Asv+wn*(A-Asv) ! [v2] nucleon box: pi0 p channel
   Asv=A
! two-photon vertex
   p1s=t
   p2s=nuclm**2
   p1p2s=nuclm**2
   m1s=mpi2
   m2s=m1s
   m3s=mN2c                  ! [v2] pi+ n channel

   C2=C0i(cc2,p1s,p2s,p1p2s,m1s,m2s,m3s)

   A(1) = A(1) + 4.D0*C2/2.D0

! two-photon vertex, isoscalar graph

   ps=t
   m1s=mpi2
   m2s=m1s

   Bf0=B0(ps,m1s,m2s)

   A(1) = A(1) - 2.D0*Bf0/2.D0

! isoscalar three-pion loop

   p1s=0.D0
   p2s=t
   p1p2s=0.D0
   m1s=mpi2
   m2s=m1s
   m3s=m1s

   C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)

   A(1) = A(1) + 4.D0*C00

   A(2) = A(2) - 4.D0*C12   
   IF(ipimode.GE.2) A=Asv+wc*(A-Asv) ! [v2] two-photon / isoscalar pi-pi graphs: charged pions

   A=gafac*A

END SUBROUTINE CovOneLoop



!
!  Inclusion of the pi-Delta loops at order p^4/(M_D-M_N)
!  the external parameter Delta the subtraction constant in dimreg
!  Delta=0 corresponds with the MSbar scheme
!

SUBROUTINE ZeroCovOneDeltaLoop(s,t,mpi2,gafac,A,Delta)
  USE trigd
  IMPLICIT NONE
   REAL*8, INTENT(in) :: s,t,mpi2,gafac,Delta
   COMPLEX*16,INTENT(out) :: A(8)
  REAL*8, PARAMETER :: SGN(8)=(/1.D0,1.D0,-1.D0,-1.D0,-1.D0,-1.D0,-1.D0,1.D0/)
   COMPLEX*16  A0(8)
  REAL*8 s0,t0,u0,nuclm,theta,omega,omegap

  nuclm=1.d0
  theta=3.D0
  omega=1.d0
  omegap=omega/(1.D0+omega*(1.D0-DCOSd(theta))/nuclm)
  s0=2.d0*omega*nuclm+nuclm**2
  t0=-2.D0*omega*omegap*(1.D0-DCOSd(theta))
  u0=2.d0*nuclm**2-s0-t0
  call CovOneDeltaLoop(s0,t0,mpi2,gafac,A0,Delta)
  A=A0
  call CovOneDeltaLoop(u0,t0,mpi2,gafac,A0,Delta)
  A=A+SGN*A0
  A(3)=1.d0/(omega+omegap)*A(3)
END SUBROUTINE ZeroCovOneDeltaLoop


SUBROUTINE CovOneDeltaLoop(s,t,mpi2,gafac,A,Delta)
  USE trigd
!!   IMPLICIT COMPLEX*16 B, C, D
  implicit COMPLEX*16 (B-D)
#include "looptools.h"
   REAL*8, INTENT(in) :: s,t,mpi2,gafac,Delta
   COMPLEX*16,INTENT(out) :: A(8)
   COMPLEX*16, PARAMETER :: CI=(0.d0,1.d0)
   REAL*8 x,nuclm,deltam, MD
   REAL*8 DeltaAux, sw1
   REAL*8 m1s,m2s,m3s,m4s
   REAL*8 p1s,p2s,p2so,p3s,p4s,p1p2s,p1p2so,p2p3s,ps,pso
   INTEGER i,j
   COMPLEX*16 Sum1, Sum2, Sum3, AA(8)

   A=0.d0
   AA=0.d0

   Sum1=0.d0
   Sum2=0.d0
   Sum3=0.d0

   nuclm=1.d0   !* it's assumed that all is given in nucleon-mass units
   deltam=1.31d0
   MD=deltam
   x=(s-1.d0)/2.d0

   DeltaAux=1.d0

   !* Graph 1:
   !*  Box graph with photons attached to the pions
   !* Isospin factor = 4/3

   p1s=nuclm**2
   p2s=0.0d0
   p3s=0.0d0
   p4s=nuclm**2
   p2p3s=t
   p1p2s=s

   m1s=deltam**2
   m4s=mpi2
   m2s=mpi2
   m3s=mpi2

   CALL Setdelta(Delta)

   D00=D0i(dd00,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11=D0i(dd11,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12=D0i(dd12,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D13=D0i(dd13,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22=D0i(dd22,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D23=D0i(dd23,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D33=D0i(dd33,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D111=D0i(dd111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D112=D0i(dd112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D113=D0i(dd113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D122=D0i(dd122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D123=D0i(dd123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D133=D0i(dd133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D222=D0i(dd222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D223=D0i(dd223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D233=D0i(dd233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D333=D0i(dd333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D001=D0i(dd001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D002=D0i(dd002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D003=D0i(dd003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0000=D0i(dd0000,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0011=D0i(dd0011,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0012=D0i(dd0012,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0013=D0i(dd0013,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0022=D0i(dd0022,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0023=D0i(dd0023,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0033=D0i(dd0033,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1111=D0i(dd1111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1112=D0i(dd1112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1113=D0i(dd1113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1122=D0i(dd1122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1123=D0i(dd1123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1133=D0i(dd1133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1222=D0i(dd1222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1223=D0i(dd1223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1233=D0i(dd1233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1333=D0i(dd1333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D2222=D0i(dd2222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D2223=D0i(dd2223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D2233=D0i(dd2233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D2333=D0i(dd2333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D3333=D0i(dd3333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00001=D0i(dd00001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00002=D0i(dd00002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00003=D0i(dd00003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00111=D0i(dd00111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00112=D0i(dd00112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00113=D0i(dd00113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00122=D0i(dd00122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00123=D0i(dd00123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00133=D0i(dd00133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00222=D0i(dd00222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00223=D0i(dd00223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00233=D0i(dd00233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00333=D0i(dd00333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11111=D0i(dd11111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11112=D0i(dd11112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11113=D0i(dd11113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11122=D0i(dd11122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11123=D0i(dd11123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11133=D0i(dd11133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11222=D0i(dd11222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11223=D0i(dd11223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11233=D0i(dd11233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11333=D0i(dd11333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12222=D0i(dd12222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12223=D0i(dd12223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12233=D0i(dd12233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12333=D0i(dd12333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D13333=D0i(dd13333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22222=D0i(dd22222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22223=D0i(dd22223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22233=D0i(dd22233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22333=D0i(dd22333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D23333=D0i(dd23333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D33333=D0i(dd33333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)


   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   D0000D=D0i(dd0000,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00001D=D0i(dd00001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00002D=D0i(dd00002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00003D=D0i(dd00003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)

   DED0000=2.D0*(D0000-D0000D)/sw1
   DED00001=2.D0*(D00001-D00001D)/sw1
   DED00002=2.D0*(D00002-D00002D)/sw1
   DED00003=2.D0*(D00003-D00003D)/sw1

   A(1) = A(1) &
	+DED00002*(44.D0/27.D0)+DED00001*(88.D0/27.D0)+D00002*(-40.D0/9.D0)+D00001*(-80.D0/9.D0)+MD*DED0000*(-44.D0/27.D0)&
	+MD*D0000*(40.D0/9.D0)+t*DED00002*(-13.D0/27.D0)+t*DED00001*(-16.D0/9.D0)+t*D00002*(14.D0/9.D0)+t*D00001*(16.D0/3.D0)&
	+t*MD*DED0000*(20.D0/27.D0)+t*MD*D0022*(-4.D0/9.D0)+t*MD*D0013*(-8.D0/9.D0)+t*MD*D0012*(-8.D0/9.D0)+t*MD*D0000*(-28.D0/9.D0)&
	+t**2*DED00002*(-4.D0/27.D0)+t**2*D00002*(11.D0/9.D0)+t**2*MD*D0013*(2.D0/9.D0)+t**2*MD*D0012*(-1.D0/9.D0)+x*DED00002*(44.D0/27.D0)&
	+x*D00002*(-40.D0/9.D0)+x*t*DED00002*(-16.D0/27.D0)+x*t*D00002*(44.D0/9.D0)+x*t*MD*D0022*(-10.D0/9.D0)&
	+x*t*MD*D0012*(-4.D0/3.D0)+x**2*MD*D0022*(-8.D0/9.D0)

   A(2) = A(2) &
	+DED00002*(2.D0/27.D0)+DED00001*(4.D0/27.D0)+D00222*(14.D0/9.D0)+D0022*(20.D0/9.D0)+D002*(2.D0/3.D0)&
	+D00002*(-4.D0/9.D0)+D00123*(28.D0/9.D0)+D00122*(28.D0/9.D0)+D0012*(40.D0/9.D0)+D00113*(28.D0/9.D0)&
	+D00112*(-28.D0/9.D0)+D00111*(-28.D0/9.D0)+D001*(4.D0/3.D0)+D00001*(-8.D0/9.D0)+MD*DED0000*(2.D0/27.D0)&
	+MD*D0022*(-14.D0/9.D0)+MD*D002*(-20.D0/9.D0)+MD*D0013*(-28.D0/9.D0)+MD*D0011*(28.D0/9.D0)&
	+MD*D00*(-2.D0/3.D0)+MD*D0000*(8.D0/9.D0)+t*DED00002*(1.D0/54.D0)+t*D22222*(-1.D0/9.D0)+t*D2222*(-2.D0/9.D0)&
	+t*D00222*(-1.D0/2.D0)+t*D222*(-1.D0/9.D0)+t*D0022*(-7.D0/9.D0)+t*D002*(-5.D0/18.D0)+t*D00002*(-1.D0/9.D0)&
	+t*D00123*(-5.D0/9.D0)+t*D123*(-4.D0/9.D0)+t*D1223*(-8.D0/9.D0)+t*D12223*(-2.D0/3.D0)+t*D12222*(-4.D0/9.D0)&
	+t*D1222*(-8.D0/9.D0)+t*D00122*(-16.D0/9.D0)+t*D122*(-4.D0/9.D0)+t*D0012*(-8.D0/3.D0)+t*D00113*(-4.D0/3.D0)&
	+t*D113*(-4.D0/9.D0)+t*D11233*(-2.D0/3.D0)+t*D1123*(-8.D0/9.D0)+t*D11223*(-8.D0/9.D0)+t*D1122*(-4.D0/9.D0)&
	+t*D00112*(5.D0/9.D0)+t*D112*(-2.D0/9.D0)+t*D11133*(-4.D0/9.D0)+t*D11123*(4.D0/9.D0)+t*D11122*(4.D0/9.D0)&
	+t*D11113*(4.D0/9.D0)+t*D11112*(2.D0/9.D0)+t*D00111*(4.D0/3.D0)+t*D001*(-8.D0/9.D0)+t*MD*D2222*(1.D0/9.D0)&
	+t*MD*D222*(2.D0/9.D0)+t*MD*D0022*(10.D0/9.D0)+t*MD*D22*(1.D0/9.D0)+t*MD*D002*(14.D0/9.D0)+t*MD*D0013*(4.D0/3.D0)&
	+t*MD*D13*(2.D0/9.D0)+t*MD*D123*(4.D0/9.D0)+t*MD*D1223*(4.D0/9.D0)+t*MD*D1222*(2.D0/9.D0)+t*MD*D122*(4.D0/9.D0)&
	+t*MD*D12*(2.D0/9.D0)+t*MD*D1133*(4.D0/9.D0)+t*MD*D1123*(2.D0/9.D0)+t*MD*D1122*(-2.D0/9.D0)+t*MD*D1113*(-4.D0/9.D0)&
	+t*MD*D1112*(-2.D0/9.D0)+t*MD*D0011*(-4.D0/3.D0)+t*MD*D00*(4.D0/9.D0)+t**2*D22222*(-1.D0/24.D0)+t**2*D2222*(-1.D0/12.D0)&
	+t**2*D00222*(-7.D0/18.D0)+t**2*D222*(-1.D0/24.D0)+t**2*D0022*(-11.D0/18.D0)+t**2*D002*(-2.D0/9.D0)&
	+t**2*D00123*(-11.D0/18.D0)+t**2*D123*(-1.D0/36.D0)+t**2*D1223*(-1.D0/18.D0)+t**2*D12223*(-1.D0/9.D0)&
	+t**2*D12222*(-1.D0/9.D0)+t**2*D1222*(-2.D0/9.D0)+t**2*D122*(-1.D0/9.D0)+t**2*D113*(1.D0/9.D0)+t**2*D11233*(-1.D0/36.D0)&
	+t**2*D1123*(2.D0/9.D0)+t**2*D11222*(1.D0/18.D0)+t**2*D1122*(-1.D0/18.D0)+t**2*D00112*(11.D0/18.D0)+t**2*D112*(-1.D0/36.D0)&
	+t**2*D11133*(1.D0/9.D0)+t**2*D11122*(1.D0/9.D0)+t**2*D11113*(-1.D0/9.D0)+t**2*D11112*(1.D0/36.D0)+t**2*MD*D13*(-1.D0/18.D0)&
	+t**2*MD*D123*(-1.D0/9.D0)+t**2*MD*D1223*(-1.D0/18.D0)+t**2*MD*D1222*(1.D0/36.D0)+t**2*MD*D122*(1.D0/18.D0)&
	+t**2*MD*D12*(1.D0/36.D0)+t**2*MD*D1133*(-1.D0/9.D0)+t**2*MD*D1123*(1.D0/36.D0)+t**2*MD*D1113*(1.D0/9.D0)&
	+t**2*MD*D1112*(-1.D0/36.D0)+t**3*D123*(1.D0/36.D0)+t**3*D1223*(1.D0/18.D0)+t**3*D12223*(1.D0/36.D0)+t**3*D11233*(1.D0/18.D0)&
	+t**3*D11123*(-1.D0/18.D0)+x*DED00002*(2.D0/27.D0)+x*D00222*(14.D0/9.D0)+x*D0022*(20.D0/9.D0)+x*D002*(2.D0/3.D0)&
	+x*D00002*(-4.D0/9.D0)+x*D00123*(28.D0/9.D0)+x*D00112*(-28.D0/9.D0)+x*t*D22222*(-7.D0/18.D0)+x*t*D2222*(-7.D0/9.D0)&
	+x*t*D00222*(-14.D0/9.D0)+x*t*D222*(-7.D0/18.D0)+x*t*D0022*(-22.D0/9.D0)+x*t*D002*(-8.D0/9.D0)&
	+x*t*D00123*(-22.D0/9.D0)+x*t*D123*(-5.D0/9.D0)+x*t*D1223*(-10.D0/9.D0)+x*t*D12223*(-4.D0/3.D0)+x*t*D12222*(-10.D0/9.D0)&
	+x*t*D1222*(-20.D0/9.D0)+x*t*D122*(-10.D0/9.D0)+x*t*D11233*(-7.D0/9.D0)+x*t*D11223*(-10.D0/9.D0)&
	+x*t*D11222*(4.D0/9.D0)+x*t*D1122*(-2.D0/3.D0)+x*t*D00112*(22.D0/9.D0)+x*t*D112*(-1.D0/3.D0)+x*t*D11123*(4.D0/9.D0)&
	+x*t*D11122*(10.D0/9.D0)+x*t*D11112*(1.D0/3.D0)+x*t*MD*D2222*(5.D0/18.D0)+x*t*MD*D222*(5.D0/9.D0)+x*t*MD*D22*(5.D0/18.D0)&
	+x*t*MD*D1223*(5.D0/9.D0)+x*t*MD*D1222*(1.D0/3.D0)+x*t*MD*D122*(2.D0/3.D0)+x*t*MD*D12*(1.D0/3.D0)+x*t*MD*D1123*(1.D0/3.D0)&
	+x*t*MD*D1122*(-5.D0/9.D0)+x*t*MD*D1112*(-1.D0/3.D0)+x*t**2*D22222*(-1.D0/12.D0)+x*t**2*D2222*(-1.D0/6.D0)&
	+x*t**2*D222*(-1.D0/12.D0)+x*t**2*D123*(1.D0/9.D0)+x*t**2*D1223*(2.D0/9.D0)+x*t**2*D12223*(-1.D0/18.D0)&
	+x*t**2*D12222*(-1.D0/9.D0)+x*t**2*D1222*(-2.D0/9.D0)+x*t**2*D122*(-1.D0/9.D0)+x*t**2*D11233*(2.D0/9.D0)&
	+x*t**2*D11223*(-1.D0/9.D0)+x*t**2*D11222*(1.D0/6.D0)+x*t**2*D11123*(-2.D0/9.D0)+x*t**2*D11122*(1.D0/9.D0)&
	+x**2*D22222*(-2.D0/9.D0)+x**2*D2222*(-4.D0/9.D0)+x**2*D222*(-2.D0/9.D0)+x**2*D12223*(-4.D0/9.D0)&
	+x**2*D12222*(-4.D0/9.D0)+x**2*D1222*(-8.D0/9.D0)+x**2*D122*(-4.D0/9.D0)+x**2*D11223*(-4.D0/9.D0)&
	+x**2*D11222*(4.D0/9.D0)+x**2*D11122*(4.D0/9.D0)+x**2*MD*D2222*(2.D0/9.D0)+x**2*MD*D222*(4.D0/9.D0)+x**2*MD*D22*(2.D0/9.D0)&
	+x**2*MD*D1223*(4.D0/9.D0)+x**2*MD*D1122*(-4.D0/9.D0)+x**2*t*D22222*(-7.D0/18.D0)+x**2*t*D2222*(-7.D0/9.D0)&
	+x**2*t*D222*(-7.D0/18.D0)+x**2*t*D12223*(-7.D0/9.D0)+x**2*t*D12222*(-4.D0/9.D0)+x**2*t*D1222*(-8.D0/9.D0)&
	+x**2*t*D122*(-4.D0/9.D0)+x**2*t*D11223*(-4.D0/9.D0)+x**2*t*D11222*(7.D0/9.D0)+x**2*t*D11122*(4.D0/9.D0)&
	+x**3*D22222*(-2.D0/9.D0)+x**3*D2222*(-4.D0/9.D0)+x**3*D222*(-2.D0/9.D0)+x**3*D12223*(-4.D0/9.D0)+x**3*D11222*(4.D0/9.D0)

   A(3) = A(3)

   A(4) = A(4) &
	+DED00002*(2.D0/9.D0)+D00002*(-20.D0/3.D0)+MD*D0022*(4.D0/9.D0)+MD*D0012*(8.D0/9.D0)&
	+t*DED00002*(-8.D0/27.D0)+t*D00002*(22.D0/9.D0)+t*MD*D0012*(-2.D0/9.D0)+x*MD*D0022*(4.D0/9.D0)

   A(5) = A(5) &
	+DED00002*(1.D0/9.D0)+DED0000*(5.D0/27.D0)+D00222*(-2.D0/9.D0)+D0022*(-2.D0/9.D0)+D00002*(-10.D0/3.D0)+D0013*(-4.D0/9.D0)&
	+D00123*(-4.D0/9.D0)+D00122*(-8.D0/9.D0)+D0012*(-8.D0/9.D0)+D00112*(-4.D0/9.D0)+D0011*(-4.D0/9.D0)+D0000*(-22.D0/9.D0)&
	+MD*D0022*(2.D0/9.D0)+MD*D002*(2.D0/9.D0)+MD*D0012*(4.D0/9.D0)+MD*D001*(4.D0/9.D0)+t*DED00002*(-4.D0/27.D0)&
	+t*DED0000*(-4.D0/27.D0)+t*D0022*(1.D0/6.D0)+t*D00002*(11.D0/9.D0)+t*D0013*(5.D0/9.D0)+t*D00123*(2.D0/9.D0)+t*D00122*(2.D0/9.D0)&
	+t*D0012*(4.D0/9.D0)+t*D00112*(2.D0/9.D0)+t*D0011*(1.D0/9.D0)+t*D0000*(8.D0/9.D0)+t*MD*D0012*(-1.D0/9.D0)+t*MD*D001*(-1.D0/9.D0)&
	+t**2*D0013*(-1.D0/9.D0)+t**2*D00123*(-1.D0/18.D0)+x*D00222*(-4.D0/9.D0)+x*D0022*(-4.D0/9.D0)+x*D00122*(-8.D0/9.D0)&
	+x*D0012*(-8.D0/9.D0)+x*MD*D0022*(2.D0/9.D0)+x*MD*D002*(2.D0/9.D0)+x*t*D0022*(1.D0/3.D0)&
	+x*t*D00122*(2.D0/9.D0)+x*t*D0012*(4.D0/9.D0)+x**2*D00222*(-2.D0/9.D0)+x**2*D0022*(2.D0/9.D0)

   A(6) = A(6) &
	+DED00002*(-1.D0/9.D0)+DED0000*(-5.D0/27.D0)+D00222*(2.D0/9.D0)+D0022*(2.D0/9.D0)+D00002*(10.D0/3.D0)+D0013*(4.D0/9.D0)&
	+D00123*(4.D0/9.D0)+D00122*(8.D0/9.D0)+D0012*(8.D0/9.D0)+D00112*(4.D0/9.D0)+D0011*(4.D0/9.D0)+D0000*(22.D0/9.D0)&
	+MD*D0022*(-2.D0/9.D0)+MD*D002*(-2.D0/9.D0)+MD*D0012*(-4.D0/9.D0)+MD*D001*(-4.D0/9.D0)+t*DED00002*(4.D0/27.D0)+t*DED0000*(4.D0/27.D0)&
	+t*D0022*(-1.D0/6.D0)+t*D00002*(-11.D0/9.D0)+t*D0013*(-5.D0/9.D0)+t*D00123*(-4.D0/9.D0)+t*D00122*(-2.D0/9.D0)+t*D0012*(-4.D0/9.D0)&
	+t*D0011*(-1.D0/9.D0)+t*D0000*(-8.D0/9.D0)+t*MD*D0012*(1.D0/9.D0)+t*MD*D001*(1.D0/9.D0)+t**2*D0013*(1.D0/9.D0)+t**2*D00123*(1.D0/18.D0)&
	+x*D00222*(4.D0/9.D0)+x*D0022*(4.D0/9.D0)+x*D00122*(8.D0/9.D0)+x*D0012*(8.D0/9.D0)+x*MD*D0022*(-2.D0/9.D0)+x*MD*D002*(-2.D0/9.D0)&
	+x*t*D0022*(-1.D0/3.D0)+x*t*D00122*(-2.D0/9.D0)+x*t*D0012*(-4.D0/9.D0)+x**2*D00222*(2.D0/9.D0)+x**2*D0022*(-2.D0/9.D0)

   A(7) = A(7) &
	+DED00002*(-1.D0/27.D0)+D22222*(-1.D0/9.D0)+D2222*(-2.D0/9.D0)+D00222*(-19.D0/9.D0)+D222*(-1.D0/9.D0)+D0022*(-10.D0/3.D0)&
	+D002*(-11.D0/9.D0)+D00002*(2.D0/9.D0)+D00123*(-38.D0/9.D0)+D123*(-2.D0/9.D0)+D1223*(-4.D0/9.D0)+D12223*(-4.D0/9.D0)&
	+D12222*(-4.D0/9.D0)+D1222*(-8.D0/9.D0)+D122*(-4.D0/9.D0)+D11233*(-2.D0/9.D0)+D11223*(-4.D0/9.D0)+D1122*(-4.D0/9.D0)&
	+D00112*(38.D0/9.D0)+D112*(-2.D0/9.D0)+D11122*(4.D0/9.D0)+D11112*(2.D0/9.D0)+MD*D2222*(1.D0/9.D0)+MD*D222*(2.D0/9.D0)&
	+MD*D22*(1.D0/9.D0)+MD*D1223*(2.D0/9.D0)+MD*D1222*(2.D0/9.D0)+MD*D122*(4.D0/9.D0)+MD*D12*(2.D0/9.D0)&
	+MD*D1123*(2.D0/9.D0)+MD*D1122*(-2.D0/9.D0)+MD*D1112*(-2.D0/9.D0)+t*D22222*(1.D0/12.D0)+t*D2222*(1.D0/6.D0)&
	+t*D00222*(7.D0/9.D0)+t*D222*(1.D0/12.D0)+t*D0022*(11.D0/9.D0)+t*D002*(4.D0/9.D0)+t*D00123*(11.D0/9.D0)+t*D123*(5.D0/18.D0)&
	+t*D1223*(5.D0/9.D0)+t*D12223*(4.D0/9.D0)+t*D12222*(2.D0/9.D0)+t*D1222*(4.D0/9.D0)+t*D122*(2.D0/9.D0)+t*D11233*(1.D0/2.D0)&
	+t*D11223*(2.D0/9.D0)+t*D11222*(-1.D0/9.D0)+t*D1122*(1.D0/9.D0)+t*D00112*(-11.D0/9.D0)+t*D112*(1.D0/18.D0)+t*D11123*(-4.D0/9.D0)&
	+t*D11122*(-2.D0/9.D0)+t*D11112*(-1.D0/18.D0)+t*MD*D1222*(-1.D0/18.D0)+t*MD*D122*(-1.D0/9.D0)+t*MD*D12*(-1.D0/18.D0)&
	+t*MD*D1123*(-1.D0/18.D0)+t*MD*D1112*(1.D0/18.D0)+t**2*D123*(-1.D0/18.D0)+t**2*D1223*(-1.D0/9.D0)+t**2*D12223*(-1.D0/18.D0)&
	+t**2*D11233*(-1.D0/9.D0)+t**2*D11123*(1.D0/9.D0)+x*D22222*(-2.D0/9.D0)+x*D2222*(-4.D0/9.D0)+x*D222*(-2.D0/9.D0)&
	+x*D12223*(-4.D0/9.D0)+x*D12222*(-4.D0/9.D0)+x*D1222*(-8.D0/9.D0)+x*D122*(-4.D0/9.D0)+x*D11223*(-4.D0/9.D0)+x*D11222*(4.D0/9.D0)&
	+x*D11122*(4.D0/9.D0)+x*MD*D2222*(1.D0/9.D0)+x*MD*D222*(2.D0/9.D0)+x*MD*D22*(1.D0/9.D0)+x*MD*D1223*(2.D0/9.D0)&
	+x*MD*D1122*(-2.D0/9.D0)+x*t*D22222*(1.D0/6.D0)+x*t*D2222*(1.D0/3.D0)+x*t*D222*(1.D0/6.D0)+x*t*D12223*(1.D0/3.D0)&
	+x*t*D12222*(2.D0/9.D0)+x*t*D1222*(4.D0/9.D0)+x*t*D122*(2.D0/9.D0)+x*t*D11223*(2.D0/9.D0)+x*t*D11222*(-1.D0/3.D0)&
	+x*t*D11122*(-2.D0/9.D0)+x**2*D22222*(1.D0/9.D0)+x**2*D2222*(2.D0/9.D0)+x**2*D222*(1.D0/9.D0)+x**2*D12223*(2.D0/9.D0)&
	+x**2*D11222*(-2.D0/9.D0)


   A(8) = A(8)

   !Sum1=A(1)

   !write(*,*) '1', Sum1/Sum1

   !* Graph 2:
   !*  Box graphs with one photon attached to the pion and the other attached to the piNDelta vertex
   !* Isospin factor = 4/3

   p1s=nuclm**2
   p2s=0.d0
   p1p2s=s

   m1s=deltam**2
   m2s=mpi2
   m3s=mpi2

   CALL Setdelta(Delta)

   C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C11=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C22=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C111=C0i(cc111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C112=C0i(cc112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C122=C0i(cc122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C222=C0i(cc222,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0000=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0011=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0012=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0022=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1111=C0i(cc1111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1112=C0i(cc1112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1122=C0i(cc1122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1222=C0i(cc1222,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C2222=C0i(cc2222,p1s,p2s,p1p2s,m1s,m2s,m3s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   C00D=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001D=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002D=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0000D=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0011D=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0012D=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0022D=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)

   DEC00=2.D0*(C00-C00D)/sw1
   DEC001=2.D0*(C001-C001D)/sw1
   DEC002=2.D0*(C002-C002D)/sw1
   DEC0000=2.D0*(C0000-C0000D)/sw1
   DEC0011=2.D0*(C0011-C0011D)/sw1
   DEC0012=2.D0*(C0012-C0012D)/sw1
   DEC0022=2.D0*(C0022-C0022D)/sw1

   A(1) = A(1) &
	+DEC0022*(4.D0/27.D0)+DEC0012*(8.D0/27.D0)+DEC0011*(4.D0/27.D0)+DEC0000*(4.D0/27.D0)+C0022*(-8.D0/9.D0)&
	+C0012*(-16.D0/9.D0)+C0011*(-8.D0/9.D0)+C0000*(-8.D0/9.D0)+MD*DEC002*(-4.D0/27.D0)+MD*DEC001*(-4.D0/27.D0)&
	+MD*C002*(8.D0/9.D0)+MD*C001*(8.D0/9.D0)+t*DEC0022*(1.D0/27.D0)+t*DEC0012*(1.D0/27.D0)+t*C0022*(-2.D0/9.D0)&
	+t*C0012*(-2.D0/9.D0)+t*MD*DEC002*(-1.D0/27.D0)+t*MD*C002*(-1.D0/9.D0)+x*DEC0022*(8.D0/27.D0)+x*DEC0012*(8.D0/27.D0)&
	+x*C0022*(-16.D0/9.D0)+x*C0012*(-16.D0/9.D0)+x*MD*DEC002*(-4.D0/27.D0)+x*MD*C002*(8.D0/9.D0)+x*t*DEC0022*(1.D0/27.D0)&
	+x*t*C0022*(-2.D0/9.D0)+x**2*DEC0022*(4.D0/27.D0)+x**2*C0022*(-8.D0/9.D0)

   A(2) = A(2) &
	+DEC0022*(-5.D0/9.D0)+DEC002*(-4.D0/9.D0)+DEC0012*(-26.D0/27.D0)+DEC0011*(-11.D0/27.D0)+DEC001*(-10.D0/27.D0)&
	+C2222*(4.D0/9.D0)+C222*(4.D0/9.D0)+C0022*(2.D0)+C002*(4.D0/3.D0)+C1222*(4.D0/3.D0)+C122*(8.D0/9.D0)+C0012*(28.D0/9.D0)&
	+C1122*(4.D0/3.D0)+C112*(4.D0/9.D0)+C1112*(4.D0/9.D0)+C0011*(10.D0/9.D0)+C001*(8.D0/9.D0)+MD*DEC002*(11.D0/27.D0)&
	+MD*DEC001*(11.D0/27.D0)+MD*DEC00*(10.D0/27.D0)+MD*C222*(-4.D0/9.D0)+MD*C22*(-4.D0/9.D0)+MD*C002*(-16.D0/9.D0)&
	+MD*C122*(-8.D0/9.D0)+MD*C12*(-4.D0/9.D0)+MD*C112*(-4.D0/9.D0)+MD*C001*(-10.D0/9.D0)+MD*C00*(-8.D0/9.D0)+t*DEC0022*(-1.D0/12.D0)&
	+t*DEC002*(-2.D0/27.D0)+t*DEC0012*(-7.D0/108.D0)+t*C2222*(5.D0/36.D0)+t*C222*(5.D0/36.D0)+t*C0022*(2.D0/3.D0)&
	+t*C002*(4.D0/9.D0)+t*C1222*(11.D0/36.D0)+t*C122*(1.D0/6.D0)+t*C0012*(5.D0/9.D0)+t*C1122*(7.D0/36.D0)+t*C112*(1.D0/36.D0)&
	+t*C1112*(1.D0/36.D0)+t*MD*C222*(1.D0/36.D0)+t*MD*C22*(1.D0/36.D0)+t*MD*C12*(-1.D0/36.D0)+t*MD*C112*(-1.D0/36.D0)&
	+x*DEC0022*(-1.D0/3.D0)+x*DEC002*(-8.D0/27.D0)+x*DEC0012*(-7.D0/27.D0)+x*C2222*(11.D0/9.D0)+x*C222*(11.D0/9.D0)&
	+x*C0022*(8.D0/3.D0)+x*C002*(16.D0/9.D0)+x*C1222*(25.D0/9.D0)+x*C122*(14.D0/9.D0)+x*C0012*(20.D0/9.D0)+x*C1122*(17.D0/9.D0)&
	+x*C112*(1.D0/3.D0)+x*C1112*(1.D0/3.D0)+x*MD*C222*(-7.D0/9.D0)+x*MD*C22*(-7.D0/9.D0)+x*MD*C122*(-10.D0/9.D0)&
	+x*MD*C12*(-1.D0/3.D0)+x*MD*C112*(-1.D0/3.D0)+x*t*C2222*(2.D0/9.D0)+x*t*C222*(2.D0/9.D0)+x*t*C1222*(1.D0/3.D0)&
	+x*t*C122*(1.D0/9.D0)+x*t*C1122*(1.D0/9.D0)+x**2*C2222*(8.D0/9.D0)+x**2*C222*(8.D0/9.D0)+x**2*C1222*(4.D0/3.D0)&
	+x**2*C122*(4.D0/9.D0)+x**2*C1122*(4.D0/9.D0)

   A(3) = A(3) &
	+DEC0022*(4.D0/27.D0)+DEC0012*(8.D0/27.D0)+DEC0011*(4.D0/27.D0)+DEC0000*(4.D0/9.D0)+C0022*(4.D0/9.D0)&
 	+C0012*(8.D0/9.D0)+C0011*(4.D0/9.D0)+C0000*(8.D0/3.D0)+MD*DEC002*(-4.D0/27.D0)+MD*DEC001*(-4.D0/27.D0)&
	+MD*C002*(-4.D0/9.D0)+MD*C001*(-4.D0/9.D0)+x*DEC0022*(8.D0/27.D0)+x*DEC0012*(8.D0/27.D0)+x*C0022*(8.D0/9.D0)&
	+x*C0012*(8.D0/9.D0)+x*MD*DEC002*(-4.D0/27.D0)+x*MD*C002*(-4.D0/9.D0)

   A(4) = A(4) &
	+DEC0022*(2.D0/27.D0)+DEC0012*(2.D0/27.D0)+C0022*(-4.D0/9.D0)+C0012*(-4.D0/9.D0)+MD*DEC002*(-2.D0/27.D0)&
	+MD*C002*(-2.D0/9.D0)+x*DEC0022*(2.D0/27.D0)+x*C0022*(-4.D0/9.D0)

   A(5) = A(5) &
	+DEC0022*(-1.D0/6.D0)+DEC002*(-5.D0/54.D0)+DEC0012*(-5.D0/27.D0)+DEC0011*(-1.D0/54.D0)+DEC001*(1.D0/54.D0)&
	+DEC0000*(1.D0/9.D0)+C2222*(-1.D0/9.D0)+C222*(-1.D0/9.D0)+C0022*(-4.D0/3.D0)+C002*(-7.D0/9.D0)+C1222*(-1.D0/3.D0)&
	+C122*(-2.D0/9.D0)+C0012*(-14.D0/9.D0)+C1122*(-1.D0/3.D0)+C112*(-1.D0/9.D0)+C1112*(-1.D0/9.D0)+C0011*(-2.D0/9.D0)&
	+C001*(-1.D0/9.D0)+C0000*(-4.D0/3.D0)+MD*DEC002*(1.D0/54.D0)+MD*DEC001*(1.D0/54.D0)+MD*DEC00*(-1.D0/54.D0)&
	+MD*C222*(1.D0/9.D0)+MD*C22*(1.D0/9.D0)+MD*C002*(2.D0/9.D0)+MD*C122*(2.D0/9.D0)+MD*C12*(1.D0/9.D0)&
	+MD*C112*(1.D0/9.D0)+MD*C001*(2.D0/9.D0)+MD*C00*(1.D0/9.D0)+x*DEC0022*(-5.D0/54.D0)+x*DEC002*(-1.D0/27.D0)&
	+x*DEC0012*(-7.D0/54.D0)+x*C2222*(-5.D0/18.D0)+x*C222*(-5.D0/18.D0)+x*C0022*(-10.D0/9.D0)+x*C002*(-4.D0/9.D0)&
	+x*C1222*(-11.D0/18.D0)+x*C122*(-1.D0/3.D0)+x*C0012*(-8.D0/9.D0)+x*C1122*(-7.D0/18.D0)+x*C112*(-1.D0/18.D0)&
	+x*C1112*(-1.D0/18.D0)+x*MD*C222*(1.D0/6.D0)+x*MD*C22*(1.D0/6.D0)+x*MD*C122*(2.D0/9.D0)+x*MD*C12*(1.D0/18.D0)&
	+x*MD*C112*(1.D0/18.D0)+x**2*C2222*(-1.D0/9.D0)+x**2*C222*(-1.D0/9.D0)+x**2*C1222*(-2.D0/9.D0)+x**2*C122*(-1.D0/9.D0)&
	+x**2*C1122*(-1.D0/9.D0)

   A(6) = A(6) &
	+DEC0022*(-1.D0/54.D0)+DEC002*(-1.D0/54.D0)+DEC0011*(1.D0/54.D0)+DEC001*(-1.D0/54.D0)+DEC0000*(-1.D0/9.D0)&
	+C0022*(4.D0/9.D0)+C002*(1.D0/9.D0)+C0012*(2.D0/3.D0)+C0011*(2.D0/9.D0)+C001*(1.D0/9.D0)+C0000*(4.D0/3.D0)&
	+MD*DEC002*(1.D0/54.D0)+MD*DEC001*(-1.D0/54.D0)+MD*DEC00*(1.D0/54.D0)+MD*C002*(-1.D0/9.D0)+MD*C001*(-2.D0/9.D0)&
	+MD*C00*(-1.D0/9.D0)+x*DEC0022*(5.D0/54.D0)+x*DEC002*(1.D0/27.D0)+x*DEC0012*(7.D0/54.D0)+x*C2222*(1.D0/18.D0)&
	+x*C222*(1.D0/18.D0)+x*C0022*(10.D0/9.D0)+x*C002*(4.D0/9.D0)+x*C1222*(1.D0/6.D0)+x*C122*(1.D0/9.D0)&
	+x*C0012*(8.D0/9.D0)+x*C1122*(1.D0/6.D0)+x*C112*(1.D0/18.D0)+x*C1112*(1.D0/18.D0)+x*MD*C222*(-1.D0/18.D0)&
	+x*MD*C22*(-1.D0/18.D0)+x*MD*C122*(-1.D0/9.D0)+x*MD*C12*(-1.D0/18.D0)+x*MD*C112*(-1.D0/18.D0)+x**2*C2222*(1.D0/9.D0)&
	+x**2*C222*(1.D0/9.D0)+x**2*C1222*(2.D0/9.D0)+x**2*C122*(1.D0/9.D0)+x**2*C1122*(1.D0/9.D0)

   A(7) = A(7) &
	+DEC0022*(1.D0/6.D0)+DEC002*(4.D0/27.D0)+DEC0012*(7.D0/54.D0)+C2222*(-5.D0/18.D0)+C222*(-5.D0/18.D0)&
	+C0022*(-4.D0/3.D0)+C002*(-8.D0/9.D0)+C1222*(-11.D0/18.D0)+C122*(-1.D0/3.D0)+C0012*(-10.D0/9.D0)&
	+C1122*(-7.D0/18.D0)+C112*(-1.D0/18.D0)+C1112*(-1.D0/18.D0)+MD*C222*(-1.D0/18.D0)+MD*C22*(-1.D0/18.D0)&
	+MD*C12*(1.D0/18.D0)+MD*C112*(1.D0/18.D0)+x*C2222*(-4.D0/9.D0)+x*C222*(-4.D0/9.D0)+x*C1222*(-2.D0/3.D0)&
	+x*C122*(-2.D0/9.D0)+x*C1122*(-2.D0/9.D0)

   A(8) = A(8)

   !Sum2=(A(1)-Sum1)/Sum1
   !Sum3=A(1)

   !write(*,*) '2', Sum2

   !* Graph 3:
   !*  Box graph with photons attached to the piNDelta vertices
   !* Isospin factor = 4/3

   p1s=s

   m1s=deltam**2
   m2s=mpi2

   CALL Setdelta(Delta)

   Bf00=B0i(bb00,p1s,m1s,m2s)
   Bf11=B0i(bb11,p1s,m1s,m2s)
   Bf001=B0i(bb001,p1s,m1s,m2s)
   Bf111=B0i(bb111,p1s,m1s,m2s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   Bf00D=B0i(bb00,p1s,m1s,m2s)
   Bf001D=B0i(bb001,p1s,m1s,m2s)
   Bf11D=B0i(bb11,p1s,m1s,m2s)
   Bf111D=B0i(bb111,p1s,m1s,m2s)

   DEBf00=2.D0*(Bf00-Bf00D)/sw1
   DEBf001=2.D0*(Bf001-Bf001D)/sw1
   DEBf11=2.D0*(Bf11-Bf11D)/sw1
   DEBf111=2.D0*(Bf111-Bf111D)/sw1

   A(1) = A(1) &
	+DEBf111*(1.D0/27.D0)+DEBf001*(11.D0/27.D0)+Bf111*(-2.D0/9.D0)+Bf001*(-10.D0/9.D0)+MD*DEBf11*(-1.D0/27.D0)+MD*DEBf00*(-1.D0/3.D0)&
	+MD*Bf11*(2.D0/9.D0)+MD*Bf00*(2.D0/3.D0)+t*DEBf111*(1.D0/108.D0)+t*DEBf001*(11.D0/108.D0)+t*Bf111*(-1.D0/18.D0)+t*Bf001*(-5.D0/18.D0)&
	+x*DEBf111*(1.D0/9.D0)+x*DEBf001*(11.D0/27.D0)+x*Bf111*(-2.D0/3.D0)+x*Bf001*(-10.D0/9.D0)+x*MD*DEBf11*(-2.D0/27.D0)+x*MD*Bf11*(4.D0/9.D0)&
	+x*t*DEBf111*(1.D0/54.D0)+x*t*Bf111*(-1.D0/9.D0)+x**2*DEBf111*(2.D0/27.D0)+x**2*Bf111*(-4.D0/9.D0)

   A(2) = A(2) &
	+DEBf111*(1.D0/108.D0)+Bf111*(-1.D0/18.D0)+MD*DEBf11*(-1.D0/36.D0)+t*DEBf111*(1.D0/432.D0)+t*Bf111*(-1.D0/72.D0)&
	+x*DEBf111*(1.D0/108.D0)+x*Bf111*(-1.D0/18.D0)

   A(3) = A(3) &
	+DEBf111*(1.D0/27.D0)+DEBf001*(1.D0/9.D0)+Bf111*(1.D0/9.D0)+Bf001*(2.D0/3.D0)+MD*DEBf11*(-1.D0/27.D0)+MD*DEBf00*(1.D0/27.D0)&
	+MD*Bf11*(-1.D0/9.D0)+MD*Bf00*(-2.D0/9.D0)+t*DEBf111*(1.D0/108.D0)+t*DEBf001*(1.D0/36.D0)+t*Bf111*(1.D0/36.D0)+t*Bf001*(1.D0/6.D0)&
	+x*DEBf111*(1.D0/9.D0)+x*DEBf001*(1.D0/9.D0)+x*Bf111*(1.D0/3.D0)+x*Bf001*(2.D0/3.D0)+x*MD*DEBf11*(-2.D0/27.D0)+x*MD*Bf11*(-2.D0/9.D0)&
	+x*t*DEBf111*(1.D0/54.D0)+x*t*Bf111*(1.D0/18.D0)+x**2*DEBf111*(2.D0/27.D0)+x**2*Bf111*(2.D0/9.D0)

   A(4) = A(4) &
	+DEBf111*(1.D0/54.D0)+DEBf001*(11.D0/54.D0)+Bf111*(-1.D0/9.D0)+Bf001*(-5.D0/9.D0)+x*DEBf111*(1.D0/27.D0)+x*Bf111*(-2.D0/9.D0)

   A(5) = A(5) &
	+DEBf111*(-1.D0/108.D0)+DEBf001*(-1.D0/27.D0)+Bf111*(-1.D0/36.D0)+Bf001*(-1.D0/9.D0)+MD*DEBf11*(-1.D0/108.D0)+MD*Bf11*(-1.D0/36.D0)&
	+x*DEBf111*(-1.D0/54.D0)+x*Bf111*(-1.D0/18.D0)

   A(6) = A(6) &
	+DEBf111*(1.D0/108.D0)+DEBf001*(1.D0/27.D0)+Bf111*(1.D0/36.D0)+Bf001*(1.D0/9.D0)+MD*DEBf11*(-1.D0/108.D0)+MD*Bf11*(-1.D0/36.D0)&
	+x*DEBf111*(1.D0/54.D0)+x*Bf111*(1.D0/18.D0)

   A(7) = A(7) &
	+DEBf111*(-1.D0/216.D0)+Bf111*(1.D0/36.D0)

   A(8) = A(8) &
	+DEBf111*(-1.D0/54.D0)+DEBf001*(-1.D0/18.D0)+Bf111*(-1.D0/18.D0)+Bf001*(-1.D0/3.D0)+x*DEBf111*(-1.D0/27.D0)+x*Bf111*(-1.D0/9.D0)

!    AA=A*3.D0/4.D0
!    A=0.d0

   !Sum2=(A(1)-Sum3)/Sum1
   !Sum3=A(1)

   !write(*,*) '3', Sum2

   !* Graph 4:
   !*  Reducible graphs with gammaNN vertices
   !* Isospin factor = -2/3

   p1s=nuclm**2

   m1s=deltam**2
   m2s=mpi2

   CALL Setdelta(Delta)

   Bf00o=B0i(bb00,p1s,m1s,m2s)
   Bf11o=B0i(bb11,p1s,m1s,m2s)
   Bf001o=B0i(bb001,p1s,m1s,m2s)
   Bf111o=B0i(bb111,p1s,m1s,m2s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   Bf00Do=B0i(bb00,p1s,m1s,m2s)
   Bf001Do=B0i(bb001,p1s,m1s,m2s)
   Bf11Do=B0i(bb11,p1s,m1s,m2s)
   Bf111Do=B0i(bb111,p1s,m1s,m2s)

   DEBf00o=2.D0*(Bf00o-Bf00Do)/sw1
   DEBf001o=2.D0*(Bf001o-Bf001Do)/sw1
   DEBf11o=2.D0*(Bf11o-Bf11Do)/sw1
   DEBf111o=2.D0*(Bf111o-Bf111Do)/sw1

   p1s=nuclm**2+0.0001d0
   p2s=0.d0
   p1p2s=nuclm**2

   m1s=deltam**2
   m2s=mpi2
   m3s=mpi2

   CALL Setdelta(Delta)

   C00o=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C11o=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12o=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C22o=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001o=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002o=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C111o=C0i(cc111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C112o=C0i(cc112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C122o=C0i(cc122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C222o=C0i(cc222,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0000o=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0011o=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0012o=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0022o=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1111o=C0i(cc1111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1112o=C0i(cc1112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1122o=C0i(cc1122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1222o=C0i(cc1222,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C2222o=C0i(cc2222,p1s,p2s,p1p2s,m1s,m2s,m3s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   C00Do=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001Do=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002Do=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0000Do=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0011Do=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0012Do=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0022Do=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)

   DEC00o=2.D0*(C00o-C00Do)/sw1
   DEC001o=2.D0*(C001o-C001Do)/sw1
   DEC002o=2.D0*(C002o-C002Do)/sw1
   DEC0000o=2.D0*(C0000o-C0000Do)/sw1
   DEC0011o=2.D0*(C0011o-C0011Do)/sw1
   DEC0012o=2.D0*(C0012o-C0012Do)/sw1
   DEC0022o=2.D0*(C0022o-C0022Do)/sw1

   A(1) = A(1) &
	+x**(-1)*t*DEBf001o*(-1.D0/108.D0)+x**(-1)*t*DEBf001*(1.D0/108.D0)+x**(-1)*t*Bf001o*(1.D0/18.D0)+x**(-1)*t*Bf001*(-1.D0/18.D0)&
	+x**(-1)*t*DEC0022o*(1.D0/27.D0)+x**(-1)*t*DEC0012o*(2.D0/27.D0)+x**(-1)*t*DEC0011o*(1.D0/27.D0)+x**(-1)*t*DEC0000o*(-5.D0/54.D0)&
	+x**(-1)*t*DEC0022*(-1.D0/27.D0)+x**(-1)*t*DEC0012*(-2.D0/27.D0)+x**(-1)*t*DEC0011*(-1.D0/27.D0)+x**(-1)*t*DEC0000*(5.D0/54.D0)&
	+x**(-1)*t*C0022o*(1.D0/9.D0)+x**(-1)*t*C0012o*(2.D0/9.D0)+x**(-1)*t*C0011o*(1.D0/9.D0)+x**(-1)*t*C0000o*(11.D0/9.D0)&
	+x**(-1)*t*C0022*(-1.D0/9.D0)+x**(-1)*t*C0012*(-2.D0/9.D0)+x**(-1)*t*C0011*(-1.D0/9.D0)+x**(-1)*t*C0000*(-11.D0/9.D0)&
	+x**(-1)*t*MD*DEBf00o*(1.D0/108.D0)+x**(-1)*t*MD*DEBf00*(-1.D0/108.D0)+x**(-1)*t*MD*Bf00o*(-1.D0/18.D0)+x**(-1)*t*MD*Bf00*(1.D0/18.D0)&
	+x**(-1)*t*MD*DEC002o*(-1.D0/27.D0)+x**(-1)*t*MD*DEC001o*(-1.D0/27.D0)+x**(-1)*t*MD*DEC002*(1.D0/27.D0)+x**(-1)*t*MD*DEC001*(1.D0/27.D0)&
	+x**(-1)*t*MD*C002o*(-1.D0/9.D0)+x**(-1)*t*MD*C001o*(-1.D0/9.D0)+x**(-1)*t*MD*C002*(1.D0/9.D0)+x**(-1)*t*MD*C001*(1.D0/9.D0)&
	+DEBf001o*(19.D0/54.D0)+DEBf001*(1.D0/54.D0)+Bf001o*(-7.D0/9.D0)+Bf001*(-1.D0/9.D0)+DEC0022o*(-1.D0/3.D0)+DEC002o*(-1.D0/3.D0)&
	+DEC0012o*(-2.D0/3.D0)+DEC0011o*(-1.D0/3.D0)+DEC001o*(-1.D0/3.D0)+DEC0000o*(-1.D0/3.D0)+DEC0000*(11.D0/27.D0)+C0022o*(2.D0/3.D0)&
	+C002o*(2.D0/3.D0)+C0012o*(4.D0/3.D0)+C0011o*(2.D0/3.D0)+C001o*(2.D0/3.D0)+C0000o*(2.D0/3.D0)+C0000*(-10.D0/9.D0)+MD*DEBf00o*(-19.D0/54.D0)&
	+MD*DEBf00*(1.D0/54.D0)+MD*Bf00o*(7.D0/9.D0)+MD*Bf00*(-1.D0/9.D0)+MD*DEC002o*(1.D0/3.D0)+MD*DEC001o*(1.D0/3.D0)+MD*DEC00o*(1.D0/3.D0)&
	+MD*C002o*(-2.D0/3.D0)+MD*C001o*(-2.D0/3.D0)+MD*C00o*(-2.D0/3.D0)+t*DEBf001*(1.D0/108.D0)+t*Bf001*(-1.D0/18.D0)+t*DEC0022*(-5.D0/108.D0)&
	+t*DEC0012*(-1.D0/18.D0)+t*DEC0011*(-1.D0/108.D0)+t*DEC0000*(2.D0/27.D0)+t*C0022*(-5.D0/36.D0)+t*C0012*(-1.D0/6.D0)&
	+t*C0011*(-1.D0/36.D0)+t*C0000*(-4.D0/9.D0)+t*MD*DEC002*(1.D0/36.D0)+t*MD*DEC001*(1.D0/108.D0)+t*MD*C002*(1.D0/12.D0)&
	+t*MD*C001*(1.D0/36.D0)+x*DEBf001*(1.D0/27.D0)+x*Bf001*(-2.D0/9.D0)+x*DEC0022*(-1.D0/27.D0)+x*DEC0012*(-2.D0/27.D0)&
	+x*DEC0011*(-1.D0/27.D0)+x*DEC0000*(8.D0/27.D0)+x*C0022*(-1.D0/9.D0)+x*C0012*(-2.D0/9.D0)+x*C0011*(-1.D0/9.D0)+x*C0000*(-16.D0/9.D0)&
	+x*MD*DEC002*(1.D0/27.D0)+x*MD*DEC001*(1.D0/27.D0)+x*MD*C002*(1.D0/9.D0)+x*MD*C001*(1.D0/9.D0)+x*t*DEC0022*(-1.D0/54.D0)&
	+x*t*DEC0012*(-1.D0/108.D0)+x*t*C0022*(-1.D0/18.D0)+x*t*C0012*(-1.D0/9.D0)+x**2*DEC0022*(-2.D0/27.D0)+x**2*DEC0012*(-1.D0/27.D0)&
	+x**2*C0022*(-2.D0/9.D0)+x**2*C0012*(-4.D0/9.D0)

   A(2) = A(2) &
	+x**(-1)*DEBf001o*(5.D0/27.D0)+x**(-1)*DEBf001*(-5.D0/27.D0)+x**(-1)*Bf001o*(-4.D0/9.D0)+x**(-1)*Bf001*(4.D0/9.D0)&
	+x**(-1)*DEC0022o*(-11.D0/27.D0)+x**(-1)*DEC002o*(-1.D0/3.D0)+x**(-1)*DEC0012o*(-22.D0/27.D0)+x**(-1)*DEC0011o*(-11.D0/27.D0)&
	+x**(-1)*DEC001o*(-1.D0/3.D0)+x**(-1)*DEC0000o*(-4.D0/27.D0)+x**(-1)*DEC0022*(11.D0/27.D0)+x**(-1)*DEC002*(1.D0/3.D0)&
	+x**(-1)*DEC0012*(22.D0/27.D0)+x**(-1)*DEC0011*(11.D0/27.D0)+x**(-1)*DEC001*(1.D0/3.D0)+x**(-1)*DEC0000*(4.D0/27.D0)&
	+x**(-1)*C0022o*(4.D0/9.D0)+x**(-1)*C002o*(2.D0/3.D0)+x**(-1)*C0012o*(8.D0/9.D0)+x**(-1)*C0011o*(4.D0/9.D0)&
	+x**(-1)*C001o*(2.D0/3.D0)+x**(-1)*C0000o*(-16.D0/9.D0)+x**(-1)*C0022*(-4.D0/9.D0)+x**(-1)*C002*(-2.D0/3.D0)&
	+x**(-1)*C0012*(-8.D0/9.D0)+x**(-1)*C0011*(-4.D0/9.D0)+x**(-1)*C001*(-2.D0/3.D0)+x**(-1)*C0000*(16.D0/9.D0)&
	+x**(-1)*MD*DEBf00o*(-5.D0/27.D0)+x**(-1)*MD*DEBf00*(5.D0/27.D0)+x**(-1)*MD*Bf00o*(4.D0/9.D0)+x**(-1)*MD*Bf00*(-4.D0/9.D0)&
	+x**(-1)*MD*DEC002o*(11.D0/27.D0)+x**(-1)*MD*DEC001o*(11.D0/27.D0)+x**(-1)*MD*DEC00o*(1.D0/3.D0)+x**(-1)*MD*DEC002*(-11.D0/27.D0)&
	+x**(-1)*MD*DEC001*(-11.D0/27.D0)+x**(-1)*MD*DEC00*(-1.D0/3.D0)+x**(-1)*MD*C002o*(-4.D0/9.D0)+x**(-1)*MD*C001o*(-4.D0/9.D0)&
	+x**(-1)*MD*C00o*(-2.D0/3.D0)+x**(-1)*MD*C002*(4.D0/9.D0)+x**(-1)*MD*C001*(4.D0/9.D0)+x**(-1)*MD*C00*(2.D0/3.D0)&
	+DEC0022*(11.D0/54.D0)+DEC002*(4.D0/27.D0)+DEC0012*(11.D0/27.D0)+DEC0011*(11.D0/54.D0)+DEC001*(5.D0/27.D0)+C0022*(1.D0/9.D0)&
	+C002*(-2.D0/9.D0)+C0012*(-4.D0/9.D0)+C0011*(-5.D0/9.D0)+C001*(-4.D0/9.D0)+MD*DEC002*(-11.D0/54.D0)+MD*DEC001*(-11.D0/54.D0)&
	+MD*DEC00*(-5.D0/27.D0)+MD*C002*(2.D0/9.D0)+MD*C001*(5.D0/9.D0)+MD*C00*(4.D0/9.D0)+x*C1222*(-1.D0/9.D0)+x*C122*(-1.D0/9.D0)&
	+x*C1122*(-2.D0/9.D0)+x*C112*(-1.D0/9.D0)+x*C1112*(-1.D0/9.D0)+x*MD*C122*(1.D0/9.D0)+x*MD*C12*(1.D0/9.D0)+x*MD*C112*(1.D0/9.D0)

   A(3) = A(3) &
	+x**(-1)*t*DEBf001o*(-1.D0/108.D0)+x**(-1)*t*DEBf001*(1.D0/108.D0)+x**(-1)*t*Bf001o*(1.D0/18.D0)+x**(-1)*t*Bf001*(-1.D0/18.D0)&
	+x**(-1)*t*DEC0022o*(1.D0/27.D0)+x**(-1)*t*DEC0012o*(2.D0/27.D0)+x**(-1)*t*DEC0011o*(1.D0/27.D0)+x**(-1)*t*DEC0000o*(-5.D0/54.D0)&
	+x**(-1)*t*DEC0022*(-1.D0/27.D0)+x**(-1)*t*DEC0012*(-2.D0/27.D0)+x**(-1)*t*DEC0011*(-1.D0/27.D0)+x**(-1)*t*DEC0000*(5.D0/54.D0)&
	+x**(-1)*t*C0022o*(1.D0/9.D0)+x**(-1)*t*C0012o*(2.D0/9.D0)+x**(-1)*t*C0011o*(1.D0/9.D0)+x**(-1)*t*C0000o*(11.D0/9.D0)&
	+x**(-1)*t*C0022*(-1.D0/9.D0)+x**(-1)*t*C0012*(-2.D0/9.D0)+x**(-1)*t*C0011*(-1.D0/9.D0)+x**(-1)*t*C0000*(-11.D0/9.D0)&
	+x**(-1)*t*MD*DEBf00o*(1.D0/108.D0)+x**(-1)*t*MD*DEBf00*(-1.D0/108.D0)+x**(-1)*t*MD*Bf00o*(-1.D0/18.D0)+x**(-1)*t*MD*Bf00*(1.D0/18.D0)&
	+x**(-1)*t*MD*DEC002o*(-1.D0/27.D0)+x**(-1)*t*MD*DEC001o*(-1.D0/27.D0)+x**(-1)*t*MD*DEC002*(1.D0/27.D0)+x**(-1)*t*MD*DEC001*(1.D0/27.D0)&
	+x**(-1)*t*MD*C002o*(-1.D0/9.D0)+x**(-1)*t*MD*C001o*(-1.D0/9.D0)+x**(-1)*t*MD*C002*(1.D0/9.D0)+x**(-1)*t*MD*C001*(1.D0/9.D0)&
	+DEBf001o*(19.D0/54.D0)+DEBf001*(1.D0/54.D0)+Bf001o*(-7.D0/9.D0)+Bf001*(-1.D0/9.D0)+DEC0022o*(-1.D0/3.D0)+DEC002o*(-1.D0/3.D0)&
	+DEC0012o*(-2.D0/3.D0)+DEC0011o*(-1.D0/3.D0)+DEC001o*(-1.D0/3.D0)+DEC0000o*(-1.D0/3.D0)+DEC0000*(11.D0/27.D0)+C0022o*(2.D0/3.D0)&
	+C002o*(2.D0/3.D0)+C0012o*(4.D0/3.D0)+C0011o*(2.D0/3.D0)+C001o*(2.D0/3.D0)+C0000o*(2.D0/3.D0)+C0000*(-10.D0/9.D0)&
	+MD*DEBf00o*(-19.D0/54.D0)+MD*DEBf00*(1.D0/54.D0)+MD*Bf00o*(7.D0/9.D0)+MD*Bf00*(-1.D0/9.D0)+MD*DEC002o*(1.D0/3.D0)&
	+MD*DEC001o*(1.D0/3.D0)+MD*DEC00o*(1.D0/3.D0)+MD*C002o*(-2.D0/3.D0)+MD*C001o*(-2.D0/3.D0)+MD*C00o*(-2.D0/3.D0)&
	+t*DEBf001*(1.D0/108.D0)+t*Bf001*(-1.D0/18.D0)+t*DEC0022*(-5.D0/108.D0)+t*DEC0012*(-1.D0/18.D0)+t*DEC0011*(-1.D0/108.D0)&
	+t*DEC0000*(2.D0/27.D0)+t*C0022*(-5.D0/36.D0)+t*C0012*(-1.D0/6.D0)+t*C0011*(-1.D0/36.D0)+t*C0000*(-4.D0/9.D0)&
	+t*MD*DEC002*(1.D0/36.D0)+t*MD*DEC001*(1.D0/108.D0)+t*MD*C002*(1.D0/12.D0)+t*MD*C001*(1.D0/36.D0)&
	+x*DEBf001*(1.D0/27.D0)+x*Bf001*(-2.D0/9.D0)+x*DEC0022*(-1.D0/27.D0)+x*DEC0012*(-2.D0/27.D0)+x*DEC0011*(-1.D0/27.D0)&
	+x*DEC0000*(8.D0/27.D0)+x*C0022*(-1.D0/9.D0)+x*C0012*(-2.D0/9.D0)+x*C0011*(-1.D0/9.D0)+x*C0000*(-16.D0/9.D0)&
	+x*MD*DEC002*(1.D0/27.D0)+x*MD*DEC001*(1.D0/27.D0)+x*MD*C002*(1.D0/9.D0)+x*MD*C001*(1.D0/9.D0)+x*t*DEC0022*(-1.D0/54.D0)&
	+x*t*DEC0012*(-1.D0/108.D0)+x*t*C0022*(-1.D0/18.D0)+x*t*C0012*(-1.D0/9.D0)+x**2*DEC0022*(-2.D0/27.D0)&
	+x**2*DEC0012*(-1.D0/27.D0)+x**2*C0022*(-2.D0/9.D0)+x**2*C0012*(-4.D0/9.D0)

   A(4) = A(4) &
	+x**(-1)*DEBf001o*(-1.D0/54.D0)+x**(-1)*DEBf001*(1.D0/54.D0)+x**(-1)*Bf001o*(1.D0/9.D0)+x**(-1)*Bf001*(-1.D0/9.D0)&
	+x**(-1)*DEC0022o*(2.D0/27.D0)+x**(-1)*DEC0012o*(4.D0/27.D0)+x**(-1)*DEC0011o*(2.D0/27.D0)+x**(-1)*DEC0000o*(-5.D0/27.D0)&
	+x**(-1)*DEC0022*(-2.D0/27.D0)+x**(-1)*DEC0012*(-4.D0/27.D0)+x**(-1)*DEC0011*(-2.D0/27.D0)+x**(-1)*DEC0000*(5.D0/27.D0)&
	+x**(-1)*C0022o*(2.D0/9.D0)+x**(-1)*C0012o*(4.D0/9.D0)+x**(-1)*C0011o*(2.D0/9.D0)+x**(-1)*C0000o*(22.D0/9.D0)&
	+x**(-1)*C0022*(-2.D0/9.D0)+x**(-1)*C0012*(-4.D0/9.D0)+x**(-1)*C0011*(-2.D0/9.D0)+x**(-1)*C0000*(-22.D0/9.D0)&
	+x**(-1)*MD*DEBf00o*(1.D0/54.D0)+x**(-1)*MD*DEBf00*(-1.D0/54.D0)+x**(-1)*MD*Bf00o*(-1.D0/9.D0)+x**(-1)*MD*Bf00*(1.D0/9.D0)&
	+x**(-1)*MD*DEC002o*(-2.D0/27.D0)+x**(-1)*MD*DEC001o*(-2.D0/27.D0)+x**(-1)*MD*DEC002*(2.D0/27.D0)+x**(-1)*MD*DEC001*(2.D0/27.D0)&
	+x**(-1)*MD*C002o*(-2.D0/9.D0)+x**(-1)*MD*C001o*(-2.D0/9.D0)+x**(-1)*MD*C002*(2.D0/9.D0)+x**(-1)*MD*C001*(2.D0/9.D0)&
	+DEBf001*(1.D0/54.D0)+Bf001*(-1.D0/9.D0)+DEC0022*(-5.D0/54.D0)+DEC0012*(-1.D0/9.D0)+DEC0011*(-1.D0/54.D0)+DEC0000*(4.D0/27.D0)&
	+C0022*(-5.D0/18.D0)+C0012*(-1.D0/3.D0)+C0011*(-1.D0/18.D0)+C0000*(-8.D0/9.D0)+MD*DEC002*(1.D0/18.D0)+MD*DEC001*(1.D0/54.D0)&
	+MD*C002*(1.D0/6.D0)+MD*C001*(1.D0/18.D0)+x*DEC0022*(-1.D0/27.D0)+x*DEC0012*(-1.D0/54.D0)+x*C0022*(-1.D0/9.D0)&
	+x*C0012*(-2.D0/9.D0)

   A(5) = A(5) &
	+x**(-1)*DEBf001o*(1.D0/54.D0)+x**(-1)*DEBf001*(-1.D0/54.D0)+x**(-1)*Bf001o*(-1.D0/9.D0)+x**(-1)*Bf001*(1.D0/9.D0)&
	+x**(-1)*DEC0022o*(-2.D0/27.D0)+x**(-1)*DEC0012o*(-4.D0/27.D0)+x**(-1)*DEC0011o*(-2.D0/27.D0)+x**(-1)*DEC0000o*(5.D0/27.D0)&
	+x**(-1)*DEC0022*(2.D0/27.D0)+x**(-1)*DEC0012*(4.D0/27.D0)+x**(-1)*DEC0011*(2.D0/27.D0)+x**(-1)*DEC0000*(-5.D0/27.D0)&
	+x**(-1)*C0022o*(-2.D0/9.D0)+x**(-1)*C0012o*(-4.D0/9.D0)+x**(-1)*C0011o*(-2.D0/9.D0)+x**(-1)*C0000o*(-22.D0/9.D0)&
	+x**(-1)*C0022*(2.D0/9.D0)+x**(-1)*C0012*(4.D0/9.D0)+x**(-1)*C0011*(2.D0/9.D0)+x**(-1)*C0000*(22.D0/9.D0)&
	+x**(-1)*MD*DEBf00o*(-1.D0/54.D0)+x**(-1)*MD*DEBf00*(1.D0/54.D0)+x**(-1)*MD*Bf00o*(1.D0/9.D0)+x**(-1)*MD*Bf00*(-1.D0/9.D0)&
	+x**(-1)*MD*DEC002o*(2.D0/27.D0)+x**(-1)*MD*DEC001o*(2.D0/27.D0)+x**(-1)*MD*DEC002*(-2.D0/27.D0)+x**(-1)*MD*DEC001*(-2.D0/27.D0)&
	+x**(-1)*MD*C002o*(2.D0/9.D0)+x**(-1)*MD*C001o*(2.D0/9.D0)+x**(-1)*MD*C002*(-2.D0/9.D0)+x**(-1)*MD*C001*(-2.D0/9.D0)&
	+DEBf001*(1.D0/36.D0)+DEC0022*(-1.D0/108.D0)+DEC002*(-5.D0/54.D0)+DEC0012*(1.D0/108.D0)+DEC0011*(1.D0/54.D0)&
	+DEC001*(1.D0/108.D0)+DEC0000*(-1.D0/6.D0)+C0022*(17.D0/36.D0)+C002*(2.D0/9.D0)+C0012*(4.D0/9.D0)+C0011*(-1.D0/36.D0)&
	+C001*(-1.D0/18.D0)+C0000*(2.D0/3.D0)+MD*DEC002*(-1.D0/18.D0)+MD*DEC001*(-1.D0/54.D0)+MD*DEC00*(-1.D0/108.D0)&
	+MD*C002*(-1.D0/12.D0)+MD*C001*(1.D0/36.D0)+MD*C00*(1.D0/18.D0)+x*DEC0022*(-7.D0/108.D0)+x*DEC002*(-5.D0/54.D0)&
	+x*DEC0012*(-11.D0/108.D0)+x*C0022*(2.D0/9.D0)+x*C002*(2.D0/9.D0)+x*C0012*(4.D0/9.D0)+x**2*C1222*(1.D0/18.D0)+x**2*C122*(1.D0/18.D0)&
	+x**2*C1122*(1.D0/18.D0)

   A(6) = A(6) &
	+x**(-1)*DEBf001o*(2.D0/27.D0)+x**(-1)*DEBf001*(-2.D0/27.D0)+x**(-1)*Bf001o*(-1.D0/9.D0)+x**(-1)*Bf001*(1.D0/9.D0)&
	+x**(-1)*DEC0022o*(-7.D0/54.D0)+x**(-1)*DEC002o*(-1.D0/6.D0)+x**(-1)*DEC0012o*(-7.D0/27.D0)+x**(-1)*DEC0011o*(-7.D0/54.D0)&
	+x**(-1)*DEC001o*(-1.D0/6.D0)+x**(-1)*DEC0000o*(-7.D0/27.D0)+x**(-1)*DEC0022*(7.D0/54.D0)+x**(-1)*DEC002*(1.D0/6.D0)&
	+x**(-1)*DEC0012*(7.D0/27.D0)+x**(-1)*DEC0011*(7.D0/54.D0)+x**(-1)*DEC001*(1.D0/6.D0)+x**(-1)*DEC0000*(7.D0/27.D0)&
	+x**(-1)*C0022o*(4.D0/9.D0)+x**(-1)*C002o*(1.D0/3.D0)+x**(-1)*C0012o*(8.D0/9.D0)+x**(-1)*C0011o*(4.D0/9.D0)&
	+x**(-1)*C001o*(1.D0/3.D0)+x**(-1)*C0000o*(14.D0/9.D0)+x**(-1)*C0022*(-4.D0/9.D0)+x**(-1)*C002*(-1.D0/3.D0)&
	+x**(-1)*C0012*(-8.D0/9.D0)+x**(-1)*C0011*(-4.D0/9.D0)+x**(-1)*C001*(-1.D0/3.D0)+x**(-1)*C0000*(-14.D0/9.D0)&
	+x**(-1)*MD*DEBf00o*(-2.D0/27.D0)+x**(-1)*MD*DEBf00*(2.D0/27.D0)+x**(-1)*MD*Bf00o*(1.D0/9.D0)+x**(-1)*MD*Bf00*(-1.D0/9.D0)&
	+x**(-1)*MD*DEC002o*(7.D0/54.D0)+x**(-1)*MD*DEC001o*(7.D0/54.D0)+x**(-1)*MD*DEC00o*(1.D0/6.D0)+x**(-1)*MD*DEC002*(-7.D0/54.D0)&
	+x**(-1)*MD*DEC001*(-7.D0/54.D0)+x**(-1)*MD*DEC00*(-1.D0/6.D0)+x**(-1)*MD*C002o*(-4.D0/9.D0)+x**(-1)*MD*C001o*(-4.D0/9.D0)&
	+x**(-1)*MD*C00o*(-1.D0/3.D0)+x**(-1)*MD*C002*(4.D0/9.D0)+x**(-1)*MD*C001*(4.D0/9.D0)+x**(-1)*MD*C00*(1.D0/3.D0)&
	+DEBf001*(-1.D0/36.D0)+DEC0022*(1.D0/9.D0)+DEC002*(1.D0/6.D0)+DEC0012*(7.D0/36.D0)+DEC0011*(1.D0/12.D0)+DEC001*(1.D0/12.D0)&
	+DEC0000*(1.D0/6.D0)+C0022*(-5.D0/12.D0)+C002*(-1.D0/3.D0)+C0012*(-2.D0/3.D0)+C0011*(-1.D0/4.D0)+C001*(-1.D0/6.D0)&
	+C0000*(-2.D0/3.D0)+MD*DEC002*(-5.D0/108.D0)+MD*DEC001*(-1.D0/12.D0)+MD*DEC00*(-1.D0/12.D0)+MD*C002*(7.D0/36.D0)&
	+MD*C001*(1.D0/4.D0)+MD*C00*(1.D0/6.D0)+x*DEC0022*(7.D0/108.D0)+x*DEC002*(5.D0/54.D0)+x*DEC0012*(11.D0/108.D0)+x*C0022*(-2.D0/9.D0)&
	+x*C002*(-2.D0/9.D0)+x*C1222*(-1.D0/18.D0)+x*C122*(-1.D0/18.D0)+x*C0012*(-4.D0/9.D0)+x*C1122*(-1.D0/9.D0)+x*C112*(-1.D0/18.D0)&
	+x*C1112*(-1.D0/18.D0)+x*MD*C122*(1.D0/18.D0)+x*MD*C12*(1.D0/18.D0)+x*MD*C112*(1.D0/18.D0)+x**2*C1222*(-1.D0/18.D0)&
	+x**2*C122*(-1.D0/18.D0)+x**2*C1122*(-1.D0/18.D0)

   A(7) = A(7)

   A(8) = A(8) &
	+x**(-1)*DEBf001o*(1.D0/54.D0)+x**(-1)*DEBf001*(-1.D0/54.D0)+x**(-1)*Bf001o*(-1.D0/9.D0)+x**(-1)*Bf001*(1.D0/9.D0)&
	+x**(-1)*DEC0022o*(-2.D0/27.D0)+x**(-1)*DEC0012o*(-4.D0/27.D0)+x**(-1)*DEC0011o*(-2.D0/27.D0)+x**(-1)*DEC0000o*(5.D0/27.D0)&
	+x**(-1)*DEC0022*(2.D0/27.D0)+x**(-1)*DEC0012*(4.D0/27.D0)+x**(-1)*DEC0011*(2.D0/27.D0)+x**(-1)*DEC0000*(-5.D0/27.D0)&
	+x**(-1)*C0022o*(-2.D0/9.D0)+x**(-1)*C0012o*(-4.D0/9.D0)+x**(-1)*C0011o*(-2.D0/9.D0)+x**(-1)*C0000o*(-22.D0/9.D0)&
	+x**(-1)*C0022*(2.D0/9.D0)+x**(-1)*C0012*(4.D0/9.D0)+x**(-1)*C0011*(2.D0/9.D0)+x**(-1)*C0000*(22.D0/9.D0)&
	+x**(-1)*MD*DEBf00o*(-1.D0/54.D0)+x**(-1)*MD*DEBf00*(1.D0/54.D0)+x**(-1)*MD*Bf00o*(1.D0/9.D0)+x**(-1)*MD*Bf00*(-1.D0/9.D0)&
	+x**(-1)*MD*DEC002o*(2.D0/27.D0)+x**(-1)*MD*DEC001o*(2.D0/27.D0)+x**(-1)*MD*DEC002*(-2.D0/27.D0)+x**(-1)*MD*DEC001*(-2.D0/27.D0)&
	+x**(-1)*MD*C002o*(2.D0/9.D0)+x**(-1)*MD*C001o*(2.D0/9.D0)+x**(-1)*MD*C002*(-2.D0/9.D0)+x**(-1)*MD*C001*(-2.D0/9.D0)&
	+DEBf001*(-1.D0/54.D0)+Bf001*(1.D0/9.D0)+DEC0022*(5.D0/54.D0)+DEC0012*(1.D0/9.D0)+DEC0011*(1.D0/54.D0)+DEC0000*(-4.D0/27.D0)&
	+C0022*(5.D0/18.D0)+C0012*(1.D0/3.D0)+C0011*(1.D0/18.D0)+C0000*(8.D0/9.D0)+MD*DEC002*(-1.D0/18.D0)+MD*DEC001*(-1.D0/54.D0)&
	+MD*C002*(-1.D0/6.D0)+MD*C001*(-1.D0/18.D0)+x*DEC0022*(1.D0/27.D0)+x*DEC0012*(1.D0/54.D0)+x*C0022*(1.D0/9.D0)+x*C0012*(2.D0/9.D0)

!    AA=AA-3.D0*A/2.D0
!    A=0.D0

   !Sum2=(A(1)-Sum3)/Sum1
   !Sum3=A(1)

   !write(*,*) '4', Sum2

   !* Graph 5:
   !*  Self-energy graph
   !* Isospin factor = 2

   A(1) = A(1) &
	+x**(-2)*t*DEBf001o*(1.D0/4.D0)+x**(-2)*t*DEBf001*(-1.D0/4.D0)+x**(-2)*t*Bf001o*(-1.D0/2.D0)+x**(-2)*t*Bf001*(1.D0/2.D0)&
	+x**(-2)*t*MD*DEBf00o*(-1.D0/4.D0)+x**(-2)*t*MD*DEBf00*(1.D0/4.D0)+x**(-2)*t*MD*Bf00o*(1.D0/2.D0)+x**(-2)*t*MD*Bf00*(-1.D0/2.D0)&
	+x**(-1)*t*DEC0022o*(-1.D0/4.D0)+x**(-1)*t*DEC002o*(-1.D0/4.D0)+x**(-1)*t*DEC0012o*(-1.D0/2.D0)+x**(-1)*t*DEC0011o*(-1.D0/4.D0)&
	+x**(-1)*t*DEC001o*(-1.D0/4.D0)+x**(-1)*t*DEC0000o*(-1.D0/4.D0)+x**(-1)*t*C0022o*(1.D0/2.D0)+x**(-1)*t*C002o*(1.D0/2.D0)&
	+x**(-1)*t*C0012o*(1.D0)+x**(-1)*t*C0011o*(1.D0/2.D0)+x**(-1)*t*C001o*(1.D0/2.D0)+x**(-1)*t*C0000o*(1.D0/2.D0)&
	+x**(-1)*t*DEBf001o*(1.D0/4.D0)+x**(-1)*t*DEBf001*(-3.D0/8.D0)+x**(-1)*t*Bf001o*(-1.D0/2.D0)+x**(-1)*t*Bf001*(3.D0/4.D0)&
	+x**(-1)*t*MD*DEC002o*(1.D0/4.D0)+x**(-1)*t*MD*DEC001o*(1.D0/4.D0)+x**(-1)*t*MD*DEC00o*(1.D0/4.D0)+x**(-1)*t*MD*C002o*(-1.D0/2.D0)&
	+x**(-1)*t*MD*C001o*(-1.D0/2.D0)+x**(-1)*t*MD*C00o*(-1.D0/2.D0)+x**(-1)*t*MD*DEBf00o*(-1.D0/4.D0)+x**(-1)*t*MD*DEBf00*(1.D0/4.D0)&
	+x**(-1)*t*MD*Bf00o*(1.D0/2.D0)+x**(-1)*t*MD*Bf00*(-1.D0/2.D0)+DEC0022o*(-1.D0/2.D0)+DEC002o*(-1.D0/2.D0)+DEC0012o*(-1.D0)&
	+DEC0011o*(-1.D0/2.D0)+DEC001o*(-1.D0/2.D0)+DEC0000o*(-1.D0/2.D0)+C0022o*(1.D0)+C002o*(1.D0)+C0012o*(2.D0)+C0011o*(1.D0)&
	+C001o*(1.D0)+C0000o*(1.D0)+DEBf001o*(1.D0/2.D0)+DEBf001*(-1.D0/4.D0)+Bf001o*(-1.D0)+Bf001*(1.D0/2.D0)+MD*DEC002o*(1.D0/2.D0)&
	+MD*DEC001o*(1.D0/2.D0)+MD*DEC00o*(1.D0/2.D0)+MD*C002o*(-1.D0)+MD*C001o*(-1.D0)+MD*C00o*(-1.D0)+MD*DEBf00o*(-1.D0/2.D0)&
	+MD*Bf00o*(1.D0)+t*DEBf001*(-1.D0/8.D0)+t*Bf001*(1.D0/4.D0)+x*DEBf001*(-1.D0/2.D0)+x*Bf001*(1.D0)

   A(2) = A(2)

   A(3) = A(3) &
	+x**(-2)*t*DEBf001o*(1.D0/4.D0)+x**(-2)*t*DEBf001*(-1.D0/4.D0)+x**(-2)*t*Bf001o*(-1.D0/2.D0)+x**(-2)*t*Bf001*(1.D0/2.D0)&
	+x**(-2)*t*MD*DEBf00o*(-1.D0/4.D0)+x**(-2)*t*MD*DEBf00*(1.D0/4.D0)+x**(-2)*t*MD*Bf00o*(1.D0/2.D0)+x**(-2)*t*MD*Bf00*(-1.D0/2.D0)&
	+x**(-1)*t*DEC0022o*(-1.D0/4.D0)+x**(-1)*t*DEC002o*(-1.D0/4.D0)+x**(-1)*t*DEC0012o*(-1.D0/2.D0)+x**(-1)*t*DEC0011o*(-1.D0/4.D0)&
	+x**(-1)*t*DEC001o*(-1.D0/4.D0)+x**(-1)*t*DEC0000o*(-1.D0/4.D0)+x**(-1)*t*C0022o*(1.D0/2.D0)+x**(-1)*t*C002o*(1.D0/2.D0)&
	+x**(-1)*t*C0012o*(1.D0)+x**(-1)*t*C0011o*(1.D0/2.D0)+x**(-1)*t*C001o*(1.D0/2.D0)+x**(-1)*t*C0000o*(1.D0/2.D0)&
	+x**(-1)*t*DEBf001o*(1.D0/4.D0)+x**(-1)*t*DEBf001*(-3.D0/8.D0)+x**(-1)*t*Bf001o*(-1.D0/2.D0)+x**(-1)*t*Bf001*(3.D0/4.D0)&
	+x**(-1)*t*MD*DEC002o*(1.D0/4.D0)+x**(-1)*t*MD*DEC001o*(1.D0/4.D0)+x**(-1)*t*MD*DEC00o*(1.D0/4.D0)+x**(-1)*t*MD*C002o*(-1.D0/2.D0)&
	+x**(-1)*t*MD*C001o*(-1.D0/2.D0)+x**(-1)*t*MD*C00o*(-1.D0/2.D0)+x**(-1)*t*MD*DEBf00o*(-1.D0/4.D0)+x**(-1)*t*MD*DEBf00*(1.D0/4.D0)&
	+x**(-1)*t*MD*Bf00o*(1.D0/2.D0)+x**(-1)*t*MD*Bf00*(-1.D0/2.D0)+DEC0022o*(-1.D0/2.D0)+DEC002o*(-1.D0/2.D0)+DEC0012o*(-1.D0)&
	+DEC0011o*(-1.D0/2.D0)+DEC001o*(-1.D0/2.D0)+DEC0000o*(-1.D0/2.D0)+C0022o*(1.D0)+C002o*(1.D0)+C0012o*(2.D0)+C0011o*(1.D0)&
	+C001o*(1.D0)+C0000o*(1.D0)+DEBf001o*(1.D0/2.D0)+DEBf001*(-1.D0/4.D0)+Bf001o*(-1.D0)+Bf001*(1.D0/2.D0)+MD*DEC002o*(1.D0/2.D0)&
	+MD*DEC001o*(1.D0/2.D0)+MD*DEC00o*(1.D0/2.D0)+MD*C002o*(-1.D0)+MD*C001o*(-1.D0)+MD*C00o*(-1.D0)+MD*DEBf00o*(-1.D0/2.D0)&
	+MD*Bf00o*(1.D0)+t*DEBf001*(-1.D0/8.D0)+t*Bf001*(1.D0/4.D0)+x*DEBf001*(-1.D0/2.D0)+x*Bf001*(1.D0)

   A(4) = A(4) &
	+x**(-2)*DEBf001o*(1.D0/2.D0)+x**(-2)*DEBf001*(-1.D0/2.D0)+x**(-2)*Bf001o*(-1.D0)+x**(-2)*Bf001*(1.D0)+x**(-2)*MD*DEBf00o*(-1.D0/2.D0)&
	+x**(-2)*MD*DEBf00*(1.D0/2.D0)+x**(-2)*MD*Bf00o*(1.D0)+x**(-2)*MD*Bf00*(-1.D0)+x**(-1)*DEC0022o*(-1.D0/2.D0)+x**(-1)*DEC002o*(-1.D0/2.D0)&
	+x**(-1)*DEC0012o*(-1.D0)+x**(-1)*DEC0011o*(-1.D0/2.D0)+x**(-1)*DEC001o*(-1.D0/2.D0)+x**(-1)*DEC0000o*(-1.D0/2.D0)&
	+x**(-1)*C0022o*(1.D0)+x**(-1)*C002o*(1.D0)+x**(-1)*C0012o*(2.D0)+x**(-1)*C0011o*(1.D0)+x**(-1)*C001o*(1.D0)&
	+x**(-1)*C0000o*(1.D0)+x**(-1)*DEBf001o*(1.D0/2.D0)+x**(-1)*DEBf001*(-3.D0/4.D0)+x**(-1)*Bf001o*(-1.D0)+x**(-1)*Bf001*(3.D0/2.D0)&
	+x**(-1)*MD*DEC002o*(1.D0/2.D0)+x**(-1)*MD*DEC001o*(1.D0/2.D0)+x**(-1)*MD*DEC00o*(1.D0/2.D0)+x**(-1)*MD*C002o*(-1.D0)&
	+x**(-1)*MD*C001o*(-1.D0)+x**(-1)*MD*C00o*(-1.D0)+x**(-1)*MD*DEBf00o*(-1.D0/2.D0)+x**(-1)*MD*DEBf00*(1.D0/2.D0)&
	+x**(-1)*MD*Bf00o*(1.D0)+x**(-1)*MD*Bf00*(-1.D0)+DEBf001*(-1.D0/4.D0)+Bf001*(1.D0/2.D0)

   A(5) = A(5) &
	+x**(-2)*DEBf001o*(-1.D0/2.D0)+x**(-2)*DEBf001*(1.D0/2.D0)+x**(-2)*Bf001o*(1.D0)+x**(-2)*Bf001*(-1.D0)+x**(-2)*MD*DEBf00o*(1.D0/2.D0)&
	+x**(-2)*MD*DEBf00*(-1.D0/2.D0)+x**(-2)*MD*Bf00o*(-1.D0)+x**(-2)*MD*Bf00*(1.D0)+x**(-1)*DEC0022o*(1.D0/2.D0)+x**(-1)*DEC002o*(1.D0/2.D0)&
	+x**(-1)*DEC0012o*(1.D0)+x**(-1)*DEC0011o*(1.D0/2.D0)+x**(-1)*DEC001o*(1.D0/2.D0)+x**(-1)*DEC0000o*(1.D0/2.D0)+x**(-1)*C0022o*(-1.D0)&
	+x**(-1)*C002o*(-1.D0)+x**(-1)*C0012o*(-2.D0)+x**(-1)*C0011o*(-1.D0)+x**(-1)*C001o*(-1.D0)+x**(-1)*C0000o*(-1.D0)&
	+x**(-1)*DEBf001o*(-1.D0/2.D0)+x**(-1)*DEBf001*(3.D0/4.D0)+x**(-1)*Bf001o*(1.D0)+x**(-1)*Bf001*(-3.D0/2.D0)+x**(-1)*MD*DEC002o*(-1.D0/2.D0)&
	+x**(-1)*MD*DEC001o*(-1.D0/2.D0)+x**(-1)*MD*DEC00o*(-1.D0/2.D0)+x**(-1)*MD*C002o*(1.D0)+x**(-1)*MD*C001o*(1.D0)&
	+x**(-1)*MD*C00o*(1.D0)+x**(-1)*MD*DEBf00o*(1.D0/2.D0)+x**(-1)*MD*DEBf00*(-1.D0/2.D0)+x**(-1)*MD*Bf00o*(-1.D0)&
	+x**(-1)*MD*Bf00*(1.D0)+DEBf001*(1.D0/4.D0)+Bf001*(-1.D0/2.D0)

   A(6) = A(6) &
	+x**(-2)*DEBf001o*(1.D0/2.D0)+x**(-2)*DEBf001*(-1.D0/2.D0)+x**(-2)*Bf001o*(-1.D0)+x**(-2)*Bf001*(1.D0)+x**(-2)*MD*DEBf00o*(-1.D0/2.D0)&
	+x**(-2)*MD*DEBf00*(1.D0/2.D0)+x**(-2)*MD*Bf00o*(1.D0)+x**(-2)*MD*Bf00*(-1.D0)+x**(-1)*DEC0022o*(-1.D0/2.D0)&
	+x**(-1)*DEC002o*(-1.D0/2.D0)+x**(-1)*DEC0012o*(-1.D0)+x**(-1)*DEC0011o*(-1.D0/2.D0)+x**(-1)*DEC001o*(-1.D0/2.D0)&
	+x**(-1)*DEC0000o*(-1.D0/2.D0)+x**(-1)*C0022o*(1.D0)+x**(-1)*C002o*(1.D0)+x**(-1)*C0012o*(2.D0)+x**(-1)*C0011o*(1.D0)&
	+x**(-1)*C001o*(1.D0)+x**(-1)*C0000o*(1.D0)+x**(-1)*DEBf001o*(1.D0/2.D0)+x**(-1)*DEBf001*(-3.D0/4.D0)+x**(-1)*Bf001o*(-1.D0)&
	+x**(-1)*Bf001*(3.D0/2.D0)+x**(-1)*MD*DEC002o*(1.D0/2.D0)+x**(-1)*MD*DEC001o*(1.D0/2.D0)+x**(-1)*MD*DEC00o*(1.D0/2.D0)&
	+x**(-1)*MD*C002o*(-1.D0)+x**(-1)*MD*C001o*(-1.D0)+x**(-1)*MD*C00o*(-1.D0)+x**(-1)*MD*DEBf00o*(-1.D0/2.D0)&
	+x**(-1)*MD*DEBf00*(1.D0/2.D0)+x**(-1)*MD*Bf00o*(1.D0)+x**(-1)*MD*Bf00*(-1.D0)+DEBf001*(-1.D0/4.D0)+Bf001*(1.D0/2.D0)

   A(7) = A(7)

   A(8) = A(8) &
	+x**(-2)*DEBf001o*(-1.D0/2.D0)+x**(-2)*DEBf001*(1.D0/2.D0)+x**(-2)*Bf001o*(1.D0)+x**(-2)*Bf001*(-1.D0)+x**(-2)*MD*DEBf00o*(1.D0/2.D0)&
	+x**(-2)*MD*DEBf00*(-1.D0/2.D0)+x**(-2)*MD*Bf00o*(-1.D0)+x**(-2)*MD*Bf00*(1.D0)+x**(-1)*DEC0022o*(1.D0/2.D0)+x**(-1)*DEC002o*(1.D0/2.D0)&
	+x**(-1)*DEC0012o*(1.D0)+x**(-1)*DEC0011o*(1.D0/2.D0)+x**(-1)*DEC001o*(1.D0/2.D0)+x**(-1)*DEC0000o*(1.D0/2.D0)&
	+x**(-1)*C0022o*(-1.D0)+x**(-1)*C002o*(-1.D0)+x**(-1)*C0012o*(-2.D0)+x**(-1)*C0011o*(-1.D0)+x**(-1)*C001o*(-1.D0)&
	+x**(-1)*C0000o*(-1.D0)+x**(-1)*DEBf001o*(-1.D0/2.D0)+x**(-1)*DEBf001*(3.D0/4.D0)+x**(-1)*Bf001o*(1.D0)+x**(-1)*Bf001*(-3.D0/2.D0)&
	+x**(-1)*MD*DEC002o*(-1.D0/2.D0)+x**(-1)*MD*DEC001o*(-1.D0/2.D0)+x**(-1)*MD*DEC00o*(-1.D0/2.D0)+x**(-1)*MD*C002o*(1.D0)&
	+x**(-1)*MD*C001o*(1.D0)+x**(-1)*MD*C00o*(1.D0)+x**(-1)*MD*DEBf00o*(1.D0/2.D0)+x**(-1)*MD*DEBf00*(-1.D0/2.D0)+x**(-1)*MD*Bf00o*(-1.D0)&
	+x**(-1)*MD*Bf00*(1.D0)+DEBf001*(1.D0/4.D0)+Bf001*(-1.D0/2.D0)

!    AA=AA+A/2.D0
!    A=0.d0

   !Sum2=(A(1)-Sum3)/Sum1
   !Sum3=A(1)

   !write(*,*) '5', Sum2


   !* Graph 6:
   !*  Two-photon vertex graph
   !* Isospin factor = 4/3

   p1s=nuclm**2
   p2s=t
   p1p2s=nuclm**2

   m1s=deltam**2
   m2s=mpi2
   m3s=mpi2

   CALL Setdelta(Delta)

   C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C11=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C22=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C111=C0i(cc111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C112=C0i(cc112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C122=C0i(cc122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C222=C0i(cc222,p1s,p2s,p1p2s,m1s,m2s,m3s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   C00D=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001D=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002D=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)

   DEC00=2.D0*(C00-C00D)/sw1
   DEC001=2.D0*(C001-C001D)/sw1
   DEC002=2.D0*(C002-C002D)/sw1

   A(1) = A(1) &
	+(DEC002*(-2.D0/3.D0)+DEC001*(-2.D0/3.D0)+C002*(4.D0/3.D0)+C001*(4.D0/3.D0)+MD*DEC00*(2.D0/3.D0)&
	+MD*C00*(-4.D0/3.D0)+t*DEC002*(10.D0/27.D0)+t*DEC001*(10.D0/27.D0)+t*C002*(-8.D0/9.D0)+t*C122*(-4.D0/9.D0)&
	+t*C112*(-4.D0/9.D0)+t*C001*(-8.D0/9.D0)+t*MD*DEC00*(-10.D0/27.D0)+t*MD*C12*(4.D0/9.D0)+t*MD*C00*(8.D0/9.D0)&
	+t**2*C122*(1.D0/9.D0)+t**2*C112*(1.D0/9.D0)+t**2*MD*C12*(-1.D0/9.D0))/2.D0

!    Sum1=-1.D0/3.D0*DEC0022o-1.D0/3.D0*DEC002o-2.D0/3.D0*DEC0012o-1.D0/3.D0*DEC0011o-1.D0/3.D0*DEC001o-1.D0/3.D0*DEC0000o&
! 	+2.D0/3.D0*C0022o+2.D0/3.D0*C002o+4.D0/3.D0*C0012o+2.D0/3.D0*C0011o+2.D0/3.D0*C001o+2.D0/3.D0*C0000o+1.D0/3.D0*MD*DEC002o+1.D0/3.D0*MD*DEC001o&
! 	+1.D0/3.D0*MD*DEC00o-2.D0/3.D0*MD*C002o-2.D0/3.D0*MD*C001o-2.D0/3.D0*MD*C00o+1.D0/6.D0*DEB001o-1.D0/3.D0*B001o&
! 	-1.D0/6.D0*MD*DEB00o+1.D0/3.D0*MD*B00o
! 
!    Sum2=-DEC0022o-DEC002o-2.D0*DEC0012o-DEC0011o-DEC001o-8.D0*DEC0000o+2.D0*C0022o+2.D0*C002o+4.D0*C0012o+2.D0*C0011o+2.D0*C001o&
!   	+12.D0*C0000o+MD*DEC002o+MD*DEC001o+MD*DEC00o-2.D0*MD*C002o-2.D0*MD*C001o-2.D0*MD*C00o+MD*DEB00o-2.D0*MD*B00o

!   Sum1=(+2.D0*C0022o+2.D0*C002o+4.D0*C0012o+2.D0*C0011o+2.D0*C001o+2.D0*C0000o-2.D0*MD*C002o-2.D0*MD*C001o-2.D0*MD*C00o-2.D0*B001o+2.D0*MD*B00o)/3.D0

!   Sum2= -2.D0*C0022o-2.D0*C002o-4.D0*C0012o-2.D0*C0011o-2.D0*C001o-2.D0*C0000o+2.D0*MD*C002o+2.D0*MD*C001o+2.D0*MD*C00o+2.D0*B001o-2.D0*MD*B00o

!   write(*,*) Sum1, Sum2, Sum2/Sum1, Sum2-Sum1
!   stop

!    A=AA+A*3.D0/4.D0

   !Sum2=(A(1)-Sum3)/Sum1
   !Sum3=A(1)

   !write(*,*) '6', Sum2

   A=gafac*A


 END SUBROUTINE CovOneDeltaLoop

!
!

!
! Isospin factors of 1PR graphs are fudged to provide current conservation
!
SUBROUTINE CorrIsoCovOneDeltaLoop(s,t,mpi2,gafac,A,Delta)
  USE trigd
!!  IMPLICIT COMPLEX*16 B, C, D
  implicit COMPLEX*16 (B-D)
#include "looptools.h"
  REAL*8, INTENT(in) :: s,t,mpi2,gafac,Delta
  COMPLEX*16,INTENT(out) :: A(8)
  COMPLEX*16, PARAMETER :: CI=(0.d0,1.d0)
  REAL*8 x,nuclm,deltam, MD
  REAL*8 DeltaAux, sw1,facv,face
  !* [v2] Delta mass inside the loops taken from the main program (baryon.dat), audit S1;
  !*      LoopTools UV parameter restored on exit, audit B5
  REAL*8 dmloop, DeltaIn
  COMMON /v2loopmass/ dmloop
  REAL*8 m1s,m2s,m3s,m4s
  REAL*8 p1s,p2s,p2so,p3s,p4s,p1p2s,p1p2so,p2p3s,ps,pso
  INTEGER i,j
  COMPLEX*16 Sum1, Sum2, Sum3, AA(8)

  A=0.d0
  AA=0.d0

  Sum1=0.d0
  Sum2=0.d0
  Sum3=0.d0

  facv=(-3.D0/2.D0)*(4.D0/3.D0)
  face=(1.D0/2.D0)*(4.D0/3.D0)   

!  facv=1.D0
!  face=1.D0

  nuclm=1.d0   !* it's assumed that all is given in nucleon-mass units
!  deltam=1.31d0
  deltam = 1.3130548775939D0 !MD = 1232
  IF(dmloop.GT.0.d0) deltam = dmloop   ! [v2] = bary(3)%rm (1232 MeV/M_N by default)
  MD=deltam
  x=(s-1.d0)/2.d0
  DeltaIn = getdelta()                  ! [v2] restored at exit

  DeltaAux=1.d0

  !* Graph 1:
  !*  Box graph with photons attached to the pions
  !* Isospin factor = 4/3

  p1s=nuclm**2
  p2s=0.0d0
  p3s=0.0d0
  p4s=nuclm**2
  p2p3s=t
  p1p2s=s

  m1s=deltam**2
  m4s=mpi2
  m2s=mpi2
  m3s=mpi2

  CALL Setdelta(Delta)

  D00=D0i(dd00,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11=D0i(dd11,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D12=D0i(dd12,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D13=D0i(dd13,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D22=D0i(dd22,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D23=D0i(dd23,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D33=D0i(dd33,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D111=D0i(dd111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D112=D0i(dd112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D113=D0i(dd113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D122=D0i(dd122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D123=D0i(dd123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D133=D0i(dd133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D222=D0i(dd222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D223=D0i(dd223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D233=D0i(dd233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D333=D0i(dd333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D001=D0i(dd001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D002=D0i(dd002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D003=D0i(dd003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D0000=D0i(dd0000,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D0011=D0i(dd0011,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D0012=D0i(dd0012,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D0013=D0i(dd0013,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D0022=D0i(dd0022,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D0023=D0i(dd0023,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D0033=D0i(dd0033,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1111=D0i(dd1111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1112=D0i(dd1112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1113=D0i(dd1113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1122=D0i(dd1122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1123=D0i(dd1123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1133=D0i(dd1133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1222=D0i(dd1222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1223=D0i(dd1223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1233=D0i(dd1233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D1333=D0i(dd1333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D2222=D0i(dd2222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D2223=D0i(dd2223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D2233=D0i(dd2233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D2333=D0i(dd2333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D3333=D0i(dd3333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00001=D0i(dd00001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00002=D0i(dd00002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00003=D0i(dd00003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00111=D0i(dd00111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00112=D0i(dd00112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00113=D0i(dd00113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00122=D0i(dd00122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00123=D0i(dd00123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00133=D0i(dd00133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00222=D0i(dd00222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00223=D0i(dd00223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00233=D0i(dd00233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00333=D0i(dd00333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11111=D0i(dd11111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11112=D0i(dd11112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11113=D0i(dd11113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11122=D0i(dd11122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11123=D0i(dd11123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11133=D0i(dd11133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11222=D0i(dd11222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11223=D0i(dd11223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11233=D0i(dd11233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D11333=D0i(dd11333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D12222=D0i(dd12222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D12223=D0i(dd12223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D12233=D0i(dd12233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D12333=D0i(dd12333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D13333=D0i(dd13333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D22222=D0i(dd22222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D22223=D0i(dd22223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D22233=D0i(dd22233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D22333=D0i(dd22333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D23333=D0i(dd23333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D33333=D0i(dd33333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)


  if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
  else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
  endif
  
  D0000D=D0i(dd0000,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00001D=D0i(dd00001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00002D=D0i(dd00002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
  D00003D=D0i(dd00003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)

  DED0000=2.D0*(D0000-D0000D)/sw1
  DED00001=2.D0*(D00001-D00001D)/sw1
  DED00002=2.D0*(D00002-D00002D)/sw1
  DED00003=2.D0*(D00003-D00003D)/sw1

  !* Graph 2:
  !*  Box graphs with one photon attached to the pion and the other attached to the piNDelta vertex
  !* Isospin factor = 4/3

  p1s=nuclm**2
  p2s=0.d0
  p1p2s=s

  m1s=deltam**2
  m2s=mpi2
  m3s=mpi2

  CALL Setdelta(Delta)

  C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C11=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C22=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C001=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C002=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C111=C0i(cc111,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C112=C0i(cc112,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C122=C0i(cc122,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C222=C0i(cc222,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C0000=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C0011=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C0012=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C0022=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C1111=C0i(cc1111,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C1112=C0i(cc1112,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C1122=C0i(cc1122,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C1222=C0i(cc1222,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C2222=C0i(cc2222,p1s,p2s,p1p2s,m1s,m2s,m3s)

  if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
  else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
  endif

  C00D=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C001D=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C002D=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C0000D=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C0011D=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C0012D=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C0022D=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)

  DEC00=2.D0*(C00-C00D)/sw1
  DEC001=2.D0*(C001-C001D)/sw1
  DEC002=2.D0*(C002-C002D)/sw1
  DEC0000=2.D0*(C0000-C0000D)/sw1
  DEC0011=2.D0*(C0011-C0011D)/sw1
  DEC0012=2.D0*(C0012-C0012D)/sw1
  DEC0022=2.D0*(C0022-C0022D)/sw1

  !write(*,*) '2', Sum2

  !* Graph 3:
  !*  Box graph with photons attached to the piNDelta vertices
  !* Isospin factor = 4/3

  p1s=s

  m1s=deltam**2
  m2s=mpi2

  CALL Setdelta(Delta)

  Bf00=B0i(bb00,p1s,m1s,m2s)
  Bf11=B0i(bb11,p1s,m1s,m2s)
  Bf001=B0i(bb001,p1s,m1s,m2s)
  Bf111=B0i(bb111,p1s,m1s,m2s)

  if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
  else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
  endif

  Bf00D=B0i(bb00,p1s,m1s,m2s)
  Bf001D=B0i(bb001,p1s,m1s,m2s)
  Bf11D=B0i(bb11,p1s,m1s,m2s)
  Bf111D=B0i(bb111,p1s,m1s,m2s)

  DEBf00=2.D0*(Bf00-Bf00D)/sw1
  DEBf001=2.D0*(Bf001-Bf001D)/sw1
  DEBf11=2.D0*(Bf11-Bf11D)/sw1
  DEBf111=2.D0*(Bf111-Bf111D)/sw1


  AA=A!*3.D0/4.D0
  A=0.d0

  !* Graph 4:
  !*  Reducible graphs with gammaNN vertices
  !* Isospin factor = -2/3
  !*
  !* [v2, 2026-10] The on-shell (p^2=M^2, real photon) subtraction constants of this
  !* graph are kinematics-independent. v1 evaluated the C-functions at
  !* p1p2s = M^2 + 1e-4 (one-sided offset, needed because the LoopTools reduction is
  !* singular at q=0), which leaves an O(1e-4) UV-dependent remainder and violates the
  !* Thomson limit (A1(0,0) = 1.4e-5, audit B1). They are now obtained from the
  !* symmetric average of p1p2s = M^2 +/- eps (error O(eps^2)) and cached
  !* (SUBROUTINE piD_onshell_consts below), which also saves ~40 LoopTools calls
  !* and two cache flushes per call.
  CALL piD_onshell_consts(mpi2,deltam,Delta, &
       Bf00o,Bf001o,DEBf00o,DEBf001o, &
       C00o,C001o,C002o,C0000o,C0011o,C0012o,C0022o, &
       DEC00o,DEC001o,DEC002o,DEC0000o,DEC0011o,DEC0012o,DEC0022o)

  !* Graph 5:
  !*  Self-energy graph
  !* Isospin factor = 2

  !Summing up

  !Box

  A(1) = A(1) &
       +DED00002*(44.D0/27.D0)+DED00001*(88.D0/27.D0)+D00002*(-40.D0/9.D0)+D00001*(-80.D0/9.D0)&
       +MD*DED0000*(-44.D0/27.D0)+MD*D0000*(40.D0/9.D0)+t*DED00002*(-13.D0/27.D0)+t*DED00001*(-16.D0/9.D0)&
       +t*D00222*(4.D0/9.D0)+t*D00002*(14.D0/9.D0)+t*D00123*(16.D0/9.D0)+t*D00122*(16.D0/9.D0)&
       +t*D00113*(16.D0/9.D0)+t*D00112*(8.D0/9.D0)+t*D00001*(16.D0/3.D0)+t*MD*DED0000*(20.D0/27.D0)&
       +t*MD*D0022*(-4.D0/9.D0)+t*MD*D0013*(-8.D0/9.D0)+t*MD*D0012*(-8.D0/9.D0)+t*MD*D0000*(-28.D0/9.D0)&
       +t**2*DED00002*(-4.D0/27.D0)+t**2*D00222*(1.D0/6.D0)+t**2*D00002*(11.D0/9.D0)+t**2*D00123*(1.D0/9.D0)&
       +t**2*D00122*(4.D0/9.D0)+t**2*D00113*(-4.D0/9.D0)+t**2*D00112*(1.D0/9.D0)+t**2*MD*D0013*(2.D0/9.D0)&
       +t**2*MD*D0012*(-1.D0/9.D0)+t**3*D00123*(-1.D0/9.D0)+x*DED00002*(44.D0/27.D0)+x*D00002*(-40.D0/9.D0)&
       +x*t*DED00002*(-16.D0/27.D0)+x*t*D00222*(14.D0/9.D0)+x*t*D00002*(44.D0/9.D0)+x*t*D00123*(20.D0/9.D0)&
       +x*t*D00122*(40.D0/9.D0)+x*t*D00112*(4.D0/3.D0)+x*t*MD*D0022*(-10.D0/9.D0)+x*t*MD*D0012*(-4.D0/3.D0)&
       +x*t**2*D00222*(1.D0/3.D0)+x*t**2*D00123*(-4.D0/9.D0)+x*t**2*D00122*(4.D0/9.D0)+x**2*D00222*(8.D0/9.D0)&
       +x**2*D00122*(16.D0/9.D0)+x**2*MD*D0022*(-8.D0/9.D0)+x**2*t*D00222*(14.D0/9.D0)+x**2*t*D00122*(16.D0/9.D0)&
       +x**3*D00222*(8.D0/9.D0)

  A(2) = A(2) &
       +DED00002*(2.D0/27.D0)+DED00001*(4.D0/27.D0)+D00222*(14.D0/9.D0)+D0022*(20.D0/9.D0)+D002*(2.D0/3.D0)&
       +D00002*(-4.D0/9.D0)+D00123*(28.D0/9.D0)+D00122*(28.D0/9.D0)+D0012*(40.D0/9.D0)+D00113*(28.D0/9.D0)&
       +D00112*(-28.D0/9.D0)+D00111*(-28.D0/9.D0)+D001*(4.D0/3.D0)+D00001*(-8.D0/9.D0)+MD*DED0000*(2.D0/27.D0)&
       +MD*D0022*(-14.D0/9.D0)+MD*D002*(-20.D0/9.D0)+MD*D0013*(-28.D0/9.D0)+MD*D0011*(28.D0/9.D0)&
       +MD*D00*(-2.D0/3.D0)+MD*D0000*(8.D0/9.D0)+t*DED00002*(1.D0/54.D0)+t*D22222*(-1.D0/9.D0)+t*D2222*(-2.D0/9.D0)&
       +t*D00222*(-1.D0/2.D0)+t*D222*(-1.D0/9.D0)+t*D0022*(-7.D0/9.D0)+t*D002*(-5.D0/18.D0)+t*D00002*(-1.D0/9.D0)&
       +t*D00123*(-5.D0/9.D0)+t*D123*(-4.D0/9.D0)+t*D1223*(-8.D0/9.D0)+t*D12223*(-2.D0/3.D0)+t*D12222*(-4.D0/9.D0)&
       +t*D1222*(-8.D0/9.D0)+t*D00122*(-16.D0/9.D0)+t*D122*(-4.D0/9.D0)+t*D0012*(-8.D0/3.D0)+t*D00113*(-4.D0/3.D0)&
       +t*D113*(-4.D0/9.D0)+t*D11233*(-2.D0/3.D0)+t*D1123*(-8.D0/9.D0)+t*D11223*(-8.D0/9.D0)+t*D1122*(-4.D0/9.D0)&
       +t*D00112*(5.D0/9.D0)+t*D112*(-2.D0/9.D0)+t*D11133*(-4.D0/9.D0)+t*D11123*(4.D0/9.D0)+t*D11122*(4.D0/9.D0)&
       +t*D11113*(4.D0/9.D0)+t*D11112*(2.D0/9.D0)+t*D00111*(4.D0/3.D0)+t*D001*(-8.D0/9.D0)+t*MD*D2222*(1.D0/9.D0)&
       +t*MD*D222*(2.D0/9.D0)+t*MD*D0022*(10.D0/9.D0)+t*MD*D22*(1.D0/9.D0)+t*MD*D002*(14.D0/9.D0)+t*MD*D0013*(4.D0/3.D0)&
       +t*MD*D13*(2.D0/9.D0)+t*MD*D123*(4.D0/9.D0)+t*MD*D1223*(4.D0/9.D0)+t*MD*D1222*(2.D0/9.D0)+t*MD*D122*(4.D0/9.D0)&
       +t*MD*D12*(2.D0/9.D0)+t*MD*D1133*(4.D0/9.D0)+t*MD*D1123*(2.D0/9.D0)+t*MD*D1122*(-2.D0/9.D0)+t*MD*D1113*(-4.D0/9.D0)&
       +t*MD*D1112*(-2.D0/9.D0)+t*MD*D0011*(-4.D0/3.D0)+t*MD*D00*(4.D0/9.D0)+t**2*D22222*(-1.D0/24.D0)&
       +t**2*D2222*(-1.D0/12.D0)+t**2*D00222*(-7.D0/18.D0)+t**2*D222*(-1.D0/24.D0)+t**2*D0022*(-11.D0/18.D0)&
       +t**2*D002*(-2.D0/9.D0)+t**2*D00123*(-2.D0/3.D0)+t**2*D123*(-1.D0/36.D0)+t**2*D1223*(-1.D0/18.D0)&
       +t**2*D12223*(-1.D0/9.D0)+t**2*D12222*(-1.D0/9.D0)+t**2*D1222*(-2.D0/9.D0)+t**2*D122*(-1.D0/9.D0)&
       +t**2*D113*(1.D0/9.D0)+t**2*D11233*(-1.D0/36.D0)+t**2*D1123*(2.D0/9.D0)+t**2*D11222*(1.D0/18.D0)&
       +t**2*D1122*(-1.D0/18.D0)+t**2*D00112*(2.D0/3.D0)+t**2*D112*(-1.D0/36.D0)+t**2*D11133*(1.D0/9.D0)&
       +t**2*D11122*(1.D0/9.D0)+t**2*D11113*(-1.D0/9.D0)+t**2*D11112*(1.D0/36.D0)+t**2*MD*D13*(-1.D0/18.D0)&
       +t**2*MD*D123*(-1.D0/9.D0)+t**2*MD*D1223*(-1.D0/18.D0)+t**2*MD*D1222*(1.D0/36.D0)+t**2*MD*D122*(1.D0/18.D0)&
       +t**2*MD*D12*(1.D0/36.D0)+t**2*MD*D1133*(-1.D0/9.D0)+t**2*MD*D1123*(1.D0/36.D0)+t**2*MD*D1113*(1.D0/9.D0)&
       +t**2*MD*D1112*(-1.D0/36.D0)+t**3*D123*(1.D0/36.D0)+t**3*D1223*(1.D0/18.D0)+t**3*D12223*(1.D0/36.D0)&
       +t**3*D11233*(1.D0/18.D0)+t**3*D11123*(-1.D0/18.D0)+x*DED00002*(2.D0/27.D0)+x*D00222*(14.D0/9.D0)&
       +x*D0022*(20.D0/9.D0)+x*D002*(2.D0/3.D0)+x*D00002*(-4.D0/9.D0)+x*D00123*(28.D0/9.D0)&
       +x*D00112*(-28.D0/9.D0)+x*t*D22222*(-7.D0/18.D0)+x*t*D2222*(-7.D0/9.D0)+x*t*D00222*(-14.D0/9.D0)&
       +x*t*D222*(-7.D0/18.D0)+x*t*D0022*(-22.D0/9.D0)+x*t*D002*(-8.D0/9.D0)+x*t*D00123*(-8.D0/3.D0)&
       +x*t*D123*(-5.D0/9.D0)+x*t*D1223*(-10.D0/9.D0)+x*t*D12223*(-4.D0/3.D0)+x*t*D12222*(-10.D0/9.D0)&
       +x*t*D1222*(-20.D0/9.D0)+x*t*D122*(-10.D0/9.D0)+x*t*D11233*(-7.D0/9.D0)+x*t*D11223*(-10.D0/9.D0)&
       +x*t*D11222*(4.D0/9.D0)+x*t*D1122*(-2.D0/3.D0)+x*t*D00112*(8.D0/3.D0)+x*t*D112*(-1.D0/3.D0)&
       +x*t*D11123*(4.D0/9.D0)+x*t*D11122*(10.D0/9.D0)+x*t*D11112*(1.D0/3.D0)+x*t*MD*D2222*(5.D0/18.D0)&
       +x*t*MD*D222*(5.D0/9.D0)+x*t*MD*D22*(5.D0/18.D0)+x*t*MD*D1223*(5.D0/9.D0)+x*t*MD*D1222*(1.D0/3.D0)&
       +x*t*MD*D122*(2.D0/3.D0)+x*t*MD*D12*(1.D0/3.D0)+x*t*MD*D1123*(1.D0/3.D0)+x*t*MD*D1122*(-5.D0/9.D0)&
       +x*t*MD*D1112*(-1.D0/3.D0)+x*t**2*D22222*(-1.D0/12.D0)+x*t**2*D2222*(-1.D0/6.D0)+x*t**2*D222*(-1.D0/12.D0)&
       +x*t**2*D123*(1.D0/9.D0)+x*t**2*D1223*(2.D0/9.D0)+x*t**2*D12223*(-1.D0/18.D0)+x*t**2*D12222*(-1.D0/9.D0)&
       +x*t**2*D1222*(-2.D0/9.D0)+x*t**2*D122*(-1.D0/9.D0)+x*t**2*D11233*(2.D0/9.D0)+x*t**2*D11223*(-1.D0/9.D0)&
       +x*t**2*D11222*(1.D0/6.D0)+x*t**2*D11123*(-2.D0/9.D0)+x*t**2*D11122*(1.D0/9.D0)+x**2*D22222*(-2.D0/9.D0)&
       +x**2*D2222*(-4.D0/9.D0)+x**2*D222*(-2.D0/9.D0)+x**2*D12223*(-4.D0/9.D0)+x**2*D12222*(-4.D0/9.D0)&
       +x**2*D1222*(-8.D0/9.D0)+x**2*D122*(-4.D0/9.D0)+x**2*D11223*(-4.D0/9.D0)+x**2*D11222*(4.D0/9.D0)&
       +x**2*D11122*(4.D0/9.D0)+x**2*MD*D2222*(2.D0/9.D0)+x**2*MD*D222*(4.D0/9.D0)+x**2*MD*D22*(2.D0/9.D0)&
       +x**2*MD*D1223*(4.D0/9.D0)+x**2*MD*D1122*(-4.D0/9.D0)+x**2*t*D22222*(-7.D0/18.D0)+x**2*t*D2222*(-7.D0/9.D0)&
       +x**2*t*D222*(-7.D0/18.D0)+x**2*t*D12223*(-7.D0/9.D0)+x**2*t*D12222*(-4.D0/9.D0)+x**2*t*D1222*(-8.D0/9.D0)&
       +x**2*t*D122*(-4.D0/9.D0)+x**2*t*D11223*(-4.D0/9.D0)+x**2*t*D11222*(7.D0/9.D0)+x**2*t*D11122*(4.D0/9.D0)&
       +x**3*D22222*(-2.D0/9.D0)+x**3*D2222*(-4.D0/9.D0)+x**3*D222*(-2.D0/9.D0)+x**3*D12223*(-4.D0/9.D0)&
       +x**3*D11222*(4.D0/9.D0)

  A(3) = A(3)

  A(4) = A(4) &
       +DED00002*(2.D0/9.D0)+D00222*(-4.D0/9.D0)+D00002*(-20.D0/3.D0)+D00123*(-8.D0/9.D0)+D00122*(-16.D0/9.D0)&
       +D00112*(-8.D0/9.D0)+MD*D0022*(4.D0/9.D0)+MD*D0012*(8.D0/9.D0)+t*DED00002*(-8.D0/27.D0)+t*D00222*(1.D0/3.D0)&
       +t*D00002*(22.D0/9.D0)+t*D00123*(10.D0/9.D0)+t*D00122*(8.D0/9.D0)+t*D00112*(2.D0/9.D0)+t*MD*D0012*(-2.D0/9.D0)&
       +t**2*D00123*(-2.D0/9.D0)+x*D00222*(-8.D0/9.D0)+x*D00122*(-16.D0/9.D0)+x*MD*D0022*(4.D0/9.D0)+x*t*D00222*(2.D0/3.D0)&
       +x*t*D00122*(8.D0/9.D0)+x**2*D00222*(4.D0/9.D0)

  A(5) = A(5) &
       +DED00002*(1.D0/9.D0)+DED0000*(5.D0/27.D0)+D00222*(-2.D0/9.D0)+D0022*(-2.D0/9.D0)+D00002*(-10.D0/3.D0)&
       +D0013*(-4.D0/9.D0)+D00123*(-4.D0/9.D0)+D00122*(-8.D0/9.D0)+D0012*(-8.D0/9.D0)+D00112*(-4.D0/9.D0)&
       +D0011*(-4.D0/9.D0)+D0000*(-22.D0/9.D0)+MD*D0022*(2.D0/9.D0)+MD*D002*(2.D0/9.D0)+MD*D0012*(4.D0/9.D0)&
       +MD*D001*(4.D0/9.D0)+t*DED00002*(-4.D0/27.D0)+t*DED0000*(-4.D0/27.D0)+t*D00222*(1.D0/6.D0)+t*D0022*(1.D0/6.D0)&
       +t*D00002*(11.D0/9.D0)+t*D0013*(5.D0/9.D0)+t*D00123*(5.D0/9.D0)+t*D00122*(4.D0/9.D0)+t*D0012*(4.D0/9.D0)&
       +t*D00112*(1.D0/9.D0)+t*D0011*(1.D0/9.D0)+t*D0000*(8.D0/9.D0)+t*MD*D0012*(-1.D0/9.D0)+t*MD*D001*(-1.D0/9.D0)&
       +t**2*D0013*(-1.D0/9.D0)+t**2*D00123*(-1.D0/9.D0)+x*D00222*(-4.D0/9.D0)+x*D0022*(-4.D0/9.D0)&
       +x*D00122*(-8.D0/9.D0)+x*D0012*(-8.D0/9.D0)+x*MD*D0022*(2.D0/9.D0)+x*MD*D002*(2.D0/9.D0)&
       +x*t*D00222*(1.D0/3.D0)+x*t*D0022*(1.D0/3.D0)+x*t*D00122*(4.D0/9.D0)+x*t*D0012*(4.D0/9.D0)&
       +x**2*D00222*(2.D0/9.D0)+x**2*D0022*(2.D0/9.D0)

  A(6) = A(6) &
       +DED00002*(-1.D0/9.D0)+DED0000*(-5.D0/27.D0)+D00222*(2.D0/9.D0)+D0022*(2.D0/9.D0)+D00002*(10.D0/3.D0)&
       +D0013*(4.D0/9.D0)+D00123*(4.D0/9.D0)+D00122*(8.D0/9.D0)+D0012*(8.D0/9.D0)+D00112*(4.D0/9.D0)+D0011*(4.D0/9.D0)&
       +D0000*(22.D0/9.D0)+MD*D0022*(-2.D0/9.D0)+MD*D002*(-2.D0/9.D0)+MD*D0012*(-4.D0/9.D0)+MD*D001*(-4.D0/9.D0)&
       +t*DED00002*(4.D0/27.D0)+t*DED0000*(4.D0/27.D0)+t*D00222*(-1.D0/6.D0)+t*D0022*(-1.D0/6.D0)+t*D00002*(-11.D0/9.D0)&
       +t*D0013*(-5.D0/9.D0)+t*D00123*(-5.D0/9.D0)+t*D00122*(-4.D0/9.D0)+t*D0012*(-4.D0/9.D0)+t*D00112*(-1.D0/9.D0)&
       +t*D0011*(-1.D0/9.D0)+t*D0000*(-8.D0/9.D0)+t*MD*D0012*(1.D0/9.D0)+t*MD*D001*(1.D0/9.D0)+t**2*D0013*(1.D0/9.D0)&
       +t**2*D00123*(1.D0/9.D0)+x*D00222*(4.D0/9.D0)+x*D0022*(4.D0/9.D0)+x*D00122*(8.D0/9.D0)+x*D0012*(8.D0/9.D0)&
       +x*MD*D0022*(-2.D0/9.D0)+x*MD*D002*(-2.D0/9.D0)+x*t*D00222*(-1.D0/3.D0)+x*t*D0022*(-1.D0/3.D0)+x*t*D00122*(-4.D0/9.D0)&
       +x*t*D0012*(-4.D0/9.D0)+x**2*D00222*(-2.D0/9.D0)+x**2*D0022*(-2.D0/9.D0)

  A(7) = A(7) &
       +DED00002*(-1.D0/27.D0)+D22222*(-1.D0/9.D0)+D2222*(-2.D0/9.D0)+D00222*(-19.D0/9.D0)+D222*(-1.D0/9.D0)&
       +D0022*(-10.D0/3.D0)+D002*(-11.D0/9.D0)+D00002*(2.D0/9.D0)+D00123*(-38.D0/9.D0)+D123*(-2.D0/9.D0)&
       +D1223*(-4.D0/9.D0)+D12223*(-4.D0/9.D0)+D12222*(-4.D0/9.D0)+D1222*(-8.D0/9.D0)+D122*(-4.D0/9.D0)&
       +D11233*(-2.D0/9.D0)+D11223*(-4.D0/9.D0)+D1122*(-4.D0/9.D0)+D00112*(38.D0/9.D0)+D112*(-2.D0/9.D0)&
       +D11122*(4.D0/9.D0)+D11112*(2.D0/9.D0)+MD*D2222*(1.D0/9.D0)+MD*D222*(2.D0/9.D0)+MD*D22*(1.D0/9.D0)&
       +MD*D1223*(2.D0/9.D0)+MD*D1222*(2.D0/9.D0)+MD*D122*(4.D0/9.D0)+MD*D12*(2.D0/9.D0)+MD*D1123*(2.D0/9.D0)&
       +MD*D1122*(-2.D0/9.D0)+MD*D1112*(-2.D0/9.D0)+t*D22222*(1.D0/12.D0)+t*D2222*(1.D0/6.D0)+t*D00222*(7.D0/9.D0)&
       +t*D222*(1.D0/12.D0)+t*D0022*(11.D0/9.D0)+t*D002*(4.D0/9.D0)+t*D00123*(4.D0/3.D0)+t*D123*(5.D0/18.D0)&
       +t*D1223*(5.D0/9.D0)+t*D12223*(4.D0/9.D0)+t*D12222*(2.D0/9.D0)+t*D1222*(4.D0/9.D0)+t*D122*(2.D0/9.D0)&
       +t*D11233*(1.D0/2.D0)+t*D11223*(2.D0/9.D0)+t*D11222*(-1.D0/9.D0)+t*D1122*(1.D0/9.D0)+t*D00112*(-4.D0/3.D0)&
       +t*D112*(1.D0/18.D0)+t*D11123*(-4.D0/9.D0)+t*D11122*(-2.D0/9.D0)+t*D11112*(-1.D0/18.D0)+t*MD*D1222*(-1.D0/18.D0)&
       +t*MD*D122*(-1.D0/9.D0)+t*MD*D12*(-1.D0/18.D0)+t*MD*D1123*(-1.D0/18.D0)+t*MD*D1112*(1.D0/18.D0)+t**2*D123*(-1.D0/18.D0)&
       +t**2*D1223*(-1.D0/9.D0)+t**2*D12223*(-1.D0/18.D0)+t**2*D11233*(-1.D0/9.D0)+t**2*D11123*(1.D0/9.D0)+x*D22222*(-2.D0/9.D0)&
       +x*D2222*(-4.D0/9.D0)+x*D222*(-2.D0/9.D0)+x*D12223*(-4.D0/9.D0)+x*D12222*(-4.D0/9.D0)+x*D1222*(-8.D0/9.D0)&
       +x*D122*(-4.D0/9.D0)+x*D11223*(-4.D0/9.D0)+x*D11222*(4.D0/9.D0)+x*D11122*(4.D0/9.D0)+x*MD*D2222*(1.D0/9.D0)&
       +x*MD*D222*(2.D0/9.D0)+x*MD*D22*(1.D0/9.D0)+x*MD*D1223*(2.D0/9.D0)+x*MD*D1122*(-2.D0/9.D0)+x*t*D22222*(1.D0/6.D0)&
       +x*t*D2222*(1.D0/3.D0)+x*t*D222*(1.D0/6.D0)+x*t*D12223*(1.D0/3.D0)+x*t*D12222*(2.D0/9.D0)+x*t*D1222*(4.D0/9.D0)&
       +x*t*D122*(2.D0/9.D0)+x*t*D11223*(2.D0/9.D0)+x*t*D11222*(-1.D0/3.D0)+x*t*D11122*(-2.D0/9.D0)+x**2*D22222*(1.D0/9.D0)&
       +x**2*D2222*(2.D0/9.D0)+x**2*D222*(1.D0/9.D0)+x**2*D12223*(2.D0/9.D0)+x**2*D11222*(-2.D0/9.D0)

  A(8) = A(8)

  !KR 

  A(1) = A(1) &
       +DEBf111*(1.D0/27.D0)+DEBf001*(11.D0/27.D0)+Bf111*(-2.D0/9.D0)+Bf001*(-10.D0/9.D0)+MD*DEBf11*(-1.D0/27.D0)&
       +MD*DEBf00*(-1.D0/3.D0)+MD*Bf11*(2.D0/9.D0)+MD*Bf00*(2.D0/3.D0)+t*DEBf111*(1.D0/108.D0)+t*DEBf001*(11.D0/108.D0)&
       +t*Bf111*(-1.D0/18.D0)+t*Bf001*(-5.D0/18.D0)+x*DEBf111*(1.D0/9.D0)+x*DEBf001*(11.D0/27.D0)+x*Bf111*(-2.D0/3.D0)&
       +x*Bf001*(-10.D0/9.D0)+x*MD*DEBf11*(-2.D0/27.D0)+x*MD*Bf11*(4.D0/9.D0)+x*t*DEBf111*(1.D0/54.D0)+x*t*Bf111*(-1.D0/9.D0)&
       +x**2*DEBf111*(2.D0/27.D0)+x**2*Bf111*(-4.D0/9.D0)

  A(2) = A(2) &
       +DEBf111*(1.D0/108.D0)+Bf111*(-1.D0/18.D0)+MD*DEBf11*(-1.D0/36.D0)+t*DEBf111*(1.D0/432.D0)&
       +t*Bf111*(-1.D0/72.D0)+x*DEBf111*(1.D0/108.D0)+x*Bf111*(-1.D0/18.D0)

  A(3) = A(3) &
       +DEBf111*(1.D0/27.D0)+DEBf001*(1.D0/9.D0)+Bf111*(1.D0/9.D0)+Bf001*(2.D0/3.D0)+MD*DEBf11*(-1.D0/27.D0)&
       +MD*DEBf00*(1.D0/27.D0)+MD*Bf11*(-1.D0/9.D0)+MD*Bf00*(-2.D0/9.D0)+t*DEBf111*(1.D0/108.D0)+t*DEBf001*(1.D0/36.D0)&
       +t*Bf111*(1.D0/36.D0)+t*Bf001*(1.D0/6.D0)+x*DEBf111*(1.D0/9.D0)+x*DEBf001*(1.D0/9.D0)+x*Bf111*(1.D0/3.D0)&
       +x*Bf001*(2.D0/3.D0)+x*MD*DEBf11*(-2.D0/27.D0)+x*MD*Bf11*(-2.D0/9.D0)+x*t*DEBf111*(1.D0/54.D0)+x*t*Bf111*(1.D0/18.D0)&
       +x**2*DEBf111*(2.D0/27.D0)+x**2*Bf111*(2.D0/9.D0)

  A(4) = A(4) &
       +DEBf111*(1.D0/54.D0)+DEBf001*(11.D0/54.D0)+Bf111*(-1.D0/9.D0)+Bf001*(-5.D0/9.D0)+x*DEBf111*(1.D0/27.D0)+x*Bf111*(-2.D0/9.D0)

  A(5) = A(5) &
       +DEBf111*(-1.D0/108.D0)+DEBf001*(-1.D0/27.D0)+Bf111*(-1.D0/36.D0)+Bf001*(-1.D0/9.D0)+MD*DEBf11*(-1.D0/108.D0)&
       +MD*Bf11*(-1.D0/36.D0)+x*DEBf111*(-1.D0/54.D0)+x*Bf111*(-1.D0/18.D0)

  A(6) = A(6) &
       +DEBf111*(1.D0/108.D0)+DEBf001*(1.D0/27.D0)+Bf111*(1.D0/36.D0)+Bf001*(1.D0/9.D0)+MD*DEBf11*(-1.D0/108.D0)&
       +MD*Bf11*(-1.D0/36.D0)+x*DEBf111*(1.D0/54.D0)+x*Bf111*(1.D0/18.D0)

  A(7) = A(7) &
       +DEBf111*(-1.D0/216.D0)+Bf111*(1.D0/36.D0)

  A(8) = A(8) &
       +DEBf111*(-1.D0/54.D0)+DEBf001*(-1.D0/18.D0)+Bf111*(-1.D0/18.D0)+Bf001*(-1.D0/3.D0)+x*DEBf111*(-1.D0/27.D0)+x*Bf111*(-1.D0/9.D0)

  !KP

  A(1) = A(1) &
       +DEC0022*(4.D0/27.D0)+DEC0012*(8.D0/27.D0)+DEC0011*(4.D0/27.D0)+DEC0000*(4.D0/27.D0)+C0022*(-8.D0/9.D0)&
       +C0012*(-16.D0/9.D0)+C0011*(-8.D0/9.D0)+C0000*(-8.D0/9.D0)+MD*DEC002*(-4.D0/27.D0)+MD*DEC001*(-4.D0/27.D0)&
       +MD*C002*(8.D0/9.D0)+MD*C001*(8.D0/9.D0)+t*DEC0022*(1.D0/27.D0)+t*DEC0012*(1.D0/27.D0)+t*C0022*(-2.D0/9.D0)&
       +t*C0012*(-2.D0/9.D0)+t*MD*DEC002*(-1.D0/27.D0)+t*MD*C002*(-1.D0/9.D0)+x*DEC0022*(8.D0/27.D0)&
       +x*DEC0012*(8.D0/27.D0)+x*C0022*(-16.D0/9.D0)+x*C0012*(-16.D0/9.D0)+x*MD*DEC002*(-4.D0/27.D0)&
       +x*MD*C002*(8.D0/9.D0)+x*t*DEC0022*(1.D0/27.D0)+x*t*C0022*(-2.D0/9.D0)+x**2*DEC0022*(4.D0/27.D0)&
       +x**2*C0022*(-8.D0/9.D0)

  A(2) = A(2) &
       +DEC0022*(-5.D0/9.D0)+DEC002*(-4.D0/9.D0)+DEC0012*(-26.D0/27.D0)+DEC0011*(-11.D0/27.D0)+DEC001*(-10.D0/27.D0)&
       +C2222*(4.D0/9.D0)+C222*(4.D0/9.D0)+C0022*(2.D0)+C002*(4.D0/3.D0)+C1222*(4.D0/3.D0)+C122*(8.D0/9.D0)&
       +C0012*(28.D0/9.D0)+C1122*(4.D0/3.D0)+C112*(4.D0/9.D0)+C1112*(4.D0/9.D0)+C0011*(10.D0/9.D0)&
       +C001*(8.D0/9.D0)+MD*DEC002*(11.D0/27.D0)+MD*DEC001*(11.D0/27.D0)+MD*DEC00*(10.D0/27.D0)+MD*C222*(-4.D0/9.D0)&
       +MD*C22*(-4.D0/9.D0)+MD*C002*(-16.D0/9.D0)+MD*C122*(-8.D0/9.D0)+MD*C12*(-4.D0/9.D0)+MD*C112*(-4.D0/9.D0)&
       +MD*C001*(-10.D0/9.D0)+MD*C00*(-8.D0/9.D0)+t*DEC0022*(-1.D0/12.D0)+t*DEC002*(-2.D0/27.D0)+t*DEC0012*(-7.D0/108.D0)&
       +t*C2222*(5.D0/36.D0)+t*C222*(5.D0/36.D0)+t*C0022*(2.D0/3.D0)+t*C002*(4.D0/9.D0)+t*C1222*(11.D0/36.D0)&
       +t*C122*(1.D0/6.D0)+t*C0012*(5.D0/9.D0)+t*C1122*(7.D0/36.D0)+t*C112*(1.D0/36.D0)+t*C1112*(1.D0/36.D0)&
       +t*MD*C222*(1.D0/36.D0)+t*MD*C22*(1.D0/36.D0)+t*MD*C12*(-1.D0/36.D0)+t*MD*C112*(-1.D0/36.D0)+x*DEC0022*(-1.D0/3.D0)&
       +x*DEC002*(-8.D0/27.D0)+x*DEC0012*(-7.D0/27.D0)+x*C2222*(11.D0/9.D0)+x*C222*(11.D0/9.D0)+x*C0022*(8.D0/3.D0)&
       +x*C002*(16.D0/9.D0)+x*C1222*(25.D0/9.D0)+x*C122*(14.D0/9.D0)+x*C0012*(20.D0/9.D0)+x*C1122*(17.D0/9.D0)&
       +x*C112*(1.D0/3.D0)+x*C1112*(1.D0/3.D0)+x*MD*C222*(-7.D0/9.D0)+x*MD*C22*(-7.D0/9.D0)+x*MD*C122*(-10.D0/9.D0)&
       +x*MD*C12*(-1.D0/3.D0)+x*MD*C112*(-1.D0/3.D0)+x*t*C2222*(2.D0/9.D0)+x*t*C222*(2.D0/9.D0)+x*t*C1222*(1.D0/3.D0)&
       +x*t*C122*(1.D0/9.D0)+x*t*C1122*(1.D0/9.D0)+x**2*C2222*(8.D0/9.D0)+x**2*C222*(8.D0/9.D0)+x**2*C1222*(4.D0/3.D0)&
       +x**2*C122*(4.D0/9.D0)+x**2*C1122*(4.D0/9.D0)

  A(3) = A(3) &
       +DEC0022*(4.D0/27.D0)+DEC0012*(8.D0/27.D0)+DEC0011*(4.D0/27.D0)+DEC0000*(4.D0/9.D0)+C0022*(4.D0/9.D0)&
       +C0012*(8.D0/9.D0)+C0011*(4.D0/9.D0)+C0000*(8.D0/3.D0)+MD*DEC002*(-4.D0/27.D0)+MD*DEC001*(-4.D0/27.D0)&
       +MD*C002*(-4.D0/9.D0)+MD*C001*(-4.D0/9.D0)+x*DEC0022*(8.D0/27.D0)+x*DEC0012*(8.D0/27.D0)+x*C0022*(8.D0/9.D0)&
       +x*C0012*(8.D0/9.D0)+x*MD*DEC002*(-4.D0/27.D0)+x*MD*C002*(-4.D0/9.D0)

  A(4) = A(4) &
       +DEC0022*(2.D0/27.D0)+DEC0012*(2.D0/27.D0)+C0022*(-4.D0/9.D0)+C0012*(-4.D0/9.D0)+MD*DEC002*(-2.D0/27.D0)&
       +MD*C002*(-2.D0/9.D0)+x*DEC0022*(2.D0/27.D0)+x*C0022*(-4.D0/9.D0)

  A(5) = A(5) &
       +DEC0022*(-1.D0/6.D0)+DEC002*(-5.D0/54.D0)+DEC0012*(-5.D0/27.D0)+DEC0011*(-1.D0/54.D0)+DEC001*(1.D0/54.D0)&
       +DEC0000*(1.D0/9.D0)+C2222*(-1.D0/9.D0)+C222*(-1.D0/9.D0)+C0022*(-4.D0/3.D0)+C002*(-7.D0/9.D0)+C1222*(-1.D0/3.D0)&
       +C122*(-2.D0/9.D0)+C0012*(-14.D0/9.D0)+C1122*(-1.D0/3.D0)+C112*(-1.D0/9.D0)+C1112*(-1.D0/9.D0)+C0011*(-2.D0/9.D0)&
       +C001*(-1.D0/9.D0)+C0000*(-4.D0/3.D0)+MD*DEC002*(1.D0/54.D0)+MD*DEC001*(1.D0/54.D0)+MD*DEC00*(-1.D0/54.D0)&
       +MD*C222*(1.D0/9.D0)+MD*C22*(1.D0/9.D0)+MD*C002*(2.D0/9.D0)+MD*C122*(2.D0/9.D0)+MD*C12*(1.D0/9.D0)&
       +MD*C112*(1.D0/9.D0)+MD*C001*(2.D0/9.D0)+MD*C00*(1.D0/9.D0)+x*DEC0022*(-5.D0/54.D0)+x*DEC002*(-1.D0/27.D0)&
       +x*DEC0012*(-7.D0/54.D0)+x*C2222*(-5.D0/18.D0)+x*C222*(-5.D0/18.D0)+x*C0022*(-10.D0/9.D0)+x*C002*(-4.D0/9.D0)&
       +x*C1222*(-11.D0/18.D0)+x*C122*(-1.D0/3.D0)+x*C0012*(-8.D0/9.D0)+x*C1122*(-7.D0/18.D0)+x*C112*(-1.D0/18.D0)&
       +x*C1112*(-1.D0/18.D0)+x*MD*C222*(1.D0/6.D0)+x*MD*C22*(1.D0/6.D0)+x*MD*C122*(2.D0/9.D0)+x*MD*C12*(1.D0/18.D0)&
       +x*MD*C112*(1.D0/18.D0)+x**2*C2222*(-1.D0/9.D0)+x**2*C222*(-1.D0/9.D0)+x**2*C1222*(-2.D0/9.D0)&
       +x**2*C122*(-1.D0/9.D0)+x**2*C1122*(-1.D0/9.D0)

  A(6) = A(6) &
       +DEC0022*(-1.D0/54.D0)+DEC002*(-1.D0/54.D0)+DEC0011*(1.D0/54.D0)+DEC001*(-1.D0/54.D0)+DEC0000*(-1.D0/9.D0)&
       +C0022*(4.D0/9.D0)+C002*(1.D0/9.D0)+C0012*(2.D0/3.D0)+C0011*(2.D0/9.D0)+C001*(1.D0/9.D0)+C0000*(4.D0/3.D0)&
       +MD*DEC002*(1.D0/54.D0)+MD*DEC001*(-1.D0/54.D0)+MD*DEC00*(1.D0/54.D0)+MD*C002*(-1.D0/9.D0)+MD*C001*(-2.D0/9.D0)&
       +MD*C00*(-1.D0/9.D0)+x*DEC0022*(5.D0/54.D0)+x*DEC002*(1.D0/27.D0)+x*DEC0012*(7.D0/54.D0)+x*C2222*(1.D0/18.D0)&
       +x*C222*(1.D0/18.D0)+x*C0022*(10.D0/9.D0)+x*C002*(4.D0/9.D0)+x*C1222*(1.D0/6.D0)+x*C122*(1.D0/9.D0)&
       +x*C0012*(8.D0/9.D0)+x*C1122*(1.D0/6.D0)+x*C112*(1.D0/18.D0)+x*C1112*(1.D0/18.D0)+x*MD*C222*(-1.D0/18.D0)&
       +x*MD*C22*(-1.D0/18.D0)+x*MD*C122*(-1.D0/9.D0)+x*MD*C12*(-1.D0/18.D0)+x*MD*C112*(-1.D0/18.D0)&
       +x**2*C2222*(1.D0/9.D0)+x**2*C222*(1.D0/9.D0)+x**2*C1222*(2.D0/9.D0)+x**2*C122*(1.D0/9.D0)+x**2*C1122*(1.D0/9.D0)

  A(7) = A(7) &
       +DEC0022*(1.D0/6.D0)+DEC002*(4.D0/27.D0)+DEC0012*(7.D0/54.D0)+C2222*(-5.D0/18.D0)+C222*(-5.D0/18.D0)&
       +C0022*(-4.D0/3.D0)+C002*(-8.D0/9.D0)+C1222*(-11.D0/18.D0)+C122*(-1.D0/3.D0)+C0012*(-10.D0/9.D0)&
       +C1122*(-7.D0/18.D0)+C112*(-1.D0/18.D0)+C1112*(-1.D0/18.D0)+MD*C222*(-1.D0/18.D0)+MD*C22*(-1.D0/18.D0)&
       +MD*C12*(1.D0/18.D0)+MD*C112*(1.D0/18.D0)+x*C2222*(-4.D0/9.D0)+x*C222*(-4.D0/9.D0)+x*C1222*(-2.D0/3.D0)&
       +x*C122*(-2.D0/9.D0)+x*C1122*(-2.D0/9.D0)

  A(8) = A(8)

  !VR

  A(1) = A(1) &
       +facv*( &
       +x**(-1)*t*DEBf001o*(-1.D0/216.D0)+x**(-1)*t*DEBf001*(1.D0/216.D0)+x**(-1)*t*Bf001o*(1.D0/36.D0)&
       +x**(-1)*t*Bf001*(-1.D0/36.D0)+x**(-1)*t*DEC0022o*(1.D0/54.D0)+x**(-1)*t*DEC0012o*(1.D0/27.D0)&
       +x**(-1)*t*DEC0011o*(1.D0/54.D0)+x**(-1)*t*DEC0000o*(-5.D0/108.D0)+x**(-1)*t*DEC0022*(-1.D0/54.D0)&
       +x**(-1)*t*DEC0012*(-1.D0/27.D0)+x**(-1)*t*DEC0011*(-1.D0/54.D0)+x**(-1)*t*DEC0000*(5.D0/108.D0)&
       +x**(-1)*t*C0022o*(1.D0/18.D0)+x**(-1)*t*C0012o*(1.D0/9.D0)+x**(-1)*t*C0011o*(1.D0/18.D0)&
       +x**(-1)*t*C0000o*(11.D0/18.D0)+x**(-1)*t*C0022*(-1.D0/18.D0)+x**(-1)*t*C0012*(-1.D0/9.D0)&
       +x**(-1)*t*C0011*(-1.D0/18.D0)+x**(-1)*t*C0000*(-11.D0/18.D0)+x**(-1)*t*MD*DEBf00o*(1.D0/216.D0)&
       +x**(-1)*t*MD*DEBf00*(-1.D0/216.D0)+x**(-1)*t*MD*Bf00o*(-1.D0/36.D0)+x**(-1)*t*MD*Bf00*(1.D0/36.D0)&
       +x**(-1)*t*MD*DEC002o*(-1.D0/54.D0)+x**(-1)*t*MD*DEC001o*(-1.D0/54.D0)+x**(-1)*t*MD*DEC002*(1.D0/54.D0)&
       +x**(-1)*t*MD*DEC001*(1.D0/54.D0)+x**(-1)*t*MD*C002o*(-1.D0/18.D0)+x**(-1)*t*MD*C001o*(-1.D0/18.D0)&
       +x**(-1)*t*MD*C002*(1.D0/18.D0)+x**(-1)*t*MD*C001*(1.D0/18.D0)&
       +DEBf001o*(19.D0/54.D0)+DEBf001*(1.D0/54.D0)+Bf001o*(-7.D0/9.D0)+Bf001*(-1.D0/9.D0)+DEC0022o*(-1.D0/3.D0)&
       +DEC002o*(-1.D0/3.D0)+DEC0012o*(-2.D0/3.D0)+DEC0011o*(-1.D0/3.D0)+DEC001o*(-1.D0/3.D0)+DEC0000o*(-1.D0/3.D0)&
       +DEC0000*(11.D0/27.D0)+C0022o*(2.D0/3.D0)+C002o*(2.D0/3.D0)+C0012o*(4.D0/3.D0)+C0011o*(2.D0/3.D0)&
       +C001o*(2.D0/3.D0)+C0000o*(2.D0/3.D0)+C0000*(-10.D0/9.D0)+MD*DEBf00o*(-19.D0/54.D0)+MD*DEBf00*(1.D0/54.D0)&
       +MD*Bf00o*(7.D0/9.D0)+MD*Bf00*(-1.D0/9.D0)+MD*DEC002o*(1.D0/3.D0)+MD*DEC001o*(1.D0/3.D0)&
       +MD*DEC00o*(1.D0/3.D0)+MD*C002o*(-2.D0/3.D0)+MD*C001o*(-2.D0/3.D0)+MD*C00o*(-2.D0/3.D0)+t*DEBf001*(1.D0/108.D0)&
       +t*Bf001*(-1.D0/18.D0)+t*DEC0022*(-5.D0/108.D0)+t*DEC0012*(-1.D0/18.D0)+t*DEC0011*(-1.D0/108.D0)&
       +t*DEC0000*(2.D0/27.D0)+t*C0022*(-5.D0/36.D0)+t*C0012*(-1.D0/6.D0)+t*C0011*(-1.D0/36.D0)&
       +t*C0000*(-4.D0/9.D0)+t*MD*DEC002*(1.D0/36.D0)+t*MD*DEC001*(1.D0/108.D0)+t*MD*C002*(1.D0/12.D0)&
       +t*MD*C001*(1.D0/36.D0)+x*DEBf001*(1.D0/27.D0)+x*Bf001*(-2.D0/9.D0)+x*DEC0022*(-1.D0/27.D0)&
       +x*DEC0012*(-2.D0/27.D0)+x*DEC0011*(-1.D0/27.D0)+x*DEC0000*(8.D0/27.D0)+x*C0022*(-1.D0/9.D0)&
       +x*C0012*(-2.D0/9.D0)+x*C0011*(-1.D0/9.D0)+x*C0000*(-16.D0/9.D0)+x*MD*DEC002*(1.D0/27.D0)&
       +x*MD*DEC001*(1.D0/27.D0)+x*MD*C002*(1.D0/9.D0)+x*MD*C001*(1.D0/9.D0)+x*t*DEC0022*(-1.D0/54.D0)&
       +x*t*DEC0012*(-1.D0/108.D0)+x*t*C0022*(-1.D0/18.D0)+x*t*C0012*(-1.D0/9.D0)+x**2*DEC0022*(-2.D0/27.D0)&
       +x**2*DEC0012*(-1.D0/27.D0)+x**2*C0022*(-2.D0/9.D0)+x**2*C0012*(-4.D0/9.D0)&
       )

  A(2) = A(2) &
       +facv*( &
       +x**(-1)*DEBf001o*(5.D0/54.D0)+x**(-1)*DEBf001*(-5.D0/54.D0)+x**(-1)*Bf001o*(-2.D0/9.D0)&
       +x**(-1)*Bf001*(2.D0/9.D0)+x**(-1)*DEC0022o*(-11.D0/54.D0)+x**(-1)*DEC002o*(-1.D0/6.D0)&
       +x**(-1)*DEC0012o*(-11.D0/27.D0)+x**(-1)*DEC0011o*(-11.D0/54.D0)+x**(-1)*DEC001o*(-1.D0/6.D0)&
       +x**(-1)*DEC0000o*(-2.D0/27.D0)+x**(-1)*DEC0022*(11.D0/54.D0)+x**(-1)*DEC002*(1.D0/6.D0)&
       +x**(-1)*DEC0012*(11.D0/27.D0)+x**(-1)*DEC0011*(11.D0/54.D0)+x**(-1)*DEC001*(1.D0/6.D0)&
       +x**(-1)*DEC0000*(2.D0/27.D0)+x**(-1)*C0022o*(2.D0/9.D0)+x**(-1)*C002o*(1.D0/3.D0)&
       +x**(-1)*C0012o*(4.D0/9.D0)+x**(-1)*C0011o*(2.D0/9.D0)+x**(-1)*C001o*(1.D0/3.D0)&
       +x**(-1)*C0000o*(-8.D0/9.D0)+x**(-1)*C0022*(-2.D0/9.D0)+x**(-1)*C002*(-1.D0/3.D0)&
       +x**(-1)*C0012*(-4.D0/9.D0)+x**(-1)*C0011*(-2.D0/9.D0)+x**(-1)*C001*(-1.D0/3.D0)&
       +x**(-1)*C0000*(8.D0/9.D0)+x**(-1)*MD*DEBf00o*(-5.D0/54.D0)+x**(-1)*MD*DEBf00*(5.D0/54.D0)&
       +x**(-1)*MD*Bf00o*(2.D0/9.D0)+x**(-1)*MD*Bf00*(-2.D0/9.D0)+x**(-1)*MD*DEC002o*(11.D0/54.D0)&
       +x**(-1)*MD*DEC001o*(11.D0/54.D0)+x**(-1)*MD*DEC00o*(1.D0/6.D0)+x**(-1)*MD*DEC002*(-11.D0/54.D0)&
       +x**(-1)*MD*DEC001*(-11.D0/54.D0)+x**(-1)*MD*DEC00*(-1.D0/6.D0)+x**(-1)*MD*C002o*(-2.D0/9.D0)&
       +x**(-1)*MD*C001o*(-2.D0/9.D0)+x**(-1)*MD*C00o*(-1.D0/3.D0)+x**(-1)*MD*C002*(2.D0/9.D0)&
       +x**(-1)*MD*C001*(2.D0/9.D0)+x**(-1)*MD*C00*(1.D0/3.D0)&
       +DEC0022*(11.D0/54.D0)+DEC002*(4.D0/27.D0)+DEC0012*(11.D0/27.D0)+DEC0011*(11.D0/54.D0)&
       +DEC001*(5.D0/27.D0)+C0022*(1.D0/9.D0)+C002*(-2.D0/9.D0)+C0012*(-4.D0/9.D0)+C0011*(-5.D0/9.D0)&
       +C001*(-4.D0/9.D0)+MD*DEC002*(-11.D0/54.D0)+MD*DEC001*(-11.D0/54.D0)+MD*DEC00*(-5.D0/27.D0)&
       +MD*C002*(2.D0/9.D0)+MD*C001*(5.D0/9.D0)+MD*C00*(4.D0/9.D0)+x*C1222*(-1.D0/9.D0)+x*C122*(-1.D0/9.D0)&
       +x*C1122*(-2.D0/9.D0)+x*C112*(-1.D0/9.D0)+x*C1112*(-1.D0/9.D0)+x*MD*C122*(1.D0/9.D0)&
       +x*MD*C12*(1.D0/9.D0)+x*MD*C112*(1.D0/9.D0)&
       )

  A(3) = A(3) &
       +facv*( &
       +x**(-1)*t*DEBf001o*(-1.D0/216.D0)+x**(-1)*t*DEBf001*(1.D0/216.D0)+x**(-1)*t*Bf001o*(1.D0/36.D0)&
       +x**(-1)*t*Bf001*(-1.D0/36.D0)+x**(-1)*t*DEC0022o*(1.D0/54.D0)+x**(-1)*t*DEC0012o*(1.D0/27.D0)&
       +x**(-1)*t*DEC0011o*(1.D0/54.D0)+x**(-1)*t*DEC0000o*(-5.D0/108.D0)+x**(-1)*t*DEC0022*(-1.D0/54.D0)&
       +x**(-1)*t*DEC0012*(-1.D0/27.D0)+x**(-1)*t*DEC0011*(-1.D0/54.D0)+x**(-1)*t*DEC0000*(5.D0/108.D0)&
       +x**(-1)*t*C0022o*(1.D0/18.D0)+x**(-1)*t*C0012o*(1.D0/9.D0)+x**(-1)*t*C0011o*(1.D0/18.D0)&
       +x**(-1)*t*C0000o*(11.D0/18.D0)+x**(-1)*t*C0022*(-1.D0/18.D0)+x**(-1)*t*C0012*(-1.D0/9.D0)&
       +x**(-1)*t*C0011*(-1.D0/18.D0)+x**(-1)*t*C0000*(-11.D0/18.D0)+x**(-1)*t*MD*DEBf00o*(1.D0/216.D0)&
       +x**(-1)*t*MD*DEBf00*(-1.D0/216.D0)+x**(-1)*t*MD*Bf00o*(-1.D0/36.D0)+x**(-1)*t*MD*Bf00*(1.D0/36.D0)&
       +x**(-1)*t*MD*DEC002o*(-1.D0/54.D0)+x**(-1)*t*MD*DEC001o*(-1.D0/54.D0)+x**(-1)*t*MD*DEC002*(1.D0/54.D0)&
       +x**(-1)*t*MD*DEC001*(1.D0/54.D0)+x**(-1)*t*MD*C002o*(-1.D0/18.D0)+x**(-1)*t*MD*C001o*(-1.D0/18.D0)&
       +x**(-1)*t*MD*C002*(1.D0/18.D0)+x**(-1)*t*MD*C001*(1.D0/18.D0)&
       +DEBf001o*(19.D0/54.D0)+DEBf001*(1.D0/54.D0)+Bf001o*(-7.D0/9.D0)+Bf001*(-1.D0/9.D0)+DEC0022o*(-1.D0/3.D0)&
       +DEC002o*(-1.D0/3.D0)+DEC0012o*(-2.D0/3.D0)+DEC0011o*(-1.D0/3.D0)+DEC001o*(-1.D0/3.D0)&
       +DEC0000o*(-1.D0/3.D0)+DEC0000*(11.D0/27.D0)+C0022o*(2.D0/3.D0)+C002o*(2.D0/3.D0)+C0012o*(4.D0/3.D0)&
       +C0011o*(2.D0/3.D0)+C001o*(2.D0/3.D0)+C0000o*(2.D0/3.D0)+C0000*(-10.D0/9.D0)+MD*DEBf00o*(-19.D0/54.D0)&
       +MD*DEBf00*(1.D0/54.D0)+MD*Bf00o*(7.D0/9.D0)+MD*Bf00*(-1.D0/9.D0)+MD*DEC002o*(1.D0/3.D0)&
       +MD*DEC001o*(1.D0/3.D0)+MD*DEC00o*(1.D0/3.D0)+MD*C002o*(-2.D0/3.D0)+MD*C001o*(-2.D0/3.D0)&
       +MD*C00o*(-2.D0/3.D0)+t*DEBf001*(1.D0/108.D0)+t*Bf001*(-1.D0/18.D0)+t*DEC0022*(-5.D0/108.D0)&
       +t*DEC0012*(-1.D0/18.D0)+t*DEC0011*(-1.D0/108.D0)+t*DEC0000*(2.D0/27.D0)+t*C0022*(-5.D0/36.D0)&
       +t*C0012*(-1.D0/6.D0)+t*C0011*(-1.D0/36.D0)+t*C0000*(-4.D0/9.D0)+t*MD*DEC002*(1.D0/36.D0)&
       +t*MD*DEC001*(1.D0/108.D0)+t*MD*C002*(1.D0/12.D0)+t*MD*C001*(1.D0/36.D0)+x*DEBf001*(1.D0/27.D0)&
       +x*Bf001*(-2.D0/9.D0)+x*DEC0022*(-1.D0/27.D0)+x*DEC0012*(-2.D0/27.D0)+x*DEC0011*(-1.D0/27.D0)&
       +x*DEC0000*(8.D0/27.D0)+x*C0022*(-1.D0/9.D0)+x*C0012*(-2.D0/9.D0)+x*C0011*(-1.D0/9.D0)&
       +x*C0000*(-16.D0/9.D0)+x*MD*DEC002*(1.D0/27.D0)+x*MD*DEC001*(1.D0/27.D0)+x*MD*C002*(1.D0/9.D0)&
       +x*MD*C001*(1.D0/9.D0)+x*t*DEC0022*(-1.D0/54.D0)+x*t*DEC0012*(-1.D0/108.D0)+x*t*C0022*(-1.D0/18.D0)&
       +x*t*C0012*(-1.D0/9.D0)+x**2*DEC0022*(-2.D0/27.D0)+x**2*DEC0012*(-1.D0/27.D0)+x**2*C0022*(-2.D0/9.D0)&
       +x**2*C0012*(-4.D0/9.D0)&
       )

  A(4) = A(4) &
       +facv*( &
       +x**(-1)*DEBf001o*(-1.D0/108.D0)+x**(-1)*DEBf001*(1.D0/108.D0)+x**(-1)*Bf001o*(1.D0/18.D0)&
       +x**(-1)*Bf001*(-1.D0/18.D0)+x**(-1)*DEC0022o*(1.D0/27.D0)+x**(-1)*DEC0012o*(2.D0/27.D0)&
       +x**(-1)*DEC0011o*(1.D0/27.D0)+x**(-1)*DEC0000o*(-5.D0/54.D0)+x**(-1)*DEC0022*(-1.D0/27.D0)&
       +x**(-1)*DEC0012*(-2.D0/27.D0)+x**(-1)*DEC0011*(-1.D0/27.D0)+x**(-1)*DEC0000*(5.D0/54.D0)&
       +x**(-1)*C0022o*(1.D0/9.D0)+x**(-1)*C0012o*(2.D0/9.D0)+x**(-1)*C0011o*(1.D0/9.D0)&
       +x**(-1)*C0000o*(11.D0/9.D0)+x**(-1)*C0022*(-1.D0/9.D0)+x**(-1)*C0012*(-2.D0/9.D0)&
       +x**(-1)*C0011*(-1.D0/9.D0)+x**(-1)*C0000*(-11.D0/9.D0)+x**(-1)*MD*DEBf00o*(1.D0/108.D0)&
       +x**(-1)*MD*DEBf00*(-1.D0/108.D0)+x**(-1)*MD*Bf00o*(-1.D0/18.D0)+x**(-1)*MD*Bf00*(1.D0/18.D0)&
       +x**(-1)*MD*DEC002o*(-1.D0/27.D0)+x**(-1)*MD*DEC001o*(-1.D0/27.D0)+x**(-1)*MD*DEC002*(1.D0/27.D0)&
       +x**(-1)*MD*DEC001*(1.D0/27.D0)+x**(-1)*MD*C002o*(-1.D0/9.D0)+x**(-1)*MD*C001o*(-1.D0/9.D0)&
       +x**(-1)*MD*C002*(1.D0/9.D0)+x**(-1)*MD*C001*(1.D0/9.D0)&
       +DEBf001*(1.D0/54.D0)+Bf001*(-1.D0/9.D0)+DEC0022*(-5.D0/54.D0)+DEC0012*(-1.D0/9.D0)&
       +DEC0011*(-1.D0/54.D0)+DEC0000*(4.D0/27.D0)+C0022*(-5.D0/18.D0)+C0012*(-1.D0/3.D0)&
       +C0011*(-1.D0/18.D0)+C0000*(-8.D0/9.D0)+MD*DEC002*(1.D0/18.D0)+MD*DEC001*(1.D0/54.D0)&
       +MD*C002*(1.D0/6.D0)+MD*C001*(1.D0/18.D0)+x*DEC0022*(-1.D0/27.D0)+x*DEC0012*(-1.D0/54.D0)&
       +x*C0022*(-1.D0/9.D0)+x*C0012*(-2.D0/9.D0)&
       )

  A(5) = A(5) &
       +facv*( &
       +x**(-1)*DEBf001o*(1.D0/108.D0)+x**(-1)*DEBf001*(-1.D0/108.D0)+x**(-1)*Bf001o*(-1.D0/18.D0)&
       +x**(-1)*Bf001*(1.D0/18.D0)+x**(-1)*DEC0022o*(-1.D0/27.D0)+x**(-1)*DEC0012o*(-2.D0/27.D0)&
       +x**(-1)*DEC0011o*(-1.D0/27.D0)+x**(-1)*DEC0000o*(5.D0/54.D0)+x**(-1)*DEC0022*(1.D0/27.D0)&
       +x**(-1)*DEC0012*(2.D0/27.D0)+x**(-1)*DEC0011*(1.D0/27.D0)+x**(-1)*DEC0000*(-5.D0/54.D0)&
       +x**(-1)*C0022o*(-1.D0/9.D0)+x**(-1)*C0012o*(-2.D0/9.D0)+x**(-1)*C0011o*(-1.D0/9.D0)&
       +x**(-1)*C0000o*(-11.D0/9.D0)+x**(-1)*C0022*(1.D0/9.D0)+x**(-1)*C0012*(2.D0/9.D0)&
       +x**(-1)*C0011*(1.D0/9.D0)+x**(-1)*C0000*(11.D0/9.D0)+x**(-1)*MD*DEBf00o*(-1.D0/108.D0)&
       +x**(-1)*MD*DEBf00*(1.D0/108.D0)+x**(-1)*MD*Bf00o*(1.D0/18.D0)+x**(-1)*MD*Bf00*(-1.D0/18.D0)&
       +x**(-1)*MD*DEC002o*(1.D0/27.D0)+x**(-1)*MD*DEC001o*(1.D0/27.D0)+x**(-1)*MD*DEC002*(-1.D0/27.D0)&
       +x**(-1)*MD*DEC001*(-1.D0/27.D0)+x**(-1)*MD*C002o*(1.D0/9.D0)+x**(-1)*MD*C001o*(1.D0/9.D0)&
       +x**(-1)*MD*C002*(-1.D0/9.D0)+x**(-1)*MD*C001*(-1.D0/9.D0)&
       +DEBf001*(1.D0/36.D0)+DEC0022*(-1.D0/108.D0)+DEC002*(-5.D0/54.D0)+DEC0012*(1.D0/108.D0)&
       +DEC0011*(1.D0/54.D0)+DEC001*(1.D0/108.D0)+DEC0000*(-1.D0/6.D0)+C0022*(17.D0/36.D0)&
       +C002*(2.D0/9.D0)+C0012*(4.D0/9.D0)+C0011*(-1.D0/36.D0)+C001*(-1.D0/18.D0)+C0000*(2.D0/3.D0)&
       +MD*DEC002*(-1.D0/18.D0)+MD*DEC001*(-1.D0/54.D0)+MD*DEC00*(-1.D0/108.D0)+MD*C002*(-1.D0/12.D0)&
       +MD*C001*(1.D0/36.D0)+MD*C00*(1.D0/18.D0)+x*DEC0022*(-7.D0/108.D0)+x*DEC002*(-5.D0/54.D0)&
       +x*DEC0012*(-11.D0/108.D0)+x*C0022*(2.D0/9.D0)+x*C002*(2.D0/9.D0)+x*C0012*(4.D0/9.D0)+x**2*C1222*(1.D0/18.D0)&
       +x**2*C122*(1.D0/18.D0)+x**2*C1122*(1.D0/18.D0)&
       )

  A(6) = A(6) &
       +facv*( &
       +x**(-1)*DEBf001o*(1.D0/27.D0)+x**(-1)*DEBf001*(-1.D0/27.D0)+x**(-1)*Bf001o*(-1.D0/18.D0)&
       +x**(-1)*Bf001*(1.D0/18.D0)+x**(-1)*DEC0022o*(-7.D0/108.D0)+x**(-1)*DEC002o*(-1.D0/12.D0)&
       +x**(-1)*DEC0012o*(-7.D0/54.D0)+x**(-1)*DEC0011o*(-7.D0/108.D0)+x**(-1)*DEC001o*(-1.D0/12.D0)&
       +x**(-1)*DEC0000o*(-7.D0/54.D0)+x**(-1)*DEC0022*(7.D0/108.D0)+x**(-1)*DEC002*(1.D0/12.D0)&
       +x**(-1)*DEC0012*(7.D0/54.D0)+x**(-1)*DEC0011*(7.D0/108.D0)+x**(-1)*DEC001*(1.D0/12.D0)&
       +x**(-1)*DEC0000*(7.D0/54.D0)+x**(-1)*C0022o*(2.D0/9.D0)+x**(-1)*C002o*(1.D0/6.D0)&
       +x**(-1)*C0012o*(4.D0/9.D0)+x**(-1)*C0011o*(2.D0/9.D0)+x**(-1)*C001o*(1.D0/6.D0)&
       +x**(-1)*C0000o*(7.D0/9.D0)+x**(-1)*C0022*(-2.D0/9.D0)+x**(-1)*C002*(-1.D0/6.D0)&
       +x**(-1)*C0012*(-4.D0/9.D0)+x**(-1)*C0011*(-2.D0/9.D0)+x**(-1)*C001*(-1.D0/6.D0)&
       +x**(-1)*C0000*(-7.D0/9.D0)+x**(-1)*MD*DEBf00o*(-1.D0/27.D0)+x**(-1)*MD*DEBf00*(1.D0/27.D0)&
       +x**(-1)*MD*Bf00o*(1.D0/18.D0)+x**(-1)*MD*Bf00*(-1.D0/18.D0)+x**(-1)*MD*DEC002o*(7.D0/108.D0)&
       +x**(-1)*MD*DEC001o*(7.D0/108.D0)+x**(-1)*MD*DEC00o*(1.D0/12.D0)+x**(-1)*MD*DEC002*(-7.D0/108.D0)&
       +x**(-1)*MD*DEC001*(-7.D0/108.D0)+x**(-1)*MD*DEC00*(-1.D0/12.D0)+x**(-1)*MD*C002o*(-2.D0/9.D0)&
       +x**(-1)*MD*C001o*(-2.D0/9.D0)+x**(-1)*MD*C00o*(-1.D0/6.D0)+x**(-1)*MD*C002*(2.D0/9.D0)&
       +x**(-1)*MD*C001*(2.D0/9.D0)+x**(-1)*MD*C00*(1.D0/6.D0)&
       +DEBf001*(-1.D0/36.D0)+DEC0022*(1.D0/9.D0)+DEC002*(1.D0/6.D0)+DEC0012*(7.D0/36.D0)&
       +DEC0011*(1.D0/12.D0)+DEC001*(1.D0/12.D0)+DEC0000*(1.D0/6.D0)+C0022*(-5.D0/12.D0)&
       +C002*(-1.D0/3.D0)+C0012*(-2.D0/3.D0)+C0011*(-1.D0/4.D0)+C001*(-1.D0/6.D0)+C0000*(-2.D0/3.D0)&
       +MD*DEC002*(-5.D0/108.D0)+MD*DEC001*(-1.D0/12.D0)+MD*DEC00*(-1.D0/12.D0)+MD*C002*(7.D0/36.D0)&
       +MD*C001*(1.D0/4.D0)+MD*C00*(1.D0/6.D0)+x*DEC0022*(7.D0/108.D0)+x*DEC002*(5.D0/54.D0)&
       +x*DEC0012*(11.D0/108.D0)+x*C0022*(-2.D0/9.D0)+x*C002*(-2.D0/9.D0)+x*C1222*(-1.D0/18.D0)&
       +x*C122*(-1.D0/18.D0)+x*C0012*(-4.D0/9.D0)+x*C1122*(-1.D0/9.D0)+x*C112*(-1.D0/18.D0)&
       +x*C1112*(-1.D0/18.D0)+x*MD*C122*(1.D0/18.D0)+x*MD*C12*(1.D0/18.D0)+x*MD*C112*(1.D0/18.D0)&
       +x**2*C1222*(-1.D0/18.D0)+x**2*C122*(-1.D0/18.D0)+x**2*C1122*(-1.D0/18.D0)&
       )

  A(7) = A(7)

  A(8) = A(8) &
       +facv*( &
       +x**(-1)*DEBf001o*(1.D0/108.D0)+x**(-1)*DEBf001*(-1.D0/108.D0)+x**(-1)*Bf001o*(-1.D0/18.D0)&
       +x**(-1)*Bf001*(1.D0/18.D0)+x**(-1)*DEC0022o*(-1.D0/27.D0)+x**(-1)*DEC0012o*(-2.D0/27.D0)&
       +x**(-1)*DEC0011o*(-1.D0/27.D0)+x**(-1)*DEC0000o*(5.D0/54.D0)+x**(-1)*DEC0022*(1.D0/27.D0)&
       +x**(-1)*DEC0012*(2.D0/27.D0)+x**(-1)*DEC0011*(1.D0/27.D0)+x**(-1)*DEC0000*(-5.D0/54.D0)&
       +x**(-1)*C0022o*(-1.D0/9.D0)+x**(-1)*C0012o*(-2.D0/9.D0)+x**(-1)*C0011o*(-1.D0/9.D0)&
       +x**(-1)*C0000o*(-11.D0/9.D0)+x**(-1)*C0022*(1.D0/9.D0)+x**(-1)*C0012*(2.D0/9.D0)&
       +x**(-1)*C0011*(1.D0/9.D0)+x**(-1)*C0000*(11.D0/9.D0)+x**(-1)*MD*DEBf00o*(-1.D0/108.D0)&
       +x**(-1)*MD*DEBf00*(1.D0/108.D0)+x**(-1)*MD*Bf00o*(1.D0/18.D0)+x**(-1)*MD*Bf00*(-1.D0/18.D0)&
       +x**(-1)*MD*DEC002o*(1.D0/27.D0)+x**(-1)*MD*DEC001o*(1.D0/27.D0)+x**(-1)*MD*DEC002*(-1.D0/27.D0)&
       +x**(-1)*MD*DEC001*(-1.D0/27.D0)+x**(-1)*MD*C002o*(1.D0/9.D0)+x**(-1)*MD*C001o*(1.D0/9.D0)&
       +x**(-1)*MD*C002*(-1.D0/9.D0)+x**(-1)*MD*C001*(-1.D0/9.D0)&
       +DEBf001*(-1.D0/54.D0)+Bf001*(1.D0/9.D0)+DEC0022*(5.D0/54.D0)+DEC0012*(1.D0/9.D0)&
       +DEC0011*(1.D0/54.D0)+DEC0000*(-4.D0/27.D0)+C0022*(5.D0/18.D0)+C0012*(1.D0/3.D0)&
       +C0011*(1.D0/18.D0)+C0000*(8.D0/9.D0)+MD*DEC002*(-1.D0/18.D0)+MD*DEC001*(-1.D0/54.D0)&
       +MD*C002*(-1.D0/6.D0)+MD*C001*(-1.D0/18.D0)+x*DEC0022*(1.D0/27.D0)+x*DEC0012*(1.D0/54.D0)&
       +x*C0022*(1.D0/9.D0)+x*C0012*(2.D0/9.D0)&
       )

  !SE

  A(1) = A(1) &
       +face*( &
       +x**(-2)*t*DEBf001o*(1.D0/16.D0)+x**(-2)*t*DEBf001*(-1.D0/16.D0)+x**(-2)*t*Bf001o*(-1.D0/8.D0)&
       +x**(-2)*t*Bf001*(1.D0/8.D0)+x**(-2)*t*MD*DEBf00o*(-1.D0/16.D0)+x**(-2)*t*MD*DEBf00*(1.D0/16.D0)&
       +x**(-2)*t*MD*Bf00o*(1.D0/8.D0)+x**(-2)*t*MD*Bf00*(-1.D0/8.D0)&
       +x**(-1)*t*DEC0022o*(-1.D0/8.D0)+x**(-1)*t*DEC002o*(-1.D0/8.D0)+x**(-1)*t*DEC0012o*(-1.D0/4.D0)&
       +x**(-1)*t*DEC0011o*(-1.D0/8.D0)+x**(-1)*t*DEC001o*(-1.D0/8.D0)+x**(-1)*t*DEC0000o*(-1.D0/8.D0)&
       +x**(-1)*t*C0022o*(1.D0/4.D0)+x**(-1)*t*C002o*(1.D0/4.D0)+x**(-1)*t*C0012o*(1.D0/2.D0)&
       +x**(-1)*t*C0011o*(1.D0/4.D0)+x**(-1)*t*C001o*(1.D0/4.D0)+x**(-1)*t*C0000o*(1.D0/4.D0)&
       +x**(-1)*t*DEBf001o*(1.D0/8.D0)+x**(-1)*t*DEBf001*(-3.D0/16.D0)+x**(-1)*t*Bf001o*(-1.D0/4.D0)&
       +x**(-1)*t*Bf001*(3.D0/8.D0)+x**(-1)*t*MD*DEC002o*(1.D0/8.D0)+x**(-1)*t*MD*DEC001o*(1.D0/8.D0)&
       +x**(-1)*t*MD*DEC00o*(1.D0/8.D0)+x**(-1)*t*MD*C002o*(-1.D0/4.D0)+x**(-1)*t*MD*C001o*(-1.D0/4.D0)&
       +x**(-1)*t*MD*C00o*(-1.D0/4.D0)+x**(-1)*t*MD*DEBf00o*(-1.D0/8.D0)+x**(-1)*t*MD*DEBf00*(1.D0/8.D0)&
       +x**(-1)*t*MD*Bf00o*(1.D0/4.D0)+x**(-1)*t*MD*Bf00*(-1.D0/4.D0)&
       +DEC0022o*(-1.D0/2.D0)+DEC002o*(-1.D0/2.D0)+DEC0012o*(-1.D0)+DEC0011o*(-1.D0/2.D0)+DEC001o*(-1.D0/2.D0)&
       +DEC0000o*(-1.D0/2.D0)+C0022o*(1.D0)+C002o*(1.D0)+C0012o*(2.D0)+C0011o*(1.D0)+C001o*(1.D0)+C0000o*(1.D0)&
       +DEBf001o*(1.D0/2.D0)+DEBf001*(-1.D0/4.D0)+Bf001o*(-1.D0)+Bf001*(1.D0/2.D0)+MD*DEC002o*(1.D0/2.D0)&
       +MD*DEC001o*(1.D0/2.D0)+MD*DEC00o*(1.D0/2.D0)+MD*C002o*(-1.D0)+MD*C001o*(-1.D0)+MD*C00o*(-1.D0)&
       +MD*DEBf00o*(-1.D0/2.D0)+MD*Bf00o*(1.D0)+t*DEBf001*(-1.D0/8.D0)+t*Bf001*(1.D0/4.D0)+x*DEBf001*(-1.D0/2.D0)&
       +x*Bf001*(1.D0)&
       )

  A(2) = A(2)

  A(3) = A(3) &
       +face*( &
       +x**(-2)*t*DEBf001o*(1.D0/16.D0)+x**(-2)*t*DEBf001*(-1.D0/16.D0)+x**(-2)*t*Bf001o*(-1.D0/8.D0)&
       +x**(-2)*t*Bf001*(1.D0/8.D0)+x**(-2)*t*MD*DEBf00o*(-1.D0/16.D0)+x**(-2)*t*MD*DEBf00*(1.D0/16.D0)&
       +x**(-2)*t*MD*Bf00o*(1.D0/8.D0)+x**(-2)*t*MD*Bf00*(-1.D0/8.D0)&
       +x**(-1)*t*DEC0022o*(-1.D0/8.D0)+x**(-1)*t*DEC002o*(-1.D0/8.D0)+x**(-1)*t*DEC0012o*(-1.D0/4.D0)&
       +x**(-1)*t*DEC0011o*(-1.D0/8.D0)+x**(-1)*t*DEC001o*(-1.D0/8.D0)+x**(-1)*t*DEC0000o*(-1.D0/8.D0)&
       +x**(-1)*t*C0022o*(1.D0/4.D0)+x**(-1)*t*C002o*(1.D0/4.D0)+x**(-1)*t*C0012o*(1.D0/2.D0)&
       +x**(-1)*t*C0011o*(1.D0/4.D0)+x**(-1)*t*C001o*(1.D0/4.D0)+x**(-1)*t*C0000o*(1.D0/4.D0)&
       +x**(-1)*t*DEBf001o*(1.D0/8.D0)+x**(-1)*t*DEBf001*(-3.D0/16.D0)+x**(-1)*t*Bf001o*(-1.D0/4.D0)&
       +x**(-1)*t*Bf001*(3.D0/8.D0)+x**(-1)*t*MD*DEC002o*(1.D0/8.D0)+x**(-1)*t*MD*DEC001o*(1.D0/8.D0)&
       +x**(-1)*t*MD*DEC00o*(1.D0/8.D0)+x**(-1)*t*MD*C002o*(-1.D0/4.D0)+x**(-1)*t*MD*C001o*(-1.D0/4.D0)&
       +x**(-1)*t*MD*C00o*(-1.D0/4.D0)+x**(-1)*t*MD*DEBf00o*(-1.D0/8.D0)+x**(-1)*t*MD*DEBf00*(1.D0/8.D0)&
       +x**(-1)*t*MD*Bf00o*(1.D0/4.D0)+x**(-1)*t*MD*Bf00*(-1.D0/4.D0)&
       +DEC0022o*(-1.D0/2.D0)+DEC002o*(-1.D0/2.D0)+DEC0012o*(-1.D0)+DEC0011o*(-1.D0/2.D0)&
       +DEC001o*(-1.D0/2.D0)+DEC0000o*(-1.D0/2.D0)+C0022o*(1.D0)+C002o*(1.D0)+C0012o*(2.D0)+C0011o*(1.D0)&
       +C001o*(1.D0)+C0000o*(1.D0)+DEBf001o*(1.D0/2.D0)+DEBf001*(-1.D0/4.D0)+Bf001o*(-1.D0)+Bf001*(1.D0/2.D0)&
       +MD*DEC002o*(1.D0/2.D0)+MD*DEC001o*(1.D0/2.D0)+MD*DEC00o*(1.D0/2.D0)+MD*C002o*(-1.D0)+MD*C001o*(-1.D0)&
       +MD*C00o*(-1.D0)+MD*DEBf00o*(-1.D0/2.D0)+MD*Bf00o*(1.D0)+t*DEBf001*(-1.D0/8.D0)+t*Bf001*(1.D0/4.D0)&
       +x*DEBf001*(-1.D0/2.D0)+x*Bf001*(1.D0)&
       )

  A(4) = A(4) &
       +face*( &
       +x**(-2)*DEBf001o*(1.D0/8.D0)+x**(-2)*DEBf001*(-1.D0/8.D0)+x**(-2)*Bf001o*(-1.D0/4.D0)&
       +x**(-2)*Bf001*(1.D0/4.D0)+x**(-2)*MD*DEBf00o*(-1.D0/8.D0)+x**(-2)*MD*DEBf00*(1.D0/8.D0)&
       +x**(-2)*MD*Bf00o*(1.D0/4.D0)+x**(-2)*MD*Bf00*(-1.D0/4.D0)&
       +x**(-1)*DEC0022o*(-1.D0/4.D0)+x**(-1)*DEC002o*(-1.D0/4.D0)+x**(-1)*DEC0012o*(-1.D0/2.D0)&
       +x**(-1)*DEC0011o*(-1.D0/4.D0)+x**(-1)*DEC001o*(-1.D0/4.D0)+x**(-1)*DEC0000o*(-1.D0/4.D0)&
       +x**(-1)*C0022o*(1.D0/2.D0)+x**(-1)*C002o*(1.D0/2.D0)+x**(-1)*C0012o*(1.D0)&
       +x**(-1)*C0011o*(1.D0/2.D0)+x**(-1)*C001o*(1.D0/2.D0)+x**(-1)*C0000o*(1.D0/2.D0)&
       +x**(-1)*DEBf001o*(1.D0/4.D0)+x**(-1)*DEBf001*(-3.D0/8.D0)+x**(-1)*Bf001o*(-1.D0/2.D0)&
       +x**(-1)*Bf001*(3.D0/4.D0)+x**(-1)*MD*DEC002o*(1.D0/4.D0)+x**(-1)*MD*DEC001o*(1.D0/4.D0)&
       +x**(-1)*MD*DEC00o*(1.D0/4.D0)+x**(-1)*MD*C002o*(-1.D0/2.D0)+x**(-1)*MD*C001o*(-1.D0/2.D0)&
       +x**(-1)*MD*C00o*(-1.D0/2.D0)+x**(-1)*MD*DEBf00o*(-1.D0/4.D0)+x**(-1)*MD*DEBf00*(1.D0/4.D0)&
       +x**(-1)*MD*Bf00o*(1.D0/2.D0)+x**(-1)*MD*Bf00*(-1.D0/2.D0)&
       +DEBf001*(-1.D0/4.D0)+Bf001*(1.D0/2.D0)&
       )

  A(5) = A(5) &
       +face*( &
       +x**(-2)*DEBf001o*(-1.D0/8.D0)+x**(-2)*DEBf001*(1.D0/8.D0)+x**(-2)*Bf001o*(1.D0/4.D0)&
       +x**(-2)*Bf001*(-1.D0/4.D0)+x**(-2)*MD*DEBf00o*(1.D0/8.D0)+x**(-2)*MD*DEBf00*(-1.D0/8.D0)&
       +x**(-2)*MD*Bf00o*(-1.D0/4.D0)+x**(-2)*MD*Bf00*(1.D0/4.D0)&
       +x**(-1)*DEC0022o*(1.D0/4.D0)+x**(-1)*DEC002o*(1.D0/4.D0)+x**(-1)*DEC0012o*(1.D0/2.D0)&
       +x**(-1)*DEC0011o*(1.D0/4.D0)+x**(-1)*DEC001o*(1.D0/4.D0)+x**(-1)*DEC0000o*(1.D0/4.D0)&
       +x**(-1)*C0022o*(-1.D0/2.D0)+x**(-1)*C002o*(-1.D0/2.D0)+x**(-1)*C0012o*(-1.D0)&
       +x**(-1)*C0011o*(-1.D0/2.D0)+x**(-1)*C001o*(-1.D0/2.D0)+x**(-1)*C0000o*(-1.D0/2.D0)&
       +x**(-1)*DEBf001o*(-1.D0/4.D0)+x**(-1)*DEBf001*(3.D0/8.D0)+x**(-1)*Bf001o*(1.D0/2.D0)&
       +x**(-1)*Bf001*(-3.D0/4.D0)+x**(-1)*MD*DEC002o*(-1.D0/4.D0)+x**(-1)*MD*DEC001o*(-1.D0/4.D0)&
       +x**(-1)*MD*DEC00o*(-1.D0/4.D0)+x**(-1)*MD*C002o*(1.D0/2.D0)+x**(-1)*MD*C001o*(1.D0/2.D0)&
       +x**(-1)*MD*C00o*(1.D0/2.D0)+x**(-1)*MD*DEBf00o*(1.D0/4.D0)+x**(-1)*MD*DEBf00*(-1.D0/4.D0)&
       +x**(-1)*MD*Bf00o*(-1.D0/2.D0)+x**(-1)*MD*Bf00*(1.D0/2.D0)&
       +DEBf001*(1.D0/4.D0)+Bf001*(-1.D0/2.D0)&
       )

  A(6) = A(6) &
       +face*( &
       +x**(-2)*DEBf001o*(1.D0/8.D0)+x**(-2)*DEBf001*(-1.D0/8.D0)+x**(-2)*Bf001o*(-1.D0/4.D0)&
       +x**(-2)*Bf001*(1.D0/4.D0)+x**(-2)*MD*DEBf00o*(-1.D0/8.D0)+x**(-2)*MD*DEBf00*(1.D0/8.D0)&
       +x**(-2)*MD*Bf00o*(1.D0/4.D0)+x**(-2)*MD*Bf00*(-1.D0/4.D0)&
       +x**(-1)*DEC0022o*(-1.D0/4.D0)+x**(-1)*DEC002o*(-1.D0/4.D0)+x**(-1)*DEC0012o*(-1.D0/2.D0)&
       +x**(-1)*DEC0011o*(-1.D0/4.D0)+x**(-1)*DEC001o*(-1.D0/4.D0)+x**(-1)*DEC0000o*(-1.D0/4.D0)&
       +x**(-1)*C0022o*(1.D0/2.D0)+x**(-1)*C002o*(1.D0/2.D0)+x**(-1)*C0012o*(1.D0)&
       +x**(-1)*C0011o*(1.D0/2.D0)+x**(-1)*C001o*(1.D0/2.D0)+x**(-1)*C0000o*(1.D0/2.D0)&
       +x**(-1)*DEBf001o*(1.D0/4.D0)+x**(-1)*DEBf001*(-3.D0/8.D0)+x**(-1)*Bf001o*(-1.D0/2.D0)&
       +x**(-1)*Bf001*(3.D0/4.D0)+x**(-1)*MD*DEC002o*(1.D0/4.D0)+x**(-1)*MD*DEC001o*(1.D0/4.D0)&
       +x**(-1)*MD*DEC00o*(1.D0/4.D0)+x**(-1)*MD*C002o*(-1.D0/2.D0)+x**(-1)*MD*C001o*(-1.D0/2.D0)&
       +x**(-1)*MD*C00o*(-1.D0/2.D0)+x**(-1)*MD*DEBf00o*(-1.D0/4.D0)+x**(-1)*MD*DEBf00*(1.D0/4.D0)&
       +x**(-1)*MD*Bf00o*(1.D0/2.D0)+x**(-1)*MD*Bf00*(-1.D0/2.D0)&
       +DEBf001*(-1.D0/4.D0)+Bf001*(1.D0/2.D0)&
       )

  A(7) = A(7)

  A(8) = A(8) &
       +face*( &
       +x**(-2)*DEBf001o*(-1.D0/8.D0)+x**(-2)*DEBf001*(1.D0/8.D0)+x**(-2)*Bf001o*(1.D0/4.D0)&
       +x**(-2)*Bf001*(-1.D0/4.D0)+x**(-2)*MD*DEBf00o*(1.D0/8.D0)+x**(-2)*MD*DEBf00*(-1.D0/8.D0)&
       +x**(-2)*MD*Bf00o*(-1.D0/4.D0)+x**(-2)*MD*Bf00*(1.D0/4.D0)&
       +x**(-1)*DEC0022o*(1.D0/4.D0)+x**(-1)*DEC002o*(1.D0/4.D0)+x**(-1)*DEC0012o*(1.D0/2.D0)&
       +x**(-1)*DEC0011o*(1.D0/4.D0)+x**(-1)*DEC001o*(1.D0/4.D0)+x**(-1)*DEC0000o*(1.D0/4.D0)&
       +x**(-1)*C0022o*(-1.D0/2.D0)+x**(-1)*C002o*(-1.D0/2.D0)+x**(-1)*C0012o*(-1.D0)&
       +x**(-1)*C0011o*(-1.D0/2.D0)+x**(-1)*C001o*(-1.D0/2.D0)+x**(-1)*C0000o*(-1.D0/2.D0)&
       +x**(-1)*DEBf001o*(-1.D0/4.D0)+x**(-1)*DEBf001*(3.D0/8.D0)+x**(-1)*Bf001o*(1.D0/2.D0)&
       +x**(-1)*Bf001*(-3.D0/4.D0)+x**(-1)*MD*DEC002o*(-1.D0/4.D0)+x**(-1)*MD*DEC001o*(-1.D0/4.D0)&
       +x**(-1)*MD*DEC00o*(-1.D0/4.D0)+x**(-1)*MD*C002o*(1.D0/2.D0)+x**(-1)*MD*C001o*(1.D0/2.D0)&
       +x**(-1)*MD*C00o*(1.D0/2.D0)+x**(-1)*MD*DEBf00o*(1.D0/4.D0)+x**(-1)*MD*DEBf00*(-1.D0/4.D0)&
       +x**(-1)*MD*Bf00o*(-1.D0/2.D0)+x**(-1)*MD*Bf00*(1.D0/2.D0)&
       +DEBf001*(1.D0/4.D0)+Bf001*(-1.D0/2.D0)&
       )

  !* Graph 6:
  !*  Two-photon vertex graph
  !* Isospin factor = 4/3

  p1s=nuclm**2
  p2s=t
  p1p2s=nuclm**2

  m1s=deltam**2
  m2s=mpi2
  m3s=mpi2

  CALL Setdelta(Delta)

  C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C11=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C22=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C001=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C002=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C111=C0i(cc111,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C112=C0i(cc112,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C122=C0i(cc122,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C222=C0i(cc222,p1s,p2s,p1p2s,m1s,m2s,m3s)

  if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
  else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
  endif

  C00D=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C001D=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
  C002D=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)

  DEC00=2.D0*(C00-C00D)/sw1
  DEC001=2.D0*(C001-C001D)/sw1
  DEC002=2.D0*(C002-C002D)/sw1

  !TP

  A(1) = A(1) &
       +(DEC002*(-2.D0/3.D0)+DEC001*(-2.D0/3.D0)+C002*(4.D0/3.D0)+C001*(4.D0/3.D0)+MD*DEC00*(2.D0/3.D0)&
       +MD*C00*(-4.D0/3.D0)+t*DEC002*(10.D0/27.D0)+t*DEC001*(10.D0/27.D0)+t*C002*(-8.D0/9.D0)&
       +t*C122*(-4.D0/9.D0)+t*C112*(-4.D0/9.D0)+t*C001*(-8.D0/9.D0)+t*MD*DEC00*(-10.D0/27.D0)&
       +t*MD*C12*(4.D0/9.D0)+t*MD*C00*(8.D0/9.D0)+t**2*C122*(1.D0/9.D0)+t**2*C112*(1.D0/9.D0)&
       +t**2*MD*C12*(-1.D0/9.D0))/2.D0

  A=gafac*A

  CALL Setdelta(DeltaIn)                ! [v2] no LoopTools state leak (audit B5)

END SUBROUTINE CorrIsoCovOneDeltaLoop

!=============================================================================
! [v2, 2026-10] On-shell subtraction constants of the reducible gammaNN-vertex
! graphs of CorrIsoCovOneDeltaLoop (graph 4): B- and C-functions at p^2 = M^2
! with a real photon. LoopTools' C-reduction is singular for p1p2s = p1s, q^2=0
! (zero Gram determinant), so v1 used p1p2s = M^2 + 1e-4 (one-sided).  Here the
! C-functions are the average over p1p2s = M^2 + eps and M^2 - eps (exact up to
! O(eps^2); the functions are analytic there, far below all thresholds).
! The results depend only on (mpi2, deltam, Delta, mu) and are cached.
! Mode/offset can be changed via COMMON /v2piDos/ (diagnostics only):
!   iosmode = 0 : symmetric average (default);  1 : one-sided v1 prescription.
!=============================================================================
SUBROUTINE piD_onshell_consts(mpi2,deltam,Delta, &
     Bf00o,Bf001o,DEBf00o,DEBf001o, &
     C00o,C001o,C002o,C0000o,C0011o,C0012o,C0022o, &
     DEC00o,DEC001o,DEC002o,DEC0000o,DEC0011o,DEC0012o,DEC0022o)
  implicit none
#include "looptools.h"
  REAL*8, INTENT(in) :: mpi2,deltam,Delta
  COMPLEX*16, INTENT(out) :: Bf00o,Bf001o,DEBf00o,DEBf001o
  COMPLEX*16, INTENT(out) :: C00o,C001o,C002o,C0000o,C0011o,C0012o,C0022o
  COMPLEX*16, INTENT(out) :: DEC00o,DEC001o,DEC002o,DEC0000o,DEC0011o,DEC0012o,DEC0022o
  REAL*8 epsos
  INTEGER iosmode
  COMMON /v2piDos/ epsos, iosmode
  COMPLEX*16, SAVE :: v2cache(18)
  REAL*8, SAVE :: v2key(6) = -1.d300
  REAL*8 DeltaAux, sw1, p1s, p2s, p1p2s, m1s, m2s, m3s, DeltaIn, eps, wgt
  COMPLEX*16 Bv(2), BvD(2), Cv(7), CvD(7)
  COMPLEX*16 Bacc(2), DBacc(2), Cacc(7), DCacc(7)
  INTEGER k, npt
  REAL*8, PARAMETER :: nuclm=1.d0

  eps=epsos
  IF(eps.LE.0.d0) eps=1.d-4
  IF(v2key(1).EQ.mpi2 .AND. v2key(2).EQ.deltam .AND. v2key(3).EQ.Delta .AND. &
       v2key(4).EQ.getmudim() .AND. v2key(5).EQ.eps .AND. v2key(6).EQ.DBLE(iosmode)) THEN
     Bf00o=v2cache(1); Bf001o=v2cache(2); DEBf00o=v2cache(3); DEBf001o=v2cache(4)
     C00o=v2cache(5); C001o=v2cache(6); C002o=v2cache(7); C0000o=v2cache(8)
     C0011o=v2cache(9); C0012o=v2cache(10); C0022o=v2cache(11)
     DEC00o=v2cache(12); DEC001o=v2cache(13); DEC002o=v2cache(14); DEC0000o=v2cache(15)
     DEC0011o=v2cache(16); DEC0012o=v2cache(17); DEC0022o=v2cache(18)
     RETURN
  ENDIF
  DeltaIn=getdelta()
  DeltaAux=1.d0
  IF(.not.Delta.eq.0.D0) THEN
     sw1=Delta
  ELSE
     sw1=-DeltaAux
  ENDIF

  !* B functions at p^2 = M^2 (no offset needed)
  p1s=nuclm**2
  m1s=deltam**2
  m2s=mpi2
  CALL Setdelta(Delta)
  Bv(1)=B0i(bb00,p1s,m1s,m2s)
  Bv(2)=B0i(bb001,p1s,m1s,m2s)
  IF(.not.Delta.eq.0.D0) THEN
     CALL Setdelta(0.D0)
  ELSE
     CALL Setdelta(DeltaAux)
  ENDIF
  BvD(1)=B0i(bb00,p1s,m1s,m2s)
  BvD(2)=B0i(bb001,p1s,m1s,m2s)
  Bf00o=Bv(1); Bf001o=Bv(2)
  DEBf00o=2.D0*(Bv(1)-BvD(1))/sw1
  DEBf001o=2.D0*(Bv(2)-BvD(2))/sw1

  !* C functions at p1p2s -> M^2
  p1s=nuclm**2
  p2s=0.d0
  m1s=deltam**2
  m2s=mpi2
  m3s=mpi2
  Cacc=0.d0; DCacc=0.d0
  IF(iosmode.EQ.1) THEN
     npt=1
  ELSE
     npt=2
  ENDIF
  wgt=1.d0/DBLE(npt)
  DO k=1,npt
     IF(k.EQ.1) THEN
        p1p2s=nuclm**2+eps
     ELSE
        p1p2s=nuclm**2-eps
     ENDIF
     CALL Setdelta(Delta)
     Cv(1)=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
     Cv(2)=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
     Cv(3)=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
     Cv(4)=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
     Cv(5)=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
     Cv(6)=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
     Cv(7)=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)
     IF(.not.Delta.eq.0.D0) THEN
        CALL Setdelta(0.D0)
     ELSE
        CALL Setdelta(DeltaAux)
     ENDIF
     CvD(1)=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
     CvD(2)=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
     CvD(3)=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
     CvD(4)=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
     CvD(5)=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
     CvD(6)=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
     CvD(7)=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)
     Cacc=Cacc+wgt*Cv
     DCacc=DCacc+wgt*2.D0*(Cv-CvD)/sw1
  ENDDO
  C00o=Cacc(1); C001o=Cacc(2); C002o=Cacc(3); C0000o=Cacc(4)
  C0011o=Cacc(5); C0012o=Cacc(6); C0022o=Cacc(7)
  DEC00o=DCacc(1); DEC001o=DCacc(2); DEC002o=DCacc(3); DEC0000o=DCacc(4)
  DEC0011o=DCacc(5); DEC0012o=DCacc(6); DEC0022o=DCacc(7)
  CALL Setdelta(DeltaIn)
  v2cache=(/Bf00o,Bf001o,DEBf00o,DEBf001o,C00o,C001o,C002o,C0000o,C0011o,C0012o,C0022o, &
       DEC00o,DEC001o,DEC002o,DEC0000o,DEC0011o,DEC0012o,DEC0022o/)
  v2key=(/mpi2,deltam,Delta,getmudim(),eps,DBLE(iosmode)/)
END SUBROUTINE piD_onshell_consts
!
! Older version of the pi-Delta loop contribution
! Isospin factors of 1PR graphs are fudged to provide current conservation
!
SUBROUTINE OldCorrIsoCovOneDeltaLoop(s,t,mpi2,gafac,A,Delta)
  USE trigd
  implicit COMPLEX*16 (B-D)

#include "looptools.h"

   REAL*8, INTENT(in) :: s,t,mpi2,gafac,Delta
   COMPLEX*16,INTENT(out) :: A(8)
   COMPLEX*16, PARAMETER :: CI=(0.d0,1.d0)
   REAL*8 x,nuclm,deltam, MD
   REAL*8 DeltaAux, sw1
   REAL*8 m1s,m2s,m3s,m4s
   REAL*8 p1s,p2s,p2so,p3s,p4s,p1p2s,p1p2so,p2p3s,ps,pso
   INTEGER i,j
   COMPLEX*16 Sum1, Sum2, Sum3, AA(8)

   A=0.d0
   AA=0.d0

   Sum1=0.d0
   Sum2=0.d0
   Sum3=0.d0

   nuclm=1.d0   !* it's assumed that all is given in nucleon-mass units
   deltam=1.31d0
   MD=deltam
   x=(s-1.d0)/2.d0

   DeltaAux=1.d0

   !* Graph 1:
   !*  Box graph with photons attached to the pions
   !* Isospin factor = 4/3

   p1s=nuclm**2
   p2s=0.0d0
   p3s=0.0d0
   p4s=nuclm**2
   p2p3s=t
   p1p2s=s

   m1s=deltam**2
   m4s=mpi2
   m2s=mpi2
   m3s=mpi2

   CALL Setdelta(Delta)

   D00=D0i(dd00,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11=D0i(dd11,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12=D0i(dd12,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D13=D0i(dd13,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22=D0i(dd22,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D23=D0i(dd23,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D33=D0i(dd33,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D111=D0i(dd111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D112=D0i(dd112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D113=D0i(dd113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D122=D0i(dd122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D123=D0i(dd123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D133=D0i(dd133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D222=D0i(dd222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D223=D0i(dd223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D233=D0i(dd233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D333=D0i(dd333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D001=D0i(dd001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D002=D0i(dd002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D003=D0i(dd003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0000=D0i(dd0000,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0011=D0i(dd0011,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0012=D0i(dd0012,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0013=D0i(dd0013,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0022=D0i(dd0022,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0023=D0i(dd0023,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D0033=D0i(dd0033,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1111=D0i(dd1111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1112=D0i(dd1112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1113=D0i(dd1113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1122=D0i(dd1122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1123=D0i(dd1123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1133=D0i(dd1133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1222=D0i(dd1222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1223=D0i(dd1223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1233=D0i(dd1233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D1333=D0i(dd1333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D2222=D0i(dd2222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D2223=D0i(dd2223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D2233=D0i(dd2233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D2333=D0i(dd2333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D3333=D0i(dd3333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00001=D0i(dd00001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00002=D0i(dd00002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00003=D0i(dd00003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00111=D0i(dd00111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00112=D0i(dd00112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00113=D0i(dd00113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00122=D0i(dd00122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00123=D0i(dd00123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00133=D0i(dd00133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00222=D0i(dd00222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00223=D0i(dd00223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00233=D0i(dd00233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00333=D0i(dd00333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11111=D0i(dd11111,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11112=D0i(dd11112,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11113=D0i(dd11113,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11122=D0i(dd11122,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11123=D0i(dd11123,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11133=D0i(dd11133,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11222=D0i(dd11222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11223=D0i(dd11223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11233=D0i(dd11233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D11333=D0i(dd11333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12222=D0i(dd12222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12223=D0i(dd12223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12233=D0i(dd12233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D12333=D0i(dd12333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D13333=D0i(dd13333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22222=D0i(dd22222,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22223=D0i(dd22223,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22233=D0i(dd22233,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D22333=D0i(dd22333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D23333=D0i(dd23333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D33333=D0i(dd33333,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)


   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   D0000D=D0i(dd0000,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00001D=D0i(dd00001,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00002D=D0i(dd00002,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)
   D00003D=D0i(dd00003,p1s,p2s,p3s,p4s,s,t,m1s,m2s,m3s,m4s)

   DED0000=2.D0*(D0000-D0000D)/sw1
   DED00001=2.D0*(D00001-D00001D)/sw1
   DED00002=2.D0*(D00002-D00002D)/sw1
   DED00003=2.D0*(D00003-D00003D)/sw1

   !* Graph 2:
   !*  Box graphs with one photon attached to the pion and the other attached to the piNDelta vertex
   !* Isospin factor = 4/3

   p1s=nuclm**2
   p2s=0.d0
   p1p2s=s

   m1s=deltam**2
   m2s=mpi2
   m3s=mpi2

   CALL Setdelta(Delta)

   C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C11=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C22=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C111=C0i(cc111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C112=C0i(cc112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C122=C0i(cc122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C222=C0i(cc222,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0000=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0011=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0012=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0022=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1111=C0i(cc1111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1112=C0i(cc1112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1122=C0i(cc1122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1222=C0i(cc1222,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C2222=C0i(cc2222,p1s,p2s,p1p2s,m1s,m2s,m3s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   C00D=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001D=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002D=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0000D=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0011D=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0012D=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0022D=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)

   DEC00=2.D0*(C00-C00D)/sw1
   DEC001=2.D0*(C001-C001D)/sw1
   DEC002=2.D0*(C002-C002D)/sw1
   DEC0000=2.D0*(C0000-C0000D)/sw1
   DEC0011=2.D0*(C0011-C0011D)/sw1
   DEC0012=2.D0*(C0012-C0012D)/sw1
   DEC0022=2.D0*(C0022-C0022D)/sw1

   !write(*,*) '2', Sum2

   !* Graph 3:
   !*  Box graph with photons attached to the piNDelta vertices
   !* Isospin factor = 4/3

   p1s=s

   m1s=deltam**2
   m2s=mpi2

   CALL Setdelta(Delta)

   Bf00=B0i(bb00,p1s,m1s,m2s)
   Bf11=B0i(bb11,p1s,m1s,m2s)
   Bf001=B0i(bb001,p1s,m1s,m2s)
   Bf111=B0i(bb111,p1s,m1s,m2s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   Bf00D=B0i(bb00,p1s,m1s,m2s)
   Bf001D=B0i(bb001,p1s,m1s,m2s)
   Bf11D=B0i(bb11,p1s,m1s,m2s)
   Bf111D=B0i(bb111,p1s,m1s,m2s)

   DEBf00=2.D0*(Bf00-Bf00D)/sw1
   DEBf001=2.D0*(Bf001-Bf001D)/sw1
   DEBf11=2.D0*(Bf11-Bf11D)/sw1
   DEBf111=2.D0*(Bf111-Bf111D)/sw1


    AA=A!*3.D0/4.D0
    A=0.d0

   !* Graph 4:
   !*  Reducible graphs with gammaNN vertices
   !* Isospin factor = -2/3

   p1s=nuclm**2

   m1s=deltam**2
   m2s=mpi2

   CALL Setdelta(Delta)

   Bf00o=B0i(bb00,p1s,m1s,m2s)
   Bf11o=B0i(bb11,p1s,m1s,m2s)
   Bf001o=B0i(bb001,p1s,m1s,m2s)
   Bf111o=B0i(bb111,p1s,m1s,m2s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   Bf00Do=B0i(bb00,p1s,m1s,m2s)
   Bf001Do=B0i(bb001,p1s,m1s,m2s)
   Bf11Do=B0i(bb11,p1s,m1s,m2s)
   Bf111Do=B0i(bb111,p1s,m1s,m2s)

   DEBf00o=2.D0*(Bf00o-Bf00Do)/sw1
   DEBf001o=2.D0*(Bf001o-Bf001Do)/sw1
   DEBf11o=2.D0*(Bf11o-Bf11Do)/sw1
   DEBf111o=2.D0*(Bf111o-Bf111Do)/sw1

   p1s=nuclm**2
   p2s=0.d0
   p1p2s=nuclm**2+0.0001d0

   m1s=deltam**2
   m2s=mpi2
   m3s=mpi2

   CALL Setdelta(Delta)

   C00o=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C11o=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12o=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C22o=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001o=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002o=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C111o=C0i(cc111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C112o=C0i(cc112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C122o=C0i(cc122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C222o=C0i(cc222,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0000o=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0011o=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0012o=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0022o=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1111o=C0i(cc1111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1112o=C0i(cc1112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1122o=C0i(cc1122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C1222o=C0i(cc1222,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C2222o=C0i(cc2222,p1s,p2s,p1p2s,m1s,m2s,m3s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   C00Do=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001Do=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002Do=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0000Do=C0i(cc0000,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0011Do=C0i(cc0011,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0012Do=C0i(cc0012,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C0022Do=C0i(cc0022,p1s,p2s,p1p2s,m1s,m2s,m3s)

   DEC00o=2.D0*(C00o-C00Do)/sw1
   DEC001o=2.D0*(C001o-C001Do)/sw1
   DEC002o=2.D0*(C002o-C002Do)/sw1
   DEC0000o=2.D0*(C0000o-C0000Do)/sw1
   DEC0011o=2.D0*(C0011o-C0011Do)/sw1
   DEC0012o=2.D0*(C0012o-C0012Do)/sw1
   DEC0022o=2.D0*(C0022o-C0022Do)/sw1


   !* Graph 5:
   !*  Self-energy graph
   !* Isospin factor = 2

   !Summing up

   A(1) = A(1) &
	+(+DED00002*(44.D0/27.D0)+DED00001*(88.D0/27.D0)+D00002*(-40.D0/9.D0)+D00001*(-80.D0/9.D0)+MD*DED0000*(-44.D0/27.D0)&
	+MD*D0000*(40.D0/9.D0)+t*DED00002*(-13.D0/27.D0)+t*DED00001*(-16.D0/9.D0)+t*D00002*(14.D0/9.D0)+t*D00001*(16.D0/3.D0)&
	+t*MD*DED0000*(20.D0/27.D0)+t*MD*D0022*(-4.D0/9.D0)+t*MD*D0013*(-8.D0/9.D0)+t*MD*D0012*(-8.D0/9.D0)+t*MD*D0000*(-28.D0/9.D0)&
	+t**2*DED00002*(-4.D0/27.D0)+t**2*D00002*(11.D0/9.D0)+t**2*MD*D0013*(2.D0/9.D0)+t**2*MD*D0012*(-1.D0/9.D0)+x*DED00002*(44.D0/27.D0)&
	+x*D00002*(-40.D0/9.D0)+x*t*DED00002*(-16.D0/27.D0)+x*t*D00002*(44.D0/9.D0)+x*t*MD*D0022*(-10.D0/9.D0)&
	+x*t*MD*D0012*(-4.D0/3.D0)+x**2*MD*D0022*(-8.D0/9.D0)&
!!$
        +DEC0022*(4.D0/27.D0)+DEC0012*(8.D0/27.D0)+DEC0011*(4.D0/27.D0)+DEC0000*(4.D0/27.D0)+C0022*(-8.D0/9.D0)&
	+C0012*(-16.D0/9.D0)+C0011*(-8.D0/9.D0)+C0000*(-8.D0/9.D0)+MD*DEC002*(-4.D0/27.D0)+MD*DEC001*(-4.D0/27.D0)&
	+MD*C002*(8.D0/9.D0)+MD*C001*(8.D0/9.D0)+t*DEC0022*(1.D0/27.D0)+t*DEC0012*(1.D0/27.D0)+t*C0022*(-2.D0/9.D0)&
	+t*C0012*(-2.D0/9.D0)+t*MD*DEC002*(-1.D0/27.D0)+t*MD*C002*(-1.D0/9.D0)+x*DEC0022*(8.D0/27.D0)+x*DEC0012*(8.D0/27.D0)&
	+x*C0022*(-16.D0/9.D0)+x*C0012*(-16.D0/9.D0)+x*MD*DEC002*(-4.D0/27.D0)+x*MD*C002*(8.D0/9.D0)+x*t*DEC0022*(1.D0/27.D0)&
	+x*t*C0022*(-2.D0/9.D0)+x**2*DEC0022*(4.D0/27.D0)+x**2*C0022*(-8.D0/9.D0)&
!!$
        +DEBf111*(1.D0/27.D0)+DEBf001*(11.D0/27.D0)+Bf111*(-2.D0/9.D0)+Bf001*(-10.D0/9.D0)+MD*DEBf11*(-1.D0/27.D0)+MD*DEBf00*(-1.D0/3.D0)&
	+MD*Bf11*(2.D0/9.D0)+MD*Bf00*(2.D0/3.D0)+t*DEBf111*(1.D0/108.D0)+t*DEBf001*(11.D0/108.D0)+t*Bf111*(-1.D0/18.D0)+t*Bf001*(-5.D0/18.D0)&
	+x*DEBf111*(1.D0/9.D0)+x*DEBf001*(11.D0/27.D0)+x*Bf111*(-2.D0/3.D0)+x*Bf001*(-10.D0/9.D0)+x*MD*DEBf11*(-2.D0/27.D0)+x*MD*Bf11*(4.D0/9.D0)&
	+x*t*DEBf111*(1.D0/54.D0)+x*t*Bf111*(-1.D0/9.D0)+x**2*DEBf111*(2.D0/27.D0)+x**2*Bf111*(-4.D0/9.D0))*(2.D0*x)**2 &
!!$
        +(+t*DEBf001o*(-1.D0/108.D0)+t*DEBf001*(1.D0/108.D0)+t*Bf001o*(1.D0/18.D0)+t*Bf001*(-1.D0/18.D0)+t*DEC0022o*(1.D0/27.D0)&
        +t*DEC0012o*(2.D0/27.D0)+t*DEC0011o*(1.D0/27.D0)+t*DEC0000o*(-5.D0/54.D0)+t*DEC0022*(-1.D0/27.D0)+t*DEC0012*(-2.D0/27.D0)&
        +t*DEC0011*(-1.D0/27.D0)+t*DEC0000*(5.D0/54.D0)+t*C0022o*(1.D0/9.D0)+t*C0012o*(2.D0/9.D0)+t*C0011o*(1.D0/9.D0)&
        +t*C0000o*(11.D0/9.D0)+t*C0022*(-1.D0/9.D0)+t*C0012*(-2.D0/9.D0)+t*C0011*(-1.D0/9.D0)+t*C0000*(-11.D0/9.D0)&
        +t*MD*DEBf00o*(1.D0/108.D0)+t*MD*DEBf00*(-1.D0/108.D0)+t*MD*Bf00o*(-1.D0/18.D0)+t*MD*Bf00*(1.D0/18.D0)&
        +t*MD*DEC002o*(-1.D0/27.D0)+t*MD*DEC001o*(-1.D0/27.D0)+t*MD*DEC002*(1.D0/27.D0)+t*MD*DEC001*(1.D0/27.D0)&
        +t*MD*C002o*(-1.D0/9.D0)+t*MD*C001o*(-1.D0/9.D0)+t*MD*C002*(1.D0/9.D0)+t*MD*C001*(1.D0/9.D0)+x*DEBf001o*(19.D0/27.D0)&
        +x*DEBf001*(1.D0/27.D0)+x*Bf001o*(-14.D0/9.D0)+x*Bf001*(-2.D0/9.D0)+x*DEC0022o*(-2.D0/3.D0)+x*DEC002o*(-2.D0/3.D0)&
        +x*DEC0012o*(-4.D0/3.D0)+x*DEC0011o*(-2.D0/3.D0)+x*DEC001o*(-2.D0/3.D0)+x*DEC0000o*(-2.D0/3.D0)+x*DEC0000*(22.D0/27.D0)&
        +x*C0022o*(4.D0/3.D0)+x*C002o*(4.D0/3.D0)+x*C0012o*(8.D0/3.D0)+x*C0011o*(4.D0/3.D0)+x*C001o*(4.D0/3.D0)&
        +x*C0000o*(4.D0/3.D0)+x*C0000*(-20.D0/9.D0)+x*MD*DEBf00o*(-19.D0/27.D0)+x*MD*DEBf00*(1.D0/27.D0)+x*MD*Bf00o*(14.D0/9.D0)&
        +x*MD*Bf00*(-2.D0/9.D0)+x*MD*DEC002o*(2.D0/3.D0)+x*MD*DEC001o*(2.D0/3.D0)+x*MD*DEC00o*(2.D0/3.D0)+x*MD*C002o*(-4.D0/3.D0)&
        +x*MD*C001o*(-4.D0/3.D0)+x*MD*C00o*(-4.D0/3.D0)+x*t*DEBf001*(1.D0/54.D0)+x*t*Bf001*(-1.D0/9.D0)+x*t*DEC0022*(-5.D0/54.D0)&
        +x*t*DEC0012*(-1.D0/9.D0)+x*t*DEC0011*(-1.D0/54.D0)+x*t*DEC0000*(4.D0/27.D0)+x*t*C0022*(-5.D0/18.D0)+x*t*C0012*(-1.D0/3.D0)&
        +x*t*C0011*(-1.D0/18.D0)+x*t*C0000*(-8.D0/9.D0)+x*t*MD*DEC002*(1.D0/18.D0)+x*t*MD*DEC001*(1.D0/54.D0)+x*t*MD*C002*(1.D0/6.D0)&
        +x*t*MD*C001*(1.D0/18.D0)+x**2*DEBf001*(2.D0/27.D0)+x**2*Bf001*(-4.D0/9.D0)+x**2*DEC0022*(-2.D0/27.D0)+x**2*DEC0012*(-4.D0/27.D0)&
        +x**2*DEC0011*(-2.D0/27.D0)+x**2*DEC0000*(16.D0/27.D0)+x**2*C0022*(-2.D0/9.D0)+x**2*C0012*(-4.D0/9.D0)+x**2*C0011*(-2.D0/9.D0)&
        +x**2*C0000*(-32.D0/9.D0)+x**2*MD*DEC002*(2.D0/27.D0)+x**2*MD*DEC001*(2.D0/27.D0)+x**2*MD*C002*(2.D0/9.D0)&
        +x**2*MD*C001*(2.D0/9.D0)+x**2*t*DEC0022*(-1.D0/27.D0)+x**2*t*DEC0012*(-1.D0/54.D0)+x**2*t*C0022*(-1.D0/9.D0)&
        +x**2*t*C0012*(-2.D0/9.D0)+x**3*DEC0022*(-4.D0/27.D0)+x**3*DEC0012*(-2.D0/27.D0)+x**3*C0022*(-4.D0/9.D0)&
        +x**3*C0012*(-8.D0/9.D0))*(2.D0*x)*(-3.D0/2.D0)*4.D0/3.D0&
!!$
        +(+t*DEBf001o*(1.D0/4.D0)+t*DEBf001*(-1.D0/4.D0)+t*Bf001o*(-1.D0/2.D0)+t*Bf001*(1.D0/2.D0)+t*MD*DEBf00o*(-1.D0/4.D0)&
        +t*MD*DEBf00*(1.D0/4.D0)+t*MD*Bf00o*(1.D0/2.D0)+t*MD*Bf00*(-1.D0/2.D0)+x*t*DEC0022o*(-1.D0/2.D0)+x*t*DEC002o*(-1.D0/2.D0)&
        +x*t*DEC0012o*(-1.D0)+x*t*DEC0011o*(-1.D0/2.D0)+x*t*DEC001o*(-1.D0/2.D0)+x*t*DEC0000o*(-1.D0/2.D0)+x*t*C0022o*(1.D0)&
        +x*t*C002o*(1.D0)+x*t*C0012o*(2.D0)+x*t*C0011o*(1.D0)+x*t*C001o*(1.D0)+x*t*C0000o*(1.D0)+x*t*DEBf001o*(1.D0/2.D0)&
        +x*t*DEBf001*(-3.D0/4.D0)+x*t*Bf001o*(-1.D0)+x*t*Bf001*(3.D0/2.D0)+x*t*MD*DEC002o*(1.D0/2.D0)+x*t*MD*DEC001o*(1.D0/2.D0)&
        +x*t*MD*DEC00o*(1.D0/2.D0)+x*t*MD*C002o*(-1.D0)+x*t*MD*C001o*(-1.D0)+x*t*MD*C00o*(-1.D0)+x*t*MD*DEBf00o*(-1.D0/2.D0)&
        +x*t*MD*DEBf00*(1.D0/2.D0)+x*t*MD*Bf00o*(1.D0)+x*t*MD*Bf00*(-1.D0)+x**2*DEC0022o*(-2.D0)+x**2*DEC002o*(-2.D0)&
        +x**2*DEC0012o*(-4.D0)+x**2*DEC0011o*(-2.D0)+x**2*DEC001o*(-2.D0)+x**2*DEC0000o*(-2.D0)+x**2*C0022o*(4.D0)&
        +x**2*C002o*(4.D0)+x**2*C0012o*(8.D0)+x**2*C0011o*(4.D0)+x**2*C001o*(4.D0)+x**2*C0000o*(4.D0)+x**2*DEBf001o*(2.D0)&
        +x**2*DEBf001*(-1.D0)+x**2*Bf001o*(-4.D0)+x**2*Bf001*(2.D0)+x**2*MD*DEC002o*(2.D0)+x**2*MD*DEC001o*(2.D0)+x**2*MD*DEC00o*(2.D0)&
        +x**2*MD*C002o*(-4.D0)+x**2*MD*C001o*(-4.D0)+x**2*MD*C00o*(-4.D0)+x**2*MD*DEBf00o*(-2.D0)+x**2*MD*Bf00o*(4.D0)&
        +x**2*t*DEBf001*(-1.D0/2.D0)+x**2*t*Bf001*(1.D0)+x**3*DEBf001*(-2.D0)+x**3*Bf001*(4.D0))*(1.D0/2.D0)*(4.D0/3.D0)
        


   A(2) = A(2) &
	+(+DED00002*(2.D0/27.D0)+DED00001*(4.D0/27.D0)+D00222*(14.D0/9.D0)+D0022*(20.D0/9.D0)+D002*(2.D0/3.D0)&
	+D00002*(-4.D0/9.D0)+D00123*(28.D0/9.D0)+D00122*(28.D0/9.D0)+D0012*(40.D0/9.D0)+D00113*(28.D0/9.D0)&
	+D00112*(-28.D0/9.D0)+D00111*(-28.D0/9.D0)+D001*(4.D0/3.D0)+D00001*(-8.D0/9.D0)+MD*DED0000*(2.D0/27.D0)&
	+MD*D0022*(-14.D0/9.D0)+MD*D002*(-20.D0/9.D0)+MD*D0013*(-28.D0/9.D0)+MD*D0011*(28.D0/9.D0)&
	+MD*D00*(-2.D0/3.D0)+MD*D0000*(8.D0/9.D0)+t*DED00002*(1.D0/54.D0)+t*D22222*(-1.D0/9.D0)+t*D2222*(-2.D0/9.D0)&
	+t*D00222*(-1.D0/2.D0)+t*D222*(-1.D0/9.D0)+t*D0022*(-7.D0/9.D0)+t*D002*(-5.D0/18.D0)+t*D00002*(-1.D0/9.D0)&
	+t*D00123*(-5.D0/9.D0)+t*D123*(-4.D0/9.D0)+t*D1223*(-8.D0/9.D0)+t*D12223*(-2.D0/3.D0)+t*D12222*(-4.D0/9.D0)&
	+t*D1222*(-8.D0/9.D0)+t*D00122*(-16.D0/9.D0)+t*D122*(-4.D0/9.D0)+t*D0012*(-8.D0/3.D0)+t*D00113*(-4.D0/3.D0)&
	+t*D113*(-4.D0/9.D0)+t*D11233*(-2.D0/3.D0)+t*D1123*(-8.D0/9.D0)+t*D11223*(-8.D0/9.D0)+t*D1122*(-4.D0/9.D0)&
	+t*D00112*(5.D0/9.D0)+t*D112*(-2.D0/9.D0)+t*D11133*(-4.D0/9.D0)+t*D11123*(4.D0/9.D0)+t*D11122*(4.D0/9.D0)&
	+t*D11113*(4.D0/9.D0)+t*D11112*(2.D0/9.D0)+t*D00111*(4.D0/3.D0)+t*D001*(-8.D0/9.D0)+t*MD*D2222*(1.D0/9.D0)&
	+t*MD*D222*(2.D0/9.D0)+t*MD*D0022*(10.D0/9.D0)+t*MD*D22*(1.D0/9.D0)+t*MD*D002*(14.D0/9.D0)+t*MD*D0013*(4.D0/3.D0)&
	+t*MD*D13*(2.D0/9.D0)+t*MD*D123*(4.D0/9.D0)+t*MD*D1223*(4.D0/9.D0)+t*MD*D1222*(2.D0/9.D0)+t*MD*D122*(4.D0/9.D0)&
	+t*MD*D12*(2.D0/9.D0)+t*MD*D1133*(4.D0/9.D0)+t*MD*D1123*(2.D0/9.D0)+t*MD*D1122*(-2.D0/9.D0)+t*MD*D1113*(-4.D0/9.D0)&
	+t*MD*D1112*(-2.D0/9.D0)+t*MD*D0011*(-4.D0/3.D0)+t*MD*D00*(4.D0/9.D0)+t**2*D22222*(-1.D0/24.D0)+t**2*D2222*(-1.D0/12.D0)&
	+t**2*D00222*(-7.D0/18.D0)+t**2*D222*(-1.D0/24.D0)+t**2*D0022*(-11.D0/18.D0)+t**2*D002*(-2.D0/9.D0)&
	+t**2*D00123*(-11.D0/18.D0)+t**2*D123*(-1.D0/36.D0)+t**2*D1223*(-1.D0/18.D0)+t**2*D12223*(-1.D0/9.D0)&
	+t**2*D12222*(-1.D0/9.D0)+t**2*D1222*(-2.D0/9.D0)+t**2*D122*(-1.D0/9.D0)+t**2*D113*(1.D0/9.D0)+t**2*D11233*(-1.D0/36.D0)&
	+t**2*D1123*(2.D0/9.D0)+t**2*D11222*(1.D0/18.D0)+t**2*D1122*(-1.D0/18.D0)+t**2*D00112*(11.D0/18.D0)+t**2*D112*(-1.D0/36.D0)&
	+t**2*D11133*(1.D0/9.D0)+t**2*D11122*(1.D0/9.D0)+t**2*D11113*(-1.D0/9.D0)+t**2*D11112*(1.D0/36.D0)+t**2*MD*D13*(-1.D0/18.D0)&
	+t**2*MD*D123*(-1.D0/9.D0)+t**2*MD*D1223*(-1.D0/18.D0)+t**2*MD*D1222*(1.D0/36.D0)+t**2*MD*D122*(1.D0/18.D0)&
	+t**2*MD*D12*(1.D0/36.D0)+t**2*MD*D1133*(-1.D0/9.D0)+t**2*MD*D1123*(1.D0/36.D0)+t**2*MD*D1113*(1.D0/9.D0)&
	+t**2*MD*D1112*(-1.D0/36.D0)+t**3*D123*(1.D0/36.D0)+t**3*D1223*(1.D0/18.D0)+t**3*D12223*(1.D0/36.D0)+t**3*D11233*(1.D0/18.D0)&
	+t**3*D11123*(-1.D0/18.D0)+x*DED00002*(2.D0/27.D0)+x*D00222*(14.D0/9.D0)+x*D0022*(20.D0/9.D0)+x*D002*(2.D0/3.D0)&
	+x*D00002*(-4.D0/9.D0)+x*D00123*(28.D0/9.D0)+x*D00112*(-28.D0/9.D0)+x*t*D22222*(-7.D0/18.D0)+x*t*D2222*(-7.D0/9.D0)&
	+x*t*D00222*(-14.D0/9.D0)+x*t*D222*(-7.D0/18.D0)+x*t*D0022*(-22.D0/9.D0)+x*t*D002*(-8.D0/9.D0)&
	+x*t*D00123*(-22.D0/9.D0)+x*t*D123*(-5.D0/9.D0)+x*t*D1223*(-10.D0/9.D0)+x*t*D12223*(-4.D0/3.D0)+x*t*D12222*(-10.D0/9.D0)&
	+x*t*D1222*(-20.D0/9.D0)+x*t*D122*(-10.D0/9.D0)+x*t*D11233*(-7.D0/9.D0)+x*t*D11223*(-10.D0/9.D0)&
	+x*t*D11222*(4.D0/9.D0)+x*t*D1122*(-2.D0/3.D0)+x*t*D00112*(22.D0/9.D0)+x*t*D112*(-1.D0/3.D0)+x*t*D11123*(4.D0/9.D0)&
	+x*t*D11122*(10.D0/9.D0)+x*t*D11112*(1.D0/3.D0)+x*t*MD*D2222*(5.D0/18.D0)+x*t*MD*D222*(5.D0/9.D0)+x*t*MD*D22*(5.D0/18.D0)&
	+x*t*MD*D1223*(5.D0/9.D0)+x*t*MD*D1222*(1.D0/3.D0)+x*t*MD*D122*(2.D0/3.D0)+x*t*MD*D12*(1.D0/3.D0)+x*t*MD*D1123*(1.D0/3.D0)&
	+x*t*MD*D1122*(-5.D0/9.D0)+x*t*MD*D1112*(-1.D0/3.D0)+x*t**2*D22222*(-1.D0/12.D0)+x*t**2*D2222*(-1.D0/6.D0)&
	+x*t**2*D222*(-1.D0/12.D0)+x*t**2*D123*(1.D0/9.D0)+x*t**2*D1223*(2.D0/9.D0)+x*t**2*D12223*(-1.D0/18.D0)&
	+x*t**2*D12222*(-1.D0/9.D0)+x*t**2*D1222*(-2.D0/9.D0)+x*t**2*D122*(-1.D0/9.D0)+x*t**2*D11233*(2.D0/9.D0)&
	+x*t**2*D11223*(-1.D0/9.D0)+x*t**2*D11222*(1.D0/6.D0)+x*t**2*D11123*(-2.D0/9.D0)+x*t**2*D11122*(1.D0/9.D0)&
	+x**2*D22222*(-2.D0/9.D0)+x**2*D2222*(-4.D0/9.D0)+x**2*D222*(-2.D0/9.D0)+x**2*D12223*(-4.D0/9.D0)&
	+x**2*D12222*(-4.D0/9.D0)+x**2*D1222*(-8.D0/9.D0)+x**2*D122*(-4.D0/9.D0)+x**2*D11223*(-4.D0/9.D0)&
	+x**2*D11222*(4.D0/9.D0)+x**2*D11122*(4.D0/9.D0)+x**2*MD*D2222*(2.D0/9.D0)+x**2*MD*D222*(4.D0/9.D0)+x**2*MD*D22*(2.D0/9.D0)&
	+x**2*MD*D1223*(4.D0/9.D0)+x**2*MD*D1122*(-4.D0/9.D0)+x**2*t*D22222*(-7.D0/18.D0)+x**2*t*D2222*(-7.D0/9.D0)&
	+x**2*t*D222*(-7.D0/18.D0)+x**2*t*D12223*(-7.D0/9.D0)+x**2*t*D12222*(-4.D0/9.D0)+x**2*t*D1222*(-8.D0/9.D0)&
	+x**2*t*D122*(-4.D0/9.D0)+x**2*t*D11223*(-4.D0/9.D0)+x**2*t*D11222*(7.D0/9.D0)+x**2*t*D11122*(4.D0/9.D0)&
	+x**3*D22222*(-2.D0/9.D0)+x**3*D2222*(-4.D0/9.D0)+x**3*D222*(-2.D0/9.D0)+x**3*D12223*(-4.D0/9.D0)+x**3*D11222*(4.D0/9.D0)&
!!$
        +DEC0022*(-5.D0/9.D0)+DEC002*(-4.D0/9.D0)+DEC0012*(-26.D0/27.D0)+DEC0011*(-11.D0/27.D0)+DEC001*(-10.D0/27.D0)&
	+C2222*(4.D0/9.D0)+C222*(4.D0/9.D0)+C0022*(2.D0)+C002*(4.D0/3.D0)+C1222*(4.D0/3.D0)+C122*(8.D0/9.D0)+C0012*(28.D0/9.D0)&
	+C1122*(4.D0/3.D0)+C112*(4.D0/9.D0)+C1112*(4.D0/9.D0)+C0011*(10.D0/9.D0)+C001*(8.D0/9.D0)+MD*DEC002*(11.D0/27.D0)&
	+MD*DEC001*(11.D0/27.D0)+MD*DEC00*(10.D0/27.D0)+MD*C222*(-4.D0/9.D0)+MD*C22*(-4.D0/9.D0)+MD*C002*(-16.D0/9.D0)&
	+MD*C122*(-8.D0/9.D0)+MD*C12*(-4.D0/9.D0)+MD*C112*(-4.D0/9.D0)+MD*C001*(-10.D0/9.D0)+MD*C00*(-8.D0/9.D0)+t*DEC0022*(-1.D0/12.D0)&
	+t*DEC002*(-2.D0/27.D0)+t*DEC0012*(-7.D0/108.D0)+t*C2222*(5.D0/36.D0)+t*C222*(5.D0/36.D0)+t*C0022*(2.D0/3.D0)&
	+t*C002*(4.D0/9.D0)+t*C1222*(11.D0/36.D0)+t*C122*(1.D0/6.D0)+t*C0012*(5.D0/9.D0)+t*C1122*(7.D0/36.D0)+t*C112*(1.D0/36.D0)&
	+t*C1112*(1.D0/36.D0)+t*MD*C222*(1.D0/36.D0)+t*MD*C22*(1.D0/36.D0)+t*MD*C12*(-1.D0/36.D0)+t*MD*C112*(-1.D0/36.D0)&
	+x*DEC0022*(-1.D0/3.D0)+x*DEC002*(-8.D0/27.D0)+x*DEC0012*(-7.D0/27.D0)+x*C2222*(11.D0/9.D0)+x*C222*(11.D0/9.D0)&
	+x*C0022*(8.D0/3.D0)+x*C002*(16.D0/9.D0)+x*C1222*(25.D0/9.D0)+x*C122*(14.D0/9.D0)+x*C0012*(20.D0/9.D0)+x*C1122*(17.D0/9.D0)&
	+x*C112*(1.D0/3.D0)+x*C1112*(1.D0/3.D0)+x*MD*C222*(-7.D0/9.D0)+x*MD*C22*(-7.D0/9.D0)+x*MD*C122*(-10.D0/9.D0)&
	+x*MD*C12*(-1.D0/3.D0)+x*MD*C112*(-1.D0/3.D0)+x*t*C2222*(2.D0/9.D0)+x*t*C222*(2.D0/9.D0)+x*t*C1222*(1.D0/3.D0)&
	+x*t*C122*(1.D0/9.D0)+x*t*C1122*(1.D0/9.D0)+x**2*C2222*(8.D0/9.D0)+x**2*C222*(8.D0/9.D0)+x**2*C1222*(4.D0/3.D0)&
	+x**2*C122*(4.D0/9.D0)+x**2*C1122*(4.D0/9.D0)&
!!$
        +DEBf111*(1.D0/108.D0)+Bf111*(-1.D0/18.D0)+MD*DEBf11*(-1.D0/36.D0)+t*DEBf111*(1.D0/432.D0)+t*Bf111*(-1.D0/72.D0)&
	+x*DEBf111*(1.D0/108.D0)+x*Bf111*(-1.D0/18.D0))*(2.D0*x)**2&
!!$
        +(+DEBf001o*(5.D0/27.D0)+DEBf001*(-5.D0/27.D0)+Bf001o*(-4.D0/9.D0)+Bf001*(4.D0/9.D0)+DEC0022o*(-11.D0/27.D0)+DEC002o*(-1.D0/3.D0)&
        +DEC0012o*(-22.D0/27.D0)+DEC0011o*(-11.D0/27.D0)+DEC001o*(-1.D0/3.D0)+DEC0000o*(-4.D0/27.D0)+DEC0022*(11.D0/27.D0)+DEC002*(1.D0/3.D0)&
        +DEC0012*(22.D0/27.D0)+DEC0011*(11.D0/27.D0)+DEC001*(1.D0/3.D0)+DEC0000*(4.D0/27.D0)+C0022o*(4.D0/9.D0)+C002o*(2.D0/3.D0)&
        +C0012o*(8.D0/9.D0)+C0011o*(4.D0/9.D0)+C001o*(2.D0/3.D0)+C0000o*(-16.D0/9.D0)+C0022*(-4.D0/9.D0)+C002*(-2.D0/3.D0)+C0012*(-8.D0/9.D0)&
        +C0011*(-4.D0/9.D0)+C001*(-2.D0/3.D0)+C0000*(16.D0/9.D0)+MD*DEBf00o*(-5.D0/27.D0)+MD*DEBf00*(5.D0/27.D0)+MD*Bf00o*(4.D0/9.D0)&
        +MD*Bf00*(-4.D0/9.D0)+MD*DEC002o*(11.D0/27.D0)+MD*DEC001o*(11.D0/27.D0)+MD*DEC00o*(1.D0/3.D0)+MD*DEC002*(-11.D0/27.D0)&
        +MD*DEC001*(-11.D0/27.D0)+MD*DEC00*(-1.D0/3.D0)+MD*C002o*(-4.D0/9.D0)+MD*C001o*(-4.D0/9.D0)+MD*C00o*(-2.D0/3.D0)+MD*C002*(4.D0/9.D0)&
        +MD*C001*(4.D0/9.D0)+MD*C00*(2.D0/3.D0)+x*DEC0022*(11.D0/27.D0)+x*DEC002*(8.D0/27.D0)+x*DEC0012*(22.D0/27.D0)+x*DEC0011*(11.D0/27.D0)&
        +x*DEC001*(10.D0/27.D0)+x*C0022*(2.D0/9.D0)+x*C002*(-4.D0/9.D0)+x*C0012*(-8.D0/9.D0)+x*C0011*(-10.D0/9.D0)+x*C001*(-8.D0/9.D0)&
        +x*MD*DEC002*(-11.D0/27.D0)+x*MD*DEC001*(-11.D0/27.D0)+x*MD*DEC00*(-10.D0/27.D0)+x*MD*C002*(4.D0/9.D0)+x*MD*C001*(10.D0/9.D0)&
        +x*MD*C00*(8.D0/9.D0)+x**2*C1222*(-2.D0/9.D0)+x**2*C122*(-2.D0/9.D0)+x**2*C1122*(-4.D0/9.D0)+x**2*C112*(-2.D0/9.D0)+x**2*C1112*(-2.D0/9.D0)&
        +x**2*MD*C122*(2.D0/9.D0)+x**2*MD*C12*(2.D0/9.D0)+x**2*MD*C112*(2.D0/9.D0))*(2.D0*x)*(-3.D0/2.D0)*4.D0/3.D0&
!!$
        +(0.D0)*(1.D0/2.D0)*(4.D0/3.D0)


   A(3) = A(3) &
        +(0.D0&
!!$
        +DEC0022*(4.D0/27.D0)+DEC0012*(8.D0/27.D0)+DEC0011*(4.D0/27.D0)+DEC0000*(4.D0/9.D0)+C0022*(4.D0/9.D0)&
 	+C0012*(8.D0/9.D0)+C0011*(4.D0/9.D0)+C0000*(8.D0/3.D0)+MD*DEC002*(-4.D0/27.D0)+MD*DEC001*(-4.D0/27.D0)&
	+MD*C002*(-4.D0/9.D0)+MD*C001*(-4.D0/9.D0)+x*DEC0022*(8.D0/27.D0)+x*DEC0012*(8.D0/27.D0)+x*C0022*(8.D0/9.D0)&
	+x*C0012*(8.D0/9.D0)+x*MD*DEC002*(-4.D0/27.D0)+x*MD*C002*(-4.D0/9.D0)&
!!$
	+DEBf111*(1.D0/27.D0)+DEBf001*(1.D0/9.D0)+Bf111*(1.D0/9.D0)+Bf001*(2.D0/3.D0)+MD*DEBf11*(-1.D0/27.D0)+MD*DEBf00*(1.D0/27.D0)&
	+MD*Bf11*(-1.D0/9.D0)+MD*Bf00*(-2.D0/9.D0)+t*DEBf111*(1.D0/108.D0)+t*DEBf001*(1.D0/36.D0)+t*Bf111*(1.D0/36.D0)+t*Bf001*(1.D0/6.D0)&
	+x*DEBf111*(1.D0/9.D0)+x*DEBf001*(1.D0/9.D0)+x*Bf111*(1.D0/3.D0)+x*Bf001*(2.D0/3.D0)+x*MD*DEBf11*(-2.D0/27.D0)+x*MD*Bf11*(-2.D0/9.D0)&
	+x*t*DEBf111*(1.D0/54.D0)+x*t*Bf111*(1.D0/18.D0)+x**2*DEBf111*(2.D0/27.D0)+x**2*Bf111*(2.D0/9.D0))*(2.D0*x)**2&
!!$
        +(+t*DEBf001o*(-1.D0/108.D0)+t*DEBf001*(1.D0/108.D0)+t*Bf001o*(1.D0/18.D0)+t*Bf001*(-1.D0/18.D0)+t*DEC0022o*(1.D0/27.D0)&
        +t*DEC0012o*(2.D0/27.D0)+t*DEC0011o*(1.D0/27.D0)+t*DEC0000o*(-5.D0/54.D0)+t*DEC0022*(-1.D0/27.D0)+t*DEC0012*(-2.D0/27.D0)&
        +t*DEC0011*(-1.D0/27.D0)+t*DEC0000*(5.D0/54.D0)+t*C0022o*(1.D0/9.D0)+t*C0012o*(2.D0/9.D0)+t*C0011o*(1.D0/9.D0)&
        +t*C0000o*(11.D0/9.D0)+t*C0022*(-1.D0/9.D0)+t*C0012*(-2.D0/9.D0)+t*C0011*(-1.D0/9.D0)+t*C0000*(-11.D0/9.D0)&
        +t*MD*DEBf00o*(1.D0/108.D0)+t*MD*DEBf00*(-1.D0/108.D0)+t*MD*Bf00o*(-1.D0/18.D0)+t*MD*Bf00*(1.D0/18.D0)&
        +t*MD*DEC002o*(-1.D0/27.D0)+t*MD*DEC001o*(-1.D0/27.D0)+t*MD*DEC002*(1.D0/27.D0)+t*MD*DEC001*(1.D0/27.D0)&
        +t*MD*C002o*(-1.D0/9.D0)+t*MD*C001o*(-1.D0/9.D0)+t*MD*C002*(1.D0/9.D0)+t*MD*C001*(1.D0/9.D0)+x*DEBf001o*(19.D0/27.D0)&
        +x*DEBf001*(1.D0/27.D0)+x*Bf001o*(-14.D0/9.D0)+x*Bf001*(-2.D0/9.D0)+x*DEC0022o*(-2.D0/3.D0)+x*DEC002o*(-2.D0/3.D0)&
        +x*DEC0012o*(-4.D0/3.D0)+x*DEC0011o*(-2.D0/3.D0)+x*DEC001o*(-2.D0/3.D0)+x*DEC0000o*(-2.D0/3.D0)+x*DEC0000*(22.D0/27.D0)&
        +x*C0022o*(4.D0/3.D0)+x*C002o*(4.D0/3.D0)+x*C0012o*(8.D0/3.D0)+x*C0011o*(4.D0/3.D0)+x*C001o*(4.D0/3.D0)&
        +x*C0000o*(4.D0/3.D0)+x*C0000*(-20.D0/9.D0)+x*MD*DEBf00o*(-19.D0/27.D0)+x*MD*DEBf00*(1.D0/27.D0)+x*MD*Bf00o*(14.D0/9.D0)&
        +x*MD*Bf00*(-2.D0/9.D0)+x*MD*DEC002o*(2.D0/3.D0)+x*MD*DEC001o*(2.D0/3.D0)+x*MD*DEC00o*(2.D0/3.D0)+x*MD*C002o*(-4.D0/3.D0)&
        +x*MD*C001o*(-4.D0/3.D0)+x*MD*C00o*(-4.D0/3.D0)+x*t*DEBf001*(1.D0/54.D0)+x*t*Bf001*(-1.D0/9.D0)+x*t*DEC0022*(-5.D0/54.D0)&
        +x*t*DEC0012*(-1.D0/9.D0)+x*t*DEC0011*(-1.D0/54.D0)+x*t*DEC0000*(4.D0/27.D0)+x*t*C0022*(-5.D0/18.D0)+x*t*C0012*(-1.D0/3.D0)&
        +x*t*C0011*(-1.D0/18.D0)+x*t*C0000*(-8.D0/9.D0)+x*t*MD*DEC002*(1.D0/18.D0)+x*t*MD*DEC001*(1.D0/54.D0)+x*t*MD*C002*(1.D0/6.D0)&
        +x*t*MD*C001*(1.D0/18.D0)+x**2*DEBf001*(2.D0/27.D0)+x**2*Bf001*(-4.D0/9.D0)+x**2*DEC0022*(-2.D0/27.D0)+x**2*DEC0012*(-4.D0/27.D0)&
        +x**2*DEC0011*(-2.D0/27.D0)+x**2*DEC0000*(16.D0/27.D0)+x**2*C0022*(-2.D0/9.D0)+x**2*C0012*(-4.D0/9.D0)+x**2*C0011*(-2.D0/9.D0)&
        +x**2*C0000*(-32.D0/9.D0)+x**2*MD*DEC002*(2.D0/27.D0)+x**2*MD*DEC001*(2.D0/27.D0)+x**2*MD*C002*(2.D0/9.D0)&
        +x**2*MD*C001*(2.D0/9.D0)+x**2*t*DEC0022*(-1.D0/27.D0)+x**2*t*DEC0012*(-1.D0/54.D0)+x**2*t*C0022*(-1.D0/9.D0)&
        +x**2*t*C0012*(-2.D0/9.D0)+x**3*DEC0022*(-4.D0/27.D0)+x**3*DEC0012*(-2.D0/27.D0)+x**3*C0022*(-4.D0/9.D0)&
        +x**3*C0012*(-8.D0/9.D0))*(2.D0*x)*(-3.D0/2.D0)*4.D0/3.D0&
!!$
        +(+t*DEBf001o*(1.D0/4.D0)+t*DEBf001*(-1.D0/4.D0)+t*Bf001o*(-1.D0/2.D0)+t*Bf001*(1.D0/2.D0)+t*MD*DEBf00o*(-1.D0/4.D0)&
        +t*MD*DEBf00*(1.D0/4.D0)+t*MD*Bf00o*(1.D0/2.D0)+t*MD*Bf00*(-1.D0/2.D0)+x*t*DEC0022o*(-1.D0/2.D0)+x*t*DEC002o*(-1.D0/2.D0)&
        +x*t*DEC0012o*(-1.D0)+x*t*DEC0011o*(-1.D0/2.D0)+x*t*DEC001o*(-1.D0/2.D0)+x*t*DEC0000o*(-1.D0/2.D0)+x*t*C0022o*(1.D0)&
        +x*t*C002o*(1.D0)+x*t*C0012o*(2.D0)+x*t*C0011o*(1.D0)+x*t*C001o*(1.D0)+x*t*C0000o*(1.D0)+x*t*DEBf001o*(1.D0/2.D0)&
        +x*t*DEBf001*(-3.D0/4.D0)+x*t*Bf001o*(-1.D0)+x*t*Bf001*(3.D0/2.D0)+x*t*MD*DEC002o*(1.D0/2.D0)+x*t*MD*DEC001o*(1.D0/2.D0)&
        +x*t*MD*DEC00o*(1.D0/2.D0)+x*t*MD*C002o*(-1.D0)+x*t*MD*C001o*(-1.D0)+x*t*MD*C00o*(-1.D0)+x*t*MD*DEBf00o*(-1.D0/2.D0)&
        +x*t*MD*DEBf00*(1.D0/2.D0)+x*t*MD*Bf00o*(1.D0)+x*t*MD*Bf00*(-1.D0)+x**2*DEC0022o*(-2.D0)+x**2*DEC002o*(-2.D0)+x**2*DEC0012o*(-4.D0)&
        +x**2*DEC0011o*(-2.D0)+x**2*DEC001o*(-2.D0)+x**2*DEC0000o*(-2.D0)+x**2*C0022o*(4.D0)+x**2*C002o*(4.D0)+x**2*C0012o*(8.D0)&
        +x**2*C0011o*(4.D0)+x**2*C001o*(4.D0)+x**2*C0000o*(4.D0)+x**2*DEBf001o*(2.D0)+x**2*DEBf001*(-1.D0)+x**2*Bf001o*(-4.D0)&
        +x**2*Bf001*(2.D0)+x**2*MD*DEC002o*(2.D0)+x**2*MD*DEC001o*(2.D0)+x**2*MD*DEC00o*(2.D0)+x**2*MD*C002o*(-4.D0)&
        +x**2*MD*C001o*(-4.D0)+x**2*MD*C00o*(-4.D0)+x**2*MD*DEBf00o*(-2.D0)+x**2*MD*Bf00o*(4.D0)+x**2*t*DEBf001*(-1.D0/2.D0)&
        +x**2*t*Bf001*(1.D0)+x**3*DEBf001*(-2.D0)+x**3*Bf001*(4.D0))*(1.D0/2.D0)*(4.D0/3.D0)

   A(4) = A(4) &
	+(+DED00002*(2.D0/9.D0)+D00002*(-20.D0/3.D0)+MD*D0022*(4.D0/9.D0)+MD*D0012*(8.D0/9.D0)&
	+t*DED00002*(-8.D0/27.D0)+t*D00002*(22.D0/9.D0)+t*MD*D0012*(-2.D0/9.D0)+x*MD*D0022*(4.D0/9.D0)&
!!$
        +DEC0022*(2.D0/27.D0)+DEC0012*(2.D0/27.D0)+C0022*(-4.D0/9.D0)+C0012*(-4.D0/9.D0)+MD*DEC002*(-2.D0/27.D0)&
	+MD*C002*(-2.D0/9.D0)+x*DEC0022*(2.D0/27.D0)+x*C0022*(-4.D0/9.D0)&
!!$
        +DEBf111*(1.D0/54.D0)+DEBf001*(11.D0/54.D0)+Bf111*(-1.D0/9.D0)+Bf001*(-5.D0/9.D0)+x*DEBf111*(1.D0/27.D0)&
        +x*Bf111*(-2.D0/9.D0))*(2.D0*x)**2&
!!$
        +(+DEBf001o*(-1.D0/54.D0)+DEBf001*(1.D0/54.D0)+Bf001o*(1.D0/9.D0)+Bf001*(-1.D0/9.D0)+DEC0022o*(2.D0/27.D0)+DEC0012o*(4.D0/27.D0)&
        +DEC0011o*(2.D0/27.D0)+DEC0000o*(-5.D0/27.D0)+DEC0022*(-2.D0/27.D0)+DEC0012*(-4.D0/27.D0)+DEC0011*(-2.D0/27.D0)+DEC0000*(5.D0/27.D0)&
        +C0022o*(2.D0/9.D0)+C0012o*(4.D0/9.D0)+C0011o*(2.D0/9.D0)+C0000o*(22.D0/9.D0)+C0022*(-2.D0/9.D0)+C0012*(-4.D0/9.D0)+C0011*(-2.D0/9.D0)&
        +C0000*(-22.D0/9.D0)+MD*DEBf00o*(1.D0/54.D0)+MD*DEBf00*(-1.D0/54.D0)+MD*Bf00o*(-1.D0/9.D0)+MD*Bf00*(1.D0/9.D0)+MD*DEC002o*(-2.D0/27.D0)&
        +MD*DEC001o*(-2.D0/27.D0)+MD*DEC002*(2.D0/27.D0)+MD*DEC001*(2.D0/27.D0)+MD*C002o*(-2.D0/9.D0)+MD*C001o*(-2.D0/9.D0)+MD*C002*(2.D0/9.D0)&
        +MD*C001*(2.D0/9.D0)+x*DEBf001*(1.D0/27.D0)+x*Bf001*(-2.D0/9.D0)+x*DEC0022*(-5.D0/27.D0)+x*DEC0012*(-2.D0/9.D0)+x*DEC0011*(-1.D0/27.D0)&
        +x*DEC0000*(8.D0/27.D0)+x*C0022*(-5.D0/9.D0)+x*C0012*(-2.D0/3.D0)+x*C0011*(-1.D0/9.D0)+x*C0000*(-16.D0/9.D0)+x*MD*DEC002*(1.D0/9.D0)&
        +x*MD*DEC001*(1.D0/27.D0)+x*MD*C002*(1.D0/3.D0)+x*MD*C001*(1.D0/9.D0)+x**2*DEC0022*(-2.D0/27.D0)+x**2*DEC0012*(-1.D0/27.D0)&
        +x**2*C0022*(-2.D0/9.D0)+x**2*C0012*(-4.D0/9.D0))*(2.D0*x)*(-3.D0/2.D0)*4.D0/3.D0&
!!$
        +(+DEBf001o*(1.D0/2.D0)+DEBf001*(-1.D0/2.D0)+Bf001o*(-1.D0)+Bf001*(1.D0)+MD*DEBf00o*(-1.D0/2.D0)+MD*DEBf00*(1.D0/2.D0)&
        +MD*Bf00o*(1.D0)+MD*Bf00*(-1.D0)+x*DEC0022o*(-1.D0)+x*DEC002o*(-1.D0)+x*DEC0012o*(-2.D0)+x*DEC0011o*(-1.D0)&
        +x*DEC001o*(-1.D0)+x*DEC0000o*(-1.D0)+x*C0022o*(2.D0)+x*C002o*(2.D0)+x*C0012o*(4.D0)+x*C0011o*(2.D0)+x*C001o*(2.D0)&
        +x*C0000o*(2.D0)+x*DEBf001o*(1.D0)+x*DEBf001*(-3.D0/2.D0)+x*Bf001o*(-2.D0)+x*Bf001*(3.D0)+x*MD*DEC002o*(1.D0)&
        +x*MD*DEC001o*(1.D0)+x*MD*DEC00o*(1.D0)+x*MD*C002o*(-2.D0)+x*MD*C001o*(-2.D0)+x*MD*C00o*(-2.D0)+x*MD*DEBf00o*(-1.D0)&
        +x*MD*DEBf00*(1.D0)+x*MD*Bf00o*(2.D0)+x*MD*Bf00*(-2.D0)+x**2*DEBf001*(-1.D0)+x**2*Bf001*(2.D0))*(1.D0/2.D0)*(4.D0/3.D0)

   A(5) = A(5) &
	+(+DED00002*(1.D0/9.D0)+DED0000*(5.D0/27.D0)+D00222*(-2.D0/9.D0)+D0022*(-2.D0/9.D0)+D00002*(-10.D0/3.D0)+D0013*(-4.D0/9.D0)&
	+D00123*(-4.D0/9.D0)+D00122*(-8.D0/9.D0)+D0012*(-8.D0/9.D0)+D00112*(-4.D0/9.D0)+D0011*(-4.D0/9.D0)+D0000*(-22.D0/9.D0)&
	+MD*D0022*(2.D0/9.D0)+MD*D002*(2.D0/9.D0)+MD*D0012*(4.D0/9.D0)+MD*D001*(4.D0/9.D0)+t*DED00002*(-4.D0/27.D0)&
	+t*DED0000*(-4.D0/27.D0)+t*D0022*(1.D0/6.D0)+t*D00002*(11.D0/9.D0)+t*D0013*(5.D0/9.D0)+t*D00123*(2.D0/9.D0)+t*D00122*(2.D0/9.D0)&
	+t*D0012*(4.D0/9.D0)+t*D00112*(2.D0/9.D0)+t*D0011*(1.D0/9.D0)+t*D0000*(8.D0/9.D0)+t*MD*D0012*(-1.D0/9.D0)+t*MD*D001*(-1.D0/9.D0)&
	+t**2*D0013*(-1.D0/9.D0)+t**2*D00123*(-1.D0/18.D0)+x*D00222*(-4.D0/9.D0)+x*D0022*(-4.D0/9.D0)+x*D00122*(-8.D0/9.D0)&
	+x*D0012*(-8.D0/9.D0)+x*MD*D0022*(2.D0/9.D0)+x*MD*D002*(2.D0/9.D0)+x*t*D0022*(1.D0/3.D0)&
	+x*t*D00122*(2.D0/9.D0)+x*t*D0012*(4.D0/9.D0)+x**2*D00222*(-2.D0/9.D0)+x**2*D0022*(2.D0/9.D0)&
!!$
        +DEC0022*(-1.D0/6.D0)+DEC002*(-5.D0/54.D0)+DEC0012*(-5.D0/27.D0)+DEC0011*(-1.D0/54.D0)+DEC001*(1.D0/54.D0)&
	+DEC0000*(1.D0/9.D0)+C2222*(-1.D0/9.D0)+C222*(-1.D0/9.D0)+C0022*(-4.D0/3.D0)+C002*(-7.D0/9.D0)+C1222*(-1.D0/3.D0)&
	+C122*(-2.D0/9.D0)+C0012*(-14.D0/9.D0)+C1122*(-1.D0/3.D0)+C112*(-1.D0/9.D0)+C1112*(-1.D0/9.D0)+C0011*(-2.D0/9.D0)&
	+C001*(-1.D0/9.D0)+C0000*(-4.D0/3.D0)+MD*DEC002*(1.D0/54.D0)+MD*DEC001*(1.D0/54.D0)+MD*DEC00*(-1.D0/54.D0)&
	+MD*C222*(1.D0/9.D0)+MD*C22*(1.D0/9.D0)+MD*C002*(2.D0/9.D0)+MD*C122*(2.D0/9.D0)+MD*C12*(1.D0/9.D0)&
	+MD*C112*(1.D0/9.D0)+MD*C001*(2.D0/9.D0)+MD*C00*(1.D0/9.D0)+x*DEC0022*(-5.D0/54.D0)+x*DEC002*(-1.D0/27.D0)&
	+x*DEC0012*(-7.D0/54.D0)+x*C2222*(-5.D0/18.D0)+x*C222*(-5.D0/18.D0)+x*C0022*(-10.D0/9.D0)+x*C002*(-4.D0/9.D0)&
	+x*C1222*(-11.D0/18.D0)+x*C122*(-1.D0/3.D0)+x*C0012*(-8.D0/9.D0)+x*C1122*(-7.D0/18.D0)+x*C112*(-1.D0/18.D0)&
	+x*C1112*(-1.D0/18.D0)+x*MD*C222*(1.D0/6.D0)+x*MD*C22*(1.D0/6.D0)+x*MD*C122*(2.D0/9.D0)+x*MD*C12*(1.D0/18.D0)&
	+x*MD*C112*(1.D0/18.D0)+x**2*C2222*(-1.D0/9.D0)+x**2*C222*(-1.D0/9.D0)+x**2*C1222*(-2.D0/9.D0)+x**2*C122*(-1.D0/9.D0)&
	+x**2*C1122*(-1.D0/9.D0)&
!!$
        +DEBf111*(-1.D0/108.D0)+DEBf001*(-1.D0/27.D0)+Bf111*(-1.D0/36.D0)+Bf001*(-1.D0/9.D0)+MD*DEBf11*(-1.D0/108.D0)+MD*Bf11*(-1.D0/36.D0)&
	+x*DEBf111*(-1.D0/54.D0)+x*Bf111*(-1.D0/18.D0))*(2.D0*x)**2&
!!$
        +(+DEBf001o*(1.D0/54.D0)+DEBf001*(-1.D0/54.D0)+Bf001o*(-1.D0/9.D0)+Bf001*(1.D0/9.D0)+DEC0022o*(-2.D0/27.D0)+DEC0012o*(-4.D0/27.D0)&
        +DEC0011o*(-2.D0/27.D0)+DEC0000o*(5.D0/27.D0)+DEC0022*(2.D0/27.D0)+DEC0012*(4.D0/27.D0)+DEC0011*(2.D0/27.D0)+DEC0000*(-5.D0/27.D0)&
        +C0022o*(-2.D0/9.D0)+C0012o*(-4.D0/9.D0)+C0011o*(-2.D0/9.D0)+C0000o*(-22.D0/9.D0)+C0022*(2.D0/9.D0)+C0012*(4.D0/9.D0)+C0011*(2.D0/9.D0)&
        +C0000*(22.D0/9.D0)+MD*DEBf00o*(-1.D0/54.D0)+MD*DEBf00*(1.D0/54.D0)+MD*Bf00o*(1.D0/9.D0)+MD*Bf00*(-1.D0/9.D0)+MD*DEC002o*(2.D0/27.D0)&
        +MD*DEC001o*(2.D0/27.D0)+MD*DEC002*(-2.D0/27.D0)+MD*DEC001*(-2.D0/27.D0)+MD*C002o*(2.D0/9.D0)+MD*C001o*(2.D0/9.D0)+MD*C002*(-2.D0/9.D0)&
        +MD*C001*(-2.D0/9.D0)+x*DEBf001*(1.D0/18.D0)+x*DEC0022*(-1.D0/54.D0)+x*DEC002*(-5.D0/27.D0)+x*DEC0012*(1.D0/54.D0)+x*DEC0011*(1.D0/27.D0)&
        +x*DEC001*(1.D0/54.D0)+x*DEC0000*(-1.D0/3.D0)+x*C0022*(17.D0/18.D0)+x*C002*(4.D0/9.D0)+x*C0012*(8.D0/9.D0)+x*C0011*(-1.D0/18.D0)&
        +x*C001*(-1.D0/9.D0)+x*C0000*(4.D0/3.D0)+x*MD*DEC002*(-1.D0/9.D0)+x*MD*DEC001*(-1.D0/27.D0)+x*MD*DEC00*(-1.D0/54.D0)&
        +x*MD*C002*(-1.D0/6.D0)+x*MD*C001*(1.D0/18.D0)+x*MD*C00*(1.D0/9.D0)+x**2*DEC0022*(-7.D0/54.D0)+x**2*DEC002*(-5.D0/27.D0)&
        +x**2*DEC0012*(-11.D0/54.D0)+x**2*C0022*(4.D0/9.D0)+x**2*C002*(4.D0/9.D0)+x**2*C0012*(8.D0/9.D0)+x**3*C1222*(1.D0/9.D0)&
        +x**3*C122*(1.D0/9.D0)+x**3*C1122*(1.D0/9.D0))*(2.D0*x)*(-3.D0/2.D0)*4.D0/3.D0&
!!$
        +(+DEBf001o*(-1.D0/2.D0)+DEBf001*(1.D0/2.D0)+Bf001o*(1.D0)+Bf001*(-1.D0)+MD*DEBf00o*(1.D0/2.D0)+MD*DEBf00*(-1.D0/2.D0)&
        +MD*Bf00o*(-1.D0)+MD*Bf00*(1.D0)+x*DEC0022o*(1.D0)+x*DEC002o*(1.D0)+x*DEC0012o*(2.D0)+x*DEC0011o*(1.D0)+x*DEC001o*(1.D0)&
        +x*DEC0000o*(1.D0)+x*C0022o*(-2.D0)+x*C002o*(-2.D0)+x*C0012o*(-4.D0)+x*C0011o*(-2.D0)+x*C001o*(-2.D0)+x*C0000o*(-2.D0)&
        +x*DEBf001o*(-1.D0)+x*DEBf001*(3.D0/2.D0)+x*Bf001o*(2.D0)+x*Bf001*(-3.D0)+x*MD*DEC002o*(-1.D0)+x*MD*DEC001o*(-1.D0)&
        +x*MD*DEC00o*(-1.D0)+x*MD*C002o*(2.D0)+x*MD*C001o*(2.D0)+x*MD*C00o*(2.D0)+x*MD*DEBf00o*(1.D0)+x*MD*DEBf00*(-1.D0)&
        +x*MD*Bf00o*(-2.D0)+x*MD*Bf00*(2.D0)+x**2*DEBf001*(1.D0)+x**2*Bf001*(-2.D0))*(1.D0/2.D0)*(4.D0/3.D0)

   A(6) = A(6) &
	+(+DED00002*(-1.D0/9.D0)+DED0000*(-5.D0/27.D0)+D00222*(2.D0/9.D0)+D0022*(2.D0/9.D0)+D00002*(10.D0/3.D0)+D0013*(4.D0/9.D0)&
	+D00123*(4.D0/9.D0)+D00122*(8.D0/9.D0)+D0012*(8.D0/9.D0)+D00112*(4.D0/9.D0)+D0011*(4.D0/9.D0)+D0000*(22.D0/9.D0)&
	+MD*D0022*(-2.D0/9.D0)+MD*D002*(-2.D0/9.D0)+MD*D0012*(-4.D0/9.D0)+MD*D001*(-4.D0/9.D0)+t*DED00002*(4.D0/27.D0)+t*DED0000*(4.D0/27.D0)&
	+t*D0022*(-1.D0/6.D0)+t*D00002*(-11.D0/9.D0)+t*D0013*(-5.D0/9.D0)+t*D00123*(-4.D0/9.D0)+t*D00122*(-2.D0/9.D0)+t*D0012*(-4.D0/9.D0)&
	+t*D0011*(-1.D0/9.D0)+t*D0000*(-8.D0/9.D0)+t*MD*D0012*(1.D0/9.D0)+t*MD*D001*(1.D0/9.D0)+t**2*D0013*(1.D0/9.D0)+t**2*D00123*(1.D0/18.D0)&
	+x*D00222*(4.D0/9.D0)+x*D0022*(4.D0/9.D0)+x*D00122*(8.D0/9.D0)+x*D0012*(8.D0/9.D0)+x*MD*D0022*(-2.D0/9.D0)+x*MD*D002*(-2.D0/9.D0)&
	+x*t*D0022*(-1.D0/3.D0)+x*t*D00122*(-2.D0/9.D0)+x*t*D0012*(-4.D0/9.D0)+x**2*D00222*(2.D0/9.D0)+x**2*D0022*(-2.D0/9.D0)&
!!$
        +DEC0022*(-1.D0/54.D0)+DEC002*(-1.D0/54.D0)+DEC0011*(1.D0/54.D0)+DEC001*(-1.D0/54.D0)+DEC0000*(-1.D0/9.D0)&
	+C0022*(4.D0/9.D0)+C002*(1.D0/9.D0)+C0012*(2.D0/3.D0)+C0011*(2.D0/9.D0)+C001*(1.D0/9.D0)+C0000*(4.D0/3.D0)&
	+MD*DEC002*(1.D0/54.D0)+MD*DEC001*(-1.D0/54.D0)+MD*DEC00*(1.D0/54.D0)+MD*C002*(-1.D0/9.D0)+MD*C001*(-2.D0/9.D0)&
	+MD*C00*(-1.D0/9.D0)+x*DEC0022*(5.D0/54.D0)+x*DEC002*(1.D0/27.D0)+x*DEC0012*(7.D0/54.D0)+x*C2222*(1.D0/18.D0)&
	+x*C222*(1.D0/18.D0)+x*C0022*(10.D0/9.D0)+x*C002*(4.D0/9.D0)+x*C1222*(1.D0/6.D0)+x*C122*(1.D0/9.D0)&
	+x*C0012*(8.D0/9.D0)+x*C1122*(1.D0/6.D0)+x*C112*(1.D0/18.D0)+x*C1112*(1.D0/18.D0)+x*MD*C222*(-1.D0/18.D0)&
	+x*MD*C22*(-1.D0/18.D0)+x*MD*C122*(-1.D0/9.D0)+x*MD*C12*(-1.D0/18.D0)+x*MD*C112*(-1.D0/18.D0)+x**2*C2222*(1.D0/9.D0)&
	+x**2*C222*(1.D0/9.D0)+x**2*C1222*(2.D0/9.D0)+x**2*C122*(1.D0/9.D0)+x**2*C1122*(1.D0/9.D0)&
!!$
        +DEBf111*(1.D0/108.D0)+DEBf001*(1.D0/27.D0)+Bf111*(1.D0/36.D0)+Bf001*(1.D0/9.D0)+MD*DEBf11*(-1.D0/108.D0)+MD*Bf11*(-1.D0/36.D0)&
	+x*DEBf111*(1.D0/54.D0)+x*Bf111*(1.D0/18.D0))*(2.D0*x)**2&
!!$
        +(+DEBf001o*(2.D0/27.D0)+DEBf001*(-2.D0/27.D0)+Bf001o*(-1.D0/9.D0)+Bf001*(1.D0/9.D0)+DEC0022o*(-7.D0/54.D0)+DEC002o*(-1.D0/6.D0)&
        +DEC0012o*(-7.D0/27.D0)+DEC0011o*(-7.D0/54.D0)+DEC001o*(-1.D0/6.D0)+DEC0000o*(-7.D0/27.D0)+DEC0022*(7.D0/54.D0)+DEC002*(1.D0/6.D0)&
        +DEC0012*(7.D0/27.D0)+DEC0011*(7.D0/54.D0)+DEC001*(1.D0/6.D0)+DEC0000*(7.D0/27.D0)+C0022o*(4.D0/9.D0)+C002o*(1.D0/3.D0)&
        +C0012o*(8.D0/9.D0)+C0011o*(4.D0/9.D0)+C001o*(1.D0/3.D0)+C0000o*(14.D0/9.D0)+C0022*(-4.D0/9.D0)+C002*(-1.D0/3.D0)&
        +C0012*(-8.D0/9.D0)+C0011*(-4.D0/9.D0)+C001*(-1.D0/3.D0)+C0000*(-14.D0/9.D0)+MD*DEBf00o*(-2.D0/27.D0)+MD*DEBf00*(2.D0/27.D0)&
        +MD*Bf00o*(1.D0/9.D0)+MD*Bf00*(-1.D0/9.D0)+MD*DEC002o*(7.D0/54.D0)+MD*DEC001o*(7.D0/54.D0)+MD*DEC00o*(1.D0/6.D0)+MD*DEC002*(-7.D0/54.D0)&
        +MD*DEC001*(-7.D0/54.D0)+MD*DEC00*(-1.D0/6.D0)+MD*C002o*(-4.D0/9.D0)+MD*C001o*(-4.D0/9.D0)+MD*C00o*(-1.D0/3.D0)+MD*C002*(4.D0/9.D0)&
        +MD*C001*(4.D0/9.D0)+MD*C00*(1.D0/3.D0)+x*DEBf001*(-1.D0/18.D0)+x*DEC0022*(2.D0/9.D0)+x*DEC002*(1.D0/3.D0)+x*DEC0012*(7.D0/18.D0)&
        +x*DEC0011*(1.D0/6.D0)+x*DEC001*(1.D0/6.D0)+x*DEC0000*(1.D0/3.D0)+x*C0022*(-5.D0/6.D0)+x*C002*(-2.D0/3.D0)+x*C0012*(-4.D0/3.D0)&
        +x*C0011*(-1.D0/2.D0)+x*C001*(-1.D0/3.D0)+x*C0000*(-4.D0/3.D0)+x*MD*DEC002*(-5.D0/54.D0)+x*MD*DEC001*(-1.D0/6.D0)+x*MD*DEC00*(-1.D0/6.D0)&
        +x*MD*C002*(7.D0/18.D0)+x*MD*C001*(1.D0/2.D0)+x*MD*C00*(1.D0/3.D0)+x**2*DEC0022*(7.D0/54.D0)+x**2*DEC002*(5.D0/27.D0)+x**2*DEC0012*(11.D0/54.D0)&
        +x**2*C0022*(-4.D0/9.D0)+x**2*C002*(-4.D0/9.D0)+x**2*C1222*(-1.D0/9.D0)+x**2*C122*(-1.D0/9.D0)+x**2*C0012*(-8.D0/9.D0)+x**2*C1122*(-2.D0/9.D0)&
        +x**2*C112*(-1.D0/9.D0)+x**2*C1112*(-1.D0/9.D0)+x**2*MD*C122*(1.D0/9.D0)+x**2*MD*C12*(1.D0/9.D0)+x**2*MD*C112*(1.D0/9.D0)+x**3*C1222*(-1.D0/9.D0)&
        +x**3*C122*(-1.D0/9.D0)+x**3*C1122*(-1.D0/9.D0))*(2.D0*x)*(-3.D0/2.D0)*4.D0/3.D0&
!!$
        +(+DEBf001o*(1.D0/2.D0)+DEBf001*(-1.D0/2.D0)+Bf001o*(-1.D0)+Bf001*(1.D0)+MD*DEBf00o*(-1.D0/2.D0)+MD*DEBf00*(1.D0/2.D0)&
        +MD*Bf00o*(1.D0)+MD*Bf00*(-1.D0)+x*DEC0022o*(-1.D0)+x*DEC002o*(-1.D0)+x*DEC0012o*(-2.D0)+x*DEC0011o*(-1.D0)+x*DEC001o*(-1.D0)&
        +x*DEC0000o*(-1.D0)+x*C0022o*(2.D0)+x*C002o*(2.D0)+x*C0012o*(4.D0)+x*C0011o*(2.D0)+x*C001o*(2.D0)+x*C0000o*(2.D0)&
        +x*DEBf001o*(1.D0)+x*DEBf001*(-3.D0/2.D0)+x*Bf001o*(-2.D0)+x*Bf001*(3.D0)+x*MD*DEC002o*(1.D0)+x*MD*DEC001o*(1.D0)&
        +x*MD*DEC00o*(1.D0)+x*MD*C002o*(-2.D0)+x*MD*C001o*(-2.D0)+x*MD*C00o*(-2.D0)+x*MD*DEBf00o*(-1.D0)+x*MD*DEBf00*(1.D0)&
        +x*MD*Bf00o*(2.D0)+x*MD*Bf00*(-2.D0)+x**2*DEBf001*(-1.D0)+x**2*Bf001*(2.D0))*(1.D0/2.D0)*(4.D0/3.D0)

   A(7) = A(7) &
	+(+DED00002*(-1.D0/27.D0)+D22222*(-1.D0/9.D0)+D2222*(-2.D0/9.D0)+D00222*(-19.D0/9.D0)+D222*(-1.D0/9.D0)+D0022*(-10.D0/3.D0)&
	+D002*(-11.D0/9.D0)+D00002*(2.D0/9.D0)+D00123*(-38.D0/9.D0)+D123*(-2.D0/9.D0)+D1223*(-4.D0/9.D0)+D12223*(-4.D0/9.D0)&
	+D12222*(-4.D0/9.D0)+D1222*(-8.D0/9.D0)+D122*(-4.D0/9.D0)+D11233*(-2.D0/9.D0)+D11223*(-4.D0/9.D0)+D1122*(-4.D0/9.D0)&
	+D00112*(38.D0/9.D0)+D112*(-2.D0/9.D0)+D11122*(4.D0/9.D0)+D11112*(2.D0/9.D0)+MD*D2222*(1.D0/9.D0)+MD*D222*(2.D0/9.D0)&
	+MD*D22*(1.D0/9.D0)+MD*D1223*(2.D0/9.D0)+MD*D1222*(2.D0/9.D0)+MD*D122*(4.D0/9.D0)+MD*D12*(2.D0/9.D0)&
	+MD*D1123*(2.D0/9.D0)+MD*D1122*(-2.D0/9.D0)+MD*D1112*(-2.D0/9.D0)+t*D22222*(1.D0/12.D0)+t*D2222*(1.D0/6.D0)&
	+t*D00222*(7.D0/9.D0)+t*D222*(1.D0/12.D0)+t*D0022*(11.D0/9.D0)+t*D002*(4.D0/9.D0)+t*D00123*(11.D0/9.D0)+t*D123*(5.D0/18.D0)&
	+t*D1223*(5.D0/9.D0)+t*D12223*(4.D0/9.D0)+t*D12222*(2.D0/9.D0)+t*D1222*(4.D0/9.D0)+t*D122*(2.D0/9.D0)+t*D11233*(1.D0/2.D0)&
	+t*D11223*(2.D0/9.D0)+t*D11222*(-1.D0/9.D0)+t*D1122*(1.D0/9.D0)+t*D00112*(-11.D0/9.D0)+t*D112*(1.D0/18.D0)+t*D11123*(-4.D0/9.D0)&
	+t*D11122*(-2.D0/9.D0)+t*D11112*(-1.D0/18.D0)+t*MD*D1222*(-1.D0/18.D0)+t*MD*D122*(-1.D0/9.D0)+t*MD*D12*(-1.D0/18.D0)&
	+t*MD*D1123*(-1.D0/18.D0)+t*MD*D1112*(1.D0/18.D0)+t**2*D123*(-1.D0/18.D0)+t**2*D1223*(-1.D0/9.D0)+t**2*D12223*(-1.D0/18.D0)&
	+t**2*D11233*(-1.D0/9.D0)+t**2*D11123*(1.D0/9.D0)+x*D22222*(-2.D0/9.D0)+x*D2222*(-4.D0/9.D0)+x*D222*(-2.D0/9.D0)&
	+x*D12223*(-4.D0/9.D0)+x*D12222*(-4.D0/9.D0)+x*D1222*(-8.D0/9.D0)+x*D122*(-4.D0/9.D0)+x*D11223*(-4.D0/9.D0)+x*D11222*(4.D0/9.D0)&
	+x*D11122*(4.D0/9.D0)+x*MD*D2222*(1.D0/9.D0)+x*MD*D222*(2.D0/9.D0)+x*MD*D22*(1.D0/9.D0)+x*MD*D1223*(2.D0/9.D0)&
	+x*MD*D1122*(-2.D0/9.D0)+x*t*D22222*(1.D0/6.D0)+x*t*D2222*(1.D0/3.D0)+x*t*D222*(1.D0/6.D0)+x*t*D12223*(1.D0/3.D0)&
	+x*t*D12222*(2.D0/9.D0)+x*t*D1222*(4.D0/9.D0)+x*t*D122*(2.D0/9.D0)+x*t*D11223*(2.D0/9.D0)+x*t*D11222*(-1.D0/3.D0)&
	+x*t*D11122*(-2.D0/9.D0)+x**2*D22222*(1.D0/9.D0)+x**2*D2222*(2.D0/9.D0)+x**2*D222*(1.D0/9.D0)+x**2*D12223*(2.D0/9.D0)&
	+x**2*D11222*(-2.D0/9.D0)&
!!$
        +DEC0022*(1.D0/6.D0)+DEC002*(4.D0/27.D0)+DEC0012*(7.D0/54.D0)+C2222*(-5.D0/18.D0)+C222*(-5.D0/18.D0)&
	+C0022*(-4.D0/3.D0)+C002*(-8.D0/9.D0)+C1222*(-11.D0/18.D0)+C122*(-1.D0/3.D0)+C0012*(-10.D0/9.D0)&
	+C1122*(-7.D0/18.D0)+C112*(-1.D0/18.D0)+C1112*(-1.D0/18.D0)+MD*C222*(-1.D0/18.D0)+MD*C22*(-1.D0/18.D0)&
	+MD*C12*(1.D0/18.D0)+MD*C112*(1.D0/18.D0)+x*C2222*(-4.D0/9.D0)+x*C222*(-4.D0/9.D0)+x*C1222*(-2.D0/3.D0)&
	+x*C122*(-2.D0/9.D0)+x*C1122*(-2.D0/9.D0)&
!!$
        +DEBf111*(-1.D0/216.D0)+Bf111*(1.D0/36.D0))*(2.D0*x)**2&
!!$
        +(0.D0)*(2.D0*x)*(-3.D0/2.D0)*4.D0/3.D0&
!!$
        +(0.D0)*(1.D0/2.D0)*(4.D0/3.D0)


   A(8) = A(8) &
        +(0.D0&
!!$
        +0.D0&
!!$
        +0.D0&
!!$
        +DEBf111*(-1.D0/54.D0)+DEBf001*(-1.D0/18.D0)+Bf111*(-1.D0/18.D0)+Bf001*(-1.D0/3.D0)+x*DEBf111*(-1.D0/27.D0)+x*Bf111*(-1.D0/9.D0))*(2.D0*x)**2&
!!$
        +(+DEBf001o*(1.D0/54.D0)+DEBf001*(-1.D0/54.D0)+Bf001o*(-1.D0/9.D0)+Bf001*(1.D0/9.D0)+DEC0022o*(-2.D0/27.D0)+DEC0012o*(-4.D0/27.D0)&
        +DEC0011o*(-2.D0/27.D0)+DEC0000o*(5.D0/27.D0)+DEC0022*(2.D0/27.D0)+DEC0012*(4.D0/27.D0)+DEC0011*(2.D0/27.D0)+DEC0000*(-5.D0/27.D0)&
        +C0022o*(-2.D0/9.D0)+C0012o*(-4.D0/9.D0)+C0011o*(-2.D0/9.D0)+C0000o*(-22.D0/9.D0)+C0022*(2.D0/9.D0)+C0012*(4.D0/9.D0)+C0011*(2.D0/9.D0)&
        +C0000*(22.D0/9.D0)+MD*DEBf00o*(-1.D0/54.D0)+MD*DEBf00*(1.D0/54.D0)+MD*Bf00o*(1.D0/9.D0)+MD*Bf00*(-1.D0/9.D0)+MD*DEC002o*(2.D0/27.D0)&
        +MD*DEC001o*(2.D0/27.D0)+MD*DEC002*(-2.D0/27.D0)+MD*DEC001*(-2.D0/27.D0)+MD*C002o*(2.D0/9.D0)+MD*C001o*(2.D0/9.D0)+MD*C002*(-2.D0/9.D0)&
        +MD*C001*(-2.D0/9.D0)+x*DEBf001*(-1.D0/27.D0)+x*Bf001*(2.D0/9.D0)+x*DEC0022*(5.D0/27.D0)+x*DEC0012*(2.D0/9.D0)+x*DEC0011*(1.D0/27.D0)&
        +x*DEC0000*(-8.D0/27.D0)+x*C0022*(5.D0/9.D0)+x*C0012*(2.D0/3.D0)+x*C0011*(1.D0/9.D0)+x*C0000*(16.D0/9.D0)+x*MD*DEC002*(-1.D0/9.D0)&
        +x*MD*DEC001*(-1.D0/27.D0)+x*MD*C002*(-1.D0/3.D0)+x*MD*C001*(-1.D0/9.D0)+x**2*DEC0022*(2.D0/27.D0)+x**2*DEC0012*(1.D0/27.D0)&
        +x**2*C0022*(2.D0/9.D0)+x**2*C0012*(4.D0/9.D0))*(2.D0*x)*(-3.D0/2.D0)*4.D0/3.D0&
!!$
        +(+DEBf001o*(-1.D0/2.D0)+DEBf001*(1.D0/2.D0)+Bf001o*(1.D0)+Bf001*(-1.D0)+MD*DEBf00o*(1.D0/2.D0)+MD*DEBf00*(-1.D0/2.D0)&
        +MD*Bf00o*(-1.D0)+MD*Bf00*(1.D0)+x*DEC0022o*(1.D0)+x*DEC002o*(1.D0)+x*DEC0012o*(2.D0)+x*DEC0011o*(1.D0)+x*DEC001o*(1.D0)&
        +x*DEC0000o*(1.D0)+x*C0022o*(-2.D0)+x*C002o*(-2.D0)+x*C0012o*(-4.D0)+x*C0011o*(-2.D0)+x*C001o*(-2.D0)+x*C0000o*(-2.D0)&
        +x*DEBf001o*(-1.D0)+x*DEBf001*(3.D0/2.D0)+x*Bf001o*(2.D0)+x*Bf001*(-3.D0)+x*MD*DEC002o*(-1.D0)+x*MD*DEC001o*(-1.D0)&
        +x*MD*DEC00o*(-1.D0)+x*MD*C002o*(2.D0)+x*MD*C001o*(2.D0)+x*MD*C00o*(2.D0)+x*MD*DEBf00o*(1.D0)+x*MD*DEBf00*(-1.D0)&
        +x*MD*Bf00o*(-2.D0)+x*MD*Bf00*(2.D0)+x**2*DEBf001*(1.D0)+x**2*Bf001*(-2.D0))*(1.D0/2.D0)*(4.D0/3.D0)

	
   A=A/(2.D0*x)**2

   !* Graph 6:
   !*  Two-photon vertex graph
   !* Isospin factor = 4/3

   p1s=nuclm**2
   p2s=t
   p1p2s=nuclm**2

   m1s=deltam**2
   m2s=mpi2
   m3s=mpi2

   CALL Setdelta(Delta)

   C00=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C11=C0i(cc11,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C12=C0i(cc12,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C22=C0i(cc22,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C111=C0i(cc111,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C112=C0i(cc112,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C122=C0i(cc122,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C222=C0i(cc222,p1s,p2s,p1p2s,m1s,m2s,m3s)

   if(.not.Delta.eq.0.D0) then
     CALL Setdelta(0.D0)
     sw1=Delta
   else
     CALL Setdelta(DeltaAux)
     sw1=-DeltaAux
   endif

   C00D=C0i(cc00,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C001D=C0i(cc001,p1s,p2s,p1p2s,m1s,m2s,m3s)
   C002D=C0i(cc002,p1s,p2s,p1p2s,m1s,m2s,m3s)

   DEC00=2.D0*(C00-C00D)/sw1
   DEC001=2.D0*(C001-C001D)/sw1
   DEC002=2.D0*(C002-C002D)/sw1

   A(1) = A(1) &
	+(DEC002*(-2.D0/3.D0)+DEC001*(-2.D0/3.D0)+C002*(4.D0/3.D0)+C001*(4.D0/3.D0)+MD*DEC00*(2.D0/3.D0)&
	+MD*C00*(-4.D0/3.D0)+t*DEC002*(10.D0/27.D0)+t*DEC001*(10.D0/27.D0)+t*C002*(-8.D0/9.D0)+t*C122*(-4.D0/9.D0)&
	+t*C112*(-4.D0/9.D0)+t*C001*(-8.D0/9.D0)+t*MD*DEC00*(-10.D0/27.D0)+t*MD*C12*(4.D0/9.D0)+t*MD*C00*(8.D0/9.D0)&
	+t**2*C122*(1.D0/9.D0)+t**2*C112*(1.D0/9.D0)+t**2*MD*C12*(-1.D0/9.D0))/2.D0

   A=gafac*A


END SUBROUTINE OldCorrIsoCovOneDeltaLoop

SUBROUTINE s_nucl_sub(s,t,Q1,amm1,A)
  implicit none
  REAL*8 s,x,t,kappa,Q,Q1,amm,amm1
  COMPLEX*16 A(8)

  Q=1.D0
  amm=1.79D0

  x=(s-1.d0)/2.d0

  A=0.d0
 
  A(1) = A(1) &
       +x**(-1)*t*amm*amm1*(-1.D0/4.D0)+x**(-1)*t*amm*Q1*(-1.D0/4.D0)+x**(-1)*t*Q*amm1*(-1.D0/4.D0)&
       +x**(-1)*t*Q*Q1*(-1.D0/4.D0)+Q*Q1*(-1.D0)+t*amm*amm1*(-1.D0/8.D0)+x*amm*amm1*(-1.D0)

  A(2) = A(2) &
       +x**(-1)*amm*amm1*(1.D0/2.D0)+x**(-1)*amm*Q1*(1.D0/2.D0)+x**(-1)*Q*amm1*(1.D0/2.D0)&
       +x**(-1)*t*amm*amm1*(-1.D0/16.D0)+amm*amm1*(-1.D0/4.D0)

  A(3) = A(3) &
       +x**(-1)*t*amm*amm1*(-1.D0/4.D0)+x**(-1)*t*amm*Q1*(-1.D0/4.D0)+x**(-1)*t*Q*amm1*(-1.D0/4.D0)&
       +x**(-1)*t*Q*Q1*(-1.D0/4.D0)+Q*Q1*(-1.D0)+t*amm*amm1*(-1.D0/8.D0)+x*amm*amm1*(-1.D0)

  A(4) = A(4) &
       +x**(-1)*amm*amm1*(-1.D0/2.D0)+x**(-1)*amm*Q1*(-1.D0/2.D0)+x**(-1)*Q*amm1*(-1.D0/2.D0)&
       +x**(-1)*Q*Q1*(-1.D0/2.D0)+amm*amm1*(-1.D0/4.D0)

  A(5) = A(5) &
       +x**(-1)*amm*amm1*(1.D0/2.D0)+x**(-1)*amm*Q1*(1.D0/2.D0)+x**(-1)*Q*amm1*(1.D0/2.D0)+x**(-1)*Q*Q1*(1.D0/2.D0)

  A(6) = A(6) &
       +x**(-1)*amm*Q1*(-1.D0/4.D0)+x**(-1)*Q*amm1*(-1.D0/4.D0)+x**(-1)*Q*Q1*(-1.D0/2.D0)

  A(7) = A(7) &
       +x**(-1)*amm*amm1*(1.D0/8.D0)

  A(8) = A(8) &
       +x**(-1)*amm*amm1*(1.D0/2.D0)+x**(-1)*amm*Q1*(1.D0/2.D0)+x**(-1)*Q*amm1*(1.D0/2.D0)&
       +x**(-1)*Q*Q1*(1.D0/2.D0)+amm*amm1*(1.D0/4.D0)

end subroutine s_nucl_sub

 
 
