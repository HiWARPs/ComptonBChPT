      SUBROUTINE s_nucl(x,t,Q,kappa,A)
        implicit none
        REAL*8 x,t,kappa,Q
        COMPLEX*16 A(8)
        REAL*8 kapa2,Qk2,denom
        REAL*8 a1,a2,a3,a4,a5,a6,a7,a8

      Qk2=(Q+kappa)**2
      kapa2=kappa**2

      a1 = x*t*kapa2 * ( 1.D0/8.D0 )
      a1 = a1 + x*Q**2
      a1 = a1 + x**2*kapa2 * ( 1.D0/2.D0 )
      a1 = a1 + t*Qk2 * ( 1.D0/4.D0 )

      a2 = x*kapa2 * (  1.D0/4.D0 ) + t*kapa2 * ( 1.D0/16.D0 )
      a2 = a2 + Q**2 * ( 1.D0/2.D0 )
      a2 = a2 + Qk2 * (  - 1.D0/2.D0 )

      a3 = a1

      a4 = x*kapa2 * ( 1.D0/4.D0 )
      a4 = a4 + Qk2 * ( 1.D0/2.D0 )

      a5 = Qk2 * (  - 1.D0/2.D0 )

      a6 = Q**2 * ( 1.D0/4.D0 )
      a6 = a6 + Qk2 * ( 1.D0/4.D0 )
      a6 = a6 + kapa2 * (  - 1.D0/4.D0 )

      a7 = kapa2 * (  - 1.D0/8.D0 )

      a8 = -a4

      denom=-1./(2.*x)
    
      A(1) = denom*a1      
      A(2) = denom*a2
      A(3) = denom*a3      
      A(4) = denom*a4      
      A(5) = denom*a5
      A(6) = denom*a6      
      A(7) = denom*a7      
      A(8) = denom*a8

    end subroutine s_nucl


      SUBROUTINE s_Delta(x,t,Md,gg,A)
        implicit none
        REAL*8, INTENT(in) :: x,t,Md,gg
        COMPLEX*16, INTENT(out) :: A(8)
        COMPLEX*16 denom, CI
        REAL*8 a1,a2,a3,a4,a5,a6,a7,a8,W,ggsq,edem

      a1 = x*t * ( 7.D0/6.D0 + 5.D0/6.D0*Md )*0.d0
      a1 = a1 + x*t**2 * ( 1.D0/4.D0 )*0.d0
      a1 = a1 + x**2*t * ( 7.D0/6.D0 )*0.d0
      a1 = a1 + x**2 * ( 2.D0/3.D0 + 2.D0/3.D0*Md )
      a1 = a1 + x**3 * ( 2.D0/3.D0 )*0.d0
      a1 = a1 + t * ( 1.D0/3.D0 + 1.D0/3.D0*Md )*0.d0
      a1 = a1 + t**2 * ( 1.D0/8.D0 )*0.d0

      a2 = x*t * (  - 1.D0/6.D0 )*0.d0
      a2 = a2 + x * (  - 5.D0/3.D0 - Md )*0.d0
      a2 = a2 + x**2 * (  - 7.D0/6.D0 )*0.d0
      a2 = a2 + t * (  - 1.D0/6.D0 + 1.D0/12.D0*Md )*0.d0
      a2 = a2 + t**2 * ( 1.D0/32.D0 )*0.d0
      a2 = a2 - 2.D0/3.D0 - 2.D0/3.D0*Md

      a3 =x*t * (  - 1.D0/2.D0 - 1.D0/3.D0*Md )
      a3 = a3 + x*t**2 * (  - 1.D0/12.D0 )
      a3 = a3 + x**2*t * (  - 5.D0/12.D0 )
      a3 = a3 + x**2 * (  - 1.D0/3.D0 - 1.D0/3.D0*Md )
      a3 = a3 + x**3 * (  - 1.D0/3.D0 )
      a3 = a3 + t * (  - 1.D0/6.D0 - 1.D0/6.D0*Md )
      a3 = a3 + t**2 * (  - 1.D0/24.D0 )

      a4 = x*t * ( 1.D0/2.D0 )
      a4 = a4 + x * (  - 2.D0/3.D0 - 1.D0/3.D0*Md )
      a4 = a4 + x**2 * ( 1.D0/3.D0 )
      a4 = a4 + t * ( 1.D0/4.D0 )
      a4 = a4 - 1.D0/3.D0 - 1.D0/3.D0*Md

      a5 =x*t * ( 1.D0/12.D0 )
      a5 = a5 + x * ( 2.D0/3.D0 + 1.D0/3.D0*Md )
      a5 = a5 + x**2 * ( 1.D0/3.D0 )
      a5 = a5 + t * ( 1.D0/24.D0 - 1.D0/24.D0*Md )
      a5 = a5 + 1.D0/3.D0 + 1.D0/3.D0*Md

      a6 = x*t * (  - 1.D0/12.D0 )
      a6 = a6 + x * (  - 1.D0/6.D0 - 1.D0/6.D0*Md )
      a6 = a6 + x**2 * (  - 1.D0/3.D0 )
      a6 = a6 + t * (  - 1.D0/24.D0 - 1.D0/24.D0*Md )

      a7 = x * ( 7.D0/12.D0 )
      a7 = a7 + t * (  - 1.D0/16.D0 )
      a7 = a7 + 7.D0/12.D0 - 1.D0/12.D0*Md

      a8 = x*t * ( 1.D0/6.D0 )
      a8 = a8 + x**2 * ( 1.D0/6.D0 )
      a8 = a8 + t * ( 1.D0/12.D0 )

      W=1.d0+2.d0*x
      if(W.gt.1.d0) then
         CI=(0.d0,1.d0)
         W=dsqrt(W)
         !denom=1./(2.d0*x-(Md-0.5d0*CI*edem(W,Md))**2+1.d0)
         denom=1./(2.d0*x-Md**2+CI*Md*edem(W,Md)+1.d0)
      else
         denom=1./(2.d0*x-Md**2+1.d0)
      endif

      ggsq=-gg*gg/4.d0
      ggsq=(2.d0/3.d0)*ggsq ! isospin factor T_3*T_3
      ! Jones & Scadron convention:
      ggsq=(3.d0/(1.d0+Md))**2*ggsq
      A(1)=ggsq*denom*a1      
      A(2)=ggsq*denom*a2
      A(3)=ggsq*denom*a3      
      A(4)=ggsq*denom*a4      
      A(5)=ggsq*denom*a5
      A(6)=ggsq*denom*a6      
      A(7)=ggsq*denom*a7      
      A(8)=ggsq*denom*a8

    end subroutine s_Delta

    SUBROUTINE s_Delta2(x,t,Md,mpi,G1,G2,A)
      implicit none
      REAL*8, INTENT(in) :: x,t,Md,mpi,G1,G2
      COMPLEX*16, INTENT(out) :: A(8)
      COMPLEX*16 denom, CI
      REAL*8 W,s,gg,GM,ggsq,edem,ImGE
      Real*8 ImSigmaPrime, Width
      COMPLEX*16 a1,a2,a3,a4,a5,a6,a7,a8,GE
      

      GM=(3.d0/2.d0)/(1.d0+Md)*G1
      GE=(3.d0/2.d0)/(1.d0+Md)*G2
      !GM=1.d0/(2.d0*Md)*G1
      !GE=1.d0/(2.d0*Md)*G2

      s=1.d0+2.d0*x
      if(s.gt.(1.d0+mpi)**2) then
         CI=(0.d0,1.d0)
        W=dsqrt(s)
         !denom=1./(2.d0*x-(Md-0.5d0*CI*edem(W,Md))**2+1.d0)
          denom=1./(s-Md**2+CI*Md*edem(W,Md))  !; denom=1./(s-Md**2)
!!        Width=0.115
!!	ImSigmaPrime= -0.655
!!        denom = 1.d0/( (s-Md**2)*(1.d0 - CI*ImSigmaPrime) + CI*Md*Width)
!!$         !GE= GE + CI*ImGE(W,Md)
         write(*,*) "width =", edem(Md,Md)*0.9383
!!$         write(*,*) "GE = ", GE/(3.d0/2.d0)*(1.d0+Md)
      else
         denom=1./(s-Md**2)
      endif

      
      a1 = x*t*GM*GE * (  - 1.D0/6.D0 - 1.D0/6.D0*Md )
      a1 = a1 + x*t*GM**2 * ( 7.D0/6.D0 + 5.D0/6.D0*Md )
      a1 = a1 + x*t**2*GM**2 * ( 1.D0/4.D0 )
      a1 = a1 + x**2*t*GM*GE * (  - 1.D0/6.D0 )
      a1 = a1 + x**2*t*GM**2 * ( 7.D0/6.D0 )
      a1 = a1 + x**2*t*GE**2 * ( 1.D0/6.D0 )
      a1 = a1 + x**2*GM**2 * ( 2.D0/3.D0 + 2.D0/3.D0*Md )
      a1 = a1 + x**2*GE**2 * ( 2.D0/3.D0 - 2.D0/3.D0*Md )
      a1 = a1 + x**3*GM*GE * (  - 2.D0/3.D0 )
      a1 = a1 + x**3*GM**2 * ( 2.D0/3.D0 )
      a1 = a1 + x**3*GE**2 * ( 2.D0/3.D0 )
      a1 = a1 + t*GM**2 * ( 1.D0/3.D0 + 1.D0/3.D0*Md )
      a1 = a1 + t**2*GM**2 * ( 1.D0/8.D0 )


      a2 = x*t*GM*GE * ( 5.D0/12.D0 )
      a2 = a2 + x*t*GM**2 * (  - 1.D0/6.D0 )
      a2 = a2 + x*t*GE**2 * ( 1.D0/3.D0 )
      a2 = a2 + x*GM*GE * ( 1.D0/3.D0 + 1.D0/3.D0*Md )
      a2 = a2 + x*GM**2 * (  - 5.D0/3.D0 - Md )
      a2 = a2 + x*GE**2 * ( 2.D0/3.D0 - 2.D0/3.D0*Md )
      a2 = a2 + x**2*GM*GE * ( 2.D0/3.D0 )
      a2 = a2 + x**2*GM**2 * (  - 7.D0/6.D0 )
      a2 = a2 + x**2*GE**2 * ( 5.D0/6.D0 )
      a2 = a2 + t*GM**2 * (  - 1.D0/6.D0 + 1.D0/12.D0*Md )
      a2 = a2 + t*GE**2 * ( 1.D0/12.D0 - 1.D0/12.D0*Md )
      a2 = a2 + t**2*GM*GE * ( 1.D0/16.D0 )
      a2 = a2 + t**2*GM**2 * ( 1.D0/32.D0 )
      a2 = a2 + t**2*GE**2 * ( 1.D0/32.D0 )
      a2 = a2 + GM**2 * (  - 2.D0/3.D0 - 2.D0/3.D0*Md )


      a3 = x*t*GM*GE * ( 1.D0/3.D0 + 1.D0/3.D0*Md )
      a3 = a3 + x*t*GM**2 * ( 1.D0/6.D0 + 1.D0/3.D0*Md )
      a3 = a3 + x*t**2*GM**2 * ( 1.D0/4.D0 )
      a3 = a3 + x**2*t*GM*GE * ( 1.D0/3.D0 )
      a3 = a3 + x**2*t*GM**2 * ( 11.D0/12.D0 )
      a3 = a3 + x**2*t*GE**2 * (  - 1.D0/12.D0 )
      a3 = a3 + x**2*GM**2 * (  - 1.D0/3.D0 - 1.D0/3.D0*Md )
      a3 = a3 + x**2*GE**2 * (  - 1.D0/3.D0 + 1.D0/3.D0*Md )
      a3 = a3 + x**3*GM*GE * ( 4.D0/3.D0 )
      a3 = a3 + x**3*GM**2 * (  - 1.D0/3.D0 )
      a3 = a3 + x**3*GE**2 * (  - 1.D0/3.D0 )
      a3 = a3 + t*GM**2 * (  - 1.D0/6.D0 - 1.D0/6.D0*Md )
      a3 = a3 + t**2*GM**2 * ( 1.D0/8.D0 )


      a4 = x*t*GM**2 * ( 1.D0/2.D0 )
      a4 = a4 + x*GM*GE * (  - 1.D0/3.D0 - 1.D0/3.D0*Md )
      a4 = a4 + x*GM**2 * (  - 2.D0/3.D0 - 1.D0/3.D0*Md )
      a4 = a4 + x**2*GM*GE * (  - 1.D0/3.D0 )
      a4 = a4 + x**2*GM**2 * ( 1.D0/3.D0 )
      a4 = a4 + x**2*GE**2 * ( 1.D0/3.D0 )
      a4 = a4 + t*GM**2 * ( 1.D0/4.D0 )
      a4 = a4 + GM**2 * (  - 1.D0/3.D0 - 1.D0/3.D0*Md )


      a5 = x*t*GM*GE * ( 1.D0/4.D0 )
      a5 = a5 + x*t*GM**2 * (  - 1.D0/4.D0 )
      a5 = a5 + x*GM*GE * (  - 1.D0/6.D0 - 1.D0/6.D0*Md )
      a5 = a5 + x*GM**2 * (  - 1.D0/3.D0*Md )
      a5 = a5 + x*GE**2 * (  - 1.D0/6.D0 + 1.D0/6.D0*Md )
      a5 = a5 + x**2*GM*GE * ( 1.D0 )
      a5 = a5 + x**2*GM**2 * (  - 1.D0 )
      a5 = a5 + t*GM*GE * ( 1.D0/8.D0 + 1.D0/8.D0*Md )
      a5 = a5 + t*GM**2 * (  - 1.D0/8.D0 + 1.D0/8.D0*Md )
      a5 = a5 + GM**2 * ( 1.D0/3.D0 + 1.D0/3.D0*Md )

      a6 = x*t*GM*GE * (  - 1.D0/4.D0 )
      a6 = a6 + x*t*GM**2 * ( 1.D0/4.D0 )
      a6 = a6 + x*GM*GE * (  - 1.D0/2.D0 + 1.D0/2.D0*Md )
      a6 = a6 + x*GM**2 * ( 1.D0/2.D0 + 1.D0/2.D0*Md )
      a6 = a6 + x**2*GM*GE * (  - 1.D0 )
      a6 = a6 + x**2*GM**2 * ( 1.D0 )
      a6 = a6 + t*GM*GE * (  - 1.D0/8.D0 + 1.D0/8.D0*Md )
      a6 = a6 + t*GM**2 * ( 1.D0/8.D0 + 1.D0/8.D0*Md )

      a7 = x*GM*GE * (  - 1.D0/3.D0 )
      a7 = a7 + x*GM**2 * ( 7.D0/12.D0 )
      a7 = a7 + x*GE**2 * (  - 5.D0/12.D0 )
      a7 = a7 + t*GM*GE * (  - 1.D0/8.D0 )
      a7 = a7 + t*GM**2 * (  - 1.D0/16.D0 )
      a7 = a7 + t*GE**2 * (  - 1.D0/16.D0 )
      a7 = a7 + GM*GE * ( 1.D0/2.D0 + 1.D0/2.D0*Md )
      a7 = a7 + GM**2 * ( 7.D0/12.D0 + 7.D0/12.D0*Md )
      a7 = a7 + GE**2 * ( 1.D0/12.D0 - 1.D0/12.D0*Md )

      a8 = x*t*GM**2 * (  - 1.D0/2.D0 )
      a8 = a8 + x*GM*GE * (  - 2.D0/3.D0 - 2.D0/3.D0*Md )
      a8 = a8 + x*GM**2 * ( 8.D0/3.D0 + 4.D0/3.D0*Md )
      a8 = a8 + x**2*GM*GE * (  - 2.D0/3.D0 )
      a8 = a8 + x**2*GM**2 * ( 1.D0/6.D0 )
      a8 = a8 + x**2*GE**2 * ( 1.D0/6.D0 )
      a8 = a8 + t*GM**2 * (  - 1.D0/4.D0 )
      a8 = a8 + GM**2 * ( 4.D0/3.D0 + 4.D0/3.D0*Md )


      
      gg=1.d0
      ggsq=-gg*gg
      ggsq=(2.d0/3.d0)*ggsq ! isospin factor T_3*T_3
   
      A(1)=ggsq*denom*a1      
      A(2)=ggsq*denom*a2
      A(3)=ggsq*denom*a3      
      A(4)=ggsq*denom*a4      
      A(5)=ggsq*denom*a5
      A(6)=ggsq*denom*a6      
      A(7)=ggsq*denom*a7      
      A(8)=ggsq*denom*a8

    end subroutine s_Delta2


    REAL*8 FUNCTION EDEM(W,Md)
      implicit none
      REAL*8, intent(in) :: W,Md
      REAL*8 gA,hA,fpi,pi0mass,pionm
      COMMON/Consts/hA,gA,fpi,pi0mass,pionm
      REAL*8 what2,Ehat,phat,K0,Csq
      REAL*8 PI,MPI,MPI2

      MPI=pionm
      MPI2 = MPI*MPI
      IF(W.LT.(1.+MPI)) THEN
         EDEM=0.D0
         RETURN
      ELSE
         PI=4.d0*datan(1.d0)
         Csq=(hA*Md/(8*PI*fpi))**2
         what2 = ( (W*W - 1.d0+MPI2)/(2.D0*W) )**2
         Ehat=(W*W +1.d0-MPI2)/(2.D0*W)
         phat=dsqrt(what2-MPI2)         
         K0 = (Ehat+1.d0)*phat**3/W
         EDEM =4.d0*PI/3.d0*Csq*K0/Md**2
         EDEM=(W/Md)**2*EDEM
         RETURN
      endif
    END FUNCTION EDEM

    REAL*8 FUNCTION ImGE(W,Md)
      implicit none
      REAL*8, intent(in) :: W,Md
      REAL*8 gA,hA,fpi,pi0mass,pionm
      COMMON/Consts/hA,gA,fpi,pi0mass,pionm
      REAL*8 what,Ehat,phat,qhat,K0,s,gg
      REAL*8 PI,MPI,MPI2,rQfun

      MPI=pionm
      MPI2 = MPI*MPI
      IF(W.LT.(1.d0+MPI)) THEN
         ImGE=0.D0
         RETURN
      ELSE
         PI=4.d0*datan(1.d0)
         gg=4.d0*gA*hA/(2.*fpi)**2/(16.d0*PI)
         MPI2 = MPI*MPI
         s=W*W
         what =  (s - 1.d0+MPI2)/(2.D0*W) 
         Ehat=(s +1.d0-MPI2)/(2.D0*W)
         phat=dsqrt(what*what-MPI2)  
         qhat = (s-1.d0)/(2.*W)
         K0 = Ehat*phat**2/W  
         ImGE = gg*K0/(qhat*W) *(1./3.)/what*(rQfun(0,what/phat) - rQfun(2,what/phat) ) 
         ImGE=(W/Md)*ImGE
         !write(*,*) qhat*Mn,phat*Mn
         RETURN
      endif
    END FUNCTION ImGE

RECURSIVE FUNCTION rQfun(l,x)  RESULT(Qfun) ! Legendre function of the 2d kind
  REAL*8, INTENT(in):: x                  ! with real argument: x
  INTEGER, INTENT(in):: l                  ! and  integer index: l
  REAL*8 Qfun
  REAL*8, PARAMETER :: EPS=1.d-15
  REAL*8 rlg
  
  IF(dabs(x).LT.(1.d0+EPS).AND.dabs(x).GT.(1.d0-EPS)) THEN
     write(*,*) 'Qfun: Hitting the log singularity, return Qfun=0.'
     Qfun = 0.
     RETURN
  ENDIF
  rlg= dlog( dabs((1.d0+x)/(1.d0-x)) ) 
  SELECT CASE(l)
     CASE(0); Qfun = 0.5d0*rlg
     CASE(1); Qfun = 0.5d0*x*rlg - 1.d0
     CASE(2); Qfun = 0.25d0*(3.d0*x*x - 1.d0)*rlg - 1.5d0*x
     CASE(3); Qfun = 0.25d0*(5.d0*x*x - 3.d0)*x*rlg-2.5*x*x+2.d0/3.d0  ! l=3
     CASE(4:);Qfun = ( (2*l-1)*x*rQfun(l-1,x) - (l-1)*rQfun(l-2,x) )/(1.d0*l)
     CASE DEFAULT; write(*,*) 'Qfun: l is negative, out of range'
   END SELECT 
  RETURN
END FUNCTION rQfun
