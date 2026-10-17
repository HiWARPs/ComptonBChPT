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
      REAL*8 W,s,gg,ggsq,edem
      Real*8 ImSigmaPrime, Width
      COMPLEX*16 a1,a2,a3,a4,a5,a6,a7,a8,GE,GM
      !* [v2] Delta-region options (set by the main program from namelist &v2):
      !*   ldvertex : imaginary parts of the one-loop gamma N Delta vertex (LMP2015 Fig. 4,
      !*              PV2006 App. B), s-channel only (they vanish below the pi N threshold)
      !*   ldrunning: running Re gM(s), Re gE(s), M_Delta(s), 1+Re tau_R(s) (LMP2015 App. D); the
      !*              subtraction scheme (irun = 1, 2, 3) is passed to delta_running via COMMON /v2drun/
      LOGICAL ldvertex, ldrunning
      COMMON /v2delta/ ldvertex, ldrunning
      REAL*8 vImgM, vImgE, ReG1, ReG2, Mrun, ztau, Mp
      CI=(0.d0,1.d0)
      s=1.d0+2.d0*x
      ReG1=G1; ReG2=G2; Mrun=Md; ztau=1.d0
      IF(ldrunning) CALL delta_running(s,Md,G1,G2,ReG1,ReG2,Mrun,ztau)
      GM=(3.d0/2.d0)/(1.d0+Md)*ReG1
      GE=(3.d0/2.d0)/(1.d0+Md)*ReG2
      !GM=1.d0/(2.d0*Md)*G1
      !GE=1.d0/(2.d0*Md)*G2

      if(s.gt.(1.d0+mpi)**2) then
        W=dsqrt(s)
         !denom=1./(2.d0*x-(Md-0.5d0*CI*edem(W,Md))**2+1.d0)
          denom=1./(s-Mrun**2+CI*Md*edem(W,Md))  !; denom=1./(s-Md**2)
!!        Width=0.115
!!	ImSigmaPrime= -0.655
!!        denom = 1.d0/( (s-Md**2)*(1.d0 - CI*ImSigmaPrime) + CI*Md*Width)
!!$         !GE= GE + CI*ImGE(W,Md)
         !* [v2] per-call debug print removed: write(*,*) "width =", edem(Md,Md)*0.9383
         IF(ldvertex) THEN
            CALL gnd_vertex_im(W,Md,vImgM,vImgE)
            GM=GM+CI*(3.d0/2.d0)/(1.d0+Md)*vImgM
            GE=GE+CI*(3.d0/2.d0)/(1.d0+Md)*vImgE
         ENDIF
!!$         write(*,*) "GE = ", GE/(3.d0/2.d0)*(1.d0+Md)
      else
         denom=1./(s-Mrun**2)
      endif
      denom=denom*ztau
      Mp=Md
      IF(ldrunning) Mp=Mrun      ! numerator (pslash + M_Delta(s)) with the running mass

      
      a1 = x*t*GM*GE * (  - 1.D0/6.D0 - 1.D0/6.D0*Mp )
      a1 = a1 + x*t*GM**2 * ( 7.D0/6.D0 + 5.D0/6.D0*Mp )
      a1 = a1 + x*t**2*GM**2 * ( 1.D0/4.D0 )
      a1 = a1 + x**2*t*GM*GE * (  - 1.D0/6.D0 )
      a1 = a1 + x**2*t*GM**2 * ( 7.D0/6.D0 )
      a1 = a1 + x**2*t*GE**2 * ( 1.D0/6.D0 )
      a1 = a1 + x**2*GM**2 * ( 2.D0/3.D0 + 2.D0/3.D0*Mp )
      a1 = a1 + x**2*GE**2 * ( 2.D0/3.D0 - 2.D0/3.D0*Mp )
      a1 = a1 + x**3*GM*GE * (  - 2.D0/3.D0 )
      a1 = a1 + x**3*GM**2 * ( 2.D0/3.D0 )
      a1 = a1 + x**3*GE**2 * ( 2.D0/3.D0 )
      a1 = a1 + t*GM**2 * ( 1.D0/3.D0 + 1.D0/3.D0*Mp )
      a1 = a1 + t**2*GM**2 * ( 1.D0/8.D0 )


      a2 = x*t*GM*GE * ( 5.D0/12.D0 )
      a2 = a2 + x*t*GM**2 * (  - 1.D0/6.D0 )
      a2 = a2 + x*t*GE**2 * ( 1.D0/3.D0 )
      a2 = a2 + x*GM*GE * ( 1.D0/3.D0 + 1.D0/3.D0*Mp )
      a2 = a2 + x*GM**2 * (  - 5.D0/3.D0 - Mp )
      a2 = a2 + x*GE**2 * ( 2.D0/3.D0 - 2.D0/3.D0*Mp )
      a2 = a2 + x**2*GM*GE * ( 2.D0/3.D0 )
      a2 = a2 + x**2*GM**2 * (  - 7.D0/6.D0 )
      a2 = a2 + x**2*GE**2 * ( 5.D0/6.D0 )
      a2 = a2 + t*GM**2 * (  - 1.D0/6.D0 + 1.D0/12.D0*Mp )
      a2 = a2 + t*GE**2 * ( 1.D0/12.D0 - 1.D0/12.D0*Mp )
      a2 = a2 + t**2*GM*GE * ( 1.D0/16.D0 )
      a2 = a2 + t**2*GM**2 * ( 1.D0/32.D0 )
      a2 = a2 + t**2*GE**2 * ( 1.D0/32.D0 )
      a2 = a2 + GM**2 * (  - 2.D0/3.D0 - 2.D0/3.D0*Mp )


      a3 = x*t*GM*GE * ( 1.D0/3.D0 + 1.D0/3.D0*Mp )
      a3 = a3 + x*t*GM**2 * ( 1.D0/6.D0 + 1.D0/3.D0*Mp )
      a3 = a3 + x*t**2*GM**2 * ( 1.D0/4.D0 )
      a3 = a3 + x**2*t*GM*GE * ( 1.D0/3.D0 )
      a3 = a3 + x**2*t*GM**2 * ( 11.D0/12.D0 )
      a3 = a3 + x**2*t*GE**2 * (  - 1.D0/12.D0 )
      a3 = a3 + x**2*GM**2 * (  - 1.D0/3.D0 - 1.D0/3.D0*Mp )
      a3 = a3 + x**2*GE**2 * (  - 1.D0/3.D0 + 1.D0/3.D0*Mp )
      a3 = a3 + x**3*GM*GE * ( 4.D0/3.D0 )
      a3 = a3 + x**3*GM**2 * (  - 1.D0/3.D0 )
      a3 = a3 + x**3*GE**2 * (  - 1.D0/3.D0 )
      a3 = a3 + t*GM**2 * (  - 1.D0/6.D0 - 1.D0/6.D0*Mp )
      a3 = a3 + t**2*GM**2 * ( 1.D0/8.D0 )


      a4 = x*t*GM**2 * ( 1.D0/2.D0 )
      a4 = a4 + x*GM*GE * (  - 1.D0/3.D0 - 1.D0/3.D0*Mp )
      a4 = a4 + x*GM**2 * (  - 2.D0/3.D0 - 1.D0/3.D0*Mp )
      a4 = a4 + x**2*GM*GE * (  - 1.D0/3.D0 )
      a4 = a4 + x**2*GM**2 * ( 1.D0/3.D0 )
      a4 = a4 + x**2*GE**2 * ( 1.D0/3.D0 )
      a4 = a4 + t*GM**2 * ( 1.D0/4.D0 )
      a4 = a4 + GM**2 * (  - 1.D0/3.D0 - 1.D0/3.D0*Mp )


      a5 = x*t*GM*GE * ( 1.D0/4.D0 )
      a5 = a5 + x*t*GM**2 * (  - 1.D0/4.D0 )
      a5 = a5 + x*GM*GE * (  - 1.D0/6.D0 - 1.D0/6.D0*Mp )
      a5 = a5 + x*GM**2 * (  - 1.D0/3.D0*Mp )
      a5 = a5 + x*GE**2 * (  - 1.D0/6.D0 + 1.D0/6.D0*Mp )
      a5 = a5 + x**2*GM*GE * ( 1.D0 )
      a5 = a5 + x**2*GM**2 * (  - 1.D0 )
      a5 = a5 + t*GM*GE * ( 1.D0/8.D0 + 1.D0/8.D0*Mp )
      a5 = a5 + t*GM**2 * (  - 1.D0/8.D0 + 1.D0/8.D0*Mp )
      a5 = a5 + GM**2 * ( 1.D0/3.D0 + 1.D0/3.D0*Mp )

      a6 = x*t*GM*GE * (  - 1.D0/4.D0 )
      a6 = a6 + x*t*GM**2 * ( 1.D0/4.D0 )
      a6 = a6 + x*GM*GE * (  - 1.D0/2.D0 + 1.D0/2.D0*Mp )
      a6 = a6 + x*GM**2 * ( 1.D0/2.D0 + 1.D0/2.D0*Mp )
      a6 = a6 + x**2*GM*GE * (  - 1.D0 )
      a6 = a6 + x**2*GM**2 * ( 1.D0 )
      a6 = a6 + t*GM*GE * (  - 1.D0/8.D0 + 1.D0/8.D0*Mp )
      a6 = a6 + t*GM**2 * ( 1.D0/8.D0 + 1.D0/8.D0*Mp )

      a7 = x*GM*GE * (  - 1.D0/3.D0 )
      a7 = a7 + x*GM**2 * ( 7.D0/12.D0 )
      a7 = a7 + x*GE**2 * (  - 5.D0/12.D0 )
      a7 = a7 + t*GM*GE * (  - 1.D0/8.D0 )
      a7 = a7 + t*GM**2 * (  - 1.D0/16.D0 )
      a7 = a7 + t*GE**2 * (  - 1.D0/16.D0 )
      a7 = a7 + GM*GE * ( 1.D0/2.D0 + 1.D0/2.D0*Mp )
      a7 = a7 + GM**2 * ( 7.D0/12.D0 + 7.D0/12.D0*Mp )
      a7 = a7 + GE**2 * ( 1.D0/12.D0 - 1.D0/12.D0*Mp )

      a8 = x*t*GM**2 * (  - 1.D0/2.D0 )
      a8 = a8 + x*GM*GE * (  - 2.D0/3.D0 - 2.D0/3.D0*Mp )
      a8 = a8 + x*GM**2 * ( 8.D0/3.D0 + 4.D0/3.D0*Mp )
      a8 = a8 + x**2*GM*GE * (  - 2.D0/3.D0 )
      a8 = a8 + x**2*GM**2 * ( 1.D0/6.D0 )
      a8 = a8 + x**2*GE**2 * ( 1.D0/6.D0 )
      a8 = a8 + t*GM**2 * (  - 1.D0/4.D0 )
      a8 = a8 + GM**2 * ( 4.D0/3.D0 + 4.D0/3.D0*Mp )


      
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

!=============================================================================
! [v2, 2026-10] Imaginary parts of the one-loop (pi N) gamma N Delta vertex for real photons,
! as functions of the Delta invariant mass W (LMP2015 Sect. 2.2, Fig. 4: only the imaginary
! parts are of NLO in the Delta region). Literal transcription of Pascalutsa & Vanderhaeghen,
! PRD 73 (2006) 034003, App. B [Eqs. for Im g_i^+ of graphs (a) photon-on-pion and (b)
! photon-on-nucleon], evaluated in quadruple precision at Q^2/s = 1e-14 (the coefficient
! functions contain 1/Q^2 terms that cancel; residual error ~1e-13), and mapped onto the
! Lagrangian couplings with the on-shell relations (PV2006 Eq. gmecgi, Q^2=0):
!     Im gM = -(2/3) M (M_D+M)/M_D Im g1 ,  Im gE = (2/3) M (M_D+M)/M_D (Im g1 - M_D Im g2).
! Cross-checks (work/code_v2/tools): identical to a 60-digit mpmath evaluation; at W=M_D equal
! to the unitarity (Watson) construction Im G = M_{1+},E_{1+}(Born, I=3/2) x sqrt(k M_D Gamma/F)
! from the code's own Delta-pole partial waves (Im gM 0.2166, Im gE -2.201 for LMP2015 constants).
! Units: W, Md in nucleon-mass units; hA, gA, fpi, pionm from COMMON /Consts/ (M units).
!=============================================================================
    SUBROUTINE gnd_vertex_im(Win,Md,ImgM,ImgE)
      implicit none
      REAL*8, INTENT(in) :: Win, Md
      REAL*8, INTENT(out) :: ImgM, ImgE
      REAL*8 gA,hA,fpi,pi0mass,pionm
      COMMON/Consts/hA,gA,fpi,pi0mass,pionm
      REAL*16 W,rt,mut,rD,lam,bet,alp,kap,bg,wpi,wN,Lp,LN,f,Q,PI16,K,c
      REAL*16 c1p_pi,d1p_pi,d2_pi,d2p_pi
      REAL*16 a1N,a2N,a3N,b1N,b2N,b3N,c1N,c3N,c1pN,c2N,d1N,d1pN,d2N,d2pN
      REAL*16 g1a,g2a,g1b,g2b,g1,g2,MD16
      REAL*8, SAVE :: lastW=-1.d0, lastMd=-1.d0, lastk(5)=0.d0, lastM=0.d0, lastE=0.d0
      INTEGER ivxmode
      COMMON /v2vx/ ivxmode
      REAL*8 aem_c
      COMMON /v2aem/ aem_c
      REAL*8 w_Wm, w_kpi, w_Eqp, w_om, w_aW, w_kap, w_MB, w_EB, w_sq, EDEM, w_ae
      EXTERNAL EDEM
      ImgM=0.d0; ImgE=0.d0
      IF(Win.LE.1.d0+pionm*(1.d0+1.d-9)) RETURN
      IF(ivxmode.EQ.2) THEN
         !* Watson/unitarity construction (exact Watson's theorem for the code's own off-shell Delta
         !* pole): the s-channel Delta-pole Compton multipole is f_MM^{1+} = -a(W) G_M*(W)^2/D(W) with
         !* a(W) = alpha w^2 (W+M)^2 (W+M_D)/(8(M+M_D)^2) (determined from the code's partial waves,
         !* exact) and G_M* = gM - kappa(W) gE, kappa = (W-M)/(2(W+M)); unitarity in P33 then requires
         !*   Im G_M* = M1+^B sqrt(k M_D Gamma/a),  Im gE = (E1+^B/kappa) sqrt(k M_D Gamma/a),
         !* with the LO Born multipoles M1+^B, E1+^B (I=3/2, CG-normalised) of SUBROUTINE born_p33.
         w_ae=aem_c
         IF(w_ae.LE.0.d0) w_ae=1.d0/137.036d0
         w_Wm=Win
         w_Eqp=(w_Wm*w_Wm+pionm**2-1.d0)/(2.d0*w_Wm)
         w_kpi=dsqrt(MAX(w_Eqp**2-pionm**2,0.d0))
         w_om=(w_Wm*w_Wm-1.d0)/(2.d0*w_Wm)
         w_aW=w_ae*w_om**2*(w_Wm+1.d0)**2*(w_Wm+Md)/(8.d0*(1.d0+Md)**2)
         w_kap=(w_Wm-1.d0)/(2.d0*(w_Wm+1.d0))
         CALL born_p33(w_Wm,pionm,gA/fpi,dsqrt(4.d0*4.d0*datan(1.d0)*w_ae),w_MB,w_EB)
         w_sq=dsqrt(w_kpi*Md*EDEM(w_Wm,Md)/w_aW)
         ImgE=w_EB/w_kap*w_sq
         ImgM=w_MB*w_sq+w_kap*ImgE
         RETURN
      ENDIF
      IF(Win.EQ.lastW.AND.Md.EQ.lastMd.AND.lastk(1).EQ.hA.AND.lastk(2).EQ.gA.AND.lastk(3).EQ.fpi &
           .AND.lastk(4).EQ.pionm.AND.lastk(5).EQ.DBLE(ivxmode)) THEN
         ImgM=lastM; ImgE=lastE
         RETURN
      ENDIF
      PI16=4.q0*ATAN(1.q0)
      W=Win
      MD16=Md
      rt=1.q0/W; mut=pionm/W; rD=MD16/W
      Q=1.q-14
      lam=SQRT(((1+rt)**2-mut**2)*((1-rt)**2-mut**2))/2
      bet=(1-rt**2+mut**2)/2
      alp=1-bet
      kap=SQRT(((1+rt)**2+Q)*((1-rt)**2+Q))/2
      bg=(1-rt**2-Q)/2
      wpi=2*kap*lam/(Q+2*bet*bg)
      wN=2*kap*lam/(Q+2*alp*bg)
      Lp=LOG((1+wpi)/(1-wpi))
      LN=LOG((1+wN)/(1-wN))
      f=1/(4*PI16)
      !* graph (a): photon couples to the pion
      c1p_pi=f/(4*kap**2)*(lam*((1-3*Q/kap**2)*(2*bet*bg+Q)+2*Q) &
           +1/(2*kap)*Lp*(lam**2*(2*bg**2*(1-1/wpi**2)+4*Q/wpi**2)-Q*(2*bet*bg+Q)))
      d1p_pi=f*bg/(8*kap**4)*(lam*(2*lam**2*(2.q0/3*(bg**2-Q)+5*Q/wpi**2)-Q*(4*bet**2+2*Q+2*bet*bg)) &
           +Q/(2*kap)*Lp*(lam**2*(bg*(2*bet-bg)*(3-5/wpi**2)+kap**2*(1-1/wpi**2)) &
           +bet*(2*bet-bg)*(2*bet*bg+Q)))
      d2_pi=f*lam**3/(2*kap**2)*(2.q0/3-1/wpi**2+1/(2*wpi)*Lp*(-1+1/wpi**2))
      d2p_pi=f/(8*kap**4)*(lam*(2*lam**2*(bg**2*(-4.q0/3+2/wpi**2)-3*Q/wpi**2)+1.5q0*Q*(2*bet*bg+Q)) &
           +lam**2/(2*kap)*Lp*(4*kap**2*bet*bg*(1-1/wpi**2)-bg*Q*(2*bet-bg)*(3-5/wpi**2)))
      !* graph (b): photon couples to the nucleon
      a1N=f*lam*alp
      a2N=f*(-lam**3/3)
      a3N=f*lam*(rt**2+4*lam**2/3)
      b1N=f/(2*kap)*LN
      b2N=f/(4*kap**2)*(4*lam*bg+Q*(1+rt**2+Q+2*(rt**2-mut**2))/(2*kap)*LN)
      b3N=f/(2*kap**2)*(-2*lam+(2*alp*bg+Q)/(2*kap)*LN)
      c1N=1/(4*kap**2)*(-2*bg*a1N+2*(bg**2*rt**2-Q*lam**2)*b1N+Q*bg*b2N+Q*(kap**2+bg**2)*b3N)
      c3N=1/kap**2*(3*c1N+lam**2*b1N)
      c1pN=-(1-3*Q/kap**2)*c1N+Q*lam**2/kap**2*b1N-Q/2*b3N
      c2N=-2*c1N-Q*c3N+rt**2*b1N+Q*b3N
      d1N=1/(2*kap**2)*(-bg*(2*a2N+a3N)+(rt**2*(kap**2-Q)-Q*lam**2)*b2N &
           +Q/(2*bg)*(2*kap**2-Q)*c1pN+Q/2*bg*(c2N+2*c1N))
      d1pN=-2.5q0*d1N+rt**2/2*b2N-(1-rt**2)/(4*kap**2)*a3N &
           +Q/(2*kap**2)*(5*d1N+a3N/2-rt**2*b2N+bg/2*(c1pN+c2N))
      d2N=bg*d1N/Q+c1N/2-a2N/Q
      d2pN=bg*(d1N+2*d1pN)/(2*Q)+c1pN/4
      !* Im g1, Im g2 (nucleon-mass units; M_N = 1)
      K=REAL(gA,16)*REAL(hA,16)/(4*REAL(fpi,16)**2)
      g1a=K*d2_pi
      g2a=-K/W*((rD+rt)/bg*(2*d2p_pi-c1p_pi)+2*rD/bg*d1p_pi)
      g1b=-K*(d2N+rt*(rD-rt)/2*b3N+(1-rt*rD)/(2*bg)*c1pN-Q/2*c3N)
      g2b=K/W*(-rt*b3N-(rD-rt)/bg*c1pN+2*rD/bg*d1pN+2*(rD+rt)/bg*d2pN)
      g1=g1a+g1b
      g2=g2a+g2b
      !* mapping g1,g2 -> gM,gE: on-shell (M_D) relations [ivxmode=0, default] or with M_D -> W
      !* [ivxmode=1]; the two differ only off the resonance peak (by up to ~10% in Im gE)
      IF(ivxmode.EQ.1) MD16=W
      c=2.q0/3*(MD16+1)/MD16
      ImgM=REAL(-c*g1,8)
      ImgE=REAL(c*(g1-MD16*g2),8)
      lastW=Win; lastMd=Md; lastk=(/hA,gA,fpi,pionm,DBLE(ivxmode)/); lastM=ImgM; lastE=ImgE
    END SUBROUTINE gnd_vertex_im

!=============================================================================
! [v2] running couplings / mass / field renormalisation of the Delta pole (&v2 drunning, irun)
!=============================================================================
    SUBROUTINE delta_running(s,Md,G1,G2,ReG1,ReG2,Mrun,ztau)
      !* [v2] higher-order Delta-pole option (&v2 drunning / irun), LMP2015 Sect. 2.3 and App. D:
      !*  Re gM(s), Re gE(s): dispersion relation of Im gM, Im gE (gnd_vertex_im), principal values by
      !*     singularity subtraction; nodes s' = s_th + c (y/(1-y))^2, 400 Gauss points (cached).
      !*  M_D(s) = M_D + Re sigma_R(s),  ztau = 1 + Re tau_R(s): piN (C_N = hA^2) and piDelta
      !*     (C_D = (5/3 H_A)^2, H_A = 9/5 gA) self-energy loops, Eqs. (D1)-(D4), Xi = 0, scale M_N.
      !* Subtraction schemes (COMMON /v2drun/ irun_c; main program: &v2 irun, or drunning=t -> irun=1):
      !*  irun=1  once subtracted at s = M_D^2 (LMP2015 Eq. D5 and D4): Re g = G + F1(s),
      !*          F1(s) = (s-M_D^2)/pi P int ds' Im g(s')/((s'-M_D^2)(s'-s)); sigma_R = sigma - Re sigma(M_D^2)
      !*          - (s-M_D^2) Re sigma'(M_D^2); tau_R = tau - Re tau(M_D^2).  Changes the static polarisabilities
      !*          strongly (Delta-pole beta_M1 7.15 -> 9.35); for studies only.
      !*  irun=2  doubly subtracted: additionally subtracted at the static point s0 = M_N^2 (= 1 here), so that
      !*          Re g(s0) = G, sigma_R(s0) = 0, tau_R(s0) = 0 and alpha_E1, beta_M1 are unchanged:
      !*            Re g(s) = G + (s-M_D^2)(s-s0)/pi P int ds' Im g(s')/((s'-s)(s'-M_D^2)(s'-s0))
      !*                    = G + F1(s) - (s-M_D^2)/(s0-M_D^2) F1(s0)             [partial fractions in s';
      !*            the twice-subtracted integral vanishes at both subtraction points, hence no further linear term]
      !*            sigma_R(s) = sigma(s) - [a + b (s-M_D^2) + c (s-M_D^2)^2], a, b as for irun=1 (value and slope
      !*            at the pole), c = sigma_R^(1)(s0)/(s0-M_D^2)^2 ;  tau_R(s) = tau(s) - Re tau(M_D^2) - c' (s-M_D^2),
      !*            c' = tau_R^(1)(s0)/(s0-M_D^2).
      !*          The subtracted pieces are polynomials in s (analytic): they correspond to counterterms of
      !*          higher order, the leading ones being the O(p^4) contact terms (delta alpha, delta beta, delta gamma)
      !*          of the analysis, which absorb them.  The spin polarisabilities still depend on the slopes of the
      !*          running functions at s0 (they enter the amplitude at O(nu^3)).
      !*  irun=3  (exploratory) as irun=2 with, in addition, vanishing slopes at s0 (Re g'(s0) = 0, sigma_R'(s0) = 0,
      !*          tau_R'(s0) = 0), i.e. one more polynomial term each: -(s-M_D^2)(s-s0) b, b = (F1'(s0)-a)/(s0-M_D^2),
      !*          a = F1(s0)/(s0-M_D^2), F1'(s0) = 1/pi int ds' Im g(s')/(s'-s0)^2; -d (s-M_D^2)^2 (s-s0),
      !*          d = sigma_R^(2)'(s0)/(s0-M_D^2)^2; -d' (s-M_D^2)(s-s0), d' = tau_R^(2)'(s0)/(s0-M_D^2).  Then all
      !*          six static dipole polarisabilities are unchanged.
      !* Units: nucleon-mass units. Verified against tools/running_delta.py (schemes 1-3, 1e-6) and, for irun=1,
      !* LMP2015 Fig. 5 (M_D(s)/M_N and 1+Re tau_R reproduced; gM(s) to 0.02).
      implicit none
      REAL*8, INTENT(in) :: s,Md,G1,G2
      REAL*8, INTENT(out) :: ReG1,ReG2,Mrun,ztau
      REAL*8 gA,hA,fpi,pi0mass,pionm
      COMMON/Consts/hA,gA,fpi,pi0mass,pionm
      INTEGER, PARAMETER :: NQ=400
      REAL*8, PARAMETER :: s0=1.d0        ! static point s0 = M_N^2 (nucleon-mass units)
      REAL*8, SAVE :: yq(NQ), wq(NQ), sp(NQ), jac(NQ), imgq(NQ,2), imgD(2), key(6)=-1.d0
      REAL*8, SAVE :: s0N(3), s0D(3), F0(2), dF0(2), sR0, tR0, dsR0, dtR0
      REAL*8 sth, cc, PI, ims(2), G(2), F(2), wS, wM, MD2, sig(2), tau(2), dsg(2), dtg(2), sR, tR
      REAL*8 dd, ca(2), cb(2), cm, cmp, dm, dmp
      INTEGER k
      INTEGER ivxmode
      COMMON /v2vx/ ivxmode
      INTEGER irun_c
      COMMON /v2drun/ irun_c
      PI=4.d0*datan(1.d0)
      sth=(1.d0+pionm)**2
      cc=1.d0
      MD2=Md**2
      IF(key(1).NE.Md.OR.key(2).NE.hA.OR.key(3).NE.gA.OR.key(4).NE.fpi.OR.key(5).NE.pionm &
           .OR.key(6).NE.DBLE(ivxmode)) THEN
         CALL gauleg_v2(0.d0,1.d0,yq,wq,NQ)
         DO k=1,NQ
            sp(k)=sth+cc*(yq(k)/(1.d0-yq(k)))**2
            jac(k)=cc*2.d0*yq(k)/(1.d0-yq(k))**3*wq(k)
            CALL gnd_vertex_im(dsqrt(sp(k)),Md,imgq(k,1),imgq(k,2))
         ENDDO
         CALL gnd_vertex_im(Md,Md,imgD(1),imgD(2))
         CALL delta_selfenergy(MD2,Md,1,s0N(1),s0N(2),s0N(3),dtg(1))
         CALL delta_selfenergy(MD2,Md,2,s0D(1),s0D(2),s0D(3),dtg(2))
         !* once-subtracted integral F1 and its slope at the static point s0 (below threshold: no principal value)
         F0=0.d0; dF0=0.d0
         DO k=1,NQ
            wM=(MD2-sth+cc)/(sp(k)-sth+cc)
            G=imgq(k,:)/(sp(k)-s0)-(imgq(k,:)-imgD*wM)/(sp(k)-MD2)
            F0=F0+G*jac(k)
            dF0=dF0+imgq(k,:)/(sp(k)-s0)**2*jac(k)
         ENDDO
         F0=(F0+imgD*dlog((MD2-sth)/cc))/PI
         dF0=dF0/PI
         !* once-subtracted sigma_R, tau_R and their slopes at s0
         CALL delta_selfenergy(s0,Md,1,sig(1),tau(1),dsg(1),dtg(1))
         CALL delta_selfenergy(s0,Md,2,sig(2),tau(2),dsg(2),dtg(2))
         sR0=(sig(1)-s0N(1)-(s0-MD2)*s0N(3))+(sig(2)-s0D(1)-(s0-MD2)*s0D(3))
         tR0=(tau(1)-s0N(2))+(tau(2)-s0D(2))
         dsR0=(dsg(1)-s0N(3))+(dsg(2)-s0D(3))
         dtR0=dtg(1)+dtg(2)
         key=(/Md,hA,gA,fpi,pionm,DBLE(ivxmode)/)
      ENDIF
      !* running couplings: once-subtracted integral F1(s) (principal value by singularity subtraction with the
      !* weight w_x(s') = (x-s_th+c)/(s'-s_th+c), P int w_x/(s'-x) = -log((x-s_th)/c))
      ims=0.d0
      IF(s.GT.sth) CALL gnd_vertex_im(dsqrt(s),Md,ims(1),ims(2))
      F=0.d0
      DO k=1,NQ
         wM=(MD2-sth+cc)/(sp(k)-sth+cc)
         IF(s.GT.sth) THEN
            wS=(s-sth+cc)/(sp(k)-sth+cc)
            G=(imgq(k,:)-ims*wS)/(sp(k)-s)
         ELSE
            G=imgq(k,:)/(sp(k)-s)
         ENDIF
         G=G-(imgq(k,:)-imgD*wM)/(sp(k)-MD2)
         F=F+G*jac(k)
      ENDDO
      IF(s.GT.sth) F=F-ims*dlog((s-sth)/cc)
      F=F+imgD*dlog((MD2-sth)/cc)
      F=F/PI
      !* running mass and field renormalisation: sigma_R = sigma - Re sigma(MD^2) - (s-MD^2) Re sigma'(MD^2)
      CALL delta_selfenergy(s,Md,1,sig(1),tau(1),dsg(1),dtg(1))
      CALL delta_selfenergy(s,Md,2,sig(2),tau(2),dsg(2),dtg(2))
      sR=(sig(1)-s0N(1)-(s-MD2)*s0N(3))+(sig(2)-s0D(1)-(s-MD2)*s0D(3))
      tR=(tau(1)-s0N(2))+(tau(2)-s0D(2))
      IF(irun_c.GE.2) THEN
         !* second subtraction at the static point s0 (polynomial pieces = higher-order counterterms, see above)
         dd=s0-MD2
         ca=F0/dd
         F=F-(s-MD2)*ca
         cm=sR0/dd**2; cmp=tR0/dd
         sR=sR-cm*(s-MD2)**2; tR=tR-cmp*(s-MD2)
         IF(irun_c.GE.3) THEN
            !* third condition: vanishing slopes at s0
            cb=(dF0-ca)/dd
            F=F-(s-MD2)*(s-s0)*cb
            dm=(dsR0-2.d0*cm*dd)/dd**2; dmp=(dtR0-cmp)/dd
            sR=sR-dm*(s-MD2)**2*(s-s0); tR=tR-dmp*(s-MD2)*(s-s0)
         ENDIF
      ENDIF
      ReG1=G1+F(1); ReG2=G2+F(2)
      Mrun=Md+sR
      ztau=1.d0+tR
    END SUBROUTINE delta_running

    SUBROUTINE delta_selfenergy(s,Md,which,ReSig,ReTau,ReDSig,ReDTau)
      !* Re sigma(s), Re tau(s), Re dsigma/ds, Re dtau/ds of LMP2015 Eq. (D2) for the piN (which=1) or piDelta
      !* (which=2) loop; Xi=0, log scale M_N=1; Feynman-parameter integral split at the zeros of M^2(x)
      implicit none
      REAL*8, INTENT(in) :: s, Md
      INTEGER, INTENT(in) :: which
      REAL*8, INTENT(out) :: ReSig, ReTau, ReDSig, ReDTau
      REAL*8 gA,hA,fpi,pi0mass,pionm
      COMMON/Consts/hA,gA,fpi,pi0mass,pionm
      INTEGER, PARAMETER :: NX=96
      REAL*8, SAVE :: xg(NX), wg(NX)
      LOGICAL, SAVE :: init=.true.
      REAL*8 C, Mx, pref, a, b, c2, disc, r1, r2, pts(4), x, w, M2, L, PI, xa, xb
      INTEGER np, i, k
      PI=4.d0*datan(1.d0)
      IF(init) THEN
         CALL gauleg_v2(-1.d0,1.d0,xg,wg,NX); init=.false.
      ENDIF
      IF(which.EQ.1) THEN
         C=hA**2; Mx=1.d0
      ELSE
         C=(5.d0/3.d0*9.d0/5.d0*gA)**2; Mx=Md
      ENDIF
      pref=-C/(2.d0*(8.d0*PI*fpi)**2)
      a=s; b=-(s-pionm**2+Mx**2); c2=Mx**2
      disc=b*b-4.d0*a*c2
      np=2; pts(1)=0.d0; pts(2)=1.d0
      IF(disc.GT.0.d0.AND.a.NE.0.d0) THEN
         r1=(-b-dsqrt(disc))/(2.d0*a); r2=(-b+dsqrt(disc))/(2.d0*a)
         IF(r1.GT.0.d0.AND.r1.LT.1.d0) THEN; np=np+1; pts(np)=r1; ENDIF
         IF(r2.GT.0.d0.AND.r2.LT.1.d0) THEN; np=np+1; pts(np)=r2; ENDIF
      ENDIF
      DO i=1,np-1                     ! sort
         DO k=i+1,np
            IF(pts(k).LT.pts(i)) THEN; x=pts(i); pts(i)=pts(k); pts(k)=x; ENDIF
         ENDDO
      ENDDO
      ReSig=0.d0; ReTau=0.d0; ReDSig=0.d0; ReDTau=0.d0
      DO i=1,np-1
         xa=pts(i); xb=pts(i+1)
         IF(xb-xa.LT.1.d-13) CYCLE     ! coincident roots (exactly at threshold): 0 * log(0)
         DO k=1,NX
            !* cubic endpoint map x = mid + half (3u - u^3)/2 (vanishing Jacobian at the ends): the log
            !* singularities of the integrand at the zeros of M^2(x) then converge like a smooth integrand
            !* (sigma'(M_D^2) accurate to 1e-8 with 96 nodes; plain Gauss-Legendre: 2e-5)
            x=0.5d0*(xb+xa)+0.5d0*(xb-xa)*(3.d0*xg(k)-xg(k)**3)/2.d0
            w=0.5d0*(xb-xa)*1.5d0*(1.d0-xg(k)**2)*wg(k)
            M2=pionm**2*x+Mx**2*(1.d0-x)-s*x*(1.d0-x)
            L=dlog(dabs(M2))                   ! real part of log(M^2 - i0)
            ReSig=ReSig+w*(x*Md+Mx)*M2*(-1.d0+L)
            ReTau=ReTau+w*x*M2*(-1.d0+L)
            ReDSig=ReDSig+w*(x*Md+Mx)*(-x*(1.d0-x))*L      ! d/ds [M2 (-1+L)] = dM2/ds L, dM2/ds = -x(1-x)
            ReDTau=ReDTau+w*x*(-x*(1.d0-x))*L
         ENDDO
      ENDDO
      ReSig=pref*ReSig; ReTau=pref*ReTau; ReDSig=pref*ReDSig; ReDTau=pref*ReDTau
    END SUBROUTINE delta_selfenergy

    SUBROUTINE gauleg_v2(x1,x2,x,w,n)
      !* Gauss-Legendre nodes/weights (Numerical Recipes algorithm), independent of gauss.f90
      implicit none
      INTEGER, INTENT(in) :: n
      REAL*8, INTENT(in) :: x1,x2
      REAL*8, INTENT(out) :: x(n),w(n)
      INTEGER i,j,m
      REAL*8 p1,p2,p3,pp,xl,xm,z,z1,PI
      PI=4.d0*datan(1.d0)
      m=(n+1)/2; xm=0.5d0*(x2+x1); xl=0.5d0*(x2-x1)
      DO i=1,m
         z=dcos(PI*(i-0.25d0)/(n+0.5d0))
         DO
            p1=1.d0; p2=0.d0
            DO j=1,n
               p3=p2; p2=p1; p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
            ENDDO
            pp=n*(z*p1-p2)/(z*z-1.d0)
            z1=z; z=z1-p1/pp
            IF(dabs(z-z1).LT.1.d-15) EXIT
         ENDDO
         x(i)=xm-xl*z; x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1.d0-z*z)*pp*pp); w(n+1-i)=w(i)
      ENDDO
    END SUBROUTINE gauleg_v2

!=============================================================================
! [v2] Tree-level (LO: charge couplings, pseudoscalar = pseudovector piNN at kappa=0) pion
! photoproduction multipoles M1+ and E1+ in isospin 3/2, CG-normalised (= sqrt(2/3) A^{3/2},
! A^{3/2} = A(pi0 p) - A(pi+ n)/sqrt2, MAID/CGLN phase convention: E0+(pi+ n) > 0), in nucleon-
! mass units. Explicit Dirac algebra, CGLN decomposition by traces, Gauss-Legendre projection.
! Fortran port of tools/born_multipoles.py (verified there: gauge invariance, LETs, unitarity of
! the code's piN loops). Used by the Watson construction of the vertex Im parts (vxmode=2).
!   T(pi0 p) = -i e [ V S(p+k) gamma.eps + gamma.eps S(p-q) V ]
!   T(pi+ n) = -i e sqrt2 [ V S(p+k) gamma.eps + V (2q-k).eps / ((q-k)^2 - m^2) ],  V = g gamma5
!=============================================================================
    SUBROUTINE born_p33(W,mpi,gpiN,e,M1p,E1p)
      implicit none
      REAL*8, INTENT(in) :: W, mpi, gpiN, e
      REAL*8, INTENT(out) :: M1p, E1p
      INTEGER, PARAMETER :: NG=32
      COMPLEX*16 G(4,4,0:3), G5(4,4), Id(4,4), Ss(4,4), Su(4,4), V(4,4), Op(4,4,2,2), T(2,2,2,2)
      COMPLEX*16 uf(4,2), ui(4,2), ubf(4,2), F(4,2), Fx(2,2), Fy(2,2), CI, tr1, tr2, tr3, tr4
      REAL*8 xg(NG), wg(NG), s, k, Eq, q, c, sn, kv(0:3), qv(0:3), pv(0:3), ppv(0:3), tv(0:3), pole
      REAL*8 PI, P0, P1, P2, P3, Ml(2), El(2), mt(0:3)
      INTEGER ig, ch, ie, mu, i, j
      LOGICAL, SAVE :: init=.true.
      REAL*8, SAVE :: xs(NG), ws(NG)
      CI=(0.d0,1.d0); PI=4.d0*datan(1.d0)
      IF(init) THEN
         CALL gauleg_v2(-1.d0,1.d0,xs,ws,NG); init=.false.
      ENDIF
      mt=(/1.d0,-1.d0,-1.d0,-1.d0/)
      CALL dirac_bd(G,G5)
      Id=0.d0
      DO i=1,4
         Id(i,i)=1.d0
      ENDDO
      s=W*W
      k=(s-1.d0)/(2.d0*W)
      Eq=(s+mpi**2-1.d0)/(2.d0*W)
      q=dsqrt(MAX(Eq*Eq-mpi**2,0.d0))
      Ml=0.d0; El=0.d0
      V=gpiN*G5
      DO ig=1,NG
         c=xs(ig); sn=dsqrt(1.d0-c*c)
         kv=(/k,0.d0,0.d0,k/); qv=(/Eq,q*sn,0.d0,q*c/)
         pv=(/dsqrt(k*k+1.d0),0.d0,0.d0,-k/); ppv=(/dsqrt(q*q+1.d0),-q*sn,0.d0,-q*c/)
         CALL prop_inv(pv+kv,G,Id,Ss)
         CALL prop_inv(pv-qv,G,Id,Su)
         tv=qv-kv
         pole=1.d0/(tv(0)**2-tv(1)**2-tv(2)**2-tv(3)**2-mpi**2)
         CALL spinors_bd(pv,ui); CALL spinors_bd(ppv,uf)
         DO j=1,2                       ! ubar = u^dagger gamma0
            DO i=1,4
               ubf(i,j)=DCONJG(uf(i,j))*G(i,i,0)
            ENDDO
         ENDDO
         DO ch=1,2                      ! 1 = pi0 p, 2 = pi+ n
            DO ie=1,2                   ! eps = x, y  (eps_mu = -eps^i)
               mu=ie
               IF(ch.EQ.1) THEN
                  Op(:,:,ie,ch)=MATMUL(MATMUL(V,Ss),G(:,:,mu))+MATMUL(MATMUL(G(:,:,mu),Su),V)
               ELSE
                  Op(:,:,ie,ch)=dsqrt(2.d0)*(MATMUL(MATMUL(V,Ss),G(:,:,mu))+V*(2.d0*qv(mu)-kv(mu))*pole)
               ENDIF
               Op(:,:,ie,ch)=-CI*e*(-1.d0)*Op(:,:,ie,ch)        ! overall -i e, eps_mu = -eps^mu
               DO i=1,2
                  DO j=1,2
                     T(i,j,ie,ch)=DOT_PRODUCT(DCONJG(ubf(:,i)),MATMUL(Op(:,:,ie,ch),ui(:,j)))/(8.d0*PI*W)
                  ENDDO
               ENDDO
            ENDDO
            Fx=T(:,:,1,ch); Fy=T(:,:,2,ch)
            ! traces with sigma matrices: Tr(sx M)=M21+M12, Tr(sy M)=i(M12-M21)... sy=[[0,-i],[i,0]]
            tr1=Fy(1,1)+Fy(2,2)                          ! Tr(Fy)
            tr2=-CI*Fy(2,1)+CI*Fy(1,2)                   ! Tr(sy Fy) = sum_ij sy_ij Fy_ji
            tr3=Fx(2,1)+Fx(1,2)                          ! Tr(sx Fx)
            tr4=Fx(1,1)-Fx(2,2)                          ! Tr(sz Fx)
            F(2,ch)=-tr1/(2.d0*sn)
            F(1,ch)=tr2/(2.d0*CI)+c*F(2,ch)
            F(4,ch)=(tr3/(2.d0*CI)-(F(1,ch)-c*F(2,ch)))/sn**2
            F(3,ch)=tr4/(2.d0*CI*sn)-F(2,ch)-c*F(4,ch)
            P0=1.d0; P1=c; P2=0.5d0*(3.d0*c*c-1.d0); P3=0.5d0*(5.d0*c**3-3.d0*c)
            Ml(ch)=Ml(ch)+ws(ig)*DREAL(P1*F(1,ch)-P2*F(2,ch)-(P0-P2)/3.d0*F(3,ch))/4.d0
            El(ch)=El(ch)+ws(ig)*DREAL(P1*F(1,ch)-P2*F(2,ch)+(P0-P2)/3.d0*F(3,ch) &
                 +2.d0/5.d0*(P1-P3)*F(4,ch))/4.d0
         ENDDO
      ENDDO
      M1p=dsqrt(2.d0/3.d0)*(Ml(1)-Ml(2)/dsqrt(2.d0))
      E1p=dsqrt(2.d0/3.d0)*(El(1)-El(2)/dsqrt(2.d0))
    END SUBROUTINE born_p33

    SUBROUTINE dirac_bd(G,G5)
      implicit none
      COMPLEX*16, INTENT(out) :: G(4,4,0:3), G5(4,4)
      COMPLEX*16 sig(2,2,3), CI
      INTEGER k
      CI=(0.d0,1.d0)
      sig=0.d0
      sig(1,2,1)=1.d0; sig(2,1,1)=1.d0
      sig(1,2,2)=-CI; sig(2,1,2)=CI
      sig(1,1,3)=1.d0; sig(2,2,3)=-1.d0
      G=0.d0; G5=0.d0
      G(1,1,0)=1.d0; G(2,2,0)=1.d0; G(3,3,0)=-1.d0; G(4,4,0)=-1.d0
      DO k=1,3
         G(1:2,3:4,k)=sig(:,:,k)
         G(3:4,1:2,k)=-sig(:,:,k)
      ENDDO
      G5(1,3)=1.d0; G5(2,4)=1.d0; G5(3,1)=1.d0; G5(4,2)=1.d0
    END SUBROUTINE dirac_bd

    SUBROUTINE prop_inv(p,G,Id,S)
      !* S = (pslash - 1)^-1 = (pslash + 1)/(p^2 - 1)  (nucleon mass 1)
      implicit none
      REAL*8, INTENT(in) :: p(0:3)
      COMPLEX*16, INTENT(in) :: G(4,4,0:3), Id(4,4)
      COMPLEX*16, INTENT(out) :: S(4,4)
      REAL*8 p2
      p2=p(0)**2-p(1)**2-p(2)**2-p(3)**2
      S=(p(0)*G(:,:,0)-p(1)*G(:,:,1)-p(2)*G(:,:,2)-p(3)*G(:,:,3)+Id)/(p2-1.d0)
    END SUBROUTINE prop_inv

    SUBROUTINE spinors_bd(p,u)
      !* u(p,chi) for chi=(1,0),(0,1), ubar u = 2M (M=1)
      implicit none
      REAL*8, INTENT(in) :: p(0:3)
      COMPLEX*16, INTENT(out) :: u(4,2)
      REAL*8 a
      COMPLEX*16 CI
      CI=(0.d0,1.d0)
      a=dsqrt(p(0)+1.d0)
      u=0.d0
      u(1,1)=a; u(2,2)=a
      ! lower components sigma.p chi/(E+M)
      u(3,1)=p(3)/a;            u(4,1)=(p(1)+CI*p(2))/a
      u(3,2)=(p(1)-CI*p(2))/a;  u(4,2)=-p(3)/a
    END SUBROUTINE spinors_bd
