!!****************************************************************************
!!* Module: gauss.mod               Version: 1.2     Update: 05/22/98
!!*
!!*___________________________________________________________________________

        MODULE gauss
          
  !! np = total # of spatial integration points
  !! np0 = total # of temporal integration points
  !! ncosx = # of the angle-integration points
  !! cp, mapp -  parameters of mapping the infty to unit interval

      REAL*8, ALLOCATABLE :: eq(:),p(:),dp(:)
      REAL*8, ALLOCATABLE :: k0(:),dk0(:)
      REAL*8 cp
      INTEGER npa,nph,np,np0,mapp,ncosx
      REAL*8, ALLOCATABLE :: cosx(:),dcosx(:),phi(:),dphi(:)
      LOGICAL bnd 
       
       CONTAINS

      SUBROUTINE space_points(phat)
      IMPLICIT NONE
      REAL*8, INTENT(in) :: phat   ! singularity position
      REAL*8, ALLOCATABLE :: x(:),dx(:)
      REAL*8 zero,one,tphat,xv,xp,xpm1,denom
      REAL*8, SAVE :: lastphat
      INTEGER npend,i,ip,nt
      nt=nph+npa
      IF(.NOT.ALLOCATED(p)) THEN
         ALLOCATE(eq(0:nt),p(0:nt),dp(0:nt))
         lastphat=phat
      ELSE
         IF(phat.EQ.lastphat) RETURN
      ENDIF
      lastphat=phat
      npend=npa
      np=npend
      zero=0.d0
      IF(phat.GT.0.d0) THEN
         !!       symmetric integration points on [zero, 2*phat]
         tphat=2.d0*phat 
         CALL gauleg(zero,tphat,p(1),dp(1),nph)
         np=nt
      ENDIF
      IF(npend.NE.0) THEN
         !!       integration points on [tphat, infinity)
         one=1.d0
         ALLOCATE(x(npend),dx(npend))
         CALL gauleg(zero,one,x,dx,npend)
         !!       mapping from [tphat,infty) to [0,1] 
         IF (mapp.EQ.0) THEN
            ip=np+1
            DO i=1,npend
               ip=ip-1
               p(ip)=tphat-cp*dlog(x(i))
               dp(ip)=cp*dx(i)/x(i)
            ENDDO
         ELSE
            DO i=1,npend
               xv=x(i)
               xp=x(i)**mapp
               xpm1=xp/xv
               denom=1.d0/(1.d0-xv)
               ip = i+nph
               p(ip)=cp*xp*denom**mapp + tphat
               dp(ip)=cp*mapp*xpm1*dx(i)*denom**(mapp+1)
            ENDDO
         ENDIF
         DEALLOCATE(x,dx)
      ENDIF
      !!     on-shell point
      dp(0)=1.d0
      p(0)=phat
      !!     energy of a particle with unit mass
      DO i=0,np 
         eq(i)=dsqrt(1.d0+p(i)**2)
      ENDDO
      END SUBROUTINE space_points
      
      SUBROUTINE time_points
        IMPLICIT NONE
        REAL*8, ALLOCATABLE :: x(:),dx(:)
        REAL*8 xv,xp,xpm1,denom
        INTEGER npend,i,ip
        
        npend=np0/2
        IF(.NOT.ALLOCATED(k0)) ALLOCATE(k0(0:np0),dk0(0:np0),x(npend),dx(npend))
        CALL gauleg(0.d0,1.d0,x,dx,npend)
        
        DO ip=1,npend
           xv=x(ip)
           xp=x(ip)**mapp
           xpm1=xp/xv
           denom=1.d0/(1.d0-xv)
           k0(npend-ip+1)=-(cp*xp*denom**mapp)
           dk0(npend-ip+1)=cp*mapp*xpm1*dx(ip)*denom**(mapp+1)
           k0(ip+npend)=-k0(npend-ip+1)
           dk0(ip+npend)=dk0(npend-ip+1)
        ENDDO
        DEALLOCATE(x,dx)
        k0(0)=0.d0
        dk0(0)=1.d0
      END SUBROUTINE time_points
      
      SUBROUTINE angular_points
        IMPLICIT NONE 
        REAL*8 Pi,twoPi
        COMMON/numbers/Pi
              IF(.NOT.ALLOCATED(cosx))   ALLOCATE(cosx(ncosx),dcosx(ncosx),phi(ncosx),dphi(ncosx))

        CALL gauleg(-1.d0,1.d0,cosx,dcosx,ncosx)
        twoPi=2.d0*Pi      
        CALL gauleg(0.d0,twoPi,phi,dphi,ncosx)
      END SUBROUTINE angular_points
      
      END MODULE gauss
       

      SUBROUTINE gauleg(x1,x2,x,w,n)
        INTEGER, INTENT(in) :: n
        REAL*8, INTENT(in) :: x1,x2
        REAL*8, INTENT(out) :: x(n),w(n)
        DOUBLE PRECISION EPS
        PARAMETER (EPS=3.d-14)
        INTEGER i,j,m
        DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
        m=(n+1)/2
        xm=0.5d0*(x2+x1)
        xl=0.5d0*(x2-x1)
        DO i=1,m
           z=dcos(3.141592654d0*(i-.25d0)/(n+.5d0))
1          CONTINUE
           p1=1.d0
           p2=0.d0
           DO j=1,n
              p3=p2
              p2=p1
              p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
           ENDDO
           pp=n*(z*p1-p2)/(z*z-1.d0)
           z1=z
           z=z1-p1/pp
           IF(dabs(z-z1).GT.EPS) GOTO 1
           x(i)=xm-xl*z
           x(n+1-i)=xm+xl*z
           w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
           w(n+1-i)=w(i)
        ENDDO
        RETURN
      END SUBROUTINE gauleg
