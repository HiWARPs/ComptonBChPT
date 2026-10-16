!!****************************************************************************
!!* Module: dirac.f90               Version: 1.3     Update: 08/01/98
!!*
!!*___________________________________________________________________________

        MODULE dirac

        TYPE matrix4
        PRIVATE 
        COMPLEX*16 a(4,4)
        END TYPE matrix4
        
        TYPE spinor
        PRIVATE 
        COMPLEX*16 sp(4)
        END TYPE spinor
        
        TYPE (matrix4) gzero, gone, gamm(0:5), gammunu(0:3,0:3)
        INTEGER epsilon4(0:3,0:3,0:3,0:3),kronikers(4,4),D_(0:3,0:3)
        REAL*8 gmunu(0:3)
        COMPLEX*16 veps(0:3,3), vepsp(0:3,3)
        COMPLEX*16, PARAMETER :: ci=(0.d0,1.d0)


        INTERFACE OPERATOR (+)
           MODULE PROCEDURE msum
           MODULE PROCEDURE mr_sum
           MODULE PROCEDURE mc_sum
           MODULE PROCEDURE rm_sum
           MODULE PROCEDURE cm_sum
           MODULE PROCEDURE spsum
        END INTERFACE 
        INTERFACE OPERATOR (-)
           MODULE PROCEDURE msub
           MODULE PROCEDURE mr_sub
           MODULE PROCEDURE mc_sub
           MODULE PROCEDURE rm_sub
           MODULE PROCEDURE cm_sub
           MODULE PROCEDURE spsub
        END INTERFACE 

        INTERFACE OPERATOR (*)
           MODULE PROCEDURE mm_prod
           MODULE PROCEDURE rm_prod
           MODULE PROCEDURE cm_prod
           MODULE PROCEDURE mr_prod
           MODULE PROCEDURE mc_prod
           MODULE PROCEDURE im_prod
           MODULE PROCEDURE mi_prod
           MODULE PROCEDURE rsp_prod
           MODULE PROCEDURE csp_prod
           MODULE PROCEDURE spr_prod
           MODULE PROCEDURE spc_prod
           MODULE PROCEDURE spmtx_prod
           MODULE PROCEDURE mtxsp_prod
       END INTERFACE

       INTERFACE dot
          MODULE PROCEDURE rrdot
           MODULE PROCEDURE rcdot
           MODULE PROCEDURE crdot
           MODULE PROCEDURE ccdot
        END INTERFACE

       INTERFACE ini_polar
          MODULE PROCEDURE ini_polar_virtual
          MODULE PROCEDURE ini_polar_real
        END INTERFACE

        INTERFACE OPERATOR (/)
           MODULE PROCEDURE mr_div
           MODULE PROCEDURE mc_div
        END INTERFACE 


        INTERFACE slash
           MODULE PROCEDURE rslash
           MODULE PROCEDURE cslash
        END INTERFACE

        INTERFACE g_
           MODULE PROCEDURE rslash
           MODULE PROCEDURE cslash
        END INTERFACE

        INTERFACE ag_
           MODULE PROCEDURE rrag2
           MODULE PROCEDURE rcag2
           MODULE PROCEDURE crag2
           MODULE PROCEDURE ccag2
        END INTERFACE

         INTERFACE eps4_
           MODULE PROCEDURE rrrreps4
           MODULE PROCEDURE rrrceps4
           MODULE PROCEDURE rrcceps4
           MODULE PROCEDURE rcrceps4
           MODULE PROCEDURE crcreps4
           MODULE PROCEDURE crrceps4
           MODULE PROCEDURE ccrreps4
           MODULE PROCEDURE rccreps4
        END INTERFACE

       INTERFACE contract
           MODULE PROCEDURE contract0
           MODULE PROCEDURE contract1in
           MODULE PROCEDURE contract1out
           MODULE PROCEDURE contract2
        END INTERFACE
      CONTAINS

        !!--matrix operations:
        !!--
        FUNCTION msum(m1,m2) RESULT(m3)
          TYPE(matrix4), INTENT(in) :: m1,m2
          TYPE(matrix4) :: m3
          m3%a=m1%a+m2%a
        END FUNCTION msum

        FUNCTION mc_SUM(m1,x1) RESULT(m3)
          COMPLEX*16, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m1
          TYPE(matrix4) :: m3
          INTEGER i
          DO i=1,4
             m3%a(i,i)=m1%a(i,i)+x1
          ENDDO
        END FUNCTION mc_sum

        FUNCTION cm_SUM(x1,m1) RESULT(m3)
          COMPLEX*16, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m1
          TYPE(matrix4) :: m3
          INTEGER i
          DO i=1,4
             m3%a(i,i)=m1%a(i,i)+x1
          ENDDO
        END FUNCTION cm_sum

        FUNCTION mr_SUM(m1,x1) RESULT(m3)
          REAL*8, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m1
          TYPE(matrix4) :: m3
          INTEGER i
          DO i=1,4
             m3%a(i,i)=m1%a(i,i)+x1
          ENDDO
        END FUNCTION mr_sum

        FUNCTION rm_SUM(x1,m1) RESULT(m3)
          REAL*8, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m1
          TYPE(matrix4) :: m3
          INTEGER i
          DO i=1,4
             m3%a(i,i)=m1%a(i,i)+x1
          ENDDO
        END FUNCTION rm_sum

        FUNCTION msub(m1,m2) RESULT(m3)
          TYPE(matrix4), INTENT(in) :: m1,m2
          TYPE(matrix4) :: m3
          m3%a=m1%a-m2%a
        END FUNCTION msub

        FUNCTION mc_sub(m1,x1) RESULT(m3)
          COMPLEX*16, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m1
          TYPE(matrix4) :: m3
          INTEGER i
          DO i=1,4
             m3%a(i,i)=m1%a(i,i)-x1
          ENDDO
        END FUNCTION mc_sub

        FUNCTION cm_sub(x1,m1) RESULT(m3)
          COMPLEX*16, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m1
          TYPE(matrix4) :: m3
          m3=x1*gone-m1
        END FUNCTION cm_sub

        FUNCTION mr_sub(m2,x1) RESULT(m3)
          REAL*8, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m2
          TYPE(matrix4) :: m3
          INTEGER i
          DO i=1,4
             m3%a(i,i)=m2%a(i,i)-x1
          ENDDO
        END FUNCTION mr_sub

        FUNCTION rm_sub(x1,m1) RESULT(m3)
          REAL*8, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m1
          TYPE(matrix4) :: m3
          m3=x1*gone-m1
        END FUNCTION rm_sub

        FUNCTION mm_prod(m1,m2) RESULT(m3)
        TYPE(matrix4), INTENT(in) :: m1,m2
        TYPE(matrix4) :: m3
        INTEGER i,j,k
        COMPLEX*16 x
        DO i=1,4
        DO j=1,4
           x=0.d0
           DO k=1,4
              x=x+m1%a(i,k)*m2%a(k,j)
           ENDDO
           m3%a(i,j)=x
        ENDDO
        ENDDO
        END FUNCTION mm_prod

        FUNCTION cm_prod(x1,m2) RESULT(m3)
        COMPLEX*16, INTENT(in) :: x1
        TYPE(matrix4), INTENT(in) :: m2
        TYPE(matrix4) :: m3
        INTEGER i,j
        DO i=1,4
        DO j=1,4
           m3%a(i,j)=x1*m2%a(i,j)
        ENDDO
        ENDDO
        END FUNCTION cm_prod

        FUNCTION mc_prod(m2,x1) RESULT(m3)
        COMPLEX*16, INTENT(in) :: x1
        TYPE(matrix4), INTENT(in) :: m2
        TYPE(matrix4) :: m3
        INTEGER i,j
        DO i=1,4
        DO j=1,4
           m3%a(i,j)=x1*m2%a(i,j)
        ENDDO
        ENDDO
        END FUNCTION mc_prod

        FUNCTION rm_prod(x1,m2) RESULT(m3)
        REAL*8, INTENT(in) :: x1
        TYPE(matrix4), INTENT(in) :: m2
        TYPE(matrix4) :: m3
        INTEGER i,j
        DO i=1,4
        DO j=1,4
           m3%a(i,j)=x1*m2%a(i,j)
        ENDDO
        ENDDO
        END FUNCTION rm_prod

        FUNCTION mr_prod(m2,x1) RESULT(m3)
        REAL*8, INTENT(in) :: x1
        TYPE(matrix4), INTENT(in) :: m2
        TYPE(matrix4) :: m3
        INTEGER i,j
        DO i=1,4
        DO j=1,4
           m3%a(i,j)=x1*m2%a(i,j)
        ENDDO
        ENDDO
        END FUNCTION mr_prod

        FUNCTION im_prod(x1,m2) RESULT(m3)
        INTEGER, INTENT(in) :: x1
        TYPE(matrix4), INTENT(in) :: m2
        TYPE(matrix4) :: m3
        INTEGER i,j
        DO i=1,4
        DO j=1,4
           m3%a(i,j)=x1*m2%a(i,j)
        ENDDO
        ENDDO
        END FUNCTION im_prod

        FUNCTION mi_prod(m2,x1) RESULT(m3)
        INTEGER, INTENT(in) :: x1
        TYPE(matrix4), INTENT(in) :: m2
        TYPE(matrix4) :: m3
        INTEGER i,j
        DO i=1,4
        DO j=1,4
           m3%a(i,j)=x1*m2%a(i,j)
        ENDDO
        ENDDO
        END FUNCTION mi_prod

        FUNCTION mr_div(m2,x1) RESULT(m3)
          REAL*8, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m2
          TYPE(matrix4) :: m3
          m3=m2*(1.d0/x1)
        END FUNCTION mr_div

        FUNCTION mc_div(m2,x1) RESULT(m3)
          COMPLEX*16, INTENT(in) :: x1
          TYPE(matrix4), INTENT(in) :: m2
          TYPE(matrix4) :: m3
          m3=m2*(1.d0/x1)
        END FUNCTION mc_div
 
       !!---- spinor operations
       !!----
        FUNCTION spsum(s1,s2) RESULT(s3)
          TYPE(spinor), INTENT(in) :: s1,s2
          TYPE(spinor) :: s3
          s3%sp=s1%sp+s2%sp
        END FUNCTION spsum
 
        FUNCTION spsub(s1,s2) RESULT(s3)
          TYPE(spinor), INTENT(in) :: s1,s2
          TYPE(spinor) :: s3
          s3%sp=s1%sp-s2%sp
        END FUNCTION spsub

        FUNCTION csp_prod(x1,s2) RESULT(s3)
        COMPLEX*16, INTENT(in) :: x1
        TYPE(spinor), INTENT(in) :: s2
        TYPE(spinor) :: s3
        INTEGER i
        DO i=1,4
           s3%sp(i)=x1*s2%sp(i)
        ENDDO
        END FUNCTION csp_prod

        FUNCTION spc_prod(s2,x1) RESULT(s3)
        COMPLEX*16, INTENT(in) :: x1
        TYPE(spinor), INTENT(in) :: s2
        TYPE(spinor) :: s3
        INTEGER i
        DO i=1,4
           s3%sp(i)=x1*s2%sp(i)
        ENDDO
        END FUNCTION spc_prod

        FUNCTION rsp_prod(x1,s2) RESULT(s3)
        REAL*8, INTENT(in) :: x1
        TYPE(spinor), INTENT(in) :: s2
        TYPE(spinor) :: s3
        INTEGER i
        DO i=1,4
           s3%sp(i)=x1*s2%sp(i)
        ENDDO
        END FUNCTION rsp_prod

        FUNCTION spr_prod(s2,x1) RESULT(s3)
        REAL*8, INTENT(in) :: x1
        TYPE(spinor), INTENT(in) :: s2
        TYPE(spinor) :: s3
        INTEGER i
        DO i=1,4
           s3%sp(i)=x1*s2%sp(i)
        ENDDO
        END FUNCTION spr_prod

        FUNCTION spmtx_prod(s2,m4) RESULT(s3)
        TYPE(matrix4), INTENT(in) :: m4
        TYPE(spinor), INTENT(in) :: s2
        TYPE(spinor) :: s3
        INTEGER i,j
        DO i=1,4
           s3%sp(i)=0.d0
           DO j=1,4
              s3%sp(i)=s3%sp(i)+s2%sp(j)*m4%a(j,i)
           ENDDO
        ENDDO
        END FUNCTION spmtx_prod 

        FUNCTION mtxsp_prod(m4,s2) RESULT(s3)
        TYPE(matrix4), INTENT(in) :: m4
        TYPE(spinor), INTENT(in) :: s2
        TYPE(spinor) :: s3
        INTEGER i,j
        DO i=1,4
           s3%sp(i)=0.d0
           DO j=1,4
              s3%sp(i)=s3%sp(i)+s2%sp(j)*m4%a(i,j)
           ENDDO
        ENDDO
        END FUNCTION mtxsp_prod 
        
        INTEGER FUNCTION lct4(L,I,J,K)
          !! 4-dim anti-symm Levi-Cevita tensor
          INTEGER, INTENT(in) :: L,I,J,K
          INTEGER IP,JP,KP,M

          lct4=0
          IF( (L-I)*(L-J)*(L-K)*(I-J)*(I-K)*(J-K).EQ.0 ) RETURN
          lct4=1-2*MOD(L,2)    ! perm. phase corresponding to a cyclic perm.
          M=4-L
          IP=MOD(I+M,4)          ! Do cyclic permutation, first index=0
          JP=MOD(J+M,4)
          KP=MOD(K+M,4)
          M=MAX(JP-IP,KP-JP)
          IF(M.NE.1) lct4=-lct4  ! M=1 Corresponds to cycl. perm. (1,2,3)
          RETURN
        END FUNCTION lct4

        SUBROUTINE ini_gamma
          INTEGER i,j,mu,nu,al,be 
          DO i=1,4
             DO j=1,4
                gone%a(I,J)=0.d0
                gzero%a(I,J)=0.d0
                kronikers(i,j)=0
                DO mu=0,5
                   gamm(mu)%a(i,j)=(0.d0,0.d0)
                ENDDO
             ENDDO
             kronikers(i,i)=1             
          ENDDO
          gamm(0)%a(1,1)=1.d0
          gamm(0)%a(2,2)=1.d0
          gamm(0)%a(3,3)=-1.d0
          gamm(0)%a(4,4)=-1.d0
          gamm(1)%a(1,4)=1.d0
          gamm(1)%a(2,3)=1.d0
          gamm(1)%a(3,2)=-1.d0
          gamm(1)%a(4,1)=-1.d0
          gamm(2)%a(1,4)=-ci
          gamm(2)%a(2,3)=ci
          gamm(2)%a(3,2)=ci
          gamm(2)%a(4,1)=-ci
          gamm(3)%a(1,3)=1.d0
          gamm(3)%a(2,4)=-1.d0
          gamm(3)%a(3,1)=-1.d0
          gamm(3)%a(4,2)=1.d0
          gamm(5)%a(1,3)=1.d0
          gamm(5)%a(2,4)=1.d0
          gamm(5)%a(3,1)=1.d0
          gamm(5)%a(4,2)=1.d0
          ! unit matrix
          gone%a(1,1)=1.d0
          gone%a(2,2)=1.d0
          gone%a(3,3)=1.d0
          gone%a(4,4)=1.d0
          ! metric
          gmunu(0)=1.d0
          gmunu(1)=-1.d0
          gmunu(2)=-1.d0
          gmunu(3)=-1.d0
          DO mu=0,3
             DO nu=0,3
                gammunu(mu,nu)=kronikers(mu+1,nu+1)*gmunu(mu)*gone
                D_(mu,nu)=kronikers(mu+1,nu+1)*gmunu(mu)
                DO al=0,3
                   DO be=0,3
                      epsilon4(mu,nu,al,be)=lct4(mu,nu,al,be)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          RETURN
        END SUBROUTINE ini_gamma

        SUBROUTINE prnmatx(m1)
        TYPE(matrix4),INTENT(in) :: m1
        INTEGER i,j
        WRITE(*,*)
        DO i=1,4
        WRITE(*,'(1x,4(f7.4,1x,f5.1,4x))') (CMPLX(m1%a(i,j)),j=1,4)
        ENDDO
        END SUBROUTINE prnmatx
        
        FUNCTION tracematx(m1) RESULT(tr)
          TYPE(matrix4),INTENT(in) :: m1
          COMPLEX*16 tr
          INTEGER i
          tr=0.d0
          DO i=1,4
             tr=tr+m1%a(i,i)
          ENDDO
        END FUNCTION tracematx

         FUNCTION u_(p,theta,phi,hel,irho) RESULT(u)
           REAL*8, INTENT(in) :: p,theta,phi,hel
           INTEGER, INTENT(in) :: irho
           REAL*8 ep,comp,RHO1,RHO2
           COMPLEX*16 ch1,ch2
           TYPE(spinor) u
           IF (hel.GT.0.d0) THEN
              ihel=1
              ch1=dcos(theta/2.d0)
              ch2=cdexp(ci*phi)*dsin(theta/2.d0)
           ELSE
              ihel=-1
              ch1=-cdexp(-ci*phi)*dsin(theta/2.d0)
              ch2=dcos(theta/2.d0)
           ENDIF
           ep=dsqrt(P*P+1.d0)
           comp=dsqrt((ep+1.d0)/(2.d0*ep))
           IF(irho.EQ.1) THEN
              RHO1=comp
              RHO2=ihel*p*comp/(ep+1.d0)
           ELSE
              RHO2=comp
              RHO1=-ihel*p*comp/(ep+1.d0)
           ENDIF
           u%sp(1)=RHO1*ch1
           u%sp(2)=RHO1*ch2
           u%sp(3)=RHO2*ch1
           u%sp(4)=RHO2*ch2
           RETURN
         END FUNCTION u_

         FUNCTION ubar_(p,theta,phi,hel,irho) RESULT(ubar)
           REAL*8, INTENT(in) :: p,theta,phi,hel
           INTEGER, INTENT(in) :: irho
           TYPE(spinor) ubar,uu
           uu=u_(p,theta,phi,hel,irho)
           ubar%sp(1) = CONJG( uu%sp(1) )
           ubar%sp(2) = CONJG( uu%sp(2) )
           ubar%sp(3) = -CONJG( uu%sp(3) )
           ubar%sp(4) = -CONJG( uu%sp(4) )
           RETURN
         END FUNCTION ubar_

         FUNCTION u0_(dum) RESULT(u)
           REAL*8, intent(in):: dum
           TYPE(spinor) u
           u%sp(1) = dum
           u%sp(2) = dum
           u%sp(3) = dum
           u%sp(4) = dum
           RETURN
         END FUNCTION u0_

         FUNCTION u_ubar_(p,theta,phi,hel,irho) RESULT(uubar)
           REAL*8, INTENT(in) :: p,theta,phi,hel
           INTEGER, INTENT(in) :: irho
           INTEGER i,j
           TYPE(matrix4) uubar
           TYPE(spinor) ubar,uu
           uu=u_(p,theta,phi,hel,irho)
           ubar=ubar_(p,theta,phi,hel,irho)

           do i=1,4
              do j=1,4
                 uubar%a(i,j)=uu%sp(i)*ubar%sp(j)
              enddo
           enddo
           RETURN
         END FUNCTION u_ubar_

         FUNCTION um_(p,theta,phi,hel,irho,rm) RESULT(u)
           REAL*8, INTENT(in) :: p,theta,phi,hel,rm
           INTEGER, INTENT(in) :: irho
           REAL*8 ep,RHO1,RHO2
           COMPLEX*16 ch1,ch2
           TYPE(spinor) u
           IF (hel.GT.0.d0) THEN
              ihel=1
              ch1=dcos(theta/2.d0)
              ch2=cdexp(ci*phi)*dsin(theta/2.d0)
           ELSE
              ihel=-1
              ch1=-cdexp(-ci*phi)*dsin(theta/2.d0)
              ch2=dcos(theta/2.d0)
           ENDIF
           ep=dsqrt(P*P+rm*rm)
           IF(irho.GT.0) THEN
              RHO1=dsqrt((ep+rm)/(2.d0*ep))
              RHO2=ihel*dsqrt((ep-rm)/(2.d0*ep))
           ELSE
              RHO2=dsqrt((ep+rm)/(2.d0*ep))
              RHO1=-ihel*dsqrt((ep-rm)/(2.d0*ep))
           ENDIF
           u%sp(1)=RHO1*ch1
           u%sp(2)=RHO1*ch2
           u%sp(3)=RHO2*ch1
           u%sp(4)=RHO2*ch2
           RETURN
         END FUNCTION um_

         FUNCTION umbar_(p,theta,phi,hel,irho,rm) RESULT(ubar)
           REAL*8, INTENT(in) :: p,theta,phi,hel,rm
           INTEGER, INTENT(in) :: irho
           TYPE(spinor) ubar,uu
           uu=um_(p,theta,phi,hel,irho,rm)
           ubar%sp(1) = CONJG( uu%sp(1) )
           ubar%sp(2) = CONJG( uu%sp(2) )
           ubar%sp(3) = -CONJG( uu%sp(3) )
           ubar%sp(4) = -CONJG( uu%sp(4) )
           RETURN
         END FUNCTION umbar_

      FUNCTION rsu_(p,theta,rm) RESULT(u)
       IMPLICIT NONE 
       REAL*8, INTENT(in) :: p,theta,rm
       TYPE(spinor) u(4,0:3)

       REAL*8 hnuc(2)
       DATA hnuc/0.5d0,-0.5d0/
       INTEGER hindx,mu
       REAL*8 sq2,sq2m1,sq3m1
       REAL*8 X,Z,h,PI
       complex*16 eps(0:3,3)
       TYPE(spinor) uu,uu1
       sq2=dsqrt(2.d0)
       sq2m1=1.d0/sq2
       sq3m1=1.d0/dsqrt(3.d0)
       X= dsin(PI-theta)
       Z= dcos(PI-theta)
       eps(0,1)= 0.d0   ! Theta1=Theta_gamma + 90, Phi1=Phi_gamma=0.
       eps(1,1)= -sq2m1*Z
       eps(2,1)= -ci*sq2m1
       eps(3,1)= sq2m1*X
       eps(0,2)= 0.d0 
       eps(1,2)= sq2m1*Z
       eps(2,2)= -ci*sq2m1
       eps(3,2)= -sq2m1*X
       eps(0,3)= p/rm
       !eps(0,3)= p/dsqrt(rm**2+p**2)
       eps(1,3)= 0.d0
       eps(2,3)= 0.d0
       eps(3,3)= -dsqrt(rm**2+p**2)/rm
       !eps(3,3)=-1.d0
       
       !write(*,*) 'p = ',p
       DO hindx=1,2
          h=hnuc(hindx)
          uu=um_(p,theta,0.d0,h,1,rm)
          uu1=um_(p,theta,0.d0,-h,1,rm)
          DO mu=0,3
             u(hindx,mu)=eps(mu,hindx)*uu
             u(hindx+2,mu)=(eps(mu,hindx)*uu1 &
                  +  sq2*eps(mu,3)*uu)*sq3m1
          ENDDO
       enddo
       RETURN
   END FUNCTION rsu_

     FUNCTION rsubar_(p,theta,rm) RESULT(ubar)
       IMPLICIT NONE 
       REAL*8, INTENT(in) :: p,theta,rm
       TYPE(spinor) ubar(4,0:3)
       INTEGER hindx,mu
       TYPE(spinor) uu(4,0:3)
           uu=rsu_(p,theta,rm)
           do mu=0,3
              DO hindx=1,4
                 ubar(hindx,mu)%sp(1) = CONJG( uu(hindx,mu)%sp(1) )
                 ubar(hindx,mu)%sp(2) = CONJG( uu(hindx,mu)%sp(2) )
                 ubar(hindx,mu)%sp(3) = -CONJG( uu(hindx,mu)%sp(3) )
                 ubar(hindx,mu)%sp(4) = -CONJG( uu(hindx,mu)%sp(4) )
              ENDDO
           enddo
       RETURN
     END FUNCTION rsubar_


         FUNCTION rsu_rsubar_(p,theta,rm) RESULT(rsuubar)
           REAL*8, INTENT(in) :: p,theta,rm
           INTEGER i,j,k,mu,nu
           TYPE(matrix4) rsuubar(4,0:3,0:3)
           TYPE(spinor) ubar(4,0:3),uu(4,0:3)
           uu=rsu_(p,theta,rm)
           ubar=rsubar_(p,theta,rm)

           do 21 mu=0,3
           do 21 nu=0,3
           do 21 k=1,4
           do i=1,4
              do j=1,4
                 rsuubar(k,mu,nu)%a(i,j)=uu(k,mu)%sp(i)*ubar(k,nu)%sp(j)
              enddo
           enddo
21       continue
           RETURN
         END FUNCTION rsu_rsubar_

        SUBROUTINE prnspinor(u)
        TYPE(spinor),INTENT(in) :: u
        INTEGER i
        WRITE(*,*)
        WRITE(*,*) (CMPLX(u%sp(i)),i=1,4)
        END SUBROUTINE prnspinor

         SUBROUTINE ini_polar_real(theta)
           REAL*8, INTENT(in) :: theta
           REAL*8 X,Y,Z,sq2
           INTEGER i,mu
           X= dsin(theta)
           Y=0.d0
           Z= dcos(theta)
           DO i=1,3
              DO mu=0,3
                 vEPS(mu,i) = (0.d0,0.d0)
                 vEPSP(mu,i) =(0.d0,0.d0)
              ENDDO
           ENDDO
           sq2=1.d0/dsqrt(2.d0)
           vEPS(1,1)= -sq2
           vEPS(2,1)= -ci*sq2
           vEPS(1,2)= SQ2
           vEPS(2,2)= -ci*sq2
           veps(3,3)=1.d0
           vEPSP(0,1)= 0.d0   ! Theta1=Theta_gamma + 90, Phi1=Phi_gamma=0.
           vEPSP(1,1)= -sq2*Z
           vEPSP(2,1)= ci*sq2
           vEPSP(3,1)= sq2*X
           vEPSP(1,2)= sq2*Z
           vEPSP(2,2)= ci*sq2
           vEPSP(3,2)= -sq2*X
           RETURN
         END SUBROUTINE ini_polar_real

         SUBROUTINE ini_polar_virtual(theta,Qsq,om)
           REAL*8, INTENT(in) :: theta,Qsq,om
           REAL*8 X,Y,Z,sq2,dQsq,qvec
           INTEGER i,mu
           if(Qsq.eq.0.d0) then
           dQsq=om
           qvec=om
           else
           dQsq=dsqrt(Qsq)
           qvec=dsqrt(Qsq+om**2)
           endif
           X= dsin(theta)
           Y=0.d0
           Z= dcos(theta)
           DO i=1,3
              DO mu=0,3
                 vEPS(mu,i) = (0.d0,0.d0)
                 vEPSP(mu,i) =(0.d0,0.d0)
              ENDDO
           ENDDO
           sq2=1.d0/dsqrt(2.d0)
           vEPS(1,1)= -sq2
           vEPS(2,1)= -ci*sq2
           vEPS(1,2)= SQ2
           vEPS(2,2)= -ci*sq2
           veps(0,3)=1.d0 
         !  veps(0,3)=qvec/dQsq
           veps(3,3)=-1.d0 
         !  veps(3,3)=-om/dQsq
         !  vEPSP(0,1)= 0.d0   ! Theta1=Theta_gamma + 90, Phi1=Phi_gamma=0.
         !  vEPSP(1,1)= -sq2*Z
         !  vEPSP(2,1)= ci*sq2
         !  vEPSP(3,1)= sq2*X
         !  vEPSP(1,2)= sq2*Z
         !  vEPSP(2,2)= ci*sq2
         !  vEPSP(3,2)= -sq2*X
           RETURN
         END SUBROUTINE ini_polar_virtual


         FUNCTION rslash(p) RESULT(slashp)
           REAL*8, INTENT(in) :: p(0:3)
           TYPE(matrix4) slashp
           slashp=p(0)*gamm(0)-p(1)*gamm(1)-p(2)*gamm(2)-p(3)*gamm(3)
           RETURN
         END FUNCTION rslash

         FUNCTION cslash(p) RESULT(slashp)
           COMPLEX*16, INTENT(in) :: p(0:3)
           TYPE(matrix4) slashp
           slashp=p(0)*gamm(0)-p(1)*gamm(1)-p(2)*gamm(2)-p(3)*gamm(3)
           RETURN
         END FUNCTION cslash
  

        FUNCTION braket(ub,g,u) RESULT(s)
          TYPE(spinor), INTENT(in) :: ub,u
          TYPE(matrix4), INTENT(in) :: g
          INTEGER i,j
          COMPLEX*16 s,x
          s=0.d0
          DO i=1,4
             x=0.d0
             DO j=1,4
                x=x+g%a(i,j)*u%sp(j)
             ENDDO
             s=s+ub%sp(i)*x
          ENDDO
          RETURN
        END FUNCTION braket

        FUNCTION contract0(verout,tprop,verin) RESULT(ample)
          IMPLICIT NONE
          TYPE(matrix4), INTENT(in) :: verin(0:3),verout(0:3),tprop(0:3,0:3)
          TYPE(matrix4) ample
          INTEGER mu,nu
          ample=gzero
          DO mu=0,3
             DO nu=0,3
                ample=ample+gmunu(mu)*gmunu(nu)*verout(mu) &
                     *tprop(mu,nu)*verin(nu)
             ENDDO
          ENDDO
          RETURN
        END FUNCTION contract0

        FUNCTION contract1in(verout,tprop,verin) RESULT(ample)
          IMPLICIT NONE
          TYPE(matrix4), INTENT(in) :: verin(0:3,0:3),verout(0:3), &
               tprop(0:3,0:3)
          TYPE(matrix4) ample(0:3)
          INTEGER mu,nu,al

          DO al=0,3
             ample(al)=gzero
             DO mu=0,3
                DO nu=0,3
                   ample(al)=ample(al)+gmunu(mu)*gmunu(nu)*verout(mu) &
                        *tprop(mu,nu)*verin(nu,al)
                ENDDO
             ENDDO
          ENDDO
          RETURN
        END FUNCTION contract1in

        FUNCTION contract1out(verout,tprop,verin) RESULT(ample)
          IMPLICIT NONE
          TYPE(matrix4), INTENT(in) :: verin(0:3),verout(0:3,0:3), &
               tprop(0:3,0:3)
          TYPE(matrix4) ample(0:3)
          INTEGER mu,nu,al

          DO al=0,3
             ample(al)=gzero
             DO mu=0,3
                DO nu=0,3
                   ample(al)=ample(al)+gmunu(mu)*gmunu(nu)*verout(mu,al) &
                        *tprop(mu,nu)*verin(nu)
                ENDDO
             ENDDO
          ENDDO
          RETURN
        END FUNCTION contract1out

        FUNCTION contract2(verout,tprop,verin) RESULT(ample)
          IMPLICIT NONE
          TYPE(matrix4), INTENT(in) :: verin(0:3,0:3),verout(0:3,0:3), &
               tprop(0:3,0:3)
          TYPE(matrix4) ample(0:3,0:3)
          INTEGER mu,nu,al,be

          DO al=0,3
             DO be=0,3
             ample(al,be)=gzero
             DO mu=0,3
                DO nu=0,3
                   ample(al,be)=ample(al,be)+gmunu(mu)*gmunu(nu) &
                        *verout(mu,al)*tprop(mu,nu)*verin(nu,be)
                ENDDO
             ENDDO
             ENDDO
          ENDDO
          RETURN
        END FUNCTION contract2

        FUNCTION rrdot(p,q)
          REAL*8, INTENT(in) :: p(0:3),q(0:3)
          REAL*8 rrdot
          rrdot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
        END FUNCTION rrdot

        FUNCTION rcdot(p,q)
          REAL*8, INTENT(in) :: p(0:3)
          COMPLEX*16, INTENT(in) :: q(0:3)
          COMPLEX*16 rcdot
          rcdot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
        END FUNCTION rcdot

        FUNCTION crdot(p,q)
          REAL*8, INTENT(in) :: q(0:3)
          COMPLEX*16, INTENT(in) :: p(0:3)
          COMPLEX*16 crdot
          crdot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
        END FUNCTION crdot

        FUNCTION ccdot(p,q)
          COMPLEX*16, INTENT(in) :: p(0:3),q(0:3)
          COMPLEX*16 ccdot
          ccdot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
        END FUNCTION ccdot

        FUNCTION crag2(p,q)
          REAL*8, INTENT(in) :: q(0:3)
          COMPLEX*16, INTENT(in) :: p(0:3)
          TYPE(matrix4) crag2
          crag2=gzero
          DO mu=0,3
          DO nu=0,3
             crag2=crag2+gmunu(mu)*gmunu(nu) &
                        *(gamm(mu)*gamm(nu)-gammunu(mu,nu))*p(mu)*q(nu)
          ENDDO
          ENDDO
        END FUNCTION crag2

        FUNCTION rcag2(p,q)
          REAL*8, INTENT(in) :: p(0:3)
          COMPLEX*16, INTENT(in) :: q(0:3)
          TYPE(matrix4) rcag2
          rcag2=gzero
          DO mu=0,3
          DO nu=0,3
             rcag2=rcag2+gmunu(mu)*gmunu(nu) &
                        *(gamm(mu)*gamm(nu)-gammunu(mu,nu))*p(mu)*q(nu)
          ENDDO
          ENDDO
        END FUNCTION rcag2

        FUNCTION rrag2(p,q)
          REAL*8, INTENT(in) :: p(0:3), q(0:3)
          TYPE(matrix4) rrag2
          rrag2=gzero
          DO mu=0,3
          DO nu=0,3
             rrag2=rrag2+gmunu(mu)*gmunu(nu) &
                        *(gamm(mu)*gamm(nu)-gammunu(mu,nu))*p(mu)*q(nu)
          ENDDO
          ENDDO
        END FUNCTION rrag2

       FUNCTION ccag2(p,q)
          COMPLEX*16, INTENT(in) :: p(0:3),q(0:3)
          TYPE(matrix4) ccag2
          ccag2=gzero
          DO mu=0,3
          DO nu=0,3
             ccag2=ccag2+gmunu(mu)*gmunu(nu) &
                        *(gamm(mu)*gamm(nu)-D_(mu,nu)*gone)*p(mu)*q(nu)
          ENDDO
          ENDDO
        END FUNCTION ccag2

        FUNCTION rrrreps4(v1,v2,v3,v4)
          REAL*8, INTENT(in) :: v1(0:3),v2(0:3),v3(0:3),v4(0:3)
          REAL*8 rrrreps4
          INTEGER mu,nu,al,be
          rrrreps4=0.d0
          DO mu=0,3
          DO nu=0,3
             DO al=0,3
                DO be=0,3
                   rrrreps4=rrrreps4+gmunu(mu)*gmunu(nu)*gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(mu)*v2(nu)*v3(al)*v4(be)
                ENDDO
             ENDDO
          ENDDO
          ENDDO
        end function rrrreps4

       FUNCTION rrrceps4(v1,v2,v3,v4)
          REAL*8, INTENT(in) :: v1(0:3),v2(0:3),v3(0:3)
          Complex*16, INTENT(in) :: v4(0:3)
          complex*16 rrrceps4
          INTEGER mu,nu,al,be
          rrrceps4=0.d0
          DO mu=0,3
          DO nu=0,3
             DO al=0,3
                DO be=0,3
                   rrrceps4=rrrceps4+gmunu(mu)*gmunu(nu)*gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(mu)*v2(nu)*v3(al)*v4(be)
                ENDDO
             ENDDO
          ENDDO
          ENDDO
        end function rrrceps4

       FUNCTION rrcceps4(v1,v2,v3,v4)
          REAL*8, INTENT(in) :: v1(0:3),v2(0:3)
          Complex*16, INTENT(in) :: v4(0:3),v3(0:3)
          complex*16 rrcceps4
          INTEGER mu,nu,al,be
          rrcceps4=0.d0
          DO mu=0,3
          DO nu=0,3
             DO al=0,3
                DO be=0,3
                   rrcceps4=rrcceps4+gmunu(mu)*gmunu(nu)*gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(mu)*v2(nu)*v3(al)*v4(be)
                ENDDO
             ENDDO
          ENDDO
          ENDDO
        end function rrcceps4

       FUNCTION rcrceps4(v1,v2,v3,v4)
          REAL*8, INTENT(in) :: v1(0:3),v3(0:3)
          Complex*16, INTENT(in) :: v4(0:3),v2(0:3)
          complex*16 rcrceps4
          INTEGER mu,nu,al,be
          rcrceps4=0.d0
          DO mu=0,3
          DO nu=0,3
             DO al=0,3
                DO be=0,3
                   rcrceps4=rcrceps4+gmunu(mu)*gmunu(nu)*gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(mu)*v2(nu)*v3(al)*v4(be)
                ENDDO
             ENDDO
          ENDDO
          ENDDO
        end function rcrceps4

      FUNCTION crcreps4(v1,v2,v3,v4)
          REAL*8, INTENT(in) :: v2(0:3),v4(0:3)
          Complex*16, INTENT(in) :: v1(0:3),v3(0:3)
          COMPLEX*16 crcreps4
          INTEGER mu,nu,al,be
          crcreps4=0.d0
          DO mu=0,3
          DO nu=0,3
             DO al=0,3
                DO be=0,3
                   crcreps4=crcreps4+gmunu(mu)*gmunu(nu)*gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(mu)*v2(nu)*v3(al)*v4(be)
                ENDDO
             ENDDO
          ENDDO
          ENDDO
        end function crcreps4

       FUNCTION crrceps4(v1,v2,v3,v4)
          REAL*8, INTENT(in) :: v3(0:3),v2(0:3)
          Complex*16, INTENT(in) :: v4(0:3),v1(0:3)
          complex*16 crrceps4
          INTEGER mu,nu,al,be
          crrceps4=0.d0
          DO mu=0,3
          DO nu=0,3
             DO al=0,3
                DO be=0,3
                   crrceps4=crrceps4+gmunu(mu)*gmunu(nu)*gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(mu)*v2(nu)*v3(al)*v4(be)
                ENDDO
             ENDDO
          ENDDO
          ENDDO
        end function crrceps4

      FUNCTION ccrreps4(v1,v2,v3,v4)
          REAL*8, INTENT(in) :: v3(0:3),v4(0:3)
          Complex*16, INTENT(in) :: v2(0:3),v1(0:3)
          complex*16 ccrreps4
          INTEGER mu,nu,al,be
          ccrreps4=0.d0
          DO mu=0,3
          DO nu=0,3
             DO al=0,3
                DO be=0,3
                   ccrreps4=ccrreps4+gmunu(mu)*gmunu(nu)*gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(mu)*v2(nu)*v3(al)*v4(be)
                ENDDO
             ENDDO
          ENDDO
          ENDDO
        end function ccrreps4

     FUNCTION rccreps4(v1,v2,v3,v4)
          REAL*8, INTENT(in) :: v1(0:3),v4(0:3)
          Complex*16, INTENT(in) :: v2(0:3),v3(0:3)
          complex*16 rccreps4
          INTEGER mu,nu,al,be
          rccreps4=0.d0
          DO mu=0,3
          DO nu=0,3
             DO al=0,3
                DO be=0,3
                   rccreps4=rccreps4+gmunu(mu)*gmunu(nu)*gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(mu)*v2(nu)*v3(al)*v4(be)
                ENDDO
             ENDDO
          ENDDO
          ENDDO
        end function rccreps4
      END MODULE dirac

 
