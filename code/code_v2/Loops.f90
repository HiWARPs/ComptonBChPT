 SUBROUTINE OneLoop(omega,t,gafac,mpi,A,xp,dxp,nop)
   USE gauss
   IMPLICIT NONE
   INTEGER, INTENT(in) :: nop
   REAL*8, INTENT(in) :: t,gafac,mpi,omega,xp(nop),dxp(nop)
   COMPLEX*16,INTENT(inout) :: A(8)
   COMPLEX*16, PARAMETER :: CI=(0.d0,1.d0)
   REAL*8 sqt,mpi2,PI,omega2,costh,Rxz,z0,dzeta
   REAL*8 aa,aa2,bb,bb2, albe,beta,inv_nmass_2fm
   COMPLEX*16 a1,a2,a3,a4,a5,a6
   COMPLEX*16 sumAz1,sumAz2,sumAz3,sumAz4,sumAz5,Wxz,arcsinxz
   INTEGER i,j,np1,np2
   REAL*8, ALLOCATABLE :: z1(:),dz1(:),z2(:),dz2(:)

   PI=4.d0*datan(1.d0)
   mpi2=mpi*mpi
   omega2=omega*omega
   costh=t/(2.d0*omega2)+1.d0
   !write(*,*) ' costh =', costh, -t-4.*omega2
   sqt=dsqrt(-t)

   sumAz1=0.d0;sumAz2=0.d0;sumAz3=0.d0;sumAz4=0.d0;sumAz5=0.d0
   IF(omega < mpi) THEN ! below the threshold
      DO i=1,nop
         sumAz1 = sumAz1 + dxp(i)*Az1(xp(i))
         sumAz2 = sumAz2 + dxp(i)*Az2(xp(i))
         CALL integrate(xp,dxp,nop)
      ENDDO
   ELSE ! above threshold
      aa=sqt/(2.d0*omega)
      bb=mpi/omega
      call getAZ12(aa,bb,sumAz1,sumAz2)
      call getAZ345(sumAZ3,sumAZ4,sumAZ5)
   ENDIF
      aa=sqt/(2.d0*omega)
      bb=mpi/omega
      !call getAZ12(aa,bb,sumAz1,sumAz2)
      WRITE(8,*) 1./bb, sqt/mpi
      WRITE(8,*) 'sumAZ1 =',sumAZ1
      write(8,*) 'sumAZ2 =',sumAZ2*sqt
      !WRITE(8,*) 'sumAZ3=',sumAZ3
      !write(8,*) 'sumAZ4 =',sumAZ4
      !WRITE(8,*) 'sumAZ5=',sumAZ5
   IF(sqt > 0.d0) THEN
      a1 = gafac*(mpi - Esqrt(mpi2 - omega2) &
           + (2.0d0*mpi2 - t)/sqt &
           * (datan(sqt/(2.d0*mpi))/2.0d0 - sumAZ1))
      a2 = gafac*(2.0d0*mpi2/t-1.d0)*(sumAZ1/sqt - sumAZ2) ! "+"is CORRECT
   ELSE ! forward scattering only
      a1 = gafac*( (3./2.d0)*mpi - Esqrt(mpi2 - omega2) &
           - (mpi2/omega)*( ArcSin(omega/mpi) &
           - (mpi - Esqrt(mpi2 - omega2))/omega ))
      a2=-1.d0
   ENDIF
   
   a3 = gafac/pi*(mpi2/omega*(ArcSin(omega/mpi)**2) - omega) &
        + 2.0d0*gafac/pi*(omega**4)*(1.d0-costh**2)*sumAZ3
   a4 = 2.0d0*gafac/pi*sumAZ4
   a5 = gafac/pi*(-sumAZ5 + 2.d0*omega2*costh*sumAZ3)
   a6 = gafac/pi*(sumAZ5 - 2.d0*omega2*sumAZ3)
   
   A(1)=A(1)+a1
   A(2)=A(2)+a2
   A(3)=A(3)+a3      
   A(4)=A(4)+a4      
   A(5)=A(5)+a5
   A(6)=A(6)+a6
   

                albe=dreal(a1+0.5d0*t*a2)
                beta=-dreal(a2)
                write(11,*) omega, t
                inv_nmass_2fm=(1.d0/0.9383)*0.19733
                write(11,*) 'alpha+beta= ', (1./137.)*albe/(omega**2)*inv_nmass_2fm**3*1.d4
                write(11,*) 'beta = ', (1./137.)*beta*inv_nmass_2fm**3*1.d4 
	write(11,*)
   
 CONTAINS
   
   SUBROUTINE integrate(z,dz,n)
     IMPLICIT NONE
     INTEGER, INTENT(in) :: n
     INTEGER k
     REAL*8, INTENT(in) :: z(n),dz(n)


     DO k=1,n
        Wxz=W(xp(i),z(k))
        Rxz=R(xp(i),z(k))
        arcsinxz = ArcSin(omega*z(k)/Rxz)
        sumAZ3 = sumAZ3 + dxp(i)*dz(k)*AZ3(xp(i),z(k))
        sumAZ4 = sumAZ4 + dxp(i)*dz(k)*AZ4(xp(i),z(k))
        sumAZ5 = sumAZ5 + dxp(i)*dz(k)*AZ5(xp(i),z(k))
        !write(*,*) 'W, R= ',k,Wxz,Rxz
        !write(*,*) 'A3..A5= ', &
             !     AZ3(xp(i),z(k)),AZ4(xp(i),z(k)),AZ5(xp(i),z(k))
     END DO
   END SUBROUTINE integrate



   COMPLEX*16 FUNCTION W(x,z)
     REAL*8 x,z,arg
     arg = mpi2 - (omega*z)**2 + t*(1.0d0 - z)**2*x*(x - 1.0d0)
     W = Esqrt(arg)
   END FUNCTION W

   FUNCTION R(x,z)
     REAL*8 R,x,z
     R= dsqrt(mpi2 + t*(1.0d0-z)**2*x*(x - 1.0d0))
   END FUNCTION R

   COMPLEX*16 FUNCTION AZ1(z)
     REAL*8 z,arg,y
     REAL*8, PARAMETER ::  EPS=1.e-5
     COMPLEX*16, PARAMETER :: ci=(0.d0,1.d0)
     arg = mpi2 - (omega*z)**2
     IF(arg >= 0.d0) THEN
        AZ1=datan((1.0d0 - z)*sqt/(2.0d0*dsqrt(arg)))
     ELSE
        y=(1.0d0 - z)*sqt/(2.0d0*dsqrt(-arg))
        AZ1=-0.5d0*ci*dlog(  (1.d0+y)/(1.d0-y)  )
     ENDIF
     !write(*,*) arg,y,AZ1
     !stop
   END FUNCTION AZ1

   COMPLEX*16 FUNCTION AZ2(z)
     REAL*8 z,arg
     arg = mpi2 - (omega*z)**2
    AZ2 = 2.0d0*(1.d0-z)*Esqrt(arg) &
          /(4.0d0*arg - t*(1.d0 - z)**2)
   END FUNCTION AZ2

   COMPLEX*16 FUNCTION AZ3(x,z)
     REAL*8, INTENT(in) :: x,z
     AZ3= x*(1.0d0-x)*z*(1.d0-z)**3/(Wxz*Wxz*Wxz) &
          *(arcsinxz + omega*z*Wxz/(Rxz**2))
   END FUNCTION AZ3

   SUBROUTINE getAZ345(sAZ3,sAZ4,sAZ5)
     Complex*16, INTENT(inout) :: sAZ3,sAZ4,sAZ5
     REAL*8 z,dzeta,x,z0,z1
     integer i,j

     sAZ3=0.d0;sAZ4=0.d0;sAZ5=0.d0
      DO i=1,np
         ! determine singularity in W(x,z): z0
         dzeta=t*xp(i)*(xp(i)-1.d0)/omega2 ! defined on [0,1]
         IF(dzeta < 1.d0) THEN
            z0 = (dsqrt(dzeta+(1.d0-dzeta)*mpi2/omega2)-dzeta) &
                 /(1.d0-dzeta)
         ELSE ! dzeta =1
            z0 = 0.5d0*(1.d0+mpi2/omega2)
         ENDIF
         z1=(dsqrt(dzeta+(1.d0-dzeta)*mpi2/omega2)+dzeta) &
                 /(1.d0-dzeta)
         x=xp(i)
         !write(*,*)dzeta,z0,z1
         IF(z0 > 0.d0.AND.z0 <1.d0) THEN 
            DO j=1,np
               z=xp(j)
               sAZ3=sAZ3+ (fAZ3a(z,x)-fAZ3a(z0,x))/((z-z0)*(z+z1)) &
                    *dxp(i)*dxp(j)/Esqrt((z-z0)*(z+z1)) 

               sAZ3=sAZ3+(fAZ3b(z,x)-fAZ3b(z0,x))/((z-z0)*(z+z1))*dxp(i)*dxp(j)
               
               sAZ4=sAZ4+ (fAZ4(z,x)-fAZ4(z0,x)) &
                    *dxp(i)*dxp(j)/Esqrt((z-z0)*(z+z1)) 
               sAZ5=sAZ5+ (fAZ5(z,x)-fAZ5(z0,x)) &
                    *dxp(i)*dxp(j)/Esqrt((z-z0)*(z+z1)) 
            enddo
            sAZ3=sAZ3+ (fAZ3a(z0,x)*(intAZ3a(1.d0,z0,z1)-intAZ3a(0.d0,z0,z1)) &
                 + fAZ3b(z0,x)/(z0+z1) &
                 *cdlog((1.d0,0.d0)*((1.d0-z0)*z1)/(-z0*(1.d0+z1))) )*dxp(i)
            sAZ4=sAZ4 + fAZ4(z0,x)*dxp(i) &
                 *(intAZ4(1.d0,z0,z1)-intAZ4(0.d0,z0,z1))
            sAZ5=sAZ5 + fAZ5(z0,x)*dxp(i) &
                 *(intAZ4(1.d0,z0,z1)-intAZ4(0.d0,z0,z1))
         else
            DO j=1,np
               Wxz=W(xp(i),xp(j))
               Rxz=R(xp(i),xp(j))
               arcsinxz = ArcSin(omega*xp(j)/Rxz)
               sAZ3 = sAZ3 + dxp(i)*dxp(j)*AZ3(xp(i),xp(j))
               sAZ4 = sAZ4 + dxp(i)*dxp(j)*AZ4(xp(i),xp(j))
               sAZ5 = sAZ5 + dxp(i)*dxp(j)*AZ5(xp(i),xp(j))
           enddo
         endif
      enddo
   END SUBROUTINE getAZ345

   COMPLEX*16 FUNCTION intAZ3a(z,z0,z1)
     REAL*8, INTENT(in) :: z,z0,z1
     intAZ3a = -2.d0*(2.*z-z0+z1)/(z1+z0)**2 &
          *1.d0/Esqrt((z-z0)*(z+z1))
   end function intAZ3a

   COMPLEX*16 FUNCTION fAZ3a(z,x)
     REAL*8, INTENT(in) :: z,x
     arcsinxz = ArcSin(omega*z/R(x,z))
     fAZ3a= x*(1.0d0-x)*z*(1.d0-z)**3*arcsinxz
   end function fAZ3a


   COMPLEX*16 FUNCTION fAZ3b(z,x)
     REAL*8, INTENT(in) :: z,x
     fAZ3b= x*(1.0d0-x)*z*(1.d0-z)**3*omega*z/(R(x,z)**2)
   end function fAZ3b

   COMPLEX*16 FUNCTION intAZ4(z,z0,z1)
     REAL*8, INTENT(in) :: z,z0,z1
     intAZ4 = cdlog(z+(-z0+z1)/2.d0 + Esqrt((z-z0)*(z+z1)))
   end function intAZ4

   COMPLEX*16 FUNCTION fAZ4(x,z) 
     REAL*8, INTENT(in) :: x,z
     arcsinxz = ArcSin(omega*z/R(x,z))
     fAZ4= z*(1.d0-z)*arcsinxz
   END FUNCTION fAZ4

   COMPLEX*16 FUNCTION fAZ5(x,z) 
     REAL*8, INTENT(in) :: x,z
     arcsinxz = ArcSin(omega*z/R(x,z))
     fAZ5= (1.d0-z)**2*arcsinxz
   END FUNCTION fAZ5

   COMPLEX*16 FUNCTION AZ4(x,z) 
     REAL*8, INTENT(in) :: x,z
     AZ4= z*(1.d0-z)/Wxz*arcsinxz
   END FUNCTION AZ4

   COMPLEX*16 FUNCTION AZ5(x,z) 
     REAL*8, INTENT(in) :: x,z
     AZ5= (1.0d0-z)**2/Wxz*arcsinxz
   END FUNCTION AZ5

   COMPLEX*16 FUNCTION ArcTan(z)
     COMPLEX*16, INTENT(in) :: z
     COMPLEX*16, PARAMETER ::  ci=(0.d0,1.d0)
     ArcTan= -0.5d0*ci*cdlog( (1.d0+ci*z)/(1.d0-ci*z) )
     RETURN
   END FUNCTION ArcTan

   COMPLEX*16 FUNCTION ArcSin(z)
     REAL*8, INTENT(in) :: z
     COMPLEX*16, PARAMETER ::  ci=(0.d0,1.d0)
     IF(dabs(z).LE.1.d0) THEN
        ArcSin=dasin(z)
     ELSE
        ArcSin= -ci*cdlog(ci*(z + dsqrt(z*z-1.d0)) )
     ENDIF
   END FUNCTION ArcSin

   COMPLEX*16 FUNCTION Esqrt(z)
     REAL*8, INTENT(in) :: z
     COMPLEX*16, PARAMETER ::  ci=(0.d0,1.d0)
     IF(z.GE.0.d0) THEN
        Esqrt=dsqrt(z)
     ELSE
        Esqrt= ci*dsqrt(-z)
     ENDIF
   END FUNCTION Esqrt


   Subroutine getAZ12(a,b,sAZ1,sAZ2)
     REAL*8, intent(in) :: a,b
     real*8 a2,b2,z1,z2
     complex*16, intent(out) :: sAZ1,sAZ2
     write(*,*)
     write(*,*) 'a,b = ', a,b
     a2=a*a; b2=b*b
     if(a.ne.1.d0) then
        z1=-(a2+dsqrt(a2+b2-a2*b2))/(1.d0-a2)
        z2=-(a2-dsqrt(a2+b2-a2*b2))/(1.d0-a2)
        sAZ1= a/(1.d0-a2)*(NintAz1(z1,b)-NintAz1(z2,b))/(z1-z2)
        !write(*,*) 'AZ1(1) = ', AZ1(1.d0)
        !write(*,*) 'z1 = ', z1
        !write(*,*) intAz1(z1,b),NintAz1(z1,b),NNintAz1(z1,b)
        !write(*,*) intAz2(z1,b) ,NintAz2(z1,b)
        !write(*,*) 'z2 = ', z2
        !write(*,*) intAz1(z2,b),NintAz1(z2,b),NNintAz1(z2,b)
        !write(*,*) intAz2(z2,b),NintAz2(z2,b)
        
           sAZ2=-1.d0/(1.d0-a2)*(NintAz2(z1,b)-NintAz2(z2,b))/(z1-z2)

     else
        z1=(b2+1.d0)/2.d0
        sAZ2=-0.5d0*(NintAz2(z1,b))
        sAZ1= 0.5d0*NintAz1(z1,b)
        !Write(*,*) 'a=1 ', a
     endif
      sAz2=sAz2/(2.d0*omega)
   END subroutine getAZ12

   COMPLEX*16 FUNCTION intAz2(x0,b)
     REAL*8, intent(in):: x0,b
     complex*16 sum
     sum=0.d0

     do i=1,np
        sum = sum + dxp(i)*( - Esqrt(b**2-xp(i)**2) + (1.d0-x0) &
             * (Esqrt(b**2-xp(i)**2)-Esqrt(b**2-x0**2))/(xp(i)-x0))
     ENDDO
     sum=sum+Esqrt(b**2-x0**2)*dlog( dabs((1.d0-x0)/(-x0)) )*(1.d0-x0)
     

     intAz2 = sum
     !write(*,*) x,x0
     !write(*,*) intAz2
   end function intAz2

   FUNCTION NNintAz1(x0,b) RESULT(res)
     real*8 x0,b
     COMPLEX*16 res
     complex*16 sum
     sum=0.d0

     do i=1,np
        sum = sum + dxp(i)*(rfun(xp(i),b)-rfun(x0,b))/(xp(i)-x0)
     ENDDO
     sum=sum+rfun(x0,b)*dlog( dabs( (1.d0-x0)/(-x0)) )
     res=sum

   end function NNintAz1
   COMPLEX*16 FUNCTION rfun(z,b)
     REAL*8 z,b
     rfun = b**2-z**2
   END FUNCTION rfun

   FUNCTION NintAz2(x0,b) RESULT(res)
     real*8 x0,b
     COMPLEX*16 res
     complex*16 sum
     sum=0.d0

     do i=1,np
        sum = sum + dxp(i)*(ffun(xp(i),b)-ffun(x0,b))/(xp(i)-x0)
     ENDDO
     sum=sum+ffun(x0,b)*dlog( dabs((1.d0-x0)/(-x0)) )
     res=sum

   end function NintAz2



   COMPLEX*16 FUNCTION ffun(z,b)
     REAL*8 z,b
     ffun = (1.d0-z)*Esqrt(b**2-z**2)
   END FUNCTION ffun

   FUNCTION NintAz1(x0,b) RESULT(res)
     real*8, intent(in) :: x0,b
     real*8, parameter :: eps=0.01
     COMPLEX*16 res
     complex*16 sum, test
     sum=0.d0

     !if(dabs(b-x0).lt.eps) then
     if(x0.gt.0.d0.and.x0.lt.1.d0) then
        write(*,*)' b-x0 = ', b-x0
        do i=1,np
           sum = sum + dxp(i)*(gfun2(xp(i),b)-gfun2(x0,b))/(xp(i)-x0) 
        ENDDO
        sum=sum+gfun2(x0,b)*dlog( dabs( (1.d0-x0)/(-x0)) )
        write(*,*) 'sum = ', sum
        sum=sum - b*(b-1.d0-Esqrt(b**2-1.d0))
     else
        do i=1,np
           sum = sum + dxp(i)*(gfun(xp(i),x0,b)-gfun(x0,x0,b))/(xp(i)-x0)
           !&
           !     +b**2*Esqrt((b-xp(i))/(b+xp(i))))
        ENDDO
        sum=sum+gfun(x0,x0,b)*dlog( dabs( (1.d0-x0)/(-x0)) )
        sum=sum*(1.d0-x0)
        sum=sum+b**2*(1.d0-b)*ArcSin(1.d0/b)
        sum=sum/(b-x0)


        !if(x0.gt.0.d0.and.x0.lt.1.d0) sum=0.d0
        do i=1,np
           sum = sum - dxp(i) &
                *(Esqrt(b**2-xp(i)**2)-Esqrt(b**2-x0**2))/(xp(i)-x0)
        ENDDO
        sum=sum-Esqrt(b**2-x0**2)*dlog( dabs((1.d0-x0)/(-x0)) )

     endif
 

     res=(1.d0,0.d0)*dreal(sum) + CI*dimag(intAz1(x0,b))
    ! res=sum

   end function NintAz1

    COMPLEX*16 FUNCTION gfun(z,x0,b)
     REAL*8 z,b,x0
     gfun = b**2*Esqrt((b-z)/(b+z))
   END FUNCTION gfun

    COMPLEX*16 FUNCTION gfun2(z,b)
     REAL*8 z,b
     gfun2 = -Esqrt(b**2-z**2)+b**2*Esqrt((b-z)/(b+z))
   END FUNCTION gfun2

   COMPLEX*16 FUNCTION intAz1(x0,b)
     REAL*8, intent(in):: x0,b
     intAz1 = -Esqrt(b**2-1.d0)*(1.d0-x0)/dabs(1.d0-x0)- b*x0/dabs(x0) &
          +(b**2-x0)*(-ArcSin(1.d0/b) &
          - x0/Esqrt(x0**2-b**2)*( ArcSin((b**2-x0)/(b*(1.d0-x0))) &
          + ArcSin(b/x0) ) )
     !if(-1.d0/x0*1.d0/(1.d0-x0) < 0.d0) then
     !   intAz1 = intAz1 + (b**2-x0)*Pi
     !endif
   end function intAz1

   COMPLEX*16 FUNCTION FintAz1(x0,b)
     REAL*8, intent(in):: x0,b
     FintAz1 = fAz1(1.d0/(1.d0-x0),x0,b) - fAz1(1.d0/(-x0),x0,b)
   end function FintAz1

   COMPLEX*16 FUNCTION fAz1(t,x0,b)
     REAL*8, intent(in):: t,x0,b
     fAz1 = -Esqrt( ((b-x0)*t-1.d0)*((b+x0)*t+1.d0))/t &
          + (b**2-x0)*( ArcSin( (-1.d0-x0*t)/(b*t)) &
          -x0/Esqrt(x0**2-b**2)*ArcSin( (-x0+(b**2-x0**2)*t)/b) )
   end function fAz1

END SUBROUTINE OneLoop
