!! check_pi0ff.F90 -- cross-check of the closed-form pi0 -> gamma gamma transition form factor of the proton
!! triangle loop used by ampl.f90 (&v2 pi0ff, FUNCTION pi0_formfactor) against LoopTools' scalar C0:
!!     F(t) = -2 M^2 C0(0, 0, t; M^2, M^2, M^2)     (both photons real, pion momentum^2 = t; F(0) = 1)
!!     t < 0:      F = (4 M^2/(-t)) arsinh^2( sqrt(-t)/(2M) )
!!     0 < t < 4M^2: F = (4 M^2/t)  arcsin^2( sqrt(t)/(2M) )
!! Build and run (from work/code_v2):
!!     LT=$HOME/.local/looptools-gfortran-native
!!     SDKROOT=$(xcrun --sdk macosx --show-sdk-path) /opt/homebrew/bin/gfortran-16 -O2 -cpp -ffree-form \
!!        -ffree-line-length-none -std=legacy -I$LT/include tools/check_pi0ff.F90 $LT/lib/libooptools.a -o /tmp/check_pi0ff
!!     /tmp/check_pi0ff
!! Output: t [GeV^2], E_lab/theta of the point (where applicable), F(LoopTools), F(closed form), difference, F-1.
PROGRAM check_pi0ff
  IMPLICIT NONE
#include "looptools.h"
  REAL*8, PARAMETER :: M=0.93827d0
  REAL*8 t, Flt, Fcf, r, s, wcm, tmin
  COMPLEX*16 c0v
  INTEGER i
  REAL*8 tlist(9), elist(8)
  CHARACTER(len=40) lab(9)
  CALL ltini
  CALL setdelta(0.d0)
  tlist=(/ -1.d-4, -1.d-2, -0.05d0, -0.1d0, -0.2d0, -0.5d0, -2.d0, 0.1d0, 1.d0 /)
  lab=' '; lab(1)='small t'; lab(8)='timelike, below p pbar'; lab(9)='timelike, below p pbar'
  WRITE(*,'(a)') '#  t [GeV^2]        F(LoopTools)          F(closed form)        diff        F-1'
  DO i=1,9
     t=tlist(i)
     c0v=C0(0.d0,0.d0,t,M**2,M**2,M**2)
     Flt=-2.d0*M**2*dreal(c0v)
     Fcf=closed(t)
     WRITE(*,'(f12.6,2f22.15,es12.2,f12.6,2x,a)') t,Flt,Fcf,Flt-Fcf,Fcf-1.d0,TRIM(lab(i))
     IF(dabs(dimag(c0v)).GT.1.d-12) WRITE(*,*) '   (Im C0 =',dimag(c0v),')'
  ENDDO
  WRITE(*,'(a)') '# backward Compton kinematics (theta_cm = 180 deg, t_min = -4 w_cm^2):'
  WRITE(*,'(a)') '#  E_lab [MeV]   t_min [GeV^2]   F(LoopTools)       F(closed form)      F-1'
  elist=(/0.05d0,0.1d0,0.14d0,0.15d0,0.2d0,0.25d0,0.3d0,0.35d0/)
  DO i=1,8
     s=M**2+2.d0*M*elist(i)
     wcm=(s-M**2)/(2.d0*dsqrt(s))
     tmin=-4.d0*wcm**2
     c0v=C0(0.d0,0.d0,tmin,M**2,M**2,M**2)
     Flt=-2.d0*M**2*dreal(c0v)
     WRITE(*,'(f12.1,f16.6,2f20.12,f12.6)') 1.d3*elist(i),tmin,Flt,closed(tmin),closed(tmin)-1.d0
  ENDDO
CONTAINS
  REAL*8 FUNCTION closed(t)
    REAL*8, INTENT(in) :: t
    REAL*8 r
    IF(t.LT.0.d0) THEN
       r=dsqrt(-t)/(2.d0*M); closed=(asinh(r)/r)**2
    ELSE
       r=dsqrt(t)/(2.d0*M); closed=(dasin(r)/r)**2
    ENDIF
  END FUNCTION closed
END PROGRAM check_pi0ff
