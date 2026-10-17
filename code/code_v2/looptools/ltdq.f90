!* ltdq.f90 -- double-precision front end of LoopTools-dq: the usual LoopTools calls with double
!* arguments, evaluated by a quad-precision (REAL*16) LoopTools core whose public names carry the prefix q
!* (built by build_looptools_dq.sh). Only the routines below are wrapped; any other LoopTools routine is
!* absent from the library and gives a link error, never a silent precision mix.
!* shortcut: no Xget/Xval cache access, no complex-mass (..C) or five-point routines; add wrappers when needed.

double complex function A0(m)
  double precision m
  complex(16), external :: qA0
  A0 = qA0(real(m,16))
end function

double complex function B0(p, m1, m2)
  double precision p, m1, m2
  complex(16), external :: qB0
  B0 = qB0(real(p,16), real(m1,16), real(m2,16))
end function

double complex function B1(p, m1, m2)
  double precision p, m1, m2
  complex(16), external :: qB1
  B1 = qB1(real(p,16), real(m1,16), real(m2,16))
end function

double complex function DB1(p, m1, m2)
  double precision p, m1, m2
  complex(16), external :: qDB1
  DB1 = qDB1(real(p,16), real(m1,16), real(m2,16))
end function

double complex function B0i(id, p, m1, m2)
  integer id
  double precision p, m1, m2
  complex(16), external :: qB0i
  B0i = qB0i(id, real(p,16), real(m1,16), real(m2,16))
end function

double complex function C0(p1, p2, p1p2, m1, m2, m3)
  double precision p1, p2, p1p2, m1, m2, m3
  complex(16), external :: qC0
  C0 = qC0(real(p1,16), real(p2,16), real(p1p2,16), real(m1,16), real(m2,16), real(m3,16))
end function

double complex function C0i(id, p1, p2, p1p2, m1, m2, m3)
  integer id
  double precision p1, p2, p1p2, m1, m2, m3
  complex(16), external :: qC0i
  C0i = qC0i(id, real(p1,16), real(p2,16), real(p1p2,16), real(m1,16), real(m2,16), real(m3,16))
end function

double complex function D0(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
  double precision p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4
  complex(16), external :: qD0
  D0 = qD0(real(p1,16), real(p2,16), real(p3,16), real(p4,16), real(p1p2,16), real(p2p3,16), &
           real(m1,16), real(m2,16), real(m3,16), real(m4,16))
end function

double complex function D0i(id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
  integer id
  double precision p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4
  complex(16), external :: qD0i
  D0i = qD0i(id, real(p1,16), real(p2,16), real(p3,16), real(p4,16), real(p1p2,16), real(p2p3,16), &
             real(m1,16), real(m2,16), real(m3,16), real(m4,16))
end function

subroutine ltini()
  call qltini()
end subroutine

subroutine ltexi()
  call qltexi()
end subroutine

subroutine clearcache()
  call qclearcache()
end subroutine

subroutine setdelta(x)
  double precision x
  call qsetdelta(real(x,16))
end subroutine

double precision function getdelta()
  real(16), external :: qgetdelta
  getdelta = qgetdelta()
end function

subroutine setmudim(x)
  double precision x
  call qsetmudim(real(x,16))
end subroutine

double precision function getmudim()
  real(16), external :: qgetmudim
  getmudim = qgetmudim()
end function

subroutine setlambda(x)
  double precision x
  call qsetlambda(real(x,16))
end subroutine

double precision function getlambda()
  real(16), external :: qgetlambda
  getlambda = qgetlambda()
end function

subroutine setuvdiv(x)
  double precision x
  call qsetuvdiv(real(x,16))
end subroutine

double precision function getuvdiv()
  real(16), external :: qgetuvdiv
  getuvdiv = qgetuvdiv()
end function

subroutine setminmass(x)
  double precision x
  call qsetminmass(real(x,16))
end subroutine

double precision function getminmass()
  real(16), external :: qgetminmass
  getminmass = qgetminmass()
end function

subroutine setversionkey(k)
  integer k
  call qsetversionkey(k)
end subroutine

integer function getversionkey()
  integer, external :: qgetversionkey
  getversionkey = qgetversionkey()
end function
