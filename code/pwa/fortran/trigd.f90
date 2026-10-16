MODULE trigd
REAL(8), PARAMETER, PRIVATE :: PI = 3.141592653d0
CONTAINS
REAL *8 FUNCTION dcosd(x)
  REAL*8 x
  dcosd=cos(x*PI/180.d0)
END FUNCTION dcosd
REAL *8 FUNCTION  dsind(x)
  REAL*8 x
  dsind=sin(x*PI/180.)
END FUNCTION dsind
REAL *8 FUNCTION  dtand(x)
  REAL*8 x
  dtand=tan(x*PI/180.d0)
END FUNCTION dtand
REAL *8 FUNCTION  dacosd(x)
  REAL*8 x
  dacosd=180.d0/PI*acos(x)
END FUNCTION dacosd
END MODULE trigd
