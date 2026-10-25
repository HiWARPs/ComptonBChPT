!* p4_glue.f90 -- hand-written interface to the generated O(p^4) module (../n3lo/py/build_fortran.py)
SUBROUTINE p4_set_part(k)
  USE p4mod
  IMPLICIT NONE
  INTEGER, INTENT(in) :: k
  p4_part=k
END SUBROUTINE p4_set_part
