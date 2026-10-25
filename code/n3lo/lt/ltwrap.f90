! ltwrap.f90 -- bind(C) access to the quad-precision LoopTools-dq core for Python tests (work/n3lo)
! fin(k): eps^0 coefficient (with the current Delta), uv(k): coefficient of 1/eps; order as in coeffs.txt
subroutine lt_init() bind(C,name='lt_init')
  call qltini()
end subroutine
subroutine lt_setmudim(m2) bind(C,name='lt_setmudim')
  real*8, intent(in) :: m2
  call qsetmudim(real(m2,16))
end subroutine
subroutine lt_setdelta(d) bind(C,name='lt_setdelta')
  real*8, intent(in) :: d
  call qsetdelta(real(d,16))
end subroutine
subroutine lt_clear() bind(C,name='lt_clear')
  call qclearcache()
end subroutine
subroutine lt_bget(p2,m1,m2,fin,uv) bind(C,name='lt_bget')
  real*8, intent(in) :: p2,m1,m2
  complex*16, intent(out) :: fin(11), uv(11)
  integer, parameter :: id(11) = (/1,4,7,10,13,16,19,22,25,28,31/)
  complex(16), external :: qB0i
  integer k
  do k=1,11
     fin(k)=qB0i(id(k),real(p2,16),real(m1,16),real(m2,16))
     uv(k)=qB0i(id(k)+1,real(p2,16),real(m1,16),real(m2,16))
  enddo
end subroutine
subroutine lt_cget(p1,p2,p3,m1,m2,m3,fin,uv) bind(C,name='lt_cget')
  real*8, intent(in) :: p1,p2,p3,m1,m2,m3
  complex*16, intent(out) :: fin(22), uv(22)
  integer, parameter :: id(22) = (/1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64/)
  complex(16), external :: qC0i
  integer k
  do k=1,22
     fin(k)=qC0i(id(k),real(p1,16),real(p2,16),real(p3,16),real(m1,16),real(m2,16),real(m3,16))
     uv(k)=qC0i(id(k)+1,real(p1,16),real(p2,16),real(p3,16),real(m1,16),real(m2,16),real(m3,16))
  enddo
end subroutine
subroutine lt_dget(p1,p2,p3,p4,p12,p23,m1,m2,m3,m4,fin,uv) bind(C,name='lt_dget')
  real*8, intent(in) :: p1,p2,p3,p4,p12,p23,m1,m2,m3,m4
  complex*16, intent(out) :: fin(80), uv(80)
  integer, parameter :: id(80) = (/1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91,94,97,100,103,106,109,112,115,118,121,124,127,130,133,136,139,142,145,148,151,154,157,160,163,166,169,172,175,178,181,184,187,190,193,196,199,202,205,208,211,214,217,220,223,226,229,232,235,238/)
  complex(16), external :: qD0i
  integer k
  do k=1,80
     fin(k)=qD0i(id(k),real(p1,16),real(p2,16),real(p3,16),real(p4,16),real(p12,16),real(p23,16),real(m1,16),real(m2,16),real(m3,16),real(m4,16))
     uv(k)=qD0i(id(k)+1,real(p1,16),real(p2,16),real(p3,16),real(p4,16),real(p12,16),real(p23,16),real(m1,16),real(m2,16),real(m3,16),real(m4,16))
  enddo
end subroutine
subroutine lt_setversion(k) bind(C,name='lt_setversion')
  integer, intent(in) :: k
  call qsetversionkey(k)
end subroutine
