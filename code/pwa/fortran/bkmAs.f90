SUBROUTINE BornCS(x,t,q,kappa,A)
  REAL*8 x,t,kappa,q,omega,costh
  COMPLEX*16, intent(inout) :: A(8)
  REAL*8 a1,a2,a3,a4,a5,a6
  
  !kappa=0.d0
  omega = x
  costh = t/(2.d0*omega**2)+1.d0
  a1 = -q**2
  a2 = 0.d0
  a3 = omega/(2.d0)*((q + 2.d0*kappa)*q - (q + kappa)**2*costh )
  a4 = -(q + kappa)**2/(2.d0)/omega
  a5 = (q + kappa)**2/(2.d0)/omega
  a6 = -q*(q + kappa)/(2.d0)/omega

  A(1)=A(1)+a1
  A(2)=A(2)+a2
  A(3)=A(3)+a3
  A(4)=A(4)+a4      
  A(5)=A(5)+a5
  A(6)=A(6)+a6

end subroutine BornCS


