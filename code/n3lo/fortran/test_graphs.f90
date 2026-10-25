program tg
  implicit none
  real*8 x,t,mpi2,MD
  complex*16 W(6), Z0
  integer ios,k
  mpi2=(139.57d0/938.27d0)**2; MD=1232.d0/938.27d0; Z0=(0.d0,0.d0)
  call ltini
  call setversionkey(64)
  do
    read(*,*,iostat=ios) x,t,k
    if(ios/=0) exit
    if(k==1) call dlg_d_piAB_D(x,t,0.d0,0.d0,mpi2,MD,Z0,Z0,Z0,Z0,W)
    if(k==2) call dlg_d_intA_piB_DD_ch(x,t,0.d0,0.d0,mpi2,MD,Z0,Z0,Z0,Z0,W)
    if(k==3) call dlg_d_intAB_DDD_chch(x,t,0.d0,0.d0,mpi2,MD,Z0,Z0,Z0,Z0,W)
    if(k==4) call dlg_d_V1AcD_V2BcD_D(x,t,0.d0,0.d0,mpi2,MD,Z0,Z0,Z0,Z0,W)
    if(k==5) call dlg_d_inA_piB_ND_ch(x,t,0.d0,0.d0,mpi2,MD,Z0,Z0,Z0,Z0,W)
    if(k==6) call dlg_d_intA_V1Bcpi_DD_ch(x,t,0.d0,0.d0,mpi2,MD,Z0,Z0,Z0,Z0,W)
    if(k==7) call dlg_d_sg_D(x,t,0.d0,0.d0,mpi2,MD,Z0,Z0,Z0,Z0,W)
    if(k==8) call dlg_d_inA_intB_NDD_chch(x,t,0.d0,0.d0,mpi2,MD,Z0,Z0,Z0,Z0,W)
    write(*,'(12es24.15)') (dble(W(ios)),dimag(W(ios)),ios=1,6)
  end do
end program
