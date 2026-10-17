program running_table
  !* [v2] tabulate the running Delta parameters of As.f90 (SUBROUTINE delta_running) for the subtraction schemes
  !* irun = 1, 2, 3: Re gM(s), Re gE(s), M_Delta(s)/M_N, 1 + Re tau_R(s), and the vertex Im parts.
  !* Build (as in tests/regression.py):  gfortran-16 -O2 -ffree-line-length-none -std=legacy running_table.f90
  !*        build-gfortran-native/As.o -o running_table
  !* usage: running_table [vxmode [W1 W2 ...]]  (vxmode default 2 = Watson; 0 = PV2006 on-shell mapping, as
  !*        tools/running_delta.py; optional list of W values [GeV] replaces the default grid, up to 64 values)
  !* output columns: W[M_N units]  irun  ReG1  ReG2  Mrun[M_N]  ztau  ImgM  ImgE   for W on the grid below
  implicit none
  REAL*8 gA,hA,fpi,pi0mass,pionm
  COMMON/Consts/hA,gA,fpi,pi0mass,pionm
  INTEGER ivxmode
  COMMON /v2vx/ ivxmode
  REAL*8 aem_c
  COMMON /v2aem/ aem_c
  INTEGER irun_c
  COMMON /v2drun/ irun_c
  REAL*8, PARAMETER :: MN=0.93827d0
  REAL*8 Md, W(64), ReG1, ReG2, Mrun, ztau, ImgM, ImgE, s
  INTEGER i, ir, narg, nW
  CHARACTER(len=32) arg
  nW=8
  W(1:8)=(/0.93827d0, 0.93827d0+0.13957d0, 1.1d0, 1.15d0, 1.2d0, 1.232d0, 1.3d0, 1.4d0/)
  hA=2.85d0; gA=1.27d0; fpi=0.09221d0/MN; pionm=0.13957d0/MN; pi0mass=0.13498d0/MN
  Md=1.232d0/MN; ivxmode=2; aem_c=1.d0/137.036d0
  narg=command_argument_count()
  IF(narg.GE.1) THEN
     CALL get_command_argument(1,arg); READ(arg,*) ivxmode
  ENDIF
  IF(narg.GE.2) THEN
     nW=MIN(narg-1,64)
     DO i=1,nW
        CALL get_command_argument(i+1,arg); READ(arg,*) W(i)
     ENDDO
  ENDIF
  DO ir=1,3
     irun_c=ir
     DO i=1,nW
        s=(W(i)/MN)**2
        CALL delta_running(s,Md,2.97d0,-1.d0,ReG1,ReG2,Mrun,ztau)
        ImgM=0.d0; ImgE=0.d0
        IF(W(i)/MN.GT.1.d0+pionm) CALL gnd_vertex_im(W(i)/MN,Md,ImgM,ImgE)
        WRITE(*,'(f12.8,i3,6es24.15)') W(i)/MN, ir, ReG1, ReG2, Mrun, ztau, ImgM, ImgE
     ENDDO
  ENDDO
end program
