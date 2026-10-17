!!****************************************************************************
!!* Module: ampl.f90                   
!!* Created:  25/08/02             Update: Jan 06, Jan 08, Sep 08, Jan 12, Mar 13 
!!*
!!* Author: Vladimir Pascalutsa
!!* based on work in collaboration with 
!!* Olaf Scholten, Daniel R. Phillips, Vadim Lensky, Nadya Krupina)
!!*
!!* Calculates: Nucleon Compton scattering amplitudes and diff. cross-sections
!!* and other observables
!!* using the A1 ... A8 representation of the amplitude
!!*___________________________________________________________________________

      INCLUDE 'gauss.f90'
      INCLUDE 'dirac.f90'
      INCLUDE 'frules.f90'
      PROGRAM compton 

        USE dirac
        USE gauss
        USE trigd
        IMPLICIT NONE       
        REAL*8, PARAMETER :: smm=0.93827d0 ! protons mass in GeV
       !  REAL*8, PARAMETER :: hbarc=0.197327d0
       !  REAL*8, PARAMETER :: kapgp=1.7928474d0,kapgn=-1.913043d0
       ! INTEGER, PARAMETER :: Nbrn=3
       !  REAL*8, PARAMETER :: e2=1.d0/137.036d0
       !  REAL*8, PARAMETER :: alphaQED=1.d0/137.036d0
       REAL*8, PARAMETER :: hbarc=0.19733d0
        REAL*8, PARAMETER :: kapgp=1.793d0,kapgn=-1.913043d0
       INTEGER, PARAMETER :: Nbrn=3
        REAL*8, PARAMETER :: e2=1.d0/137.036d0
        REAL*8, PARAMETER :: alphaQED=1.d0/137.036d0
        REAL*8 :: hA=2.85
        REAL*8, PARAMETER :: alphaNNLOfm=11.2461d0,betaNNLOfm=3.87054d0, &
           gamEEfmNNLO=-3.34375,  gamMMfmNNLO=2.95934, gamEMfmNNLO=0.185504, gamMEfmNNLO=1.14328
        INTEGER  loopoints ! # of int point for HB calculation
        REAL*8, ALLOCATABLE :: xp(:),dxp(:)
       
        TYPE Baryon
           LOGICAL su(3)        ! include to s- (su1=T), u- (su2=T) channel;
           INTEGER ip, j, iso   ! parity, spin and isospin
           REAL*8 uiso          ! isospin factor in u-channel
           REAL*8 rm, rm2, wid  ! mass, mass^2, width
           REAL*8 cutoff        ! cutoff mass^2
           REAL*8 gnpi, gnpi2   ! gnpi=B-pi-N coupling, gnpi2=gnpi^2
           REAL*8 gngam, gngam2 ! gngam=N-B-photon coupling, gngam2=gngam^2
           LOGICAL rproc        ! T if B is to be renormalized
           REAL*8 z2,z1,rm0     ! w-f, c-c and mass renormalization
           REAL*8 aux1
        END TYPE Baryon

        TYPE(Baryon) bary(Nbrn)
        CHARACTER*1 signc
        INTEGER ih,i,irho,irhop
        INTEGER kvar,lwave
        REAL*8  m(4), msq(4), PI, GeV2fm


        REAL*8 Emin,Emax,Edel, Eminmem,thminmem
        NAMELIST/energy/ Emin,Emax,Edel,kvar

        REAL*8 s,t,u,W,Ephat,wphat,Eqhat,wqhat,Tplab,Tglab,crossec,cserror,diffcrossec
        REAL*8 ccnst,ccnst0,phat,phat2,sthr,esq,cutN,cutpi,qhat,qhat2
        REAL*8 inv_nmass_2fm,polza2LECs
        REAL*8, DIMENSION(0:3) :: vps
        REAL*8 ph,Tglabmev
        COMMON/pinphase/ph
        COMPLEX*16 dummy,f(2,2),ee(2,2),mm(2,2),ee0(2,2),mm0(2,2)
        COMPLEX*16 csampl(2,2,2,2),csampl_err(2,2,2,2)

        LOGICAL exch(3)
        INTEGER nq(3)
        REAL*8 rm,wid,cutf,gstr,gem
        NAMELIST/hadr/ exch,nq,rm,wid,cutf,gstr,gem
        LOGICAL anomaly,oneNloop,oneDloop,relchpt
        NAMELIST /incl/anomaly,oneNloop,oneDloop,relchpt,loopoints
        LOGICAL gaugeinv
        REAL*8 dng1,dng2,dg1,dg2,dz1,dz2,REM,Del
        NAMELIST /Delta/ gaugeinv,dng1,dng2,dg1,dg2,dz1,dz2
        COMMON /pardnga/dng1,dng2

        REAL*8 thmin,thmax,thdel,thetad,th,costh,thLAB
        LOGICAL cmsangle
        NAMELIST/angles/thmin,thmax,thdel,cmsangle
        REAL*8 gA,fpi,pi0mass,pionm2,pionm,exchm
        NAMELIST/constants/gA,fpi,pi0mass,pionm
        COMMON/Consts/hA,gA,fpi,pi0mass,pionm
        LOGICAL hbchpt,born,nuclamm,onlyLEX
        NAMELIST /nucleon/ hbchpt,born,nuclamm,onlyLEX
        LOGICAL pwa,proton,feynrule
        INTEGER jmax
        NAMELIST /multipoles/ pwa,jmax
        NAMELIST /waytocompute/ feynrule
        NAMELIST /isospin/ proton 
        REAL*8 alphaE1,betaM1,gamEE,gamEM,gamME,gamMM,gam0,gampi
        REAL*8 calpha,cbeta,cgamEE,cgamEM,cgamME,cgamMM,cgam0,cgampi
        REAL*8 alphaE1fm,betaM1fm,gamEEfm,gamEMfm,gamMEfm,gamMMfm
		LOGICAL varyscalarPs
        NAMELIST /scalarPs/ alphaE1fm,betaM1fm,varyscalarPs
	    NAMELIST /spinPs/gamEEfm,gamEMfm,gamMEfm,gamMMfm
		INTEGER iscalarPs,NscalarPs
		Real xscalarPs,SPinterval,SPdel

        Real, dimension(4) :: plzty_loop,plzty_delta, plzty_wzw

        !* [2026-10] optional fit-driver mode (see driver/README.md).
        !* If the namelist &fitdriver is absent from in.dat, nothing changes.
        LOGICAL pointsmode, selftest, staticmode, helmode
        CHARACTER(len=256) pointsfile, outfile
        REAL*8 refalphaE1, refbetaM1, refgamEE, refgamMM, refgamEM, refgamME
        INTEGER ios_fd
        NAMELIST /fitdriver/ pointsmode, pointsfile, outfile, selftest, staticmode, helmode, &
             refalphaE1, refbetaM1, refgamEE, refgamMM, refgamEM, refgamME

        PI=4.d0*datan(1.d0)

		

        OPEN (55,file='baryon.dat',status='old')


        DO ih=1,Nbrn
           READ (55,hadr)
           DO i=1,3
              bary(ih)%su(i) = exch(i)
           ENDDO
           bary(ih)%ip = nq(1)
           bary(ih)%j = nq(2)
           bary(ih)%iso = nq(3)
           bary(ih)%rm = rm/smm
           !!           IF(ih.EQ.3) bary(ih)%rm=(1.232d0/smm)
           bary(ih)%rm2= bary(ih)%rm**2
           bary(ih)%wid = wid/smm
           bary(ih)%cutoff = (cutf/smm)**2
           bary(ih)%gnpi2 = gstr
           IF(bary(ih)%j.GE.3) &
                bary(ih)%gnpi2=gstr*2.d0/3.d0
           bary(ih)%gnpi = dsqrt(gstr)
           bary(ih)%rproc=.FALSE.
           bary(ih)%aux1= gem
        ENDDO
        CLOSE(55)

        OPEN (5,file='in.dat',status='old')
      

        ncosx=16 ! set the number of angular-integration Gaussian points
        READ (5,energy)
        READ (5,angles)
        READ (5,incl)
        READ (5,Delta)
        READ (5,constants)
        READ (5,nucleon)
        READ (5,multipoles)
        READ (5,waytocompute)
        READ (5,isospin)
        READ (5,scalarPs)
        READ (5,spinPs)
        !* [2026-10] optional fit-driver namelist (defaults = off)
        pointsmode=.false.; selftest=.false.; staticmode=.false.; helmode=.false.
        pointsfile='points.txt'; outfile='points.out'
        refalphaE1=alphaNNLOfm; refbetaM1=betaNNLOfm
        refgamEE=gamEEfmNNLO; refgamMM=gamMMfmNNLO
        refgamEM=gamEMfmNNLO; refgamME=gamMEfmNNLO
        READ (5,fitdriver,IOSTAT=ios_fd)
        IF(ios_fd.GT.0) THEN
           WRITE(*,*) 'ERROR reading namelist &fitdriver in in.dat, iostat=',ios_fd
           STOP
        ENDIF
        IF(ios_fd.LT.0) THEN
           pointsmode=.false.; staticmode=.false.
        ENDIF

! set the polarizabilities
! go from fm units to nucleon mass units for polarizabilities
	    inv_nmass_2fm=(1.d0/smm)*hbarc
		polza2LECs=1.d0/(alphaQED*inv_nmass_2fm**3*1.d4)
		
		alphaE1= alphaE1fm*polza2LECs
        betaM1= betaM1fm*polza2LECs
        calpha= (alphaE1fm-alphaNNLOfm)*polza2LECs
        cbeta=(betaM1fm-betaNNLOfm )*polza2LECs

       gamEE=gamEEfm*polza2LECs/inv_nmass_2fm
       cgamEE=(gamEEfm-gamEEfmNNLO)*polza2LECs/inv_nmass_2fm
       gamMM=gamMMfm*polza2LECs/inv_nmass_2fm
       cgamMM=(gamMMfm-gamMMfmNNLO)*polza2LECs/inv_nmass_2fm
       gamEM=gamEMfm*polza2LECs/inv_nmass_2fm
       cgamEM=(gamEMfm-gamEMfmNNLO)*polza2LECs/inv_nmass_2fm
       gamME=gamMEfm*polza2LECs/inv_nmass_2fm
       cgamME=(gamMEfm-gamMEfmNNLO)*polza2LECs/inv_nmass_2fm
  
       gam0=-(gamEEfm+gamEMfm+gamMEfm+gamMMfm)
       gampi=-(gamEEfm+gamEMfm-gamMEfm-gamMMfm)
	          write(*,*) "forward and backward SPs", real(gam0), real(gampi)

	   cgam0= (gam0 +(gamEEfmNNLO+gamEMfmNNLO+gamMEfmNNLO+gamMMfmNNLO))*polza2LECs/inv_nmass_2fm
       cgampi= (gampi +(gamEEfmNNLO+gamEMfmNNLO-gamMEfmNNLO-gamMMfmNNLO))*polza2LECs/inv_nmass_2fm
        gam0=gam0*polza2LECs/inv_nmass_2fm
        gampi=gampi*polza2LECs/inv_nmass_2fm
		!write(*,*) 'LECs = ', calpha,cbeta,cgamEE,cgamEM,cgamME,cgamMM,cgam0,cgampi
	   If(onlyLEX) then
		write(*,*) 'Option onlyLEX is ON: Running only the low-energy expansion = Born+anomaly+polarizabilties'
		write(*,*) 'See in.dat for polarizability values and further info'
		bary(3)%su(1)=.false.
		bary(3)%su(2)=.false.
		oneNloop=.false.
		oneDloop=.false.
	    calpha=alphaE1
		cbeta=betaM1
		cgamEE=gamEE
		cgamEM=gamEM
		cgamME=gamME
		cgamMM=gamMM
		cgam0=gam0
		cgampi=gampi
			!write(*,*) calpha,cbeta,cgamEE,cgamEM,cgamME,cgamMM
	   endif
	   If(bary(1)%su(1).and.bary(1)%su(2).and.bary(3)%su(1).and.bary(3)%su(2).and.anomaly.and.oneNloop.and.oneDloop.and.relchpt) then
	   write(*,*) "Complete NNLO BChPT calculation running..."
       !write(*,*) "forward and backward SPs", real(gam0), real(gampi)
	   else
	   write(*,*) "Warning: running INCOMPLETE NNLO BChPT calculation, setting all polarizability LECs = 0"
	     if(.not.bary(1)%su(1)) write(*,*) "Missing the nucleon s-channel graph"
		 if(.not.bary(1)%su(2)) write(*,*) "Missing the nucleon u-channel graph"
		 if(.not.anomaly) write(*,*) "Missing the pi0 anomaly graph"	
		 if(.not.bary(3)%su(1)) write(*,*) "Missing the Delta(1232) s-channel graph"
		 if(.not.bary(3)%su(2)) write(*,*) "Missing the Delta(1232) u-channel graph"
		 if(.not.oneNloop) write(*,*) "Missing p-qube piN loops"
		 if(.not.oneDloop) write(*,*) "Missing p-forth/Delta piDelta loops"
		 if(.not.relchpt) write(*,*) "Invoking the Heavy-Baryon expansion: incorrect above the pion-productio threshold!"
		calpha=0.
		cbeta=0.
		cgamEE=0.
		cgamEM=0.
		cgamME=0.
		cgamMM=0.
		cgam0=0.
		cgampi=0.
	   endif


        IF(oneNloop.OR.oneDloop) THEN
           if(relchpt) then 
			  write(*,*) 'Initiating LoopTools...'
              call ltini    !* Initiate LoopTools
           else
              IF(.NOT.ALLOCATED(xp)) ALLOCATE(xp(1:loopoints),dxp(1:loopoints))
              CALL gauleg(0.d0,1.d0,xp(1),dxp(1),loopoints)
           endif
        ENDIF
 
        DO ih=1,Nbrn
           bary(ih)%z1=1.d0
           bary(ih)%z2=1.d0
           bary(ih)%rm0=bary(ih)%rm
        ENDDO
 
        
        Emin=Emin/smm
        Emax=Emax/smm
        Edel=Edel/smm
        fpi=fpi/smm
        pi0mass=pi0mass/smm
        m(1)=1.d0
        m(2)=cutN/smm
        m(3)=pionm/smm
        m(4)=cutPi/smm
        msq(1:4)=m(1:4)**2
        
        pionm2=msq(3)
        pionm=m(3)
 
        !!
        !! now everything is in the nucleon mass units
        !!
        WRITE(*,*) 'masses=',(m(i),i=1,4)

        ccnst=m(3)*dsqrt(PI)

        !dng1=dng1*dsqrt(2.d0/3.d0)
        !dng2=dng2*dsqrt(2.d0/3.d0)
        
        ! ratio E2/M1 from the gND coupling
        Del=bary(3)%rm-1.d0 ! Delta-Nulceon mass-difference
        REM = (Del/(bary(3)%rm+1.d0))**2*dng2/dng1
        !REM = Del*(dng2-Del*dng1)/( 2.d0*(bary(3)%rm+1.d0)*dng1-(1.d0+Del)*Del*dng2)
        write(*,*) 'REM = ', real(REM*100.), ' %'

        ! polarizabilities: 1-electric, 2-magnetic, 3-forward, 4-backward
        GeV2fm= hbarc/smm
        plzty_loop(1) = e2*5.*Pi/(6.*pionm)*(gA/(4.*Pi*fpi))**2*(hbarc/smm)**3*1.d4
        plzty_loop(2) = plzty_loop(1)/10.d0
        plzty_loop(3) = e2*2./(3.*pionm2)*(gA/(4.*Pi*fpi))**2*(hbarc/smm)**4*1.d4
        plzty_loop(4) = plzty_loop(3)
        plzty_wzw=0.
        plzty_wzw(4)= -e2*2.*gA/(2.*Pi*fpi)**2/pi0mass**2*(hbarc/smm)**4*1.d4
        plzty_delta(2) =2.*e2*dng1**2/(2.d0+Del)**2/Del*(hbarc/smm)**3*1.d4
        plzty_delta(1) =-2.*e2*dng2**2/(2.d0+Del)**3*(hbarc/smm)**3*1.d4
        plzty_delta(3) =-2.*e2/(2.d0+Del)**3/Del*(dng1**2/Del-dng2**2/(2.+Del)) &
            *(hbarc/smm)**4*1.d4
         plzty_delta(4) =2.*e2/(2.d0+Del)**3/Del*(5.*dng1**2/Del+dng2**2/(2.+Del)) &
            *(hbarc/smm)**4*1.d4
         OPEN (57,file='xpolarizab.out')
		write(57,*) "Heavy-Baryon Polarizabilities"
        do i=1,4
           write(57,*) plzty_loop(i),plzty_delta(i),plzty_wzw(i) &
                ,plzty_loop(i)+plzty_delta(i)+plzty_wzw(i)
        enddo
        CLOSE(57)

        if(pwa) then
           OPEN (501,file='aE1.out')
           OPEN (502,file='bM1.out')
           OPEN (503,file='aE2.out')
           OPEN (504,file='bM2.out')
           OPEN (505,file='gaEE.out')
           OPEN (506,file='gaMM.out')
           OPEN (507,file='gaEM.out')
           OPEN (508,file='gaME.out')
		   cmsangle=.true.  !integration over lab angle is unstable, switch to cms angle
        endif

        CALL angular_points
        CALL ini_gamma

        !* [2026-10] fit-driver / static-polarisability modes; the original
        !* main loops below are skipped only if one of these modes is requested.
        IF(pointsmode.OR.staticmode.OR.helmode) THEN
           IF(staticmode) CALL run_static()
           IF(pointsmode) CALL run_points()
           IF(helmode) CALL run_helpoints()
           if((oneNloop.or.oneDloop).and.relchpt) call ltexi
           STOP
        ENDIF

        if(varyscalarPs) then
			write(*,*) "Due to option varyscalarPs=t, computing only the scalar polarizability dependence for fixed E and Theta" 
            !Eminmem=Emin
			thminmem=thmin	
			thmax=thmin+thdel*0.1
			Emax=Emin+Edel*0.1
			OPEN (420,file='xBa_Alpha_dependence.output')
			OPEN (430,file='xBa_Beta_dependence.output')
			OPEN (320,file='xcs_Alpha_dependence.output')
			OPEN (330,file='xcs_Beta_dependence.output')
			OPEN (370,file='xcs_polar_Alpha_dependence.output')
			OPEN (380,file='xcs_polar_Beta_dependence.output')
			OPEN (390,file='xcs_beta_Alpha_dependence.output')
			OPEN (400,file='xcs_beta_Beta_dependence.output')
			OPEN (350,file='xBaLEX_Alpha_dependence.output')
			OPEN (360,file='xBaLEX_Beta_dependence.output')
			OPEN (325,file='xcsLEX_Alpha_dependence.output')
			OPEN (335,file='xcsLEX_Beta_dependence.output')
			OPEN (375,file='xcsLEX_polar_Alpha_dependence.output')
			OPEN (385,file='xcsLEX_polar_Beta_dependence.output')
			OPEN (395,file='xcsLEX_beta_Alpha_dependence.output')
			OPEN (405,file='xcsLEX_beta_Beta_dependence.output')

			NscalarPs=2
			SPinterval=3.
			SPdel=0.1
			write(*,*) 'Dependence on alphaE1 computed around',real(alphaE1fm),'+/-',real(SPinterval) 
			write(*,*) 'Dependence on betaM1 computed around',real(betaM1fm),'+/-',real(SPinterval) 
		else
			NscalarPs=1
			SPinterval=0.
			SPdel=0.1
			OPEN (31,file='xA_lab_cms.output')
			OPEN (32,file='xcs_Edependence.output')
			OPEN (33,file='xcs_Adistribution.output')
			OPEN (37,file='xcs_polar_Edependence.output')
			OPEN (38,file='xcs_polar_Adistribution.output')
			OPEN (39,file='xcs_beta_Edependence.output')
			OPEN (40,file='xcs_beta_Adistribution.output')
			
			OPEN (42,file='xBa_Edependence.output')
			OPEN (43,file='xBa_Adistribution.output')
			OPEN (52,file='xTa_Edependence.output')
			OPEN (53,file='xTa_Adistribution.output')
			OPEN (62,file='xGa_Edependence.output')
			OPEN (63,file='xGa_Adistribution.output')
			OPEN (72,file='xSigma2x_Edependence.output')
			OPEN (73,file='xSigma2x_Adistribution.output')
			OPEN (82,file='xSigma2z_Edependence.output')
			OPEN (83,file='xSigma2z_Adistribution.output')
			OPEN (44,file='xcs_Phi_beampolariz.output')
			OPEN (34,file='xtotal.output')
			OPEN (340,file='xLEXtotal.output')
			!OPEN (35,file='xcs_Klein-Nishina_lab.output')
			!OPEN (35,file='xcs_PowellMaximon_Edependence.output')
			!OPEN (36,file='xcs_PowellMaximon_Adistribution.output')
			OPEN (351,file='xbaLEX_Edependence.output')
			OPEN (361,file='xbaLEX_Adistribution.output')
			OPEN (651,file='xSig2z_LEX_Edependence.output')
			OPEN (661,file='xSig2z_LEX_Adistribution.output')
			OPEN (671,file='xSig2x_LEX_Edependence.output')
			OPEN (681,file='xSig2x_LEX_Adistribution.output')
			OPEN (352,file='xcsLEX_Edependence.output')
			OPEN (362,file='xcsLEX_Adistribution.output')
			OPEN (353,file='xcsLEX_polar_Edependence.output')
			OPEN (363,file='xcsLEX_polar_Adistribution.output')
			OPEN (354,file='xcsLEX_beta_Edependence.output')
			OPEN (364,file='xcsLEX_beta_Adistribution.output')
		endif

        Eminmem=Emin	
		Do 777 iscalarPs=1,NscalarPs
		Do 777 xscalarPs=-SPinterval,SPinterval,SPdel
			!Write(*,*) ' xscalarPs =', xscalarPs
		if(varyscalarPs) then
			thmin=thminmem
			SELECT CASE(iscalarPs)
			CASE(1)
				alphaE1= (alphaE1fm+xscalarPs)*polza2LECs
				betaM1= betaM1fm*polza2LECs
			CASE(2)
				alphaE1= alphaE1fm*polza2LECs
				betaM1= (betaM1fm+xscalarPs)*polza2LECs
			CASE default
              WRITE(*,*)'Debug: NscalarPs is out of range'
              STOP
			END SELECT
			If(onlyLEX) then
				calpha=alphaE1
				cbeta=betaM1
			else
				calpha= alphaE1-alphaNNLOfm*polza2LECs
				cbeta=betaM1-betaNNLOfm*polza2LECs
			endif
		endif		

        DO WHILE (thmin.LE.thmax)
           WRITE(*,*)
           IF(cmsangle) THEN
              WRITE(*,*) '--- THETA = ', thmin
              thetad=thmin
              th=PI*thetad/180.d0
              costh=dcos(th)
              CALL ini_polar(th)
           ENDIF
        
        Emin=Eminmem
        DO WHILE (Emin.LE.Emax)
           SELECT CASE(kvar)
           CASE(1)
              Tplab = Emin
              s = msq(1) + 2.d0*m(1)*(m(3)+Tplab) + msq(3)
              W = dsqrt(s)         
           CASE(2)
              W = Emin
              s = W*W
              Tplab=(s-msq(1)-msq(3))/(2.*m(1)) - m(3)
           CASE(3)
              Tglab = Emin
              s = 2.*m(1)*Tglab + msq(1)
              W = dsqrt(s)   
              Tplab=(s-msq(1)-msq(3))/(2.d0*m(1)) - m(3)
           CASE(4)
              W=exchm
              s = W*W
              Tplab=(s-msq(1)-msq(3))/(2.d0*m(1)) - m(3)
              m(2)=Emin
              msq(2)=m(2)*m(2)
           CASE default
              WRITE(*,*)'ERROR: kvar is out of range'
              STOP
           END SELECT
          
           Tglab = (s-msq(1))/(2.d0*m(1))     ! lab photon energy
           !WRITE(*,*)'Tplab=',REAL(Tplab*1.e3*smm), &
           !     '  Tglab=',REAL(Tglab*1.e3*smm),'  W=',W*1.e3*smm,' [MeV]'
           IF(.NOT.cmsangle) THEN
              thetad=thmin
              thLAB=PI*thetad/180.d0
			  !write(*,*) 'printout 1', s,thLAB, th,costh
              !th=PI*dacosd(1.d0-(1.d0+2.d0*Tglab)*(1.d0-dcosd(thetad)) &
              !     /(1.d0+Tglab*(1.d0-dcosd(thetad)) ))/180.d0
			   costh=1.d0-(1.d0+2.d0*Tglab)*(1.d0-dcos(thLAB)) &
                   /(1.d0+Tglab*(1.d0-dcos(thLAB)) )
              th=acos(real(costh))
			  !write(*,*) 'printout 2', s,thLAB,th,costh,dacos(-1.d0)
              CALL ini_polar(th)
           ELSE
			thLAB=dacos(1.d0-(1.d0-costh)/(1.d0+Tglab*(1.d0+costh)) )
		   ENDIF
           !WRITE(*,*) 'thetaLAB = ', 180.d0/PI*thLAB,'  degrees'
           WRITE(31,*) 180.d0/PI*thLAB,180.d0/PI*th

           Eqhat=(s+1.d0)/(2.d0*W)
           wqhat=(s-1.d0)/(2.d0*W)
           qhat2=Eqhat**2-1.d0
           IF(qhat2.GE.1.d-14) THEN
              qhat=dsqrt(qhat2)
           ELSE
              qhat=0.d0
           ENDIF           
		   t=-2.d0*wqhat**2*(1.d0-costh)
           u=2.d0-t-s
			write(*,*) 's,t,u =',s*smm**2,t*smm**2,u*smm**2
           Ephat=(s+1.d0-pionm2)/(2.d0*W)
           wphat=(s-1.d0+pionm2)/(2.d0*W)
           esq=Ephat**2
           phat2=esq-1.d0
           sthr=(1.d0+m(3))**2
           bnd=.TRUE.
           IF(phat2.GE.1.d-14) THEN
              phat=dsqrt(phat2)
              IF(s.GT.sthr) bnd=.FALSE.
           ELSE
              phat=0.d0
           ENDIF
           vps(0)=W
           DO i=1,3
              vps(i)=0.d0  ! c.m.s
           ENDDO
           ccnst0=Ephat/W
           Tglabmev=Tglab*1.d3*smm
           if(feynrule) then
              ! CALL fr_gaga_ampl(qhat,th,csampl)
              write(*,*) 'Warning: Option feynrule=true has no effect' 
           endif
		   CALL gaga_ampl(qhat,th,csampl,csampl_err)
           CALL observables(csampl,csampl_err,crossec,cserror,diffcrossec,.true.)
		   CALL LEXobservables(.true.)
           if (pwa) then 
              CALL pwa_gaga(qhat)
           endif
           Emin=Emin+Edel
        ENDDO 
           thmin=thmin+thdel
        ENDDO 
777     CONTINUE
		
        if((oneNloop.or.oneDloop).and.relchpt) call ltexi    !* Close LoopTools

      CONTAINS


        SUBROUTINE gaga_ampl(qq,th,mgg,mgg_err) 
          !USE frules
          IMPLICIT NONE
          REAL*8, INTENT(in) :: qq,th
          COMPLEX*16, INTENT(out) :: mgg(2,2,2,2), mgg_err(2,2,2,2) 
          COMPLEX*16 sum
          INTEGER ih,irho,irhop
          INTEGER jb,isob,ita,hindx,hpindx,mu,nu
          INTEGER i,j
          REAL*8 rm,rm2,gg1,gg2
          REAL*8 rt,ru,costh,sinth,x,Pdotq,pnucl
          REAL*8, DIMENSION(0:3) :: P,vp,vq,vpp,vqp,vqf,vpt,vpu
          REAL*8 hel,help
          COMPLEX*16, DIMENSION(0:3,2) :: ep,e
          TYPE(spinor) ubar,uu
          TYPE(matrix4) curnt(0:3),curntp(0:3)
          REAL*8 gg,gpin
          TYPE(matrix4) amplgg(2,2), m_vqp_vq
          COMPLEX*16 A(8),scalar
          common/temp/ A

         
          !! c.m.s. kinematics
          vq(0)=qhat
          vq(1)=0.d0
          vq(2)=0.d0
          vq(3)=qhat
          vp(0)=Eqhat
          vp(1)=0.d0
          vp(2)=0.d0
          vp(3)=-qhat

          costh=dcos(th)
          sinth=dsin(th)
          vqp(0)=qhat
          vqp(1)=qhat*sinth
          vqp(2)=0.d0
          vqp(3)=qhat*costh
          vpp(0)=Eqhat
          vpp(1)=-qhat*sinth
          vpp(2)=0.d0
          vpp(3)=-qhat*costh
          vqf=-vqp

          P=vpp+vp
          vpu=vp-vqp
          vpt=vq-vqp
          ru=dot(vpu,vpu)
          rt=dot(vpt,vpt)
          x=(s-1.d0)*0.5d0
          Pdotq=(s-ru)*0.5d0

          pnucl=qhat
          A=0.d0
          IF (hbchpt) THEN
             CALL comp_bkmAs(qhat,rt,A)
             pnucl=0.d0
          ELSE             
             CALL computeAs(x,rt,A)
          ENDIF

            ! for testing purposes
            !A=0.d0; P=0.d0
            !!A(1)=-1.d0
            !A(1)=costh*2.
            !A(6)=1./qhat**2
             !A(6)=10.*costh/qhat**2
           !A(3)=costh*2.
            ! end testing text

          A=alphaQED*A

          DO mu=0,3
          DO i=1,2
             ep(mu,i)=vepsp(mu,i)- dot(P,vepsp(:,i))/(Pdotq) *vqp(mu)
             e(mu,i) =veps(mu,i) - dot(P,veps(:,i))/(Pdotq) *vq(mu)
          ENDDO
          ENDDO


          m_vqp_vq = ag_(vqp,vq)
          DO i=1,2
             DO j=1,2
                amplgg(i,j)= (-1.d0)*A(1)*dot(ep(:,i),e(:,j))*gone &
                     + A(2)*dot(vq,ep(:,i))*dot(vqp,e(:,j))*gone &
                     + (-1.d0)*A(3)*ag_(ep(:,i),e(:,j)) &
                     + A(4)*dot(ep(:,i),e(:,j))*m_vqp_vq &
                     - A(5)*( dot(vqp,e(:,j))*ag_(vq,ep(:,i)) &
                                  - dot(vq,ep(:,i))*ag_(vqp,e(:,j)) ) &
                     - A(6)*( dot(vqp,e(:,j))*ag_(vqp,ep(:,i)) &
                                  - dot(vq,ep(:,i))*ag_(vq,e(:,j)) ) &
                     + A(7)*dot(ep(:,i),vq)*dot(vqp,e(:,j))*m_vqp_vq &           
                     - A(8)*(0.d0,1.d0)*eps4_(ep(:,i),e(:,j),vq,vqp)*gamm(5)
             ENDDO
          ENDDO

! Intel compiler does not like this:
!!$         DO mu=0,3
!!$          DO i=1,2
!!$             ep(mu,i)=vepsp(mu,i)- (P*vepsp(:,i))/(Pdotq) *vqp(mu)
!!$             e(mu,i) =veps(mu,i) - (P*veps(:,i))/(Pdotq) *vq(mu)
!!$          ENDDO
!!$          ENDDO
!!$
!!$
!!$          m_vqp_vq = ag_(vqp,vq)
!!$          DO i=1,2
!!$             DO j=1,2
!!$                amplgg(i,j)= (-1.d0)*A(1)*(ep(:,i)*e(:,j))*gone &
!!$                     + A(2)*(vq*ep(:,i))*(vqp*e(:,j))*gone &
!!$                     + (-1.d0)*A(3)*ag_(ep(:,i),e(:,j)) &
!!$                     + A(4)*(ep(:,i)*e(:,j))*m_vqp_vq &
!!$                     - A(5)*( (vqp*e(:,j))*ag_(vq,ep(:,i)) &
!!$                                  - (vq*ep(:,i))*ag_(vqp,e(:,j)) ) &
!!$                     - A(6)*( (vqp*e(:,j))*ag_(vqp,ep(:,i)) &
!!$                                  - (vq*ep(:,i))*ag_(vq,e(:,j)) ) &
!!$                     + A(7)*(ep(:,i)*vq)*(vqp*e(:,j))*m_vqp_vq &           
!!$                     - A(8)*(0.d0,1.d0)*eps4_(ep(:,i),e(:,j),vq,vqp)*gamm(5)
!!$             ENDDO
!!$          ENDDO

          DO hindx=1,2
             hel=-(hindx-1.5d0)
             uu=u_(pnucl,0.d0+PI,0.d0,hel,1)
               !call prnspinor(uu)
             DO hpindx=1,2
                help=-(hpindx-1.5d0)
                ubar=ubar_(pnucl,th+PI,0.d0,help,1)
               !call prnspinor(ubar)
               !write(9,*) hpindx,hindx,braket(ubar,gamm(5)*gamm(0)*(gamm(1)+ci*gamm(2)),uu)
               !write(10,*) hpindx,hindx,braket(ubar,gamm(5)*gamm(0)*(gamm(1)-ci*gamm(2)),uu)
               DO i=1,2
                   DO j=1,2
                      sum=braket(ubar,amplgg(i,j),uu)
                     ! IF(dabs(dreal(sum)).LT.1.e-10) sum=(0.d0,1.d0)*dimag(sum)
                      mgg(i,j,hpindx,hindx)=sum
                       if (x.gt.0.2) then
                          scalar=-dot(ep(:,i),e(:,j))*x*alphaQED
                       else
                          scalar=-dot(ep(:,i),e(:,j))*alphaQED*(x**2+x**3/Del)/2.d0
                      endif
                     mgg_err(i,j,hpindx,hindx)=braket(ubar,gone*scalar,uu)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          
          !write(*,*) 
          RETURN
        END SUBROUTINE gaga_ampl
        


       SUBROUTINE comp_bkmAs(x,t,A)
          REAL*8, INTENT(in):: x,t
          COMPLEX*16, INTENT(inout):: A(8)
          REAL*8 q,gg,gafac,kappa
       
         if(proton) then
             Q=1.d0 ! 1-proton, 0-neutron
             kappa=kapgp
          else
             Q=0.d0
             kappa=kapgn
          endif
          IF(born) THEN
          if(.NOT.nuclamm) kappa=0.d0
            ! Born diagrams (nucleon s and u)
             CALL BornCS(x,t,q,kappa,A)
          ENDIF
          IF(anomaly) THEN
             gg=(2.d0*q - 1.d0)*gA/(2.d0*Pi*fpi)**2
             gg=-gg/(pi0mass**2-t)
             A(3)=A(3)+gg*t*x
             A(5)=A(5)+gg*x
             A(6)=A(6)-gg*x            
          ENDIF
          IF (oneNloop) THEN
             gafac=gA**2/(8.d0*Pi*fpi**2)
             CALL OneLoop(x,t,gafac,pionm,A,xp,dxp,loopoints)
          ENDIF
    
        END SUBROUTINE comp_bkmAs

        SUBROUTINE computeAs(x,t,A)
          REAL*8, INTENT(in):: x,t
          COMPLEX*16, INTENT(inout):: A(8)
          COMPLEX*16 A0(8),A00(8)
          REAL*8 y,gg,Md,Q,gafac,cutf2,kappa
          REAL*8 albe,beta,gam0,dimregDelta,xi
          INTEGER i,ih
          Logical, parameter :: id=.false.

          if(proton) then
             Q=1.d0 ! 1-proton, 0-neutron
             kappa=kapgp
          else
             Q=0.d0
             kappa=kapgn
          endif
         if(.NOT.nuclamm) kappa=0.d0


          y=-x-t*0.5d0
		  xi=0.5d0*(x-y)             
          ih=1 ! the nucleon pole diagrams
            !WRITE(*,*) 'x =', REAL(x)
            !WRITE(*,*) 't = ', REAL(t)

          A0=0.d0
          IF(bary(ih)%su(1)) THEN
             CALL s_nucl(x,t,Q,kappa,A0)
             A=A+A0
          ENDIF
          IF(bary(ih)%su(2)) THEN
             !CALL u_nucl(y,t,Q,kapgp,A)
             !breit frame:
             CALL s_nucl(y,t,Q,kappa,A0)
             A(1)=A(1)+A0(1)
             A(2)=A(2)+A0(2)
             DO i=3,7
                A(i)=A(i)-A0(i)
             ENDDO             
             A(8)=A(8)+A0(8)
          ENDIF

          ih=3 ! the Delta pole diagram
          Md=bary(ih)%rm
          gg=dng1
          ! Delta s-exchange
          IF(bary(ih)%su(1)) THEN
             !CALL s_Delta(x,t,Md,gg,A0)
             CALL s_Delta2(x,t,Md,pionm,dng1,dng2,A0)
             if(id) then
                ! identity of rewriting 8 into 6 check
                A0(3)=A0(3)-t*(x**2+t/2.)*A0(7)-t*(x/2.)*A0(8)
                A0(5)=A0(5)+(x**2+t/2.)*A0(7)+(x/2.)*A0(8)
                A0(6)=A0(6)-x**2*A0(7)-(x/2.)*A0(8)
                A0(7)=(0.d0,0.d0)
                A0(8)=(0.d0,0.d0)
              endif
              !A0(3:)=0.d0
            A=A+A0
          ENDIF
          IF(bary(ih)%su(2)) THEN
             !breit frame:
             !CALL s_Delta(y,t,Md,gg,A0)
             CALL s_Delta2(y,t,Md,pionm,dng1,dng2,A0)
             A(1)=A(1)+A0(1)
             A(2)=A(2)+A0(2)
              !A0(3:)=0.d0
             DO i=3,7
                A(i)=A(i)-A0(i)
             ENDDO             
             A(8)=A(8)+A0(8)
          ENDIF

          IF(anomaly) THEN
             gg=(2.d0*Q - 1.d0)*gA/(2.d0*Pi*fpi)**2
             !cutf2=(0.7/0.9383)**2
             !gg=gg*(cutf2)/(cutf2-t)
             !write(*,*) 'anomaly gg = ', gg
             A(8)=A(8) - gg/(pi0mass**2-t)
          ENDIF

          IF (oneNloop) THEN
             IF(relchpt) THEN
                gafac=gA**2/(4.d0*Pi*fpi)**2
                CALL CovOneLoop(1.d0+2.*x,t,pionm**2,gafac,A0)
                A=A+A0
                albe=dreal(A0(1)+0.5d0*t*A0(2))
                beta=-dreal(A0(2))
                gam0=dreal(A0(3))
                CALL CovOneLoop(1.d0+2.*y,t,pionm**2,gafac,A0)
                A(1)=A(1)+A0(1)
                A(2)=A(2)+A0(2)
                DO i=3,7
                   A(i)=A(i)-A0(i)
                ENDDO
                A(8)=A(8)+A0(8)
                albe=albe+dreal(A0(1)+0.5d0*t*A0(2))
                beta=beta-dreal(A0(2))
                gam0=gam0-dreal(A0(3))
                !write(11,*) x, t
                !write(11,*) 'alpha+beta= ', (1./137.)*albe/(-x*y)*inv_nmass_2fm**3*1.d4 
                !write(11,*) 'beta = ', (1./137.)*beta*inv_nmass_2fm**3*1.d4 
                !write(11,*) 'alpha-beta = ', (1./137.)*dreal(A(1))/(-x*y)*inv_nmass_2fm**3*1.d4 
	        !write(11,*) 'gamma0 = ', (1./137.)*(gam0/x**3)*inv_nmass_2fm**4*1.d4 
             ELSE
                gafac=gA**2/(8.d0*Pi*fpi**2)
              !CALL OneLoop(W-1.d0,t,gafac,pionm,A,xp,dxp,loopoints) !* prescription
              CALL OneLoop(qhat,t,gafac,pionm,A,xp,dxp,loopoints)  !* true heavy-baryon result
             ENDIF
          ENDIF
    
          IF (oneDloop) THEN
             IF(relchpt) THEN
                gafac=hA**2/(4.d0*Md*Pi*fpi)**2
                dimregDelta=0.d0
                CALL CorrIsoCovOneDeltaLoop(1.d0+2.*x,t,pionm**2,gafac,A0,dimregDelta)
                A=A+A0
                albe=dreal(A0(1)+0.5d0*t*A0(2))
                beta=-dreal(A0(2))
                gam0=dreal(A0(3))
                CALL CorrIsoCovOneDeltaLoop(1.d0+2.*y,t,pionm**2,gafac,A0,dimregDelta)
                A(1)=A(1)+A0(1)
                A(2)=A(2)+A0(2)
                DO i=3,7
                   A(i)=A(i)-A0(i)
                ENDDO
                A(8)=A(8)+A0(8)
                !CALL ZeroCovOneDeltaLoop(s,t,pionm**2,gafac,A00,dimregDelta)
                !A(1)=A(1)-A00(  A(5)=A(5) + 0.5*(x-y)*cgamMM1)
                !A(3)=A(3)-(x-y)*A00(3)

                albe=albe+dreal(A0(1)+0.5d0*t*A0(2))
                beta=beta-dreal(A0(2))
                gam0=gam0-dreal(A0(3))
               ! write(13,*) x, t
               ! write(13,*) 'D alpha+beta= ', (1./137.)*albe/(-x*y)*inv_nmass_2fm**3*1.d4 
               ! write(13,*) 'D beta = ', (1./137.)*beta*inv_nmass_2fm**3*1.d4 
                !write(13,*) 'D alpha-beta = ', (1./137.)*dreal(A0(1))/(-x*y)*inv_nmass_2fm**3*1.d4 
	        !write(13,*) 'D gamma0 = ', (1./137.)*(gam0/x**3)*inv_nmass_2fm**4*1.d4 
             ELSE
             call piDloop(qhat,t, A)
          ENDIF
          endif
          !albe=dreal(A(1)+0.5d0*t*A(2)) !+1.d0 ! -Born contribution
          !beta=-dreal(A(2))
          !write(12,*) x, t
          !write(12,*) 'alpha+beta= ', (1./137.)*albe/(qhat**2)*inv_nmass_2fm**3*1.d4
          !write(12,*) 'beta = ', (1./137.)*beta*inv_nmass_2fm**3*1.d4

          !* Add the O(p4) LECs contribution (scalar polarizabilities):
	    !write(*,*) 'LECs = ', calpha,cbeta,cgamEE,cgamEM,cgamME,cgamMM,cgam0,cgampi
          !* [2026-10] the contact-term formulas were moved verbatim into
          !* SUBROUTINE add_LECs below, so that the fit driver (run_points) uses
          !* exactly the same expressions (single source of truth).
          CALL add_LECs(x,t,A,calpha,cbeta,cgamEM,cgamME,cgamMM,cgam0)
          !* Add the O(p5) LECs contribution:
          !polza2LECs=1.d0/(alphaQED*inv_nmass_2fm**4*1.d4)
          !A(3)=A(3)+(calpha+cbeta)*x*(-y)+0.5d0*t*cbeta
          !A(2)=A(2)-cbeta
          
          !write(*,*) ' added ', calpha*(1./137.)*inv_nmass_2fm**3*1.d4,' [10^-4 fm^3] to alpha' 
          !write(*,*) ' added ', cbeta*(1./137.)*inv_nmass_2fm**3*1.d4,' [10^-4 fm^3] to beta'
        END SUBROUTINE computeAs

        SUBROUTINE add_LECs(x,t,A,ca,cb,cgEM,cgME,cgMM,cg0)
          !* O(p^4) contact (polarizability-LEC) contributions to the
          !* invariant amplitudes A1..A8 (nucleon-mass units, without alpha_em).
          !* Moved verbatim from computeAs [2026-10]; arithmetic unchanged.
          !*  ca, cb   : shifts of alpha_E1, beta_M1   (in units alpha_em/M^3)
          !*  cgEM,cgME,cgMM : shifts of gamma_E1M2, gamma_M1E2, gamma_M1M1
          !*  cg0      : shift of gamma_0 = -(gEE+gEM+gME+gMM)  (units alpha_em/M^4)
          REAL*8, INTENT(in) :: x,t,ca,cb,cgEM,cgME,cgMM,cg0
          COMPLEX*16, INTENT(inout) :: A(8)
          REAL*8 y,xi
          y=-x-t*0.5d0
          xi=0.5d0*(x-y)
          A(1)= A(1) + (0.5d0*(ca+cb)*(x**2+(-y)**2)+0.5d0*t*cb) &
                    + 0.5d0*xi**2*t*(cgMM-cgME)
          A(2)=A(2) - cb  - (1.d0/8.d0)*(ca+cb)*t &
                    - xi**2*(cgMM+cgEM)
          A(3)=A(3) + xi**3*cg0 - 0.5d0*xi*t*(cgMM+cgME)
          A(4)=A(4) + (1.d0-0.25d0*t)*xi*(-cgMM+cgME)
          A(5)=A(5) + xi*cgMM + (1.d0/8.d0)*xi*t*(cgEM+cgME)
          A(6)=A(6) + xi*cgEM - (1.d0/8.d0)*xi*t*(cgEM+cgME)
          A(7)=A(7) + 0.5d0*xi*(cgMM+cgEM)
          A(8)=A(8) - 0.5d0*xi**2*cg0
        END SUBROUTINE add_LECs

!!=========================================================================
!!  [2026-10]  FIT-DRIVER SUPPORT ROUTINES
!!  run_points : reads a list of kinematic points, writes for every point the
!!               observable as an exact quadratic form in the shifts of the six
!!               dipole polarisabilities (see driver/README.md)
!!  run_static : dumps the invariant amplitudes A1..A8, split by mechanism, on a
!!               small (nu,t) grid for the extraction of static polarisabilities
!!=========================================================================

        SUBROUTINE add_crossed(Ad,A0)
          !* adds the u-channel (crossed) partner exactly as done in computeAs
          COMPLEX*16, INTENT(inout) :: Ad(8)
          COMPLEX*16, INTENT(in) :: A0(8)
          Ad(1)=Ad(1)+A0(1)
          Ad(2)=Ad(2)+A0(2)
          Ad(3:7)=Ad(3:7)-A0(3:7)
          Ad(8)=Ad(8)+A0(8)
        END SUBROUTINE add_crossed

        SUBROUTINE components_As(x,t,Ac)
          !* Same physics as computeAs (covariant branch), split into
          !* 1 Born(N s+u), 2 pi0 anomaly, 3 piN loops, 4 Delta pole, 5 piDelta loops,
          !* 6 contact terms with the current LEC values. Sum_c Ac(:,c) = computeAs.
          REAL*8, INTENT(in) :: x,t
          COMPLEX*16, INTENT(out) :: Ac(8,6)
          COMPLEX*16 A0(8)
          REAL*8 y,Q,kappa,gg,gafac,Md
          Ac=(0.d0,0.d0)
          if(proton) then
             Q=1.d0; kappa=kapgp
          else
             Q=0.d0; kappa=kapgn
          endif
          if(.NOT.nuclamm) kappa=0.d0
          y=-x-t*0.5d0
          IF(bary(1)%su(1)) THEN
             A0=0.d0; CALL s_nucl(x,t,Q,kappa,A0); Ac(:,1)=Ac(:,1)+A0
          ENDIF
          IF(bary(1)%su(2)) THEN
             A0=0.d0; CALL s_nucl(y,t,Q,kappa,A0); CALL add_crossed(Ac(:,1),A0)
          ENDIF
          IF(anomaly) THEN
             gg=(2.d0*Q - 1.d0)*gA/(2.d0*Pi*fpi)**2
             Ac(8,2)=Ac(8,2) - gg/(pi0mass**2-t)
          ENDIF
          IF(oneNloop.AND.relchpt) THEN
             gafac=gA**2/(4.d0*Pi*fpi)**2
             CALL CovOneLoop(1.d0+2.*x,t,pionm**2,gafac,A0); Ac(:,3)=Ac(:,3)+A0
             CALL CovOneLoop(1.d0+2.*y,t,pionm**2,gafac,A0); CALL add_crossed(Ac(:,3),A0)
          ENDIF
          Md=bary(3)%rm
          IF(bary(3)%su(1)) THEN
             CALL s_Delta2(x,t,Md,pionm,dng1,dng2,A0); Ac(:,4)=Ac(:,4)+A0
          ENDIF
          IF(bary(3)%su(2)) THEN
             CALL s_Delta2(y,t,Md,pionm,dng1,dng2,A0); CALL add_crossed(Ac(:,4),A0)
          ENDIF
          IF(oneDloop.AND.relchpt) THEN
             gafac=hA**2/(4.d0*Md*Pi*fpi)**2
             CALL CorrIsoCovOneDeltaLoop(1.d0+2.*x,t,pionm**2,gafac,A0,0.d0); Ac(:,5)=Ac(:,5)+A0
             CALL CorrIsoCovOneDeltaLoop(1.d0+2.*y,t,pionm**2,gafac,A0,0.d0); CALL add_crossed(Ac(:,5),A0)
          ENDIF
          CALL add_LECs(x,t,Ac(:,6),calpha,cbeta,cgamEM,cgamME,cgamMM,cgam0)
        END SUBROUTINE components_As

        SUBROUTINE hel_from_A(thx,Ain,mgg)
          !* c.m. helicity amplitudes mgg(lambda',lambda,h',h) from invariant amplitudes
          !* Ain(1:8) (ALREADY multiplied by alphaQED). Verbatim copy of the projection
          !* in gaga_ampl (covariant branch); host s,Eqhat,qhat and ini_polar(thx) must be set.
          IMPLICIT NONE
          REAL*8, INTENT(in) :: thx
          COMPLEX*16, INTENT(in) :: Ain(8)
          COMPLEX*16, INTENT(out) :: mgg(2,2,2,2)
          REAL*8, DIMENSION(0:3) :: P,vp,vq,vpp,vqp,vpu
          REAL*8 cst,snt,Pdotq,pnucl,hel,help,ru
          COMPLEX*16, DIMENSION(0:3,2) :: ep,e
          TYPE(spinor) ubar,uu
          TYPE(matrix4) amplgg(2,2), m_vqp_vq
          INTEGER mu,i,j,hindx,hpindx
          vq(0)=qhat; vq(1)=0.d0; vq(2)=0.d0; vq(3)=qhat
          vp(0)=Eqhat; vp(1)=0.d0; vp(2)=0.d0; vp(3)=-qhat
          cst=dcos(thx); snt=dsin(thx)
          vqp(0)=qhat; vqp(1)=qhat*snt; vqp(2)=0.d0; vqp(3)=qhat*cst
          vpp(0)=Eqhat; vpp(1)=-qhat*snt; vpp(2)=0.d0; vpp(3)=-qhat*cst
          P=vpp+vp
          vpu=vp-vqp
          ru=dot(vpu,vpu)
          Pdotq=(s-ru)*0.5d0
          pnucl=qhat
          DO mu=0,3
          DO i=1,2
             ep(mu,i)=vepsp(mu,i)- dot(P,vepsp(:,i))/(Pdotq) *vqp(mu)
             e(mu,i) =veps(mu,i) - dot(P,veps(:,i))/(Pdotq) *vq(mu)
          ENDDO
          ENDDO
          m_vqp_vq = ag_(vqp,vq)
          DO i=1,2
             DO j=1,2
                amplgg(i,j)= (-1.d0)*Ain(1)*dot(ep(:,i),e(:,j))*gone &
                     + Ain(2)*dot(vq,ep(:,i))*dot(vqp,e(:,j))*gone &
                     + (-1.d0)*Ain(3)*ag_(ep(:,i),e(:,j)) &
                     + Ain(4)*dot(ep(:,i),e(:,j))*m_vqp_vq &
                     - Ain(5)*( dot(vqp,e(:,j))*ag_(vq,ep(:,i)) &
                                  - dot(vq,ep(:,i))*ag_(vqp,e(:,j)) ) &
                     - Ain(6)*( dot(vqp,e(:,j))*ag_(vqp,ep(:,i)) &
                                  - dot(vq,ep(:,i))*ag_(vq,e(:,j)) ) &
                     + Ain(7)*dot(ep(:,i),vq)*dot(vqp,e(:,j))*m_vqp_vq &
                     - Ain(8)*(0.d0,1.d0)*eps4_(ep(:,i),e(:,j),vq,vqp)*gamm(5)
             ENDDO
          ENDDO
          DO hindx=1,2
             hel=-(hindx-1.5d0)
             uu=u_(pnucl,0.d0+PI,0.d0,hel,1)
             DO hpindx=1,2
                help=-(hpindx-1.5d0)
                ubar=ubar_(pnucl,thx+PI,0.d0,help,1)
                DO i=1,2
                   DO j=1,2
                      mgg(i,j,hpindx,hindx)=braket(ubar,amplgg(i,j),uu)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
        END SUBROUTINE hel_from_A

        SUBROUTINE bilin_obs(Ma,Mb,scl,ob)
          !* Real symmetric bilinear forms B(Ma,Mb) such that B(M,M) reproduces the
          !* observables of SUBROUTINE observables (same conventions):
          !*  ob(1) = unpolarised sigma0        = SCALE*sum|M|^2/4
          !*  ob(2) = sigma_par  (SI_PAR*SCALE) , ob(3) = sigma_perp (SI_PERP*SCALE)
          !*  ob(4) = Sigma2x numerator = -SCALE*SI_Xm/2 ,  ob(5) = Sigma2z numerator = -SCALE*SI_Zm/2
          !*  ob(6) = SCALE*SI_Xp/2 , ob(7) = SCALE*SI_Zp/2  (= sigma0 by parity; checks only)
          COMPLEX*16, INTENT(in) :: Ma(2,2,2,2), Mb(2,2,2,2)
          REAL*8, INTENT(in) :: scl
          REAL*8, INTENT(out) :: ob(7)
          COMPLEX*16 xa,xb,ya,yb,sa,sb
          REAL*8 sq2
          INTEGER h,hp,i,j
          sq2=1.d0/dsqrt(2.d0)
          ob=0.d0
          DO h=1,2
          DO hp=1,2
          DO i=1,2
             DO j=1,2
                ob(1)=ob(1)+dreal(Ma(i,j,hp,h)*CONJG(Mb(i,j,hp,h)))/4.d0
             ENDDO
             xa=-sq2*(Ma(i,1,hp,h)-Ma(i,2,hp,h)); xb=-sq2*(Mb(i,1,hp,h)-Mb(i,2,hp,h))
             ya=(0.d0,-1.d0)*sq2*(Ma(i,1,hp,h)+Ma(i,2,hp,h))
             yb=(0.d0,-1.d0)*sq2*(Mb(i,1,hp,h)+Mb(i,2,hp,h))
             ob(2)=ob(2)+dreal(xa*CONJG(xb))/2.d0
             ob(3)=ob(3)+dreal(ya*CONJG(yb))/2.d0
          ENDDO
          ENDDO
          ENDDO
          DO hp=1,2
          DO i=1,2
          DO j=1,2
             sa=Ma(i,j,hp,1)+Ma(i,j,hp,2); sb=Mb(i,j,hp,1)+Mb(i,j,hp,2)
             ob(4)=ob(4)-(3-2*j)*dreal(sa*CONJG(sb))/4.d0
             ob(5)=ob(5)-(3-2*j)*dreal(Ma(i,j,hp,1)*CONJG(Mb(i,j,hp,1)))/2.d0
             ob(6)=ob(6)+dreal(sa*CONJG(sb))/4.d0
             ob(7)=ob(7)+dreal(Ma(i,j,hp,1)*CONJG(Mb(i,j,hp,1)))/2.d0
          ENDDO
          ENDDO
          ENDDO
          ob=ob*scl
        END SUBROUTINE bilin_obs

        SUBROUTINE set_kinematics(Elab_MeV,thdeg,labframe)
          !* sets all host kinematic variables exactly as the main loop does for
          !* kvar=3 (photon lab energy); angle either lab (default) or c.m.
          !* NOTE: uses DOUBLE precision acos (main loop uses acos(real(costh))).
          REAL*8, INTENT(in) :: Elab_MeV,thdeg
          LOGICAL, INTENT(in) :: labframe
          INTEGER i
          Tglab = Elab_MeV*1.d-3/smm
          s = 2.d0*m(1)*Tglab + msq(1)
          W = dsqrt(s)
          Tplab=(s-msq(1)-msq(3))/(2.d0*m(1)) - m(3)
          Tglab = (s-msq(1))/(2.d0*m(1))
          IF(labframe) THEN
             thLAB=PI*thdeg/180.d0
             costh=1.d0-(1.d0+2.d0*Tglab)*(1.d0-dcos(thLAB)) &
                  /(1.d0+Tglab*(1.d0-dcos(thLAB)) )
             costh=MAX(-1.d0,MIN(1.d0,costh))
             th=dacos(costh)
          ELSE
             th=PI*thdeg/180.d0
             costh=dcos(th)
             thLAB=dacos(1.d0-(1.d0-costh)/(1.d0+Tglab*(1.d0+costh)) )
          ENDIF
          thetad=thdeg
          CALL ini_polar(th)
          Eqhat=(s+1.d0)/(2.d0*W)
          wqhat=(s-1.d0)/(2.d0*W)
          qhat2=Eqhat**2-1.d0
          IF(qhat2.GE.1.d-14) THEN
             qhat=dsqrt(qhat2)
          ELSE
             qhat=0.d0
          ENDIF
          t=-2.d0*wqhat**2*(1.d0-costh)
          u=2.d0-t-s
          Ephat=(s+1.d0-pionm2)/(2.d0*W)
          wphat=(s-1.d0+pionm2)/(2.d0*W)
          esq=Ephat**2
          phat2=esq-1.d0
          sthr=(1.d0+m(3))**2
          bnd=.TRUE.
          IF(phat2.GE.1.d-14) THEN
             phat=dsqrt(phat2)
             IF(s.GT.sthr) bnd=.FALSE.
          ELSE
             phat=0.d0
          ENDIF
          vps(0)=W
          DO i=1,3
             vps(i)=0.d0
          ENDDO
          ccnst0=Ephat/W
          Tglabmev=Tglab*1.d3*smm
        END SUBROUTINE set_kinematics

        SUBROUTINE run_points()
          !* Fit driver: see driver/README.md for the file formats.
          IMPLICIT NONE
          INTEGER, PARAMETER :: NPAR=6, NQ=28
          CHARACTER(len=512) line
          CHARACTER(len=64) pid
          CHARACTER(len=8) obs, frame
          CHARACTER(len=12), DIMENSION(NPAR) :: pname
          REAL*8 Elab, thdeg, xx, rt, sclcm, scllab, omg, omgp, g3, g4, val, vdir, den
          REAL*8 vecq(0:3), vecqp(0:3), vpt(0:3)
          COMPLEX*16 Abase(8), Bk(8), Mh(2,2,2,2,0:NPAR), Mfull(2,2,2,2), Merr(2,2,2,2)
          REAL*8 ob(7), obf(7), B(7,0:NPAR,0:NPAR), cn(NQ), cd(NQ), dpar(0:NPAR)
          REAL*8 Pref(NPAR), Puser(NPAR), qv(NQ)
          REAL*8 sv_ca,sv_cb,sv_gEE,sv_gEM,sv_gME,sv_gMM,sv_g0,sv_gpi
          REAL*8 crs,crserr,dcrs, dev, maxdev(4)
          INTEGER ios, ios2, a, b2, k, npts, nbad, iq
          LOGICAL labframe, okobs, sv_cms
          pname=(/'alphaE1     ','betaM1      ','gammaE1E1   ','gammaM1M1   ', &
                  'gammaE1M2   ','gammaM1E2   '/)
          IF(hbchpt) THEN
             WRITE(*,*) 'run_points: heavy-baryon branch (hbchpt=t) has no LEC terms; not supported'
             STOP
          ENDIF
          sv_cms=cmsangle
          cmsangle=.false.
          g3=polza2LECs                   ! 1e-4 fm^3 -> internal LEC units
          g4=polza2LECs/inv_nmass_2fm     ! 1e-4 fm^4 -> internal LEC units
          IF(onlyLEX) THEN
             Pref=0.d0
          ELSE
             Pref=(/refalphaE1,refbetaM1,refgamEE,refgamMM,refgamEM,refgamME/)
          ENDIF
          Puser=(/alphaE1fm,betaM1fm,gamEEfm,gamMMfm,gamEMfm,gamMEfm/)
          dpar(0)=1.d0
          dpar(1:NPAR)=Puser-Pref
          sv_ca=calpha; sv_cb=cbeta; sv_gEE=cgamEE; sv_gEM=cgamEM
          sv_gME=cgamME; sv_gMM=cgamMM; sv_g0=cgam0; sv_gpi=cgampi
          maxdev=0.d0
          OPEN(91,file=pointsfile,status='old',iostat=ios)
          IF(ios.NE.0) THEN
             WRITE(*,*) 'run_points: cannot open points file ',TRIM(pointsfile)
             STOP
          ENDIF
          OPEN(92,file=outfile)
          WRITE(92,'(a)') '# ComptonBChPT fit-driver output (ampl.exe, &fitdriver pointsmode=t)'
          WRITE(92,'(a)') '# observable = Q_num(d)/Q_den(d),  Q(d) = sum_{a<=b} c_ab d_a d_b,  d_0=1,'
          WRITE(92,'(a)') '#   d_k = P_k - Pref_k  (k=1..6), P = (alphaE1,betaM1 [1e-4 fm^3]; gE1E1,gM1M1,gE1M2,gM1E2 [1e-4 fm^4])'
          WRITE(92,'(a)') '#   pairs (a,b) ordered (0,0),(0,1),..,(0,6),(1,1),(1,2),..,(6,6)  [28 coefficients]'
          WRITE(92,'(a,6(1x,f12.6))') '# Pref  =',Pref
          WRITE(92,'(a,6(1x,f12.6))') '# Puser =',Puser
          WRITE(92,'(a,l2,a,l2,a,l2,a,l2,a,l2)') '# flags: onlyLEX=',onlyLEX,' anomaly=',anomaly, &
               ' oneNloop=',oneNloop,' oneDloop=',oneDloop,' relchpt=',relchpt
          WRITE(92,'(a)') '# CS,SPAR,SPERP: lab-frame dsigma/dOmega_lab [nb/sr]; CSCM: c.m. [nb/sr]; S3,S2X,S2Z: asymmetries (fraction, not %)'
          WRITE(92,'(a)') '# columns: id obs E_lab[MeV] theta_lab[deg] theta_cm[deg] value(Puser) c_num(28) c_den(28)'
          IF(selftest) THEN
             OPEN(93,file=TRIM(outfile)//'.selftest')
             WRITE(93,'(a)') '# id obs value(quadratic form) value(direct full-amplitude evaluation) reldev  '// &
                  'sigma0_vs_observables()  (sperp-spar)_vs_observables()  parity(Sx,Sz)'
          ENDIF
          npts=0; nbad=0
          DO
             READ(91,'(a)',iostat=ios) line
             IF(ios.NE.0) EXIT
             line=ADJUSTL(line)
             IF(LEN_TRIM(line).EQ.0) CYCLE
             IF(line(1:1).EQ.'#'.OR.line(1:1).EQ.'!') CYCLE
             frame='LAB'
             READ(line,*,iostat=ios) pid, Elab, thdeg, obs, frame
             IF(ios.NE.0) THEN
                frame='LAB'
                READ(line,*,iostat=ios2) pid, Elab, thdeg, obs
                IF(ios2.NE.0) THEN
                   WRITE(*,*) 'run_points: cannot parse line: ',TRIM(line)
                   nbad=nbad+1
                   CYCLE
                ENDIF
             ENDIF
             CALL upcase(obs); CALL upcase(frame)
             labframe=(frame(1:2).NE.'CM')
             okobs = (obs.EQ.'CS').OR.(obs.EQ.'CSCM').OR.(obs.EQ.'S3').OR.(obs.EQ.'S2X') &
                  .OR.(obs.EQ.'S2Z').OR.(obs.EQ.'SPAR').OR.(obs.EQ.'SPERP')
             IF(.NOT.okobs) THEN
                WRITE(*,*) 'run_points: unknown observable code ',TRIM(obs),' in line: ',TRIM(line)
                nbad=nbad+1
                CYCLE
             ENDIF
             CALL set_kinematics(Elab,thdeg,labframe)
             xx=(s-1.d0)*0.5d0
             !* t exactly as in gaga_ampl (from the c.m. four-vectors)
             vecq=(/qhat,0.d0,0.d0,qhat/)
             vecqp=(/qhat,qhat*dsin(th),0.d0,qhat*dcos(th)/)
             vpt=vecq-vecqp
             rt=dot(vpt,vpt)
             !* base amplitude: all LEC shifts switched off
             calpha=0.d0; cbeta=0.d0; cgamEE=0.d0; cgamEM=0.d0
             cgamME=0.d0; cgamMM=0.d0; cgam0=0.d0; cgampi=0.d0
             Abase=(0.d0,0.d0)
             CALL setdelta(0.d0)   ! deterministic LoopTools state (LoopTools Delta)
             CALL computeAs(xx,rt,Abase)
             CALL hel_from_A(th,alphaQED*Abase,Mh(:,:,:,:,0))
             !* basis amplitudes: unit shift (1e-4 fm^3 or fm^4) of each polarisability
             DO k=1,NPAR
                Bk=(0.d0,0.d0)
                SELECT CASE(k)
                CASE(1); CALL add_LECs(xx,rt,Bk,g3,0.d0,0.d0,0.d0,0.d0,0.d0)
                CASE(2); CALL add_LECs(xx,rt,Bk,0.d0,g3,0.d0,0.d0,0.d0,0.d0)
                CASE(3); CALL add_LECs(xx,rt,Bk,0.d0,0.d0,0.d0,0.d0,0.d0,-g4)
                CASE(4); CALL add_LECs(xx,rt,Bk,0.d0,0.d0,0.d0,0.d0,g4,-g4)
                CASE(5); CALL add_LECs(xx,rt,Bk,0.d0,0.d0,g4,0.d0,0.d0,-g4)
                CASE(6); CALL add_LECs(xx,rt,Bk,0.d0,0.d0,0.d0,g4,0.d0,-g4)
                END SELECT
                CALL hel_from_A(th,alphaQED*Bk,Mh(:,:,:,:,k))
             ENDDO
             !* phase-space factors (as in SUBROUTINE observables)
             sclcm=(2.d0*Eqhat)**2/(4.d0*s) * 1.E7 *(hbarc/smm)**2
             omg=(s-1.d0)/2.d0
             omgp=(1.d0-u)/2.d0
             scllab=sclcm*(omgp/omg)**2*s
             DO a=0,NPAR
                DO b2=a,NPAR
                   CALL bilin_obs(Mh(:,:,:,:,a),Mh(:,:,:,:,b2),scllab,ob)
                   B(:,a,b2)=ob
                   B(:,b2,a)=ob
                ENDDO
             ENDDO
             !* numerator / denominator forms
             iq=0
             cd=0.d0
             DO a=0,NPAR
                DO b2=a,NPAR
                   iq=iq+1
                   SELECT CASE(TRIM(obs))
                   CASE('CS');    qv(iq)=B(1,a,b2)
                   CASE('CSCM');  qv(iq)=B(1,a,b2)*sclcm/scllab
                   CASE('SPAR');  qv(iq)=B(2,a,b2)
                   CASE('SPERP'); qv(iq)=B(3,a,b2)
                   CASE('S3');    qv(iq)=0.5d0*(B(2,a,b2)-B(3,a,b2))
                   CASE('S2X');   qv(iq)=B(4,a,b2)
                   CASE('S2Z');   qv(iq)=B(5,a,b2)
                   END SELECT
                   IF(a.NE.b2) qv(iq)=2.d0*qv(iq)
                   cn(iq)=qv(iq)
                   IF(obs(1:1).EQ.'S'.AND.obs.NE.'SPAR'.AND.obs.NE.'SPERP') THEN
                      cd(iq)=B(1,a,b2)
                      IF(a.NE.b2) cd(iq)=2.d0*cd(iq)
                   ENDIF
                ENDDO
             ENDDO
             IF(.NOT.(obs(1:1).EQ.'S'.AND.obs.NE.'SPAR'.AND.obs.NE.'SPERP')) cd(1)=1.d0
             val=qform(cn,dpar)/qform(cd,dpar)
             WRITE(92,'(a,1x,a,3(1x,f10.4),1x,es20.12,56(1x,es20.12))') TRIM(pid),TRIM(obs), &
                  Elab, thLAB*180.d0/PI, th*180.d0/PI, val, cn, cd
             npts=npts+1
             IF(selftest) THEN
                !* direct evaluation with the LECs set as the main program would set them
                calpha=dpar(1)*g3; cbeta=dpar(2)*g3
                cgamEM=dpar(5)*g4; cgamME=dpar(6)*g4; cgamMM=dpar(4)*g4
                cgamEE=dpar(3)*g4
                cgam0=-(dpar(3)+dpar(4)+dpar(5)+dpar(6))*g4
                CALL setdelta(0.d0)
                CALL gaga_ampl(qhat,th,Mfull,Merr)
                CALL bilin_obs(Mfull,Mfull,scllab,obf)
                SELECT CASE(TRIM(obs))
                CASE('CS');    vdir=obf(1)
                CASE('CSCM');  vdir=obf(1)*sclcm/scllab
                CASE('SPAR');  vdir=obf(2)
                CASE('SPERP'); vdir=obf(3)
                CASE('S3');    vdir=(obf(2)-obf(3))/(obf(2)+obf(3))
                CASE('S2X');   vdir=obf(4)/obf(6)
                CASE('S2Z');   vdir=obf(5)/obf(7)
                END SELECT
                CALL observables(Mfull,Merr,crs,crserr,dcrs,.false.)
                dev=dabs(val-vdir)/MAX(dabs(vdir),1.d-12)
                maxdev(1)=MAX(maxdev(1),dev)
                maxdev(2)=MAX(maxdev(2),dabs(crs-obf(1))/dabs(crs))
                maxdev(3)=MAX(maxdev(3),dabs(dcrs-(obf(3)-obf(2)))/MAX(dabs(dcrs),1.d-12*dabs(crs)))
                maxdev(4)=MAX(maxdev(4),dabs(obf(6)-obf(1))/obf(1),dabs(obf(7)-obf(1))/obf(1))
                WRITE(93,'(a,1x,a,2(1x,es22.14),4(1x,es10.2))') TRIM(pid),TRIM(obs),val,vdir,dev, &
                     dabs(crs-obf(1))/dabs(crs), dabs(dcrs-(obf(3)-obf(2)))/MAX(dabs(dcrs),1.d-12*dabs(crs)), &
                     MAX(dabs(obf(6)-obf(1)),dabs(obf(7)-obf(1)))/obf(1)
             ENDIF
          ENDDO
          CLOSE(91); CLOSE(92)
          calpha=sv_ca; cbeta=sv_cb; cgamEE=sv_gEE; cgamEM=sv_gEM
          cgamME=sv_gME; cgamMM=sv_gMM; cgam0=sv_g0; cgampi=sv_gpi
          cmsangle=sv_cms
          WRITE(*,'(a,i8,a,a)') ' run_points: processed ',npts,' points -> ',TRIM(outfile)
          IF(nbad.GT.0) WRITE(*,'(a,i6,a)') ' run_points: WARNING ',nbad,' input lines skipped'
          IF(selftest) THEN
             CLOSE(93)
             WRITE(*,'(a,es10.2)') ' selftest: max rel. dev. quadratic form vs direct evaluation   =',maxdev(1)
             WRITE(*,'(a,es10.2)') ' selftest: max rel. dev. sigma0 vs SUBROUTINE observables      =',maxdev(2)
             WRITE(*,'(a,es10.2)') ' selftest: max rel. dev. (sperp-spar) vs SUBROUTINE observables=',maxdev(3)
             WRITE(*,'(a,es10.2)') ' selftest: max parity check |S_x+ - s0|,|S_z+ - s0| / s0        =',maxdev(4)
          ENDIF
        END SUBROUTINE run_points

        SUBROUTINE run_helpoints()
          !* [PWA copy, 2026-10] dump c.m. helicity amplitudes mgg(i,j,hp,h) per mechanism
          !* (1 Born N s+u, 2 pi0 anomaly, 3 piN loops, 4 Delta pole, 5 piDelta loops, 6 LEC/contact)
          !* input lines: id E_lab[MeV] theta[deg] frame(LAB|CM)
          !* output: id E_lab th_lab th_cm s Eqhat qhat  then for c=1..6, (i,j,hp,h) loops (h outer):
          !*         Re,Im of mgg  (mgg already includes alphaQED; nucleon-mass units; spinors u^+u=1)
          IMPLICIT NONE
          CHARACTER(len=512) line
          CHARACTER(len=64) pid
          CHARACTER(len=8) frame
          REAL*8 Elab, thdeg, xx, rt, vecq(0:3), vecqp(0:3), vpt(0:3)
          COMPLEX*16 Ac(8,6), Mh(2,2,2,2,6)
          INTEGER ios, c, i, j, hp, h, npts
          LOGICAL labframe
          IF(hbchpt) THEN
             WRITE(*,*) 'run_helpoints: hbchpt not supported'; STOP
          ENDIF
          cmsangle=.false.
          OPEN(91,file=pointsfile,status='old',iostat=ios)
          IF(ios.NE.0) THEN
             WRITE(*,*) 'run_helpoints: cannot open ',TRIM(pointsfile); STOP
          ENDIF
          OPEN(92,file=outfile)
          WRITE(92,'(a)') '# helicity amplitudes per component; see run_helpoints in ampl.f90'
          npts=0
          DO
             READ(91,'(a)',iostat=ios) line
             IF(ios.NE.0) EXIT
             line=ADJUSTL(line)
             IF(LEN_TRIM(line).EQ.0) CYCLE
             IF(line(1:1).EQ.'#') CYCLE
             READ(line,*,iostat=ios) pid, Elab, thdeg, frame
             IF(ios.NE.0) CYCLE
             CALL upcase(frame)
             labframe=(frame(1:2).NE.'CM')
             CALL set_kinematics(Elab,thdeg,labframe)
             xx=(s-1.d0)*0.5d0
             vecq=(/qhat,0.d0,0.d0,qhat/)
             vecqp=(/qhat,qhat*dsin(th),0.d0,qhat*dcos(th)/)
             vpt=vecq-vecqp
             rt=dot(vpt,vpt)
             CALL setdelta(0.d0)
             CALL components_As(xx,rt,Ac)
             DO c=1,6
                CALL hel_from_A(th,alphaQED*Ac(:,c),Mh(:,:,:,:,c))
             ENDDO
             WRITE(92,'(a,1x,f12.6,2(1x,f14.9),3(1x,es24.16),192(1x,es24.16))') TRIM(pid), Elab, &
                  thLAB*180.d0/PI, th*180.d0/PI, s, Eqhat, qhat, &
                  ((((( dreal(Mh(i,j,hp,h,c)), dimag(Mh(i,j,hp,h,c)) ,i=1,2),j=1,2),hp=1,2),h=1,2),c=1,6)
             npts=npts+1
          ENDDO
          CLOSE(91); CLOSE(92)
          WRITE(*,'(a,i8,a,a)') ' run_helpoints: processed ',npts,' points -> ',TRIM(outfile)
        END SUBROUTINE run_helpoints

        SUBROUTINE nnlo_static_pols(Pout, Pcomp)
          !* Analytic NNLO static dipole polarisabilities of the proton (LMP 2015, App. C),
          !* evaluated with the constants actually used by the code: gA, fpi, pionm (charged),
          !* smm, alphaQED, hA, gM=dng1, gE=dng2; Delta mass: bary(3)%rm for the pole and the
          !* value hard-coded in CorrIsoCovOneDeltaLoop (deltam=1.3130548775939 M) for the loops.
          !* Output in 1e-4 fm^3 (alphaE1,betaM1) and 1e-4 fm^4 (gE1E1,gM1M1,gE1M2,gM1E2).
          !* Pcomp(:,1..3) = piN loops, Delta pole, piDelta loops.
          REAL*8, INTENT(out) :: Pout(6), Pcomp(6,3)
          INTEGER, PARAMETER :: NG=400
          REAL*8 xg(NG), wg(NG)
          REAL*8 e2,mu,sq,ac,L,cN,cg,aN,bN,g1,g2,g3,g4,d,cD,Mloop,Mp,Dl,kk,f3,f4,Fp
          REAL*8 x,DD,lg,Ia,Ib,I1,I2,I3,I4
          INTEGER i
          e2=4.d0*PI*alphaQED
          Fp=fpi                       ! already in nucleon-mass units
          mu=pionm
          sq=dsqrt(4.d0-mu**2); ac=dacos(mu**2/2.d0-1.d0); L=dlog(mu**2)
          cN=e2*gA**2/(192.d0*PI**3*Fp**2)
          aN=cN/(mu*sq**5)*(mu*(-18*mu**6+157*mu**4-407*mu**2+304)*sq + mu*(9*mu**4-20*mu**2+9)*sq**5*L &
               -(9*mu**10-110*mu**8+479*mu**6-870*mu**4+590*mu**2-80)*ac)
          bN=-cN/(mu*sq**3)*(mu*(-54*mu**4+235*mu**2-127)*sq - mu*(27*mu**4-50*mu**2+9)*sq**3*L &
               -(27*mu**8-212*mu**6+471*mu**4-246*mu**2+2)*ac)
          cg=cN
          g1=cg/(mu**2*sq**5)*(sq*(18*mu**8-173*mu**6+524*mu**4-452*mu**2-(4-mu**2)**2*(9*mu**4-28*mu**2+13)*mu**2*L+32) &
               + mu*(9*mu**10-118*mu**8+563*mu**6-1150*mu**4+848*mu**2-104)*ac)
          g2=2*cg/(mu**2*sq**5)*(2*sq*(18*mu**8-149*mu**6+348*mu**4-179*mu**2-(4-mu**2)**2*(9*mu**4-16*mu**2+4)*mu**2*L+4) &
               + mu*(18*mu**10-212*mu**8+868*mu**6-1400*mu**4+701*mu**2-44)*ac)
          g3=-cg/2/(mu**2*sq**5)*(sq*(18*mu**8-141*mu**6+284*mu**4-4*mu**2-3*(mu**2-4)**2*(3*mu**4-4*mu**2-1)*mu**2*L-16) &
               + mu*(9*mu**10-102*mu**8+387*mu**6-510*mu**4+46*mu**2+56)*ac)
          g4=cg/2/(mu**2*sq**5)*(sq*(-126*mu**8+1047*mu**6-2462*mu**4+1308*mu**2+(mu**2-4)**2*(63*mu**4-114*mu**2+29)*mu**2*L-16) &
               + (-63*mu**11+744*mu**9-3059*mu**7+4970*mu**5-2534*mu**3+152*mu)*ac)
          Pcomp(:,1)=(/aN,bN,-g1-g3,g4,g3,g2+g4/)
          !* Delta pole (gM=dng1, gE=dng2)
          Mp=1.d0+bary(3)%rm; Dl=bary(3)%rm-1.d0
          kk=alphaQED/Mp**2
          g1=kk*(dng1**2/Dl-2*dng2*dng1/Mp+dng2**2*bary(3)%rm/Mp**2)
          g2=-kk*(dng1**2*bary(3)%rm/Dl**2+dng2*dng1/Dl)
          g3=-kk*(dng1**2/(2*Dl)-3*dng2*dng1/(2*Mp))
          g4=kk*(dng1**2*bary(3)%rm/Dl**2-dng2*dng1/(2*Dl)+dng2**2/(2*Mp))
          Pcomp(:,2)=(/-alphaQED*2*dng2**2/Mp**3, alphaQED*2*dng1**2/(Mp**2*Dl),-g1-g3,g4,g3,g2+g4/)
          !* piDelta loops (MSbar, Xi=0), Delta mass as hard-coded in the loop routine
          Mloop=1.3130548775939d0
          d=Mloop-1.d0
          cD=e2*hA**2/(1728.d0*PI**3*Fp**2*Mloop**2)
          CALL gauleg(0.d0,1.d0,xg,wg,NG)
          Ia=0.d0; Ib=0.d0; I1=0.d0; I2=0.d0; I3=0.d0; I4=0.d0
          DO i=1,NG
             x=xg(i)
             DD=x**2+(1+d)**2-x*(2+2*d+d**2-mu**2)
             lg=-dlog(DD)
             Ia=Ia+wg(i)*(18*(x-1)**2*x**4*(x+d+1)/DD**2 &
                  -x**2*(120*x**3+4*(27*d-16)*x**2-(166*d+115)*x+57*(d+1))/DD &
                  +6*x*(x*(22*x+15*d-1)-12*(d+1))*lg)
             Ib=Ib+wg(i)*((24*x**2-32*x+9)*x**2*(x+d+1)/DD - 6*(5*x-4)*x*(2*x+3*d+3)*lg)
             I1=I1+wg(i)*(8*(x-1)**3*x**5*(x+d+1)/DD**3 - 2*(x-1)**2*x**3*(22*x**2+8*d*x+3*d+3)/DD**2 &
                  +(x-1)*x**2*(105*x**2-2*(5*d+38)*x+33*(d+1))/DD + 2*x*(x*(-16*x+9*d+21)-12*(d+1))*lg)
             I2=I2+wg(i)*(2*(4*x**2-7*x+3)*x**3*(x+d+1)/DD**2 &
                  -x**2*(21*x**3-4*(3*d+5)*x**2+(37*d+26)*x-24*(d+1))/DD + 12*(3*x-2)*x**2*lg)
             I3=I3+wg(i)*((x-1)*x**5*(x+d+1)/DD**2 - x**2*(24*x**3+(2*d-11)*x**2+(8*d-1)*x-9*(d+1))/(2*DD) &
                  + x*(x*(34*x-9*d-33)+12*(d+1))*lg)
             I4=I4+wg(i)*(-(x**3+5*x**2-12*x+6)*x**3*(x+d+1)/DD**2 &
                  + x**2*(48*x**3-(22*d+59)*x**2+(80*d+71)*x-57*(d+1))/(2*DD) + x*(x*(-34*x+9*d+33)-12*(d+1))*lg)
          ENDDO
          g1=cD*(2*d+7.d0/4+I1); g2=cD*(5.d0/12+I2); g3=cD*(-d-2.d0/3+I3); g4=cD*(d+2.d0/3+I4)
          Pcomp(:,3)=(/cD*(8*d+35.d0/3+Ia), cD*(-8*d-65.d0/6+Ib),-g1-g3,g4,g3,g2+g4/)
          f3=inv_nmass_2fm**3*1.d4; f4=inv_nmass_2fm**4*1.d4
          DO i=1,3
             Pcomp(1:2,i)=Pcomp(1:2,i)*f3
             Pcomp(3:6,i)=Pcomp(3:6,i)*f4
          ENDDO
          Pout=Pcomp(:,1)+Pcomp(:,2)+Pcomp(:,3)
        END SUBROUTINE nnlo_static_pols

        REAL*8 FUNCTION qform(c,d)
          REAL*8, INTENT(in) :: c(28), d(0:6)
          INTEGER a,b2,iq
          qform=0.d0
          iq=0
          DO a=0,6
             DO b2=a,6
                iq=iq+1
                qform=qform+c(iq)*d(a)*d(b2)
             ENDDO
          ENDDO
        END FUNCTION qform

        SUBROUTINE upcase(str)
          CHARACTER(len=*), INTENT(inout) :: str
          INTEGER i,ic
          DO i=1,LEN(str)
             ic=ICHAR(str(i:i))
             IF(ic.GE.97.AND.ic.LE.122) str(i:i)=CHAR(ic-32)
          ENDDO
        END SUBROUTINE upcase

        SUBROUTINE run_static()
          !* Dumps Re/Im A_i(nu,t) by mechanism on a grid of small nu (=(s-u)/4M) and t,
          !* in nucleon-mass units WITHOUT alpha_em, for the static-polarisability
          !* extraction (driver/static_pols.py). Output: static_As.out
          IMPLICIT NONE
          INTEGER, PARAMETER :: NNU=14, NT=6
          REAL*8 nug(NNU), tf(NT), nuv, tv, xx
          COMPLEX*16 Ac(8,6)
          INTEGER inu, itf, c, i
          nug=(/0.003d0,0.005d0,0.0075d0,0.01d0,0.0125d0,0.015d0,0.02d0,0.025d0, &
                0.03d0,0.04d0,0.05d0,0.06d0,0.07d0,0.08d0/)
          tf=(/0.d0,-0.25d0,-0.5d0,-1.d0,-2.d0,-3.d0/)   ! t = tf * nu^2
          OPEN(95,file='static_As.out')
          WRITE(95,'(a)') '# nu t  then for comp=1..6 (Born,anomaly,piN,Dpole,piD,LEC): ReA1..ReA8 ImA1..ImA8'
          WRITE(95,'(a,6(1x,es16.8))') '# constants (M units): mpi,fpi,gA,hA,Md,pi0mass', pionm,fpi,gA,hA,bary(3)%rm,pi0mass
          DO itf=1,NT
             DO inu=1,NNU
                nuv=nug(inu)
                tv=tf(itf)*nuv**2
                xx=nuv-tv/4.d0
                CALL components_As(xx,tv,Ac)
                WRITE(95,'(2(1x,es24.16),96(1x,es24.16))') nuv,tv, &
                     ((dreal(Ac(i,c)),i=1,8),(dimag(Ac(i,c)),i=1,8),c=1,6)
             ENDDO
          ENDDO
          CLOSE(95)
          WRITE(*,*) 'run_static: wrote static_As.out'
          CALL diag_uv()
        END SUBROUTINE run_static

        SUBROUTINE diag_uv()
          !* Diagnostic: dependence of the (renormalised, supposedly UV-finite) piN loops
          !* on the LoopTools UV parameter Delta, which CorrIsoCovOneDeltaLoop leaves at 1.
          IMPLICIT NONE
          REAL*8 Eg(5), xx, tv, ww, dl(3)
          COMPLEX*16 A1(8), A2(8), A3(8)
          INTEGER ie, id
          REAL*8 r12, r13
          Eg=(/0.03d0,0.08d0,0.14d0,0.2d0,0.3d0/)
          dl=(/0.d0,1.d0,1.d3/)
          WRITE(*,'(a)') ' diag_uv: E_lab[GeV]  max_i|A_i(D=1)-A_i(D=0)|/max|A|   max_i|A_i(D=1e3)-A_i(D=0)|/max|A|'
          DO ie=1,5
             ww=Eg(ie)/smm
             xx=ww
             tv=-2.d0*ww**2/(1.d0+2.d0*ww)   ! roughly 90 deg c.m.
             CALL setdelta(dl(1)); CALL CovOneLoop(1.d0+2.d0*xx,tv,pionm**2,gA**2/(4.d0*Pi*fpi)**2,A1)
             CALL setdelta(dl(2)); CALL CovOneLoop(1.d0+2.d0*xx,tv,pionm**2,gA**2/(4.d0*Pi*fpi)**2,A2)
             CALL setdelta(dl(3)); CALL CovOneLoop(1.d0+2.d0*xx,tv,pionm**2,gA**2/(4.d0*Pi*fpi)**2,A3)
             r12=MAXVAL(ABS(A2-A1))/MAXVAL(ABS(A1))
             r13=MAXVAL(ABS(A3-A1))/MAXVAL(ABS(A1))
             WRITE(*,'(a,f8.3,2(3x,es10.2))') ' diag_uv:',Eg(ie),r12,r13
          ENDDO
          CALL setdelta(0.d0)
          !* reproducibility of repeated evaluations (LoopTools cache / state)
          CALL diag_repro()
        END SUBROUTINE diag_uv

        SUBROUTINE diag_repro()
          IMPLICIT NONE
          REAL*8 Eg(4), thl(4), xx, tv, ww, cth, d(6)
          COMPLEX*16 Ac1(8,6), Ac2(8,6), Ac3(8,6)
          INTEGER ie, c
          Eg=(/0.06d0,0.14d0,0.20d0,0.30d0/)
          thl=(/90.d0,140.d0,60.d0,120.d0/)
          WRITE(*,'(a)') ' diag_repro: E[GeV] th  max rel diff 2nd/1st call per component (Born,anom,piN,Dpole,piD,LEC); 3rd after clearcache'
          DO ie=1,4
             ww=Eg(ie)/smm
             xx=ww
             cth=1.d0-(1.d0+2.d0*ww)*(1.d0-dcos(thl(ie)*PI/180.d0))/(1.d0+ww*(1.d0-dcos(thl(ie)*PI/180.d0)))
             tv=-2.d0*(ww/dsqrt(1.d0+2.d0*ww))**2*(1.d0-cth)
             CALL setdelta(0.d0); CALL components_As(xx,tv,Ac1)
             CALL setdelta(0.d0); CALL components_As(xx,tv,Ac2)
             CALL clearcache
             CALL setdelta(0.d0); CALL components_As(xx,tv,Ac3)
             DO c=1,6
                d(c)=MAXVAL(ABS(Ac2(:,c)-Ac1(:,c)))/MAX(MAXVAL(ABS(Ac1(:,c))),1.d-300)
             ENDDO
             WRITE(*,'(a,f6.3,f6.0,6es10.2)') ' diag_repro 2nd:',Eg(ie),thl(ie),d
             DO c=1,6
                d(c)=MAXVAL(ABS(Ac3(:,c)-Ac1(:,c)))/MAX(MAXVAL(ABS(Ac1(:,c))),1.d-300)
             ENDDO
             WRITE(*,'(a,f6.3,f6.0,6es10.2)') ' diag_repro 3rd:',Eg(ie),thl(ie),d
          ENDDO
        END SUBROUTINE diag_repro

       SUBROUTINE piDloop(w,t,A)
          REAL*8, INTENT(in):: w,t
          COMPLEX*16, INTENT(inout):: A(8)
          REAL*8 costh,const,ga5,a1,a2,a3,a4,a5,a6
          Real*8 alP,beP,gaP(4),lnDloop,hA

          hA=2.7 ! pi-N-Delta coupling in the large Nc limit
          lnDloop=dlog(Del/pionm+dsqrt( (Del/pionm)**2-1.d0) )
          alP=(2./3.)/(Del**2-pionm2)*(Del+(1./9.)*lnDloop &
               * (Del**2 - 10.d0*pionm2)/dsqrt(Del**2-pionm2) ) &
               *(hA/(4.*Pi*fpi))**2
          beP=(2./3.)/dsqrt(Del**2-pionm2) *(1./9.)*lnDloop &
             *(hA/(4.*Pi*fpi))**2
          gaP(1) =- 2.*(hA/(4.*Pi*fpi))**2/(27.*Del**2)
          gaP(2) =2.*(1.d0-lnDloop)*(hA/(4.*Pi*fpi))**2/(27.*Del**2) 
          gaP(3) = (1.d0-lnDloop)*(hA/(4.*Pi*fpi))**2/(27.*Del**2) 
          gaP(4)=- gaP(3) 
     
          costh = t/(2.d0*w**2)+1.d0
          ga5=-(gaP(2)+2.d0*gaP(4))
          a1=(alP+beP)*w**2+0.5d0*t*beP
          a2= - beP
          a3= w**3*(gaP(1)+ga5*costh)
          a4= w*gaP(2)
          a5= w*gaP(4)
          a6= w*gaP(3)

          !const=4.d0*Pi
          const=1.d0
          A(1)=A(1)+const*a1
          A(2)=A(2)+const*a2
          A(3)=A(3)+const*a3      
          A(4)=A(4)+const*a4      
          A(5)=A(5)+const*a5
          A(6)=A(6)+const*a6

        END SUBROUTINE piDloop
		
        SUBROUTINE LEXobservables(OUTPUT)
          IMPLICIT NONE
		  LOGICAL OUTPUT
		  REAL*8 SI_PAR,SI_PERP,RTOT,RTOT_ERR,RPOL_ERR,RBETA_ERR,SI_X,RTOTB,totxcn,polxcn
          REAL*8 SCALE,omegap,omega,w,wp,Q,kap,zLAB,omegaB,zB,A8,A8p,mpi
		  REAL*8 basym,tasym,gasym,gasym_err,gtotal,Sigma3Born,Sigma3LEX2,Sigma3LEX4,SI_PERP_LEX2,SI_PAR_LEX2
		  REAL*8 Sigma2zBorn,Sigma2xBorn,Sigma2zLEX3,Sigma2xLEX3
		  REAL*8 RB(0:4)
		  
            !compute first the Klein-Nishina, Powell and Maximon cross-sections
            if(proton) then
             Q=1.d0 ! 1-proton, 0-neutron
             kap=kapgp
            else
             Q=0.d0
             kap=kapgn
            endif
			if(.not.nuclamm) kap=0.
            omega=(s-1.d0)/2.d0
	        omegap=(1.d0-u)/2.d0
				w=omega
				wp=omegap
			!SCALE=(alphaQED*hbarc/smm)**2*1.E7 *0.25d0*(omegap/omega)**2
			!          Phase-space factor
            SCALE=1.d0/(4.d0*s) * 1.E7 *(alphaQED*hbarc/smm)**2
          IF(.NOT.cmsangle) THEN
              SCALE=SCALE*(omegap/omega)**2*s
          ENDIF
			zLAB=dcos(thLAB)
	     ! Klein-Nishina unpolarized and polarized
		!!   RTOT=2.d0*SCALE*Q**4*(omegap/omega+omega/omegap-dsin(thLAB)**2)
	    !!   SI_PERP=2.d0*SCALE*Q**4*(omegap/omega+omega/omegap)
		!!   SI_PAR=2.d0*SCALE*Q**4*(omegap/omega+omega/omegap-2.d0*dsin(thLAB)**2)
		   omegaB=(s-1.d0+t/2.d0)/dsqrt(4.d0-t)
		   zB=1.d0+t/(2.d0*omegaB**2)
		   RTOTB=0.
			!if(.not.varyscalarPs) write(35,103) thLAB*180./PI, RTOT,SI_PERP-SI_PAR
	     ! Powell unpolarized and polarized (Maximon)
        !!         RTOT=2.d0*SCALE*(Q**4*(1.d0+zLAB**2) + omegap*omega &
        !!         *( Q**4 *(1.d0-zLAB)**2+kap*(2.*Q**3+9./2.*kap*Q**2+3.*kap**2*Q+3./4.*kap**3) &
        !!         - kap*Q*(4.*Q**2+5.*kap*Q+2.*kap**2)*zLAB &
        !!         + kap*(2.*Q**3+0.5*kap*Q**2 -kap**2*Q - 0.25*kap**3)*zLAB**2) )
	    !!   SI_PERP=SI_PERP+2.d0*SCALE*omegap*omega*(2.*kap*Q**3*(1.d0-zLAB)**2  &
        !!                     + kap**2*Q**2 *(5.*(1.d0-zLAB)+0.5*(1.d0-zLAB**2)) &
        !!                     + kap**3 *Q*(2.*(1.d0-zLAB)+2.*(1.d0-zLAB**2)) + kap**4*(2.d0-zLAB**2)/2.d0 )
	    !!   SI_PAR=SI_PAR+2.d0*SCALE*omegap*omega*(2.*kap*Q**3*(1.d0-zLAB)**2  &
        !!                     + kap**2*Q**2 *(5.*(1.d0-zLAB) - 1.5*(1.d0-zLAB**2)) &
        !!                     + kap**3 *Q*2.*(1.d0-zLAB) + kap**4/2.d0 )
		   !Sigma3Born=(-SI_PERP+SI_PAR)/(SI_PERP+SI_PAR)*100.d0
		   
		   ! Born + anomaly contribution to all observables (our calculation)
		   RB(0)=1.- (w-wp)**2 *(kap**2-1.)*(1.+2.*kap+0.5*kap**2)/(2*w*wp) + (w-wp)*3.*kap**2 *(1.+kap+kap**2/6.) + w*wp*kap**4/4.
		   RB(1)=-1.- w*wp*kap**2 *(1.+kap+kap**2/4.)
		   RB(2)=(w/wp-wp/w)*(0.5+kap) + (w/wp-w-1.)*kap**2 + (w-wp)**2 *kap**2 *(1.+4.*kap+kap**2)/(4*w*wp) &
		         -(w-wp)*kap**2 *(1.+2.*kap+0.5*kap**2)
		   RB(3)=0.5*(w/wp-w-1.)*(1.+kap) - 0.5*wp*(1.+kap)**2 - (w-wp)*kap**2 *(3.+4.*kap+kap**2)/4. &
				+ w*wp*kap**2*(1.+3.*kap+kap**2)/4.
		   RB(4)=0.5*(w+wp)*(1.+kap) + 0.5*w*kap + (w-wp)*kap**2 *(3.+kap)/4. - w*wp*kap**2*(1.+kap)/4.
		   A8=0.d0
		   if(anomaly) then ! add anomaly contribution
			A8 = - 1/(pi0mass**2-t)*(2.d0*Q - 1.d0)*gA/(2.d0*Pi*fpi)**2
			RB(0) = RB(0) + A8*(w-wp)**2 *( (1.d0-zLAB)*(1.d0+kap) + kap*(2.d0+kap) )/2.d0 + A8**2 *(w-wp)**3 /2.d0
			RB(1) = RB(1)
			RB(2) = RB(2) - A8*(1-zLAB)*( -(w+wp)*(w-wp)+ kap**2 *(-2.d0*w*(w-wp) +4.d0*w**2 *wp) &
							+ kap*((w+wp)**2 + 2.d0*w*wp*(w+wp) + 4.d0*w**2 *(-1.d0+wp)) )/2.d0 &
							+ (A8*w*wp)**2 *(1.d0-zLAB**2)
			RB(3) = RB(3) + A8*( -2.d0*wp*(w-wp) + 2.d0*kap*((w+wp)**2 + 2.d0*w*wp*(w+wp) - 4.d0*w**2) &
					+ kap**2 *w*(-4.d0*(w-wp) + wp*(w+wp) + 2.d0*w*wp) )/4.d0 &
					- (A8*w*wp)**2 *zLAB
			RB(4) = RB(4) + A8*( wp*(w-wp)/2.d0 - w*wp*(w+wp)*kap - w*wp*(3.d0*w+wp)*kap**2 /4.d0 ) + (A8*w*wp)**2
		   endif	   
		   SI_PERP = 4.d0*SCALE*RB(0)
		   SI_PAR = 4.d0*SCALE*(RB(0)+(1-zLAB**2)*RB(1))
		   RTOT=0.5d0*(SI_PAR+SI_PERP)
		   Sigma2zBorn = -2.*100.*(-RB(2) + (1-zLAB**2)*RB(4))/(2.*RB(0)+(1-zLAB**2)*RB(1))
		   Sigma2xBorn = 2.*100.*dsin(thLAB)*(RB(3) + zLAB*RB(4))/(2.*RB(0)+(1-zLAB**2)*RB(1))
		   Sigma3Born = 100.*(1-zLAB**2)*RB(1)/(2.*RB(0)+(1-zLAB**2)*RB(1))
		   !Sigma3Born= -100.d0*(1.-zB**2)*(16. + 4.*(kap**2*(2.+kap)**2-2.*(1.-zB)*zB)*omegaB**2 &
			!		-(1.-zB)*(zB**2-zB*(1.+kap)**2-kap*(2.+kap))*(zB+kap*(2.+kap))*omegaB**4) &
		!			/(16.*(1.+zB**2)+4.*(kap**4*(3.-zB**2)+2.*(1.-zB)**2*(2.+zB+zB**2)+8.*kap*(1.-zB)**2 &
		!			+2.*kap**2*(1.-zB)*(9.-zB)+4.*kap**3*(1.-zB)*(3.+zB) )*omegaB**2 &
		!			+(1.-zB)**2*( zB**4+zB**2*(1.+kap)**4+2.*kap*zB*(1.+kap)**2*(2.+kap)+kap**2*(2.+kap)**2 )*omegaB**4 )

		!Sigma2zBorn= -2.*dsqrt(2.d0*(1.d0-zB))*( 4.*(1.-zB)*((1.+kap)*(1.-zB)+2.*kap +kap**2)*omegaB &
		!			-2.*(1.+zB)*(1.+(1.+kap)**2)*((1.+kap)*(1.-zB)+kap**2)*omegaB**2 &
		!			+(1.-zB)*((1.+kap)*(1.-zB)+2.*kap +kap**2)*(zB**2-zB*(1.+kap)**2-kap*(2. +kap))*omegaB**3) &
		!			/(16.*(1.+zB**2)+8.*kap**4*omegaB**2 + (1.-zB)*omegaB**2*(4.*kap**4*(1.+zB)+8.*(1.-zB)*(2.+zB+zB**2) &
		!			+32.*(1.-zB)*kap+8.*(9.-zB)*kap**2+16.*(3.+zB)*kap**3 + (1.-zB)*omegaB**2*(zB**4+zB**2*(1.+kap)**4 &
		!			+2.*zB*kap*(1.+kap)**2*(2.+kap) +kap**2*(2.+kap)**2) ))
			write(*,*) "Born results"
			write(*,*) 'Sigma2zBorn =', Sigma2zBorn
			write(*,*) 'Sigma2xBorn =', Sigma2xBorn!,-2.*dsin(thLAB)*(RB(3) + zLAB*RB(4))/(2.*RB(0)+(1-zLAB**2)*RB(1)),dsin(thLAB)
			write(*,*) 'Sigma3Born =',Sigma3Born!,(-SI_PERP+SI_PAR)/(SI_PERP+SI_PAR)

			SI_PERP_LEX2=-SCALE*8.d0*omega*omegap*(alphaE1+betaM1*zLAB)
			SI_PAR_LEX2=-SCALE*8.d0*omega*omegap*(alphaE1*zLAB**2+betaM1*zLAB)
		   ! Sigma3LEX=Sigma3Born- 4.*(betaM1)*omegap*omega*zLAB*(1.d0-zLAB**2)/(1.d0+zLAB**2)**2
		   Sigma3LEX2=0.d0
		   Sigma3LEX4=0.d0
			if(Q.gt.0.d0) then
				Sigma3LEX2=- 4.*(betaM1)*omegaB**2*zB*(1.d0-zB**2)/(1.d0+zB**2)**2/(Q**2)
				Sigma3LEX2=Sigma3LEX2*100.d0
				
				Sigma2zLEX3=  -omegaB**3 *( &
							-alphaE1*(1.d0+zB**2)*(kap - 2.d0*zB*(1.d0+kap)**2 + zB**2 *(2.d0+3.d0*kap)) &
							+betaM1*(-(1.d0+kap)**2 + (1.d0-kap)*zB + 6.d0*(1.d0+kap)**2*zB**2 + (-6.d0-8.d0*kap)*zB**3 &
							-(1.d0+kap)**2*zB**4 + (1.d0+kap)*zB**5) &
							+2.d0*(1.d0+zB**2)*( (1.d0+zB**2)*(gamEE+zB*gamME)+2.*zB*(gamMM+zB*gamEM) ) )/(zB**2+1.d0)**2
				Sigma2zLEX3=Sigma2zLEX3*100.d0
				
				Sigma2xLEX3= -dsqrt(1.d0-zB**2)*omegaB**3 *( &
							-alphaE1*(1.d0+zB**2)*(zB+ 2.d0*zB*kap - (1.d0+kap)**2) &
							+betaM1*(kap + 3.d0*(1.d0+kap)**2*zB + (-3.d0-6.d0*kap)*zB**2 - (1.d0+kap)**2*zB**3 + (1.d0+kap)*zB**4) &
							+2.d0*(1.d0+zB**2)*( gamMM+zB*(gamEE+gamEM)+zB**2*gamME)  )/(zB**2+1.d0)**2
				Sigma2xLEX3=Sigma2xLEX3*100.d0


				!Sigma2zLEX3=  -omegaB**3 *( &
				!			-alphaE1*(kap - 2.d0*zB*(1.d0+kap)**2 + zB**2 *(2.d0+3.d0*kap))/(1.d0+zB**2) &
				!			+betaM1*zB*(2.d0 + 2.d0*zB**2 + (1.d0+kap)*(-1.d0-8.d0*zB**2+zB**4) + (1.d0+kap)**2 *zB*(5.d0-2.d0*zB**2+zB**4) )/(1.d0+zB**2) &
				!			+2.*( (1.+zB**2)*(gamEE+zB*gamME)+2.*zB*(gamMM+zB*gamEM) ) )/(zB**2+1.) &
				!			- omegaB**3*betaM1*(1.+kap)**2*(1.-zB**2)**2/(1.+zB**2)
				!Sigma2xLEX3= -dsqrt(1.d0-zB**2)*omegaB**3 *( &
				!			+alphaE1*2.d0*( (1.d0+kap)**2-(2.d0*kap+1.)*zB) &
				!			+betaM1*2.d0*(-kap**2 *zB*(-3.d0+zB**2) + zB*(3.d0-3.d0*zB-zB**2+zB**3) &
				!			+ kap*(1.d0+6.d0*zB-6.d0*zB**2-2.d0*zB**3+zB**4))/(zB**2+1.d0) &
				!			+2.*( gamMM+zB*(gamEE+gamEM)+zB**2*gamME)  )/(zB**2+1.)
				!Sigma2xLEX3= dsqrt(1.d0-zB**2)*omegaB**3 *(0.*((kap+1.)*(1.-zB)+2.*kap+kap**2)*(alphaE1+betaM1+4.*zB*betaM1/(1.+zB**2)) &
				!			+2.*( gamEE+zB*(gamMM+gamME)+zB**2*gamEM)  )/(zB**2+1.)
				!Sigma2zLEX3= omegaB**3 *(2.*alphaE1*(1.-zB)*zB*(1.+zb**2) +betaM1*(1.-zB-6*zB**2)*0.) /((zB**2+1.)**2)

			!gamEE,gamME
			
				Sigma3LEX4=	omegaB**4*((-1 + zB**2)*(alphaE1*(1 + zB**2)*(4*betaM1*Q**2*zB  &
              + (2*Q**4*(-1 + zB) + kap**4*(-1 + zB**2) - 2*kap*Q**3*(2.d0 - 4*zB + zB**2) + 4*kap**3*Q*(-1.d0 + zB + zB**2) &
			  + kap**2*Q**2*(-7.d0 + 10*zB + 2*zB**2))) + betaM1*(-2*betaM1*Q**2*(1 - 6*zB**2 + zB**4) &
			  + zB*(kap**4*(-5.d0 + 3*zB**2) + 4*kap**3*Q*(-5.d0 + 4*zB + 3*zB**2) - Q**4*(7.d0 - 8*zB + zB**4) &
			  + 2*kap*Q**3*(-8.d0 + 16*zB - 7*zB**2 + zB**4) + kap**2*Q**2*(-32.d0 + 40.*zB + zB**2 + zB**4))))) &
			  /(Q**6*(1.d0 + zB**2)**3)	
			  Sigma3LEX4=Sigma3LEX4*100.d0
			endif
			
			if (pwa) then 
			 totxcn=2*Pi*1.E4 *(alphaQED*hbarc/smm)**2*(-(betaM1*((w*(8*w**2*(1 + 2*w) &
				+ (1 + kap)**2*(1 + 5*w + (22*w**2)/3. + 2*w**3)))/(2.*(1 + 2*w)**3) + ((1 + kap)**2*Log(1/(1 + 2*w)))/4.)) &
				+ alphaE1*((2 + 10*w + 12*w**2 - 4*w**3 + (-9 + (2*kap)/3. + kap**2/3.)*w**4)/(1 + 2*w)**3 + Log(1/(1 + 2*w))/w) &
				+ ((2*w*(8 + 32*w + (36 - 8*kap - 2*kap**2 + 4*kap**3 + kap**4)*w**2 + (4 - 24*kap - 6*kap**2 + 12*kap**3 + 3*kap**4)*w**3 &
				+ 4*kap**2*(2 + kap)**2*w**4 + 2*kap**4*w**5))/(1 + 2*w)**2 + (8 + 8*w + (-4 - 8*kap - 2*kap**2 + 4*kap**3 + kap**4)*w**2) &
				*Log(1/(1 + 2*w)))/(8.*w**3) )
				if(anomaly) then
				A8p =  (gA*(-1 + 2*Q))/(4.*fpi**2*Pi**2)
				mpi=pi0mass
				totxcn= totxcn + Pi*1.E4 *(alphaQED*hbarc/smm)**2*((A8p**2*(2*w**2 - 2*mpi**2*(1 + 2*w) &
				- (mpi**4*(1 + 2*w)**2)/(4*w**2 + mpi**2*(1 + 2*w))))/(4.*(1 + 2*w)**2) &
				+ (A8p*(2*w**2 - mpi**2*(1 + 2*w) + 2*w*(1 + 2*w) + kap*(-mpi**2 + 2*w + 6*w**2 + 4*mpi**2*w**2 - 4*w**3) &
				+ kap**2*(-2*w**3 + mpi**2*w*(1 + 2*w))))/(2.*w*(1 + 2*w)**2) + (3*A8p**2*mpi**4*Log(1 + (4*w**2)/(mpi**2*(1 + 2*w))))/(16.*w**2) &
				+ (A8p*(8*(1 + kap)*w**3*Log(1/(1 + 2*w)) + (1 + kap)*mpi**6*Log(1 + (4*w**2)/(mpi**2*(1 + 2*w))) &
				- kap*(2 + kap)*mpi**4*w*(mpi**2 + 2*w)*Log(1 + (4*w**2)/(mpi**2*(1 + 2*w)))))/(8.*w**3*(mpi**2 + 2*w)))
				endif ! anomaly
			
			!	totxcn=2*Pi*1.E4 *(alphaQED*hbarc/smm)**2 *1./(2*w)*( (1.-2/w-2/w**2)*Log(1.+2*w)+1./2.+4./w-1/(2.*(1+2*w)**2) )
			!	write(*,*) kap, alphaE1,betaM1
			polxcn=-2*Pi*1.E4 *(alphaQED*hbarc/smm)**2*((-(kap**2*(2 + kap)**2*w**2*(1 + w)) + 8*(1 + 2*w))/(2.*w**2*(1 + 2*w)) &
			 + ((8 + 8*w - kap**2*(2 + kap)**2*w**2)*Log(1/(1 + 2*w)))/(4.*w**3) &
			 + 2*alphaE1*((2*(1 + w))/(1 + 2*w) + Log(1/(1 + 2*w))/w))
			write(340,'(1x,f10.5,3x,2g16.8)')  w*smm*1000., totxcn,polxcn
			ENDIF
			
			IF(.NOT.OUTPUT) RETURN
			if(.not.varyscalarPs) then
				write(352,103) Tglabmev, RTOT,0.5*(SI_PERP+SI_PAR+SI_PERP_LEX2+SI_PAR_LEX2)
				write(362,103) thetad, RTOT,0.5*(SI_PERP+SI_PAR+SI_PERP_LEX2+SI_PAR_LEX2)
				write(351,103) Tglabmev, Sigma3Born,Sigma3Born+Sigma3LEX2 !,Sigma3Born+Sigma3LEX2+Sigma3LEX4
				write(361,103) thetad, Sigma3Born,Sigma3Born+Sigma3LEX2 !,Sigma3Born+Sigma3LEX2+Sigma3LEX4
				write(651,103) Tglabmev, Sigma2zBorn,Sigma2zBorn+Sigma2zLEX3 !,Sigma3Born+Sigma3LEX2+Sigma3LEX4
				write(661,103) thetad, Sigma2zBorn,Sigma2zBorn+Sigma2zLEX3 !,Sigma3Born+Sigma3LEX2+Sigma3LEX4
				write(671,103) Tglabmev, Sigma2xBorn,Sigma2xBorn+Sigma2xLEX3 !,Sigma3Born+Sigma3LEX2+Sigma3LEX4
				write(681,103) thetad, Sigma2xBorn,Sigma2xBorn+Sigma2xLEX3 !,Sigma3Born+Sigma3LEX2+Sigma3LEX4
				write(353,103) Tglabmev,(SI_PERP-SI_PAR),SI_PERP-SI_PAR+SI_PERP_LEX2-SI_PAR_LEX2
				write(363,103) thetad,(SI_PERP-SI_PAR),SI_PERP-SI_PAR+SI_PERP_LEX2-SI_PAR_LEX2
				write(354,103) Tglabmev,(zLAB**2*SI_PERP-SI_PAR),zLAB**2*(SI_PERP+SI_PERP_LEX2)-SI_PAR-SI_PAR_LEX2
				write(364,103) thetad,(zLAB**2*SI_PERP-SI_PAR),zLAB**2*(SI_PERP+SI_PERP_LEX2)-SI_PAR-SI_PAR_LEX2
			else	
			SELECT CASE(iscalarPs)
			CASE(1)
				write(350,103) alphaE1/polza2LECs,Sigma3Born,Sigma3Born+Sigma3LEX2 !,Sigma3Born+Sigma3LEX2+Sigma3LEX4
				write(325,103) alphaE1/polza2LECs, RTOT,0.5*(SI_PERP+SI_PAR+SI_PERP_LEX2+SI_PAR_LEX2)
				write(375,103) alphaE1/polza2LECs,(SI_PERP-SI_PAR),SI_PERP-SI_PAR+SI_PERP_LEX2-SI_PAR_LEX2
				write(395,103) alphaE1/polza2LECs,(zLAB**2*SI_PERP-SI_PAR),zLAB**2*(SI_PERP+SI_PERP_LEX2)-SI_PAR-SI_PAR_LEX2
			CASE(2)
			    write(360,103) betaM1/polza2LECs,Sigma3Born,Sigma3Born+Sigma3LEX2 !,Sigma3Born+Sigma3LEX2+Sigma3LEX4
				write(325,103) betaM1/polza2LECs, RTOT,0.5*(SI_PERP+SI_PAR+SI_PERP_LEX2+SI_PAR_LEX2)
				write(385,103) betaM1/polza2LECs,(SI_PERP-SI_PAR),SI_PERP-SI_PAR+SI_PERP_LEX2-SI_PAR_LEX2
				write(405,103) betaM1/polza2LECs,(zLAB**2*SI_PERP-SI_PAR),zLAB**2*(SI_PERP+SI_PERP_LEX2)-SI_PAR-SI_PAR_LEX2
		  	END SELECT

			endif
		RETURN	
101       FORMAT(1X,4(2G16.8,3X))
102       FORMAT(1X,f6.1,3x,2G16.8)
103       FORMAT(1X,f7.2,3x,4G16.8)
		END SUBROUTINE LEXobservables
		  
        SUBROUTINE observables(MEP,MEP_ERR,RTOT,RTOT_ERR,diffRTOT,OUTPUT)
          IMPLICIT NONE
          LOGICAL OUTPUT
          COMPLEX*16, INTENT(in) :: MEP(2,2,2,2), MEP_ERR(2,2,2,2)
		  REAL*8, INTENT(out) :: RTOT,RTOT_ERR,diffRTOT
          COMPLEX*16 Mx, My,M11,M12,M21,M22
          INTEGER IH1,IH1P,I,J
          REAL*8 SI_PAR,SI_PERP,RPOL_ERR,RBETA_ERR,SI_Xp,SI_Xm,SI_Zp,SI_Zm
          REAL*8 SCALE,omegap,omega,zLAB
          REAL*8 basym,tasym,gasym,gasym_err,gtotal
		  INTEGER, PARAMETER :: nofPhipoints=5
	      INTEGER iphi
	      REAL*8 phi,cosphi,sinphi,xcsPhi(nofPhipoints),xcsPhihalfPi(nofPhipoints),Sigmaz(nofPhipoints)
          


          SI_PAR=0.d0
          SI_PERP=0.d0
          RTOT=0.d0
          RTOT_ERR=0.d0
          gasym_err=0.d0
          xcsPhi=0.d0
	      xcsPhihalfPi=0.d0
          DO 74 IH1=1,2
             DO 74 IH1P=1,2
                !WRITE(*,*) ih1p,ih1
                DO 74 I=1,2
                   !WRITE(*,*) (I,J,137.*MEP(I,J,IH1P,IH1)/smm,J=1,2)
                   Mx=-1./dsqrt(2.d0)*(MEP(I,1,IH1P,IH1)-MEP(I,2,IH1P,IH1))
                   My=(0.,-1.)/dsqrt(2.d0)*(MEP(I,1,IH1P,IH1)+MEP(I,2,IH1P,IH1))
                   SI_PAR=SI_PAR+Mx*CONJG(Mx)/2.d0
                   SI_PERP=SI_PERP+My*CONJG(My)/2.d0
                   DO iphi = 1,nofPhipoints
                         phi= (iphi-1)*180./nofPhipoints
	                 cosphi=dcosd(phi)
	                 sinphi=dsind(phi)
                         xcsPhi(iphi)=xcsPhi(iphi)+ (Mx*cosphi+My*sinphi)*CONJG(Mx*cosphi+My*sinphi)/2.d0
                         xcsPhihalfPi(iphi)=xcsPhihalfPi(iphi)+ (-Mx*sinphi+My*cosphi)*CONJG(-Mx*sinphi+My*cosphi)/2.d0
			Sigmaz(iphi)= ( xcsPhi(iphi) - xcsPhihalfPi(iphi) )/ ( xcsPhi(iphi) + xcsPhihalfPi(iphi) )
		   ENDDO
                   DO 74 J=1,2
                      RTOT=RTOT+MEP(I,J,IH1P,IH1)*CONJG(MEP(I,J,IH1P,IH1))/4.d0
                      !RTOT_ERR=RTOT_ERR+(2.d0 &
                       !    *cdabs(MEP_ERR(I,J,IH1P,IH1))*CDABS(MEP(I,J,IH1P,IH1)))**2/4.d0 
                      !gasym_err=gasym_err+(2.*cdabs(MEP_ERR(I,J,IH1P,IH1))*CDABS(MEP(I,J,IH1P,IH1)))**2/4.d0
74        CONTINUE

             SI_Xp=0.d0
			 SI_Xm=0.d0
			 SI_Zp=0.d0
			 SI_Zm=0.d0
             DO 751 IH1P=1,2
                !WRITE(*,*) ih1p,ih1
                DO 751 I=1,2
                DO 751 J=1,2
					SI_Xp=SI_Xp+(MEP(I,J,IH1P,1)+MEP(I,J,IH1P,2))*CONJG(MEP(I,J,IH1P,1)+MEP(I,J,IH1P,2))/2.d0
                    SI_Xm=SI_Xm+(3-2*J)*(MEP(I,J,IH1P,1)+MEP(I,J,IH1P,2))*CONJG(MEP(I,J,IH1P,1)+MEP(I,J,IH1P,2))/2.d0
					SI_Zp=SI_Zp + MEP(I,J,IH1P,1)*CONJG(MEP(I,J,IH1P,1))
                    SI_Zm=SI_Zm + (3-2*J)*MEP(I,J,IH1P,1)*CONJG(MEP(I,J,IH1P,1))
						!write(*,*) SI_Zp, SI_Zm
						!write(*,*) I,J, IH1P, MEP(I,J,IH1P,1)
751          CONTINUE

             gasym_err=dsqrt(gasym_err)/RTOT
!          Phase-space factor
            SCALE=(2.d0*Eqhat)**2/(4.d0*s) * 1.E7 *(hbarc/smm)**2
            omega=(s-1.d0)/2.d0
	        omegap=(1.d0-u)/2.d0
          if(.NOT.cmsangle) THEN
              SCALE=SCALE*(omegap/omega)**2*s
          ENDIF
          RTOT=RTOT*SCALE
          !RTOT_ERR=dsqrt(RTOT_ERR)*SCALE
          zLAB=dcos(thLAB)
          RTOT_ERR=4.d0*SCALE*(alphaQED**2) *omega*omegap *dsqrt( (1.d0+zLAB)**2+(2.d0*zLAB)**2 )
          RPOL_ERR=8.d0*SCALE*(alphaQED**2) *omega*omegap *dsin(thLAB)**2
          RBETA_ERR=8.d0*SCALE*(alphaQED**2) *omega*omegap *dsin(thLAB)**2*dabs(zLAB)

		  SI_PAR=SI_PAR*SCALE
          SI_PERP=SI_PERP*SCALE
		  diffRTOT=SI_PERP-SI_PAR
          IF(.NOT.OUTPUT) RETURN
          basym=100.d0*(SI_PAR-SI_PERP)/(SI_PAR+SI_PERP)
		  
		       WRITE(*,101) Tglabmev,thetad, RTOT,SI_Zm/SI_Zp,SI_Xm/SI_Xp,basym

		!write(*,*) iscalarPs,xscalarPs
         if(.not.varyscalarPs) then
          WRITE(32,101) Tglabmev,RTOT,RTOT_ERR
          WRITE(33,101) thetad,RTOT,RTOT_ERR

          WRITE(37,101) Tglabmev,-(SI_PAR-SI_PERP),RPOL_ERR
          WRITE(38,101) thetad,-(SI_PAR-SI_PERP),RPOL_ERR

          WRITE(39,101) Tglabmev,-(SI_PAR-zLAB**2*SI_PERP),RBETA_ERR
          WRITE(40,101) thetad,-(SI_PAR-zLAB**2*SI_PERP),RBETA_ERR

          write(42,101) Tglabmev,basym,100.d0*gasym_err*(0.1d0+dabs(basym))
          WRITE(43,101) thetad,basym,100.d0*gasym_err*(0.1d0+dabs(basym))
		  
		  write(72,102) Tglabmev,-SI_Xm/SI_Xp*100.d0
          WRITE(73,102) thetad,-SI_Xm/SI_Xp*100.d0
		  write(82,102) Tglabmev,-SI_Zm/SI_Zp*100.d0
          WRITE(83,102) thetad,-SI_Zm/SI_Zp*100.d0
		else
			SELECT CASE(iscalarPs)
			CASE(1)
				WRITE(320,101) alphaE1/polza2LECs,RTOT,RTOT_ERR
				WRITE(370,101) alphaE1/polza2LECs,-(SI_PAR-SI_PERP),RPOL_ERR
				WRITE(390,101) alphaE1/polza2LECs,-(SI_PAR-zLAB**2*SI_PERP),RBETA_ERR
				write(420,102) alphaE1/polza2LECs,basym,gasym_err*(0.1d0+dabs(basym))
			CASE(2)
				WRITE(330,101) betaM1/polza2LECs,RTOT,RTOT_ERR
				WRITE(380,101) betaM1/polza2LECs,-(SI_PAR-SI_PERP),RPOL_ERR
				WRITE(400,101) betaM1/polza2LECs,-(SI_PAR-zLAB**2*SI_PERP),RBETA_ERR
				WRITE(430,102) betaM1/polza2LECs,basym,gasym_err*(0.1d0+dabs(basym))
		  	END SELECT
		endif	


          WRITE(44,*)
          WRITE(44,101) Tglabmev,thetad
           DO iphi = 1,nofPhipoints
                  phi= (iphi-1)*180./nofPhipoints
	          WRITE(44,103) phi, SCALE*xcsPhi(iphi), SCALE*xcsPhihalfPi(iphi), &
         		Sigmaz(iphi), basym*dcosd(2.d0*phi)
                  !WRITE(44,102) phi, SCALE*(xcsPhi(iphi)+xcsPhihalfPi(iphi))/2.,RTOT
                  !WRITE(44,102) phi,dcosd(phi),dsind(phi)
                  !WRITE(44,*)
           ENDDO

          tasym=0.d0
          DO 75 IH1P=1,2
          DO 75 I=1,2                  
          DO 75 J=1,2
             SI_PAR=MEP(I,J,IH1P,1)*CONJG(MEP(I,J,IH1P,1))
             SI_PERP=MEP(I,J,IH1P,2)*CONJG(MEP(I,J,IH1P,2))
             tasym = tasym + (SI_PAR-SI_PERP)/(SI_PAR+SI_PERP)
 75     CONTINUE
          write(52,102) Tglabmev,tasym,gasym_err*(0.1d0+dabs(tasym))
          WRITE(53,102) thetad,tasym,gasym_err*(0.1d0+dabs(tasym))

          gasym=0.d0
	  gtotal=0.d0
          DO 76 IH1P=1,2
          DO 76 I=1,2                  
             M11=MEP(I,1,IH1P,1)
             M21=MEP(I,2,IH1P,1)
             M12=MEP(I,1,IH1P,2)
             M22=MEP(I,2,IH1P,2)
             gasym=gasym+M11*CONJG(M11)+M22*CONJG(M22)-(M12*CONJG(M12)+M21*CONJG(M21))
             gtotal=gtotal+M11*CONJG(M11)+M22*CONJG(M22)+(M12*CONJG(M12)+M21*CONJG(M21))
 76     CONTINUE
             gasym = -gasym/gtotal

          write(62,102) Tglabmev,gasym,gasym_err*(0.1d0+dabs(gasym))
          WRITE(63,102) thetad,gasym,gasym_err*(0.1d0+dabs(gasym))

101       FORMAT(1X,4(2G16.8,3X))
102       FORMAT(1X,f6.1,3x,2G16.8)
103       FORMAT(1X,f6.1,3x,4G16.8)
          RETURN
        END SUBROUTINE observables

        SUBROUTINE pwa_gaga(qq)
          REAL*8, INTENT(in) :: qq
          COMPLEX*16  mgg(2,2,2,2), mgg_err(2,2,2,2)
          REAL*8 th,dfun,tj,units,factor,xi,totalcs,total_err,totpolcs,xscn,polxscn,cserror,energy
          COMPLEX*16 tampl(6),phi(6),fmEE(3,2),fpEE(3,2),fpEM(3,2),aa,bb,cc
          COMPLEX, dimension(2) :: ab1, ab2, gaEE, gaEM
          INTEGER j,l,nparit,ix,ii,ipm(2)/1,-1/
          COMPLEX*16 A(8),fEE1m,fEE1m2
          character(LEN=2) EE(2),EM(2)
          common/temp/ A
          !INTEGER gaindex/2,4,3,1/  
          EE=(/"EE","MM"/); EM=(/"EM","ME"/) 
          ! hel=1/2 i=1   => -1/2, k=2
          ! hel=1/2 i=-1  => -3/2, k=4
          ! hel=-1/2 i=1  =>  3/2, k=3
          ! hel=-1/2 i=-1 =>  1/2, k=1
          units =1.d4*pionm
          write(*,*) ' LET fEE1- = ', -alphaQED/(3.)*units
          factor=(s+1.d0)/(2.*s)
          !factor=1.d0
          totalcs =0.d0
		  totpolcs =0.d0
          total_err = 0.d0
         do j=1,jmax,2
          tj = 0.5d0*j
          fEE1m=0.d0; fEE1m2=0.d0
          tampl=0.d0
          DO ix=1,ncosx
             th = dacos(cosx(ix))  
             CALL ini_polar(th)
             if(feynrule) then
                !CALL fr_gaga_ampl(qq,th,mgg)
             else
                CALL gaga_ampl(qq,th,mgg,mgg_err)
             endif
!             CALL gaga_ampl(qq,th,mgg,mgg_err)
             if(j.eq.1) then
                CALL observables(mgg,mgg_err,xscn,cserror,polxscn,.false.)
				!write(*,*) totalcs, totpolcs
                totalcs = totalcs + xscn*dcosx(ix)
				totpolcs = totpolcs + polxscn*dcosx(ix)
                total_err = total_err + cserror*dcosx(ix)
             endif

             !th=Pi+th
!!$             phi(1)= mgg(2,2,2,2)*factor
!!$             phi(2)= mgg(2,1,2,1)*factor
!!$             phi(3)= -mgg(2,2,1,2)*factor
!!$             phi(4)= mgg(1,2,2,2)*factor
!!$             phi(5)= mgg(1,1,2,2)*factor
!!$             phi(6)= mgg(1,2,2,1)*factor

             phi(1)= mgg(2,2,2,2)*factor
           phi(2)= mgg(2,1,2,1)*factor
           phi(4)= mgg(2,2,2,1)*factor
           phi(3)= mgg(2,1,2,2)*factor
           phi(5)= mgg(1,1,2,2)*factor
           phi(6)= mgg(1,2,2,1)*factor
           write(16,*) s, qq, -2.*qq*qq*(1.d0-cosx(ix)), phi(4)

           tampl(1)=tampl(1)+ phi(1)*dcosx(ix)*dfun(tj,-0.5d0,-0.5d0,th)
           tampl(2)=tampl(2)+ phi(2)*dcosx(ix)*dfun(tj,0.5d0,-0.5d0,th)

           xi=cosx(ix)
            if(j.gt.1) then
             tampl(3)=tampl(3)+ phi(3)*dcosx(ix)*dfun(tj,1.5d0,-0.5d0,th)
             tampl(4)=tampl(4)+ phi(4)*dcosx(ix)*dfun(tj,-1.5d0,-0.5d0,th)
             tampl(5)=tampl(5)+ phi(5)*dcosx(ix)*dfun(tj,-1.5d0,-1.5d0,th)
             tampl(6)=tampl(6)- phi(6)*dcosx(ix)*dfun(tj,1.5d0,-1.5d0,th)

            fEE1m2=fEE1m2+1./(16.*dsqrt(s))*( 2.*A(1)*(1.d0+xi**2)-A(3)*(3.-xi**2) &
                  -qhat**2 *(4.*A(6)+(2.*A(5)+A(4)+2.*A(2)+2./qhat)*xi)*(1.d0-xi**2) )*dcosx(ix)
             else
                ! fEE1- 
!!$                xi=cosx(ix)
!!$             fEE1m=fEE1m+1./(8.*dsqrt(s))*( (A(1)-A(3))*(1.d0+xi**2) &
!!$                  +qhat**2 *(2.*A(5)+A(4)-A(2)-1/qhat)*xi*(1.d0-xi**2) )*dcosx(ix)
!!$
             fEE1m2=fEE1m2+1./(8.*dsqrt(s))*( A(1)*(1.d0+xi**2)+A(3)*(3.-xi**2) &
                  +qhat**2 *(4.*A(6)+(2.*A(5)+A(4)-A(2)-1/qhat)*xi)*(1.d0-xi**2) )*dcosx(ix)

             endif
          enddo
          tampl=0.5d0*tampl ! factor 1/2 from the pw integral

          !write(*,*) tampl
          l=(j-1)/2
          write(*,*) 'l =', l
          !!write(*,*) 'fEE1-=', fEE1m*units
          !write(*,*) 'FEE1-=', fEE1m2*units

          !  fEE1- = fmEE(1,1),  fMM1- = fmEE(1,2)
          !  fEE1+ = fpEE(1,1),  fMM1+ = fmEE(1,2)
          do ii=1,2
             nparit=3-2*ii
             aa=(tampl(1)+nparit*tampl(2))
             bb=(tampl(3)+nparit*tampl(4))
             cc=(tampl(5)+nparit*tampl(6))
                 write(*,*)           
             !! write(*,'(6(1x,g12.4))') 0.5d0*aa*units/4.,nparit*dsqrt(l/(l+2.d0))*bb*units/4.,cc/8./3.*units 
             fmEE(l+1,ii)= 1./(l+1.d0)**2 *( 0.5d0*aa &
                  +nparit*dsqrt(l/(l+2.d0))*bb+ 0.5d0*l/(l+2.d0)*cc )
             write(*,*) 'f',EE(ii),'-', fmEE(l+1,ii)* units
            write(14,*) Tglabmev, 'f',EE(ii),l+1,'-', fmEE(l+1,ii)* units

             if(l>0) then
             !write(*,*)'d=', dsqrt(l*(l+2.d0))-dsqrt(3.d0)
              aa=tampl(1)-nparit*tampl(2)
                bb=(tampl(3)-nparit*tampl(4))
                cc=tampl(5)-nparit*tampl(6)
                write(*,*)
               !! write(*,'(6(1x,g12.4))') 0.5d0*aa*units/4.,3.*nparit/dsqrt(l*(l+2.d0))*bb*units/4.,3.*cc/8.*units 
 

                fpEE(l,ii) = 1./(l+1.d0)**2 *( 0.5d0*aa &
                  +nparit*(l+2.d0)/dsqrt(l*(l+2.d0))*bb +  (l+2.d0)/(2.d0*l)*cc ) 

               fpEM(l,ii) = 1./(l+1.d0)**2 *( -0.5d0*aa &
                  - nparit/dsqrt(l*(l+2.d0))*bb+ 0.5d0*cc )

                write(*,*) 'f',EE(ii),'+', fpEE(l,ii)* units
                 write(*,*)
               !!write(*,'(6(1x,g12.4))') -0.5d0*aa*units/4.,-nparit/dsqrt(l*(l+2.d0))*bb*units/4.,cc/8.*units  
               write(*,*) 'f',EM(ii),'+', fpEM(l,ii)* units
            write(15,*) Tglabmev, 'f',EE(ii),l,'+', fpEE(l,ii)* units
           write(16,*) Tglabmev, 'f',EM(ii),l,'+', fpEM(l,ii)* units
            endif
          enddo

          

       enddo ! matches do j=1,3,2

       do ii=1,2
          ab1(ii)= (2.*fpEE(1,ii)+ fmEE(1,ii))/qhat**2 *GeV2fm**3*1.d4
          ab2(ii)= 36.*(3.*fpEE(2,ii)+ 2.*fmEE(2,ii))/qhat**4 *GeV2fm**5*1.d4
          gaEE(ii)= (fpEE(1,ii) - fmEE(1,ii))/qhat**3 *GeV2fm**4*1.d4
          gaEM(ii)= 6.*fpEM(1,ii)/qhat**3 *GeV2fm**4*1.d4
       enddo
       !energy= qhat/smm*1.d3
       energy = Tglabmev
      write(501,111) energy, real(ab1(1)), aimag(ab1(1))
      write(502,111) energy, real(ab1(2)), aimag(ab1(2))
      write(503,111) energy, real(ab2(1)), aimag(ab2(1))
      write(504,111) energy, real(ab2(2)), aimag(ab2(2))
      write(505,111) energy, real(gaEE(1)), aimag(gaEE(1))
      write(506,111) energy, real(gaEE(2)), aimag(gaEE(2))
      write(507,111) energy, real(gaEM(1)), aimag(gaEM(1))
      write(508,111) energy, real(gaEM(2)), aimag(gaEM(2))
       totalcs=2.*Pi*totalcs*1.d-3
       total_err=2.*Pi*total_err*1.d-3
	   totpolcs=2.*Pi*totpolcs*1.d-3
          !WRITE(34,111) Tglabmev,totalcs,totpolcs !,total_err
		  WRITE(34,'("{",f10.1,a1,g18.10,"},")') Tglabmev,",",totalcs
 111   FORMAT(1x,f11.2,3x,2g18.10)    
    
       END SUBROUTINE pwa_gaga
!!$
!!$        SUBROUTINE fr_gaga_ampl(qq,th,mgg) 
!!$          USE frules
!!$          IMPLICIT NONE
!!$          REAL*8, INTENT(in) :: qq,th
!!$          COMPLEX*16, INTENT(out) :: mgg(2,2,2,2) 
!!$          COMPLEX*16 sum
!!$          LOGICAL proj 
!!$          INTEGER ih,irho,irhop
!!$          INTEGER jb,isob,ita,hindx,hpindx,mu,nu
!!$          INTEGER i,j
!!$          REAL*8, PARAMETER :: e2=1.d0/137.d0
!!$          REAL*8 rm,rm2,gg1,gg2
!!$          REAL*8 rt,ru,costh,sinth,const
!!$          REAL*8, DIMENSION(0:3) :: vp,vq,vpp,vqp,vqf,vpt,vpu,vcurnt
!!$          REAL*8 hel,help
!!$          TYPE(spinor) ubar,uu
!!$          TYPE(matrix4) prop,curnt(0:3),curntp(0:3)
!!$          TYPE(matrix4) rsprop(0:3,0:3)
!!$          TYPE(matrix4) fmgg(0:3,0:3),dum33(0:3,0:3)
!!$          REAL*8 pp2,pp2pi,formf
!!$          REAL*8 gg,gpin,pi0mass,epsil(0:3,0:3)
!!$
!!$
!!$          proj=.FALSE. ! ON/OFF spin-3/2 proj op.
!!$          !! c.m.s. kinematics
!!$          vq(0)=qhat
!!$          vq(1)=0.d0
!!$          vq(2)=0.d0
!!$          vq(3)=qhat
!!$          vp(0)=Eqhat
!!$          vp(1)=0.d0
!!$          vp(2)=0.d0
!!$          vp(3)=-qhat
!!$
!!$          costh=dcos(th)
!!$          sinth=dsin(th)
!!$          vqp(0)=qhat
!!$          vqp(1)=qhat*sinth
!!$          vqp(2)=0.d0
!!$          vqp(3)=qhat*costh
!!$          vpp(0)=Eqhat
!!$          vpp(1)=-qhat*sinth
!!$          vpp(2)=0.d0
!!$          vpp(3)=-qhat*costh
!!$          vqf=-vqp
!!$
!!$          vpu=vp-vqp
!!$          vpt=vqp-vq
!!$          ru=dot(vpu,vpu)
!!$          rt=dot(vpt,vpt)
!!$            !write(*,*) 'x =', real((s-1.)/2.)
!!$            !write(*,*) 't = ', real(rt)
!!$          DO mu=0,3
!!$             DO nu=0,3
!!$                fmgg(mu,nu)=gzero
!!$             ENDDO
!!$          ENDDO
!!$
!!$
!!$          ih=1 ! the nucleon pole diagrams
!!$          rm=bary(ih)%rm
!!$          rm2=bary(ih)%rm2
!!$          jb=bary(ih)%j
!!$          ita=bary(ih)%ip
!!$
!!$          curnt=v_nnvec(kapgp,vq)
!!$          curntp=v_nnvec(kapgp,vqf)
!!$          IF(bary(ih)%su(1)) THEN
!!$             prop=p_half(vps,rm)        
!!$             DO mu=0,3
!!$             DO nu=0,3
!!$                fmgg(mu,nu)=fmgg(mu,nu)-curntp(mu)*prop*curnt(nu)
!!$             ENDDO
!!$             ENDDO
!!$          ENDIF
!!$
!!$          IF(bary(ih)%su(2)) THEN
!!$                prop=p_half(vpu,rm)
!!$             DO mu=0,3
!!$             DO nu=0,3
!!$                fmgg(mu,nu)=fmgg(mu,nu)-curnt(nu)*prop*curntp(mu)
!!$             ENDDO
!!$             ENDDO
!!$          ENDIF
!!$
!!$          ih=3 ! the Delta pole diagrams
!!$          rm=bary(ih)%rm
!!$          rm2=bary(ih)%rm2
!!$          jb=bary(ih)%j
!!$          isob=bary(ih)%iso
!!$          ita=bary(ih)%ip
!!$            IF(gaugeinv) THEN
!!$                gg1=dng1*(3.d0/2.d0)/(1.d0+rm)*dsqrt(2.d0/3.d0)  ! isospin factor
!!$                gg2=dng2*(3.d0/2.d0)/(1.d0+rm)*dsqrt(2.d0/3.d0)
!!$                const=-1.d0
!!$            else
!!$                gg1=dg1
!!$                gg2=dg2
!!$               const=1.d0
!!$             endif
!!$          IF(bary(ih)%su(1)) THEN
!!$
!!$             IF(gaugeinv) THEN
!!$               dum33=contract(v_gnD(gg1,gg2,vqp,vpp),p_thalfsd(vps,rm),v_gnD(gg1,gg2,vq,vp))
!!$             ELSE
!!$                dum33=contract(v_gnD_conv(1,gg1,gg2,dz1,dz2,vqf,vpp),&
!!$                     p_thalfsd(vps,rm),&
!!$                     v_gnD_conv(0,gg1,gg2,dz1,dz2,vq,vp))
!!$             ENDIF
!!$             DO mu=0,3
!!$             DO nu=0,3
!!$                fmgg(mu,nu)=fmgg(mu,nu)+const*dum33(mu,nu)
!!$             ENDDO
!!$             ENDDO
!!$          ENDIF
!!$
!!$          IF(bary(ih)%su(2)) THEN
!!$             IF(gaugeinv) THEN
!!$                 gg1=dng1*dsqrt(2.d0/3.d0)  ! isospin factor
!!$                 gg2=dng2*dsqrt(2.d0/3.d0)
!!$               dum33=contract(v_gnD(gg1,gg2,vq,vpp),p_thalfsd(vpu,rm),v_gnD(gg1,gg2,vqp,vp))
!!$             ELSE 
!!$                gg1=dg1
!!$                gg2=dg2
!!$               dum33=contract(v_gnD_conv(1,gg1,gg2,dz1,dz2,vq,vpp),&
!!$                     p_thalfsd(vpu,rm),v_gnD_conv(0,gg1,gg2,dz1,dz2,vqf,vp))  
!!$             ENDIF
!!$             DO mu=0,3
!!$                DO nu=0,3
!!$                   fmgg(mu,nu)=fmgg(mu,nu)+const*dum33(nu,mu)
!!$                ENDDO
!!$             ENDDO
!!$          ENDIF
!!$
!!$          !! anomaly contribution
!!$          IF(anomaly) THEN
!!$             pi0mass=0.136d0/smm
!!$             gpin=dsqrt(13.8d0*4.d0*PI)
!!$             gg=gA/(2.d0*PI*fpi)**2/(pi0mass**2-rt)
!!$             !gg=gpin*(0.0375/pi0mass)/(pi0mass**2-rt)
!!$             write(*,*) 'anomaly gg = ', gg*(pi0mass**2-rt)
!!$             CALL anti2(vq,vqp,epsil)             
!!$             DO mu=0,3
!!$                DO nu=0,3       
!!$                   fmgg(mu,nu)=fmgg(mu,nu)+(0.d0,1.d0)*gg*epsil(mu,nu)*gamm(5)
!!$                ENDDO
!!$             ENDDO
!!$          ENDIF
!!$
!!$          DO mu=0,3
!!$             DO nu=0,3
!!$                fmgg(mu,nu)=fmgg(mu,nu)*e2
!!$             ENDDO
!!$          ENDDO
!!$
!!$
!!$          !! helicity amplitudes
!!$          DO i=1,2
!!$             DO j=1,2
!!$                DO hindx=1,2
!!$                   DO hpindx=1,2
!!$                      mgg(i,j,hpindx,hindx)=0.d0
!!$                  ENDDO
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$       
!!$          DO hindx=1,2
!!$             hel=-(hindx-1.5d0)
!!$             uu=u_(qhat,0.d0+PI,0.d0,hel,1)
!!$             DO hpindx=1,2
!!$                help=-(hpindx-1.5d0)
!!$                ubar=ubar_(qhat,th+PI,0.d0,help,1)
!!$                !!current conservation test
!!$                IF(hindx.EQ.hpindx) THEN
!!$                   DO mu=0,3
!!$                     ! WRITE(*,*) mu, '(qj)=',&
!!$                      !     braket(ubar,fmgg(mu,0)-fmgg(mu,3),uu)
!!$                      WRITE(7,*) mu, '(qj)=',&
!!$                           braket(ubar,fmgg(0,mu)-fmgg(3,mu),uu)
!!$                   ENDDO
!!$                ENDIF
!!$                DO i=1,2
!!$                   DO j=1,2
!!$                      sum=0.d0
!!$                      DO mu=0,3
!!$                      DO nu=0,3
!!$                         sum=sum &
!!$                              + gmunu(mu)*vepsp(mu,i)*gmunu(nu)*veps(nu,j) &
!!$                              *braket(ubar,fmgg(mu,nu),uu)
!!$                      ENDDO
!!$                      ENDDO
!!$                      if(dabs(dreal(sum)).lt.1.e-10) sum=(0.d0,1.d0)*dimag(sum)
!!$                      mgg(i,j,hpindx,hindx)=sum
!!$                   ENDDO
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$
!!$          RETURN
!!$        END SUBROUTINE fr_gaga_ampl

       END PROGRAM
        
!!________________________________________ ! VP, 17/09/95
      REAL*8 FUNCTION dfun(j,mp,m,beta)    ! reduced d-functions
        REAL*8 j,mp,m,beta                 ! notation Edmons=Varshalovich
        REAL*8 c2,s2,fact
        INTEGER l,max,min,jmp,jm,jmp1,jm1,mpm 

        IF(j.LT.0.OR.ABS(m).GT.j.OR.ABS(mp).GT.j) &
             PAUSE 'in Dfun: incorrect range of arguments' 
        c2=dcos(beta*0.5d0) 
        s2=dsin(beta*0.5d0) 
        jmp = j+mp
        jm = j-m
        jmp1 = j-mp
        jm1 = j+m
        mpm = mp-m
        min = 0
        IF(mpm.GT.0) min = mpm
        max = jmp
        IF(jm.LT.jmp) max = jm    
        dfun=0.d0 
        DO l=min,max
           dfun = dfun + (-1)**l*c2**(jm+jmp-2*l)*s2**(2*l-mpm) &
                /( fact(jmp-l)*fact(jm-l)*fact(l)*fact(l-mpm) ) 
        ENDDO
        dfun = dfun*dsqrt(fact(jmp)*fact(jmp1)*fact(jm1)*fact(jm))  
        RETURN
      END FUNCTION dfun

      REAL*8 FUNCTION fact(n)  ! factorial n
        INTEGER n,i
        IF(n.LT.0) PAUSE 'in Fact: factorial of a negative number?'
        fact = 1.d0 
        IF(n.LE.1) RETURN
        DO i=2,n
           fact=fact*i 
        ENDDO
        RETURN       
      END FUNCTION fact
