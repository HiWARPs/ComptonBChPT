!!****************************************************************************
!!* Module: frules.f90               Version: 2.2     Update: 09/13/98
!!*
!!*___________________________________________________________________________

   MODULE frules
     USE dirac 
     REAL*8, PARAMETER :: pimass=0.139d0/0.9383d0 ! pion mass     
   CONTAINS
            
     FUNCTION P_half(p,rm) RESULT(prop)
       !!    spin-1/2 propagator 
       IMPLICIT NONE
       REAL*8, INTENT(in) :: p(0:3),rm
       REAL*8 s
       TYPE(matrix4) prop
       s=dot(p,p)
       prop=(slash(p)+rm*gone)*(1.d0/(s-rm*rm))
       RETURN
     END FUNCTION P_half

     FUNCTION P_Chalf(p,rm,c2) RESULT(prop)
       !!   Composite  spin-1/2 propagator 
       IMPLICIT NONE
       REAL*8, INTENT(in) :: p(0:3),rm,c2
       REAL*8 s
       TYPE(matrix4) prop,psl
       s=dot(p,p)
       psl=slash(p)
       prop=(psl+rm*gone)*(1.d0/(s-rm*rm))
       prop=prop*(c2*psl+rm*(1.-c2)*gone)*(1/(rm*(1.d0-c2)**2-c2*c2*s/rm))
       RETURN
     END FUNCTION P_Chalf

     FUNCTION P_thalf(p,rm,sproj) RESULT(rsprop)
       !!    spin-3/2 propagator 
       IMPLICIT NONE
       REAL*8, INTENT(in) :: p(0:3),rm
       LOGICAL, INTENT(in) :: sproj
       INTEGER mu,nu
       REAL*8 s,p2
       TYPE(matrix4) slp,prop,rsprop(0:3,0:3)

       s=dot(p,p)
       p2=rm**2
       IF(sproj) p2=s
       slp=slash(p)
       prop=(slp+rm*gone)*(1.d0/(s-rm*rm))
       DO mu=0,3
          DO nu=0,3
             rsprop(mu,nu)=(-1.d0/3.d0)*prop*( gamm(mu)*gamm(nu) &
                  +(slp*gamm(mu)*p(nu)+p(mu)*gamm(nu)*slp)*(1./p2) )
          ENDDO
          rsprop(mu,mu)=rsprop(mu,mu)+prop*gmunu(mu)
       ENDDO
       RETURN
     END FUNCTION P_thalf
     

     FUNCTION P_thalfWidth(p,rm,wid,sproj) RESULT(rsprop)
       !!    spin-3/2 propagator 
       IMPLICIT NONE
       REAL*8, INTENT(in) :: p(0:3),rm,wid
       LOGICAL, INTENT(in) :: sproj
       INTEGER mu,nu
       REAL*8 s,p2
       TYPE(matrix4) slp,prop,rsprop(0:3,0:3)

       s=dot(p,p)
       p2=rm**2
       IF(sproj) p2=s
       slp=slash(p)
       prop=(slp+rm*gone)*(1.d0/(s-rm*rm+CI*rm*wid))
       DO mu=0,3
          DO nu=0,3
             rsprop(mu,nu)=(-1.d0/3.d0)*prop*( gamm(mu)*gamm(nu) &
                  +(slp*gamm(mu)*p(nu)+p(mu)*gamm(nu)*slp)*(1./p2) )
          ENDDO
          rsprop(mu,mu)=rsprop(mu,mu)+prop*gmunu(mu)
       ENDDO
       RETURN
     END FUNCTION P_thalfWidth

     FUNCTION P_thalfsd(p,rm) RESULT(rsprop)
       !!    standard spin-3/2 propagator 
       IMPLICIT NONE
       REAL*8, INTENT(in) :: p(0:3),rm
       INTEGER mu,nu
       REAL*8 s,p2, edem
       TYPE(matrix4) slp,prop,rsprop(0:3,0:3)

       s=dot(p,p)
       p2=rm*rm
       slp=slash(p)
       !slp=rm*gamm(0) ! positive energy proj
       prop=(slp+rm*gone)*(1.d0/(s-rm*rm))
       DO mu=0,3
          DO nu=0,3
             rsprop(mu,nu)=(-1.d0/3.d0)*prop*( gamm(mu)*gamm(nu) &
                  +(gamm(mu)*p(nu)-p(mu)*gamm(nu))*(1.d0/rm) &
                  + (2.d0/p2)*p(mu)*p(nu)*gone )
          ENDDO
          rsprop(mu,mu)=rsprop(mu,mu)+prop*gmunu(mu)
       ENDDO
       RETURN
     END FUNCTION P_thalfsd

     FUNCTION P_thalf_proj(p,rm) RESULT(rsprop)
       !!    standard spin-3/2 pos.energy proj. operator
       IMPLICIT NONE
       REAL*8, INTENT(in) :: p(0:3),rm
       INTEGER mu,nu
       REAL*8 s,p2, edem
       TYPE(matrix4) slp,prop,rsprop(0:3,0:3)

       s=dot(p,p)
       p2=rm*rm
       slp=slash(p)
       !slp=rm*gamm(0) ! positive energy proj
       prop=(slp+rm*gone)*(1.d0/(2.*rm))
       DO mu=0,3
          DO nu=0,3
             rsprop(mu,nu)=(-1.d0/3.d0)*prop*( gamm(mu)*gamm(nu) &
                  +(gamm(mu)*p(nu)-p(mu)*gamm(nu))*(1.d0/rm) &
                  + (2.d0/p2)*p(mu)*p(nu)*gone )
          ENDDO
          rsprop(mu,mu)=rsprop(mu,mu)+prop*gmunu(mu)
       ENDDO
       RETURN
     END FUNCTION P_thalf_proj

     FUNCTION v_nnpi(ita,vk) RESULT(nnpi)
       !!    (pseudo)vector  pi-N-N Vertex
       INTEGER, INTENT(in) :: ita          
       REAL*8, INTENT(in) :: vk(0:3)
       TYPE(matrix4) nnpi
       nnpi=0.5d0*slash(vk)
       IF(ita.EQ.-1) nnpi=gamm(5)*nnpi
       RETURN
     END FUNCTION v_nnpi
     
     FUNCTION v_nnvec(kappa,vk) RESULT(nnv)
       !! N-N-vector meson vertex
       !! mementum vk directed into vertex (otherwise change sign KAPPA)
       !! NNV(mu)=[\gamma^{mu}- KAPPA/2M_nucl \sigma^{mu,nu} q_{nu}]
       IMPLICIT NONE
       REAL*8, INTENT(in) :: kappa,vk(0:3)
       TYPE(matrix4) nnv(0:3)
       INTEGER mu
       TYPE(matrix4) slk
       slk= slash(vk)
       DO mu=0,3
          nnv(mu)=gamm(mu)-0.25d0*kappa*(gamm(mu)*slk-slk*gamm(mu))
       ENDDO
       RETURN
     END FUNCTION v_nnvec


     
     FUNCTION v_nRvec(kappa,vk) RESULT(nnv)
       !! N-Roper-vector meson vertex
       !! mementum vk directed into vertex (otherwise change sign KAPPA)
       !! NNV(mu)=[- KAPPA/2M_nucl \sigma^{mu,nu} q_{nu}]
       IMPLICIT NONE
       REAL*8, INTENT(in) :: kappa,vk(0:3)
       TYPE(matrix4) nnv(0:3)
       INTEGER mu
       TYPE(matrix4) slk
       slk= slash(vk)
       DO mu=0,3
          nnv(mu)=(-0.25d0*kappa)*(gamm(mu)*slk-slk*gamm(mu)) 
       ENDDO
       RETURN
     END FUNCTION v_nRvec

     FUNCTION v_dnpi(api,vk) RESULT(ndpi)
       !! conventional P-N->Delta Vertex
       IMPLICIT NONE
       REAL*8, INTENT(in) :: api,vk(0:3)
       TYPE(matrix4) ndpi(0:3)
       INTEGER mu
       DO mu=0,3
          ndpi(mu)=(-1.d0/pimass)*(vk(mu)*gone+api*gamm(mu)*slash(vk))
       ENDDO
       RETURN
     END FUNCTION v_dnpi
     
     FUNCTION v_pind(api,vk) RESULT(ndpi)
       !! conventional Delta->N-Pion Vertex
       IMPLICIT NONE
       REAL*8, INTENT(in) :: api,vk(0:3)
       TYPE(matrix4) ndpi(0:3)
       INTEGER mu
       DO mu=0,3
          ndpi(mu)=(-1.d0/pimass)*(vk(mu)*gone+api*slash(vk)*gamm(mu))
       ENDDO
       RETURN
     END FUNCTION v_pind
          
     FUNCTION v_pindelta(vp,vk,rm) RESULT(ndpi)
       !! gauge-invariant piNDelta Vertex
       IMPLICIT NONE
       REAL*8, INTENT(in) :: rm,vp(0:3),vk(0:3)
       TYPE(matrix4) ndpi(0:3)
       REAL*8 epsil(0:3,0:3)
       INTEGER mu,nu
       COMPLEX*16 cc

       cc=1.d0/(rm)
       CALL anti2(vp,vk,epsil)
       DO mu=0,3
          ndpi(mu)=gzero
          DO nu=0,3
             ndpi(mu)=ndpi(mu)+cc*gmunu(nu)*epsil(mu,nu)*gamm(nu)*gamm(5)
          ENDDO
       ENDDO
       RETURN
     END FUNCTION v_pindelta

     SUBROUTINE anti2(v1,v2,epsil)
       IMPLICIT NONE
       REAL*8, INTENT(in) :: v1(0:3),v2(0:3)
       REAL*8 epsil(0:3,0:3)
       INTEGER mu,nu,al,be
       DO mu=0,3
          DO nu=0,3
             epsil(mu,nu)=0.d0
             DO al=0,3
                DO be=0,3
                   epsil(mu,nu)=epsil(mu,nu)+gmunu(al)*gmunu(be) &
                        *epsilon4(mu,nu,al,be)*v1(al)*v2(be)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       RETURN
     END SUBROUTINE anti2

     FUNCTION v_gnD_conv(iflag,G1,G2,z1,z2,q,p) RESULT(dver)
       !!       see 'nath et al., phys. rev D for info on a
       IMPLICIT NONE
       INTEGER, INTENT(in) :: iflag     ! 0 for n-delta, otherwise d-n
       COMPLEX*16, INTENT(in) :: G1,G2
       REAL*8, INTENT(in) :: z1,z2,q(0:3),p(0:3)  ! q-photon, p-nucleon
       TYPE(matrix4) DVER(0:3,0:3)
       INTEGER al,be,mu
       REAL*8 dotpq
       COMPLEX*16 f
       TYPE(matrix4) O_1(0:3,0:3),O_2(0:3,0:3)
       TYPE(matrix4) gf,term1,term2,slashq
       
       slashq= slash(q)
       dotpq = dot(p,q)
       
       DO al=0,3
          DO be=0,3
             gf=gamm(al)*gamm(be)
             O_1(al,be)=(-0.5d0-z1)*gf
             O_2(al,be)=(-0.5d0-z2)*gf
          ENDDO
          gf=gmunu(al)*gone
          O_1(al,al)=gf+O_1(al,al)
          O_2(al,al)=gf+O_2(al,al)
       ENDDO
       
       IF(iflag.EQ.0) THEN
          DO al=0,3
             DO mu=0,3
                term1=gzero
                term2=gzero
                DO be=0,3
                   F=kronikers(1+mu,1+be)*gmunu(mu)
                   term1=term1 - O_1(al,be)*gmunu(be) &
                        *(F*slashq-gamm(mu)*q(be))*gamm(5)
                   term2=term2 + O_2(al,be)*gmunu(be) &
                        *(dotpq*F-p(mu)*q(be) )*gamm(5)
                ENDDO
                dver(al,mu)=0.5d0*(G1*term1+0.5d0*G2*term2)
             ENDDO
          ENDDO
       ELSE
          DO al=0,3
             DO mu=0,3
                term1=gzero
                term2=gzero
                DO be=0,3
                   F=kronikers(1+mu,1+be)*gmunu(mu)
                   term1=term1 - gamm(5) &
                        *(F*slashq-gamm(mu)*q(be))*O_1(be,al)*gmunu(be)
                   term2=term2 + gamm(5) &
                        *(dotpq*F-p(mu)*q(be) )*O_2(be,al)*gmunu(be)
                ENDDO
                dver(al,mu)=0.5d0*(G1*term1+0.5d0*G2*term2)
             ENDDO
          ENDDO
       ENDIF
      END FUNCTION v_gnD_conv

     FUNCTION v_gnD(G1,G2,q,p) RESULT(dver)
       !! Double gauge-inv. gamma-N-Delta vertex
       IMPLICIT NONE
       COMPLEX*16, INTENT(in) :: G1,G2
       REAL*8, INTENT(in) :: q(0:3),p(0:3)  ! q,mu-photon, p,al-Delta
       TYPE(matrix4) DVER(0:3,0:3)
       INTEGER al,be,mu
       REAL*8 epsil(0:3,0:3),dotpq
       TYPE(matrix4) term1,term2
       
       dotpq = dot(p,q)
       CALL anti2(p,q,epsil)
       DO al=0,3
          DO mu=0,3
             term1=epsil(al,mu)*gone
             term2=CI*(dotpq*kronikers(1+mu,1+al)*gmunu(mu)-q(al)*p(mu))*gamm(5)
             dver(al,mu)=(G1*term1 - G2*term2)
             !dver(al,mu)=0.5d0*(G1*(term1-0.5d0*term2)+0.5d0*G2*term2)
          ENDDO
       ENDDO
     END FUNCTION v_gnD

     FUNCTION v_gpivec(q,p) RESULT(dver)
       !! gauge-inv. gamma-pi-rho vertex
       IMPLICIT NONE
       REAL*8, INTENT(in) :: q(0:3),p(0:3)  ! q-photon, p-nucleon
       TYPE(matrix4) DVER(0:3,0:3)
       INTEGER al,mu
       REAL*8 epsil(0:3,0:3)
       
       CALL anti2(p,q,epsil)
       DO al=0,3
          DO mu=0,3
             dver(al,mu)=epsil(al,mu)*gone
          ENDDO
       ENDDO

     END FUNCTION v_gpivec


     
!!$    FUNCTION v_mod_nnvec(kappa,vpp,vp) RESULT(nnv)
!!$       !! N-N-gamma vertex with formfactors
!!$       !! mementum vk directed into vertex (otherwise change sign KAPPA)
!!$       !! NNV(mu)=[\gamma^{mu}- KAPPA/2M_nucl \sigma^{mu,nu} q_{nu}]
!!$       IMPLICIT NONE
!!$       REAL*8, INTENT(in) :: kappa,vpp(0:3),vp(0:3)
!!$       TYPE(matrix4) nnv(0:3)
!!$       INTEGER mu
!!$       Real*8 vk(0:3), vpl(0:3)
!!$       real*8 fp,fpp,p2,pp2,theta,fnucl
!!$       TYPE(matrix4) wti, slk
!!$       vk=vpp-vp
!!$       vpl=vpp+vp
!!$       p2=dot(vp,vp)
!!$       pp2=dot(vpp,vpp)
!!$       fp=1.d0/(fnucl(p2))
!!$       fpp=1.d0/(fnucl(pp2))
!!$       theta=(fpp-fp)/(pp2-p2)
!!$       wti=(fpp*(slash(vpp)-gone) + fp*(slash(vp)-gone))*theta
!!$       slk= slash(vk)
!!$       DO mu=0,3
!!$          nnv(mu)=fp*fpp*gamm(mu) + wti*vpl(mu) &
!!$               + 0.25d0*kappa*(slk*gamm(mu)-gamm(mu)*slk)
!!$       ENDDO
!!$       RETURN
!!$     END FUNCTION v_mod_nnvec
!!$
!!$     FUNCTION v_mod_pipig(vpp,vp) RESULT(pipiv)
!!$       !! pi-pi-gamma vertex with formfactors
!!$       IMPLICIT NONE
!!$       REAL*8, INTENT(in) :: vpp(0:3),vp(0:3)
!!$       real*8 pipiv(0:3)
!!$       INTEGER mu
!!$       Real*8 vpl(0:3)
!!$       real*8 fp,fpp,p2,pp2,theta,fpion,wti
!!$       real*8 m(4),msq(4)
!!$       common/masses/ m,msq
!!$     
!!$       vpl=vpp+vp
!!$       p2=dot(vp,vp)
!!$       pp2=dot(vpp,vpp)
!!$       fp=1.d0/(fpion(p2))
!!$       fpp=1.d0/(fpion(pp2))
!!$       theta=(fpp-fp)/(pp2-p2)
!!$       wti=fpp*fp + ( fpp*(pp2-msq(3)) + fp*(p2-msq(3)) )*theta
!!$       DO mu=0,3
!!$          pipiv(mu)=vpl(mu)*wti
!!$       ENDDO
!!$       RETURN
!!$     END FUNCTION v_mod_pipig
!!$
!!$        REAL*8  FUNCTION fpion(p2)   ! cutoff for pion
!!$          IMPLICIT NONE
!!$          REAL*8, INTENT(in):: p2  
!!$          REAL*8 m(4),msq(4)
!!$          COMMON /masses/ m,msq
!!$          REAL*8 cutpion
!!$          COMMON/cutp/cutpion
!!$          fpion=(cutpion - msq(3))/(cutpion - p2)
!!$          RETURN
!!$        END FUNCTION fpion
!!$
!!$        REAL*8  FUNCTION fnucl(p2)   ! cutoff for nucleon
!!$          IMPLICIT NONE
!!$          REAL*8, INTENT(in):: p2 
!!$          REAL*8 m(4),msq(4), f_PJ
!!$          COMMON /masses/ m,msq
!!$          fnucl=1.d0
!!$          fnucl=f_PJ(p2,1.d0,msq(2))
!!$          RETURN
!!$        END FUNCTION fnucl
!!$
!!$        REAL*8  FUNCTION f_PJ(p2,rm2,cutm)   ! Pearce&Jennings cutoff
!!$          REAL*8, INTENT(in) ::  p2,rm2,cutm 
!!$          INTEGER nal
!!$          F_PJ=1.d0
!!$          nal=2
!!$          !cutf=nal*(cutm)**2
!!$          cutf=nal*(cutm/2.d0)**2
!!$          f_PJ= cutf/(cutf+(rm2-p2)**2)
!!$          f_PJ=f_PJ**nal
!!$          RETURN
!!$        END FUNCTION f_PJ
!!$
!!$        REAL*8  FUNCTION f_OBE(p2,rm2,cutm)   ! OBE cutoff 
!!$          IMPLICIT NONE
!!$          REAL*8, INTENT(in) :: p2,rm2,cutm
!!$          f_OBE= (cutm)/(cutm - p2) 
!!$          RETURN
!!$        END FUNCTION f_OBE
    END MODULE frules

