#include common.h
Symbols Q,kap;
#procedure KSUB
#endprocedure
* s-channel Born: Gam = - V(ep,out) (ps/+M) V(ee,in) * ixs,  ps = p+q
*  V(ee,k_in) = Q ee/ - kap/2/M ag(ee,k_in);  outgoing photon: k_in = -qp
L GAM = (-ixs)*(Q*g_(1,ep)+kap/2/M*(g_(1,ep,qp)-ep.qp))*(g_(1,p)+g_(1,q)+M*gi_(1))*(Q*g_(1,ee)-kap/2/M*(g_(1,ee,q)-ee.q))
       +(-iys)*(Q*g_(1,ee)-kap/2/M*(g_(1,ee,q)-ee.q))*(g_(1,p)-g_(1,qp)+M*gi_(1))*(Q*g_(1,ep)+kap/2/M*(g_(1,ep,qp)-ep.qp));
.sort
#call Proj
#call Finish
B ixs,iys,iPq;
Format doublefortran;
Print +s;
.end
