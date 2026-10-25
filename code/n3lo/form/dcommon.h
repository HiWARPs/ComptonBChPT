* dcommon.h -- additions for graphs with Delta(1232) lines (work/n3lo/py/dgraphs.py, delta.py)
*   Scheme as common.h: external photon polarisations 4-dimensional (sums after the tensor reduction); the rank
*   reduction acts before, with l.ee, l.ep irreducible.  (PolSumD: D-dim sums, kept for tests; not consistent
*   with the 4-dim Gram matrices used for the projections.)
*   MD = M_Delta, iD1 = 1/(D-1) (D-dim RS propagator); loop-momentum reduction by topology (TOPO bitmask,
*   DD(i) = inverse propagator i, SHIFT = re-routing), tensor coefficients TCFT(mask, n00, a, b, c).
#include common.h
Symbols MD,iD1,STUCK,ixD,iyD,iMD;
Indices a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16;
CFunctions TOPO,DD,SHIFT,TCFT,UCFT;
Set DI: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16;

#procedure ContractD
* D-dim contraction of the Delta-index pairs inside the gamma string before the trace
* (g^i A g^j g_i = 2 g^j A - g^i A g_i g^j); 4x faster traces for the graphs with Delta lines, identical results
  repeat;
    id g_(1,?a,i1?DI,i1?DI,?c) = D*g_(1,?a,?c);
    id g_(1,?a,i1?DI,?b,i2?,i1?DI,?c) = 2*g_(1,?a,i2,?b,?c) - g_(1,?a,i1,?b,i1,i2,?c);
    id g_(1,?a,i1?DI,?b,w1?,i1?DI,?c) = 2*g_(1,?a,w1,?b,?c) - g_(1,?a,i1,?b,i1,w1,?c);
  endrepeat;
#endprocedure

#procedure KinV
  id p2 = p + q;
  id P = 2*p + q;
  id p.p = M^2;
  id q.q = -dl;
  id p.q = dl/2;
#endprocedure

#procedure KinS
  id p.p = ss;
#endprocedure

#procedure Dexp
  id TOPO(n1?)*TCF(?x) = TCFT(n1,?x);
  .sort
  id D = 4 + dd;
#ifdef `RS4D'
  id iD1 = 1/3;
#else
  id iD1 = 1/3 - dd/9;
#endif
  .sort
  id dd^n1?{>1} = 0;
  id dd*TCFT(?x) = -2*UCFT(?x);
  id dd = 0;
  .sort
#endprocedure

#procedure PolSumD
* P-gauge polarisation sums in D dimensions (CDR; the sums are taken before the tensor reduction):
*   sum ee(a)ee(b) = -g(a,b) + (P(a)q(b)+q(a)P(b))/P.q - P^2 q(a)q(b)/(P.q)^2,  trace 2-D
  id ee.ee = 2-D;
  id ee.w1?^2 = -w1.w1 + 2*w1.P*w1.q*iPq - PP*w1.q^2*iPq^2;
  id ee.w1?*ee.w2? = -w1.w2 + (w1.P*w2.q+w1.q*w2.P)*iPq - PP*w1.q*w2.q*iPq^2;
  #call Kin
  id ep.P = 0;
  id ep.qp = 0;
  id ep.p = -ep.q/2;
  .sort
  id ep.ep = 2-D;
  id ep.w1?^2 = -w1.w1 + 2*w1.P*w1.qp*iPq - PP*w1.qp^2*iPq^2;
  id ep.w1?*ep.w2? = -w1.w2 + (w1.P*w2.qp+w1.qp*w2.P)*iPq - PP*w1.qp*w2.qp*iPq^2;
  #call Kin
#endprocedure

#procedure FinishR
  Drop GAM;
  .sort
  #call ContractD
  .sort
  tracen,1;
  .sort
  #call Kin
  #call Gauge
  .sort
  #call REDUCE
  #call TensRed
  .sort
  #call KSUBT
  #call Kin
  #call Gauge
  .sort
  #call PolSum
  .sort
  #call Dexp
  id M = 1;
  .sort
#endprocedure

#procedure GFinishR
  Drop GAM;
  .sort
  #call ContractD
  .sort
  tracen,1;
  .sort
  #call Kin
  id ep.qp = 0;
  id ep.p = -ep.q/2;
  .sort
  #call REDUCE
  #call TensRed
  .sort
  #call KSUBT
  #call Kin
  id ep.qp = 0;
  id ep.p = -ep.q/2;
  .sort
  id ep.ep = -2;
  id ep.w1?^2 = -w1.w1 + 2*w1.P*w1.qp*iPq - PP*w1.qp^2*iPq^2;
  id ep.w1?*ep.w2? = -w1.w2 + (w1.P*w2.qp+w1.qp*w2.P)*iPq - PP*w1.qp*w2.qp*iPq^2;
  #call Kin
  .sort
  #call Dexp
  id M = 1;
  .sort
#endprocedure

#procedure VFinishR
  #call ContractD
  .sort
  tracen,1;
  .sort
  #call KinV
  .sort
  #call REDUCE
  #call TensRed
  .sort
  #call KSUBT
  #call KinV
  .sort
  id ee.ee = -4;
  id ee.w1?^2 = -w1.w1;
  id ee.w1?*ee.w2? = -w1.w2;
  #call KinV
  .sort
  #call Dexp
#endprocedure

#procedure SFinishR
  L V1 = GAM/4;
  L V2 = g_(1,p)*GAM/4;
  .sort
  Drop GAM;
  .sort
  #call ContractD
  .sort
  tracen,1;
  .sort
  #call KinS
  .sort
  #call REDUCE
  #call TensRed
  .sort
  #call KSUBT
  #call KinS
  .sort
  #call Dexp
#endprocedure
