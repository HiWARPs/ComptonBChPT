* common.h -- declarations and procedures for the one-loop Compton graphs (work/n3lo)
*
* Conventions (M kept symbolic; the generated Fortran sets M=1):
*   p (in nucleon) + q (in photon, pol ee) = pp (out nucleon) + qp (out photon, pol ep);  P = p + pp
*   x = (s-M^2)/2 = p.q,  t = (q-qp)^2,  p.qp = x + t/2,  q.qp = -t/2
*   P-gauge: ee.q = ee.P = 0, ep.qp = ep.P = 0   (LMP2015 Eq. 8), polarisations real
*   Gam = T/e^2 (Dirac matrix between ubar(pp) and u(p)), iT = sum of graphs.
*   Projections V_j = Tr[(p/+M) Sbar_j (pp/+M) Gam], summed over the two physical
*   polarisations of each photon in the P-gauge; S_j = structures of hel_from_A (ampl.f90).
*   D = 4 + dd;  dd * (tensor coefficient) -> -2 * (its 1/eps UV coefficient) [rational terms].
#-
Off Statistics;
Symbols D,dd,M,x,t,mpi2,ixs,iys,iPq,PP,sq,X1,X2,X3,dl,ss,F1s,F2s,L0s,L1s;
Symbols n1,n2,n3,n4;
Dimension D;
Vectors p,pp,p2,q,qp,l,ee,ep,k1,k2,k3,P,w1,w2;
Indices i1,i2;
Tensors LT;
CFunctions DDF,KF,PRODK,TCF,UCF;

#procedure TensRed
* loop-momentum tensor decomposition (LoopTools / Denner basis g, k1, k2, k3)
  totensor l,LT;
  id LT(?x) = distrib_(1,0,DDF,KF,?x)+distrib_(1,2,DDF,KF,?x)+distrib_(1,4,DDF,KF,?x)
             +distrib_(1,6,DDF,KF,?x)+distrib_(1,8,DDF,KF,?x);
  id DDF(?x) = dd_(?x)*sq^nargs_(?x);
  id KF(?x) = PRODK(?x);
  repeat;
    id PRODK(w1?,?x) = (k1.w1*X1 + k2.w1*X2 + k3.w1*X3)*PRODK(?x);
    id PRODK(i1?,?x) = (k1(i1)*X1 + k2(i1)*X2 + k3(i1)*X3)*PRODK(?x);
  endrepeat;
  id PRODK = 1;
  multiply TCF(0,0,0,0);
  repeat;
    id sq*TCF(n1?,n2?,n3?,n4?) = TCF(n1+1,n2,n3,n4);
    id X1*TCF(n1?,n2?,n3?,n4?) = TCF(n1,n2+1,n3,n4);
    id X2*TCF(n1?,n2?,n3?,n4?) = TCF(n1,n2,n3+1,n4);
    id X3*TCF(n1?,n2?,n3?,n4?) = TCF(n1,n2,n3,n4+1);
  endrepeat;
* first argument = number of metric indices (2 per "00")
#endprocedure

#procedure Kin
* on-shell kinematics, independent vectors p, q, qp
  id pp = p + q - qp;
  id P = 2*p + q - qp;
  id p.p = M^2;
  id q.q = 0;
  id qp.qp = 0;
  id p.q = x;
  id p.qp = x + t/2;
  id q.qp = -t/2;
#endprocedure

#procedure Gauge
* P-gauge relations for the (real) polarisation vectors
  id ee.q = 0;
  id ee.p = ee.qp/2;
  id ep.qp = 0;
  id ep.p = -ep.q/2;
#endprocedure

#procedure PolSum
* sum over the two physical polarisations, P-gauge:
*   sum ee(a)ee(b) = -g(a,b) + (P(a)q(b)+q(a)P(b))/P.q - P^2 q(a)q(b)/(P.q)^2
  id ee.ee = -2;
  id ee.w1?^2 = -w1.w1 + 2*w1.P*w1.q*iPq - PP*w1.q^2*iPq^2;
  id ee.w1?*ee.w2? = -w1.w2 + (w1.P*w2.q+w1.q*w2.P)*iPq - PP*w1.q*w2.q*iPq^2;
  #call Kin
  id ep.P = 0;
  id ep.qp = 0;
  id ep.p = -ep.q/2;
  .sort
  id ep.ep = -2;
  id ep.w1?^2 = -w1.w1 + 2*w1.P*w1.qp*iPq - PP*w1.qp^2*iPq^2;
  id ep.w1?*ep.w2? = -w1.w2 + (w1.P*w2.qp+w1.qp*w2.P)*iPq - PP*w1.qp*w2.qp*iPq^2;
  #call Kin
* P.q = P.qp = 2x + t/2 ;  P^2 = 4M^2 - t (substituted in Fortran)
#endprocedure

#procedure Proj
* the six projections of the Dirac matrix in the expression GAM (spin line 1; must be defined and sorted)
  L V1 = (g_(1,p)+M*gi_(1)) * (-ep.ee) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V2 = (g_(1,p)+M*gi_(1)) * (q.ep*qp.ee) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V3 = (g_(1,p)+M*gi_(1)) * (-(g_(1,ee,ep)-ee.ep)) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V4 = (g_(1,p)+M*gi_(1)) * (ep.ee*(g_(1,q,qp)-q.qp)) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V5 = (g_(1,p)+M*gi_(1)) * (q.ep*(g_(1,ee,qp)-ee.qp)-qp.ee*(g_(1,ep,q)-ep.q))
                            * (g_(1,pp)+M*gi_(1)) * GAM;
  L V6 = (g_(1,p)+M*gi_(1)) * (q.ep*(g_(1,ee,q)-ee.q)-qp.ee*(g_(1,ep,qp)-ep.qp))
                            * (g_(1,pp)+M*gi_(1)) * GAM;
  .sort
#endprocedure

#procedure Finish
* trace, tensor reduction, kinematics, polarisation sums, D-expansion
  Drop GAM;
  .sort
  tracen,1;
  .sort
  #call TensRed
  .sort
  #call KSUB
  #call Kin
  #call Gauge
  .sort
  #call PolSum
  .sort
  id D = 4 + dd;
  .sort
  id dd^n1?{>1} = 0;
  id dd*TCF(?x) = -2*UCF(?x);
  id dd = 0;
  id M = 1;
  .sort
#endprocedure

#procedure VFinish
* on-shell vertex: u(p) -> u(p2), photon momentum q incoming, q.q = -dl, p.q = dl/2;
* structures R1 = ee/, R2 = P.ee, R3 = q.ee, P = 2p+q; polarisation sum -g (4 dim)
  tracen,1;
  .sort
  #call TensRed
  .sort
  #call KSUB
  id p2 = p + q;
  id P = 2*p + q;
  id p.p = M^2;
  id q.q = -dl;
  id p.q = dl/2;
  .sort
  id ee.ee = -4;
  id ee.w1?^2 = -w1.w1;
  id ee.w1?*ee.w2? = -w1.w2;
  id P = 2*p + q;
  id p.p = M^2;
  id q.q = -dl;
  id p.q = dl/2;
  id D = 4 + dd;
  .sort
  id dd^n1?{>1} = 0;
  id dd*TCF(?x) = -2*UCF(?x);
  id dd = 0;
  .sort
#endprocedure

#procedure VProj
  L V1 = (g_(1,p)+M*gi_(1)) * g_(1,ee) * (g_(1,p2)+M*gi_(1)) * GAM;
  L V2 = (g_(1,p)+M*gi_(1)) * P.ee * (g_(1,p2)+M*gi_(1)) * GAM;
  L V3 = (g_(1,p)+M*gi_(1)) * q.ee * (g_(1,p2)+M*gi_(1)) * GAM;
  .sort
  Drop GAM;
  .sort
#endprocedure

#procedure SFinish
* self-energy L(p) = a p/ + b,  p.p = ss:  V1 = Tr[L]/4 = b,  V2 = Tr[p/ L]/4 = a ss
  L V1 = GAM/4;
  L V2 = g_(1,p)*GAM/4;
  .sort
  Drop GAM;
  .sort
  tracen,1;
  .sort
  #call TensRed
  .sort
  #call KSUB
  id p.p = ss;
  id D = 4 + dd;
  .sort
  id dd^n1?{>1} = 0;
  id dd*TCF(?x) = -2*UCF(?x);
  id dd = 0;
  .sort
#endprocedure

#procedure GProj
* gauge-invariance test: GAM with the incoming polarisation replaced by q; projections on
* R = ep/, ep/q/, ep/qp/, ep/q/qp/, ep.q, ep.q q/, ep.q qp/, ep.q q/qp/  (sum over ep, P-gauge)
  L V1 = (g_(1,p)+M*gi_(1)) * g_(1,ep) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V2 = (g_(1,p)+M*gi_(1)) * g_(1,q,ep) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V3 = (g_(1,p)+M*gi_(1)) * g_(1,qp,ep) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V4 = (g_(1,p)+M*gi_(1)) * g_(1,qp,q,ep) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V5 = (g_(1,p)+M*gi_(1)) * ep.q * (g_(1,pp)+M*gi_(1)) * GAM;
  L V6 = (g_(1,p)+M*gi_(1)) * ep.q*g_(1,q) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V7 = (g_(1,p)+M*gi_(1)) * ep.q*g_(1,qp) * (g_(1,pp)+M*gi_(1)) * GAM;
  L V8 = (g_(1,p)+M*gi_(1)) * ep.q*g_(1,qp,q) * (g_(1,pp)+M*gi_(1)) * GAM;
  .sort
#endprocedure

#procedure GFinish
  Drop GAM;
  .sort
  tracen,1;
  .sort
  #call TensRed
  .sort
  #call KSUB
  #call Kin
  id ep.qp = 0;
  id ep.p = -ep.q/2;
  .sort
  id ep.ep = -2;
  id ep.w1?^2 = -w1.w1 + 2*w1.P*w1.qp*iPq - PP*w1.qp^2*iPq^2;
  id ep.w1?*ep.w2? = -w1.w2 + (w1.P*w2.qp+w1.qp*w2.P)*iPq - PP*w1.qp*w2.qp*iPq^2;
  #call Kin
  id D = 4 + dd;
  .sort
  id dd^n1?{>1} = 0;
  id dd*TCF(?x) = -2*UCF(?x);
  id dd = 0;
  id M = 1;
  .sort
#endprocedure
