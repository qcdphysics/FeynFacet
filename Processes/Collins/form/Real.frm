* Real.frm: |M|^2 for e(l) q(k) -> e(lp) q(p) g(r): spin average (SPIN=U) or
* transverse spin transfer (SPIN=NN, SPIN=XX).  Gluon sum -g_{mu nu}, colour
* Tr(T^a T^a)/Nc = CF.  Kinematics as in the lepton-quark process: a = (k-r)^2,
* b = (p+r)^2, tau = a + b + t, s4 = s + t + u; lp and r have hatted parts
* (spe3); the normal component nr = nn.r = -nn.lp is reduced in the notebook.
#-
Off statistics;
#include common.h
#include amplitudes.inc
.sort

L Amp = Real1 + Real2;
L AmpBar = Real1 + Real2;
.sort
Drop Born1,Real1,Real2,Vertex1,Photon1,Photon2;
#call conjugate(AmpBar)

L Square = `MINKOWSKISIGN'*Amp*AmpBar/Nc;
.sort
Drop Amp,AmpBar;
#call sew(Square)
id T(aa?,ii?,jj?)*TBar(aa?,jj?,ii?) = CF*Nc;
id Pol(r,mu?)*PolBar(r,nu?) = -d_(mu,nu);
id Den(k-r) = 1/a;
id Den(p+r) = 1/b;
id Den(k-p-r) = 1/(a+b+t);
tracen,1;
tracen,2;
#call ProjectEpsilon
#call axes(l,k,p,s,t,u)
* First pass: D-dimensional products
#call normal(l,k,p)
id nn.r = nr; id nn.lp = -nr;
id l.l = 0; id k.k = 0; id lp.lp = 0; id p.p = 0; id r.r = 0;
id l.k = s/2; id k.p = -t/2; id l.p = -u/2; id k.r = -a/2; id p.r = b/2;
id l.lp = -(a+b+t)/2; id l.r = (s+t+u+a+b)/2; id k.lp = (s+t+a)/2;
id lp.p = -(t+u+b)/2; id lp.r = (s+t+u)/2;
.sort
contract;
#call axes(l,k,p,s,t,u)
* Second pass: four-dimensional products from the contracted epsilon tensors
#call normal(l,k,p)
id nn.r = nr; id nn.lp = -nr;
id r.r = -spe3; id lp.lp = -spe3; id lp.r = (s+t+u)/2 + spe3;
id l.l = 0; id k.k = 0; id p.p = 0;
id l.k = s/2; id k.p = -t/2; id l.p = -u/2; id k.r = -a/2; id p.r = b/2;
id l.lp = -(a+b+t)/2; id l.r = (s+t+u+a+b)/2; id k.lp = (s+t+a)/2; id lp.p = -(t+u+b)/2;
.sort

Print +s Square;
Format Mathematica;
#write <Real`SPIN'.m> "HoldComplete[%E]", Square
.end
