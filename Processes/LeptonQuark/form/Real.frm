* Real.frm: |M|^2 for e(l) q(k) -> e(lp) q(p) g(r), spin averaged (SPIN=U) or
* helicity difference (SPIN=P).  The gluon polarizations are summed with
* -g_{mu nu}; the external quark colour is averaged.
*
* Kinematics: a = (k-r)^2 = -2 k.r, b = (p+r)^2 = 2 p.r, tau = (k-p-r)^2 =
* a + b + t, s4 = (lp+r)^2 = s + t + u.  In the HVBM scheme l, k, p are
* four-dimensional while lp and r have (D-4)-dimensional parts with
* rhat.rhat = lphat.lphat = -lphat.rhat = spe3.
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
* colour: sum_a (T^a)_ij (T^a)_ji = Tr(T^a T^a) = CF Nc
id T(aa?,ii?,jj?)*TBar(aa?,jj?,ii?) = CF*Nc;
id Pol(r,mu?)*PolBar(r,nu?) = -d_(mu,nu);
id Den(k-r) = 1/a;
id Den(p+r) = 1/b;
id Den(k-p-r) = 1/(a+b+t);
tracen,1;
tracen,2;
#call ProjectEpsilon

* First pass, before the epsilon tensors are contracted: every scalar product
* produced here is D-dimensional.
id l.l = 0; id k.k = 0; id lp.lp = 0; id p.p = 0; id r.r = 0;
id l.k = s/2; id k.p = -t/2; id l.p = -u/2; id k.r = -a/2; id p.r = b/2;
id l.lp = -(a+b+t)/2; id l.r = (s+t+u+a+b)/2; id k.lp = (s+t+a)/2;
id lp.p = -(t+u+b)/2; id lp.r = (s+t+u)/2;
.sort
contract;
* Second pass: products created by contracting two four-dimensional epsilon
* tensors are four-dimensional ("barred").  For the physical momenta l, k, p
* they coincide with the D-dimensional ones; for lp and r the hatted parts are
* subtracted: rbar.rbar = r.r - rhat.rhat = -spe3, lpbar.rbar = lp.r + spe3.
id r.r = -spe3; id lp.lp = -spe3; id lp.r = (s+t+u)/2 + spe3;
id l.l = 0; id k.k = 0; id p.p = 0;
id l.k = s/2; id k.p = -t/2; id l.p = -u/2; id k.r = -a/2; id p.r = b/2;
id l.lp = -(a+b+t)/2; id l.r = (s+t+u+a+b)/2; id k.lp = (s+t+a)/2; id lp.p = -(t+u+b)/2;
.sort

Print +s Square;
Format Mathematica;
#write <Real`SPIN'.m> "HoldComplete[%E]", Square
.end
