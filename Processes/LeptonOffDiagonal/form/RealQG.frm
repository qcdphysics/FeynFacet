* RealQG.frm: |M|^2 for e(l) q(k) -> e(lp) g(p) q(r), the gluon observed,
* spin averaged (SPIN=U) or helicity difference (SPIN=P).  The gluon
* polarizations are summed with -g_{mu nu}; the quark colour is averaged:
* Tr(T^a T^a)/Nc = CF.
* Kinematics: a = (k-r)^2, b = (p+r)^2, tau = (k-p-r)^2 = a + b + t,
* s4 = (lp+r)^2 = s + t + u; lp and r have hatted parts, spe3 = rhat.rhat.
#-
Off statistics;
#include common.h
#include amplitudes.inc
.sort

L Amp = RealQG1 + RealQG2;
L AmpBar = RealQG1 + RealQG2;
.sort
Drop Born1,RealQG1,RealQG2,RealGQ1,RealGQ2,PhotonQG1,PhotonQG2,PhotonGQ1,PhotonGQ2;
#call conjugate(AmpBar)

L Square = `MINKOWSKISIGN'*Amp*AmpBar/Nc;
.sort
Drop Amp,AmpBar;
#call sew(Square)
id T(aa?,ii?,jj?)*TBar(aa?,jj?,ii?) = CF*Nc;
id Pol(p,mu?)*PolBar(p,nu?) = -d_(mu,nu);
id Den(k-p) = 1/t;
id Den(-k+p) = 1/t;
id Den(p+r) = 1/b;
id Den(k-p-r) = 1/(a+b+t);
id Den(-k+p+r) = 1/(a+b+t);
tracen,1;
tracen,2;
#call ProjectEpsilon
#call realKinematics
.sort
contract;
#call realKinematicsBarred
.sort

Print +s Square;
Format Mathematica;
#write <RealQG`SPIN'.m> "HoldComplete[%E]", Square
.end
