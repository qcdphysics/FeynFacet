* RealGQ.frm: |M|^2 for e(l) g(k) -> e(lp) q(p) qbar(r), the quark observed,
* spin averaged (SPIN=U) or helicity difference (SPIN=P).  SPIN=U: the gluon
* polarizations are summed with -g_{mu nu} (the average 1/(D-2) is applied in
* the notebook); SPIN=P: circular gluon density i e_(mu,nu,k,l)/(2 k.l).
* Colour: Tr(T^a T^a)/(Nc^2-1) = TR, the colour average included.
* Kinematics as in RealQG.frm.
#-
Off statistics;
#include common.h
#include amplitudes.inc
.sort

L Amp = RealGQ1 + RealGQ2;
L AmpBar = RealGQ1 + RealGQ2;
.sort
Drop Born1,RealQG1,RealQG2,RealGQ1,RealGQ2,PhotonQG1,PhotonQG2,PhotonGQ1,PhotonGQ2;
#call conjugate(AmpBar)

L Square = `MINKOWSKISIGN'*Amp*AmpBar;
.sort
Drop Amp,AmpBar;
#call sew(Square)
id T(aa?,ii?,jj?)*TBar(aa?,jj?,ii?) = TR;
#if "`SPIN'" == "U"
  id Pol(k,mu?)*PolBar(k,nu?) = -d_(mu,nu);
#else
  id Pol(k,mu?)*PolBar(k,nu?) = i_*e_(mu,nu,k,l)/s;
#endif
id Den(k-r) = 1/a;
id Den(-k+r) = 1/a;
id Den(k-p) = 1/t;
id Den(-k+p) = 1/t;
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
#write <RealGQ`SPIN'.m> "HoldComplete[%E]", Square
.end
