* PhotonQG.frm: |M|^2 for gamma(q) q(k) -> g(p) q(r), the photon Born process
* of the photon-in-lepton subtraction of l q -> g + X.  SPIN=U: -g_{mu nu}
* for the photon (the average 1/(D-2) is applied in the notebook) and the quark
* spin average; SPIN=P: circular photon density i e_(mu,nu,q,k)/(2 q.k) and
* quark density k-slash gamma5/2.  Invariants sh = 2 q.k, th = -2 k.p, uh = -2 q.p.
#-
Off statistics;
#include common.h
Symbols sh,th,uh;
CFunction den;
#include amplitudes.inc
.sort

L Amp = PhotonQG1 + PhotonQG2;
L AmpBar = PhotonQG1 + PhotonQG2;
.sort
Drop Born1,RealQG1,RealQG2,RealGQ1,RealGQ2,PhotonQG1,PhotonQG2,PhotonGQ1,PhotonGQ2;
#call conjugate(AmpBar)

L Square = `MINKOWSKISIGN'*Amp*AmpBar/Nc;
.sort
Drop Amp,AmpBar;
#call sew(Square)
* every propagator 1/q^2 becomes den(q.q); the kinematics inside the argument then gives den(sh), den(th) or den(uh)
id Den(q?) = den(q.q);
id T(aa?,ii?,jj?)*TBar(aa?,jj?,ii?) = CF*Nc;
id Pol(p,mu?)*PolBar(p,nu?) = -d_(mu,nu);
#if "`SPIN'" == "U"
  id Pol(q,mu?)*PolBar(q,nu?) = -d_(mu,nu);
#else
  id Pol(q,mu?)*PolBar(q,nu?) = i_*e_(mu,nu,q,k)/sh;
#endif
tracen,2;
#call ProjectEpsilon
#call photonKinematics
argument den;
  #call photonKinematics
endargument;
id den(sh) = 1/sh;
id den(th) = 1/th;
.sort
contract;
#call photonKinematics
id uh = -sh-th;
.sort

Print +s Square;
Format Mathematica;
#write <PhotonQG`SPIN'.m> "HoldComplete[%E]", Square
.end
