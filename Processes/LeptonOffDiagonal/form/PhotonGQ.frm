* PhotonGQ.frm: |M|^2 for gamma(q) g(k) -> q(p) qbar(r), the photon Born
* process of the photon-in-lepton subtraction of l g -> q + X.  SPIN=U:
* -g_{mu nu} for the photon and the gluon (the averages 1/(D-2)^2 are applied
* in the notebook); SPIN=P: circular densities i e_(mu,nu,q,k)/(2 q.k) for the
* photon and i e_(rho,sigma,k,q)/(2 k.q) for the gluon.  Colour TR as in
* RealGQ.frm.  Invariants sh = 2 q.k, th = -2 k.p, uh = -2 q.p.
#-
Off statistics;
#include common.h
Symbols sh,th,uh;
CFunction den;
#include amplitudes.inc
.sort

L Amp = PhotonGQ1 + PhotonGQ2;
L AmpBar = PhotonGQ1 + PhotonGQ2;
.sort
Drop Born1,RealQG1,RealQG2,RealGQ1,RealGQ2,PhotonQG1,PhotonQG2,PhotonGQ1,PhotonGQ2;
#call conjugate(AmpBar)

L Square = `MINKOWSKISIGN'*Amp*AmpBar;
.sort
Drop Amp,AmpBar;
#call sew(Square)
* every propagator 1/q^2 becomes den(q.q); the kinematics inside the argument then gives den(sh), den(th) or den(uh)
id Den(q?) = den(q.q);
id T(aa?,ii?,jj?)*TBar(aa?,jj?,ii?) = TR;
#if "`SPIN'" == "U"
  id Pol(q,mu?)*PolBar(q,nu?) = -d_(mu,nu);
  id Pol(k,mu?)*PolBar(k,nu?) = -d_(mu,nu);
#else
  id Pol(q,mu?)*PolBar(q,nu?) = i_*e_(mu,nu,q,k)/sh;
  id Pol(k,mu?)*PolBar(k,nu?) = i_*e_(mu,nu,k,q)/sh;
#endif
tracen,2;
#call ProjectEpsilon
#call photonKinematics
argument den;
  #call photonKinematics
endargument;
id den(th) = 1/th;
id den(uh) = 1/uh;
.sort
contract;
#call photonKinematics
id uh = -sh-th;
.sort

Print +s Square;
Format Mathematica;
#write <PhotonGQ`SPIN'.m> "HoldComplete[%E]", Square
.end
