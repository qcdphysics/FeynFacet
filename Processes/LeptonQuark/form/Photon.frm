* Photon.frm: |M|^2 for gamma(q) q(k) -> q(p) g(r), the Born process needed
* for the photon-in-lepton subtraction.  SPIN=U: average over the D-2 photon
* polarizations and the quark spin; SPIN=P: circularly polarized photon,
* density +i e_(mu,nu,q,k)/(2 q.k) (FORM's e_ follows the FeynCalc convention
* here, see the gamma5 normalization in common.h), and quark density
* k-slash gamma5/2.  Invariants sh = 2 q.k, th = -2 k.p, uh = -2 q.p.
#-
Off statistics;
#include common.h
Symbols sh,th,uh;
#include amplitudes.inc
.sort

L Amp = Photon1 + Photon2;
L AmpBar = Photon1 + Photon2;
.sort
Drop Born1,Real1,Real2,Vertex1,Photon1,Photon2;
#call conjugate(AmpBar)

L Square = `MINKOWSKISIGN'*Amp*AmpBar/Nc;
.sort
Drop Amp,AmpBar;
#call sew(Square)
id T(aa?,ii?,jj?)*TBar(aa?,jj?,ii?) = CF*Nc;
id Pol(r,mu?)*PolBar(r,nu?) = -d_(mu,nu);
id Den(p+r) = 1/sh;
id Den(k-r) = 1/uh;
#if "`SPIN'" == "U"
* the photon average 1/(D-2) is applied in the notebook (FORM does not divide by polynomials)
  id Pol(q,mu?)*PolBar(q,nu?) = -d_(mu,nu);
#else
* the quark density of sew() is k-slash gamma5/2; the photon density carries its own 1/2
  id Pol(q,mu?)*PolBar(q,nu?) = i_*e_(mu,nu,q,k)/sh;
#endif
id Den(q?) = 1/(q.q);
tracen,2;
#call ProjectEpsilon
id q.q = 0; id k.k = 0; id p.p = 0; id r.r = 0;
id q.k = sh/2; id k.p = -th/2; id q.p = -uh/2; id p.r = sh/2; id k.r = -uh/2; id q.r = -th/2;
.sort
contract;
id q.q = 0; id k.k = 0; id p.p = 0; id r.r = 0;
id q.k = sh/2; id k.p = -th/2; id q.p = -uh/2; id p.r = sh/2; id k.r = -uh/2; id q.r = -th/2;
id uh = -sh-th;
.sort

Print +s Square;
Format Mathematica;
#write <Photon`SPIN'.m> "HoldComplete[%E]", Square
.end
