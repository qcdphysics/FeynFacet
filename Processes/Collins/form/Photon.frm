* Photon.frm: |M|^2 for gamma(q) q(k) -> q(p) g(r), the photon Born process
* of the photon-in-lepton subtraction: spin average (SPIN=U) or transverse spin
* transfer (SPIN=NN, SPIN=XX); the photon is unpolarized, -g_{mu nu} (the
* average 1/(D-2) is applied in the notebook).  Invariants sh = 2 q.k,
* th = -2 k.p, uh = -2 q.p; the axes use q in place of l.
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
id Pol(q,mu?)*PolBar(q,nu?) = -d_(mu,nu);
id Den(p+r) = 1/sh;
id Den(k-r) = 1/uh;
id Den(q?) = 1/(q.q);
tracen,2;
#call ProjectEpsilon
#call axes(q,k,p,sh,th,uh)
#call normal(q,k,p)
id nn.r = 0;
id q.q = 0; id k.k = 0; id p.p = 0; id r.r = 0;
id q.k = sh/2; id k.p = -th/2; id q.p = -uh/2; id p.r = sh/2; id k.r = -uh/2; id q.r = -th/2;
.sort
contract;
#call axes(q,k,p,sh,th,uh)
#call normal(q,k,p)
id nn.r = 0;
id q.q = 0; id k.k = 0; id p.p = 0; id r.r = 0;
id q.k = sh/2; id k.p = -th/2; id q.p = -uh/2; id p.r = sh/2; id k.r = -uh/2; id q.r = -th/2;
id uh = -sh-th;
.sort

Print +s Square;
Format Mathematica;
#write <Photon`SPIN'.m> "HoldComplete[%E]", Square
.end
