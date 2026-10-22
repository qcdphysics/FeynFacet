* Born.frm: |M|^2 for e(l) q(k) -> e(lp) q(p) in D dimensions: spin average
* (SPIN=U) or transverse spin transfer (SPIN=NN, SPIN=XX).
*   form -d SPIN=NN Born.frm
#-
Off statistics;
#include common.h
#include amplitudes.inc
.sort

L Amp = Born1;
L AmpBar = Born1;
.sort
Drop Born1,Real1,Real2,Vertex1,Photon1,Photon2;
#call conjugate(AmpBar)

L Square = `MINKOWSKISIGN'*Amp*AmpBar/Nc;
.sort
Drop Amp,AmpBar;
#call sew(Square)
id dF(ii?,jj?)*dF(ii?,jj?) = Nc;
id Den(k-p) = 1/t;
tracen,1;
tracen,2;
#call ProjectEpsilon
#call axes(l,k,p,s,t,u)
#call normal(l,k,p)
id nn.lp = 0;
id l.l = 0; id k.k = 0; id lp.lp = 0; id p.p = 0;
id l.k = s/2; id lp.p = s/2; id l.lp = -t/2; id k.p = -t/2; id l.p = -u/2; id k.lp = -u/2;
.sort
contract;
#call axes(l,k,p,s,t,u)
#call normal(l,k,p)
id nn.lp = 0;
id l.l = 0; id k.k = 0; id lp.lp = 0; id p.p = 0;
id l.k = s/2; id lp.p = s/2; id l.lp = -t/2; id k.p = -t/2; id l.p = -u/2; id k.lp = -u/2;
id u = -s-t;
.sort

Print +s Square;
Format Mathematica;
#write <Born`SPIN'.m> "HoldComplete[%E]", Square
.end
