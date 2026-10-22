* Born.frm: |M|^2 for e(l) q(k) -> e(lp) q(p), spin averaged (SPIN=U) or
* helicity difference (SPIN=P), in D dimensions; the same program as in the
* lepton-quark process (it includes the gamma5 normalization check).
*   form -d SPIN=U Born.frm
#-
Off statistics;
#include common.h
#include amplitudes.inc
.sort

L Amp = Born1;
L AmpBar = Born1;
.sort
Drop Born1,RealQG1,RealQG2,RealGQ1,RealGQ2,PhotonQG1,PhotonQG2,PhotonGQ1,PhotonGQ2;
#call conjugate(AmpBar)

L Square = `MINKOWSKISIGN'*Amp*AmpBar/Nc;
.sort
Drop Amp,AmpBar;
#call sew(Square)
id dF(ii?,jj?)*dF(ii?,jj?) = Nc;
id Den(k-p) = 1/t;
id Den(-k+p) = 1/t;
tracen,1;
tracen,2;
#call ProjectEpsilon
contract;
id l.l = 0; id k.k = 0; id lp.lp = 0; id p.p = 0;
id l.k = s/2; id lp.p = s/2; id l.lp = -t/2; id k.p = -t/2; id l.p = -u/2; id k.lp = -u/2;
id u = -s-t;
.sort

#if "`SPIN'" == "P"
* Check of the gamma5 normalization: Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = -4 i e_(a,b,c,d)
L G5Check = g_(1,h1,h2,h3,h4)*(-i_/24)*e_(e1,e2,e3,e4)*g_(1,e1,e2,e3,e4) + 4*i_*e_(h1,h2,h3,h4);
tracen,1;
contract;
.sort
Print +s G5Check;
#endif

Print +s Square;
Format Mathematica;
#write <Born`SPIN'.m> "HoldComplete[%E]", Square
.end
