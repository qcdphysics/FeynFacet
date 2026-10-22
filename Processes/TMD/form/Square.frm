* Square.frm: squared matrix element of one channel of the TMD process,
*
*    Square = NORM * |AMP|^2,  summed over spins, polarizations and colours,
*
* with the open operator index closed as explained in common.h.  The channel
* is selected by preprocessor variables given on the command line (-d) by
* step 2, which also lists them:
*   AMP   the amplitudes to add, e.g. BeamQQ1+BeamQQ2 (names of amplitudes.inc)
*   NORM  spin and colour averages and the operator normalization, e.g. 1/(4*Nc*pn)
*   KIN   beam (physical polarization sums with reference vector n) or
*         soft (a single gluon, summed with -g_{mu nu})
*   OUT   name of the output file (Mathematica format)
* The invariants are left as the symbols pn = n.p, kn = n.k, kb = nb.k and
* pk = p.k; step 2 inserts the on-shell kinematics of each channel.
#-
Off statistics;
#include common.h
#include amplitudes.inc
.sort

L Amp = `AMP';
L AmpBar = `AMP';
.sort
#call conjugate(AmpBar)

L Square = (`NORM')*Amp*AmpBar;
.sort
Drop Amp,AmpBar;
#call sew(Square)
#call polarization(Square)
#call colour(Square)

Skip; NSkip Square;
* eikonal denominators 1/(n.q) with q the momentum along the Wilson line, and
* the quark propagator 1/(k-p)^2 = -1/(2 p.k)
id Eik(k,n) = 1/kn;
id Eik(-p,n) = -1/pn;
id Eik(k,nb) = 1/kb;
id Eik(-k,nb) = -1/kb;
id Den(k-p) = -1/(2*pk);
id Den(-k+p) = -1/(2*pk);
tracen,1;
contract;
id n.n = 0; id nb.nb = 0; id n.nb = 2; id p.p = 0; id k.k = 0; id p.nb = 0;
id p.n = pn; id k.n = kn; id k.nb = kb; id p.k = pk;
.sort

Print +s Square;
Format Mathematica;
#write <`OUT'> "HoldComplete[%E]", Square
.end
