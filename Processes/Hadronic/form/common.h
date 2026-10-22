* common.h: declarations and procedures shared by the FORM programs of the
* proton-proton classes (square.frm, checks.frm).
*
* Every FORM job runs in a directory prepared by step 2 (2-Algebra-FORM.wl)
* for one process.  That directory contains
*   process.h       the declarations of the process (its Lorentz and colour
*                   indices, auxiliary vectors, denominator symbols), the list
*                   of its gluons and incoming quarks, and two procedures
*                   written for that process: transversality and auxiliary;
*   amp<i>.inc      diagram i of the process, and
*   bar<j>.inc      the complex conjugate of diagram j (or of the sum of the
*                   Born diagrams, for the one-loop interference).
* The vocabulary of these files:
*   Str(q1, x1, ..., xn, q2)   an open fermion line  ubar(q1) X1 ... Xn u(q2)
*                   (or with v-spinors: for massless partons the spin sums of
*                   u and v spinors are the same, so the spinor type is not
*                   needed); each x is a Lorentz index mu (the matrix gamma^mu)
*                   or a vector q (the matrix q-slash);
*   g_(3, ...)      the Dirac trace of the closed quark loop of a one-loop graph;
*   T(i,j,a,b,..)   (T^a T^b ...)_{ij}; f(a,b,c) the structure constants;
*                   d_(i,j), d_(a,b) the colour deltas (notation of color.h);
*   d_(mu,nu), q(mu)  the metric tensor and the components of a vector;
*   den<n>, dj<n>   symbols for the inverse of a tree propagator, 1/q^2, and of
*                   the j-th loop propagator of the graph's family, 1/D_j
*                   (step 2 puts in their values; FORM handles them as
*                   ordinary symbols);
*   invlr<l>        the symbol 1/(l.r) of the polarization sum of the gluon l
*                   (r its reference vector), also put in by step 2.
* An external gluon with momentum l carries the open index mu<l> in amp<i>.inc
* and nu<l> in bar<j>.inc: its polarization vector has been stripped, and the
* polarization sum (or density) of the gluon is multiplied in after the traces.
*
* Preprocessor variables given on the command line (form -d NAME=value):
*   SPIN   U (spin average, UU) or P (helicity difference, LL);
*   I, J   the diagram numbers of amp<I>.inc and bar<J>.inc;
*   OUT    the name of the result file.

* --- colour: the colour package of J. Vermaseren (color.h) ------------------
* color.h declares T, f, the group invariants NA, NR, cA, cR, I2R, ... and
* sets the default dimension of new indices to NA.  It must come before every
* other declaration.
#include color.h

* --- declarations: Lorentz algebra in D dimensions ---------------------------
* The default dimension is set back to D here: every Lorentz index declared
* from now on (also those of process.h) is D-dimensional.  An index declared
* while Dimension NA is in force would silently be NA-dimensional.
Symbols D;
Dimension D;
Symbols Nc, gs, Nf, s, t, u, A, B, C, spe3;
* the scalar products of the loop momentum (virtual): loop.loop = L2,
* loop.pa = Lpa, loop.pb = Lpb, loop.p = Lp, and Lhat2 = loophat.loophat, the
* square of its (D-4)-dimensional part
Symbols L2, Lpa, Lpb, Lp, Lhat2;
* momenta: pa, pb incoming; p observed; k1, k2 unobserved; loop the loop
* momentum.  k1b, k2b and loopb are the four-dimensional parts of the momenta
* that have (D-4)-dimensional components in the HVBM scheme.
Vectors pa, pb, p, k1, k2, loop;
Vectors k1b, k2b, loopb;
* four-dimensional indices: h1..h3 and e1..e3 for the gamma5 of the incoming
* quarks pa and pb (see sew), z1..z4 for the indices contracted between two
* epsilon tensors (see hvbm)
Indices h1=4, h2=4, h3=4, e1=4, e2=4, e3=4, z1=4, z2=4, z3=4, z4=4;
Set fourdimensional: h1, h2, h3, e1, e2, e3, z1, z2, z3, z4;
* names used only as wildcards (ix?, qq?) in the id statements
Indices ix, ix1, ix2, ix3, ix4;
Vectors qq;
* Str: open fermion line (see above); ClosedLine: a fermion line closed by the
* spin sums, before it becomes a trace g_; Epsilon: an epsilon tensor while its
* arguments are being renamed (see hvbm); LorentzPartA, LorentzPartB: the
* placeholders of the Lorentz parts in the diagram files (see square.frm)
CFunctions Str, ClosedLine, Epsilon, LorentzPartA, LorentzPartB;

#include process.h

* --- the Minkowski sign of two contracted epsilon tensors ---------------------
* In the helicity difference every term contains exactly two epsilon tensors
* (two gamma5's, one gamma5 and one circular gluon density, or two circular
* densities; the check in hvbm makes sure of it).  A product of two epsilon
* tensors does not depend on the sign convention for e^{0123}, but on the
* metric: FORM contracts two epsilon tensors with the Euclidean sign,
* e_(a,b,c,d)*e_(a,b,c,d) = +24, while in Minkowski space the product is -24.
* The polarized squares are therefore multiplied by MINKOWSKISIGN = -1.
#if "`SPIN'" == "U"
  #define MINKOWSKISIGN "1"
#else
  #define MINKOWSKISIGN "-1"
#endif

* --- colour constants ------------------------------------------------------------
* color.h gives the colour factors in terms of the group invariants; for SU(Nc)
* with the generators normalized as Tr(T^a T^b) = delta^{ab}/2:
#procedure sunConstants
  id I2R = 1/2;
  id NA = Nc^2 - 1;
  id NR = Nc;
  id cA = Nc;
  id cR = (Nc^2 - 1)/2*Nc^-1;
  id [cR-cA/2] = -1/2*Nc^-1;
* d33(cOlpR1,cOlpR2) is color.h's name for the invariant d_R^{abc} d_R^{abc} of
* two fundamental traces, where d_R^{abc} = (symmetrized trace of three
* generators) = d^{abc}/4; with d^{abc} d^{abc} = (Nc^2-4)(Nc^2-1)/Nc:
  id d33(cOlpR1,cOlpR2) = (Nc^2-1)*(Nc^2-4)/16*Nc^-1;
#endprocedure

* --- sew: spin sums and densities of the quarks -----------------------------------
* The product of a diagram and a conjugate diagram contains open lines
* Str(q1,...,q2) from both.  Two lines that end and start on the same momentum
* q are joined by inserting the spin sum of that spinor, q-slash (the same for
* u and v spinors of massless partons), until every line is closed: then
* Str(q,...,q) is the trace Tr[... q-slash].
*
* Spin average (SPIN = U): every join is the plain spin sum q-slash; the
* average 1/2 over the spins of an incoming quark is applied by step 2.
*
* Helicity difference (SPIN = P): the density of an incoming quark q is
* sign(q) q-slash gamma5, with sign(q) = +1 for a quark and -1 for an antiquark
* (SIGNa, SIGNb in process.h; the factor 1/2 of the average is again applied
* by step 2).  In the HVBM scheme gamma5 is the four-dimensional product
*    gamma5 = -i/24 e_(a,b,c,d) gamma^a gamma^b gamma^c gamma^d
* with four-dimensional indices a..d, so that it anticommutes with the
* four-dimensional gamma matrices and commutes with the (D-4)-dimensional ones.
* The constant -i/24 gives Tr[gamma^a gamma^b gamma^c gamma^d gamma5]
* = -4 i e_(a,b,c,d), the convention of the FeynCalc and Tracer routes
* (FeynCalc's $LeviCivitaSign = -1); checks.frm verifies it.
* The incoming momentum q is four-dimensional, and the product simplifies to
*    q-slash gamma5 = -i/6 e_(q,a,b,c) gamma^a gamma^b gamma^c,
* the form used below (also verified in checks.frm).
*
* An incoming quark joins a line of the diagram (which ends with u(q), or
* starts with vbar(q) for an antiquark) with a line of the conjugate diagram,
* so its density is inserted at that join.  When both incoming partons are on
* the same fermion line (q qbar -> ...), the join of the first one closes the
* line: the density of the second one is then inserted where the line closes.
* The four statements below cover both cases.
#procedure sew
  #if "`SPIN'" == "P"
    #ifdef `QUARKa'
      id Str(?x,`QUARKa')*Str(`QUARKa',?y) = `SIGNa'*(-i_/6)*e_(`QUARKa',h1,h2,h3)*Str(?x,h1,h2,h3,?y);
    #endif
    #ifdef `QUARKb'
      id Str(?x,`QUARKb')*Str(`QUARKb',?y) = `SIGNb'*(-i_/6)*e_(`QUARKb',e1,e2,e3)*Str(?x,e1,e2,e3,?y);
      id Str(`QUARKb',?x,`QUARKb') = `SIGNb'*(-i_/6)*e_(`QUARKb',e1,e2,e3)*ClosedLine(?x,e1,e2,e3);
    #endif
    #ifdef `QUARKa'
      id Str(`QUARKa',?x,`QUARKa') = `SIGNa'*(-i_/6)*e_(`QUARKa',h1,h2,h3)*ClosedLine(?x,h1,h2,h3);
    #endif
  #endif
* all other joins: the spin sum q-slash
  repeat id Str(?x,qq?)*Str(qq?,?y) = Str(?x,qq,?y);
* every line is now closed: Str(q,...,q) = Tr[... q-slash].  There are at most
* two such traces (FORM's trace lines 1 and 2; line 3 is the quark loop of a
* one-loop graph)
  id Str(qq?,?x,qq?) = ClosedLine(?x,qq);
  id,once ClosedLine(?x) = g_(1,?x);
  id,once ClosedLine(?x) = g_(2,?x);
  if ( count(ClosedLine,1,Str,1) > 0 ) exit "sew: a fermion line was not closed";
#endprocedure

* --- polarization sums and densities of the gluons --------------------------------
* Gluon l with reference vector r = REF<l> (defined in process.h), with mu the
* index of the diagram and nu that of the conjugate diagram:
*   physical sum      -d_(mu,nu) + (l(mu) r(nu) + r(mu) l(nu))/(l.r)
*                     for outgoing gluons, and for incoming gluons when
*                     SPIN = U (the average 1/(D-2) is applied by step 2);
*   circular density  i e_(mu,nu,l,r)/(2 l.r)
*                     for incoming gluons when SPIN = P: the helicity
*                     difference eps_+(mu) eps_+^*(nu) - eps_-(mu) eps_-^*(nu)
*                     times the spin average 1/2, which is therefore included
*                     here (step 2 applies only the colour average).
* The sign of the circular density goes with the gamma5 convention above (both
* are written with the same e_); with both, the LL Born squares agree with the
* table of Processes/PartonicBorn (checked in step 2).
#procedure polarizationSums
  #do l = {`OUTGLUONS'}
    #if "`l'" != ""
      Multiply -d_(mu`l',nu`l') + (`l'(mu`l')*`REF`l''(nu`l') + `REF`l''(mu`l')*`l'(nu`l'))*invlr`l';
    #endif
  #enddo
  #do l = {`INGLUONS'}
    #if "`l'" != ""
      #if "`SPIN'" == "U"
        Multiply -d_(mu`l',nu`l') + (`l'(mu`l')*`REF`l''(nu`l') + `REF`l''(mu`l')*`l'(nu`l'))*invlr`l';
      #else
        Multiply i_/2*e_(mu`l',nu`l',`l',`REF`l'')*invlr`l';
      #endif
    #endif
  #enddo
#endprocedure

* --- hvbm: contraction of the epsilon tensors in four dimensions -------------------
* After the traces and the polarization sums every index is contracted, and
* FORM has absorbed the vectors and metric tensors into the epsilon tensors.
* The epsilon tensors are four-dimensional, so
*  (1) a momentum with (D-4)-dimensional components (HATTED in process.h:
*      k1, k2 in the real emission, loop in the one-loop graphs) stands inside
*      an epsilon tensor for its four-dimensional part k1b, k2b, loopb;
*  (2) an index shared by the two epsilon tensors runs over four dimensions:
*      it is renamed to one of the four-dimensional z1..z4 (FORM would
*      otherwise sum it over D dimensions in contract).
* Then the two tensors are contracted (with the sign MINKOWSKISIGN above,
* applied in square.frm).
#procedure hvbm
* e_ is copied to the function Epsilon, whose arguments can be renamed one by
* one (an index wildcard in an epsilon tensor also matches a vector)
  id e_(ix1?,ix2?,ix3?,ix4?) = Epsilon(ix1,ix2,ix3,ix4);
  #do q = {`HATTED'}
    #if "`q'" != ""
      repeat id Epsilon(?x,`q',?y) = Epsilon(?x,`q'b,?y);
    #endif
  #enddo
  #do z = 1, 4
    id,once Epsilon(?x,ix?!fourdimensional,?y)*Epsilon(?v,ix?,?w) = Epsilon(?x,z`z',?y)*Epsilon(?v,z`z',?w);
  #enddo
  if ( ( count(Epsilon,1) != 0 ) && ( count(Epsilon,1) != 2 ) ) exit "hvbm: a term without exactly two epsilon tensors";
  id Epsilon(?x) = e_(?x);
  contract;
#endprocedure

* --- kinematics ----------------------------------------------------------------------
* Born (and one loop): s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2, k1 = pa+pb-p,
* all four-dimensional.
#procedure bornKinematics
  id pa.pa = 0; id pb.pb = 0; id p.p = 0; id k1.k1 = 0;
  id pa.pb = s/2; id pa.p = -t/2; id pb.p = -u/2;
  id pa.k1 = -u/2; id pb.k1 = -t/2; id p.k1 = s/2;
#endprocedure
* the loop momentum: its D-dimensional products, and the products of its
* four-dimensional part loopb (loopb.loopb = loop.loop - loophat.loophat; the
* external momenta are four-dimensional, so loopb.pa = loop.pa, ...)
#procedure loopProducts
  id loop.loop = L2; id loop.pa = Lpa; id loop.pb = Lpb; id loop.p = Lp; id loop.k1 = Lpa + Lpb - Lp;
  id loopb.loopb = L2 - Lhat2; id loopb.pa = Lpa; id loopb.pb = Lpb; id loopb.p = Lp; id loopb.k1 = Lpa + Lpb - Lp;
#endprocedure
* real emission: A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2, s4 = (k1+k2)^2 =
* s+t+u, k2 = pa+pb-p-k1.  The (D-4)-dimensional parts of k1 and k2 enter the
* four-dimensional products through spe3 = k1hat^2 = k2hat^2 = -k1hat.k2hat
* (k1hat + k2hat = 0, since pa, pb, p are four-dimensional); for example
* k1b.k1b = k1.k1 - k1hat.k1hat = -spe3.
#procedure realKinematics
  id pa.pa = 0; id pb.pb = 0; id p.p = 0; id k1.k1 = 0; id k2.k2 = 0;
  id pa.pb = s/2; id pa.p = -t/2; id pb.p = -u/2;
  id pa.k1 = -A/2; id pb.k1 = -B/2; id p.k1 = C/2;
  id pa.k2 = (s+t+A)/2; id pb.k2 = (s+u+B)/2; id p.k2 = (-t-u-C)/2; id k1.k2 = (s+t+u)/2;
  id k1b.k1b = -spe3; id k2b.k2b = -spe3; id k1b.k2b = (s+t+u)/2 + spe3;
  id pa.k1b = -A/2; id pb.k1b = -B/2; id p.k1b = C/2;
  id pa.k2b = (s+t+A)/2; id pb.k2b = (s+u+B)/2; id p.k2b = (-t-u-C)/2;
#endprocedure
