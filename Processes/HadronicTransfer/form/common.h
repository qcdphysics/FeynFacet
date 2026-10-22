* common.h: declarations and procedures shared by the FORM programs of the
* transverse spin transfer classes p^ p -> h^ + X (square.frm, checks.frm).
* It is the common.h of Processes/HadronicTT/form with one spin vector on an
* incoming line and one on an outgoing line: the density pa-slash gamma5
* sa-slash of the polarized incoming parton a and the projection p-slash gamma5
* sc-slash of the observed parton c, in place of the two incoming densities of
* the A_TT classes.  Read Processes/Hadronic/form/common.h first: the
* vocabulary of the process files is the same.
*
* Every FORM job runs in a directory prepared by step 2 (2-Algebra-FORM.wl)
* for one process.  That directory contains
*   process.h       the declarations of the process (its Lorentz and colour
*                   indices, auxiliary vectors, denominator symbols), the list
*                   of its gluons with their reference vectors, and two
*                   procedures written for that process: transversality and
*                   auxiliary;
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
* polarization sum of the gluon is multiplied in after the traces.
*
* Preprocessor variables given on the command line (form -d NAME=value):
*   SPIN   the component of the transverse spin transfer: XX, NN, XN or NX
*          (the spin vector of the incoming parton a along the in-plane axis
*          eX or along the normal nn, then that of the observed parton c along
*          its in-plane axis eC or along nn);
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
* the normal components: nk1 = nn.k1 = -nn.k2 of the unobserved momenta of the
* real emission, nloop = nn.loop of the loop momentum (step 2 reduces their
* even powers to invariants and drops the odd ones)
Symbols nk1, nloop;
* momenta: pa, pb incoming; p observed; k1, k2 unobserved; loop the loop
* momentum.  k1b, k2b and loopb are the four-dimensional parts of the momenta
* that have (D-4)-dimensional components in the HVBM scheme.
Vectors pa, pb, p, k1, k2, loop;
Vectors k1b, k2b, loopb;
* the spin vectors: sa of the polarized incoming parton a (momentum pa) and
* sc of the observed parton c (momentum p), kept symbolic through the traces,
* and the axes of the scattering plane they are set to afterwards: the
* in-plane axis of the incoming partons eXs = s eX = u pa + t pb + s p
* (orthogonal to pa and pb, along the transverse momentum of p), the in-plane
* axis of the observed parton eCs = (t+u)^2 eC = -u(t+u) pa + t(t+u) pb
* - s(u-t) p (orthogonal to p and to pa + pb), both scaled to polynomial
* coefficients (FORM has no polynomial denominators; step 2 divides the XX
* component by s (t+u)^2 and by the lengths of the axes), and the unit normal
* nn (nn.nn = -1, orthogonal to pa, pb, p).  All five are four-dimensional.
Vectors sa, sc, eXs, eCs, nn;
* four-dimensional indices: h1..h3 for the gamma5 of the incoming parton pa
* and e1..e3 for the gamma5 of the observed parton p (see sew), z1..z4 for the
* indices contracted between two epsilon tensors (see hvbm)
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
* Every term of a transfer square contains exactly two epsilon tensors, one
* from the gamma5 of the density of a and one from the gamma5 of the
* projection of c (the check in hvbm makes sure of it).  A product of two
* epsilon tensors does not depend on the sign convention for e^{0123}, but on
* the metric: FORM contracts two epsilon tensors with the Euclidean sign,
* e_(a,b,c,d)*e_(a,b,c,d) = +24, while in Minkowski space the product is -24.
* Every square is therefore multiplied by MINKOWSKISIGN = -1 (as the polarized
* squares of Processes/Hadronic and Processes/HadronicTT).
#define MINKOWSKISIGN "-1"
#if "`SPIN'" == "XX"
#elseif "`SPIN'" == "NN"
#elseif "`SPIN'" == "XN"
#elseif "`SPIN'" == "NX"
#else
  #message SPIN must be XX, NN, XN or NX
  #terminate
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

* --- sew: spin sums, the density of a and the projection of c -------------------
* The product of a diagram and a conjugate diagram contains open lines
* Str(q1,...,q2) from both.  Two lines that end and start on the same momentum
* q are joined by inserting the spin sum of that spinor, q-slash (the same for
* u and v spinors of massless partons), until every line is closed: then
* Str(q,...,q) is the trace Tr[... q-slash].
*
* Two joins carry the spin projections (docs/TRANSVERSE-PP.md):
*   at pa, the polarized incoming quark or antiquark a, the density
*        pa-slash gamma5 sa-slash      (the same expression for a quark and
*        for an antiquark; its factor 1/2 is applied by step 2 together with
*        the colour average);
*   at p, the observed quark or antiquark c, the projection
*        p-slash gamma5 sc-slash       (no 1/2: the partonic quantity is the
*        difference sigma(sc) - sigma(-sc)).
* In the HVBM scheme gamma5 is the four-dimensional product
*    gamma5 = -i/24 e_(a,b,c,d) gamma^a gamma^b gamma^c gamma^d
* with four-dimensional indices a..d, so that it anticommutes with the
* four-dimensional gamma matrices and commutes with the (D-4)-dimensional ones.
* The constant -i/24 gives Tr[gamma^a gamma^b gamma^c gamma^d gamma5]
* = -4 i e_(a,b,c,d), the convention of the FeynCalc route (FeynCalc's
* $LeviCivitaSign = -1); checks.frm verifies it.  The momenta pa and p are
* four-dimensional, and the product simplifies to
*    q-slash gamma5 = -i/6 e_(q,a,b,c) gamma^a gamma^b gamma^c,
* the form used below (also verified in checks.frm); the spin vector follows
* as the matrix s-slash in the string.
*
* The incoming parton a joins a line of the diagram (which ends with u(pa),
* or starts with vbar(pa) for an antiquark) with a line of the conjugate
* diagram; the observed parton c joins a line of the diagram (which starts
* with ubar(p), or ends with v(p) for an antiquark) with a line of the
* conjugate diagram.  When a and c are on the same fermion line (the exchange
* graphs), the first join closes the line: the projection of the second one
* is then inserted where the line closes.  The four statements below cover
* both cases; the order of the matrices in the trace is that of the FeynCalc
* route (... pa-slash gamma5 sa-slash ... p-slash gamma5 sc-slash ...).
#procedure sew
  id Str(?x,pa)*Str(pa,?y) = (-i_/6)*e_(pa,h1,h2,h3)*Str(?x,h1,h2,h3,sa,?y);
  id Str(?x,p)*Str(p,?y) = (-i_/6)*e_(p,e1,e2,e3)*Str(?x,e1,e2,e3,sc,?y);
  id Str(p,?x,p) = (-i_/6)*e_(p,e1,e2,e3)*ClosedLine(?x,e1,e2,e3,sc);
  id Str(pa,?x,pa) = (-i_/6)*e_(pa,h1,h2,h3)*ClosedLine(?x,h1,h2,h3,sa);
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

* --- polarization sums of the gluons ----------------------------------------------
* Gluon l with reference vector r = REF<l> (defined in process.h), with mu the
* index of the diagram and nu that of the conjugate diagram: the physical sum
*     -d_(mu,nu) + (l(mu) r(nu) + r(mu) l(nu))/(l.r)
* for the outgoing gluons and for the incoming gluon b of the q g classes,
* which is unpolarized (its average 1/(D-2) is applied by step 2).  The
* polarized parton a is never a gluon (a gluon carries no transversity).
#procedure polarizationSums
  #do l = {`OUTGLUONS'}
    #if "`l'" != ""
      Multiply -d_(mu`l',nu`l') + (`l'(mu`l')*`REF`l''(nu`l') + `REF`l''(mu`l')*`l'(nu`l'))*invlr`l';
    #endif
  #enddo
  #do l = {`INGLUONS'}
    #if "`l'" != ""
      #if "`l'" == "pa"
        #message the polarized parton a cannot be a gluon
        #terminate
      #endif
      Multiply -d_(mu`l',nu`l') + (`l'(mu`l')*`REF`l''(nu`l') + `REF`l''(mu`l')*`l'(nu`l'))*invlr`l';
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
* The spin vectors sa, sc (later eXs, eCs, nn) are four-dimensional and need
* no renaming.  Then the two tensors are contracted (with the sign
* MINKOWSKISIGN above, applied in square.frm).
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
  if ( count(Epsilon,1) != 2 ) exit "hvbm: a term without exactly two epsilon tensors";
  id Epsilon(?x) = e_(?x);
  contract;
#endprocedure

* --- spinVectors: the component SPIN of the transverse spin transfer ---------------
* After the contraction the spin vectors occur only in scalar products.  The
* component is selected by setting sa to the in-plane axis eXs or to the
* normal nn, and sc to the in-plane axis eCs or to nn (SPIN = XX: eXs, eCs;
* NN: nn, nn; XN: sa = eXs, sc = nn; NX: sa = nn, sc = eCs).  Then
*  - the products of the axes follow from their definitions
*        eXs = u pa + t pb + s p,   eCs = -u(t+u) pa + t(t+u) pb - s(u-t) p
*    (so eXs.pa = eXs.pb = 0, eCs.p = 0, eCs.(pa+pb) = 0, nn.eXs = nn.eCs = 0;
*    their products with k1, k2 and loop become products of pa, pb, p with
*    those momenta); step 2 divides the XX component by the scales s and
*    (t+u)^2 and by the product of the lengths of eX and eC, -t u/(t+u), in all
*    by -s t u (t+u), and leaves XN and NX unnormalized (they must vanish);
*  - the normal is orthogonal to pa, pb, p, with nn.nn = -1; its products with
*    the unobserved momenta of the real emission are nn.k1 = -nn.k2 = nk1 (in
*    the Born kinematics k1 = pa + pb - p lies in the plane and nn.k1 = 0), and
*    with the loop momentum nn.loop = nloop.  The four-dimensional parts k1b,
*    k2b, loopb have the same products with the four-dimensional nn and axes.
#procedure spinVectors
  #if "`SPIN'" == "XX"
    Multiply replace_(sa,eXs,sc,eCs);
  #elseif "`SPIN'" == "NN"
    Multiply replace_(sa,nn,sc,nn);
  #elseif "`SPIN'" == "XN"
    Multiply replace_(sa,eXs,sc,nn);
  #else
    Multiply replace_(sa,nn,sc,eCs);
  #endif
  id eXs.eCs = u*pa.eCs + t*pb.eCs + s*p.eCs;
  id eXs.eXs = u*pa.eXs + t*pb.eXs + s*p.eXs;
  id eXs.qq? = u*pa.qq + t*pb.qq + s*p.qq;
  id eCs.eCs = -u*(t+u)*pa.eCs + t*(t+u)*pb.eCs - s*(u-t)*p.eCs;
  id eCs.qq? = -u*(t+u)*pa.qq + t*(t+u)*pb.qq - s*(u-t)*p.qq;
  id nn.nn = -1;
  id nn.pa = 0; id nn.pb = 0; id nn.p = 0;
  #if "`KINEMATICS'" == "Real"
    id nn.k1 = nk1; id nn.k2 = -nk1; id nn.k1b = nk1; id nn.k2b = -nk1;
  #else
    id nn.k1 = 0;
    id nn.loop = nloop; id nn.loopb = nloop;
  #endif
  if ( count(sa,1,sc,1,eXs,1,eCs,1,nn,1) > 0 ) exit "spinVectors: a spin vector was left";
#endprocedure

* --- kinematics ----------------------------------------------------------------------
* Born (and one loop): s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2, k1 = pa+pb-p,
* all four-dimensional (u = -s-t is put in by step 2, at the end).
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
