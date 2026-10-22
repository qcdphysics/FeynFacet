* common.h: declarations and procedures of Born.frm, the FORM route of the
* partonic Born table.
*
* The amplitudes come from amplitudes.inc, written by 2-Algebra-FORM.wl from
* the FeynArts amplitude records of 1-Born.wl (one local expression Amp1, Amp2, ... per channel, in the
* order of the channel list); their vocabulary is
*   Chain(Ubar(p),G(mu),Sl(q),...,U(k))  the open fermion line ubar(p) gamma^mu q-slash ... u(k);
*                                         Vbar(p) and V(k) are antiquark spinors
*   Den(q) = 1/q^2        Pol(p,mu)      polarization vector of the gluon with momentum p
*   T(i,j,a) = (T^a)_ij   f(a,b,c)       generators and structure constants of
*                                         FORM's colour library color.h
* together with FORM's own d_ (metric and colour deltas), p(mu) and i_.  The
* file declarations.inc (also written by step 2) includes color.h and declares
* the Lorentz indices Lor1, Lor2, ... with dimension D, the quark colours Col1,
* Col2, ... with dimension NR and the gluon colours Glu1, Glu2, ... with
* dimension NA.  It also defines the procedure ProjectEpsilon used in Born.frm.
*
* Kinematics: a(pa) + b(pb) -> c(p) + d(pd), c is the tagged parton,
*   s = 2 pa.pb,  t = -2 pa.p,  u = -2 pb.p,  pd = pa + pb - p,
* all four momenta massless and four-dimensional (HVBM scheme).
*
* The preprocessor variable SPIN selects the spin case of the two incoming
* partons: U sums over their spins, P takes the helicity difference.  The
* spin and colour averages are applied afterwards in the notebook, exactly as
* in the FeynCalc route.

#include declarations.inc
* Mandelstam invariants and the number of colours
Symbols s,t,u,Nc;
* four-dimensional indices for the HVBM gamma5 of the two incoming quarks
Indices h1=4,h2=4,h3=4,h4=4,e1=4,e2=4,e3=4,e4=4;
* names used as wildcards in the id statements (mu?, q?, iW?, ...)
Indices mu,nu;
Indices iW=NR,jW=NR,aW=NA;
Vectors q,q1,q2;
* intermediate objects of the procedures below: a chain with its gamma string
* reversed (conjugate), the conjugate polarization vector (conjugate), a
* closed fermion line before it becomes FORM's g_ (sew, densities), and the
* markers RhoU(q), RhoV(q) of a spin sum u(q) ubar(q) or v(q) vbar(q) (sew)
CFunctions ChainReversed,PolBar,Trace,RhoU,RhoV;

* FORM contracts two e_ tensors with the Euclidean sign, e_(a,b,c,d)*e_(a,b,c,d)
* = +24, whereas in Minkowski space eps_{abcd} eps^{abcd} = -24.  Every term of a
* helicity difference (SPIN=P) contains exactly two e_ tensors (one from each
* polarized incoming parton: the gamma5 of a quark or the circular density of
* a gluon), so the polarized squares are multiplied by MINKOWSKISIGN = -1.
#if "`SPIN'" == "U"
  #define MINKOWSKISIGN "1"
#else
  #define MINKOWSKISIGN "-1"
#endif

*--------------------------------------------------------------------------
* conjugate(E): the complex conjugate of the amplitude E.
*   - i -> -i;
*   - each chain [barred spinor(a), gamma string, spinor(b)] becomes
*     [barred spinor(b), reversed gamma string, spinor(a)]: the Dirac adjoint
*     of ubar(a) gamma^mu1 ... gamma^mun u(b) is ubar(b) gamma^mun ... gamma^mu1 u(a)
*     (no gamma5 occurs in the Born amplitudes);
*   - generators are transposed, (T^a)_ij^* = (T^a)_ji; f and the deltas are real;
*   - polarization vectors are conjugated, Pol -> PolBar;
*   - every Lorentz index, and every colour index summed inside the
*     amplitude, gets a primed name, so that the summed indices of the
*     amplitude and of its conjugate stay independent.  The external colour
*     indices keep their names: multiplying amplitude and conjugate then sums
*     over the colours of all four partons.
*--------------------------------------------------------------------------
#procedure conjugate(E)
  Skip; NSkip `E';
  id i_ = -i_;
  id Chain(?x) = ChainReversed(reverse_(?x));
  id ChainReversed(U(q1?),?x,Ubar(q2?)) = Chain(Ubar(q1),?x,U(q2));
  id ChainReversed(V(q1?),?x,Ubar(q2?)) = Chain(Vbar(q1),?x,U(q2));
  id ChainReversed(U(q1?),?x,Vbar(q2?)) = Chain(Ubar(q1),?x,V(q2));
  id ChainReversed(V(q1?),?x,Vbar(q2?)) = Chain(Vbar(q1),?x,V(q2));
  id T(iW?,jW?,aW?) = T(jW,iW,aW);
  id Pol(q?,mu?) = PolBar(q,mu);
  Multiply replace_(`LORENTZ',`COLOUR');
  .sort
#endprocedure

*--------------------------------------------------------------------------
* sew: close the fermion lines of |M|^2 into Dirac traces.
* In the product of amplitude and conjugate every external quark momentum q
* occurs twice: once as the right end of a chain, u(q) or v(q), and once as
* the left end of another chain, ubar(q) or vbar(q).  The spin sum
* u(q) ubar(q) = q-slash (v(q) vbar(q) = q-slash for massless quarks) joins
* the two chains.  The junction is first marked RhoU(q) or RhoV(q), so that
* the density matrix of a polarized incoming quark can be inserted there
* (procedure densities).  A chain whose two ends carry the same momentum is
* closed into a trace; its last element is the marker of the spin sum that
* closes it.  Depending on the diagrams, the interference has one
* trace (both quark lines run through it) or two.
*--------------------------------------------------------------------------
#procedure sew
  repeat;
    id Chain(?x,U(q?))*Chain(Ubar(q?),?y) = Chain(?x,RhoU(q),?y);
    id Chain(?x,V(q?))*Chain(Vbar(q?),?y) = Chain(?x,RhoV(q),?y);
  endrepeat;
  id Chain(Ubar(q?),?x,U(q?)) = Trace(?x,RhoU(q));
  id Chain(Vbar(q?),?x,V(q?)) = Trace(?x,RhoV(q));
  if ( count(Chain,1) > 0 ) exit "sew: a fermion chain was left open";
#endprocedure

*--------------------------------------------------------------------------
* densities: the spin density matrices of the quarks, and the gluon
* polarization sums.
*
* Quarks.  Outgoing quarks and, for SPIN=U, incoming ones: the spin sum
* q-slash.  For SPIN=P an incoming quark gets q-slash gamma5 and an incoming
* antiquark -q-slash gamma5 (FeynCalc route: v vbar for helicity + is
* (1 - gamma5) q-slash/2); the factors 1/2 are in the averages of the notebook.
* gamma5 in the HVBM scheme is
*   gamma5 = c e_(h1,h2,h3,h4) gamma^h1 gamma^h2 gamma^h3 gamma^h4
* with four-dimensional indices, so that it anticommutes with the
* four-dimensional and commutes with the (D-4)-dimensional gamma matrices.
* c = -i/24 gives Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = -4 i e_(a,b,c,d),
* the convention of FeynCalc's $LeviCivitaSign = -1 (checked by G5Check in
* Born.frm).  pa uses the indices h1..h4, pb the indices e1..e4, so that the
* two gamma5's stay distinct when both quarks run through one trace.
*
* Gluons.  Physical polarization sum with a light-like reference vector r,
*   sum eps^mu eps*^nu = -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r),
* with r = pb for l = pa, and r = pa for l = pb, p, pd (as in the FeynCalc
* route; gauge invariance makes the result independent of r).  The scalar
* products are pa.pb = s/2, p.pa = -t/2, pd.pa = -u/2, which gives the factors
* 2/s, -2/t, -2/u below.  For SPIN=P an incoming gluon l gets the circular
* density (helicity + minus helicity -, with the factor 1/2)
*   eps^mu eps*^nu -> i e_(mu,nu,l,r)/(2 l.r) = i e_(mu,nu,l,r)/s,
* mu being the index of the amplitude (Pol) and nu that of the conjugate
* (PolBar).  The polarized densities are inserted first; the id statements
* after them then act only on the remaining (unpolarized) legs.
*--------------------------------------------------------------------------
#procedure densities
  #if "`SPIN'" == "P"
    id Trace(?x,RhoU(pa),?y) =  -i_/24*e_(h1,h2,h3,h4)*Trace(?x,pa,h1,h2,h3,h4,?y);
    id Trace(?x,RhoV(pa),?y) =   i_/24*e_(h1,h2,h3,h4)*Trace(?x,pa,h1,h2,h3,h4,?y);
    id Trace(?x,RhoU(pb),?y) =  -i_/24*e_(e1,e2,e3,e4)*Trace(?x,pb,e1,e2,e3,e4,?y);
    id Trace(?x,RhoV(pb),?y) =   i_/24*e_(e1,e2,e3,e4)*Trace(?x,pb,e1,e2,e3,e4,?y);
    id Pol(pa,mu?)*PolBar(pa,nu?) = i_*e_(mu,nu,pa,pb)/s;
    id Pol(pb,mu?)*PolBar(pb,nu?) = i_*e_(mu,nu,pb,pa)/s;
  #endif
  repeat;
    id Trace(?x,RhoU(q?),?y) = Trace(?x,q,?y);
    id Trace(?x,RhoV(q?),?y) = Trace(?x,q,?y);
    id Trace(?x,G(mu?),?y) = Trace(?x,mu,?y);
    id Trace(?x,Sl(q?),?y) = Trace(?x,q,?y);
  endrepeat;
  id Pol(pa,mu?)*PolBar(pa,nu?) = -d_(mu,nu) + 2/s*(pa(mu)*pb(nu) + pb(mu)*pa(nu));
  id Pol(pb,mu?)*PolBar(pb,nu?) = -d_(mu,nu) + 2/s*(pb(mu)*pa(nu) + pa(mu)*pb(nu));
  id Pol(p,mu?)*PolBar(p,nu?)   = -d_(mu,nu) - 2/t*(p(mu)*pa(nu) + pa(mu)*p(nu));
  id Pol(pd,mu?)*PolBar(pd,nu?) = -d_(mu,nu) - 2/u*(pd(mu)*pa(nu) + pa(mu)*pd(nu));
  if ( count(Pol,1,PolBar,1) > 0 ) exit "densities: a polarization vector was left";
* the traces become FORM's gamma strings g_(line,...), one line per trace
  id once Trace(?x) = g_(1,?x);
  id once Trace(?x) = g_(2,?x);
  if ( count(Trace,1) > 0 ) exit "densities: more than two traces";
#endprocedure

*--------------------------------------------------------------------------
* kinematics: propagators and scalar products of the 2 -> 2 Born process,
* a(pa) + b(pb) -> c(p) + d(pd) with massless partons.
* propagators: Den(q) = 1/q^2; the propagator momenta are first written with
* pd = pa + pb - p; the invariants are then (pa+pb)^2 = s, (pa-p)^2 = t,
* (pb-p)^2 = u.  Any other denominator stops FORM.
* scalarproducts: p_i^2 = 0, pa.pb = s/2, pa.p = -t/2, pb.p = -u/2.
*--------------------------------------------------------------------------
#procedure propagators
  Argument Den;
    id pd = pa + pb - p;
  EndArgument;
  id Den(pa+pb) = 1/s;
  id Den(-pa-pb) = 1/s;
  id Den(pa-p) = 1/t;
  id Den(p-pa) = 1/t;
  id Den(pb-p) = 1/u;
  id Den(p-pb) = 1/u;
  if ( count(Den,1) > 0 ) exit "propagators: unknown denominator";
#endprocedure

#procedure scalarproducts
  id pd = pa + pb - p;
  id pa.pa = 0;
  id pb.pb = 0;
  id p.p = 0;
  id pa.pb = s/2;
  id pa.p = -t/2;
  id pb.p = -u/2;
#endprocedure
