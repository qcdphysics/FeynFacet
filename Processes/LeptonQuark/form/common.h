* common.h: declarations and procedures shared by the FORM programs of this
* process.  The amplitudes come from amplitudes.inc, written by step 2 from the
* FeynArts records; their vocabulary is
*   Chain(Ubar(p),G(mu),Sl(q),...,U(k))  the open fermion line ubar(p) gamma^mu q-slash ... u(k)
*   Den(q) = 1/q^2        Pol(p,mu)      polarization vector of the boson with momentum p
*   T(a,i,j) = (T^a)_ij   dF(i,j)        fundamental delta
* together with FORM's own d_(mu,nu), p(mu), e_(...) and i_.
* The preprocessor variable SPIN (U or P) selects the spin average or the
* helicity difference of the incoming particles.

#include declarations.inc
Symbols s,t,u,a,b,s4,spe3,eps,Nc,CF,dm4;
Symbols d0,d1,d2;
#if "`SPIN'" == "U"
  #define MINKOWSKISIGN "1"
#else
  #define MINKOWSKISIGN "-1"
#endif
Indices mu,nu,rho,sigma;
* Four-dimensional indices, used to define gamma5 in the HVBM scheme
Indices h1=4,h2=4,h3=4,h4=4,e1=4,e2=4,e3=4,e4=4;
* names used as wildcards (ii?, jj?, aa?) in the id statements
Indices ii,jj,aa;
CFunctions Line,LineBar,ChainBar,TBar,PolBar;

* conjugate(E): complex conjugate of the expression E.  Each chain
* Chain(Ubar(p),...,U(k)) becomes ChainBar(U(k),...,Ubar(p)) with the gamma
* string reversed; i -> -i; generators are transposed; polarization vectors
* conjugated; every Lorentz index gets a primed name so that the indices of
* the amplitude and of its conjugate stay independent.
#procedure conjugate(E)
  Skip; NSkip `E';
  id i_ = -i_;
  id Chain(?args) = ChainBar(reverse_(?args));
  id T(aa?,ii?,jj?) = TBar(aa,jj,ii);
  id Pol(p?,mu?) = PolBar(p,mu);
  Multiply replace_(`LORENTZ');
  .sort
#endprocedure

* sew(E): turn the product of a chain and a conjugate chain on the same
* fermion line into one gamma string, FORM's g_(line,...), by inserting the
* spin sum p-slash of the outgoing fermion and the density matrix of the
* incoming one: k-slash/2 for the spin average, k-slash gamma5/2 for the
* helicity difference.  Line 1 is the lepton (outgoing momentum lp), line 2
* the quark (outgoing momentum p).
*
* gamma5 in the HVBM scheme is written as
*   gamma5 = c e_(h1,h2,h3,h4) gamma^h1 gamma^h2 gamma^h3 gamma^h4
* with four-dimensional indices h1..h4, so that it anticommutes with the
* four-dimensional gamma matrices and commutes with the (D-4)-dimensional
* ones (Tracer manual, Appendix A).  The constant c = -i/24 gives Tr[gamma^a gamma^b gamma^c gamma^d gamma5]
* = -4 i e_(a,b,c,d), the convention of FeynCalc's $LeviCivitaSign = -1 that
* the other routes use; Born.frm checks it.  The density matrices carry a
* further factor 1/2 each (hence -i/48 below).
*
* After tracen the programs call ProjectEpsilon (written into declarations.inc
* by step 2): the indices left inside e_ tensors are replaced by
* four-dimensional twins, because FORM would otherwise sum them over D
* dimensions when it contracts two e_ tensors.
*
* FORM contracts two e_ tensors with the Euclidean sign, e_(a,b,c,d)*e_(a,b,c,d)
* = +24, whereas in Minkowski space the result is -24.  A product of two
* gamma5 traces therefore needs an extra factor -1: the programs multiply
* the polarized squares by MINKOWSKISIGN = -1.
#procedure sew(E)
  Skip; NSkip `E';
  id Chain(Ubar(p?),?x,U(k?)) = Line(p,?x,k);
  id ChainBar(U(k?),?x,Ubar(p?)) = LineBar(k,?x,p);
  repeat;
    id Line(?x,G(mu?),?y) = Line(?x,mu,?y);
    id Line(?x,Sl(q?),?y) = Line(?x,q,?y);
    id LineBar(?x,G(mu?),?y) = LineBar(?x,mu,?y);
    id LineBar(?x,Sl(q?),?y) = LineBar(?x,q,?y);
  endrepeat;
  #if "`SPIN'" == "U"
    id Line(lp,?x,l)*LineBar(l,?y,lp) = g_(1,lp,?x,l,?y)/2;
    id Line(p,?x,k)*LineBar(k,?y,p)   = g_(2,p,?x,k,?y)/2;
  #else
    id Line(lp,?x,l)*LineBar(l,?y,lp) = -i_/48*e_(h1,h2,h3,h4)*g_(1,lp,?x,l,h1,h2,h3,h4,?y);
    id Line(p,?x,k)*LineBar(k,?y,p)   = -i_/48*e_(e1,e2,e3,e4)*g_(2,p,?x,k,e1,e2,e3,e4,?y);
  #endif
  .sort
#endprocedure
