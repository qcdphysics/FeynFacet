* common.h: declarations and procedures shared by the FORM programs of the
* TMD process.  The amplitudes come from amplitudes.inc, written by step 2
* from the FeynArts records; their vocabulary is
*   Chain(Ubar(q),G(mu),Sl(p),...,U(p))   open fermion line ubar(q) gamma^mu p-slash ... u(p);
*                                          Ubar(q) is the operator leg (open Dirac index),
*                                          V(k) the outgoing antiquark of g -> q qbar
*   Pol(k,mu)      polarization vector of the gluon with momentum k; Pol(q,mu) is
*                  the operator leg of the gluon operator (open Lorentz index)
*   Eik(k,n)       1/(n.k), eikonal denominator of a Wilson line
*   Den(q) = 1/q^2        T(a,i,j) = (T^a)_ij     f(a,b,c)     dF(i,j), dA(a,b)
* together with FORM's own d_(mu,nu), p(mu), p.k and i_.

#include declarations.inc
Symbols x,kT2,pn,kn,kb,pk,Nc,CF,CA;
Indices mu,nu,rho,sigma;
* names used as wildcards (ii?, jj?, aa?, ...) in the id statements
Indices ii,jj,aa,bb,cc;
CFunctions Line,LineBar,ChainBar,TBar,PolBar;

* conjugate(E): complex conjugate of the expression E: i -> -i, chains
* reversed, generators transposed, polarization vectors conjugated, every
* Lorentz index renamed to its primed copy (so that the indices of the
* amplitude and of its conjugate stay independent).  f(a,b,c), dF, dA and the
* eikonal denominators are real.
#procedure conjugate(E)
  Skip; NSkip `E';
  id i_ = -i_;
  id Chain(?args) = ChainBar(reverse_(?args));
  id T(aa?,ii?,jj?) = TBar(aa,jj,ii);
  id Pol(p?,mu?) = PolBar(p,mu);
  Multiply replace_(`LORENTZ');
  .sort
#endprocedure

* sew(E): close each fermion chain with its conjugate into one gamma string,
* FORM's g_(1,...): the operator leg ubar(q) ... u(q) is replaced by n-slash
* (the gamma^+ of the operator), the incoming quark by its spin sum p-slash,
* the outgoing antiquark by k-slash.  The q -> g channel has an ordinary
* quark line ubar(k) ... u(p), closed with k-slash and p-slash.
#procedure sew(E)
  Skip; NSkip `E';
  id Chain(Ubar(p?),?x,U(k?)) = Line(p,?x,k);
  id Chain(Ubar(p?),?x,V(k?)) = Line(p,?x,k);
  id ChainBar(U(k?),?x,Ubar(p?)) = LineBar(k,?x,p);
  id ChainBar(V(k?),?x,Ubar(p?)) = LineBar(k,?x,p);
  repeat;
    id Line(?x,G(mu?),?y) = Line(?x,mu,?y);
    id Line(?x,Sl(q?),?y) = Line(?x,q,?y);
    id LineBar(?x,G(mu?),?y) = LineBar(?x,mu,?y);
    id LineBar(?x,Sl(q?),?y) = LineBar(?x,q,?y);
  endrepeat;
  id Line(q,?x,p)*LineBar(p,?y,q) = g_(1,n,?x,p,?y);
  id Line(q,?x,k)*LineBar(k,?y,q) = g_(1,n,?x,k,?y);
  id Line(k,?x,p)*LineBar(p,?y,k) = g_(1,k,?x,p,?y);
  .sort
#endprocedure

* polarization(E): sums over the polarizations of the external gluons with
* the physical (light-cone gauge) sum of reference vector n; the operator
* leg of the gluon operator is contracted with -g_{mu nu}.  For the soft
* function (a single external gluon) -g_{mu nu} is used.
#procedure polarization(E)
  Skip; NSkip `E';
  #if "`KIN'" == "beam"
    id Pol(k,mu?)*PolBar(k,nu?) = -d_(mu,nu) + (k(mu)*n(nu) + n(mu)*k(nu))/kn;
    id Pol(p,mu?)*PolBar(p,nu?) = -d_(mu,nu) + (p(mu)*n(nu) + n(mu)*p(nu))/pn;
    id Pol(q,mu?)*PolBar(q,nu?) = -d_(mu,nu);
  #else
    id Pol(k,mu?)*PolBar(k,nu?) = -d_(mu,nu);
  #endif
  .sort
#endprocedure

* colour(E): the SU(N) sums.  Tr(T^a T^a) = CF Nc, delta_ii = Nc,
* delta^aa = Nc^2-1, f^{abc} f^{abc} = CA (Nc^2-1).
#procedure colour(E)
  Skip; NSkip `E';
  id T(aa?,ii?,jj?)*TBar(aa?,jj?,ii?) = CF*Nc;
  id dF(ii?,jj?)^2 = Nc;
  id dA(aa?,bb?)^2 = Nc^2-1;
  id f(aa?,bb?,cc?)^2 = CA*(Nc^2-1);
  .sort
#endprocedure
