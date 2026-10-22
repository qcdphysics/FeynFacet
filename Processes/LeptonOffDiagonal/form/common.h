* common.h: declarations and procedures shared by the FORM programs of this
* process (the lepton-quark process has the same structure; read its common.h
* for the full explanation of the vocabulary, of the HVBM gamma5 and of the
* sign of contracted epsilon tensors).  The amplitudes come from
* amplitudes.inc, written by step 2 from the FeynArts records:
*   Chain(Ubar(p),G(mu),Sl(q),...,U(k))  the open fermion line ubar(p) gamma^mu q-slash ... u(k)
*   Chain(Ubar(p),...,V(r))              the line ubar(p) ... v(r) ending on an outgoing antiquark
*   Den(q) = 1/q^2        Pol(p,mu)      polarization vector of the boson with momentum p
*   T(a,i,j) = (T^a)_ij   dF(i,j)        fundamental delta
* The preprocessor variable SPIN (U or P) selects the spin average or the
* helicity difference of the incoming particles.

#include declarations.inc
Symbols s,t,u,a,b,s4,spe3,eps,Nc,CF,TR,dm4;
#if "`SPIN'" == "U"
  #define MINKOWSKISIGN "1"
#else
* every polarized square contains exactly two epsilon tensors (two gamma5's,
* or one gamma5 and one circular density, or two circular densities), which
* FORM contracts with the Euclidean sign: hence the factor -1
  #define MINKOWSKISIGN "-1"
#endif
Indices mu,nu,rho,sigma;
* Four-dimensional indices, used to define gamma5 in the HVBM scheme
Indices h1=4,h2=4,h3=4,h4=4,e1=4,e2=4,e3=4,e4=4;
Indices ii,jj,aa;
CFunctions Line,LineBar,ChainBar,TBar,PolBar;

* conjugate(E): complex conjugate of the expression E (see the lepton-quark process)
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
* fermion line into one gamma string g_(line,...), inserting the spin sum
* p-slash of every outgoing fermion (also for an outgoing antiquark: the sum
* over v(r) vbar(r) is r-slash) and the density matrix of an incoming one:
* k-slash/2 (spin average) or k-slash gamma5/2 (helicity difference), with
* gamma5 = -i/24 e_(h1,h2,h3,h4) gamma^h1 gamma^h2 gamma^h3 gamma^h4.
* Line 1 is the lepton (incoming l, outgoing lp); line 2 the quark line:
*   l q -> l q (Born):                 incoming k, outgoing p;
*   l q -> l g q and gamma q -> g q:  incoming k, outgoing r;
*   l g -> l q qbar and gamma g -> q qbar:  outgoing p and r, no density.
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
  #if "`SPIN'" == "U"
    id Line(lp,?x,l)*LineBar(l,?y,lp) = g_(1,lp,?x,l,?y)/2;
    id Line(p,?x,k)*LineBar(k,?y,p)   = g_(2,p,?x,k,?y)/2;
    id Line(r,?x,k)*LineBar(k,?y,r)   = g_(2,r,?x,k,?y)/2;
  #else
    id Line(lp,?x,l)*LineBar(l,?y,lp) = -i_/48*e_(h1,h2,h3,h4)*g_(1,lp,?x,l,h1,h2,h3,h4,?y);
    id Line(p,?x,k)*LineBar(k,?y,p)   = -i_/48*e_(e1,e2,e3,e4)*g_(2,p,?x,k,e1,e2,e3,e4,?y);
    id Line(r,?x,k)*LineBar(k,?y,r)   = -i_/48*e_(e1,e2,e3,e4)*g_(2,r,?x,k,e1,e2,e3,e4,?y);
  #endif
  id Line(p,?x,r)*LineBar(r,?y,p) = g_(2,p,?x,r,?y);
  .sort
#endprocedure

* kinematics of the real emission, D-dimensional products (first pass) and
* the four-dimensional ("barred") products after the contraction of two
* epsilon tensors (second pass; only lp and r have hatted parts)
#procedure realKinematics
  id l.l = 0; id k.k = 0; id lp.lp = 0; id p.p = 0; id r.r = 0;
  id l.k = s/2; id k.p = -t/2; id l.p = -u/2; id k.r = -a/2; id p.r = b/2;
  id l.lp = -(a+b+t)/2; id l.r = (s+t+u+a+b)/2; id k.lp = (s+t+a)/2;
  id lp.p = -(t+u+b)/2; id lp.r = (s+t+u)/2;
#endprocedure
#procedure realKinematicsBarred
  id r.r = -spe3; id lp.lp = -spe3; id lp.r = (s+t+u)/2 + spe3;
  id l.l = 0; id k.k = 0; id p.p = 0;
  id l.k = s/2; id k.p = -t/2; id l.p = -u/2; id k.r = -a/2; id p.r = b/2;
  id l.lp = -(a+b+t)/2; id l.r = (s+t+u+a+b)/2; id k.lp = (s+t+a)/2; id lp.p = -(t+u+b)/2;
#endprocedure
* kinematics of the photon processes (all momenta four-dimensional)
#procedure photonKinematics
  id q.q = 0; id k.k = 0; id p.p = 0; id r.r = 0;
  id q.k = sh/2; id k.p = -th/2; id q.p = -uh/2; id p.r = sh/2; id k.r = -uh/2; id q.r = -th/2;
#endprocedure
