* common.h: declarations and procedures shared by the FORM programs of this
* process.  Read Processes/LeptonQuark/form/common.h first: the vocabulary of
* amplitudes.inc, the conjugation, the HVBM gamma5 and the sign of contracted
* epsilon tensors are explained there.  The preprocessor variable SPIN selects
* the spin average (U) or the transverse spin transfer with the normal axes
* (NN) or the in-plane axes (XX).

#include declarations.inc
Symbols s,t,u,a,b,s4,spe3,eps,Nc,CF,dm4,nr;
Symbols d0,d1,d2;
* the spin vectors: nn is the normal to the scattering plane; aXs = s aX and
* bXs = (t+u)^2 bX are the in-plane axes of the FeynCalc route scaled to
* polynomial coefficients, aXs = t l + u k + s p, bXs = t(t+u) l - u(t+u) k - s(u-t) p
* (the notebook divides by the scales and by the lengths of the axes)
Vectors nn,aXs,bXs;
#if "`SPIN'" == "U"
  #define MINKOWSKISIGN "1"
#else
* two gamma5's: two epsilon tensors, contracted by FORM with the Euclidean sign
  #define MINKOWSKISIGN "-1"
#endif
Indices mu,nu,rho,sigma;
* Four-dimensional indices, used to define gamma5 in the HVBM scheme
Indices h1=4,h2=4,h3=4,h4=4,e1=4,e2=4,e3=4,e4=4;
Indices ii,jj,aa;
Vectors x;
CFunctions Line,LineBar,ChainBar,TBar,PolBar;

* conjugate(E): complex conjugate of the expression E
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
* fermion line into one gamma string g_(line,...).  Line 1 is the lepton
* (unpolarized: density l-slash/2, spin sum lp-slash).  Line 2 is the quark
* line with the outgoing quark p (r in the photon process) and the incoming
* quark k: for SPIN=U the density k-slash/2 and the spin sum p-slash; for the
* transverse spin transfer the density k-slash gamma5 sIn-slash/2 and the
* projection p-slash gamma5 sOut-slash, with sIn = sOut = nn (SPIN=NN) or
* sIn = aXs, sOut = bXs (SPIN=XX).
*
* gamma5 in the HVBM scheme is written as
*   gamma5 = -i/24 e_(h1,h2,h3,h4) gamma^h1 gamma^h2 gamma^h3 gamma^h4
* with four-dimensional indices (Born.frm of the lepton-quark process checks
* the normalization); two gamma5's give (-i/24)^2 = -1/576, and with the
* density's 1/2: -1/1152.  ProjectEpsilon (declarations.inc) and MINKOWSKISIGN
* are explained in the lepton-quark process.
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
  id Line(lp,?x,l)*LineBar(l,?y,lp) = g_(1,lp,?x,l,?y)/2;
  #if "`SPIN'" == "U"
    id Line(p,?x,k)*LineBar(k,?y,p) = g_(2,p,?x,k,?y)/2;
    id Line(r,?x,k)*LineBar(k,?y,r) = g_(2,r,?x,k,?y)/2;
  #endif
  #if "`SPIN'" == "NN"
    id Line(p,?x,k)*LineBar(k,?y,p) = -1/1152*e_(h1,h2,h3,h4)*e_(e1,e2,e3,e4)*g_(2,p,h1,h2,h3,h4,nn,?x,k,e1,e2,e3,e4,nn,?y);
    id Line(r,?x,k)*LineBar(k,?y,r) = -1/1152*e_(h1,h2,h3,h4)*e_(e1,e2,e3,e4)*g_(2,r,h1,h2,h3,h4,nn,?x,k,e1,e2,e3,e4,nn,?y);
  #endif
  #if "`SPIN'" == "XX"
    id Line(p,?x,k)*LineBar(k,?y,p) = -1/1152*e_(h1,h2,h3,h4)*e_(e1,e2,e3,e4)*g_(2,p,h1,h2,h3,h4,bXs,?x,k,e1,e2,e3,e4,aXs,?y);
    id Line(r,?x,k)*LineBar(k,?y,r) = -1/1152*e_(h1,h2,h3,h4)*e_(e1,e2,e3,e4)*g_(2,r,h1,h2,h3,h4,bXs,?x,k,e1,e2,e3,e4,aXs,?y);
  #endif
  .sort
#endprocedure

* the in-plane axes: their dot products with any vector follow from the
* definitions (call after tracen and after contract, when the axes appear only
* in dot products); the arguments are the two incoming momenta, the observed
* one and the invariants
#procedure axes(LIN,KIN,POUT,S,T,U)
  id aXs.aXs = `T'*`LIN'.aXs + `U'*`KIN'.aXs + `S'*`POUT'.aXs;
  id bXs.bXs = `T'*(`T'+`U')*`LIN'.bXs - `U'*(`T'+`U')*`KIN'.bXs - `S'*(`U'-`T')*`POUT'.bXs;
  id aXs.bXs = `T'*`LIN'.bXs + `U'*`KIN'.bXs + `S'*`POUT'.bXs;
  id aXs.x? = `T'*`LIN'.x + `U'*`KIN'.x + `S'*`POUT'.x;
  id bXs.x? = `T'*(`T'+`U')*`LIN'.x - `U'*(`T'+`U')*`KIN'.x - `S'*(`U'-`T')*`POUT'.x;
#endprocedure
* the normal: orthogonal to the given three momenta, nn.nn = -1
#procedure normal(LIN,KIN,POUT)
  id nn.`LIN' = 0; id nn.`KIN' = 0; id nn.`POUT' = 0; id nn.nn = -1;
#endprocedure
