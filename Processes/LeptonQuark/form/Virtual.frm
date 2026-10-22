* Virtual.frm: the one-loop vertex graph projected onto the Born current.
* The quark chain of the vertex graph is
*   ubar(p) gamma^a (p-slash - loop-slash) gamma^mu (k-slash - loop-slash) gamma_a u(k)
* (the gluon propagator's metric contracts the two outer indices).  After loop
* integration it equals F times the Born current gamma^mu, and the scalar
* integrand of F is
*   Tr[p-slash W k-slash gamma_mu] / Tr[p-slash gamma^mu k-slash gamma_mu] .
* Written in the inverse propagators d0 = loop^2, d1 = (loop-k)^2,
* d2 = (loop-p)^2, every monomial lowers the power of a propagator, giving the
* integral family VirtualIntegral(n0,n1,n2) = int 1/(d0^n0 d1^n1 d2^n2).
* The result does not depend on the spin choice.  The lepton chain, photon
* propagator and couplings are common to the vertex and Born graphs and are
* dropped; so is the colour factor CF, restored in the notebook.
#-
Off statistics;
#include common.h
CFunctions VirtualIntegral;
Symbols n1,n2,n3;
Indices i1,i2,i3;
Vectors q1,q2;
#include amplitudes.inc
.sort

L Numerator = Vertex1;
L BornCurrent = Born1;
.sort
Drop Born1,Real1,Real2,Vertex1,Photon1,Photon2;

* keep only the quark chain and the gluon-propagator metric; the lepton chain,
* photon propagator and couplings are common to vertex and Born graphs
Skip; NSkip Numerator,BornCurrent;
id Chain(Ubar(lp),?x,U(l)) = 1;
id Den(q?) = 1;
id dF(ii?,jj?) = 1;
id T(aa?,ii?,jj?) = 1;
id ee = 1; id eq = 1; id gs = 1;
id i_ = 1;
id Chain(Ubar(p),?x,U(k)) = Line(p,?x,k);
repeat;
  id Line(?x,G(mu?),?y) = Line(?x,mu,?y);
  id Line(?x,Sl(q?),?y) = Line(?x,q,?y);
endrepeat;
* close the quark chain with the Born current gamma_mu, where mu is the photon
* index: the middle one of the three indices of the vertex chain, the only one
* of the Born chain (the metric of the lepton line has already been absorbed
* into it, the metric of the gluon propagator contracts the two outer ones)
id Line(p,i1?,q1?,i2?,q2?,i3?,k) = g_(1,p,i1,q1,i2,q2,i3,k,i2);
id Line(p,i1?,k) = g_(1,p,i1,k,i1);
tracen,1;
contract;
id k.k = 0; id p.p = 0; id k.p = -t/2;
id loop.loop = d0;
id loop.k = (d0-d1)/2;
id loop.p = (d0-d2)/2;
.sort

* map every monomial in d0, d1, d2 onto the integral family: starting from
* VirtualIntegral(1,1,1) = int 1/(d0 d1 d2), each power of d_i cancels one
* propagator.  The division by the Born current trace is done in the notebook.
Skip; NSkip Numerator;
Multiply VirtualIntegral(1,1,1);
repeat;
  id VirtualIntegral(n1?,n2?,n3?)*d0 = VirtualIntegral(n1-1,n2,n3);
  id VirtualIntegral(n1?,n2?,n3?)*d1 = VirtualIntegral(n1,n2-1,n3);
  id VirtualIntegral(n1?,n2?,n3?)*d2 = VirtualIntegral(n1,n2,n3-1);
endrepeat;
.sort

Print +s BornCurrent;
Print +s Numerator;
Format Mathematica;
#write <Virtual.m> "HoldComplete[{%E,%E}]", Numerator,BornCurrent
.end
