* checks.frm: the forms of the HVBM gamma5, of the transverse density of the
* incoming parton and of the projection of the observed parton used in
* common.h (procedure sew).
*   (1) gamma5 = -i/24 e_(a,b,c,d) gamma^a gamma^b gamma^c gamma^d (four-
*       dimensional indices a..d) gives Tr[gamma^w gamma^x gamma^y gamma^z gamma5]
*       = -4 i e_(w,x,y,z), the convention of FeynCalc with
*       $LeviCivitaSign = -1 used by the FeynCalc route;
*   (2) for a four-dimensional vector q, q-slash gamma5 = -i/6 e_(q,a,b,c)
*       gamma^a gamma^b gamma^c; checked inside a trace that also contains
*       D-dimensional vectors k1, k2, loop and D-dimensional indices mu, nu,
*       with q = pa;
*   (3) the density pa-slash gamma5 sa-slash and the projection p-slash gamma5
*       sc-slash in their short forms (-i/6 e_(q,a,b,c) gamma^a gamma^b gamma^c
*       s-slash, as inserted by sew) against the long forms with the
*       four-matrix gamma5, inside one trace with D-dimensional objects between
*       them (the exchange graphs: a and c on one line), and
*   (4) the same with the two projections in two separate traces (the
*       annihilation-exchange interference: a and c on different lines).
* Gamma5Checks is the sum of the four differences; step 2 requires that FORM
* prints Gamma5Checks = 0.
#-
Off statistics;
Symbols D;
Dimension D;
Vectors pa, p, sa, sc, k1, k2, loop;
Indices mu, nu, w, x, y, z, h1=4, h2=4, h3=4, h4=4, e1=4, e2=4, e3=4, e4=4;
L Gamma5Checks =
* (1) the trace of four gamma matrices and gamma5, minus (-4 i e_(w,x,y,z))
    (-i_/24)*e_(h1,h2,h3,h4)*g_(1,w,x,y,z,h1,h2,h3,h4) + 4*i_*e_(w,x,y,z)
* (2) a trace with pa-slash gamma5 in the long form, minus the same in the short form
  + (-i_/24)*e_(h1,h2,h3,h4)*g_(2,k1,mu,k2,loop,nu,mu,nu,pa,h1,h2,h3,h4)
  - (-i_/6)*e_(pa,e1,e2,e3)*g_(2,k1,mu,k2,loop,nu,mu,nu,e1,e2,e3)
* (3) one trace with the density of a and the projection of c, both in the long form, minus the same with both in the short form
  + (-i_/24)^2*e_(h1,h2,h3,h4)*e_(e1,e2,e3,e4)*g_(3,k1,mu,pa,h1,h2,h3,h4,sa,nu,loop,mu,p,e1,e2,e3,e4,sc,nu)
  - (-i_/6)^2*e_(pa,h1,h2,h3)*e_(p,e1,e2,e3)*g_(3,k1,mu,h1,h2,h3,sa,nu,loop,mu,e1,e2,e3,sc,nu)
* (4) the density and the projection in two traces (lines 4 and 5), long minus short
  + (-i_/24)^2*e_(h1,h2,h3,h4)*e_(e1,e2,e3,e4)*g_(4,k1,mu,pa,h1,h2,h3,h4,sa,nu)*g_(5,loop,mu,p,e1,e2,e3,e4,sc,nu,k2,w,w)
  - (-i_/6)^2*e_(pa,h1,h2,h3)*e_(p,e1,e2,e3)*g_(4,k1,mu,h1,h2,h3,sa,nu)*g_(5,loop,mu,e1,e2,e3,sc,nu,k2,w,w);
tracen,1;
tracen,2;
tracen,3;
tracen,4;
tracen,5;
contract;
Print +s;
.end
