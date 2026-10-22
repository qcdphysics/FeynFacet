* checks.frm: the two forms of the HVBM gamma5 used in common.h (procedure sew).
*   (1) gamma5 = -i/24 e_(a,b,c,d) gamma^a gamma^b gamma^c gamma^d (four-
*       dimensional indices a..d) gives Tr[gamma^w gamma^x gamma^y gamma^z gamma5]
*       = -4 i e_(w,x,y,z), the convention of FeynCalc with
*       $LeviCivitaSign = -1 used by the other routes;
*   (2) for a four-dimensional vector q, q-slash gamma5 = -i/6 e_(q,a,b,c)
*       gamma^a gamma^b gamma^c; checked inside a trace that also contains
*       D-dimensional vectors k1, k2, loop and D-dimensional indices mu, nu,
*       with q = pa.
* Gamma5Checks is the sum of the two differences; step 2 requires that FORM
* prints Gamma5Checks = 0.
#-
Off statistics;
Symbols D;
Dimension D;
Vectors pa, k1, k2, loop;
Indices mu, nu, w, x, y, z, h1=4, h2=4, h3=4, h4=4, e1=4, e2=4, e3=4;
L Gamma5Checks =
* (1) the trace of four gamma matrices and gamma5, minus (-4 i e_(w,x,y,z))
    (-i_/24)*e_(h1,h2,h3,h4)*g_(1,w,x,y,z,h1,h2,h3,h4) + 4*i_*e_(w,x,y,z)
* (2) a trace with pa-slash gamma5 in the long form, minus the same in the short form
  + (-i_/24)*e_(h1,h2,h3,h4)*g_(2,k1,mu,k2,loop,nu,mu,nu,pa,h1,h2,h3,h4)
  - (-i_/6)*e_(pa,e1,e2,e3)*g_(2,k1,mu,k2,loop,nu,mu,nu,e1,e2,e3);
tracen,1;
tracen,2;
contract;
Print +s;
.end
