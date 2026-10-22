(* Double-longitudinal NLO hard coefficient for l q -> q + X from
   P. Hinderer, M. Schlegel, W. Vogelsang, "Double-longitudinal spin asymmetry
   in single-inclusive lepton scattering at NLO", Phys. Rev. D 96 (2017) 014002,
   arXiv:1703.10872, Eqs. (11)-(14) and Appendix A, transcribed as printed.
   Used for the final comparison only.

   The paper writes the NLO term as
     CF alphaEM^2 eq^2 alphaS/(Pi s v w) { R + B D0 + C D1 + A delta(1-w) }
   with all factorization scales equal to the renormalization scale muR.  The
   function returns that prefactor and {R, B, C, A}.  Logarithms: LV = log v,
   LVbar = log(1-v), LW = log w, LWbar = log(1-w), LD = log(1-v+vw),
   LRho = log(1-vw), LS = log(s/muR^2). *)

HSV2017Polarized[s_, v_, w_, muR_, alphaEM_, alphaS_, eq_, CF_] := Module[
  {LV = Log[v], LVbar = Log[1 - v], LW = Log[w], LWbar = Log[1 - w], LD = Log[1 - v + v w],
   LRho = Log[1 - v w], LS = Log[s/muR^2]},
  <|"Prefactor" -> CF alphaEM^2 alphaS eq^2/(Pi s v w),
    "Coefficients" -> {
 (-1 + 2 LVbar + 6 LWbar + 2 LRho + 2 LS + 2 LV - 2 LW + 2 v - v^2 - 2 LVbar v^2 - 6 LWbar v^2
   - 2 LRho v^2 - 2 LS v^2 - 2 LV v^2 + 2 LW v^2 + w + LWbar w + LS w + LV w - 3 v w + 2 LD v w
   + 8 LVbar v w + 15 LWbar v w + 7 LS v w + 7 LV v w - 2 LW v w + 3 v^2 w - 2 LD v^2 w
   - 6 LVbar v^2 w - 6 LWbar v^2 w + 6 LRho v^2 w - 6 LS v^2 w - 6 LV v^2 w - 2 LW v^2 w - v^3 w
   + 2 LVbar v^3 w + 2 LWbar v^3 w - 2 LRho v^3 w + 2 LS v^3 w + 2 LV v^3 w + v w^2 + LWbar v w^2
   + LS v w^2 + LV v w^2 - 2 v^2 w^2 + 10 LVbar v^2 w^2 + 13 LWbar v^2 w^2 - 6 LRho v^2 w^2
   + 9 LS v^2 w^2 + 9 LV v^2 w^2 + v^3 w^2 - 6 LVbar v^3 w^2 - 6 LWbar v^3 w^2 + 6 LRho v^3 w^2
   - 6 LS v^3 w^2 - 6 LV v^3 w^2 + 4 LVbar v^3 w^3 + 4 LWbar v^3 w^3 - 4 LRho v^3 w^3
   + 4 LS v^3 w^3 + 4 LV v^3 w^3)/((-1 + v) (1 - v + v w)),
 (-4 (LVbar + LS + LV) (1 + v))/(-1 + v),
 (-8 (1 + v))/(-1 + v),
 -(((-8 + 3 LVbar + 3 LS + 2 LVbar LV + 2 LS LV + LV^2) (1 + v))/(-1 + v))}|>]
