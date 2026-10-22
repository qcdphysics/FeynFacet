(* The polarized NLO off-diagonal coefficients of Hinderer, Schlegel and
   Vogelsang, Phys. Rev. D 96 (2017) 014002 [arXiv:1703.10872], Eqs. (13), (14)
   and Appendix A, and the polarized photon Born cross sections of Eq. (15),
   transcribed from the LaTeX source of the paper exactly as printed.  The
   paper writes the double-longitudinal spin-dependent cross sections as

     Delta sigma_NLO^{q->g} = alpha_em^2 e_q^2 CF/(s v w) [DC1 ln(1-v(1-w)) + DC2 ln((1-v)/((1-vw)(1-v(1-w)))) + DC3 ln(v(1-w)s/mu^2) + DC4],
     Delta sigma_NLO^{g->q} = alpha_em^2 e_q^2 TR/(s v w) [DC1 ln((1-v)w/(1-vw)) + DC2 ln(v(1-w)s/mu^2) + DC3],

   as the coefficient of alpha_s/Pi, with a single scale mu, and the
   polarized Weizsaecker-Williams (photon) Born terms

     Delta sigma_LO^{gamma q -> g} =  2 Pi CF alpha_em e_q^2/(s(1-v)) v w (2-vw)/(1-vw),
     Delta sigma_LO^{gamma g -> q} = -2 Pi TR alpha_em e_q^2/(s(1-v)) (v^2 w^2 + (1-vw)^2)/(v w (1-vw)).

   The functions return the complete NLO coefficient (the quantity multiplying
   alpha_s/Pi) and the polarized Born photon cross sections. *)

HSV2017OffDiagonal["q -> g", s_, v_, w_, mu_, alphaEM_, eq_, CF_] := Module[{C1, C2, C3, C4},
  C1 = 2 v w (1 - v (1 - 2 w)) (1 + v^2 (1 - w)^2)/((1 - v) (1 - v (1 - w))^2);
  C2 = -2 v w (1 + v (1 - 2 w))/(1 - v);
  C3 = v w/((1 - v) (1 - v w)^2 (1 - v (1 - w))^2) (1 - v (1 - 2 w) + v^2 (1 - 2 w - 3 w^2) - v^3 (1 - 7 w^2 + 4 w^3)
    + 2 v^4 w (2 - 7 w + 7 w^2 - 2 w^3) - 2 v^5 w^2 (1 - w)^2 (1 - 2 w));
  C4 = v w/((1 - v) (1 - v w)^2 (1 - v (1 - w))^2) (2 - 4 v + v^2 (4 - 13 w + 12 w^2) - v^3 (2 - 14 w + 13 w^2)
    - v^4 w (3 - 4 w + 3 w^2 - 2 w^3) + v^5 w^2 (1 - w)^2);
  alphaEM^2 eq^2 CF/(s v w) (C1 Log[1 - v (1 - w)] + C2 Log[(1 - v)/((1 - v w) (1 - v (1 - w)))] + C3 Log[v (1 - w) s/mu^2] + C4)];

HSV2017OffDiagonal["g -> q", s_, v_, w_, mu_, alphaEM_, eq_, TR_] := Module[{C1, C2, C3},
  C1 = -2 (1 + v (1 - 2 w))/(1 - v);
  C2 = -1/((1 - v) (1 - v w)^2) (2 - 2 w + 2 v (1 - 4 w + 2 w^2) - v^2 w (4 - 11 w + 2 w^2) + 3 v^3 w^2 (1 - 2 w));
  C3 = 1/((1 - v) (1 - v w)^2) (2 (1 - w) + 2 v (2 - 5 w + 2 w^2) - v^2 w (7 - 12 w + 2 w^2) + v^3 w^2 (7 - 8 w));
  alphaEM^2 eq^2 TR/(s v w) (C1 Log[(1 - v) w/(1 - v w)] + C2 Log[v (1 - w) s/mu^2] + C3)];

HSV2017PhotonBorn["gamma q -> g", s_, v_, w_, alphaEM_, eq_, CF_] := 2 Pi CF alphaEM eq^2/(s (1 - v)) v w (2 - v w)/(1 - v w);
HSV2017PhotonBorn["gamma g -> q", s_, v_, w_, alphaEM_, eq_, TR_] := -2 Pi TR alphaEM eq^2/(s (1 - v)) (v^2 w^2 + (1 - v w)^2)/(v w (1 - v w));
