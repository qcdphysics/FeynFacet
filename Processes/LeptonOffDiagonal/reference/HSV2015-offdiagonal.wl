(* The NLO off-diagonal coefficients of Hinderer, Schlegel and Vogelsang,
   Phys. Rev. D 92 (2015) 014001 [arXiv:1505.06415v3], Eqs. (27), (28) and
   Appendix A, and the photon-quark Born cross sections of Eq. (29), transcribed
   from the LaTeX source of the paper exactly as printed.  The paper writes

     sigma_NLO^{q->g} = alpha_em^2 e_q^2 CF/(s v w) [C1 ln(1-v(1-w)) + C2 ln((1-v)/((1-vw)(1-v(1-w)))) + C3 ln(v(1-w)s/mu^2) + C4],
     sigma_NLO^{g->q} = alpha_em^2 e_q^2 TR/(s v w) [C1 ln((1-v)w/(1-vw)) + C2 ln(v(1-w)s/mu^2) + C3],

   as the coefficient of alpha_s/Pi, with a single scale mu for renormalization
   and all factorizations, and the Weizsaecker-Williams (photon) Born terms

     sigma_LO^{gamma q -> g} = 2 Pi CF alpha_em e_q^2/(s(1-v)) (1+(1-vw)^2)/(1-vw),
     sigma_LO^{gamma g -> q} = 2 Pi TR alpha_em e_q^2/(s(1-v)) (v^2 w^2 + (1-vw)^2)/(v w (1-vw)).

   The functions return the complete NLO coefficient (the quantity multiplying
   alpha_s/Pi) and the Born photon cross sections. *)

HSV2015OffDiagonal["q -> g", s_, v_, w_, mu_, alphaEM_, eq_, CF_] := Module[{C1, C2, C3, C4},
  C1 = 2 v w (1 + v^2 (1 - w)^2)/((1 - v)^2 (1 - v (1 - w))^2) (1 - 2 v (1 - w) + v^2 (1 - 2 w (1 - w)));
  C2 = v w (6 - 4 v w + 2 v^2 (1 - 2 w (1 - w)))/(1 - v)^2;
  C3 = v w/((1 - v)^2 (1 - v w)^2 (1 - v (1 - w))^2) (3 - 2 v (3 + w) + v^2 (6 + 4 w - w^2) - 2 v^3 (3 - 3 w + 5 w^2 - 2 w^3)
    + v^4 (3 - 4 w + 5 w^2 - 2 w^3) - 2 v^5 w (2 - 6 w + 9 w^2 - 7 w^3 + 2 w^4) + 2 v^6 (1 - w)^2 w^2 (1 - 2 w + 2 w^2));
  C4 = v w/((1 - v)^2 (1 - v w)^2 (1 - v (1 - w))^2) (2 - 2 v (5 - 3 w) + v^2 (16 - 3 w - 11 w^2) - v^3 (10 + 15 w - 27 w^2 + 2 w^3)
    + v^4 (2 + 17 w - 23 w^2 + 7 w^3 - 3 w^4) - v^5 w (5 - 5 w - w^2 + 3 w^3 - 2 w^4) + 2 v^6 (1 - w)^2 w^2 (1 - w + w^2));
  alphaEM^2 eq^2 CF/(s v w) (C1 Log[1 - v (1 - w)] + C2 Log[(1 - v)/((1 - v w) (1 - v (1 - w)))] + C3 Log[v (1 - w) s/mu^2] + C4)];

HSV2015OffDiagonal["g -> q", s_, v_, w_, mu_, alphaEM_, eq_, TR_] := Module[{C1, C2, C3},
  C1 = 2 (1 + v (4 v w^2 - 2 w (1 + v) + v))/(1 - v)^2;
  C2 = 1/((1 - v)^2 (1 - v w)^2) (2 (1 - w + w^2) - 2 v w (3 - 2 w + 2 w^2) + v^2 (2 - 4 w + 11 w^2 - 2 w^3 + 2 w^4)
    - 4 v^3 w (1 - 2 w + 3 w^2) + 3 v^4 w^2 (1 - 2 w + 2 w^2));
  C3 = 1/((1 - v)^2 (1 - v w)^2) (1 + 4 w - 6 w^2 - 2 v (1 + 3 w + w^2 - 6 w^3) + v^2 (1 + 9 w + 4 w^2 - 8 w^3 - 6 w^4)
    - v^3 w (3 + 9 w - 4 w^2 - 6 w^3) + v^4 w^2 (1 + 4 w - 4 w^2));
  alphaEM^2 eq^2 TR/(s v w) (C1 Log[(1 - v) w/(1 - v w)] + C2 Log[v (1 - w) s/mu^2] + C3)];

HSV2015PhotonBorn["gamma q -> g", s_, v_, w_, alphaEM_, eq_, CF_] := 2 Pi CF alphaEM eq^2/(s (1 - v)) (1 + (1 - v w)^2)/(1 - v w);
HSV2015PhotonBorn["gamma g -> q", s_, v_, w_, alphaEM_, eq_, TR_] := 2 Pi TR alphaEM eq^2/(s (1 - v)) (v^2 w^2 + (1 - v w)^2)/(v w (1 - v w));
