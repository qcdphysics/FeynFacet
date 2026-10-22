(* Unpolarized NLO hard coefficient for l q -> q + X from
   P. Hinderer, M. Schlegel, W. Vogelsang, "Single-inclusive production of
   hadrons and jets in lepton-nucleon scattering at NLO", Phys. Rev. D 92 (2015)
   014001, arXiv:1505.06415, Eq. (26) and Appendix A, transcribed
   as printed.  Used for the final comparison only.

   With vb = 1 - v, rho = 1 - v w, d = 1 - v + v w, D0 = [1/(1-w)]_+ and
   D1 = [log(1-w)/(1-w)]_+ , the paper writes the NLO term as

     CF alphaEM^2 eq^2 alphaS/(Pi s v w) *
       { A0 delta(1-w) + A1 D1 + [B1 log(vb/(v d)) + B2 log d + B3 log(s v^2/muR^2)] D0
         + C1 log(v(1-w)) + C2 log(vb w/rho) + C3 log(vb/(rho d)) + C4 log(s/muR^2) + C5 }.

   The function returns the prefactor and the four coefficients {R, B, C, A}
   of {1, D0, D1, delta(1-w)} in this form; note that B and C depend on w. *)

HSV2015Unpolarized[s_, v_, w_, muR_, alphaEM_, alphaS_, eq_, CF_] := Module[
  {vb = 1 - v, rho = 1 - v w, d = 1 - v + v w, pBorn, A0, A1, B1, B2, B3, C1, C2, C3, C4, C5},
  pBorn = (1 + v^2)/vb^2;
  A0 = pBorn ((3 + 2 Log[v]) Log[s vb/muR^2] + Log[v]^2 - 8);
  A1 = 8 w pBorn;
  B1 = 4 w/vb^2 (1 - v (1 - w) + v^2 (1 - w (1 - w)));
  B2 = 2 w/(vb^2 d) (1 - 2 v (1 - w) + v^2 (1 - 2 w + 2 w^2)) (2 - 2 v (1 - w) + v^2 (1 - w)^2);
  B3 = 4 w pBorn;
  C1 = (2 - w - 2 v (1 + 4 w) + v^2 (2 + 9 w - 10 w^2 + w^3) - 2 v^3 (1 - w + w^2 - 4 w^3)
        + v^4 w (2 - 2 w - 7 w^2 + 8 w^3) - 2 v^5 w^2 (1 - 3 w + 4 w^2 - 2 w^3))/(vb^2 rho d);
  C2 = 2 (1 + v^2 (1 + 2 w^2))/vb^2;
  C3 = -2 v w/vb^2 (3 - 2 v (1 - w) + v^2 (1 - 2 w + 2 w^2));
  C4 = (2 - w - 2 v (1 + 2 w) + v^2 (2 + 5 w - 6 w^2 + w^3) - 2 v^3 (1 - w + w^2 - 2 w^3)
        + v^4 w (2 - 2 w - 3 w^2 + 4 w^3) - 2 v^5 w^2 (1 - 3 w + 4 w^2 - 2 w^3))/(vb^2 rho d);
  C5 = (2 - w - 2 v (2 - w) + v^2 (4 - w - 2 w^2 + w^3) - 2 v^3 + v^4 w (2 - w^2)
        - 2 v^5 w^2 (1 - 2 w + 2 w^2 - w^3))/(vb^2 rho d);
  <|"Prefactor" -> CF alphaEM^2 alphaS eq^2/(Pi s v w),
    "Coefficients" -> {C1 Log[v (1 - w)] + C2 Log[vb w/rho] + C3 Log[vb/(rho d)] + C4 Log[s/muR^2] + C5,
      B1 Log[vb/(v d)] + B2 Log[d] + B3 Log[s v^2/muR^2], A1, A0}|>]
