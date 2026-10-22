(* NLO TMD beam-function matching coefficients, soft function and TMD PDF
   coefficients of M.-x. Luo, T.-Z. Yang, H. X. Zhu and Y. J. Zhu,
   "Unpolarized quark and gluon TMD PDFs and FFs at N3LO", JHEP 06 (2021) 115
   [arXiv:2012.03256], exactly as given in the ancillary files of the arXiv
   submission (BeamFunction.m, softfunction.m, TMDPDF.m), truncated at order
   as^1.  The paper's conventions (its ancillary_readme.pdf):

     as   = alpha_s(mu)/(4 Pi)
     Lp   = ln(b_T^2 mu^2 / b0^2),  b0 = 2 e^(-EulerGamma)
     LQ   = 2 ln(Q/nu),  Q = x P^+ the momentum of the parton that enters the operator
     Lv   = ln(nu^2/mu^2)                      (written L_nu in the paper)
     Lh   = Lv + LQ = ln(Q^2/mu^2)
     plusD[0, 1 - x] = [1/(1-x)]_+,  delta[1 - x] = delta(1-x),  zeta[2] = Zeta[2]

   Beam["ji"] is the matching coefficient I_ji(x, b_T, mu, nu) of the beam
   function with the exponential rapidity regulator (nu = 1/tau), soft["i"] the
   soft function S_i(b_T, mu, nu), and TMDpdf["ji"] the regulator-independent
   TMD PDF coefficient f_ji = I_ji Sqrt[S_i].

   This function returns the O(as) coefficients (the as^0 terms are
   delta[1 - x] for the diagonal channels, 1 for the soft function). *)

Luo2020NLO[x_, Lp_, LQ_, Lv_, Lh_, CF_, CA_, TF_] := <|
  "Beam" -> <|
    "qq" -> 2*CF + 2*CF*Lp - 2*CF*x + 2*CF*Lp*x - 2*CF*Lp*LQ*delta[1 - x] - 4*CF*Lp*plusD[0, 1 - x],
    "qg" -> -2*Lp*TF + 4*TF*x + 4*Lp*TF*x - 4*TF*x^2 - 4*Lp*TF*x^2,
    "gq" -> 4*CF*Lp - (4*CF*Lp)/x + 2*CF*x - 2*CF*Lp*x,
    "gg" -> 8*CA*Lp - (4*CA*Lp)/x - 4*CA*Lp*x + 4*CA*Lp*x^2 - 2*CA*Lp*LQ*delta[1 - x] - 4*CA*Lp*plusD[0, 1 - x]|>,
  "Soft" -> <|
    "q" -> CF*(-2*Lp^2 - 4*Lp*Lv - 2*zeta[2]),
    "g" -> CA*(-2*Lp^2 - 4*Lp*Lv - 2*zeta[2])|>,
  "TMDPDF" -> <|
    "qq" -> 2*CF + 2*CF*Lp - 2*CF*x + 2*CF*Lp*x - 2*CF*Lh*Lp*delta[1 - x] - CF*Lp^2*delta[1 - x] - 4*CF*Lp*plusD[0, 1 - x] - CF*delta[1 - x]*zeta[2],
    "qg" -> -2*Lp*TF + 4*TF*x + 4*Lp*TF*x - 4*TF*x^2 - 4*Lp*TF*x^2,
    "gq" -> 4*CF*Lp - (4*CF*Lp)/x + 2*CF*x - 2*CF*Lp*x,
    "gg" -> 8*CA*Lp - (4*CA*Lp)/x - 4*CA*Lp*x + 4*CA*Lp*x^2 - 2*CA*Lh*Lp*delta[1 - x] - CA*Lp^2*delta[1 - x] - 4*CA*Lp*plusD[0, 1 - x] - CA*delta[1 - x]*zeta[2]|>|>
