(* The NLO transverse-spin-transfer hard coefficient of l q -> q + X from the
   accompanying paper (Section "NLO transverse-spin-transfer hard part",
   Eqs. (transverse-finite-hard), (transverse-endpoints) and Appendix
   (transverse-regular)), transcribed as written there.  Used for the final
   comparison with the paper's formula only, not as an independent validation.
   The paper writes, with LI = log(muI^2/s), LO = log(muO^2/s),
   LL = log(muL^2/s) (initial, final and lepton-photon factorization scales),
   d = 1 - v + v w, D0 = [1/(1-w)]_+ and D1 = [log(1-w)/(1-w)]_+,

     sigma_T = 4 alphaEM^2 eq^2/(s (1-v)^2) delta(1-w)
             + alphaEM^2 eq^2 CF alphaS/(Pi s v w) [ RT + BT D0 + CT D1 + AT delta(1-w) ],

   for the full outgoing-spin difference T projected as (T_XX + T_NN)/2.  The
   function returns the prefactor and the coefficients {RT, BT, CT, AT}. *)

DraftTransverse[s_, v_, w_, muI_, muO_, muL_, alphaEM_, alphaS_, eq_, CF_] := Module[
  {LI = Log[muI^2/s], LO = Log[muO^2/s], LL = Log[muL^2/s], d = 1 - v + v w, pL, pa, pb, pc, pd, RT, BT, CT, AT},
  CT = 16 v/(1 - v)^2;
  BT = 4 v/(1 - v)^2 (2 Log[v] + 2 Log[1 - v] - LI - LO);
  AT = v/(1 - v)^2 (-16 + 6 Log[1 - v] - 3 LI - 3 LO + 4 Log[v] Log[1 - v] + 2 Log[v]^2 - 4 LO Log[v]);
  pL = 1 + v^2 - 2 v w - 2 v^2 w + 2 v^2 w^2;
  pa = -1 - v - 3 v^2 + v^3 - 2 v w + 2 v^2 w - 4 v^2 w^2;
  pb = -1 + 2 v - v^2 + v w + 2 v^2 w + v^3 w - 4 v^2 w^2;
  pc = -4 - w + 4 v w + v^2 w - 2 v w^2 - 2 v^2 w^2 + 4 v^2 w^3;
  pd = -7 - 2 v + v^2 - 4 w + 7 v w + 2 v^2 w - v^3 w + v w^2 - 2 v^2 w^2 + v^3 w^2 + 4 v^2 w^3;
  RT = 2 v/(1 - v)^2 (w (1 - v)^2/(1 - v w) + 2 LI + 2 (1 + w + v w^2) LO - w pL/(1 - v w) LL
    + w (3 + v^2) Log[w] + d (1 - v)^2/v Log[d] + pa/v Log[1 - v] - pb/v Log[1 - v w] + pc/(1 - v w) Log[v] + pd/(1 - v w) Log[1 - w]);
  <|"Born" -> 4 alphaEM^2 eq^2/(s (1 - v)^2), "Prefactor" -> alphaEM^2 eq^2 CF alphaS/(Pi s v w), "Coefficients" -> {RT, BT, CT, AT}|>];
