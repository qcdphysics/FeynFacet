(* ::Title:: *)
(*Lepton-quark scattering with transverse spin transfer, step 2 (FeynCalc route): spin and colour algebra*)


(* ::Text:: *)
(*Input: Results/Amplitudes.wl from step 1 (the diagrams of l q -> q + X).  Output: Results/FeynCalc/Algebra.wl with the squared matrix elements of the Born, real-emission and photon-quark processes for the spin average ("Unpolarized") and for the transverse spin transfer ("Transverse"), and the virtual form factor family, in the same format as Processes/LeptonQuark so that step 3 is the same.*)
(**)
(*The observable.  The incoming quark k carries a transverse spin, the transverse spin of the observed quark p is measured, the lepton is unpolarized and the scattered lepton unobserved.  T denotes the full difference sigma(k up -> p up) - sigma(k up -> p down) at fixed incoming spin; by parity this is half the double difference, so T is obtained from the spin-averaged square by replacing the incoming density k-slash/2 by k-slash gamma5 a-slash/2 and the outgoing spin sum p-slash by p-slash gamma5 b-slash, where a and b are the unit transverse spin vectors of the two quarks (u ubar = p-slash (1 + gamma5 s-slash)/2 for a massless fermion).  Two directions are needed: the normal n to the scattering plane, the same for both quarks (T_NN), and the in-plane transverse directions aX, bX of the lepton-quark centre-of-mass frame (T_XX):*)
(**)
(*   n: n.l = n.k = n.p = 0, n^2 = -1;   aX = (t/s) l + (u/s) k + p, orthogonal to l and k;   bX = t/(t+u) l - u/(t+u) k - s (u-t)/(t+u)^2 p, orthogonal to p and l + k,*)
(**)
(*with aX^2 = -t u/s and bX^2 = -s t u/(t+u)^2; T_XX is divided by the product of their lengths, -t u/(t+u), and the signs of the X axes are those of the paper (Section "Observable and angular normalization"): T_XX and T_NN coincide in four dimensions.  The Collins asymmetry needs the projection T_C = (T_XX + T_NN)/2, which is what is integrated ("Transverse"); the two components are kept separately in this file for inspection, since they differ by evanescent (D-4) terms.*)
(**)
(*HVBM scheme.  The physical momenta l, k, p and the spin vectors are four-dimensional, the unobserved momenta lp and r of the real emission have (D-4)-dimensional parts with spe3 = rhat.rhat = lphat.lphat = -lphat.rhat, and gamma5 anticommutes with the four-dimensional gamma matrices only ($LeviCivitaSign = -1).  The real-emission tensor depends on the component nr = n.r of the unobserved gluon momentum along the normal; terms odd in nr vanish after the phase-space integration (they are odd under reflection of the scattering plane), and nr^2 is a rational function of the invariants: writing r = rpar + rn + rhat (parts in the plane of l, k, p, along n, and in the extra dimensions) and using r^2 = 0 gives nr^2 = rpar^2 + spe3 with rpar^2 = d.G^-1.d, G the Gram matrix of l, k, p and d = (l.r, k.r, p.r).*)
(**)
(*Conventions as in Processes/LeptonQuark: s = 2 l.k, t = -2 k.p, u = -2 l.p, a = (k-r)^2, b = (p+r)^2, tau = a + b + t, s4 = s + t + u; for the photon process sh = 2 q.k, th = -2 k.p, uh = -2 q.p.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[];
amplitudes = LoadAmplitudes[]["FeynCalc"];
Keys[amplitudes]


(* ::Section:: *)
(*Spin vectors and kinematics*)


(* ::Text:: *)
(*The spin vectors are four-dimensional vectors of their own (nn, aX, bX) whose scalar products with every momentum are set from the definitions above: defineAxes[lIn, kIn, pOut, sIn, tIn, uIn, momenta] expands aX.x and bX.x for the listed momenta x through the scalar products of lIn, kIn, pOut already in force, so that the same definitions serve the photon process (with the photon momentum q in place of l).*)


(* ::Input:: *)
FCSetDiracGammaScheme["BMHV"];
$LeviCivitaSign = -1;
couplingNames = {SMP["e"] -> ee, SMP["g_s"] -> gs};
normalProducts[momenta_List] := (SPD[nn] = -1; SP[nn] = -1; Scan[(SPD[nn, #] = 0; SP[nn, #] = 0) &, momenta]; Momentum[nn, D - 4] = 0);
defineAxes[lIn_, kIn_, pOut_, sIn_, tIn_, uIn_, momenta_List] := Module[{base = {lIn, kIn, pOut}, ca, cb, product},
  ca = {tIn/sIn, uIn/sIn, 1};                                                        (* aX = (t/s) l + (u/s) k + p *)
  cb = {tIn/(tIn + uIn), -uIn/(tIn + uIn), -sIn (uIn - tIn)/(tIn + uIn)^2};         (* bX = t/(t+u) l - u/(t+u) k - s (u-t)/(t+u)^2 p *)
  product[c1_, c2_] := Factor[Sum[c1[[i]] c2[[j]] SP[base[[i]], base[[j]]], {i, 3}, {j, 3}]];
  Scan[(Momentum[#, D - 4] = 0) &, {aX, bX}];
  Do[SP[aX, x] = SPD[aX, x] = Factor[Sum[ca[[i]] SP[base[[i]], x], {i, 3}]];
     SP[bX, x] = SPD[bX, x] = Factor[Sum[cb[[i]] SP[base[[i]], x], {i, 3}]], {x, momenta}];
  SP[aX, nn] = SPD[aX, nn] = 0; SP[bX, nn] = SPD[bX, nn] = 0;
  SP[aX] = SPD[aX] = product[ca, ca]; SP[bX] = SPD[bX] = product[cb, cb]; SP[aX, bX] = SPD[aX, bX] = product[ca, cb];
  <|"aX" -> aX, "bX" -> bX, "Lengths" -> -tIn uIn/(tIn + uIn)|>];
bornKinematics := (FCClearScalarProducts[];
  SetMandelstam[s, t, u, l, k, -lp, -p, 0, 0, 0, 0];
  Scan[(Momentum[#, D - 4] = 0) &, {l, k, lp, p}];
  normalProducts[{l, k, p, lp}]);
realKinematics := (FCClearScalarProducts[];
  Scan[(Momentum[#, D - 4] = 0) &, {l, k, p, pRec}];
  SPD[l] = 0; SPD[k] = 0; SPD[lp] = 0; SPD[p] = 0; SPD[r] = 0;
  SPD[l, k] = s/2; SPD[k, p] = -t/2; SPD[l, p] = -u/2; SPD[k, r] = -a/2; SPD[p, r] = b/2;
  SPD[l, lp] = -(a + b + t)/2; SPD[l, r] = (s + t + u + a + b)/2; SPD[k, lp] = (s + t + a)/2;
  SPD[lp, p] = -(t + u + b)/2; SPD[lp, r] = (s + t + u)/2;
  SPE[lp] = spe3; SPE[r] = spe3; SPE[lp, r] = -spe3;
  SP[l] = 0; SP[k] = 0; SP[p] = 0; SP[lp] = -spe3; SP[r] = -spe3;
  SP[l, k] = s/2; SP[k, p] = -t/2; SP[l, p] = -u/2; SP[k, r] = -a/2; SP[p, r] = b/2;
  SP[l, lp] = -(a + b + t)/2; SP[l, r] = (s + t + u + a + b)/2; SP[k, lp] = (s + t + a)/2;
  SP[lp, p] = -(t + u + b)/2; SP[lp, r] = (s + t + u)/2 + spe3;
  normalProducts[{l, k, p}]; SPD[nn, r] = nr; SP[nn, r] = nr; SPD[nn, lp] = -nr; SP[nn, lp] = -nr);
photonKinematics := (FCClearScalarProducts[];
  Scan[(Momentum[#, D - 4] = 0) &, {q, k, p, r}];
  SPD[q] = 0; SPD[k] = 0; SPD[p] = 0; SPD[r] = 0;
  SPD[q, k] = sh/2; SPD[k, p] = -th/2; SPD[q, p] = -uh/2; SPD[p, r] = sh/2; SPD[k, r] = -uh/2; SPD[q, r] = -th/2;
  SP[q] = 0; SP[k] = 0; SP[p] = 0; SP[r] = 0;
  SP[q, k] = sh/2; SP[k, p] = -th/2; SP[q, p] = -uh/2; SP[p, r] = sh/2; SP[k, r] = -uh/2; SP[q, r] = -th/2;
  normalProducts[{q, k, p, r}]);
bornKinematics;
bornAxes = defineAxes[l, k, p, s, t, u, {l, k, p, lp}];
Require[Simplify[{SP[aX, l], SP[aX, k], SP[bX, p], SP[bX, l] + SP[bX, k], SP[aX] + t u/s, SP[bX] + s t u/(t + u)^2} /. u -> -s - t] === {0, 0, 0, 0, 0, 0},
  "The in-plane axes are not what the text says"];
{SP[aX], SP[bX], SP[aX, bX]}


(* ::Section:: *)
(*Squares from explicit traces*)


(* ::Text:: *)
(*For each process the lepton trace Tr[lp-slash gamma^mu l-slash gamma^nu]/2 (or, for the photon process, the photon polarization sum) multiplies the quark trace Tr[rho_out W_i rho_in W-bar_j], where W_i are the strings of the diagrams, rho_in = k-slash/2 (average) or k-slash gamma5 a-slash/2 and rho_out = p-slash (sum) or p-slash gamma5 b-slash.  spinCase is "Unpolarized", "NN" or "XX".  The strings (photon index mu, gluon index rho, all D-dimensional):*)
(**)
(*  Born:   W = gamma^mu / t;*)
(*  Real:   W_1 = gamma^mu (k-slash - r-slash) gamma^rho / a,   W_2 = gamma^rho (p-slash + r-slash) gamma^mu / b,   times 1/tau^2, gluon sum -g;*)
(*  Photon: W_1 = gamma^rho (k-slash + q-slash) gamma^mu / sh,   W_2 = gamma^mu (k-slash - r-slash) gamma^rho / uh,   photon sum -g/(D-2), gluon sum -g.*)
(**)
(*Colour CF = Tr[T^a T^a]/Nc (Born: 1).  As in Processes/LeptonOffDiagonal, the unpolarized case built from the same strings is checked against FeynCalc's generic spin sum, and the amplitude records of step 1 show that the strings are those of the diagrams.*)


(* ::Input:: *)
reverse[string_] := (string /. {mu -> nu, rho -> sigma}) /. Dot -> (Dot @@ Reverse[{##}] &);
simplifyTrace[expr_] := Factor[Together[Simplify[Contract[ExpandScalarProduct[Contract[DiracSimplify[expr, DiracTraceEvaluate -> True]]]]]]];
density[spinCase_, aIn_, bOut_] := Switch[spinCase,
  "Unpolarized", {GS[k]/2, GS[p]},
  "NN", {GS[k].GA[5].GS[nn]/2, GS[p].GA[5].GS[nn]},
  "XX", {GS[k].GA[5].GS[aIn]/2, GS[p].GA[5].GS[bOut]}];
leptonTrace = DiracTrace[GSD[lp].GAD[mu].GSD[l].GAD[nu]]/2;
quarkTrace[strings_, spinCase_, ax_] := With[{rho = density[spinCase, ax["aX"], ax["bX"]]},
  Table[DiracTrace[rho[[2]].strings[[i]].rho[[1]].reverse[strings[[j]]]], {i, Length[strings]}, {j, Length[strings]}]];
lengths[spinCase_, ax_] := If[spinCase === "XX", ax["Lengths"], 1];


(* ::Section:: *)
(*Born: l q -> l q*)


(* ::Text:: *)
(*Expected: the unpolarized square 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - eps) as in the lepton-quark process, and the transverse tensors ee^4 eq^2 [-4 s u/t^2 -+ (D-4)] (Eq. (transverse-born-tensor) of the paper): equal in four dimensions, with opposite evanescent terms, so that the Collins projection has none.*)


(* ::Input:: *)
bornKinematics; bornAxes = defineAxes[l, k, p, s, t, u, {l, k, p, lp}];
bornGeneric = Module[{amp = Total[amplitudes["Born"]]},
  amp ComplexConjugate[amp]/SUNN // FeynAmpDenominatorExplicit // SUNSimplify[#, Explicit -> True] & //
    FermionSpinSum[#, ExtraFactor -> 1/4] & // DiracSimplify[#, DiracTraceEvaluate -> True] & // Contract //
    TrickMandelstam[#, {s, t, u, 0}] & // Factor // ReplaceAll[#, couplingNames] &];
bornByTraces[spinCase_] := ee^4 eq^2/lengths[spinCase, bornAxes] TrickMandelstam[Factor[Together[
  Total[simplifyTrace[leptonTrace #] & /@ Flatten[quarkTrace[{GAD[mu]/t}, spinCase, bornAxes]]]]], {s, t, u, 0}];
Require[Simplify[bornByTraces["Unpolarized"] - bornGeneric] === 0, "The explicit traces do not reproduce the unpolarized Born square"];
born = <|"Unpolarized" -> bornGeneric, "NN" -> bornByTraces["NN"], "XX" -> bornByTraces["XX"]|>;
born["Transverse"] = Factor[(born["NN"] + born["XX"])/2];
Require[Simplify[(born["Unpolarized"] - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];
Require[Simplify[{born["NN"] - ee^4 eq^2 (-4 s u/t^2 - (D - 4)), born["XX"] - ee^4 eq^2 (-4 s u/t^2 + (D - 4))} /. u -> -s - t] === {0, 0},
  "Transverse Born tensors differ from Eq. (transverse-born-tensor) of the paper"];
born


(* ::Section:: *)
(*Virtual: the one-loop vertex correction*)


(* ::Text:: *)
(*The vertex graph gives a vector form factor times the Born current (the massless on-shell vertex has no other structure, see Processes/LeptonQuark), so the interference 2 Re(M_vertex M_Born^*) is the same multiple of the Born square for every spin projection.  The scalar family VirtualIntegral[n0,n1,n2] is extracted exactly as in the lepton-quark process.*)


(* ::Input:: *)
FCSetDiracGammaScheme["NDR"];
FCClearScalarProducts[];
SetMandelstam[s, t, u, l, k, -lp, -p, 0, 0, 0, 0];
$KeepLogDivergentScalelessIntegrals = True;
bornAmplitude = Total[amplitudes["Born"]];
vertexAmplitude = Total[amplitudes["VirtualVertex"]] // DiracSimplify // Contract // SUNSimplify[#, SUNNToCACF -> True] &;
vertexScalar = vertexAmplitude // TID[#, loop, ToPaVe -> True] & // DiracSimplify // Contract //
  SUNSimplify[#, SUNNToCACF -> True] & // FeynAmpDenominatorExplicit;
bornReduced = bornAmplitude // DiracSimplify // Contract // FeynAmpDenominatorExplicit;
square[x_, y_] := x ComplexConjugate[y] // FermionSpinSum[#, ExtraFactor -> 1/4] & // DiracSimplify //
  SUNSimplify[#, SUNNToCACF -> True] & // TrickMandelstam[#, {s, t, u, 0}] & // Simplify;
virtualRatio = Factor[square[2 vertexScalar, bornReduced]/square[bornReduced, bornReduced]];
virtualFamily = Factor[virtualRatio/(2 CF Pi^2 SMP["g_s"]^2)] /. {
  B0[0, 0, 0] -> VirtualIntegral[1, 1, 0], B0[t, 0, 0] -> VirtualIntegral[0, 1, 1],
  C0[0, 0, t, 0, 0, 0] -> VirtualIntegral[1, 1, 1]}


(* ::Section:: *)
(*Real emission: l q -> l q g*)


(* ::Text:: *)
(*The tensors are polynomials in nr; the odd part is dropped (it integrates to zero) after checking that it is at most linear, and nr^2 is replaced as explained above.  The unpolarized square is also computed with FeynCalc's generic tools and compared.*)


(* ::Input:: *)
FCSetDiracGammaScheme["BMHV"];
realKinematics;
realAxes = defineAxes[l, k, p, s, t, u, {l, k, p, lp, r}];
realStrings = {GAD[mu].(GSD[k] - GSD[r]).GAD[rho]/a, GAD[rho].(GSD[p] + GSD[r]).GAD[mu]/b};
realByTraces[spinCase_] := CF ee^4 eq^2 gs^2/((a + b + t)^2 lengths[spinCase, realAxes]) Factor[Together[
  Total[simplifyTrace[leptonTrace # (-MTD[rho, sigma])] & /@ Flatten[quarkTrace[realStrings, spinCase, realAxes]]]]];
gram = Table[SP[x, y], {x, {l, k, p}}, {y, {l, k, p}}];
parallelSquare = Factor[{SP[l, r], SP[k, r], SP[p, r]}.Inverse[gram].{SP[l, r], SP[k, r], SP[p, r]}];       (* rpar^2 *)
normalReduction[expr_] := Module[{coefficients = CoefficientList[Expand[expr], nr]},
  Require[Length[coefficients] <= 3, "The real-emission tensor has a power of nr above two", Length[coefficients] - 1];
  Factor[coefficients[[1]] + If[Length[coefficients] == 3, coefficients[[3]] (parallelSquare + spe3), 0]]];
realUnpolarized = Module[{amp = Total[SUNSimplify[#, SUNNToCACF -> True] & /@ amplitudes["Real"]] // DiracSimplify // Contract},
  amp ComplexConjugate[amp]/SUNN // FeynAmpDenominatorExplicit // SUNSimplify[#, Explicit -> True, SUNNToCACF -> True] & //
    FermionSpinSum[#, ExtraFactor -> 1/4] & // DiracSimplify // Contract // DoPolarizationSums[#, r, 0, VirtualBoson -> True] & //
    Contract // Simplify // SUNSimplify[#, SUNNToCACF -> True] & // Factor // ReplaceAll[#, couplingNames] &];
Require[Together[(realByTraces["Unpolarized"] /. spe3 -> 0) - realUnpolarized] === 0, "The explicit traces do not reproduce the unpolarized real square"];
real = <|"Unpolarized" -> realUnpolarized, "NN" -> normalReduction[realByTraces["NN"]], "XX" -> normalReduction[realByTraces["XX"]]|>;
real["Transverse"] = Factor[(real["NN"] + real["XX"])/2];
Require[FreeQ[real, nr | nn], "The normal component was not eliminated"];
{parallelSquare, Factor[(real["NN"] - real["XX"]) /. {D -> 4, spe3 -> 0}]}


(* ::Text:: *)
(*In four dimensions and without the hatted terms the two transverse components coincide (displayed above: zero), as for the Born tensor.*)


(* ::Input:: *)
Require[Factor[(real["NN"] - real["XX"]) /. {D -> 4, spe3 -> 0}] === 0, "The two transverse components differ in four dimensions"];
Factor[real["Transverse"] /. {D -> 4, spe3 -> 0}]


(* ::Section:: *)
(*Photon-quark Born: gamma q -> q g*)


(* ::Text:: *)
(*Needed through order eps for the photon-in-lepton subtraction (the photon is unpolarized).  Axes as above with the photon momentum q in place of l.  Expected transverse projection (Eq. (transverse-photon-born) of the paper, CF ee^2 eq^2 gs^2 removed): 4 - 2 eps^2 th^2/((1-eps) sh uh).*)


(* ::Input:: *)
photonKinematics;
photonAxes = defineAxes[q, k, p, sh, th, uh, {q, k, p, r}];
photonStrings = {GAD[rho].(GSD[k] + GSD[q]).GAD[mu]/sh, GAD[mu].(GSD[k] - GSD[r]).GAD[rho]/uh};
photonByTraces[spinCase_] := CF ee^2 eq^2 gs^2/lengths[spinCase, photonAxes] Factor[TrickMandelstam[Cancel[Together[
  Total[simplifyTrace[(-MTD[mu, nu]/(D - 2)) # (-MTD[rho, sigma])] & /@ Flatten[quarkTrace[photonStrings, spinCase, photonAxes]]]]], {sh, th, uh, 0}]];
photonGeneric = Module[{amp = Total[SUNSimplify[#, SUNNToCACF -> True] & /@ amplitudes["Photon"]] // DiracSimplify // Contract, expr},
  expr = amp ComplexConjugate[amp]/SUNN // FeynAmpDenominatorExplicit // SUNSimplify[#, Explicit -> True, SUNNToCACF -> True] & //
    FermionSpinSum[#, ExtraFactor -> 1/2] & // DiracSimplify // Contract // DoPolarizationSums[#, r, 0, VirtualBoson -> True] & //
    DoPolarizationSums[#, q, 0, VirtualBoson -> True] & // Contract // Simplify // SUNSimplify[#, SUNNToCACF -> True] & // Factor;
  TrickMandelstam[expr/(D - 2), {sh, th, uh, 0}] /. couplingNames];
Require[Simplify[(photonByTraces["Unpolarized"] - photonGeneric) /. uh -> -sh - th] === 0, "The explicit traces do not reproduce the unpolarized photon square"];
photon = <|"Unpolarized" -> photonGeneric, "NN" -> photonByTraces["NN"], "XX" -> photonByTraces["XX"]|>;
photon["Transverse"] = Factor[TrickMandelstam[(photon["NN"] + photon["XX"])/2, {sh, th, uh, 0}]];
Require[Simplify[(photon["Transverse"]/(CF ee^2 eq^2 gs^2) - (4 - 2 eps^2 th^2/((1 - eps) sh uh))) /. D -> 4 - 2 eps /. uh -> -sh - th] === 0,
  "The photon transverse tensor differs from Eq. (transverse-photon-born) of the paper"];
photon


(* ::Section:: *)
(*Save*)


(* ::Text:: *)
(*Algebra.wl holds everything step 3 needs, in the format of the lepton-quark process: the squares with the keys "Unpolarized" and "Transverse" (and the two components "NN", "XX" for inspection) and the virtual family.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Photon" -> photon, "VirtualFamily" -> virtualFamily|>;
SaveResult["Algebra", algebra, "FeynCalc"];
Keys[algebra]
