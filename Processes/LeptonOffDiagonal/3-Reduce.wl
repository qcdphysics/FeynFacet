(* ::Title:: *)
(*Lepton-quark and lepton-gluon scattering with a different tagged parton, step 3: phase-space integrals*)


(* ::Text:: *)
(*Input: Results/<route>/Algebra.wl from step 2 (choose the route below).  Output: Results/<route>/Reduced.wl with the real-emission phase-space integrals of both channels evaluated as Laurent series in eps, exact in the threshold variable.  Both spin cases are treated together.*)
(**)
(*Method, as in Processes/LeptonQuark/3-Reduce.wl (read its introduction first): the unobserved two-particle phase space of lp and r is written as cut propagators, Kira reduces the integrals to the two-body volume and one angular master per family, and the masters are evaluated as Euler integrals with SubTropica.  The denominators of the squared matrix elements are, besides the photon propagator tau = (k - p - r)^2 = a + b + t,*)
(**)
(*  l q -> l g q:    t = (k - p)^2 (fixed) and b = (p + r)^2:   family real_b_tau[1,1,nb,nt] = integral of 1/(b^nb tau^nt),*)
(*  l g -> l q qbar: t and a = (k - r)^2:                           family real_a_tau[1,1,na,nt] = integral of 1/(a^na tau^nt),*)
(**)
(*(indices 1 and 2 are the cut propagators), so no partial fractions are needed: each channel lives in one family.  The invariant a (or b) that is not a denominator of the channel is written as tau - b - t (or tau - a - t).  The hatted invariant spe3 of the polarized matrix elements becomes a dimension shift as in the lepton-quark process.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";        (* "FeynCalc", "FORM" or "Tracer": which step-2 result to reduce *)
algebra = LoadResult["Algebra", route];
spins = {"Unpolarized", "Polarized"};
channels = {"QG", "GQ"};
Keys[algebra]


(* ::Section:: *)
(*Real emission: terms of the two families*)


(* ::Text:: *)
(*Remove the common factors (CF or TR, ee^4 eq^2 gs^2), write D = 4 - 2 eps and u = s4 - s - t, express the matrix element through the two denominators of its family, and read off the integrals: a term x^(-i) y^(-j) with x, y the two denominators is the family integral with indices {1, 1, i, j} (negative indices are numerators).  The coefficients are rational in s, t, s4 and eps (the gluon average 1/(D-2) of the second channel makes them rational in eps) and polynomial in spe3.*)


(* ::Input:: *)
familyOf = <|"QG" -> {KRealBTau, b, tau, a -> tau - b - t}, "GQ" -> {KRealATau, a, tau, b -> tau - a - t}|>;
colourOf = <|"QG" -> CF, "GQ" -> TR|>;
(* The two channels already have two denominators. Map each additive fraction
   directly to two inert forms; keep its scalar denominator independent. *)
toFamilyForms[expr_, {head_, x_, y_, rule_}] := Module[{mapped},
  mapped = MapIntegrationDenominators[expr, {rule}, {x, y}, Function[factor,
    Which[Expand[factor - x] === 0, fX, Expand[factor + x] === 0, -fX,
      Expand[factor - y] === 0, fY, Expand[factor + y] === 0, -fY,
      True, Require[False, "A denominator is not one of the two family forms", factor]]]];
  Require[!FailureQ[mapped], "Termwise family conversion failed", mapped];
  mapped /. {x -> fX, y -> fY}];
rowVariables = {fX, fY, spe3};
reconstructUnit[term_, rows_List] := Together[Total[
  #["Coefficient"] spe3^#["HatPower"] fX^(-#["Integral"][[3]]) fY^(-#["Integral"][[4]]) & /@ rows] - term] === 0;
makeFamilyConverter[head_] := With[{familyHead = head},
  CreateIntegralRowConverter[rowVariables, Function[unit, Module[{powers, rows},
    powers = Exponent[unit, #] & /@ rowVariables;
    Require[VectorQ[powers, IntegerQ] && Last[powers] >= 0, "Unsupported family monomial", unit];
    rows = {<|"Integral" -> familyHead[1, 1, -powers[[1]], -powers[[2]]],
      "HatPower" -> powers[[3]], "Coefficient" -> 1|>};
    Require[reconstructUnit[unit, rows], "A unit family mapping does not reconstruct", unit]; rows]], reconstructUnit]];
familyConverters = Association @ Table[channel -> makeFamilyConverter[familyOf[channel][[1]]], {channel, channels}];
realStripped = Association @ Table[channel -> Association @ Table[spin ->
  Expand[algebra["Real"][channel][spin]/(colourOf[channel] ee^4 eq^2 gs^2) /.
    {D -> 4 - 2 eps, u -> s4 - s - t}], {spin, spins}], {channel, channels}];
rowInputs = Association @ Table[channel -> Association @ Table[spin ->
  CombineLaurentTerms[toFamilyForms[realStripped[channel][spin], familyOf[channel]], rowVariables],
  {spin, spins}], {channel, channels}];
Require[FreeQ[rowInputs, _Failure], "Laurent coefficient grouping failed", rowInputs];
realRowResults = Association @ Table[channel -> Association @ Table[spin ->
  familyConverters[channel][rowInputs[channel][spin]], {spin, spins}], {channel, channels}];
Require[FreeQ[realRowResults, _Failure], "Integral-row conversion failed", realRowResults];
realTerms = Association @ Table[channel -> Association @ Table[spin ->
  realRowResults[channel][spin]["Rows"], {spin, spins}], {channel, channels}];
allRows = Flatten[Table[realTerms[channel][spin], {channel, channels}, {spin, spins}]];
indexRanges = <|"denominators" -> Max[Prepend[Total[Select[List @@ #["Integral"], # > 0 &]] & /@ allRows, 0]],
  "numerators" -> Max[Prepend[-Total[Select[List @@ #["Integral"], # < 0 &]] & /@ allRows, 0]]|>;
Require[indexRanges["denominators"] <= 7 && indexRanges["numerators"] <= 3,
  "The integrals exceed the shipped Kira index ranges", indexRanges];
Dataset[realTerms["QG"]["Polarized"]]


(* ::Text:: *)
(*Check that the terms add up to the matrix element.*)


(* ::Input:: *)
rebuild[rows_, {head_, x_, y_, rule_}] := Total[#["Coefficient"] spe3^#["HatPower"] x^(-#["Integral"][[3]]) y^(-#["Integral"][[4]]) & /@ rows];
Do[Require[Together[(realStripped[channel][spin] /. familyOf[channel][[4]]) - rebuild[realTerms[channel][spin], familyOf[channel]]] === 0,
  "The family terms do not add up to the matrix element", {channel, spin}], {channel, channels}, {spin, spins}];
"family terms verified"


(* ::Section:: *)
(*Reduction with Kira*)


(* ::Text:: *)
(*The YAML files in kira/real define the two families (the same definitions as in the lepton-quark process); the target list is written from the integrals that occur.  Kira sets s = 1 (symbol_to_replace_by_one); RestoreScale puts the scale back.*)


(* ::Input:: *)
kiraNames = {KRealATau -> "real_a_tau", KRealBTau -> "real_b_tau"};
referenceRealMasters = {KRealATau[1, 1, 0, 0], KRealBTau[1, 1, 0, 0], KRealATau[1, 1, 1, 1], KRealBTau[1, 1, 1, 1]};
realTargets = Union[Flatten[Table[Lookup[realTerms[channel][spin], "Integral"], {channel, channels}, {spin, spins}]], referenceRealMasters];
kiraReal = PrepareKira[FileNameJoin[{NotebookDirectory[], "kira", "real"}], FileNameJoin[{ResultsDirectory[route], "kira", "real"}]];
WriteKiraTargets[kiraReal, realTargets, kiraNames];
realTargets


(* ::Input:: *)
RunKira[kiraReal];
realRules = RestoreScale[Join @@ Table[KiraRules[kiraReal, name, kiraNames],
  {name, Select[{"real_a_tau", "real_b_tau"},
    FileExistsQ[FileNameJoin[{kiraReal, "results", #, "kira_targets.m"}]] &]}],
  s, {t -> t/s, s4 -> s4/s}, {KRealATau, KRealBTau}];
CheckReductionCoverage[realTargets, realRules, {KRealATau, KRealBTau}, referenceRealMasters];
Column[realRules]


(* ::Text:: *)
(*Every integral reduces to the phase-space volume (real_a_tau[1,1,0,0] or real_b_tau[1,1,0,0], the same integral) and the angular integrals real_a_tau[1,1,1,1] and real_b_tau[1,1,1,1].  A term with spe3^n uses the rules at d = 4 + 2n - 2 eps; KiraMaster[family, indices, shift] denotes a master in D + shift dimensions.*)


(* ::Input:: *)
masterHeads = {KRealATau -> "real_a_tau", KRealBTau -> "real_b_tau"};
reduce[rows_] := Collect[Expand[Total[Map[Function[row, With[{n = row["HatPower"]},
  row["Coefficient"] (4 Pi)^n (-1)^n Pochhammer[-eps, n] *
    ((row["Integral"] /. realRules /. d -> 4 + 2 n - 2 eps) /.
      (head : KRealATau | KRealBTau)[indices__] :> KiraMaster[head /. masterHeads, {indices}, 2 n])]], rows]]], _KiraMaster, Together];
realReduced = Association @ Table[channel -> Association @ Table[spin -> reduce[realTerms[channel][spin]], {spin, spins}], {channel, channels}];
Require[FreeQ[realReduced, KRealATau | KRealBTau], "Some integrals were not reduced"];
Union[Cases[realReduced, _KiraMaster, Infinity]]


(* ::Section:: *)
(*Master integrals*)


(* ::Text:: *)
(*The same masters as in the lepton-quark process, with the same checks: the two-body volume Phi2 (and its D + 2 version) and the angular integral*)
(**)
(*   Omega(z) = 2 (1-2eps) integral_0^infinity dx x^(-1-eps) (1+x)^eps/(z+x),*)
(**)
(*written with the connection formula that separates the branch z^(-1-eps) from the regular hypergeometric function F1(z) = 2F1(1,1;2+eps;z) (and F2 in D + 2 dimensions):*)
(**)
(*   Omega(z) = -2 (1-2eps)/eps [ Gamma(1-eps) Gamma(1+eps) (1-z)^eps z^(-1-eps) + Gamma(1-eps) Gamma(-1-eps)/Gamma(-eps)^2 F1(z) ],*)
(*   Omega2(z) = 2 (3-2eps)/(1-eps) [ Gamma(2-eps) Gamma(eps) (1-z)^(eps-1) z^(-eps) + Gamma(2-eps) Gamma(-eps)/Gamma(1-eps)^2 F2(z) ].*)
(**)
(*In these channels the cross ratios zAT and zBT tend to 1 at the threshold w = 1 (there is no soft-gluon singularity: the gluon is observed or absent), so no plus distributions arise; the branch (1-z)^eps combines with s4^(-eps) of the volume into logarithms of 1 - w.  F1 and F2 are evaluated as Euler integrals with SubTropica and kept as inert symbols F1[n][z], F2[n][z] (the coefficients of eps^n) until step 4.*)


(* ::Input:: *)
LoadSubTropica[]


(* ::Input:: *)
masters = <|
  "phaseSpace" -> EulerIntegral[1/Gamma[1 - eps], x^(-eps) (1 + x)^(2 eps - 2), {x}, {}, 2, eps],
  "phaseSpaceD2" -> EulerIntegral[1/Gamma[2 - eps], x^(1 - eps) (1 + x)^(2 eps - 4), {x}, {}, 1, eps],
  "F1" -> EulerIntegral[1 + eps, (1 + y)^(-1 - eps)/(1 + zb y), {y}, {zb}, 2, eps] /. zb -> 1 - z,
  "F2" -> EulerIntegral[eps, (1 + y)^(-eps)/(1 + zb y), {y}, {zb}, 2, eps] /. zb -> 1 - z|>


(* ::Input:: *)
check[x_] := Require[Simplify[x] === 0, "Master integral check failed", x];
check[Normal[Series[masters["phaseSpace"] - Gamma[1 - eps]/Gamma[2 - 2 eps], {eps, 0, 2}]]];
check[Normal[Series[masters["phaseSpaceD2"] - Gamma[2 - eps]/Gamma[4 - 2 eps], {eps, 0, 1}]]];
truncationCheck[series_, exact_, power_] := Module[{d1, d2},
  d1 = N[(exact - series) /. {eps -> -1/50, z -> 1/3}, 30];
  d2 = N[(exact - series) /. {eps -> -1/100, z -> 1/3}, 30];
  Require[Abs[d2/d1 - 2^-power] < 2^-power/10, "Master integral disagrees with its closed form", {d1, d2}]; {d1, d2}];
truncationCheck[masters["F1"], Hypergeometric2F1[1, 1, 2 + eps, z], 3];
truncationCheck[masters["F2"], Hypergeometric2F1[1, 1, 1 + eps, z], 3]


(* ::Input:: *)
angular[zz_] := -2 (1 - 2 eps)/eps (Gamma[1 - eps] Gamma[1 + eps] (1 - zz)^eps zz^(-1 - eps) +
  Gamma[1 - eps] Gamma[-1 - eps]/Gamma[-eps]^2 Sum[eps^n F1[n][zz], {n, 0, 2}]);
angularD2[zz_] := 2 (3 - 2 eps)/(1 - eps) (Gamma[2 - eps] Gamma[eps] (1 - zz)^(eps - 1) zz^(-eps) +
  Gamma[2 - eps] Gamma[-eps]/Gamma[1 - eps]^2 Sum[eps^n F2[n][zz], {n, 0, 2}]);
masterFunctions = Join[{F1[0][0] -> 1, F2[0][0] -> 1, F1[_][0] -> 0, F2[_][0] -> 0},
  Table[F1[n][zz_] :> Evaluate[SeriesCoefficient[masters["F1"], {eps, 0, n}] /. z -> zz], {n, 0, 2}],
  Table[F2[n][zz_] :> Evaluate[SeriesCoefficient[masters["F2"], {eps, 0, n}] /. z -> zz], {n, 0, 2}]];
numericCheck[expr_, integrand_, power_] := Module[{d1, d2},
  d1 = N[expr /. {eps -> -1/40, z -> 1/3}, 20] - NIntegrate[integrand /. {eps -> -1/40, z -> 1/3}, {x, 0, Infinity}, WorkingPrecision -> 30];
  d2 = N[expr /. {eps -> -1/80, z -> 1/3}, 20] - NIntegrate[integrand /. {eps -> -1/80, z -> 1/3}, {x, 0, Infinity}, WorkingPrecision -> 30];
  Require[Abs[d2/d1 - 2^-power] < 2^-power/10, "Angular master disagrees with direct integration", {d1, d2}]; {d1, d2}];
numericCheck[angular[z] /. masterFunctions, 2 (1 - 2 eps) x^(-1 - eps) (1 + x)^eps/(z + x), 3];
numericCheck[angularD2[z] /. masterFunctions, 2 (3 - 2 eps) x^(-eps) (1 + x)^(eps - 1)/(z + x), 2]


(* ::Section:: *)
(*Insert the physical masters*)


(* ::Text:: *)
(*Normalization as in the lepton-quark process: Phi2 = (4 Pi)^eps s4^(-eps)/(8 Pi) Gamma(1-eps)/Gamma(2-2eps) for the measure d^D r/(2 Pi)^(D-1) delta+(r^2) 2 Pi delta+((pRec - r)^2), Phi2' = (4 Pi)^(eps-1) s4^(1-eps)/(8 Pi) Gamma(2-eps)/Gamma(4-2eps) in D + 2 dimensions, and the angular masters Phi2 Omega(z) divided by the product of the two invariants at the symmetric point, with the cross ratios*)
(**)
(*   zAT = -t (s+t-s4)/((s+t)(s4-t)),   zBT = -s t/((s-s4)(s4-t)).*)


(* ::Input:: *)
phaseSpace = (4 Pi)^eps s4^(-eps)/(8 Pi) masters["phaseSpace"];
phaseSpaceD2 = (4 Pi)^(eps - 1) s4^(1 - eps)/(8 Pi) masters["phaseSpaceD2"];
zAT = -t (s + t - s4)/((s + t) (s4 - t));  zBT = -s t/((s - s4) (s4 - t));
physicalMasters = {
  KiraMaster["real_a_tau", {1, 1, 0, 0}, 0] -> phaseSpace,
  KiraMaster["real_b_tau", {1, 1, 0, 0}, 0] -> phaseSpace,
  KiraMaster["real_a_tau", {1, 1, 1, 1}, 0] -> phaseSpace angular[zAT]/((s + t) (s4 - t)),
  KiraMaster["real_b_tau", {1, 1, 1, 1}, 0] -> -phaseSpace angular[zBT]/((s - s4) (s4 - t)),
  KiraMaster["real_a_tau", {1, 1, 0, 0}, 2] -> phaseSpaceD2,
  KiraMaster["real_b_tau", {1, 1, 0, 0}, 2] -> phaseSpaceD2,
  KiraMaster["real_a_tau", {1, 1, 1, 1}, 2] -> phaseSpaceD2 angularD2[zAT]/((s + t) (s4 - t)),
  KiraMaster["real_b_tau", {1, 1, 1, 1}, 2] -> -phaseSpaceD2 angularD2[zBT]/((s - s4) (s4 - t))};
realIntegrated = Association @ Table[channel -> Association @ Table[spin ->
  colourOf[channel] ee^4 eq^2 gs^2 (realReduced[channel][spin] /. physicalMasters), {spin, spins}], {channel, channels}];
Require[FreeQ[realIntegrated, KiraMaster], "A master integral was not inserted"];
Keys[realIntegrated]


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
reduced = <|"Real" -> realIntegrated, "Masters" -> masters, "MasterFunctions" -> masterFunctions, "RealRules" -> realRules, "Route" -> route|>;
SaveResult["Reduced", reduced, route];
Keys[reduced]
