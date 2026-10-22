(* ::Title:: *)
(*Lepton-quark scattering, step 3: scalar integrals*)


(* ::Text:: *)
(*Input: Results/<route>/Algebra.wl from step 2 (choose the route below).  Output: Results/<route>/Reduced.wl with the real-emission phase-space integrals and the virtual loop integral evaluated as Laurent series in eps, still exact in the threshold variable.  Both spin cases are treated together; step 4 turns them into the finite hard coefficients.*)
(**)
(*Method.  The unobserved two-particle phase space of lp and r is written as cut propagators (reverse unitarity), so that integration-by-parts reduction with Kira applies.  With pRec = l + k - p = lp + r, pRec^2 = s4, the cut propagators are r^2 and (pRec - r)^2 = lp^2.  The uncut denominators of the real matrix element are a = (k - r)^2, b = (p + r)^2 and tau = (k - p - r)^2 = a + b + t.  Since tau depends on a and b, only two of the three can appear in one integral; partial fractions give three families:*)
(**)
(*   real_ab[1,1,n3,n4]    = integral of 1/(a^n3 b^n4)     over the cut,*)
(*   real_a_tau[1,1,n3,n4] = integral of 1/(a^n3 tau^n4),*)
(*   real_b_tau[1,1,n3,n4] = integral of 1/(b^n3 tau^n4),*)
(**)
(*defined in kira/real/config/integralfamilies.yaml (indices 1 and 2 are the two cut propagators and always equal 1).  The virtual family is VirtualIntegral[n0,n1,n2] with denominators loop^2, (loop-k)^2, (loop-p)^2, in kira/virtual.  Kira reduces every integral to four phase-space masters (the two-body volume and the angular integral in three kinematic configurations) and one loop master (the off-shell bubble).  Their values follow from one-dimensional Euler integrals evaluated with SubTropica.*)
(**)
(*Dimension shifts.  The polarized matrix element contains the hatted invariant spe3 = rhat.rhat.  A power spe3^n under the phase-space integral equals, up to the factor (4 Pi)^n (-1)^n (-eps)_n, the same integral in D + 2n dimensions ((x)_n is the Pochhammer symbol).  Kira's reduction holds for any dimension d, so the shifted integrals are obtained by setting d -> 4 + 2n - 2 eps in its rules; the masters in D + 2 are evaluated separately.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";        (* "FeynCalc", "FORM" or "Tracer": which step-2 result to reduce *)
algebra = LoadResult["Algebra", route];
spins = {"Unpolarized", "Polarized"};
Keys[algebra]


(* ::Section:: *)
(*Real emission: partial fractions into the three families*)


(* ::Text:: *)
(*Remove the common factor CF ee^4 eq^2 gs^2, write D = 4 + dm4, and decompose the squared matrix element into terms with at most two of the denominators a, b, tau.  The coefficients are polynomials in dm4 and spe3 and rational in s, t, s4 (u = s4 - s - t). Equal Laurent monomials are combined before assigning integral families. The partial-fraction identity is applied once to each propagator geometry, with scalar coefficients left outside.*)


(* ::Input:: *)
forms = <|fA -> a, fB -> b, fT -> a + b + t|>;
formOf[factor_] := Module[{matches = Select[Keys[forms],
    Expand[factor - forms[#]] === 0 || Expand[factor + forms[#]] === 0 &]},
  Require[Length[matches] == 1, "A denominator factor is not a, b or tau", factor];
  First[matches] If[Expand[factor - forms[First[matches]]] === 0, 1, -1]];
toForms[expr_] := Module[{mapped = MapIntegrationDenominators[expr,
    {tau -> a + b + t}, {a, b}, formOf]},
  Require[!FailureQ[mapped], "Termwise denominator conversion failed", mapped];
  mapped /. {a -> fA, b -> fB}];
tripleRule = fA^i_. fB^j_. fT^k_. /; i < 0 && j < 0 && k < 0 :>
  fA^(i + 1) fB^(j + 1) fT^(k + 1) (1/(fA fB) - 1/(fA fT) - 1/(fB fT))/t;
termList[expr_] := With[{e = Expand[expr]}, Which[e === 0, {}, Head[e] === Plus, List @@ e, True, {e}]];
negativeForms[term_] := Select[Keys[forms], Exponent[term, #] < 0 &];
(* Only the unit propagator geometry is split; its scalar coefficient is a spectator. *)
partialFractionMonomial[unit_] := Module[{split},
  split = FixedPoint[Expand[# /. tripleRule] &, unit, 20];
  Require[AllTrue[termList[split], Length[negativeForms[#]] <= 2 &], "Partial fractions did not terminate", unit];
  Require[Together[(split - unit) /. forms] === 0, "A geometric partial fraction does not reconstruct", unit];
  split];
fractionConverter = CreatePartialFractionConverter[Keys[forms], spe3, partialFractionMonomial];
realStripped = Association @ Table[spin ->
  Expand[algebra["Real"][spin]/(CF ee^4 eq^2 gs^2) /. {D -> 4 + dm4, u -> s4 - s - t}], {spin, spins}];
realFractionResults = Association @ Table[spin -> fractionConverter[toForms[realStripped[spin]]], {spin, spins}];
Require[FreeQ[realFractionResults, _Failure], "Sparse partial fractions failed", realFractionResults];
realFractionTerms = Association @ Table[spin -> realFractionResults[spin]["Terms"], {spin, spins}];
Length /@ realFractionTerms


(* ::Text:: *)
(*Check that the decomposition reproduces the matrix element, then label each term by its family and index, separating the powers of dm4 and spe3.*)


(* ::Input:: *)
familyForms = <|KRealAB -> {fA, fB}, KRealATau -> {fA, fT}, KRealBTau -> {fB, fT}|>;
familyOf[denominators_List] := Which[SubsetQ[{fA, fB}, denominators], KRealAB,
  SubsetQ[{fA, fT}, denominators], KRealATau, SubsetQ[{fB, fT}, denominators], KRealBTau,
  True, Require[False, "A term has three denominators", denominators]];
familySolutions = Association @ KeyValueMap[Function[{head, pair},
  head -> First[Solve[{X == forms[pair[[1]]], Y == forms[pair[[2]]]}, {a, b}]]], familyForms];
integrand[row_] := row["Coefficient"] dm4^row["Dm4Power"] spe3^row["HatPower"] Times @@
  MapThread[#1^(-#2) &, {forms /@ familyForms[Head[row["Integral"]]], List @@ row["Integral"][[3 ;; 4]]}];
reconstructRows[term_, rows_List] := Together[Total[integrand /@ rows] - (term /. forms)] === 0;
unitRows[unit_] := Module[{den, head, pair, powers, numerator, poly, rows},
  den = negativeForms[unit]; head = familyOf[den]; pair = familyForms[head];
  powers = Association @ Table[f -> -Exponent[unit, f], {f, den}];
  numerator = unit Times @@ Table[f^(-Exponent[unit, f]), {f, den}] /. forms;
  poly = Expand[numerator /. familySolutions[head]];
  Require[PolynomialQ[poly, {X, Y, dm4, spe3}], "A numerator is not polynomial in its family", unit];
  rows = Table[<|"Integral" -> head[1, 1, Lookup[powers, pair[[1]], 0] - rule[[1, 1]],
      Lookup[powers, pair[[2]], 0] - rule[[1, 2]]], "Dm4Power" -> rule[[1, 3]],
    "HatPower" -> rule[[1, 4]], "Coefficient" -> rule[[2]]|>,
    {rule, CoefficientRules[poly, {X, Y, dm4, spe3}]}];
  Require[reconstructRows[unit, rows], "A unit integral mapping does not reconstruct", unit]; rows];
rowVariables = Join[Keys[forms], {dm4, spe3}];
rowConverter = CreateIntegralRowConverter[rowVariables, unitRows, reconstructRows];
(* Separate evanescent numerator powers before mapping, retaining trace D = 4 + dm4. *)
rowInputs = Association @ Table[spin -> CombineLaurentTerms[Total[realFractionTerms[spin]], rowVariables], {spin, spins}];
Require[FreeQ[rowInputs, _Failure], "Laurent coefficient grouping failed", rowInputs];
realRowResults = Association @ Table[spin -> rowConverter[rowInputs[spin]], {spin, spins}];
Require[FreeQ[realRowResults, _Failure], "Integral-row conversion failed", realRowResults];
realTerms = Association @ Table[spin -> realRowResults[spin]["Rows"], {spin, spins}];
indexRanges = <|"denominators" -> Max[Prepend[Total[Select[List @@ #["Integral"], # > 0 &]] & /@ Flatten[Values[realTerms]], 0]],
  "numerators" -> Max[Prepend[-Total[Select[List @@ #["Integral"], # < 0 &]] & /@ Flatten[Values[realTerms]], 0]]|>;
Require[indexRanges["denominators"] <= 6 && indexRanges["numerators"] <= 2,
  "The integrals exceed the shipped Kira index ranges", indexRanges];
(* Together verifies this identity; it is not used to build the integration rows. *)
Do[Require[Together[Total[integrand /@ realTerms[spin]] - (realStripped[spin] /. tau -> a + b + t)] === 0,
  "The family terms do not add up to the matrix element", spin], {spin, spins}];
Dataset[realTerms["Polarized"]]


(* ::Section:: *)
(*Reduction with Kira*)


(* ::Text:: *)
(*The YAML files define the families; the target list is written from the integrals that actually occur.  Kira sets s = 1 in the real families and t = 1 in the virtual family (symbol_to_replace_by_one); RestoreScale puts the scale back.  Kira's family names real_ab, real_a_tau, real_b_tau, virtual are spelled KRealAB, KRealATau, KRealBTau, KVirtual in Mathematica.*)


(* ::Input:: *)
kiraNames = {KRealAB -> "real_ab", KRealATau -> "real_a_tau", KRealBTau -> "real_b_tau", KVirtual -> "virtual"};
virtualFamily = algebra["VirtualFamily"] /. VirtualIntegral -> KVirtual;
referenceRealMasters = {KRealAB[1, 1, 0, 0], KRealAB[1, 1, 1, 1], KRealATau[1, 1, 1, 1], KRealBTau[1, 1, 1, 1]};
referenceVirtualMasters = {KVirtual[0, 1, 1]};
realTargets = Union[Flatten[Lookup[realTerms[#], "Integral"] & /@ spins], referenceRealMasters];
virtualTargets = Union[Cases[virtualFamily, _KVirtual, Infinity], referenceVirtualMasters];
kiraReal = PrepareKira[FileNameJoin[{NotebookDirectory[], "kira", "real"}], FileNameJoin[{ResultsDirectory[route], "kira", "real"}]];
kiraVirtual = PrepareKira[FileNameJoin[{NotebookDirectory[], "kira", "virtual"}], FileNameJoin[{ResultsDirectory[route], "kira", "virtual"}]];
WriteKiraTargets[kiraReal, realTargets, kiraNames];
WriteKiraTargets[kiraVirtual, virtualTargets, kiraNames];
{realTargets, virtualTargets}


(* ::Input:: *)
RunKira[kiraReal];
RunKira[kiraVirtual];
realRules = RestoreScale[Join @@ Table[KiraRules[kiraReal, name, kiraNames],
  {name, Select[{"real_ab", "real_a_tau", "real_b_tau"},
    FileExistsQ[FileNameJoin[{kiraReal, "results", #, "kira_targets.m"}]] &]}],
  s, {t -> t/s, s4 -> s4/s}, {KRealAB, KRealATau, KRealBTau}];
virtualRules = RestoreScale[KiraRules[kiraVirtual, "virtual", kiraNames], t, {}, {KVirtual}];
CheckReductionCoverage[realTargets, realRules, {KRealAB, KRealATau, KRealBTau}, referenceRealMasters];
CheckReductionCoverage[virtualTargets, virtualRules, {KVirtual}, referenceVirtualMasters];
Column[Join[virtualRules, realRules]]


(* ::Text:: *)
(*Every real integral reduces to real_ab[1,1,0,0] (the phase-space volume) and the three angular integrals real_ab[1,1,1,1], real_a_tau[1,1,1,1], real_b_tau[1,1,1,1]; the virtual family reduces to the off-shell bubble virtual[0,1,1].  The reduction rules are in the dimension d; a term with spe3^n uses them at d = 4 + 2n - 2 eps.  KiraMaster[family, indices, shift] denotes a master in D + shift dimensions.*)


(* ::Input:: *)
masterHeads = {KRealAB -> "real_ab", KRealATau -> "real_a_tau", KRealBTau -> "real_b_tau"};
realReduced = Association @ Table[spin -> Collect[Expand[Total[Map[Function[row, With[{n = row["HatPower"]},
  row["Coefficient"] (-2 eps)^row["Dm4Power"] (4 Pi)^n (-1)^n Pochhammer[-eps, n] *
    ((row["Integral"] /. realRules /. d -> 4 + 2 n - 2 eps) /.
      (head : KRealAB | KRealATau | KRealBTau)[indices__] :> KiraMaster[head /. masterHeads, {indices}, 2 n])]],
  realTerms[spin]]]], _KiraMaster, Together], {spin, spins}];
(* Trace coefficients retain their original D; only the reduction rules use d. *)
virtualReduced = Factor[Together[virtualFamily /. D -> 4 - 2 eps /. virtualRules /. d -> 4 - 2 eps]];
Require[FreeQ[realReduced, KRealAB | KRealATau | KRealBTau] && Union[Cases[virtualReduced, _KVirtual, Infinity]] === {KVirtual[0, 1, 1]},
  "Some integrals were not reduced"];
{Union[Cases[realReduced, _KiraMaster, Infinity]], virtualReduced}


(* ::Section:: *)
(*Master integrals*)


(* ::Text:: *)
(*The phase-space masters are the two-body volume Phi2 and the angular integral*)
(**)
(*   Omega(z) = 2 (1-2eps) integral_0^infinity dx x^(-1-eps) (1+x)^eps/(z+x),*)
(**)
(*where z is the cross ratio of the two uncut invariants of each family (zAB, zAT, zBT below).  At the threshold w = 1 the cross ratio zAB vanishes, and Omega(z) behaves like z^(-1-eps): this is where the plus distributions come from, so this behaviour must be kept exact instead of being expanded in eps.  Omega is a Gauss hypergeometric function, and the standard connection formula between the arguments 1-z and z separates the two behaviours at z = 0:*)
(**)
(*   Omega(z) = -2 (1-2eps)/eps [ Gamma(1-eps) Gamma(1+eps) (1-z)^eps z^(-1-eps) + Gamma(1-eps) Gamma(-1-eps)/Gamma(-eps)^2 F1(z) ],   F1(z) = 2F1(1,1;2+eps;z),*)
(**)
(*and in D + 2 dimensions, Omega2(z) = 2 (3-2eps) integral dx x^(-eps) (1+x)^(eps-1)/(z+x),*)
(**)
(*   Omega2(z) = 2 (3-2eps)/(1-eps) [ Gamma(2-eps) Gamma(eps) (1-z)^(eps-1) z^(-eps) + Gamma(2-eps) Gamma(-eps)/Gamma(1-eps)^2 F2(z) ],   F2(z) = 2F1(1,1;1+eps;z).*)
(**)
(*The first terms are exact; F1 and F2 are regular at z = 0 and are evaluated as Euler integrals, F1(z) = (1+eps) integral dy (1+y)^(-1-eps)/(1+(1-z)y) and F2(z) = eps integral dy (1+y)^(-eps)/(1+(1-z)y) over (0, infinity) (the substitution t = y/(1+y) in Euler's integral representation of 2F1).  Phi2 and the loop bubble are Euler integrals as well.*)
(**)
(*EulerIntegral[prefactor, integrand, variables, parameters, order, eps] integrates with SubTropica and returns the expansion in eps through the given order.  The volume Phi2 is needed two orders deeper than the angular functions because the latter come with a 1/eps collinear pole and the threshold expansion adds a second pole.*)


(* ::Input:: *)
LoadSubTropica[]


(* ::Input:: *)
masters = <|
  "bubble" -> EulerIntegral[Gamma[eps], x^(-eps) (1 + x)^(2 eps - 2), {x}, {}, 1, eps],
  "phaseSpace" -> EulerIntegral[1/Gamma[1 - eps], x^(-eps) (1 + x)^(2 eps - 2), {x}, {}, 2, eps],
  "phaseSpaceD2" -> EulerIntegral[1/Gamma[2 - eps], x^(1 - eps) (1 + x)^(2 eps - 4), {x}, {}, 1, eps],
  "F1" -> EulerIntegral[1 + eps, (1 + y)^(-1 - eps)/(1 + zb y), {y}, {zb}, 2, eps] /. zb -> 1 - z,
  "F2" -> EulerIntegral[eps, (1 + y)^(-eps)/(1 + zb y), {y}, {zb}, 2, eps] /. zb -> 1 - z|>


(* ::Text:: *)
(*Checks.  The two-body volumes are Beta functions, Gamma(1-eps)/Gamma(2-2eps) and Gamma(2-eps)/Gamma(4-2eps); the massless bubble is Gamma(eps) Gamma(1-eps)^2/Gamma(2-2eps).  F1 and F2 are compared with Mathematica's Hypergeometric2F1: since the latter cannot be expanded in eps symbolically, we verify that the difference between the exact function and the truncated series scales like eps^3 (the difference at eps/2 must be one eighth of that at eps), which tests every computed order.  Finally the assembled Omega(z) and Omega2(z) are compared with direct numerical integrations in the same way (the truncation error is of order eps^3 for Omega and eps^2 for Omega2, whose regular branch carries a 1/eps prefactor).*)


(* ::Input:: *)
check[x_] := Require[Simplify[x] === 0, "Master integral check failed", x];
check[Normal[Series[masters["phaseSpace"] - Gamma[1 - eps]/Gamma[2 - 2 eps], {eps, 0, 2}]]];
check[Normal[Series[masters["phaseSpaceD2"] - Gamma[2 - eps]/Gamma[4 - 2 eps], {eps, 0, 1}]]];
check[Normal[Series[masters["bubble"] - Gamma[eps] Gamma[1 - eps]^2/Gamma[2 - 2 eps], {eps, 0, 1}]]];
truncationCheck[series_, exact_, power_] := Module[{d1, d2},
  d1 = N[(exact - series) /. {eps -> -1/50, z -> 1/3}, 30];
  d2 = N[(exact - series) /. {eps -> -1/100, z -> 1/3}, 30];
  Require[Abs[d2/d1 - 2^-power] < 2^-power/10, "Master integral disagrees with its closed form", {d1, d2}]; {d1, d2}];
truncationCheck[masters["F1"], Hypergeometric2F1[1, 1, 2 + eps, z], 3];
truncationCheck[masters["F2"], Hypergeometric2F1[1, 1, 1 + eps, z], 3]


(* ::Input:: *)
(* F1 and F2 are kept as inert symbols F1[n][z] and F2[n][z], the coefficients of
   eps^n, until the end of step 4: the expressions stay small, and at the
   threshold their values are exact, F1[0][0] = F2[0][0] = 1 and F1[n][0] =
   F2[n][0] = 0 for n > 0, because 2F1(1,1;c;0) = 1 for every c. *)
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
(*Normalization of the phase-space integrals: with measure d^D r/(2 Pi)^(D-1) delta+(r^2) 2 Pi delta+((pRec - r)^2), the volume is*)
(**)
(*   Phi2 = (4 Pi)^eps s4^(-eps)/(8 Pi) Gamma(1-eps)/Gamma(2-2eps),*)
(**)
(*and in D + 2 dimensions Phi2' = (4 Pi)^(eps-1) s4^(1-eps)/(8 Pi) Gamma(2-eps)/Gamma(4-2eps).  The angular masters are Phi2 times Omega(z) divided by the product of the two invariants at the symmetric point, with the cross ratios*)
(**)
(*   zAB = -t s4/((s+t)(s-s4)),   zAT = -t (s+t-s4)/((s+t)(s4-t)),   zBT = -s t/((s-s4)(s4-t)),*)
(**)
(*all between 0 and 1 in the physical region.  The loop master is the bubble with measure d^D loop/(i Pi^(D/2)): (-t)^(-eps) Gamma(eps) Gamma(1-eps)^2/Gamma(2-2eps).*)


(* ::Input:: *)
phaseSpace = (4 Pi)^eps s4^(-eps)/(8 Pi) masters["phaseSpace"];
phaseSpaceD2 = (4 Pi)^(eps - 1) s4^(1 - eps)/(8 Pi) masters["phaseSpaceD2"];
zAB = -t s4/((s + t) (s - s4));  zAT = -t (s + t - s4)/((s + t) (s4 - t));  zBT = -s t/((s - s4) (s4 - t));
physicalMasters = {
  KiraMaster["real_ab", {1, 1, 0, 0}, 0] -> phaseSpace,
  KiraMaster["real_ab", {1, 1, 1, 1}, 0] -> -phaseSpace angular[zAB]/((s + t) (s - s4)),
  KiraMaster["real_a_tau", {1, 1, 1, 1}, 0] -> phaseSpace angular[zAT]/((s + t) (s4 - t)),
  KiraMaster["real_b_tau", {1, 1, 1, 1}, 0] -> -phaseSpace angular[zBT]/((s - s4) (s4 - t)),
  KiraMaster["real_ab", {1, 1, 0, 0}, 2] -> phaseSpaceD2,
  KiraMaster["real_ab", {1, 1, 1, 1}, 2] -> -phaseSpaceD2 angularD2[zAB]/((s + t) (s - s4)),
  KiraMaster["real_a_tau", {1, 1, 1, 1}, 2] -> phaseSpaceD2 angularD2[zAT]/((s + t) (s4 - t)),
  KiraMaster["real_b_tau", {1, 1, 1, 1}, 2] -> -phaseSpaceD2 angularD2[zBT]/((s - s4) (s4 - t))};
realIntegrated = Association @ Table[spin -> CF ee^4 eq^2 gs^2 (realReduced[spin] /. physicalMasters), {spin, spins}];
virtualIntegrated = virtualReduced /. KVirtual[0, 1, 1] -> (-t)^(-eps) masters["bubble"];
Require[FreeQ[{realIntegrated, virtualIntegrated}, KiraMaster | KVirtual], "A master integral was not inserted"];
Series[virtualIntegrated, {eps, 0, 0}]


(* ::Text:: *)
(*The virtual form factor (times the Born square) has the familiar poles -2/eps^2 - 3/eps.  The real-emission result is kept with its exact threshold factors s4^(-eps), zAB^(-1-eps) and zAB^(-eps) and with the inert symbols F1[n][z], F2[n][z]; the replacement rules masterFunctions, saved with the result, give their explicit forms.*)


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
reduced = <|"Real" -> realIntegrated, "Virtual" -> virtualIntegrated, "Masters" -> masters,
  "MasterFunctions" -> masterFunctions, "RealRules" -> realRules, "VirtualRules" -> virtualRules, "Route" -> route|>;
SaveResult["Reduced", reduced, route];
Keys[reduced]
