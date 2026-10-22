(* ::Title:: *)
(*Proton-proton classes, step 3: scalar integrals*)


(* ::Text:: *)
(*Input: Results/<route>/class-<n>/Algebra.wl (step 2).  Output: Results/<route>/class-<n>/Reduced.wl with the real-emission phase-space integrals and the one-loop integrals evaluated as Laurent series in eps, exact in the threshold variable.  The method is that of the lepton-quark process, extended to two unobserved partons that both carry propagator momenta and to one-loop boxes.*)
(**)
(*Real emission.  The two unobserved partons k1, k2 (k1 + k2 = pRec = pa + pb - p, pRec^2 = s4 = s + t + u) are integrated with the two-body phase space written as cut propagators k1^2 and k2^2.  Every denominator of the squared matrix element is a linear form in the two independent invariants A = (k1-pa)^2 and B = (k1-pb)^2: besides A and B these are C = (k1+p)^2 = -s4 - A - B and, with k2, A' = (k2-pa)^2 = -(s+t) - A, B' = (k2-pb)^2 = -(s+u) - B, C' = (k2+p)^2 = s + A + B.  Any three of them are linearly dependent (two variables), so partial fractions reduce every term to at most two forms, and the exchange k1 <-> k2 (a symmetry of the phase space) maps the primed pairs onto unprimed ones.  Six families remain, defined in kira/real/config/integralfamilies.yaml:*)
(**)
(*   real_AB[1,1,n3,n4] = integral of 1/(A^n3 B^n4),  real_AC, real_BC, real_ABp (A and B'), real_ACp (A and C'), real_BCp (B and C').*)
(**)
(*Kira reduces them to the two-body volume and one angular master per family.  The angular master of two massless "legs" q, q' (the external momenta of the two forms) depends on the angle between them in the rest frame of pRec; with z = s4 (q.q')/(2 (q.pRec)(q'.pRec)) for two forms with the same unobserved momentum, and 1 - z when one form has k1 and the other k2,*)
(**)
(*   integral of 1/(X Y) over the phase space = sigmaX sigmaY Phi2 Omega(z)/((2 q.pRec)(2 q'.pRec)),   sigma = -1 for A, B, A', B' (= -2 q.k) and +1 for C, C' (= +2 q.k),*)
(**)
(*with Omega(z) = 2(1-2eps) integral_0^Infinity dx x^(-1-eps)(1+x)^eps/(z+x) exactly as in the lepton-quark process, whose threshold branches z^(-1-eps) and regular part F1(z) are separated by the hypergeometric connection formula (that notebook explains it); the hatted invariant spe3 shifts the dimension by two for each power, with the masters Phi2' and Omega2.*)
(**)
(*Virtual.  The one-loop families virtual_ST, virtual_SU, virtual_TU (massless boxes) are reduced by Kira to the bubbles in s, t, u and the boxes in (s,t), (s,u), (t,u).  The bubbles are Euler integrals; the boxes are evaluated from their Feynman-parameter representation with SubTropica in the Euclidean region and continued to the physical region (s > 0), where their real parts enter.  Powers of the hatted loop momentum loopHat2 shift the dimension by two (masters in D+2).*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
class = "13";
algebra = LoadResult["Algebra", route <> "/class-" <> class];
spins = {"UU", "LL"};
Keys[algebra]


(* ::Section:: *)
(*Real emission: partial fractions*)


(* ::Text:: *)
(*The squared matrix element, with D = 4 - 2 eps, is a rational function of A and B (C is eliminated) with coefficients rational in eps.  Its denominator factors are the six linear forms, which are replaced by the inert symbols fA, fB, fC, fAp, fBp, fCp; then two kinds of rules are applied repeatedly until no term contains more than two forms: for a dependent pair (A A', B B', C C', whose sums are constants) 1/(L L') = (1/L + 1/L')/(L + L'), and for any three forms with alpha L1 + beta L2 + gamma L3 = delta (constants), 1/(L1 L2 L3) = [alpha/(L2 L3) + beta/(L1 L3) + gamma/(L1 L2)]/delta.*)


(* ::Input:: *)
forms = <|fA -> A, fB -> B, fC -> -s4 - A - B, fAp -> -(s + t) - A, fBp -> -(s + u) - B, fCp -> s + A + B|> /. u -> s4 - s - t;
formOf[factor_] := Module[{f = Expand[factor], matches},
  matches = Select[Keys[forms], Expand[f - forms[#]] === 0 || Expand[f + forms[#]] === 0 &];
  Require[Length[matches] == 1, "A denominator factor is not one of the six linear forms", factor];
  First[matches] If[Expand[f - forms[First[matches]]] === 0, 1, -1]];
(* write the denominator through the inert form symbols; the numerator keeps A, B (they are removed below) *)
toForms[expr_] := Module[{mapped = MapIntegrationDenominators[expr,
    {u -> s4 - s - t, C -> -s4 - A - B}, {A, B}, formOf]},
  Require[!FailureQ[mapped], "Termwise denominator conversion failed", mapped];
  mapped];
(* numerators: A and B are themselves forms *)
numeratorForms = {A -> fA, B -> fB};
relation[names_List] := Module[{L = forms /@ names, al, be, ga, sol},
  sol = Quiet[Solve[{Coefficient[al L[[1]] + be L[[2]] + ga L[[3]], A] == 0, Coefficient[al L[[1]] + be L[[2]] + ga L[[3]], B] == 0, al == 1}, {al, be, ga}], Solve::svars];
  (* no relation with al = 1, or a one-parameter family of relations (two of the forms are a dependent pair, handled by pairRules): no triple rule *)
  If[sol === {} || !FreeQ[{al, be, ga} /. First[sol], al | be | ga], None, With[{c = {al, be, ga} /. First[sol]}, {c, Expand[c.L /. {A -> 0, B -> 0}]}]]];
pairRules = Table[With[{f1 = pair[[1]], f2 = pair[[2]], c = Expand[forms[pair[[1]]] + forms[pair[[2]]]]},
  f1^a_. f2^b_. /; a < 0 && b < 0 :> f1^(a + 1) f2^(b + 1) (1/f1 + 1/f2)/c], {pair, {{fA, fAp}, {fB, fBp}, {fC, fCp}}}];
tripleRules = Flatten[Table[With[{r = relation[triple]}, If[r === None || r[[2]] === 0, Nothing,
  With[{f1 = triple[[1]], f2 = triple[[2]], f3 = triple[[3]], c = r[[1]], de = r[[2]]},
    f1^a_. f2^b_. f3^d_. /; a < 0 && b < 0 && d < 0 :> f1^(a + 1) f2^(b + 1) f3^(d + 1) (c[[1]]/(f2 f3) + c[[2]]/(f1 f3) + c[[3]]/(f1 f2))/de]]],
  {triple, Subsets[Keys[forms], {3}]}]];
negativeForms[term_] := Select[Keys[forms], Exponent[term, #] < 0 &];
(* the additive terms of an expression: a single product is one term, not a list of its factors; zero has none *)
termList[expr_] := With[{e = Expand[expr]}, Which[e === 0, {}, Head[e] === Plus, List @@ e, True, {e}]];
(* Partial fractions of a unit propagator monomial. Scalar coefficients are
   spectators and are combined by the shared sparse converter. *)
partialFractionMonomial[expr_] := Module[{e = expr, old},
  (* a numerator power of A (B) over a denominator A' (B') is written through A' (B'): the pair is dependent *)
  e = Expand[e /. {fA^n_. fAp^m_. /; n > 0 && m < 0 :> (-(s + t) - fAp)^n fAp^m, fB^n_. fBp^m_. /; n > 0 && m < 0 :> (-(s + u) - fBp)^n fBp^m} /. u -> s4 - s - t];
  Do[old = e; e = Expand[e //. pairRules]; e = Expand[e /. tripleRules]; If[e === old, Break[]], {20}];
  Require[Max[Prepend[Length[negativeForms[#]] & /@ termList[e], 0]] <= 2, "Partial fractions did not terminate"];
  Require[Together[(e - expr) /. forms] === 0,
    "A geometric partial fraction does not reconstruct its monomial", expr];
  e];
fractionConverter = CreatePartialFractionConverter[Keys[forms], spe3, partialFractionMonomial];
realFractionResults = Association @ Table[spin -> fractionConverter[
  toForms[algebra["Real"][spin] /. D -> 4 - 2 eps] /. numeratorForms], {spin, spins}];
Require[FreeQ[realFractionResults, _Failure], "Sparse partial fractions failed", realFractionResults];
realFractionTerms = Association @ Table[spin -> realFractionResults[spin]["Terms"], {spin, spins}];
realFractions = Total /@ realFractionTerms;   (* the coefficients are rational in eps (the gluon average 1/(D-2)) *)
Length /@ realFractions


(* ::Text:: *)
(*Check: the partial fractions reproduce the matrix element when the forms are substituted back.*)


(* ::Input:: *)
(* exact rational test points: a symbolic simplification of the difference of tens of thousands of terms would take hours *)
testPoints = {{s -> 7, t -> -2, s4 -> 3/2, A -> -5/3, B -> -4/7, eps -> 1/5, spe3 -> 2/9}, {s -> 11, t -> -6, s4 -> 2/3, A -> -1/4, B -> -9/5, eps -> -3/7, spe3 -> 5/13}};
Do[Require[Together[((realFractions[spin] /. forms) - (algebra["Real"][spin] /. D -> 4 - 2 eps /. u -> s4 - s - t /. C -> -s4 - A - B)) /. point] === 0,
  "Partial fractions do not add up", {spin, point}], {spin, spins}, {point, testPoints}];
"partial fractions checked at exact test points"


(* ::Section:: *)
(*Families and Kira targets*)


(* ::Text:: *)
(*Each term, a coefficient times a product of at most two forms, is mapped onto a family and its indices.  The exchange k1 <-> k2 turns A', B', C' into A, B, C; pairs with one primed form are mapped onto the mixed families real_ABp (A with B'), real_ACp (A with C') and real_BCp (B with C').  The coefficients, functions of s, t, s4 only, are not affected.  The symbols dm4 = D - 4 and the hatted invariant spe3 are separated: a power spe3^n means the integral in D + 2n dimensions times (4 Pi)^n (-1)^n (-eps)_n (the lepton-quark notebook derives this).*)


(* ::Input:: *)
(* the families and their two forms; the exchange k1 <-> k2 maps every form onto its primed partner *)
familyForms = <|KRealAB -> {fA, fB}, KRealAC -> {fA, fC}, KRealBC -> {fB, fC}, KRealABp -> {fA, fBp}, KRealACp -> {fA, fCp}, KRealBCp -> {fB, fCp}|>;
swapForms = {fA -> fAp, fAp -> fA, fB -> fBp, fBp -> fB, fC -> fCp, fCp -> fC};
swapInvariants = {A -> -(s + t) - A, B -> -(s + u) - B} /. u -> s4 - s - t;         (* k1 <-> k2 on the invariants *)
headOf[names_List] := Which[names === {} || names === {fA} || names === {fB} || Sort[names] === {fA, fB}, KRealAB,
  names === {fC} || Sort[names] === Sort[{fA, fC}], KRealAC,   Sort[names] === Sort[{fB, fC}], KRealBC,
  Sort[names] === Sort[{fA, fBp}], KRealABp,   Sort[names] === Sort[{fA, fCp}], KRealACp,   Sort[names] === Sort[{fB, fCp}], KRealBCp,
  True, None];
(* family of a term from its denominator forms: with or without the relabelling k1 <-> k2 *)
familyOf[denominators_List] := Module[{head = headOf[denominators], swap = False},
  If[head === None, swap = True; head = headOf[denominators /. swapForms]];
  Require[head =!= None, "Unexpected pair of denominator forms", denominators];
  {head, swap}];
(* A and B through the two forms X, Y of each family, solved once *)
familySolutions = Association @ KeyValueMap[Function[{head, pair}, head -> First[Solve[{X == forms[pair[[1]]], Y == forms[pair[[2]]]}, {A, B}]]], familyForms];
termRows[expr_] := Module[{terms = CombineLaurentTerms[expr, Append[Keys[forms], spe3]]},
  Require[ListQ[terms] && FreeQ[terms, _Failure], "Laurent coefficient grouping failed", terms];
  Flatten[Table[Module[{den, head, swap, pair, powers, numerator, sol, poly, rows},
    den = negativeForms[term];
    {head, swap} = familyOf[den];
    pair = familyForms[head]; sol = familySolutions[head];
    (* denominator powers of the two family forms (after the relabelling, if any) *)
    powers = Association @ Table[(f /. If[swap, swapForms, {}]) -> -Exponent[term, f], {f, den}];
    (* the rest of the term: coefficient times positive powers of forms, written in A and B (relabelled with the term) and then in the two family forms X, Y *)
    numerator = term Times @@ Table[f^(-Exponent[term, f]), {f, den}] /. forms;
    If[swap, numerator = numerator /. swapInvariants];
    poly = Expand[numerator /. sol];
    Require[PolynomialQ[poly, {X, Y, spe3}], "A numerator is not polynomial in the family forms", term];
    rows = Table[<|"Integral" -> head[1, 1, Lookup[powers, pair[[1]], 0] - rule[[1, 1]], Lookup[powers, pair[[2]], 0] - rule[[1, 2]]],
      "Dm4Power" -> 0, "HatPower" -> rule[[1, 3]], "Coefficient" -> Factor[rule[[2]]]|>, {rule, CoefficientRules[poly, {X, Y, spe3}]}];
    (* check: the rows rebuild the (relabelled) term *)
    Require[Together[(Total[(#["Coefficient"] dm4^#["Dm4Power"] spe3^#["HatPower"] X^(-#["Integral"][[3]]) Y^(-#["Integral"][[4]]) & /@ rows)] /.
        {X -> forms[pair[[1]]], Y -> forms[pair[[2]]]}) - (term /. forms /. If[swap, swapInvariants, {}])] === 0, "The rows of a term do not rebuild it", term];
    rows], {term, terms}]]];
(* The existing mapper above acts on a unit monomial and retains its exact
   reconstruction test. After multiplication by the scalar coefficient, check
   the full contribution again, including any k1 <-> k2 relabelling. *)
reconstructRows[term_, rows_List] := Module[{head, swap, pair, rebuilt},
  {head, swap} = familyOf[negativeForms[term]]; pair = familyForms[head];
  rebuilt = Total[(#["Coefficient"] dm4^#["Dm4Power"] spe3^#["HatPower"]
    forms[pair[[1]]]^(-#["Integral"][[3]]) forms[pair[[2]]]^(-#["Integral"][[4]])) & /@ rows];
  Together[rebuilt - (term /. forms /. If[swap, swapInvariants, {}])] === 0];
rowConverter = CreateIntegralRowConverter[Append[Keys[forms], spe3], termRows, reconstructRows];
realRowResults = Association @ Table[spin -> rowConverter[realFractionTerms[spin]], {spin, spins}];
Require[FreeQ[realRowResults, _Failure], "Sparse integral-row conversion failed", realRowResults];
realTerms = Association @ Table[spin -> realRowResults[spin]["Rows"], {spin, spins}];
indexRanges = <|"denominators" -> Max[Total[Select[List @@ #["Integral"], # > 0 &]] & /@ Flatten[Values[realTerms]]],
  "numerators" -> Max[-Total[Select[List @@ #["Integral"], # < 0 &]] & /@ Flatten[Values[realTerms]]]|>;
Require[indexRanges["denominators"] <= 8 && indexRanges["numerators"] <= 8, "The integrals exceed the ranges r = 8, s = 8 of kira/real/jobs.yaml", indexRanges];
indexRanges


(* ::Input:: *)
Dataset[realTerms["UU"][[1 ;; Min[8, Length[realTerms["UU"]]]]]]


(* ::Section:: *)
(*Reduction with Kira*)


(* ::Text:: *)
(*Real families real_AB ... real_BCp (kira/real) and the one-loop families (kira/virtual); Kira sets s = 1, RestoreScale puts it back.  The virtual integrals come from step 2 already labelled (KVirtualST etc.).*)


(* ::Input:: *)
realHeads = {KRealAB -> "real_AB", KRealAC -> "real_AC", KRealBC -> "real_BC", KRealABp -> "real_ABp", KRealACp -> "real_ACp", KRealBCp -> "real_BCp"};
virtualHeads = {KVirtualST -> "virtual_ST", KVirtualSU -> "virtual_SU", KVirtualTU -> "virtual_TU"};
referenceRealMasters = Join[{KRealAB[1, 1, 0, 0]}, #[1, 1, 1, 1] & /@ Keys[realHeads]];
referenceVirtualMasters = Flatten[Table[head @@ indices, {head, Keys[virtualHeads]},
  {indices, Select[Tuples[{0, 1}, 4], MemberQ[{2, 4}, Total[#]] &]}]];
realTargets = Union[Flatten[Lookup[realTerms[#], "Integral"] & /@ spins], referenceRealMasters];
virtualTargets = Union[Cases[algebra["Virtual"], (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]];
resultsDir = ResultsDirectory[route <> "/class-" <> class];
kiraReal = PrepareKira[FileNameJoin[{NotebookDirectory[], "kira", "real"}], FileNameJoin[{resultsDir, "kira", "real"}]];
WriteKiraTargets[kiraReal, realTargets, realHeads];
If[virtualTargets =!= {},
  kiraVirtual = PrepareKira[FileNameJoin[{NotebookDirectory[], "kira", "virtual"}], FileNameJoin[{resultsDir, "kira", "virtual"}]];
  WriteKiraTargets[kiraVirtual, virtualTargets, virtualHeads]];
{Length[realTargets], Length[virtualTargets]}


(* ::Input:: *)
RunKira[kiraReal];
realRules = RestoreScale[Join @@ Table[KiraRules[kiraReal, name, realHeads], {name, Select[Values[realHeads], FileExistsQ[FileNameJoin[{kiraReal, "results", #, "kira_targets.m"}]] &]}],
  s, {t -> t/s, s4 -> s4/s}, Keys[realHeads]];
virtualRules = If[virtualTargets === {}, {}, RunKira[kiraVirtual];
  RestoreScale[Join @@ Table[KiraRules[kiraVirtual, name, virtualHeads], {name, Select[Values[virtualHeads], FileExistsQ[FileNameJoin[{kiraVirtual, "results", #, "kira_targets.m"}]] &]}],
    s, {t -> t/s}, Keys[virtualHeads]]];
realMasters = Union[Cases[realRules[[All, 2]], (Alternatives @@ Keys[realHeads])[__], Infinity]];
virtualMasters = Union[Cases[virtualRules[[All, 2]], (Alternatives @@ Keys[virtualHeads])[__], Infinity]];
CheckReductionCoverage[realTargets, realRules, Keys[realHeads], referenceRealMasters];
CheckReductionCoverage[virtualTargets, virtualRules, Keys[virtualHeads], referenceVirtualMasters];
{realMasters, virtualMasters}


(* ::Text:: *)
(*Expected masters: the volume real_AB[1,1,0,0] and the angular integrals [1,1,1,1] of each family; for the loop: the bubbles [0,1,0,1] (s channel), [0,1,1,0] / [1,0,0,1] type (t, u) and the boxes [1,1,1,1].  A term with spe3^n or loopHat2^n uses the rules at d -> 4 + 2n - 2 eps; KiraMaster[family, indices, 2n] denotes a master in D + 2n dimensions.*)


(* ::Input:: *)
label[head_, heads_] := head /. heads;
(* the rows are first grouped by integral and hat power (their coefficients added), so that the reduction rules are applied once per distinct integral *)
realGrouped = Association @ Table[spin -> GroupBy[realTerms[spin], {#["Integral"], #["HatPower"]} &, Total[#["Coefficient"] (-2 eps)^#["Dm4Power"] & /@ #] &], {spin, spins}];
realDispatch = Dispatch[realRules];
realReduced = Association @ Table[spin -> Collect[Total[KeyValueMap[Function[{key, coefficient}, With[{n = key[[2]]},
  coefficient (4 Pi)^n (-1)^n Pochhammer[-eps, n] *
    ((key[[1]] /. realDispatch /. d -> 4 + 2 n - 2 eps) /.
      (head : (Alternatives @@ Keys[realHeads]))[indices__] :> KiraMaster[label[head, realHeads], {indices}, 2 n])]],
  realGrouped[spin]]], _KiraMaster, Together], {spin, spins}];
(* the trace was evaluated in D = 4 - 2 eps dimensions, and its coefficients keep that D: only the integration dimension of a hatted term is shifted
   (the Kira rules are written in their own dimension symbol d, set below to 4 + 2n - 2 eps).  Replacing the trace's D by d would shift the coefficients as well and is wrong.
   A power loopHat2^n of the hatted loop momentum under the loop integral d^D loop/(i Pi^(D/2)) equals (-1)^n (-eps)_n times the integral in D + 2n dimensions
   (the (D-4)-dimensional angular average; the measure d^D k/(2 Pi)^(D-1) of the real emission gives the additional (4 Pi)^n), see the introduction *)
virtualReduced = Association @ Table[spin -> Collect[Expand[Total[Map[Function[term, With[{n = Exponent[term, HatPower[_]] /. 0 -> 0},
  Module[{hat = FirstCase[term, HatPower[m_] :> m, 0]},
    (-1)^hat Pochhammer[-eps, hat] (term /. HatPower[_] -> 1 /. D -> 4 - 2 eps /. virtualRules /. d -> 4 + 2 hat - 2 eps) /.
      (head : (Alternatives @@ Keys[virtualHeads]))[indices__] :> KiraMaster[label[head, virtualHeads], {indices}, 2 hat]]]],
  If[algebra["Virtual"][spin] === 0, {}, termList[algebra["Virtual"][spin]]]]]], _KiraMaster, Together], {spin, spins}];
Require[FreeQ[{realReduced, virtualReduced}, Alternatives @@ Join[Keys[realHeads], Keys[virtualHeads]]], "Some integrals were not reduced"];
(* the box masters (four propagators) are evaluated through order eps^0 only: their coefficients must be regular at eps = 0 *)
boxMasters = Union[Cases[virtualReduced, m_KiraMaster /; StringStartsQ[m[[1]], "virtual"] && Count[m[[2]], 1] == 4, Infinity]];
Require[And @@ Flatten[Table[FreeQ[Normal[Series[Coefficient[virtualReduced[spin], m], {eps, 0, -1}]], eps], {spin, spins}, {m, boxMasters}]],
  "A box coefficient has a pole at eps = 0: the box masters would be needed beyond order eps^0"];
Union[Cases[{realReduced, virtualReduced}, _KiraMaster, Infinity]]


(* ::Section:: *)
(*Phase-space masters*)


(* ::Text:: *)
(*As in the lepton-quark process: the two-body volume Phi2 (and Phi2' in D+2), the angular function Omega(z) with its threshold branch z^(-1-eps) separated by the connection formula (F1, F2 inert), and for the mixed families the regular function Omega(1-z) = 2(1-2eps) integral_0^Infinity dy (1+y)^eps/(1+(1-z)y) (its D+2 analogue with (1+y)^(eps-1) and 2(3-2eps)), which is analytic at z = 0 and kept as the inert functions F3, F4.*)


(* ::Input:: *)
LoadSubTropica[]


(* ::Input:: *)
masters = <|
  "phaseSpace" -> EulerIntegral[1/Gamma[1 - eps], x^(-eps) (1 + x)^(2 eps - 2), {x}, {}, 2, eps],
  "phaseSpaceD2" -> EulerIntegral[1/Gamma[2 - eps], x^(1 - eps) (1 + x)^(2 eps - 4), {x}, {}, 1, eps],
  "F1" -> EulerIntegral[1 + eps, (1 + y)^(-1 - eps)/(1 + zb y), {y}, {zb}, 2, eps] /. zb -> 1 - z,
  "F2" -> EulerIntegral[eps, (1 + y)^(-eps)/(1 + zb y), {y}, {zb}, 2, eps] /. zb -> 1 - z,
  "F3" -> EulerIntegral[2 (1 - 2 eps), (1 + y)^eps/(1 + zb y), {y}, {zb}, 2, eps] /. zb -> 1 - z,
  "F4" -> EulerIntegral[2 (3 - 2 eps), (1 + y)^(eps - 1)/(1 + zb y), {y}, {zb}, 1, eps] /. zb -> 1 - z|>;
check[x_] := Require[Simplify[x] === 0, "Master integral check failed", x];
check[Normal[Series[masters["phaseSpace"] - Gamma[1 - eps]/Gamma[2 - 2 eps], {eps, 0, 2}]]];
check[Normal[Series[masters["phaseSpaceD2"] - Gamma[2 - eps]/Gamma[4 - 2 eps], {eps, 0, 1}]]];
truncationCheck[series_, exact_, power_] := Module[{d1, d2},
  d1 = N[(exact - series) /. {eps -> -1/50, z -> 1/3}, 30];
  d2 = N[(exact - series) /. {eps -> -1/100, z -> 1/3}, 30];
  Require[Abs[d2/d1 - 2^-power] < 2^-power/10, "Master integral disagrees with its closed form", {d1, d2}]; {d1, d2}];
truncationCheck[masters["F1"], Hypergeometric2F1[1, 1, 2 + eps, z], 3];
truncationCheck[masters["F2"], Hypergeometric2F1[1, 1, 1 + eps, z], 3];
truncationCheck[masters["F3"], -2 (1 - 2 eps)/eps Hypergeometric2F1[1, 1, 1 - eps, z], 3];
truncationCheck[masters["F4"], 2 (3 - 2 eps)/(1 - eps) Hypergeometric2F1[1, 1, 2 - eps, z], 2]


(* ::Input:: *)
angular[zz_] := -2 (1 - 2 eps)/eps (Gamma[1 - eps] Gamma[1 + eps] (1 - zz)^eps zz^(-1 - eps) +
  Gamma[1 - eps] Gamma[-1 - eps]/Gamma[-eps]^2 Sum[eps^n F1[n][zz], {n, 0, 2}]);
angularD2[zz_] := 2 (3 - 2 eps)/(1 - eps) (Gamma[2 - eps] Gamma[eps] (1 - zz)^(eps - 1) zz^(-eps) +
  Gamma[2 - eps] Gamma[-eps]/Gamma[1 - eps]^2 Sum[eps^n F2[n][zz], {n, 0, 2}]);
angularMixed[zz_] := Sum[eps^n F3[n][zz], {n, -1, 2}];                 (* Omega(1-z): regular at z = 0 *)
angularMixedD2[zz_] := Sum[eps^n F4[n][zz], {n, 0, 1}];
(* at the threshold the mixed masters are elementary: Omega(1 - z) at z = 0 is 2 (1-2eps) B(-eps, 1) = -2/eps + 4, and the D+2 one 2 (3-2eps) B(1-eps, 1) = 2 (3-2eps)/(1-eps);
   the explicit functions are regular at z = 0 and must reproduce these values (checked below) *)
thresholdValues = Join[Table[F3[n][0] -> SeriesCoefficient[-2 (1 - 2 eps)/eps, {eps, 0, n}], {n, -1, 2}], Table[F4[n][0] -> SeriesCoefficient[2 (3 - 2 eps)/(1 - eps), {eps, 0, n}], {n, 0, 1}]];
thresholdValues = Join[{F1[0][0] -> 1, F2[0][0] -> 1, F1[_][0] -> 0, F2[_][0] -> 0}, thresholdValues];   (* 2F1(1,1;c;0) = 1 for the unmixed masters *)
masterFunctions = Join[thresholdValues,
  Table[F1[n][zz_] :> Evaluate[SeriesCoefficient[masters["F1"], {eps, 0, n}] /. z -> zz], {n, 0, 2}],
  Table[F2[n][zz_] :> Evaluate[SeriesCoefficient[masters["F2"], {eps, 0, n}] /. z -> zz], {n, 0, 2}],
  Table[F3[n][zz_] :> Evaluate[SeriesCoefficient[masters["F3"], {eps, 0, n}] /. z -> zz], {n, -1, 2}],
  Table[F4[n][zz_] :> Evaluate[SeriesCoefficient[masters["F4"], {eps, 0, n}] /. z -> zz], {n, 0, 1}]];
numericCheck[expr_, integrand_, power_] := Module[{d1, d2},
  d1 = N[expr /. {eps -> -1/40, z -> 1/3}, 20] - NIntegrate[integrand /. {eps -> -1/40, z -> 1/3}, {x, 0, Infinity}, WorkingPrecision -> 30];
  d2 = N[expr /. {eps -> -1/80, z -> 1/3}, 20] - NIntegrate[integrand /. {eps -> -1/80, z -> 1/3}, {x, 0, Infinity}, WorkingPrecision -> 30];
  Require[Abs[d2/d1 - 2^-power] < 2^-power/10, "Angular master disagrees with direct integration", {d1, d2}]; {d1, d2}];
Require[And @@ Table[Simplify[Limit[F3[n][z] /. masterFunctions, z -> 0] - (F3[n][0] /. thresholdValues)] === 0, {n, -1, 2}] &&
  And @@ Table[Simplify[Limit[F4[n][z] /. masterFunctions, z -> 0] - (F4[n][0] /. thresholdValues)] === 0, {n, 0, 1}], "The mixed masters do not reach their threshold values"];
numericCheck[angular[z] /. masterFunctions, 2 (1 - 2 eps) x^(-1 - eps) (1 + x)^eps/(z + x), 3];
numericCheck[angularD2[z] /. masterFunctions, 2 (3 - 2 eps) x^(-eps) (1 + x)^(eps - 1)/(z + x), 2];
numericCheck[angularMixed[z] /. masterFunctions, 2 (1 - 2 eps) x^(-1 - eps) (1 + x)^eps/(1 - z + x), 3];
numericCheck[angularMixedD2[z] /. masterFunctions, 2 (3 - 2 eps) x^(-eps) (1 + x)^(eps - 1)/(1 - z + x), 2]


(* ::Text:: *)
(*The physical masters.  Phi2 = (4 Pi)^eps s4^(-eps)/(8 Pi) Gamma(1-eps)/Gamma(2-2eps), Phi2' = (4 Pi)^(eps-1) s4^(1-eps)/(8 Pi) Gamma(2-eps)/Gamma(4-2eps); the angular masters with the cross ratios of the pairs (pa, pb), (pa, p), (pb, p),*)
(**)
(*   zAB = s s4/((s+t)(s+u)),   zAC = -t s4/((s+t)(s-s4)),   zBC = -u s4/((s+u)(s-s4)),   (2 pa.pRec = s+t, 2 pb.pRec = s+u, 2 p.pRec = s-s4),*)
(**)
(*and the normalization sigma sigma' Phi2 Omega(z)/((2 q.pRec)(2 q'.pRec)) explained in the introduction.  The lepton-quark process checked this formula for the unmixed pairs against its literature result; the mixed one is checked here by a direct numerical angular integration in D = 4 - 2 eps dimensions at eps = -1/3 (the angular measure of the recoil rest frame is sin^(1-2eps)theta1 dtheta1 sin^(-2eps)theta2 dtheta2), against the exact Euler integral of Omega(1 - zAB).*)


(* ::Input:: *)
twoPRec = <|pa -> s + t, pb -> s + u, p -> s - s4|>;
sign = <|pa -> -1, pb -> -1, p -> 1|>;
crossRatio[q_, qp_] := s4 (2 ScalarProductValue[q, qp])/(twoPRec[q] twoPRec[qp]);        (* z = s4 (q.q')/(2 (q.pRec)(q'.pRec)) with twoPRec = 2 q.pRec *)
ScalarProductValue[pa, pb] = s/2; ScalarProductValue[pa, p] = -t/2; ScalarProductValue[pb, p] = -u/2;
zAB = Simplify[crossRatio[pa, pb]]; zAC = Simplify[crossRatio[pa, p]]; zBC = Simplify[crossRatio[pb, p]];
Require[Simplify[{zAB - s s4/((s + t) (s + u)), zAC + t s4/((s + t) (s - s4)), zBC + u s4/((s + u) (s - s4))}] === {0, 0, 0}, "cross ratios"];
phaseSpace = (4 Pi)^eps s4^(-eps)/(8 Pi) masters["phaseSpace"];
phaseSpaceD2 = (4 Pi)^(eps - 1) s4^(1 - eps)/(8 Pi) masters["phaseSpaceD2"];
pairs = <|"real_AB" -> {pa, pb, zAB, False}, "real_AC" -> {pa, p, zAC, False}, "real_BC" -> {pb, p, zBC, False},
  "real_ABp" -> {pa, pb, zAB, True}, "real_ACp" -> {pa, p, zAC, True}, "real_BCp" -> {pb, p, zBC, True}|>;
physicalMasters = Join[{KiraMaster["real_AB", {1, 1, 0, 0}, 0] -> phaseSpace, KiraMaster["real_AB", {1, 1, 0, 0}, 2] -> phaseSpaceD2},
  Flatten[Table[With[{q = pairs[name][[1]], qp = pairs[name][[2]], zz = pairs[name][[3]], mixed = pairs[name][[4]]},
    {KiraMaster[name, {1, 1, 1, 1}, 0] -> sign[q] sign[qp] phaseSpace If[mixed, angularMixed[zz], angular[zz]]/(twoPRec[q] twoPRec[qp]),
     KiraMaster[name, {1, 1, 1, 1}, 2] -> sign[q] sign[qp] phaseSpaceD2 If[mixed, angularMixedD2[zz], angularD2[zz]]/(twoPRec[q] twoPRec[qp])}],
    {name, Keys[pairs]}]]];
Length[physicalMasters]


(* ::Input:: *)
(* numerical check of the mixed-family master formula: integral over the two-body phase space of 1/(A B') at eps = -1/3, s = 1, t = -3/10, s4 = 1/5 *)
mixedCheck = Module[{epsValue = -1/3, sv = 1, tv = -3/10, s4v = 1/5, uv, Ea, Eb, Ep, ca, cb, cp, direct, formula, measure, k1, integrate},
  uv = s4v - sv - tv;
  (* recoil rest frame: pRec = (Sqrt[s4], 0); external momenta from the invariants, p along z, pa and pb in the plane phi = 0 and phi = phi_b *)
  Ea = (sv + tv)/(2 Sqrt[s4v]); Eb = (sv + uv)/(2 Sqrt[s4v]); Ep = (sv - s4v)/(2 Sqrt[s4v]);
  cp = {0, 0, 1};
  ca = Module[{c = 1 - (-tv/2)/(Ea Ep)}, {Sqrt[1 - c^2], 0, c}];
  cb = Module[{c = 1 - (-uv/2)/(Eb Ep), cab, x}, cab = 1 - (sv/2)/(Ea Eb);      (* cab = ca.cb from pa.pb = s/2 *)
    x = (cab - ca[[3]] c)/ca[[1]];                                              (* the x component of cb; the y component follows from |cb| = 1 *)
    {x, Sqrt[1 - c^2 - x^2], c}];
  k1[th1_, th2_] := {Sin[th1] Cos[th2], Sin[th1] Sin[th2], Cos[th1]};
  (* angular measure in D - 1 = 3 - 2 eps space dimensions; the integrand 1/(A B') has integrable singularities where k1 is parallel to pa or pb,
     which NIntegrate handles with its convergence warnings suppressed; the comparison tolerance is 10^-3 *)
  integrate[f_] := Quiet[NIntegrate[f, {th1, 0, Pi}, {th2, 0, Pi}, WorkingPrecision -> 20, PrecisionGoal -> 6, MaxRecursion -> 30], {NIntegrate::ncvb, NIntegrate::slwcon}];
  measure = integrate[Sin[th1]^(1 - 2 epsValue) Sin[th2]^(-2 epsValue)];
  (* 1/(A B') with A = -2 pa.k1 = -(s+t)(1 - ca.n)/2, B' = -2 pb.k2 = -(s+u)(1 + cb.n)/2 in the recoil frame *)
  direct = integrate[Sin[th1]^(1 - 2 epsValue) Sin[th2]^(-2 epsValue)/((-(sv + tv) (1 - ca.k1[th1, th2])/2) (-(sv + uv) (1 + cb.k1[th1, th2])/2))]/measure;
  (* the master Omega(1 - zAB) as its exact Euler integral (the eps series of F3 is truncated and would be off by a few percent at eps = -1/3) *)
  formula = With[{zv = zAB /. {s -> sv, t -> tv, u -> uv, s4 -> s4v}},
    N[sign[pa] sign[pb]/(twoPRec[pa] twoPRec[pb]) /. {s -> sv, t -> tv, u -> uv, s4 -> s4v}, 20] *
      NIntegrate[2 (1 - 2 epsValue) x^(-1 - epsValue) (1 + x)^epsValue/(1 - zv + x), {x, 0, Infinity}, WorkingPrecision -> 20]];
  {direct, formula}];
Require[Abs[mixedCheck[[1]]/mixedCheck[[2]] - 1] < 10^-3, "The mixed-family angular master formula fails its numerical check", mixedCheck];
mixedCheck


(* ::Section:: *)
(*Loop masters*)


(* ::Text:: *)
(*With the measure d^D loop/(i Pi^(D/2)) the massless bubble with external momentum squared P^2 < 0 is (-P^2)^(-eps) Gamma(eps) Gamma(1-eps)^2/Gamma(2-2eps) (Euler integral as in the lepton-quark process), and the massless box with invariants S, T < 0 has the Feynman-parameter representation (Cheng-Wu: x4 = 1)*)
(**)
(*   Box(S,T) = Gamma(2+eps) integral dx1 dx2 dx3 (1 + x1 + x2 + x3)^(2 eps) (-S x1 x3 - T x2)^(-2-eps),*)
(**)
(*evaluated by SubTropica as a series in eps for -S = 1 (the scale is restored by dimensional analysis, Box ~ S^(-2-eps)) and -T/-S = r.  In the physical region s > 0 the invariants become s -> s + i0 (the Feynman prescription); the real parts of the continued masters enter the cross section: (-s)^(-eps) -> s^(-eps) cos(Pi eps), log(-s) -> log(s) - i Pi inside the logarithms of the box, whose real part is checked against the standard formula Re Box(s,t) = (2 rGamma/(s t)) [ (s^(-eps) cos(Pi eps) + (-t)^(-eps))/eps^2 - log^2(s/(-t))/2 ] + O(eps) with rGamma = Gamma(1+eps) Gamma(1-eps)^2/Gamma(1-2eps).  The boxes in D+2 (needed for the hatted loop momentum, only in LL) are finite and evaluated the same way.*)


(* ::Input:: *)
loopMasters = <|
  "bubble" -> EulerIntegral[Gamma[eps], x^(-eps) (1 + x)^(2 eps - 2), {x}, {}, 2, eps],
  "bubbleD2" -> EulerIntegral[Gamma[eps - 1], x^(1 - eps) (1 + x)^(2 eps - 4), {x}, {}, 2, eps]|>;
(* SubTropica writes one constant of the order-eps^2 term as the multiple polylogarithm Mpl[{1,2},{1,1}] = Li_{1,2}(1,1) = Sum_{0<k1<k2} 1/(k1 k2^2) = zeta(2,1), which equals zeta(3) by Euler's identity *)
multiplePolylogValues = {HyperIntica`Mpl[{1, 2}, {1, 1}] -> Zeta[3]};
loopMasters = loopMasters /. multiplePolylogValues;
Require[FreeQ[loopMasters, HyperIntica`Mpl], "An unevaluated multiple polylogarithm remains in a loop master", Cases[loopMasters, _HyperIntica`Mpl, Infinity]];
(* Quiet: the order-eps^2 coefficient contains PolyGamma[2, 1] + 2 Zeta[3], which is zero but which Series cannot decide (Series::ztest1); FullSimplify in check[] does *)
check[FullSimplify[Quiet[Normal[Series[loopMasters["bubble"] - Gamma[eps] Gamma[1 - eps]^2/Gamma[2 - 2 eps], {eps, 0, 2}]], Series::ztest1]]];
check[FullSimplify[Quiet[Normal[Series[loopMasters["bubbleD2"] - Gamma[eps - 1] Gamma[2 - eps]^2/Gamma[4 - 2 eps], {eps, 0, 2}]], Series::ztest1]]];
(* The boxes are needed through order eps^0 only: the reduction coefficients of the box masters are regular at d = 4 (checked below), unlike those of
   the bubbles.  SubTropica writes the Euclidean box (r > 0) with Log[-r], Log[1/r] and PolyLog[2, 1 + r], whose imaginary parts cancel: they are
   rewritten through Log[r], Log[1 + r] and PolyLog[2, -r] (Li2(1+r) = Pi^2/6 - Li2(-r) - Log[1+r] Log[-r], Log[-r] = Log[r] + I Pi sigma with the
   branch sigma = +1 or -1 of SubTropica's continuation, whichever makes the result real), and the real result is kept. *)
realEuclidean[expr_] := Module[{candidates, real},
  candidates = Table[Expand[expr /. PolyLog[2, 1 + r] -> Pi^2/6 - PolyLog[2, -r] - Log[1 + r] (Log[r] + I Pi sigma) /. {Log[-r] -> Log[r] + I Pi sigma, Log[1/r] -> -Log[r]}], {sigma, {1, -1}}];
  imaginaryPart[e_] := Simplify[(e - (e /. Complex[x_, y_] :> Complex[x, -y]))/(2 I)];
  real = Select[candidates, imaginaryPart[#] === 0 &];
  Require[real =!= {}, "The Euclidean box is not real for either branch", candidates];
  Collect[(First[real] + (First[real] /. Complex[x_, y_] :> Complex[x, -y]))/2, eps, Simplify]];
(* Quiet: the series expansion inside EulerIntegral cannot decide the sign of combinations of Log[-r], Log[1/r] and Log[r] that vanish for r > 0 (Series::ztest1); realEuclidean resolves them *)
boxEuclid = realEuclidean[Quiet[EulerIntegral[Gamma[2 + eps], (1 + x1 + x2 + x3)^(2 eps) (x1 x3 + r x2)^(-2 - eps), {x1, x2, x3}, {r}, 0, eps], Series::ztest1]]


(* ::Input:: *)
boxD2Euclid = realEuclidean[Quiet[EulerIntegral[Gamma[1 + eps], (1 + x1 + x2 + x3)^(2 eps - 2) (x1 x3 + r x2)^(-1 - eps), {x1, x2, x3}, {r}, 0, eps], Series::ztest1]]


(* ::Text:: *)
(*Check of the Euclidean box against the literature expansion (Bern, Dixon and Kosower; here S = -1, r = T/S): Box = (2 rGamma/(S T)) [ ((-S)^(-eps) + (-T)^(-eps))/eps^2 - log^2(S/T)/2 - Pi^2/2 ] + O(eps); the finite D+2 box is checked numerically against its Feynman-parameter integral at eps = 0.*)


(* ::Input:: *)
boxLiterature = 2 Gamma[1 + eps] Gamma[1 - eps]^2/Gamma[1 - 2 eps]/r ((1 + r^(-eps))/eps^2 - Log[1/r]^2/2 - Pi^2/2);
check[Simplify[Normal[Series[boxEuclid - boxLiterature, {eps, 0, 0}]], r > 0]];
(* numerical check of the D+2 box, which is finite: at eps = 0 its order-eps^0 value is exact and is compared with the Feynman-parameter integral
   itself (the four-dimensional box has infrared poles, and its Feynman-parameter integral converges only for eps well below zero, where the
   truncated series is useless; it is checked against the literature above) *)
boxNumericCheck = Module[{rv = 3/7, direct, series},
  (* the integrand is singular at the origin (integrable); NIntegrate reaches 10^-7 but warns about its error estimate, which is suppressed *)
  direct = Quiet[NIntegrate[(1 + x1 + x2 + x3)^(-2) (x1 x3 + rv x2)^(-1), {x1, 0, Infinity}, {x2, 0, Infinity}, {x3, 0, Infinity},
    WorkingPrecision -> 20, PrecisionGoal -> 10, MaxRecursion -> 30], {NIntegrate::eincr, NIntegrate::slwcon}];
  series = N[boxD2Euclid /. {eps -> 0, r -> rv}, 20];
  Require[Abs[direct/series - 1] < 10^-6, "The D+2 box disagrees with its Feynman-parameter integral", {direct, series}]; {direct, series}];
boxNumericCheck


(* ::Text:: *)
(*Continuation to the physical region.  The reduced virtual amplitude contains, after the ordinary reduction, the bubbles in s, t, u and the boxes in (s,t), (s,u), (t,u) (and their D+2 versions).  With t, u < 0 only the s-channel objects are continued: (-s)^(-eps) -> s^(-eps) e^(i Pi eps), i.e. Re = s^(-eps) cos(Pi eps); for the boxes, the Euclidean series (a function of logs and dilogarithms of r = T/S) is written in terms of log(-S), log(-T) and polylogarithms with arguments that stay away from their cuts, and log(-s) -> log(s) - i Pi is applied; the real part is taken.  The (s,t) and (s,u) boxes: Box(s,t) = s^(-2-eps) x Box(S=-1, r = t/s) ... with the explicit scaling Box(S,T) = (-S)^(-2-eps) Box(-1, T/S); the (t,u) box has no continuation.*)


(* ::Input:: *)
(* the scaled masters: a one-loop integral with mass dimension -2-eps (box) or -eps (bubble) *)
bubbleOf[P2_] := (-P2)^(-eps) loopMasters["bubble"];
bubbleD2Of[P2_] := (-P2)^(1 - eps) loopMasters["bubbleD2"];
boxOf[S_, T_] := (-S)^(-2 - eps) (boxEuclid /. r -> T/S);
boxD2Of[S_, T_] := (-S)^(-eps) (boxD2Euclid /. r -> T/S);
(* continuation s -> s e^(i Pi) in the physical region (t, u < 0): applied to every function of -s, then the real part *)
continue[expr_] := Module[{e = expr}, 
  e = e /. (-s)^a_ :> s^a Exp[I Pi a];
  e = e /. Log[(-s)] -> Log[s] - I Pi /. Log[-s] -> Log[s] - I Pi;
  e = e /. Log[t/s] -> Log[-t/s] + I Pi /. Log[u/s] -> Log[-u/s] + I Pi /. Log[s/t] -> Log[-s/t] - I Pi /. Log[s/u] -> Log[-s/u] - I Pi;
  e];
realPart[expr_] := ComplexExpand[Re[expr], TargetFunctions -> {Re, Im}] /. Im[_] -> 0 /. Re[x_] :> x;


(* ::Text:: *)
(*Identification of the loop masters that Kira returns: a master with two propagators is the bubble of the momentum difference of the two (s, t or u channel), one with all four is the box of its family; the invariants of the family's box are S = (q3 - q1)^2 and T = (q4 - q2)^2 for the propagator momenta q1 ... q4 (the Feynman-parameter representation above has F = -S x1 x3 - T x2 x4).  Triangles do not survive the reduction (massless, they are combinations of bubbles) and would be reported.*)


(* ::Input:: *)
virtualPropagators = <|"virtual_ST" -> {loop, loop + pa, loop + pa + pb, loop + p}, "virtual_SU" -> {loop, loop + pa, loop + pa + pb, loop + pa + pb - p},
  "virtual_TU" -> {loop, loop + pa, loop + pa - p, loop + pa - p + pb}|>;
invariant[q_] := Simplify[ExpandScalarProductRules[q]];
ExpandScalarProductRules[q_] := Module[{e = Expand[q]}, (* (q)^2 for q a combination of pa, pb, p *)
  Expand[e^2] /. {pa^2 -> 0, pb^2 -> 0, p^2 -> 0, pa pb -> s/2, pa p -> -t/2, pb p -> -u/2, loop -> 0} /. u -> -s - t];
loopMasterValue[KiraMaster[family_, indices_, shift_]] := Module[{q = virtualPropagators[family], present, P2, S, T},
  present = Position[indices, 1][[All, 1]];
  Which[
    Length[present] == 2 && Count[indices, 0] == 2,
      P2 = invariant[q[[present[[1]]]] - q[[present[[2]]]]];
      Require[P2 =!= 0, "A scaleless bubble appeared as a master", {family, indices}];
      If[shift == 0, bubbleOf[P2], bubbleD2Of[P2]],
    indices === {1, 1, 1, 1},
      S = invariant[q[[3]] - q[[1]]]; T = invariant[q[[4]] - q[[2]]];
      If[shift == 0, boxOf[S, T], boxD2Of[S, T]],
    True, Require[False, "Unexpected loop master", {family, indices}]]];
loopMasterRules = Table[m -> loopMasterValue[m], {m, Union[Cases[virtualReduced, KiraMaster[f_, __] /; StringStartsQ[f, "virtual"], Infinity]]}];
Column[loopMasterRules /. (a_ -> b_) :> a -> Short[b, 3]]


(* ::Text:: *)
(*Continuation of the s-channel objects to s > 0 and the real part.  The (s,t) and (s,u) boxes are written as (-s)^(-2-eps) times the Euclidean function of r = T/S = (-t)/(-s) (or (-u)/(-s)); with s -> s e^(-i Pi) (Feynman prescription s + i0, i.e. -s -> s e^(-i Pi)) the variable r goes to rp e^(i Pi) with rp = -t/s in (0,1).  The rules below continue the functions that occur in the SubTropica result; an unknown function stops the notebook.  Then the real part is taken with ComplexExpand (eps is real).*)


(* ::Input:: *)
boxFunctions = Union[Cases[{boxEuclid, boxD2Euclid}, (Log | PolyLog)[___], Infinity]]


(* ::Input:: *)
continueBox[expr_, rp_] := Module[{e = expr},
  e = e /. {Log[r] -> Log[rp] + I Pi, Log[1 + r] -> Log[1 - rp], PolyLog[2, -r] -> PolyLog[2, rp], PolyLog[2, 1 + r] :> (Pi^2/6 - Log[rp] Log[1 - rp] - PolyLog[2, 1 - rp] + I Pi Log[1 - rp] /. Log[rp] -> Log[rp] + I Pi)};
  Require[FreeQ[Cases[e, (Log | PolyLog)[___], Infinity], r], "A function of the box was not continued", Cases[e, (Log | PolyLog)[___], Infinity]];
  e = e /. r -> -rp;                 (* the rational coefficients: r = T/S = -rp *)
  e];
(* every master is rebuilt with its continuation *)
continued[KiraMaster[family_, indices_, shift_]] := Module[{q = virtualPropagators[family], present, P2, S, T, box, rp},
  present = Position[indices, 1][[All, 1]];
  Which[
    Length[present] == 2,
      P2 = invariant[q[[present[[1]]]] - q[[present[[2]]]]];
      If[P2 === s, (* s channel: (-s)^(-eps) -> s^(-eps) cos(Pi eps) *) 
        If[shift == 0, s^(-eps) Cos[Pi eps] loopMasters["bubble"], -s^(1 - eps) Cos[Pi eps] loopMasters["bubbleD2"]],
        If[shift == 0, bubbleOf[P2], bubbleD2Of[P2]]],
    indices === {1, 1, 1, 1},
      S = invariant[q[[3]] - q[[1]]]; T = invariant[q[[4]] - q[[2]]];
      If[S === s || T === s,
        (* arrange S as the s channel *)
        If[T === s, {S, T} = {T, S}];
        rp = -T/s;
        box = If[shift == 0, boxEuclid, boxD2Euclid];
        realPart[Exp[-I Pi (If[shift == 0, -2 - eps, -eps])] s^(If[shift == 0, -2 - eps, -eps]) continueBox[box, rp]],
        If[shift == 0, boxOf[S, T], boxD2Of[S, T]]]]];
physicalLoopMasters = Table[m -> continued[m], {m, loopMasterRules[[All, 1]]}];
Column[physicalLoopMasters /. (a_ -> b_) :> a -> Short[b, 3]]


(* ::Text:: *)
(*Check of the continued (s,t) box against the literature: Re Box(s,t) = (2 rGamma/(s t)) [ (s^(-eps) cos(Pi eps) + (-t)^(-eps))/eps^2 - log^2(s/(-t))/2 ] + O(eps).*)


(* ::Input:: *)
boxSTphysical = continued[KiraMaster["virtual_ST", {1, 1, 1, 1}, 0]];
boxSTliterature = 2 Gamma[1 + eps] Gamma[1 - eps]^2/Gamma[1 - 2 eps]/(s t) ((s^(-eps) Cos[Pi eps] + (-t)^(-eps))/eps^2 - Log[s/(-t)]^2/2);
Require[Simplify[Normal[Series[boxSTphysical - boxSTliterature, {eps, 0, 0}]], s > 0 && t < 0 && t > -s] === 0, "The continued box differs from the literature formula"];
"box continuation checked"


(* ::Section:: *)
(*Insert the physical masters and save*)


(* ::Input:: *)
realIntegrated = Association @ Table[spin -> (realReduced[spin] /. physicalMasters), {spin, spins}];
virtualIntegrated = Association @ Table[spin -> (virtualReduced[spin] /. physicalLoopMasters), {spin, spins}];
Require[FreeQ[{realIntegrated, virtualIntegrated}, KiraMaster], "A master integral was not inserted"];
reduced = <|"Real" -> realIntegrated, "Virtual" -> virtualIntegrated, "Masters" -> masters, "LoopMasters" -> loopMasters,
  "MasterFunctions" -> masterFunctions, "ThresholdValues" -> thresholdValues, "RealRules" -> realRules, "VirtualRules" -> virtualRules, "Route" -> route, "Class" -> class|>;
SaveResult["ReductionMetrics", <|
  "ContributionRows" -> Map[Length, realTerms],
  "UniqueRealTargets" -> Length[realTargets], "UniqueVirtualTargets" -> Length[virtualTargets],
  "RealMastersIncludingShift" -> Length[Union[Cases[realReduced, _KiraMaster, Infinity]]],
  "VirtualMastersIncludingShift" -> Length[Union[Cases[virtualReduced, _KiraMaster, Infinity]]],
  "SubTropicaSeconds" -> Total[Lookup[$SubTropicaTimings, "Seconds", 0]]|>, route <> "/class-" <> class];
SaveResult["Reduced", reduced, route <> "/class-" <> class];
Keys[reduced]
