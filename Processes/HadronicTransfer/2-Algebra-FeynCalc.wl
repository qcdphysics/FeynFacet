(* ::Title:: *)
(*Transverse spin transfer in proton-proton collisions, step 2 (FeynCalc route): spin and colour algebra*)


(* ::Text:: *)
(*Input: the amplitudes of step 1 for the class chosen below.  Output: Results/FeynCalc/class-<n>/Algebra.wl with, for the two transverse components XX and NN defined below,*)
(*  Born:      the Born transfer square, summed over the final spins and colours of the unobserved partons, averaged over the colours of a and over the spin and colour of b, in D dimensions;*)
(*  Virtual:   the interference 2 Re(M_1loop M_Born^*) with the same projections, written in terms of scalar one-loop integrals (families VirtualST, VirtualSU, VirtualTU as in Processes/Hadronic) with D-dimensional coefficients;*)
(*  Real:      the sum over the real-emission processes of the class, each weighted as in classes.wl, of the projected |M|^2 with the unobserved momenta k1, k2 kept explicit (step 3 integrates them);*)
(*  BornTable: the D-dimensional transfer squares of the 2 -> 2 channels of classes.wl (TransferBornChannels), the lower Born terms of the collinear counterterms of step 4.*)
(**)
(*The observable (docs/TRANSVERSE-PP.md).  The incoming parton a carries a transverse spin s_a, the transverse spin s_c of the observed parton c is measured, parton b is unpolarized and the other final partons are summed.  The partonic quantity is sigma(c: s_c) - sigma(c: -s_c) at fixed s_a, averaged as ½[sigma(s_a) - sigma(-s_a)], i.e. the spin-averaged square with the incoming density pa-slash/2 of a replaced by pa-slash gamma5 s_a-slash/2 (the same expression for an incoming antiquark, class 17: v vbar = p-slash (1 + gamma5 s-slash)/2) and the final spin sum p-slash of c replaced by the projection p-slash gamma5 s_c-slash (no ½), as in Processes/Collins.  The spin vectors are written in the basis of the scattering plane of pa, pb, p: s_a = cos(phi_a) e_X + sin(phi_a) n, s_c = cos(phi_c) e_C + sin(phi_c) n, with the normal n (n.pa = n.pb = n.p = 0, n^2 = -1), the in-plane axis e_X of the incoming partons (orthogonal to pa and pb, along the transverse momentum of p: Processes/Collins' aX with l -> pb, k -> pa) and the in-plane axis e_C of the observed parton (orthogonal to p and pa + pb: Collins' bX relabelled),*)
(**)
(*   eX = (u/s) pa + (t/s) pb + p,   eX^2 = -t u/s;      eC = -u/(t+u) pa + t/(t+u) pb - s (u-t)/(t+u)^2 p,   eC^2 = -s t u/(t+u)^2,*)
(**)
(*both unnormalized (the XX component is divided by the product of their lengths, -t u/(t+u)).  The squared matrix elements are the four components XX, NN, XN, NX, the coefficients of cos phi_a cos phi_c, sin sin, cos sin, sin cos; the mixed components must vanish by parity (checked), XX and NN are both kept.*)
(**)
(*HVBM scheme.  The physical momenta pa, pb, p and the spin vectors are four-dimensional; the unobserved momenta k1, k2 and the loop momentum have (D-4)-dimensional parts (spe3 = k1hat.k1hat in the real emission, loopHat2 in the loop), and gamma5 anticommutes with the four-dimensional gamma matrices only ($LeviCivitaSign = -1).  Every projection contains two gamma5's, so the traces produce pairs of Levi-Civita tensors, which FeynCalc contracts into determinants of four-dimensional scalar products: the normal n then appears only through n.n = -1 and through its products with the momenta that have a component out of the scattering plane, n.k1 = -n.k2 (real emission) and n.loop (virtual).  Terms odd in n.k1 vanish after the phase-space integration (they are odd under the reflection of the plane), and (n.k1)^2 is a rational function of the invariants: k1 = k1par + k1n + k1hat (in the plane of pa, pb, p; along n; in the extra dimensions) with k1^2 = 0 gives (n.k1)^2 = k1par^2 + spe3, k1par^2 = d.G^-1.d with G the Gram matrix of pa, pb, p and d = (pa.k1, pb.k1, p.k1) (Processes/Collins, normalReduction).  In the loop, n.loop is reduced by the symmetric integration over the (D-3)-dimensional space orthogonal to pa, pb, p (section "Virtual").*)
(**)
(*Invariants as in Processes/Hadronic: s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2; real emission: s4 = (k1+k2)^2 = s + t + u, and the three independent angular denominators A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 with A + B + C = -s4; the denominators with k2 follow from k2 = pa + pb - p - k1 (step 3).*)


(* ::Section:: *)
(*Load*)


(* ::Text:: *)
(*Load the FeynFacet package (tool lookup, amplitude files, checks); FeynCalc and FeynArts are loaded in the next cell.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*Choices: the class; the number of parallel kernels for the squares ("auto" = one per processor core, at most 8; "1" = none); and whether the interference files of an earlier run of this same, unchanged file may be reused ("no" recomputes everything; "yes" only continues an interrupted run: a product file is reused only if its companion file <name>.input holds the fingerprint of the current inputs, the amplitudes of step 1 and this step file; a file written for other inputs, or without that companion (for example, after an interrupted write), is computed again).*)


(* ::Input:: *)
LoadFeynCalc[];
(* FeynCalc expresses this memoization threshold in MB; transient algebra can use more memory. *)
FeynCalc`$FCMemoryAvailable = Min[FeynCalc`$FCMemoryAvailable, Floor[1024^3/10^6]];
class = "1";
workers = "auto";
resume = "no";
loaded = LoadAmplitudes["class-" <> class];
amplitudes = loaded["FeynCalc"]; spec = loaded["Spec"]; weights = loaded["Weights"];
spec


(* ::Section:: *)
(*Conventions*)


(* ::Text:: *)
(*The conventions and kinematics are installed by functions, because every parallel kernel installs them too.  One unobserved momentum is not independent: inside the amplitudes k2 is replaced by pa + pb - p - k1 (real emission) and k1 by pa + pb - p (Born and one loop), except in the external spinors, whose spin sums give the slashed momentum itself.  The spin axes nn, eX, eC are four-dimensional vectors of their own whose scalar products with every momentum (including the loop momentum and the unobserved k1, k2) are set from the definitions above through the scalar products already in force (Processes/Collins' defineAxes); n.k1 = nk1 (real emission, n.k2 = -nk1) and n.loop = nl are symbols that the reductions below remove.*)


(* ::Input:: *)
conventions := (FCSetDiracGammaScheme["BMHV"]; $LeviCivitaSign = -1);
conventions;
axisCoefficients = <|eX -> {u/s, t/s, 1}, eC -> {-u/(t + u), t/(t + u), -s (u - t)/(t + u)^2}|>;    (* on the basis {pa, pb, p} *)
(* the axes: every product with the momenta of the list and among themselves; n.x = 0 for the in-plane momenta, n.n = -1 *)
defineAxes[momenta_List, normalRules_List] := Module[{base = {pa, pb, p}, product},
  product[c1_, c2_] := Factor[Sum[c1[[i]] c2[[j]] SP[base[[i]], base[[j]]], {i, 3}, {j, 3}]];
  Scan[(Momentum[#, D - 4] = 0) &, {eX, eC, nn}];
  Do[SP[axis, x] = SPD[axis, x] = Factor[Sum[axisCoefficients[axis][[i]] SP[base[[i]], x], {i, 3}]], {axis, {eX, eC}}, {x, momenta}];
  Do[SP[pair[[1]], pair[[2]]] = SPD[pair[[1]], pair[[2]]] = product[axisCoefficients[pair[[1]]], axisCoefficients[pair[[2]]]], {pair, {{eX, eX}, {eC, eC}, {eX, eC}}}];
  SPD[nn] = -1; SP[nn] = -1; SP[eX, nn] = SPD[eX, nn] = 0; SP[eC, nn] = SPD[eC, nn] = 0;
  Scan[(SPD[nn, #] = 0; SP[nn, #] = 0) &, base];
  Do[SPD[nn, rule[[1]]] = rule[[2]]; SP[nn, rule[[1]]] = rule[[2]], {rule, normalRules}]];
axisLengths = -t u/(t + u);                                                                        (* |eX| |eC| = sqrt(eX^2 eC^2), t, u < 0 *)
bornKinematics := (FCClearScalarProducts[]; SetMandelstam[s, t, u, pa, pb, -p, -k1, 0, 0, 0, 0];
  Scan[(Momentum[#, D - 4] = 0) &, {pa, pb, p, k1}];
  SP[loop, nn] = nl; SPD[loop, nn] = nl;                                                            (* the loop momentum: its products with the axes through those with pa, pb, p *)
  Do[SP[axis, loop] = SPD[axis, loop] = Sum[axisCoefficients[axis][[i]] SPD[{pa, pb, p}[[i]], loop], {i, 3}], {axis, {eX, eC}}];
  defineAxes[{pa, pb, p, k1}, {k1 -> 0}]);
realKinematics := (FCClearScalarProducts[];
  SPD[pa] = 0; SPD[pb] = 0; SPD[p] = 0; SPD[k1] = 0; SPD[k2] = 0;
  SPD[pa, pb] = s/2; SPD[pa, p] = -t/2; SPD[pb, p] = -u/2;
  SPD[pa, k1] = -A/2; SPD[pb, k1] = -B/2; SPD[p, k1] = C/2;
  SPD[pa, k2] = (s + t + A)/2; SPD[pb, k2] = (s + u + B)/2; SPD[p, k2] = (-t - u - C)/2; SPD[k1, k2] = (s + t + u)/2;
  Scan[(Momentum[#, D - 4] = 0) &, {pa, pb, p}];
  (* hatted parts of the unobserved momenta, and their four-dimensional products (SP = SPD - SPE; the physical momenta have no hatted parts) *)
  Pair[Momentum[k1, D - 4], Momentum[k1, D - 4]] = spe3; Pair[Momentum[k2, D - 4], Momentum[k2, D - 4]] = spe3;
  Pair[Momentum[k1, D - 4], Momentum[k2, D - 4]] = -spe3;
  SP[k1] = -spe3; SP[k2] = -spe3; SP[k1, k2] = (s + t + u)/2 + spe3;
  SP[pa] = 0; SP[pb] = 0; SP[p] = 0; SP[pa, pb] = s/2; SP[pa, p] = -t/2; SP[pb, p] = -u/2;
  SP[pa, k1] = -A/2; SP[pb, k1] = -B/2; SP[p, k1] = C/2;
  SP[pa, k2] = (s + t + A)/2; SP[pb, k2] = (s + u + B)/2; SP[p, k2] = (-t - u - C)/2;
  defineAxes[{pa, pb, p, k1, k2}, {k1 -> nk1, k2 -> -nk1}]);
bornConservation = k1 -> pa + pb - p;
realConservation = k2 -> pa + pb - p - k1;
incoming = spec["Born"] /. None -> spec["Real"][[1, 1]];       (* the two incoming partons *)
incoming = incoming[[1 ;; 2]];
(* a is a quark or an antiquark in every class (its density carries the ½ of the colour-and-spin average 1/(2 Nc)); b is averaged over spin and colour *)
average[in_List] := Times @@ (Switch[in[[#]], "g", 1/(D - 2)/(Nc^2 - 1), _, 1/(2 Nc)] & /@ {1, 2});
Require[incoming[[1]] =!= "g", "The polarized parton a must be a quark or an antiquark"];
colourFactors[expr_] := Factor[expr /. SUNN -> Nc] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
spinVectors = <|"XX" -> {eX, eC}, "NN" -> {nn, nn}, "XN" -> {eX, nn}, "NX" -> {nn, eC}|>;     (* {s_a, s_c} of the four components *)
spinCases = Keys[spinVectors];
outputSpins = {"XX", "NN"};
bornKinematics;
Require[Simplify[{SP[eX, pa], SP[eX, pb], SP[eC, p], SP[eC, pa] + SP[eC, pb], SP[eX] + t u/s, SP[eC] + s t u/(t + u)^2, SP[nn, eX], SP[nn, eC], SP[nn]}] === {0, 0, 0, 0, 0, 0, 0, 0, -1},
  "The axes are not what the text says"];
{incoming, average[incoming]}


(* ::Section:: *)
(*The squaring method*)


(* ::Text:: *)
(*A square sum_{i,j} M_i M_j^* of many diagrams is not formed as one product.  Five steps keep the expressions small (Processes/Hadronic):*)
(**)
(*1. Every amplitude is simplified once, before any product: Lorentz contraction, Dirac algebra on the open spinor chains, expansion of the scalar products with one independent unobserved momentum, and the conditions eps(l).l = 0 and eps(l).r = 0 for every gluon l, where r is the reference vector of the same gluon's polarization sum used later.  For the gluon whose momentum was eliminated by momentum conservation, eps(l).l = 0 is imposed in its eliminated form.*)
(*2. Colour is separated from the Lorentz algebra: each amplitude is written as sum_x c_x L_x with distinct colour structures c_x.  The colour factor of c_x c_y^* is a closed colour network; it is evaluated once and remembered; equal networks with differently named summed indices are recognized by renaming the indices canonically.*)
(*3. Only the products with i <= j are computed: for i < j, M_i M_j^* + M_j M_i^* = 2 Re(M_i M_j^*) (all invariants are real).*)
(*4. After the fermion spin sums with the transverse projections, the gluon polarization sums are contracted one gluon at a time before evaluating the Dirac traces: the terms are grouped by the factors that carry that gluon's two indices, and each distinct tensor is contracted with the polarization sum only once (and remembered).  In the real kinematics the normal component n.k1 is then reduced (normalReduction).*)
(*5. The products (i, j) are independent: they are distributed over parallel kernels, and each is written to its own file with the time it took, so that the progress is visible and an interrupted run can be continued.*)
(**)
(*The spin sums.  The spinors of a (momentum pa) and of c (momentum p) are tagged so that FermionSpinSum writes their slashed momenta as paSpin-slash and pSpin-slash, which are then replaced by the density pa-slash gamma5 s_a-slash and the projection p-slash gamma5 s_c-slash of the component (the ½ of the density is in average[]; the same replacement for a quark and an antiquark, section "Observable").  The outgoing gluons are summed with the physical sum -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r), r = pa for the observed parton and r = pa or pb for the unobserved ones, an incoming gluon b with the same sum and 1/(D-2); the results are independent of the reference vectors by gauge invariance (checked below for the real emission).*)


(* ::Input:: *)
(* ---- 1. one amplitude ---- *)
gluonsOf[expr_] := Union[Cases[expr, Polarization[l_, __] :> l, Infinity]];
transversality[references_] := Flatten[Table[{
    Pair[Momentum[Polarization[l, ___], ___], Momentum[l, ___]] -> 0,
    Pair[Momentum[Polarization[l, ___], ___], Momentum[references[l], ___]] -> 0}, {l, Keys[references]}]];
(* eps(x).x = 0 for the eliminated momentum x = sum_q c_q q: solve it for eps(x).q of the last momentum q that is not x's reference vector *)
eliminatedTransversality[rule_, references_] := Module[{x = rule[[1]], terms = termList[Expand[rule[[2]]]], momenta, coefficient, solved},
  momenta = Cases[terms, c_. q_Symbol :> q];
  solved = Last[DeleteCases[momenta, references[x]]];
  coefficient[q_] := Coefficient[rule[[2]], q];
  Pair[Momentum[Polarization[x, h___], d1___], Momentum[solved, d2___]] :>
    -Sum[coefficient[q] Pair[Momentum[Polarization[x, h], d1], Momentum[q, d2]], {q, DeleteCases[momenta, solved]}]/coefficient[solved]];
(* momentum conservation inside propagators and vertices, not in the external spinors (Spinor[...] is matched first and left unchanged) *)
conserve[expr_, rule_] := expr /. {sp_Spinor :> sp, Momentum[q_, dim___] /; FreeQ[q, Polarization] && !FreeQ[q, rule[[1]]] :>
  Momentum[Expand[q /. rule], dim]};
prepare[amp_, references_, rule_] := Module[{a = conserve[amp, rule]},
  a = FeynAmpDenominatorExplicit[a];
  a = DiracSimplify[Contract[a]];
  a = ExpandScalarProduct[a] /. eliminatedTransversality[rule, references] /. transversality[references];
  Expand[a]];

(* ---- 2. colour ---- *)
colourHead = SUNTF | SUNT | SUNF | SUNDelta | SUNFDelta | SUNTrace;
factorList[term_] := If[Head[term] === Times, List @@ term, {term}];
termList[expr_] := If[Head[expr] === Plus, List @@ expr, {expr}];
colourPart[term_] := Times @@ Select[factorList[term], !FreeQ[#, SUNIndex | SUNFIndex] &];
(* {{c_x, L_x}, ...}: the distinct colour structures of an amplitude with their Lorentz parts *)
colourGroups[amp_] := Module[{groups = GatherBy[termList[Expand[amp]], colourPart]},
  Table[{colourPart[First[g]], Total[(#/colourPart[#]) & /@ g]}, {g, groups}]];
colourCache = <||>;
colourValue[network_] := Module[{key = FCCanonicalizeDummyIndices[network, Head -> {SUNIndex, SUNFIndex}, FCI -> True,
    SUNIndexNames -> Table[Symbol["cA" <> ToString[k]], {k, 30}], SUNFIndexNames -> Table[Symbol["cF" <> ToString[k]], {k, 30}]], value},
  If[KeyExistsQ[colourCache, key], Return[colourCache[key]]];
  value = SUNSimplify[key /. SUNF[a_, b_, c_] :> -2 I (SUNTrace[SUNT[a] . SUNT[b] . SUNT[c]] - SUNTrace[SUNT[b] . SUNT[a] . SUNT[c]]),
    Explicit -> True, SUNTraceEvaluate -> True, SUNNToCACF -> False];
  value = Factor[value /. SUNN -> Nc];
  If[!FreeQ[value, colourHead | SUNIndex | SUNFIndex], Return[Failure["Colour", <|"Network" -> key, "Value" -> value|>]]];
  colourCache[key] = value];

(* ---- 4. spin sums, transverse projections and gluon polarization sums ---- *)
polarizationSum[l_, mu_, nu_, r_] := -MTD[mu, nu] + (FVD[l, mu] FVD[r, nu] + FVD[r, mu] FVD[l, nu])/SPD[l, r];
tag[l_] := Symbol[ToString[l] <> "Spin"];
indexNames[l_] := {Symbol["mu" <> ToString[l]], Symbol["nu" <> ToString[l]]};
tensorPart[term_, {mu_, nu_}] := Times @@ Select[factorList[term], !FreeQ[#, LorentzIndex[mu | nu, ___]] &];
(* the contraction uses the installed scalar products, which differ between the Born and the real kinematics: the kinematics are part of the key *)
(* Each worker reuses identical contractions across diagram pairs. The explicit
   spin projections and the original BMHV trace evaluation are unchanged. *)
currentKinematics = None;
algebraCacheScope[] := {"HadronicTransfer", class, currentKinematics, "BMHV", $LeviCivitaSign};
initializeAlgebraCaches[] := (
  traceAlgebraCache = CreateBoundedAlgebraCache[
    Function[structure, Contract[DiracSimplify[structure, DiracTraceEvaluate -> True]]], 32768, 128 2^20];
  tensorAlgebraCache = CreateBoundedAlgebraCache[
    Function[product, ExpandScalarProduct[Contract[Expand[product]]]], 8192, 128 2^20]);
initializeAlgebraCaches[];
canonicalTraceStructure[structure_] := FCCanonicalizeDummyIndices[structure,
  FCI -> True, Head -> {LorentzIndex}, DotSimplify -> False, FCTraceExpand -> False,
  LorentzIndexNames -> Table[Symbol["algebraDummy" <> ToString[j]], {j, 40}]];
groupedTrace[expr_] := Module[{rows, groups},
  (* Pure-gluon terms have no trace to evaluate; leave their scalar sums sparse. *)
  If[FreeQ[expr, _DiracTrace | _DiracGamma | _DOT], Return[expr]];
  rows = GroupIdenticalAlgebraTerms[expr, _DiracTrace | _DiracGamma | _LorentzIndex | _Eps | _DOT];
  If[FailureQ[rows], Return[rows]];
  rows = {canonicalTraceStructure[#["Structure"]], #["Coefficient"]} & /@ rows;
  groups = GatherBy[rows, First];
  Total[Table[With[{structure = First[First[group]], coefficient = Together[Total[Last /@ group]]},
    If[coefficient === 0, 0,
      coefficient traceAlgebraCache["Apply"][{algebraCacheScope[], structure}, structure]]], {group, groups}]]];
projectOnce[tensor_, density_] := tensorAlgebraCache["Apply"][
  {algebraCacheScope[], tensor, density}, tensor density];
project[expr_, l_, density_] := Module[{idx = indexNames[l], groups},
  groups = GroupIdenticalAlgebraTerms[expr, LorentzIndex[idx[[1]] | idx[[2]], ___]];
  If[FailureQ[groups], Return[groups]];
  Total[(#["Coefficient"] projectOnce[#["Structure"], density]) & /@ groups]];
(* the real-emission tensor is a polynomial in nk1 of degree at most two: the odd part integrates to zero, nk1^2 = k1par^2 + spe3 *)
gram = {{0, s/2, -t/2}, {s/2, 0, -u/2}, {-t/2, -u/2, 0}};                                            (* pa, pb, p *)
k1Parallel2 = Factor[{-A/2, -B/2, C/2} . Inverse[gram] . {-A/2, -B/2, C/2}];                          (* k1par^2 = d.G^-1.d *)
normalReduction[expr_] := Module[{coefficients = CoefficientList[Expand[expr], nk1]},
  Require[Length[coefficients] <= 3, "The real-emission tensor has a power of n.k1 above two", Length[coefficients] - 1];
  If[coefficients === {}, 0,     (* a product that vanishes identically (CoefficientList[0, nk1] is the empty list) *)
    coefficients[[1]] + If[Length[coefficients] == 3, coefficients[[3]] (k1Parallel2 + spe3), 0]]];
(* the colour-weighted product of the Lorentz parts of one pair: spin sums with the projections, traces, polarization sums *)
spinSum[product_, spinCase_, references_, gluonList_] := Module[{expr = product, sa, sc},
  {sa, sc} = spinVectors[spinCase];
  expr = expr /. (With[{l = #, lTag = tag[#]}, Spinor[c_. Momentum[l, D], m_, o_] :> Spinor[c Momentum[lTag, D], m, o]] & /@ {pa, p});
  expr = expr /. Flatten[Table[{Momentum[Polarization[l, I, ___], D] -> LorentzIndex[indexNames[l][[1]], D],
    Momentum[Polarization[l, -I, ___], D] -> LorentzIndex[indexNames[l][[2]], D]}, {l, gluonList}]];
  expr = FermionSpinSum[expr];
  expr = expr /. {DiracGamma[Momentum[paSpin, D], D] :> GS[pa] . GA[5] . GS[sa], DiracGamma[Momentum[pSpin, D], D] :> GS[p] . GA[5] . GS[sc]};
  (* Contract polarization sums before expanding the same BMHV traces. *)
  (* outgoing gluons first, then the incoming one *)
  Do[expr = project[expr, l, polarizationSum[l, Sequence @@ indexNames[l], references[l]]], {l, Complement[gluonList, {pa, pb}]}];
  Do[expr = project[expr, l, polarizationSum[l, Sequence @@ indexNames[l], references[l]]], {l, Intersection[gluonList, {pa, pb}]}];
  expr = groupedTrace[expr];
  If[!FreeQ[expr, _Failure], Return[FirstCase[expr, _Failure, Failure["Trace", <||>], Infinity]]];
  expr = ExpandScalarProduct[Contract[expr]];
  If[currentKinematics === "Real", normalReduction[expr], expr]];

(* ---- 3. one product (i, j) of colour-separated amplitudes: left = {{c_x, L_x}, ...}, right = {{c_y^*, L_y^*}, ...} ---- *)
interference[left_, right_, spinCase_, references_] := Module[{product},
  product = Total[Flatten[Table[colourValue[x[[1]] y[[1]]] x[[2]] y[[2]], {x, left}, {y, right}]]];
  If[!FreeQ[product, Failure], FirstCase[product, _Failure, Null, Infinity],
    spinSum[product, spinCase, references, gluonsOf[{left, right}]]]];
conjugateGroups[amp_] := colourGroups[ComplexConjugate[amp, FCRenameDummyIndices -> True]];
(* the real part: after expanding, every term is a number times real symbols (the denominators are real propagators), so the real part of the number is kept *)
realPartOf[expr_] := Expand[expr] /. Complex[a_, b_] :> a;
"squaring functions defined"


(* ::Text:: *)
(*Parallel kernels: each loads FeynFacet and FeynCalc, installs the conventions and receives the functions above; each job installs its kinematics (Born or real).  Every job (one product of the squaring method) is written to Results/FeynCalc/class-<n>/pairs/<name>.wl with its time.*)


(* ::Input:: *)
feynFacetFile = FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}];
nWorkers = If[workers === "auto", Min[$ProcessorCount, 8], ToExpression[workers]];   (* "auto": one per processor core, at most 8 *)
If[nWorkers > 1,
  LaunchKernels[nWorkers];
  With[{file = feynFacetFile}, ParallelEvaluate[Get[file]; LoadFeynCalc[];
    FeynCalc`$FCMemoryAvailable = Min[FeynCalc`$FCMemoryAvailable, Floor[1024^3/10^6]]]];
  DistributeDefinitions[conventions, axisCoefficients, defineAxes, bornKinematics, realKinematics, gluonsOf, transversality, eliminatedTransversality, conserve, prepare,
    colourHead, factorList, termList, colourPart, colourGroups, colourValue, polarizationSum, tag, indexNames, tensorPart,
    projectOnce, project, gram, k1Parallel2, normalReduction, spinSum, spinVectors, interference, conjugateGroups, realPartOf,
    algebraCacheScope, initializeAlgebraCaches, canonicalTraceStructure, groupedTrace];
  ParallelEvaluate[conventions; colourCache = <||>; currentKinematics = None; initializeAlgebraCaches[]]];
pairDirectory = FileNameJoin[{ResultsDirectory["FeynCalc/class-" <> class], "pairs"}];
If[!DirectoryQ[pairDirectory], CreateDirectory[pairDirectory]];
(* the fingerprint of the inputs of the product files: the amplitudes of step 1 and this step file (see the choice resume above) *)
inputFingerprint = InputFingerprint[loaded, {
  FileNameJoin[{NotebookDirectory[], "2-Algebra-FeynCalc.wl"}],
  FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "AlgebraReuse.wl"}]}];

runJob[job_] := Module[{file = FileNameJoin[{pairDirectory, job["Name"] <> ".wl"}], value, seconds},
  If[resume === "yes" && ReusableQ[file, inputFingerprint], file,
    (* Algebra caches remain local to this worker and bounded across pairs. *)
    If[job["Kinematics"] === "Born", bornKinematics, realKinematics]; currentKinematics = job["Kinematics"];
    {seconds, value} = AbsoluteTiming[interference[job["Left"], job["Right"], job["Spin"], job["References"]]];
    If[!FreeQ[value, Failure], value,
      Put[<|"Name" -> job["Name"], "Value" -> value, "Seconds" -> seconds|>, file]; MarkReusable[file, inputFingerprint]]]];
runJobs[jobs_] := Module[{files},
  DistributeDefinitions[runJob, pairDirectory, resume, inputFingerprint];
  files = If[nWorkers > 1, ParallelMap[runJob, jobs, Method -> "FinestGrained"], Map[runJob, jobs]];
  Require[FreeQ[files, Failure], "A job failed", FirstCase[files, _Failure, Null, Infinity]];
  Get /@ files];
(* the pairs i <= j of a square, and their sum: diagonal products as they are, i < j twice the real part *)
squareJobs[name_, sides_, kinematics_, references_, spinCase_] := Flatten[Table[
  <|"Name" -> name <> "-" <> spinCase <> "-" <> ToString[i] <> "-" <> ToString[j], "Kinematics" -> kinematics,
    "Left" -> sides[[1, i]], "Right" -> sides[[2, j]], "Spin" -> spinCase, "References" -> references, "Pair" -> {i, j}|>,
  {i, Length[sides[[1]]]}, {j, i, Length[sides[[1]]]}]];
(* a diagonal product |M_i|^2 is real, but its epsilon-tensor terms may leave an imaginary part that vanishes only after simplification:
   that imaginary part is checked to be exactly zero at a rational point and dropped *)
checkPoint = {s -> 1, t -> -3/10, u -> -2/5, A -> -1/7, B -> -2/9, C -> -3/10 + 1/7 + 2/9, D -> 37/10, Nc -> 3, spe3 -> 1/11, gs -> 1, Nf -> 5, nl -> 2/13};
imaginaryPartOf[expr_] := (expr - (expr /. Complex[a_, b_] :> Complex[a, -b]))/2;
square[name_, sides_, kinematics_, references_, spinCase_] := Module[{jobs = squareJobs[name, sides, kinematics, references, spinCase], results},
  results = runJobs[jobs];
  Require[And @@ MapThread[#2["Pair"][[1]] != #2["Pair"][[2]] || (imaginaryPartOf[#1["Value"]] /. checkPoint) === 0 &, {results, jobs}],
    "A diagonal product has a nonzero imaginary part", name];
  Total[MapThread[If[#2["Pair"][[1]] == #2["Pair"][[2]], realPartOf[#1["Value"]], 2 realPartOf[#1["Value"]]] &, {results, jobs}]]];
sidesOf[amps_, references_, rule_] := With[{prepared = prepare[#, references, rule] & /@ amps}, {colourGroups /@ prepared, conjugateGroups /@ prepared}];
(* the normalization of a component: XX by the product of the axis lengths, NN by nothing; the mixed components are only checked to vanish *)
normalize[expr_, spinCase_] := If[spinCase === "XX", expr/axisLengths, expr];
"parallel kernels: " <> ToString[Max[1, Length[Kernels[]]]]


(* ::Section:: *)
(*The transfer Born table*)


(* ::Text:: *)
(*The D-dimensional transfer squares of the 2 -> 2 channels a^ b -> c^ of classes.wl (first parton polarized, observed parton c, u = -s-t), divided by gs^4, with the reference vectors pa -> pb, pb -> pa, p -> pa, k1 -> pa for the gluons.  All four components are formed; XN and NX must vanish.*)


(* ::Input:: *)
bornReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa|>;
channelPartons[name_] := StringSplit[StringReplace[name, "->" -> ""]];                                (* {a, b, c} *)
bornTableRaw = Association @ KeyValueMap[Function[{name, amps}, name -> Module[{sides, in = channelPartons[name][[1 ;; 2]]},
    bornKinematics;
    sides = sidesOf[amps, bornReferences, bornConservation];
    Association @ Table[spinCase -> With[{sq = square["table-" <> StringReplace[name, {" " -> "_", "->" -> "to"}], sides, "Born", bornReferences, spinCase]},
      bornKinematics; Factor[colourFactors[Together[average[in] normalize[sq, spinCase]/gs^4]] /. u -> -s - t]], {spinCase, spinCases}]]],
  amplitudes["BornTable"]];
Require[FreeQ[bornTableRaw, Complex | Pair | Momentum | LorentzIndex | nn | eX | eC | nl], "A transfer Born square is not a real scalar", bornTableRaw];
Require[And @@ Flatten[Table[bornTableRaw[name][mixed] === 0, {name, Keys[bornTableRaw]}, {mixed, {"XN", "NX"}}]],
  "A mixed component XN or NX of a transfer Born square does not vanish", Map[KeyTake[#, {"XN", "NX"}] &, bornTableRaw]];
bornTable = Association @ Table[spinCase -> Map[#[spinCase] &, bornTableRaw], {spinCase, outputSpins}];
Dataset[Map[Simplify[# /. D -> 4] &, bornTable, {2}]]


(* ::Text:: *)
(*Checks of the table.  (1) In four dimensions the two components coincide in every channel (the evanescent terms differ).  (2) The channels whose leading-order transverse spin transfer is in the literature: the Collins hard factors of F. Yuan, Phys. Rev. Lett. 100 (2008) 032003 [arXiv:0709.3272], Eq. (6), transcribed verbatim in reference/Yuan-2008-Collins-hard-factors.wl, "exactly the same as those calculated for the transverse spin transfer in the hard partonic processes" (Stratmann and Vogelsang, Phys. Lett. B295 (1992) 277; Collins, Heppelmann and Ladinsky, Nucl. Phys. B420 (1994) 565).  They are written as H = alphaS^2 Pi/s^2 x h(s,t,u); the unpolarized hard factor in the same normalization is dsigma/dt = alphaS^2 Pi/s^2 x |M|^2/gs^4 (averaged), so h is to be compared with our table directly: the first channel, q qp -> q = -K s u/t^2 with K = (Nc^2-1)/Nc^2, fixes the normalization (no free factor), and the other channels must then agree.  The q qb channels (q qb -> q, q qb -> qb) are not listed by Yuan; they are checked against charge conjugation (qb q -> q is the conjugate of q qb -> qb) and, when Processes/Collins has been run, against its independently written resolved-photon Born table (6-ResolvedBorn.wl, Results/ResolvedBorn.wl), which was compared with the paper's table.*)


(* ::Input:: *)
fourDimensional[expr_] := Simplify[expr /. D -> 4 /. CA -> Nc /. CF -> (Nc^2 - 1)/(2 Nc) /. u -> -s - t];
Require[And @@ Table[fourDimensional[bornTable["XX"][name] - bornTable["NN"][name]] === 0, {name, Keys[bornTableRaw]}],
  "The two transverse components of a transfer Born square differ in four dimensions"];
Get[FileNameJoin[{NotebookDirectory[], "reference", "Yuan-2008-Collins-hard-factors.wl"}]];
literature = YuanCollinsHardFactors[s, t, u, Nc];                         (* h(s,t,u) for the channels of Eq. (6) *)
literatureDifferences = Association @ KeyValueMap[#1 -> fourDimensional[bornTable["XX"][#1] - #2] &, literature];
Require[And @@ (# === 0 & /@ Values[literatureDifferences]), "A transfer Born square differs from the literature", literatureDifferences];
Require[fourDimensional[bornTable["XX"]["qb q -> q"] - bornTable["XX"]["q qb -> qb"]] === 0 && Simplify[bornTable["NN"]["qb q -> q"] - bornTable["NN"]["q qb -> qb"]] === 0,
  "The charge-conjugate channels qb q -> q and q qb -> qb differ"];
collinsFile = FileNameJoin[{NotebookDirectory[], "..", "Collins", "Results", "ResolvedBorn.wl"}];
collinsDifferences = If[FileExistsQ[collinsFile], With[{resolved = Get[collinsFile]["Transfer"]},
    Association @ Table[name -> fourDimensional[bornTable["XX"][name] - resolved[name]["Collins"]], {name, Intersection[Keys[bornTableRaw], Keys[resolved]]}]],
  "Processes/Collins/Results/ResolvedBorn.wl not found: comparison skipped"];
Require[!AssociationQ[collinsDifferences] || And @@ (# === 0 & /@ Values[collinsDifferences]), "A transfer Born square differs from Processes/Collins' resolved Born table", collinsDifferences];
{literatureDifferences, collinsDifferences}


(* ::Section:: *)
(*Born*)


(* ::Text:: *)
(*The Born square of the class, in both components, with the same reference vectors; it must equal the table's entry for the class.*)


(* ::Input:: *)
bornKinematics;
bornSides = If[amplitudes["Born"] === {}, None, sidesOf[amplitudes["Born"], bornReferences, bornConservation]];
bornAll = If[bornSides === None, Association @ Table[spinCase -> 0, {spinCase, spinCases}],
  Association @ Table[spinCase -> With[{sq = square["born", bornSides, "Born", bornReferences, spinCase]},
    bornKinematics; Factor[colourFactors[Together[average[incoming] normalize[sq, spinCase]]] /. u -> -s - t]], {spinCase, spinCases}]];
Require[FreeQ[bornAll, Complex | Pair | Momentum | LorentzIndex | nn | eX | eC | nl], "The Born square is not a real scalar", bornAll];
Require[bornAll["XN"] === 0 && bornAll["NX"] === 0, "A mixed component XN or NX of the Born square does not vanish", bornAll];
born = KeyTake[bornAll, outputSpins];
If[spec["Born"] =!= None,
  Require[And @@ Table[Simplify[born[sc]/gs^4 - bornTable[sc][StringRiffle[spec["Born"][[1 ;; 2]], " "] <> " -> " <> spec["Born"][[3]]]] === 0, {sc, outputSpins}],
    "The Born square differs from the transfer Born table", born]];
born


(* ::Section:: *)
(*Virtual: the one-loop interference in terms of scalar integrals*)


(* ::Text:: *)
(*The one-loop amplitude is contracted with the conjugate Born amplitude and the spin sums as above; the loop momentum then appears in scalar products with the external momenta and in loop^2, which are combinations of the inverse propagators of one of the families*)
(**)
(*   VirtualST[n1,n2,n3,n4]: 1/[(loop^2)^n1 ((loop+pa)^2)^n2 ((loop+pa+pb)^2)^n3 ((loop+p)^2)^n4]   (external legs pa, pb, -pd, -p in cyclic order: invariants s and t),*)
(*   VirtualSU: ... ((loop+pa+pb)^2)^n3 ((loop+pa+pb-p)^2)^n4   (legs pa, pb, -p, -pd: s and u),      VirtualTU: ((loop+pa)^2)^n2 ((loop+pa-p)^2)^n3 ((loop+pa-p+pb)^2)^n4   (legs pa, -p, pb, -pd: t and u),*)
(**)
(*with negative indices for numerators (Processes/Hadronic), and in two more products: the hatted loopHat2 = loophat.loophat of the HVBM scheme, and n.loop = nl from the determinants of the Levi-Civita tensors.  FeynArts routes the loop momentum differently in every graph; for each graph the routing loop -> +-loop + shift that maps its propagators onto one family is found by trying the shifts.  Kira reduces the families in step 3.*)


(* ::Input:: *)
families = <|
  "VirtualST" -> {loop, loop + pa, loop + pa + pb, loop + p},                 (* legs pa, pb, -pd, -p in cyclic order: invariants s and t *)
  "VirtualSU" -> {loop, loop + pa, loop + pa + pb, loop + pa + pb - p},       (* legs pa, pb, -p, -pd: s and u *)
  "VirtualTU" -> {loop, loop + pa, loop + pa - p, loop + pa - p + pb}|>;      (* legs pa, -p, pb, -pd: t and u *)
familyHeads = <|"VirtualST" -> KVirtualST, "VirtualSU" -> KVirtualSU, "VirtualTU" -> KVirtualTU|>;
shifts = Flatten[Table[c1 pa + c2 pb + c3 p, {c1, -1, 1}, {c2, -1, 1}, {c3, -1, 1}]];
propagatorMomenta[amp_] := Union[Cases[amp, PropagatorDenominator[Momentum[q_, ___], 0] /; !FreeQ[q, loop] :> Expand[q /. bornConservation], Infinity]];   (* the loop propagators; tree propagators stay as they are *)
sameUpToSign[q1_, q2_] := Expand[q1 - q2] === 0 || Expand[q1 + q2] === 0;
(* find sigma, delta and the family such that every propagator of the graph, with loop -> sigma loop + delta, is one of the family's *)
routing[amp_] := Module[{props = propagatorMomenta[amp], found},
  found = Catch[Do[If[And @@ Table[Or @@ Table[sameUpToSign[q /. loop -> sigma loop + delta, f], {f, families[name]}], {q, props}],
      Throw[{name, sigma, delta}]], {name, Keys[families]}, {sigma, {1, -1}}, {delta, shifts}]; None];
  Require[found =!= None, "No family fits the propagators of a one-loop graph", props]; found];
routings = routing /@ amplitudes["Virtual"];
Counts[routings[[All, 1]]]


(* ::Text:: *)
(*The normal component of the loop momentum.  The integrand of a one-loop integral depends on the loop momentum through loop^2, loop.pa, loop.pb, loop.p and through its components orthogonal to the scattering plane: n.loop and the hatted part.  The orthogonal space has D - 3 dimensions (the direction n and the D - 4 extra dimensions, all spacelike), and the integration over its orientations is symmetric, so that with m = D - 3 and the orthogonal square loopPerp^2 = loop^2 - loopPar^2, loopPar^2 = d.G^-1.d, d = (loop.pa, loop.pb, loop.p),*)
(**)
(*   (n.loop)^(2a) loopHat2^b  ->  (-1)^a (loopPerp^2)^(a+b) (1/2)_a ((m-1)/2)_b / (m/2)_(a+b)     (odd powers of n.loop vanish),*)
(**)
(*the moments of a direction on the unit sphere of m dimensions (x_n^2 is Beta(1/2, (m-1)/2) distributed): for a = 1, b = 0 this is -loopPerp^2/(D-3), for a = 1, b = 1 it is -(D-4) (loopPerp^2)^2/((D-3)(D-1)), the contraction of the symmetric rank-four tensor with n n and the hatted metric.  Terms without n.loop keep their loopHat2 power as HatPower[b] for the dimension shift of step 3, exactly as in Processes/Hadronic (for them the two treatments are equivalent, and the dimension shift is the one that was validated against the reference codes there).*)


(* ::Input:: *)
loopProducts = {L2 -> Pair[Momentum[loop, D], Momentum[loop, D]], Lpa -> Pair[Momentum[loop, D], Momentum[pa, D]],
  Lpb -> Pair[Momentum[loop, D], Momentum[pb, D]], Lp -> Pair[Momentum[loop, D], Momentum[p, D]]};       (* the four scalar products of the loop momentum (FeynCalc's internal form of SPD) *)
inverseProps[name_] := Module[{q = families[name], dj, sol},
  dj = Array[Dj, 4];
  sol = Solve[Table[dj[[j]] == ExpandScalarProduct[SPD[q[[j]], q[[j]]]] /. Reverse /@ loopProducts, {j, 4}] /. u -> -s - t, {L2, Lpa, Lpb, Lp}];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators"]; First[sol]];
loopRules[name_] := loopRules[name] = With[{sol = inverseProps[name]}, Table[(x /. loopProducts) -> (x /. sol), {x, {L2, Lpa, Lpb, Lp}}]];
(* the orthogonal square of the loop momentum in terms of its four scalar products (Gram matrix of pa, pb, p with u = -s-t) *)
loopPerp2 = (L2 - Factor[{Lpa, Lpb, Lp} . Inverse[gram /. u -> -s - t] . {Lpa, Lpb, Lp}]) /. loopProducts;
perpMoment[a_, b_] := (-1)^a Pochhammer[1/2, a] Pochhammer[(D - 4)/2, b]/Pochhammer[(D - 3)/2, a + b];      (* m = D - 3 *)
normalLoopReduction[expr_] := Module[{terms = DeleteCases[termList[Expand[expr]], 0]},
  Total[Table[With[{a = Exponent[term, nl], b = Exponent[term, loopHat2]},
    Which[OddQ[a], 0, a == 0, term, True, (term/(nl^a loopHat2^b)) perpMoment[a/2, b] loopPerp2^(a/2 + b)]], {term, terms}]]];
(* the shifted graph with its loop propagators as inert symbols DjOf[q] *)
shiftedGraph[amp_, {name_, sigma_, delta_}] := Module[{shifted},
  shifted = amp /. Momentum[q_, dim_] /; !FreeQ[q, loop] :> Momentum[Expand[q /. loop -> sigma loop + delta], dim] /.
    PropagatorDenominator[c_. Momentum[q_, dim_], 0] :> PropagatorDenominator[Momentum[Expand[c q /. bornConservation], dim], 0];    (* c: FeynCalc writes Momentum[-loop] as -Momentum[loop] *)
  shifted /. FeynAmpDenominator[ds__] :> Times @@ (Replace[#, {PropagatorDenominator[c_. Momentum[q_, ___], 0] /; FreeQ[q, loop] :> 1/ExpandScalarProduct[SPD[q, q]],
    PropagatorDenominator[c_. Momentum[q_, ___], 0] :> 1/DjOf[c q]}] & /@ {ds})];
(* the contracted interference of one graph -> sum of coefficient x family[indices] x HatPower[n] *)
toFamilies[expr0_, name_] := Module[{expr, fourDimensionalProducts},
  (* the four-dimensional products of the HVBM scheme equal the D-dimensional ones for the four-dimensional external momenta, and loop4.loop4 = loop.loop - loophat.loophat *)
  fourDimensionalProducts = {Pair[Momentum[loop], Momentum[loop]] -> Pair[Momentum[loop, D], Momentum[loop, D]] - loopHat2,
    Pair[Momentum[loop], Momentum[x_]] :> Pair[Momentum[loop, D], Momentum[x, D]], Pair[Momentum[x_], Momentum[loop]] :> Pair[Momentum[x, D], Momentum[loop, D]]};
  (* the spin sum of an outgoing quark k1 gives k1-slash: loop.k1 by momentum conservation *)
  expr = ExpandScalarProduct[expr0] /. fourDimensionalProducts /. Pair[Momentum[k1, D], Momentum[loop, D]] -> ExpandScalarProduct[SPD[pa + pb - p, loop]] /.
    Pair[Momentum[loop, D - 4], Momentum[loop, D - 4]] -> loopHat2;
  expr = normalLoopReduction[expr] /. loopRules[name];
  Require[FreeQ[expr, nl], "The normal component of the loop momentum was not reduced"];
  Require[FreeQ[expr, Momentum[loop, ___]], "A scalar product of the loop momentum was not expressed through the propagators", Short[Union[Cases[expr, h_[a___] /; !FreeQ[{a}, Momentum[loop, ___]] && FreeQ[h, Times | Plus | Power], Infinity]], 8]];
  expr = Expand[expr /. DjOf[q_] :> Dj @@ (Position[families[name], f_ /; sameUpToSign[f, q]][[1]])];
  (* every term is a coefficient times powers of the Dj (negative for propagators, positive for numerators) and of loopHat2 *)
  Total[Table[With[{powers = Table[Exponent[term, Dj[j]], {j, 4}], hat = Exponent[term, loopHat2]},
    (term/(Times @@ Table[Dj[j]^powers[[j]], {j, 4}] loopHat2^hat)) familyHeads[name] @@ (-powers) HatPower[hat]], {term, DeleteCases[termList[expr], 0]}]]];   (* a graph can vanish exactly once eps(l).l = eps(l).r = 0 are imposed *)
bornKinematics;
virtualLeft = MapThread[colourGroups[prepare[shiftedGraph[#1, #2], bornReferences, bornConservation]] &, {amplitudes["Virtual"], routings}];
bornSumRight = If[bornSides === None, {}, conjugateGroups[Total[prepare[#, bornReferences, bornConservation] & /@ amplitudes["Born"]]]];
Length[virtualLeft]


(* ::Text:: *)
(*The interference 2 Re(M1 M0*).  FeynArts writes the one-loop amplitude with the measure d^D loop/(2 Pi)^D; the scalar integrals are defined with the measure d^D loop/(i Pi^(D/2)), so that M1 = i/(4 Pi)^(D/2) x (sum of coefficients times integrals).  The products (coefficient x M0*) formed by the spin sums are real, and*)
(**)
(*   2 Re(M1 M0*) = (2/(4 Pi)^(D/2)) Sum_j [i c_j M0*] Re(I_j),*)
(**)
(*with the real parts of the integrals in the physical region (s > 0) taken in step 3.  The association virtual below holds Sum_j [i c_j M0*] I_j for the two components, with the loop-momentum hat powers as HatPower[n]; the factor 2/(4 Pi)^(D/2) is included in step 4.  The mixed components are formed as well (before the family labelling, at an exact point) and must vanish.*)


(* ::Input:: *)
virtualJobs[spinCase_] := Table[<|"Name" -> "virtual-" <> spinCase <> "-" <> ToString[g], "Kinematics" -> "Born", "Left" -> virtualLeft[[g]],
  "Right" -> bornSumRight, "Spin" -> spinCase, "References" -> bornReferences|>, {g, Length[virtualLeft]}];
virtualAll = If[amplitudes["Virtual"] === {}, Association @ Table[spinCase -> 0, {spinCase, spinCases}],
  Association @ Table[spinCase -> Module[{results = runJobs[virtualJobs[spinCase]]},
    bornKinematics;
    Collect[I average[incoming] normalize[Total[MapThread[toFamilies[#1["Value"], #2[[1]]] &, {results, routings}]], spinCase] /. u -> -s - t,
      {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, colourFactors]], {spinCase, spinCases}]];
Require[FreeQ[virtualAll, Complex], "The virtual interference is not real", virtualAll];
Require[virtualAll["XN"] === 0 && virtualAll["NX"] === 0, "A mixed component XN or NX of the virtual interference does not vanish", KeyTake[virtualAll, {"XN", "NX"}]];
virtual = KeyTake[virtualAll, outputSpins];
{Length[virtual["XX"]], Length[virtual["NN"]], Union[Cases[virtual, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]] // Length, Union[Cases[virtual, _HatPower, Infinity]]}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*For every real-emission process of the class, the square of the sum of its diagrams with the four projections, in terms of the invariants s, t, u, A, B, C and the hatted invariant spe3 (the normal component n.k1 reduced in every product); the processes are added with their weights (classes.wl).  Reference vectors of the physical polarization sums: pa for the gluons p and k1, pb for k2 (and the other incoming momentum for an incoming gluon).  The mixed components must vanish; they are checked at exact rational points (the symbolic sum of thousands of terms is not brought to a common denominator).*)


(* ::Input:: *)
realReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa, k2 -> pb|>;
realKinematics;
realSides = Association @ KeyValueMap[#1 -> sidesOf[#2, realReferences, realConservation] &, amplitudes["Real"]];
realSquares = Association @ Table[spinCase -> Association @ KeyValueMap[
    #1 -> normalize[square["real-" <> StringReplace[#1, " " -> "_"], #2, "Real", realReferences, spinCase], spinCase] &, realSides], {spinCase, spinCases}];
realKinematics;
(* the sum is not brought to a common denominator (step 3 expands it into partial fractions anyway); the colour factors are polynomials in Nc = CA *)
realAll = Association @ Table[spinCase -> (average[incoming] Total[KeyValueMap[weights[#1] #2 &, realSquares[spinCase]]] /. Nc -> CA),
  {spinCase, spinCases}];
Require[FreeQ[realAll, Complex], "The real-emission square is not real"];
Require[FreeQ[realAll, Pair | Momentum | LorentzIndex | nn | eX | eC | nk1], "A scalar product or the normal component was left in the real-emission square", Union[Cases[realAll, _Pair, Infinity]]];
mixedPoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], CA -> 2 + 10 v[[6]], spe3 -> v[[7]]/5}], {3}];
mixedValues = Table[Together[realAll[mixed] /. C -> -(s + t + u) - A - B /. point], {mixed, {"XN", "NX"}}, {point, mixedPoints}];
Require[mixedValues === {{0, 0, 0}, {0, 0, 0}}, "A mixed component XN or NX of the real-emission square does not vanish", mixedValues];
real = KeyTake[realAll, outputSpins];
Require[Max[Exponent[#, spe3] & /@ Join[termList[Expand[real["XX"]]], termList[Expand[real["NN"]]]]] <= 1, "The real-emission square has spe3^2 terms"];
LeafCount /@ real


(* ::Text:: *)
(*Check of the gauge invariance of the real-emission square: for the first process of the class and the XX component, the physical polarization sums with other reference vectors (pb for p and k1, pa for k2) give the same result.  The amplitudes are prepared again, since the conditions eps(l).r = 0 depend on the reference vectors.*)


(* ::Input:: *)
alternativeReferences = <|pa -> pb, pb -> pa, p -> pb, k1 -> pb, k2 -> pa|>;
firstProcess = First[Keys[amplitudes["Real"]]];
realKinematics;
alternative = normalize[square["gauge-" <> StringReplace[firstProcess, " " -> "_"],
  sidesOf[amplitudes["Real"][firstProcess], alternativeReferences, realConservation], "Real", alternativeReferences, "XX"], "XX"];
(* the difference is evaluated exactly at random rational points (all invariants, D and Nc), with C = -s4 - A - B *)
gaugePoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], Nc -> 2 + 10 v[[6]], spe3 -> v[[7]]/5}], {3}];
referenceDifference = Table[Together[(alternative - realSquares["XX"][firstProcess]) /. C -> -(s + t + u) - A - B /. point], {point, gaugePoints}];
Require[referenceDifference === {0, 0, 0}, "The real-emission square depends on the gluon reference vectors", referenceDifference];
"gauge invariance of the real emission checked"


(* ::Section:: *)
(*Timing and save*)


(* ::Text:: *)
(*The number of products computed, their total time and the longest one, in seconds (read from the product files).*)


(* ::Input:: *)
pairTimes = Association @ Table[FileBaseName[f] -> Get[f]["Seconds"], {f, FileNames["*.wl", pairDirectory]}];
workerReuseStatistics = If[nWorkers > 1,
  ParallelEvaluate[<|"Trace" -> traceAlgebraCache["Statistics"][],
    "Polarization" -> tensorAlgebraCache["Statistics"][], "MemoryBytes" -> MemoryInUse[]|>],
  {<|"Trace" -> traceAlgebraCache["Statistics"][],
    "Polarization" -> tensorAlgebraCache["Statistics"][], "MemoryBytes" -> MemoryInUse[]|>}];
SaveResult["AlgebraReuseMetrics", <|"Workers" -> workerReuseStatistics,
  "PairJobs" -> Length[pairTimes], "TotalWorkerSeconds" -> Total[pairTimes]|>, "FeynCalc/class-" <> class];
{Length[pairTimes], Round[Total[pairTimes]], Round[Max[pairTimes]]}


(* ::Text:: *)
(*Save the Born square, the virtual interference, the real-emission square and the transfer Born table for steps 3 and 4 (the format of Processes/Hadronic, with the components XX and NN in place of UU and LL), and close the parallel kernels.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Virtual" -> virtual, "Real" -> real, "BornTable" -> bornTable, "Class" -> class, "Spec" -> spec, "Spins" -> outputSpins,
  "Families" -> families, "Variables" -> "s, t, u = -s-t; A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 = -s4-A-B; spe3; loopHat2 via HatPower[n]; XX divided by |eX||eC| = -t u/(t+u)"|>;
SaveResult["Algebra", algebra, "FeynCalc/class-" <> class];
CloseKernels[];
Keys[algebra]
