(* ::Title:: *)
(*Transversely polarized proton-proton classes (A_TT), step 2 (FeynCalc route): spin and colour algebra*)


(* ::Text:: *)
(*Input: the amplitudes of step 1 for the class chosen below.  Output: Results/FeynCalc/class-<n>/Algebra.wl with*)
(*  Born:       the Born square for the two transverse spin projections defined below (TT and Symmetric), summed over final spins and colours, averaged over the initial colours, in D dimensions;*)
(*  Virtual:    the interference 2 Re(M_1loop M_Born^*) for the same projections, written in terms of scalar one-loop integrals (families VirtualST, VirtualSU, VirtualTU as in Processes/Hadronic) with D-dimensional coefficients;*)
(*  Real:       the weighted sum over the real-emission processes of the class of |M|^2 for the same projections, with the unobserved momenta k1, k2 explicit (step 3 integrates them);*)
(*  BornTable:  the D-dimensional transversely polarized Born squares |M|^2/gs^4 of every 2 -> 2 channel needed by the collinear counterterms of the five classes (q qb -> q' qb', q q -> q q, q qb -> q qb, q qb -> g g, with their a <-> b relabellings), in the format of Processes/PartonicBorn, checked against the lowest-order cross sections of Soffer, Stratmann and Vogelsang (reference/SSV-born.wl).*)
(*The squaring method, the colour treatment, the gluon polarization sums, the scalar-integral families and the parallel organization are those of Processes/Hadronic/2-Algebra-FeynCalc.wl (read that file first); only the spin projections of the two incoming partons differ.*)
(**)
(*Spin densities (docs/TRANSVERSE-PP.md).  Both incoming partons are quarks or antiquarks with a transverse spin: u ubar -> pa-slash gamma5 sa-slash/2 and v vbar -> the same expression (for a massless antiquark v(p,s) vbar(p,s) = p-slash (1 + gamma5 s-slash)/2 as for the quark); these are the densities whose partonic cross section multiplies the transversity distribution delta f = f(up,up) - f(down,up) of each proton, i.e. the partonic quantity computed here is (1/2)[sigma(up,up) - sigma(up,down)].  The spin vectors are four-dimensional, orthogonal to pa and pb, and are written in the basis of the scattering plane of pa, pb, p: the in-plane axis eX, along the transverse momentum of the observed parton,*)
(**)
(*   eX = Sqrt[s/(t u)] [ p + (u/s) pa + (t/s) pb ],   eX.pa = eX.pb = 0,  eX^2 = -1,  eX.p = -Sqrt[t u/s] = -pT,*)
(**)
(*and the normal n ~ eps^{mu pa pb p}, n^2 = -1, n.pa = n.pb = n.p = 0, so that sa = cos(phia) eX + sin(phia) n and sb = cos(phib) eX + sin(phib) n.  The squared matrix element is then T = XX cos(phia) cos(phib) + NN sin(phia) sin(phib) + XN cos(phia) sin(phib) + NX sin(phia) cos(phib), and the four components are computed separately (XX from the unnormalized axis aX = p + (u/s) pa + (t/s) pb, aX^2 = -t u/s, divided by -aX^2; XN and NX up to the positive factor Sqrt[s/(t u)]).  Parity (reflection through the scattering plane) demands XN = NX = 0: for the real emission these components are odd in the normal component of the unobserved momentum and vanish after the phase-space integration.*)
(**)
(*The observable.  With both proton spins along one direction x and the observed parton at azimuth Phi with respect to it, T = (XX+NN)/2 + (XX-NN)/2 cos(2 Phi).  The physical cross section has no Phi-independent part (the two helicity flips carry two units of angular momentum along the beam: XX + NN = 0 in four dimensions, which is checked below), and the hard part of A_TT is the coefficient of cos(2 Phi),*)
(**)
(*   delta sigma_TT = (XX - NN)/2   ("TT"),   which equals the projection with F(p, sa, sb) = (s/(t u)) [2 (p.sa)(p.sb) + (t u/s) sa.sb] of Mukherjee, Stratmann and Vogelsang (F = +1 for sa = sb = eX, -1 for sa = sb = n).*)
(**)
(*In the HVBM scheme the symmetric part (XX + NN)/2 ("Symmetric") does not vanish in D dimensions: the two hatted gamma matrices that connect the two densities give a term proportional to (D - 4) s^2 (sa.sb) at Born level, and hatted momenta give similar terms at NLO.  This evanescent, Phi-independent part is kept as a separate projection through step 4 (where its poles must still cancel) and never enters TT, which is why TT is defined as the antisymmetric combination and not as XX alone.*)
(**)
(*Normal components.  n.k1 = -n.k2 = nr appears in NN, XN and NX of the real emission; its odd powers integrate to zero and nr^2 is a rational function of the invariants: with k1 = k1par + c n + k1hat (parts in the plane of pa, pb, p, along n, and in the D-4 extra dimensions), n.k1 = -c, and k1^2 = 0 gives nr^2 = k1par^2 + spe3, k1par^2 = d.G^-1.d with G the Gram matrix of pa, pb, p and d = (pa.k1, pb.k1, p.k1) (Processes/Collins, normalReduction).  The loop momentum has a normal component too: the decomposition is the same, (n.loop)^2 = looppar^2 - loop^2 + loopHat2, with looppar^2 written through the propagators, and the odd powers vanish because the tensor loop integrals are built from pa, pb, p only; aX.loop is a combination of pa.loop, pb.loop, p.loop.*)
(**)
(*Invariants as in Processes/Hadronic: s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2; real emission: s4 = (k1+k2)^2 = s + t + u, A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 with A + B + C = -s4; spe3 = k1hat.k1hat.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*Choices: the class; the number of parallel kernels for the squares ("auto" = one per processor core, at most 8; "1" = none); and whether the product files of an earlier run of this same, unchanged file may be reused ("no" recomputes everything; "yes" only continues an interrupted run: a product file is reused only if its companion file <name>.input holds the fingerprint of the current inputs, the amplitudes of step 1 and this step file; a file written for other inputs, or without that companion (for example, after an interrupted write), is computed again).*)


(* ::Input:: *)
LoadFeynCalc[];
(* FeynCalc expresses this memoization threshold in MB; transient algebra can use more memory. *)
FeynCalc`$FCMemoryAvailable = Min[FeynCalc`$FCMemoryAvailable, Floor[1024^3/10^6]];
class = "5";
workers = "auto";
resume = "no";
loaded = LoadAmplitudes["class-" <> class];
amplitudes = loaded["FeynCalc"]; spec = loaded["Spec"]; weights = loaded["Weights"];
spec


(* ::Section:: *)
(*Conventions*)


(* ::Text:: *)
(*The conventions, the kinematics and the spin axes are installed by functions, because every parallel kernel installs them too.  The axes aX and nn are vectors of their own whose scalar products with every momentum are set from the definitions above (defineAxes), as in Processes/Collins; nn.k1 = nr and nn.k2 = -nr in the real kinematics, nn.k1 = 0 in the Born kinematics (k1 = pa + pb - p), and nn.loop stays symbolic until the virtual section reduces it.  One unobserved momentum is not independent: inside the amplitudes k2 is replaced by pa + pb - p - k1 (real emission) and k1 by pa + pb - p (Born and one loop), except in the external spinors, whose spin sums give the slashed momentum itself.*)


(* ::Input:: *)
conventions := (FCSetDiracGammaScheme["BMHV"]; $LeviCivitaSign = -1);
conventions;
axisCoefficients = {u/s, t/s, 1};                      (* aX = (u/s) pa + (t/s) pb + p *)
defineAxes[momenta_List, normalComponents_List] := Module[{base = {pa, pb, p}, ca = axisCoefficients},
  Scan[(Momentum[#, D - 4] = 0) &, {aX, nn}];
  Do[SP[aX, x] = SPD[aX, x] = Factor[Sum[ca[[i]] SP[base[[i]], x], {i, 3}]], {x, momenta}];
  SP[aX] = SPD[aX] = Factor[Sum[ca[[i]] ca[[j]] SP[base[[i]], base[[j]]], {i, 3}, {j, 3}]];
  SP[nn] = SPD[nn] = -1; SP[aX, nn] = SPD[aX, nn] = 0;
  Scan[(SP[nn, #] = SPD[nn, #] = 0) &, {pa, pb, p}];
  Do[SP[nn, rule[[1]]] = SPD[nn, rule[[1]]] = rule[[2]], {rule, normalComponents}]];
bornKinematics := (FCClearScalarProducts[]; SetMandelstam[s, t, u, pa, pb, -p, -k1, 0, 0, 0, 0];
  Scan[(Momentum[#, D - 4] = 0) &, {pa, pb, p, k1}];
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
  defineAxes[{pa, pb, p, k1, k2}, {k1 -> nr, k2 -> -nr}]);
bornConservation = k1 -> pa + pb - p;
realConservation = k2 -> pa + pb - p - k1;
incoming = spec["Born"] /. None -> spec["Real"][[1, 1]];       (* the two incoming partons *)
incoming = incoming[[1 ;; 2]];
Require[FreeQ[incoming, "g"], "The transversely polarized process needs two incoming quarks or antiquarks", incoming];
(* the spin axes of pa and pb for each component, and the normalization that turns the axis aX into the unit vector eX (XN, NX: up to Sqrt[s/(t u)]) *)
components = {"XX", "NN", "XN", "NX"};
spinAxes = <|"XX" -> <|pa -> aX, pb -> aX|>, "NN" -> <|pa -> nn, pb -> nn|>, "XN" -> <|pa -> aX, pb -> nn|>, "NX" -> <|pa -> nn, pb -> aX|>|>;
normalization[spinCase_] := If[spinCase === "XX", s/(t u), 1];
average[spinCase_] := 1/(2 Nc)^2;                                   (* 1/(2 Nc) per incoming quark: the 1/2 of the density (1 +- gamma5 s-slash)/2 and the colour average *)
colourFactors[expr_] := Factor[expr /. SUNN -> Nc] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
bornKinematics;
Require[Simplify[{SP[aX, pa], SP[aX, pb], SP[aX] + t u/s, SP[aX, p] + t u/s, SP[aX, k1] - t u/s, SP[nn, aX], SP[nn] + 1} /. u -> -s - t] === {0, 0, 0, 0, 0, 0, 0},
  "The in-plane axis is not what the text says"];
{incoming, average["XX"], SP[aX]}


(* ::Section:: *)
(*The squaring method*)


(* ::Text:: *)
(*As in Processes/Hadronic: (1) every amplitude is simplified once (Lorentz contraction, Dirac algebra on the open chains, one independent unobserved momentum, transversality eps(l).l = eps(l).r = 0 of every gluon with its reference vector r); (2) colour is separated and each colour network evaluated once; (3) only the products i <= j are formed; (4) the gluon polarization sums are applied one gluon at a time on grouped tensors; (5) the products are distributed over parallel kernels, one file per product.  The only change is in the fermion spin sums (step 4 of the method): after FermionSpinSum, the slashed momenta of the two incoming partons are replaced by the densities pa-slash gamma5 sa-slash and pb-slash gamma5 sb-slash with the axes of the component being computed (the factor 1/2 of each density is in average).*)


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

(* ---- 4. spin sums and gluon polarization sums ---- *)
polarizationSum[l_, mu_, nu_, r_] := -MTD[mu, nu] + (FVD[l, mu] FVD[r, nu] + FVD[r, mu] FVD[l, nu])/SPD[l, r];
tag[l_] := Symbol[ToString[l] <> "Spin"];
indexNames[l_] := {Symbol["mu" <> ToString[l]], Symbol["nu" <> ToString[l]]};
tensorPart[term_, {mu_, nu_}] := Times @@ Select[factorList[term], !FreeQ[#, LorentzIndex[mu | nu, ___]] &];
(* the contraction uses the installed scalar products, which differ between the Born and the real kinematics: the kinematics are part of the key *)
(* Each worker reuses identical contractions across diagram pairs. The explicit
   spin projections and the original BMHV trace evaluation are unchanged. *)
currentKinematics = None;
algebraCacheScope[] := {"HadronicTT", class, currentKinematics, "BMHV", $LeviCivitaSign};
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
(* the colour-weighted product of the Lorentz parts of one pair: spin sums with the transverse densities of pa and pb, traces, polarization sums
   (every external gluon of these classes is outgoing) *)
spinSum[product_, spinCase_, references_, gluonList_] := Module[{expr = product, axes = spinAxes[spinCase]},
  expr = expr /. (With[{l = #, lTag = tag[#]}, Spinor[c_. Momentum[l, D], m_, o_] :> Spinor[c Momentum[lTag, D], m, o]] & /@ {pa, pb});
  expr = expr /. Flatten[Table[{Momentum[Polarization[l, I, ___], D] -> LorentzIndex[indexNames[l][[1]], D],
    Momentum[Polarization[l, -I, ___], D] -> LorentzIndex[indexNames[l][[2]], D]}, {l, gluonList}]];
  expr = FermionSpinSum[expr];
  expr = expr /. (With[{l = #, lTag = tag[#], ax = axes[#]}, DiracGamma[Momentum[lTag, D], D] :> GSD[l] . GA[5] . GS[ax]] & /@ {pa, pb});
  If[!FreeQ[expr, tag[pa] | tag[pb]], Return[Failure["Density", <|"Message" -> "a tagged spinor momentum survived the density insertion"|>]]];
  (* Contract the same gluon polarization sums while their Lorentz indices
     are still carried by compact traces; avoid a large expanded tensor. *)
  Do[expr = project[expr, l, polarizationSum[l, Sequence @@ indexNames[l], references[l]]], {l, gluonList}];
  expr = groupedTrace[expr];
  If[!FreeQ[expr, _Failure], Return[FirstCase[expr, _Failure, Failure["Trace", <||>], Infinity]]];
  ExpandScalarProduct[Contract[expr]]];

(* ---- 3. one product (i, j) of colour-separated amplitudes: left = {{c_x, L_x}, ...}, right = {{c_y^*, L_y^*}, ...} ---- *)
interference[left_, right_, spinCase_, references_] := Module[{product},
  product = Total[Flatten[Table[colourValue[x[[1]] y[[1]]] x[[2]] y[[2]], {x, left}, {y, right}]]];
  If[!FreeQ[product, Failure], FirstCase[product, _Failure, Null, Infinity],
    spinSum[product, spinCase, references, gluonsOf[{left, right}]]]];
conjugateGroups[amp_] := colourGroups[ComplexConjugate[amp, FCRenameDummyIndices -> True]];
(* the real part: after expanding, every term is a number times real symbols (the denominators are real propagators), so the real part of the number is kept *)
realPartOf[expr_] := Expand[expr] /. Complex[a_, b_] :> a;

(* ---- the normal component of an unobserved or loop momentum x: polynomial in nx = n.x; odd powers integrate to zero, nx^2 = xpar^2 - x4.x4
   with xpar^2 = d.G^-1.d (d the products of x with pa, pb, p; G their Gram matrix) and x4.x4 the four-dimensional square of x ---- *)
gramMatrix[products_] := Table[products[x, y], {x, {pa, pb, p}}, {y, {pa, pb, p}}];
parallelSquare[d_List, gram_] := Factor[d . Inverse[gram] . d];
normalReduction[expr_, nx_, parallel_, fourSquare_] := Module[{coefficients = CoefficientList[Expand[expr], nx]},
  If[coefficients === {}, coefficients = {0}];                                   (* a vanishing expression (a graph can vanish exactly) *)
  Require[Length[coefficients] <= 3, "A power of the normal component above two appeared", Length[coefficients] - 1];
  <|"Even" -> coefficients[[1]] + If[Length[coefficients] == 3, coefficients[[3]] (parallel - fourSquare), 0],
    "Odd" -> If[Length[coefficients] >= 2, coefficients[[2]], 0]|>];
"squaring functions defined"


(* ::Text:: *)
(*Parallel kernels: each loads FeynFacet and FeynCalc, installs the conventions and receives the functions above; each job installs its kinematics (Born or real).  Every job (one product of the squaring method, one spin component) is written to Results/FeynCalc/class-<n>/pairs/<name>.wl with its time.  A diagonal product |M_i|^2 is real; its epsilon-tensor terms may leave an imaginary part that vanishes only after simplification, which is checked to be exactly zero at a rational point and dropped.*)


(* ::Input:: *)
feynFacetFile = FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}];
nWorkers = If[workers === "auto", Min[$ProcessorCount, 8], ToExpression[workers]];   (* "auto": one per processor core, at most 8 *)
If[nWorkers > 1,
  LaunchKernels[nWorkers];
  With[{file = feynFacetFile}, ParallelEvaluate[Get[file]; LoadFeynCalc[];
    FeynCalc`$FCMemoryAvailable = Min[FeynCalc`$FCMemoryAvailable, Floor[1024^3/10^6]]]];
  DistributeDefinitions[conventions, axisCoefficients, defineAxes, bornKinematics, realKinematics, gluonsOf, transversality, eliminatedTransversality, conserve, prepare,
    colourHead, factorList, termList, colourPart, colourGroups, colourValue, polarizationSum, tag, indexNames, tensorPart,
    projectOnce, project, spinSum, interference, conjugateGroups, realPartOf, spinAxes,
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
checkPoint = {s -> 1, t -> -3/10, u -> -2/5, A -> -1/7, B -> -2/9, C -> -3/10 + 1/7 + 2/9, D -> 37/10, Nc -> 3, spe3 -> 1/11, nr -> 3/13, gs -> 1, Nf -> 5};
imaginaryPartOf[expr_] := (expr - (expr /. Complex[a_, b_] :> Complex[a, -b]))/2;
square[name_, sides_, kinematics_, references_, spinCase_] := Module[{jobs = squareJobs[name, sides, kinematics, references, spinCase], results},
  results = runJobs[jobs];
  Require[And @@ MapThread[#2["Pair"][[1]] != #2["Pair"][[2]] || (imaginaryPartOf[#1["Value"]] /. checkPoint) === 0 &, {results, jobs}],
    "A diagonal product has a nonzero imaginary part", name];
  Total[MapThread[If[#2["Pair"][[1]] == #2["Pair"][[2]], realPartOf[#1["Value"]], 2 realPartOf[#1["Value"]]] &, {results, jobs}]]];
sidesOf[amps_, references_, rule_] := With[{prepared = prepare[#, references, rule] & /@ amps}, {colourGroups /@ prepared, conjugateGroups /@ prepared}];
"parallel kernels: " <> ToString[Max[1, Length[Kernels[]]]]


(* ::Section:: *)
(*Born*)


(* ::Text:: *)
(*The Born square with the reference vectors p -> pa, k1 -> pa for the gluons, for the four components; then the checks XN = NX = 0 (exactly) and XX + NN = 0 in four dimensions, and the two projections TT = (XX - NN)/2 and Symmetric = (XX + NN)/2.  The evanescent symmetric part is displayed.*)


(* ::Input:: *)
bornReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa|>;
bornKinematics;
bornSides = If[amplitudes["Born"] === {}, None, sidesOf[amplitudes["Born"], bornReferences, bornConservation]];
bornComponents = If[bornSides === None, Association[# -> 0 & /@ components],
  Association @ Table[sc -> With[{sq = square["born", bornSides, "Born", bornReferences, sc]},
    bornKinematics; Factor[colourFactors[Together[average[sc] normalization[sc] sq]] /. u -> -s - t]], {sc, components}]];
Require[FreeQ[bornComponents, Complex | Pair | Momentum | LorentzIndex | Eps], "The Born square is not a real scalar", bornComponents];
Require[bornComponents["XN"] === 0 && bornComponents["NX"] === 0, "The mixed components of the Born square do not vanish", bornComponents];
Require[Factor[(bornComponents["XX"] + bornComponents["NN"]) /. D -> 4] === 0, "XX + NN of the Born square does not vanish in four dimensions", bornComponents];
born = <|"TT" -> Factor[(bornComponents["XX"] - bornComponents["NN"])/2], "Symmetric" -> Factor[(bornComponents["XX"] + bornComponents["NN"])/2],
  "XX" -> bornComponents["XX"], "NN" -> bornComponents["NN"]|>;
born


(* ::Section:: *)
(*The transversely polarized Born table*)


(* ::Text:: *)
(*The D-dimensional Born squares |M|^2/gs^4 (colours and spins of the unobserved parton summed, colours of a and b averaged, the two transverse densities inserted) of every 2 -> 2 channel a b -> c that the collinear counterterms of the five classes need: q qb -> q' (and qb'), q q -> q, q qb -> q (and qb), q qb -> g, their a <-> b relabellings (qb q -> ...), and the charge conjugate qb qb -> qb; plus the control channels q q' -> q, q q' -> q', q qb' -> q, which must vanish (classes 1-4: the two polarized partons on different fermion lines).  They are computed directly, as in Processes/Collins/6-ResolvedBorn.wl (FeynCalc's generic spin sum with the densities inserted, the polarization sums with reference vector pa, SUNSimplify for the colour), independently of the parallel squaring machinery above, and the class's own Born square must agree with its table entry.  Keys as in Processes/PartonicBorn ("q qb -> qp"); each entry holds XX, NN, TT and Symmetric.*)
(**)
(*Checks: XN = NX = 0 and XX + NN = 0 in four dimensions for every channel; the a <-> b relabelled channel equals the original with t <-> u (the axis aX is symmetric, the normal changes sign for both densities); the charge conjugate equals its representative; the control channels vanish; and at D = 4, Nc = 3 the TT entries reproduce Table I of Soffer, Stratmann and Vogelsang: with sigma = s E dsigma/d^3p = 2 s^2 dsigma/(dt du dPhi) and dsigma/(dt du dPhi) = dsigma/(dt dPhi) delta(s4), the cos(2 Phi) coefficient of their Eq. (6) corresponds to |M|^2_TT/gs^4 = (4/9) (t u/s^2) delta|M|^2 (the hard part is alphaS^2 |M|^2_TT/(gs^4 s v) delta(1-w), Processes/Hadronic/4-Unpolarized.wl).*)


(* ::Input:: *)
tableKinematics := (FCClearScalarProducts[]; SetMandelstam[s, t, u, pa, pb, -p, -pd, 0, 0, 0, 0];
  Scan[(Momentum[#, D - 4] = 0) &, {pa, pb, p, pd}]; defineAxes[{pa, pb, p, pd}, {pd -> 0}]);
tableChannels = {{"q", "qb", "qp"}, {"q", "qb", "qpb"}, {"q", "q", "q"}, {"q", "qb", "q"}, {"q", "qb", "qb"}, {"q", "qb", "g"},
  {"qb", "q", "qp"}, {"qb", "q", "qpb"}, {"qb", "q", "q"}, {"qb", "q", "qb"}, {"qb", "q", "g"}, {"qb", "qb", "qb"}};
controlChannels = {{"q", "qp", "q"}, {"q", "qp", "qp"}, {"q", "qpb", "q"}, {"q", "qpb", "qpb"}};
field = <|"q" -> F[3, {1}], "qb" -> -F[3, {1}], "qp" -> F[3, {2}], "qpb" -> -F[3, {2}], "g" -> V[5]|>;
channelName[{a_, b_, c_}] := a <> " " <> b <> " -> " <> c;
flavourNumber[x_] := Switch[x, "q", {1, 0}, "qb", {-1, 0}, "qp", {0, 1}, "qpb", {0, -1}, "g", {0, 0}];
partner[{a_, b_, c_}] := Module[{ok = Select[Keys[field], Total[flavourNumber /@ {a, b}] === Total[flavourNumber /@ {c, #}] &]},
  Require[Length[ok] == 1, "The unobserved parton of " <> channelName[{a, b, c}] <> " is not unique", ok]; First[ok]];
tableDiagrams[ch_] := InsertFields[CreateTopologies[0, 2 -> 2], {field[ch[[1]]], field[ch[[2]]]} -> {field[ch[[3]]], field[partner[ch]]},
  InsertionLevel -> {Classes}, Model -> "SMQCD", ExcludeParticles -> {S[_], U[_], V[1 | 2 | 3], F[1 | 2 | 4], F[3, {3}]}];
tableAmplitude[ch_] := Total[FCFAConvert[CreateFeynAmp[tableDiagrams[ch], Truncated -> False, PreFactor -> 1],
  IncomingMomenta -> {pa, pb}, OutgoingMomenta -> {p, pd}, UndoChiralSplittings -> True, ChangeDimension -> D,
  List -> True, SMP -> False, Contract -> False, DropSumOver -> True,
  FinalSubstitutions -> {FCGV["MU"] -> 0, FCGV["MC"] -> 0}] /. SMP["g_s"] -> gs /. GaugeXi[_] -> 1];
tableReferences = <|p -> pa, pd -> pa|>;
tableSquare[amp_, spinCase_] := Module[{conj, expr, gluons, indices, axes = spinAxes[spinCase]},
  tableKinematics;
  gluons = Union[Cases[amp, Polarization[l_, __] :> l, Infinity]];
  indices = Association @ Table[gluons[[i]] -> {Symbol["mu" <> ToString[i]], Symbol["nu" <> ToString[i]]}, {i, Length[gluons]}];
  expr = amp ComplexConjugate[amp];
  expr = expr /. (With[{l = #, lTag = tag[#]}, Spinor[c_. Momentum[l, D], m_, o_] :> Spinor[c Momentum[lTag, D], m, o]] & /@ {pa, pb});
  expr = expr /. Flatten[Table[{Momentum[Polarization[l, I, ___], D] -> LorentzIndex[indices[l][[1]], D],
    Momentum[Polarization[l, -I, ___], D] -> LorentzIndex[indices[l][[2]], D]}, {l, gluons}]];
  expr = expr Times @@ Table[polarizationSum[l, indices[l][[1]], indices[l][[2]], tableReferences[l]], {l, gluons}];
  expr = FermionSpinSum[FeynAmpDenominatorExplicit[expr]];
  expr = expr /. (With[{l = #, lTag = tag[#], ax = axes[#]}, DiracGamma[Momentum[lTag, D], D] :> GSD[l] . GA[5] . GS[ax]] & /@ {pa, pb});
  Require[FreeQ[expr, tag[pa] | tag[pb]], "A tagged spinor momentum survived the density insertion (table)"];
  expr = DiracSimplify[expr, DiracTraceEvaluate -> True] // Contract;
  expr = SUNSimplify[expr, Explicit -> True, SUNNToCACF -> False] // Contract // ExpandScalarProduct;
  Factor[colourFactors[Together[average[spinCase] normalization[spinCase] expr/gs^4 /. u -> -s - t]]]];
tableEntry[ch_] := Module[{amp = tableAmplitude[ch], xx, nnv, xn, nx},
  {xx, nnv, xn, nx} = tableSquare[amp, #] & /@ components;
  Require[FreeQ[{xx, nnv, xn, nx}, Complex | Pair | Momentum | LorentzIndex | Eps | SUNIndex], "A table entry is not a real scalar", channelName[ch]];
  Require[xn === 0 && nx === 0, "The mixed components of a table entry do not vanish", {channelName[ch], xn, nx}];
  Require[Factor[(xx + nnv) /. D -> 4] === 0, "XX + NN of a table entry does not vanish in four dimensions", {channelName[ch], xx, nnv}];
  <|"XX" -> xx, "NN" -> nnv, "TT" -> Factor[(xx - nnv)/2], "Symmetric" -> Factor[(xx + nnv)/2]|>];
tableEntries = Association @ Table[channelName[ch] -> tableEntry[ch], {ch, Join[tableChannels, controlChannels]}];
bornTable = Association @ Table[key -> Map[#[key] &, tableEntries], {key, {"TT", "Symmetric", "XX", "NN"}}];
TableForm[Normal[Map[Simplify[# /. D -> 4] &, bornTable["TT"]]]]


(* ::Input:: *)
tableCheck[expr_] := Factor[expr /. u -> -s - t];
Require[And @@ Table[tableCheck[bornTable["TT"][channelName[ch]]] === 0 && tableCheck[bornTable["Symmetric"][channelName[ch]]] === 0, {ch, controlChannels}],
  "A control channel of the Born table does not vanish", Map[bornTable["TT"], channelName /@ controlChannels]];
(* a <-> b: pa <-> pb exchanges t and u, i.e. t -> -s - t in the entries (u is eliminated) *)
swapped[name_] := Module[{parts = StringSplit[name, " "]}, parts[[2]] <> " " <> parts[[1]] <> " -> " <> parts[[4]]];
relabelDifferences = Association @ Flatten[Table[{key, name} -> tableCheck[bornTable[key][swapped[name]] - (bornTable[key][name] /. t -> -s - t)],
    {key, {"TT", "Symmetric"}}, {name, channelName /@ Select[tableChannels, #[[1]] === "q" && #[[2]] === "qb" &]}]];
Require[And @@ (# === 0 & /@ Values[relabelDifferences]), "An a <-> b relabelled entry is not the t <-> u image of the original", Select[relabelDifferences, # =!= 0 &]];
Require[tableCheck[bornTable["TT"]["qb qb -> qb"] - bornTable["TT"]["q q -> q"]] === 0 && tableCheck[bornTable["Symmetric"]["qb qb -> qb"] - bornTable["Symmetric"]["q q -> q"]] === 0,
  "The charge-conjugate channel differs from its representative"];
(* the class's own Born square (the parallel machinery) against the table (the direct computation) *)
If[spec["Born"] =!= None,
  Require[And @@ Table[tableCheck[born[key]/gs^4 - bornTable[key][StringRiffle[spec["Born"][[1 ;; 2]], " "] <> " -> " <> spec["Born"][[3]]]] === 0, {key, {"TT", "Symmetric", "XX", "NN"}}],
    "The Born square of the class differs from its table entry", born]];
(* the literature: Table I of Soffer, Stratmann and Vogelsang at D = 4, Nc = 3 *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "SSV-born.wl"}]];
literatureDifferences = Association @ Table[name -> Factor[(bornTable["TT"][name] /. {D -> 4, CA -> 3, CF -> 4/3} /. u -> -s - t) - ((4/9) (t u/s^2) SSVMatrixElement[name] /. u -> -s - t)],
  {name, Keys[SSVMatrixElement]}];
Require[And @@ (# === 0 & /@ Values[literatureDifferences]), "A transversely polarized Born square differs from Soffer, Stratmann and Vogelsang", literatureDifferences];
"Born table checked: control channels vanish, a <-> b relabelling, charge conjugation, the class Born, and the literature at D = 4, Nc = 3"


(* ::Section:: *)
(*Virtual: the one-loop interference in terms of scalar integrals*)


(* ::Text:: *)
(*As in Processes/Hadronic: the loop propagators of every graph are mapped by a shift of the loop momentum onto one of the families VirtualST, VirtualSU, VirtualTU, the graph is contracted with the conjugate Born amplitude and the spin sums above, and the scalar products of the loop momentum are written through the inverse propagators Dj.  New here: the products aX.loop and nn.loop of the spin axes with the loop momentum.  aX.loop is the combination (u/s) pa.loop + (t/s) pb.loop + p.loop; nn.loop is reduced as explained in the introduction, (nn.loop)^2 = looppar^2 - loop^2 + loopHat2 with looppar^2 = d.G^-1.d, d = (loop.pa, loop.pb, loop.p), and its odd powers dropped (they integrate to zero).  The hatted loop.loop is kept as loopHat2 (HatPower[n]), shifting the dimension of the integral in step 3.*)


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


(* ::Input:: *)
loopProducts = {L2 -> Pair[Momentum[loop, D], Momentum[loop, D]], Lpa -> Pair[Momentum[loop, D], Momentum[pa, D]],
  Lpb -> Pair[Momentum[loop, D], Momentum[pb, D]], Lp -> Pair[Momentum[loop, D], Momentum[p, D]]};       (* the four scalar products of the loop momentum (FeynCalc's internal form of SPD) *)
inverseProps[name_] := Module[{q = families[name], dj, sol},
  dj = Array[Dj, 4];
  sol = Solve[Table[dj[[j]] == ExpandScalarProduct[SPD[q[[j]], q[[j]]]] /. Reverse /@ loopProducts, {j, 4}] /. u -> -s - t, {L2, Lpa, Lpb, Lp}];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators"]; First[sol]];
loopRules[name_] := loopRules[name] = With[{sol = inverseProps[name]}, Table[(x /. loopProducts) -> (x /. sol), {x, {L2, Lpa, Lpb, Lp}}]];
(* the shifted graph with its loop propagators as inert symbols DjOf[q] *)
shiftedGraph[amp_, {name_, sigma_, delta_}] := Module[{shifted},
  shifted = amp /. Momentum[q_, dim_] /; !FreeQ[q, loop] :> Momentum[Expand[q /. loop -> sigma loop + delta], dim] /.
    PropagatorDenominator[c_. Momentum[q_, dim_], 0] :> PropagatorDenominator[Momentum[Expand[c q /. bornConservation], dim], 0];    (* c: FeynCalc writes Momentum[-loop] as -Momentum[loop] *)
  shifted /. FeynAmpDenominator[ds__] :> Times @@ (Replace[#, {PropagatorDenominator[c_. Momentum[q_, ___], 0] /; FreeQ[q, loop] :> 1/ExpandScalarProduct[SPD[q, q]],
    PropagatorDenominator[c_. Momentum[q_, ___], 0] :> 1/DjOf[c q]}] & /@ {ds})];
(* the products of the spin axes with the loop momentum *)
axisLoopRules = {
  Pair[Momentum[aX, ___], Momentum[loop, ___]] :> axisCoefficients . {Pair[Momentum[loop, D], Momentum[pa, D]], Pair[Momentum[loop, D], Momentum[pb, D]], Pair[Momentum[loop, D], Momentum[p, D]]},
  Pair[Momentum[loop, ___], Momentum[aX, ___]] :> axisCoefficients . {Pair[Momentum[loop, D], Momentum[pa, D]], Pair[Momentum[loop, D], Momentum[pb, D]], Pair[Momentum[loop, D], Momentum[p, D]]},
  Pair[Momentum[nn, ___], Momentum[loop, ___]] -> nLoop, Pair[Momentum[loop, ___], Momentum[nn, ___]] -> nLoop};
loopParallelSquare := parallelSquare[{Pair[Momentum[loop, D], Momentum[pa, D]], Pair[Momentum[loop, D], Momentum[pb, D]], Pair[Momentum[loop, D], Momentum[p, D]]}, gramMatrix[SP]];
(* the contracted interference of one graph -> sum of coefficient x family[indices] x HatPower[n]; the odd part in nn.loop is returned separately (it must vanish for XN, NX and is dropped) *)
toFamilies[expr0_, name_] := Module[{expr, fourDimensionalProducts, reduced},
  (* the four-dimensional products of the HVBM scheme equal the D-dimensional ones for the four-dimensional external momenta, and loop4.loop4 = loop.loop - loophat.loophat *)
  fourDimensionalProducts = {Pair[Momentum[loop], Momentum[loop]] -> Pair[Momentum[loop, D], Momentum[loop, D]] - loopHat2,
    Pair[Momentum[loop], Momentum[x_]] /; FreeQ[x, aX | nn] :> Pair[Momentum[loop, D], Momentum[x, D]], Pair[Momentum[x_], Momentum[loop]] /; FreeQ[x, aX | nn] :> Pair[Momentum[x, D], Momentum[loop, D]]};
  (* the spin sum of an outgoing quark k1 gives k1-slash: loop.k1 by momentum conservation *)
  expr = ExpandScalarProduct[expr0] /. axisLoopRules /. fourDimensionalProducts /. Pair[Momentum[k1, D], Momentum[loop, D]] -> ExpandScalarProduct[SPD[pa + pb - p, loop]] /.
    Pair[Momentum[loop, D - 4], Momentum[loop, D - 4]] -> loopHat2;
  reduced = normalReduction[expr, nLoop, loopParallelSquare, Pair[Momentum[loop, D], Momentum[loop, D]] - loopHat2];
  expr = reduced["Even"] /. loopRules[name];
  Require[FreeQ[expr, Momentum[loop, ___]], "A scalar product of the loop momentum was not expressed through the propagators", Short[Union[Cases[expr, h_[a___] /; !FreeQ[{a}, Momentum[loop, ___]] && FreeQ[h, Times | Plus | Power], Infinity]], 8]];
  expr = Expand[expr /. DjOf[q_] :> Dj @@ (Position[families[name], f_ /; sameUpToSign[f, q]][[1]])];
  (* every term is a coefficient times powers of the Dj (negative for propagators, positive for numerators) and of loopHat2 *)
  <|"Even" -> Total[Table[With[{powers = Table[Exponent[term, Dj[j]], {j, 4}], hat = Exponent[term, loopHat2]},
      (term/(Times @@ Table[Dj[j]^powers[[j]], {j, 4}] loopHat2^hat)) familyHeads[name] @@ (-powers) HatPower[hat]], {term, DeleteCases[termList[expr], 0]}]],   (* a graph can vanish exactly once eps(l).l = eps(l).r = 0 are imposed *)
    "OddLeafCount" -> LeafCount[reduced["Odd"]]|>];
bornKinematics;
virtualLeft = MapThread[colourGroups[prepare[shiftedGraph[#1, #2], bornReferences, bornConservation]] &, {amplitudes["Virtual"], routings}];
bornSumRight = If[bornSides === None, {}, conjugateGroups[Total[prepare[#, bornReferences, bornConservation] & /@ amplitudes["Born"]]]];
Length[virtualLeft]


(* ::Text:: *)
(*The interference 2 Re(M1 M0*) with the measure conventions of Processes/Hadronic: FeynArts writes the one-loop amplitude with d^D loop/(2 Pi)^D, the scalar integrals are defined with d^D loop/(i Pi^(D/2)), and the association virtual holds Sum_j [i c_j M0*] I_j per component, with the hat powers as HatPower[n]; the factor 2/(4 Pi)^(D/2) is included in step 4.  Checks: XN and NX vanish (their even parts in nn.loop; the odd parts, which integrate to zero, are dropped), and XX + NN vanishes in four dimensions once the hatted terms (which integrate to order eps) are dropped.*)


(* ::Input:: *)
virtualJobs[spinCase_] := Table[<|"Name" -> "virtual-" <> spinCase <> "-" <> ToString[g], "Kinematics" -> "Born", "Left" -> virtualLeft[[g]],
  "Right" -> bornSumRight, "Spin" -> spinCase, "References" -> bornReferences|>, {g, Length[virtualLeft]}];
virtualComponents = If[amplitudes["Virtual"] === {}, Association[# -> 0 & /@ components],
  Association @ Table[spinCase -> Module[{results = runJobs[virtualJobs[spinCase]], graphs},
    bornKinematics;
    graphs = MapThread[toFamilies[#1["Value"], #2[[1]]] &, {results, routings}];
    oddLeafCounts[spinCase] = Total[#["OddLeafCount"] & /@ graphs];
    Collect[I average[spinCase] normalization[spinCase] Total[#["Even"] & /@ graphs] /. u -> -s - t,
      {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, colourFactors]], {spinCase, components}]];
Require[FreeQ[virtualComponents, Complex | Eps | nLoop], "The virtual interference is not real or still contains nn.loop", Short[virtualComponents, 5]];
Require[virtualComponents["XN"] === 0 && virtualComponents["NX"] === 0, "The mixed components of the virtual interference do not vanish", {virtualComponents["XN"], virtualComponents["NX"]}];
Require[Collect[(virtualComponents["XX"] + virtualComponents["NN"]) /. HatPower[n_ /; n > 0] -> 0 /. HatPower[0] -> 1 /. D -> 4, {_KVirtualST, _KVirtualSU, _KVirtualTU}, Factor] === 0,
  "XX + NN of the virtual interference does not vanish in four dimensions"];
virtual = <|"TT" -> Collect[(virtualComponents["XX"] - virtualComponents["NN"])/2, {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, Factor],
  "Symmetric" -> Collect[(virtualComponents["XX"] + virtualComponents["NN"])/2, {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, Factor]|>;
{Length[virtual["TT"]], Union[Cases[virtual, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]] // Length, Map[Union[Cases[#, _HatPower, Infinity]] &, virtual],
  "leaf count of the dropped odd parts in nn.loop" -> (oddLeafCounts /@ components)}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*For every real-emission process of the class, the square of the sum of its diagrams for the four components, in terms of the invariants s, t, u, A, B, C, the hatted invariant spe3 and the normal component nr; the processes are added with their weights (classes.wl).  Reference vectors of the physical polarization sums: pa for the gluons p and k1, pb for k2.  Then the normal reduction: the parts odd in nr are dropped (they integrate to zero), nr^2 is replaced by k1par^2 + spe3.  Checks: the even parts of XN and NX vanish, and XX + NN vanishes in four dimensions without the hatted terms; both are evaluated exactly at random rational points of all invariants (the expressions have tens of thousands of terms).*)


(* ::Input:: *)
realReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa, k2 -> pb|>;
realKinematics;
realSides = Association @ KeyValueMap[#1 -> sidesOf[#2, realReferences, realConservation] &, amplitudes["Real"]];
realSquares = Association @ Table[spinCase -> Association @ KeyValueMap[
    #1 -> square["real-" <> StringReplace[#1, " " -> "_"], #2, "Real", realReferences, spinCase] &, realSides], {spinCase, components}];
realKinematics;
realParallelSquare = parallelSquare[{SP[pa, k1], SP[pb, k1], SP[p, k1]}, gramMatrix[SP]];      (* k1par^2 = d.G^-1.d *)
realReduced = Association @ Table[spinCase -> Association @ KeyValueMap[
    #1 -> normalReduction[#2, nr, realParallelSquare, SP[k1]] &, realSquares[spinCase]], {spinCase, components}];       (* SP[k1] = -spe3: nr^2 = k1par^2 + spe3 *)
(* the weighted sums per component (not brought to a common denominator: step 3 expands into partial fractions anyway); the colour factors are polynomials in Nc = CA *)
realComponents = Association @ Table[spinCase -> (average[spinCase] normalization[spinCase] Total[KeyValueMap[weights[#1] #2["Even"] &, realReduced[spinCase]]] /. Nc -> CA), {spinCase, components}];
Require[FreeQ[realComponents, Complex | nr | nn | aX | Eps], "The real-emission square is not real or still contains the normal component"];
Require[FreeQ[realComponents, Pair | Momentum | LorentzIndex], "A scalar product was left unevaluated in the real-emission square", Union[Cases[realComponents, _Pair, Infinity]]];
{realParallelSquare, Map[LeafCount, realComponents]}


(* ::Input:: *)
(* exact rational points of all invariants (C = -s4 - A - B), D, Nc, spe3 *)
realPoints = Table[With[{r = RandomInteger[{1, 97}, 8]/RandomInteger[{98, 197}, 8]},
  {s -> 1, t -> -r[[1]], u -> -r[[2]] (1 - r[[1]]), A -> -r[[3]], B -> -r[[4]], D -> 4 - r[[5]], CA -> 2 + 10 r[[6]], CF -> 1 + r[[7]], spe3 -> r[[8]], Nf -> 5}], {3}];
atPoint[expr_, point_] := expr /. C -> -(s + t + u) - A - B /. point;
Require[And @@ Table[atPoint[realComponents["XN"], point] === 0 && atPoint[realComponents["NX"], point] === 0, {point, realPoints}],
  "The even parts of the mixed components of the real emission do not vanish"];
Require[And @@ Table[atPoint[(realComponents["XX"] + realComponents["NN"]) /. {D -> 4, spe3 -> 0}, point] === 0, {point, realPoints}],
  "XX + NN of the real emission does not vanish in four dimensions"];
real = <|"TT" -> (realComponents["XX"] - realComponents["NN"])/2, "Symmetric" -> (realComponents["XX"] + realComponents["NN"])/2|>;
Require[Max[Exponent[#, spe3] & /@ termList[Expand[real["TT"] + real["Symmetric"]]]] <= 1, "The real-emission square has spe3^2 terms (step 3 has the masters in D + 2 only)"];
{"odd parts dropped (leaf counts)" -> Map[Total[KeyValueMap[LeafCount[#2["Odd"]] &, #]] &, realReduced], LeafCount /@ real}


(* ::Text:: *)
(*Check of the gauge invariance of the real-emission square: for the first process of the class and the XX component (the one without normal components), the physical polarization sums with other reference vectors (pb for p and k1, pa for k2) give the same result.  The amplitudes are prepared again, since the conditions eps(l).r = 0 depend on the reference vectors.*)


(* ::Input:: *)
alternativeReferences = <|pa -> pb, pb -> pa, p -> pb, k1 -> pb, k2 -> pa|>;
firstProcess = First[Keys[amplitudes["Real"]]];
realKinematics;
alternative = square["gauge-" <> StringReplace[firstProcess, " " -> "_"],
  sidesOf[amplitudes["Real"][firstProcess], alternativeReferences, realConservation], "Real", alternativeReferences, "XX"];
(* the difference is evaluated exactly at random rational points (all invariants, D and Nc), with C = -s4 - A - B *)
gaugePoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], Nc -> 2 + 10 v[[6]], spe3 -> v[[7]]}], {3}];
referenceDifference = Table[(alternative - realSquares["XX"][firstProcess]) /. C -> -(s + t + u) - A - B /. point, {point, gaugePoints}];
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
(*Save the Born square, the virtual interference, the real-emission square (projections TT and Symmetric) and the Born table for steps 3 and 4, and close the parallel kernels.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Virtual" -> virtual, "Real" -> real, "BornTable" -> bornTable, "Class" -> class, "Spec" -> spec,
  "Families" -> families, "Projections" -> "TT = (XX - NN)/2 (the cos(2 Phi) coefficient, the hard part of A_TT); Symmetric = (XX + NN)/2 (evanescent in the HVBM scheme)",
  "Variables" -> "s, t, u = -s-t; A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 = -s4-A-B; spe3; loopHat2 via HatPower[n]; Born table: |M|^2/gs^4 in s, t, D, CF, CA"|>;
SaveResult["Algebra", algebra, "FeynCalc/class-" <> class];
CloseKernels[];
Keys[algebra]
