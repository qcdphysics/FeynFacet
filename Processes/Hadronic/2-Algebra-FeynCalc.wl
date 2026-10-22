(* ::Title:: *)
(*Proton-proton classes, step 2 (FeynCalc route): spin and colour algebra*)


(* ::Text:: *)
(*Input: the amplitudes of step 1 for the class chosen below.  Output: Results/FeynCalc/class-<n>/Algebra.wl with*)
(*  Born:     the Born square |M|^2 (UU and LL), summed over final spins and colours, averaged over the initial ones, in D dimensions;*)
(*  Virtual:  the interference 2 Re(M_1loop M_Born^*) (UU and LL), averaged and summed in the same way, written in terms of scalar one-loop integrals (families VirtualST, VirtualSU, VirtualTU below) with D-dimensional coefficients;*)
(*  Real:     the sum over the real-emission processes of the class, each weighted as in classes.wl, of |M|^2 (UU and LL) with the unobserved partons' momenta k1, k2 kept explicit (step 3 integrates them).*)
(*The FORM and Tracer routes write the same file.*)
(**)
(*Spin sums (HVBM scheme for LL): the physical external momenta pa, pb, p are four-dimensional; the unobserved momenta k1, k2 and the loop momentum have (D-4)-dimensional parts, which appear through the invariant spe3 = k1hat.k1hat = k2hat.k2hat = -k1hat.k2hat in the real emission.  Incoming quark: p-slash/2 (UU) or +-p-slash gamma5/2 (LL, minus for an antiquark); incoming gluon: physical polarization sum averaged with 1/(D-2) (UU) or the circular density i eps^{mu nu l r}/(2 l.r) (LL; r is the other incoming momentum).  Outgoing gluons are summed with the physical sum -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r), with r = pa for the observed parton and r = pa or pb for the unobserved ones; the results are independent of these choices by gauge invariance (the Born squares were checked against the textbook with this prescription in Processes/PartonicBorn, and the real-emission squares are checked below against an alternative choice of the reference vectors).*)
(**)
(*Invariants: s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2; real emission: s4 = (k1+k2)^2 = s + t + u, and the three independent angular denominators A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 with A + B + C = -s4; the denominators with k2 follow from k2 = pa + pb - p - k1 (step 3).*)


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
class = "13";
workers = "auto";
resume = "no";
loaded = LoadAmplitudes["class-" <> class];
amplitudes = loaded["FeynCalc"]; spec = loaded["Spec"]; weights = loaded["Weights"];
spec


(* ::Section:: *)
(*Conventions*)


(* ::Text:: *)
(*The conventions and kinematics are installed by functions, because every parallel kernel installs them too.  One unobserved momentum is not independent: inside the amplitudes k2 is replaced by pa + pb - p - k1 (real emission) and k1 by pa + pb - p (Born and one loop), except in the external spinors, whose spin sums give the slashed momentum itself.*)


(* ::Input:: *)
conventions := (FCSetDiracGammaScheme["BMHV"]; $LeviCivitaSign = -1);
conventions;
bornKinematics := (FCClearScalarProducts[]; SetMandelstam[s, t, u, pa, pb, -p, -k1, 0, 0, 0, 0];
  Scan[(Momentum[#, D - 4] = 0) &, {pa, pb, p, k1}]);
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
  SP[pa, k2] = (s + t + A)/2; SP[pb, k2] = (s + u + B)/2; SP[p, k2] = (-t - u - C)/2);
bornConservation = k1 -> pa + pb - p;
realConservation = k2 -> pa + pb - p - k1;
incoming = spec["Born"] /. None -> spec["Real"][[1, 1]];       (* the two incoming partons *)
incoming = incoming[[1 ;; 2]];
polarizedQuarkSign = <|pa -> If[StringEndsQ[incoming[[1]], "b"], -1, 1], pb -> If[StringEndsQ[incoming[[2]], "b"], -1, 1]|>;
incomingGluons = Select[{pa, pb}, incoming[[Position[{pa, pb}, #][[1, 1]]]] === "g" &];
incomingQuarks = Complement[{pa, pb}, incomingGluons];
average[spinCase_] := Times @@ (Switch[incoming[[#]], "g", If[spinCase === "LL", 1, 1/(D - 2)]/(Nc^2 - 1), _, 1/(2 Nc)] & /@ {1, 2});
colourFactors[expr_] := Factor[expr /. SUNN -> Nc] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
{incomingQuarks, incomingGluons, average["UU"]}


(* ::Section:: *)
(*The squaring method*)


(* ::Text:: *)
(*A square sum_{i,j} M_i M_j^* of many diagrams is not formed as one product.  Five steps keep the expressions small:*)
(**)
(*1. Every amplitude is simplified once, before any product: Lorentz contraction, Dirac algebra on the open spinor chains, expansion of the scalar products with one independent unobserved momentum, and the conditions eps(l).l = 0 and eps(l).r = 0 for every gluon l, where r is the reference vector of the same gluon's polarization sum used later (the physical sum and the circular density are both transverse to l and to r).  For the gluon whose momentum was eliminated by momentum conservation, eps(l).l = 0 is imposed in its eliminated form.*)
(*2. Colour is separated from the Lorentz algebra: each amplitude is written as sum_x c_x L_x with distinct colour structures c_x (products of generators and structure constants).  The colour factor of c_x c_y^* is a closed colour network; it is evaluated once (structure constants written through fundamental traces, f^{abc} = -2 i (Tr[T^a T^b T^c] - Tr[T^b T^a T^c])) and remembered; equal networks with differently named summed indices are recognized by renaming the indices canonically.*)
(*3. Only the products with i <= j are computed: for i < j, M_i M_j^* + M_j M_i^* = 2 Re(M_i M_j^*) (all invariants are real).*)
(*4. After the fermion spin sums, the gluon polarization sums are contracted one gluon at a time before evaluating the Dirac traces: the terms are grouped by the factors that carry that gluon's two indices, and each distinct tensor is contracted with the polarization sum only once (and remembered).*)
(*5. The products (i, j) are independent: they are distributed over parallel kernels, and each is written to its own file with the time it took, so that the progress is visible and an interrupted run can be continued.*)


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
circular[l_, mu_, nu_, r_] := Contract[I LC[mu, nu, al, be] FV[l, al] FV[r, be]/(2 SP[l, r])];    (* four-dimensional *)
tag[l_] := Symbol[ToString[l] <> "Spin"];
indexNames[l_] := {Symbol["mu" <> ToString[l]], Symbol["nu" <> ToString[l]]};
tensorPart[term_, {mu_, nu_}] := Times @@ Select[factorList[term], !FreeQ[#, LorentzIndex[mu | nu, ___]] &];
(* the contraction uses the installed scalar products, which differ between the Born and the real kinematics: the kinematics are part of the key *)
(* Each worker reuses identical contractions across diagram pairs. The explicit
   spin projections and the original BMHV trace evaluation are unchanged. *)
currentKinematics = None;
algebraCacheScope[] := {"Hadronic", class, currentKinematics, "BMHV", $LeviCivitaSign};
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
(* the colour-weighted product of the Lorentz parts of one pair: spin sums, traces, polarization sums *)
spinSum[product_, spinCase_, references_, gluonList_] := Module[{expr = product, polarizedQuarks, polarizedGluons},
  polarizedQuarks = If[spinCase === "LL", incomingQuarks, {}];
  polarizedGluons = If[spinCase === "LL", incomingGluons, {}];
  expr = expr /. (With[{l = #, lTag = tag[#]}, Spinor[c_. Momentum[l, D], m_, o_] :> Spinor[c Momentum[lTag, D], m, o]] & /@ polarizedQuarks);
  expr = expr /. Flatten[Table[{Momentum[Polarization[l, I, ___], D] -> LorentzIndex[indexNames[l][[1]], D],
    Momentum[Polarization[l, -I, ___], D] -> LorentzIndex[indexNames[l][[2]], D]}, {l, gluonList}]];
  expr = FermionSpinSum[expr];
  expr = expr /. (With[{l = #, lTag = tag[#], sign = polarizedQuarkSign[#]}, DiracGamma[Momentum[lTag, D], D] :> sign GSD[l] . GA[5]] & /@ polarizedQuarks);
  (* Contract polarization sums before expanding the same BMHV traces. *)
  (* outgoing gluons first, then the incoming ones *)
  Do[expr = project[expr, l, polarizationSum[l, Sequence @@ indexNames[l], references[l]]], {l, Complement[gluonList, {pa, pb}]}];
  Do[expr = project[expr, l, If[MemberQ[polarizedGluons, l], circular[l, Sequence @@ indexNames[l], references[l]],
      polarizationSum[l, Sequence @@ indexNames[l], references[l]]]], {l, Intersection[gluonList, {pa, pb}]}];
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
"squaring functions defined"


(* ::Text:: *)
(*Parallel kernels: each loads FeynFacet and FeynCalc, installs the conventions and receives the functions above; each job installs its kinematics (Born or real).  Every job (one product of the squaring method) is written to Results/FeynCalc/class-<n>/pairs/<name>.wl with its time.*)


(* ::Input:: *)
feynFacetFile = FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}];
nWorkers = If[workers === "auto", Min[$ProcessorCount, 8], ToExpression[workers]];   (* "auto": one per processor core, at most 8 *)
If[nWorkers > 1,
  LaunchKernels[nWorkers];
  With[{file = feynFacetFile}, ParallelEvaluate[Get[file]; LoadFeynCalc[]]];
  DistributeDefinitions[conventions, bornKinematics, realKinematics, gluonsOf, transversality, eliminatedTransversality, conserve, prepare,
    colourHead, factorList, termList, colourPart, colourGroups, colourValue, polarizationSum, circular, tag, indexNames, tensorPart,
    projectOnce, project, spinSum, interference, conjugateGroups, realPartOf, incomingQuarks, incomingGluons, polarizedQuarkSign,
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
checkPoint = {s -> 1, t -> -3/10, u -> -2/5, A -> -1/7, B -> -2/9, C -> -3/10 + 1/7 + 2/9, D -> 37/10, Nc -> 3, spe3 -> 1/11, gs -> 1, Nf -> 5};
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
(*The Born square with the reference vectors pa -> pb, pb -> pa (incoming gluons), p -> pa, k1 -> pa.*)


(* ::Input:: *)
bornReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa|>;
bornKinematics;
bornSides = If[amplitudes["Born"] === {}, None, sidesOf[amplitudes["Born"], bornReferences, bornConservation]];
born = If[bornSides === None, <|"UU" -> 0, "LL" -> 0|>,
  Association @ Table[spinCase -> With[{sq = square["born", bornSides, "Born", bornReferences, spinCase]},
    bornKinematics; Factor[colourFactors[Together[average[spinCase] sq]] /. u -> -s - t]], {spinCase, {"UU", "LL"}}]];
Require[FreeQ[born, Complex | Pair | Momentum | LorentzIndex], "The Born square is not a real scalar", born];
born


(* ::Text:: *)
(*The Born squares must agree with the table of Processes/PartonicBorn (which was checked against the textbook).*)


(* ::Input:: *)
(* the table is read only for a class with a Born term (the classes without one do not need Processes/PartonicBorn here) *)
If[spec["Born"] =!= None,
  bornTable = Get[FileNameJoin[{NotebookDirectory[], "..", "PartonicBorn", "Results", "Born.wl"}]];
  Require[And @@ Table[Simplify[born[sc]/gs^4 - bornTable["Born"][sc][StringRiffle[spec["Born"][[1 ;; 2]], " "] <> " -> " <> spec["Born"][[3]]]] === 0, {sc, {"UU", "LL"}}],
    "The Born square differs from Processes/PartonicBorn", born]];
"Born checked"


(* ::Section:: *)
(*Virtual: the one-loop interference in terms of scalar integrals*)


(* ::Text:: *)
(*The one-loop amplitude is contracted with the conjugate Born amplitude and the spin sums as above; the loop momentum then appears only in scalar products with the external momenta and in loop^2, all of which are combinations of the inverse propagators, and every term becomes a scalar integral of one of the families*)
(**)
(*   VirtualST[n1,n2,n3,n4]: 1/[(loop^2)^n1 ((loop+pa)^2)^n2 ((loop+pa+pb)^2)^n3 ((loop+p)^2)^n4]   (external legs pa, pb, -pd, -p in cyclic order: invariants s and t),*)
(*   VirtualSU: ... ((loop+pa+pb)^2)^n3 ((loop+pa+pb-p)^2)^n4   (legs pa, pb, -p, -pd: s and u),      VirtualTU: ((loop+pa)^2)^n2 ((loop+pa-p)^2)^n3 ((loop+pa-p+pb)^2)^n4   (legs pa, -p, pb, -pd: t and u),*)
(**)
(*with negative indices for numerators.  FeynArts routes the loop momentum differently in every graph; for each graph the routing loop -> +-loop + shift that maps its propagators onto one family is found by trying the shifts.  Kira reduces the families in step 3.*)


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
(*For each graph: shift the loop momentum and replace each loop propagator by an inert symbol DjOf[q] (a tree propagator next to the loop is evaluated); simplify the graph once and separate its colour as described above, and form its product with the conjugate of the summed Born amplitude (one job of the squaring method per graph and spin case).  Then write the loop momentum's scalar products through the inverse propagators Dj of its family (the four equations loop^2 = D1, 2 loop.q_j = D_j - D_1 - q_j^2 are solved for loop^2, loop.pa, loop.pb, loop.p), and read off each term's family indices.  In the HVBM scheme the four-dimensional loop.loop is loop^2 - loopHat2; loopHat2 is kept as a separate variable (it shifts the dimension of the integral in step 3, not the dimension D of the coefficients).*)


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
(* the contracted interference of one graph -> sum of coefficient x family[indices] x HatPower[n] *)
toFamilies[expr0_, name_] := Module[{expr, fourDimensionalProducts},
  (* the four-dimensional products of the HVBM scheme equal the D-dimensional ones for the four-dimensional external momenta, and loop4.loop4 = loop.loop - loophat.loophat *)
  fourDimensionalProducts = {Pair[Momentum[loop], Momentum[loop]] -> Pair[Momentum[loop, D], Momentum[loop, D]] - loopHat2,
    Pair[Momentum[loop], Momentum[x_]] :> Pair[Momentum[loop, D], Momentum[x, D]], Pair[Momentum[x_], Momentum[loop]] :> Pair[Momentum[x, D], Momentum[loop, D]]};
  (* the spin sum of an outgoing quark k1 gives k1-slash: loop.k1 by momentum conservation *)
  expr = ExpandScalarProduct[expr0] /. fourDimensionalProducts /. Pair[Momentum[k1, D], Momentum[loop, D]] -> ExpandScalarProduct[SPD[pa + pb - p, loop]] /. loopRules[name] /. Pair[Momentum[loop, D - 4], Momentum[loop, D - 4]] -> loopHat2;
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
(*with the real parts of the integrals in the physical region (s > 0) taken in step 3.  The association virtual below holds Sum_j [i c_j M0*] I_j for both spin cases, with the loop-momentum hat powers as HatPower[n]; the factor 2/(4 Pi)^(D/2) is included in step 4.*)


(* ::Input:: *)
virtualJobs[spinCase_] := Table[<|"Name" -> "virtual-" <> spinCase <> "-" <> ToString[g], "Kinematics" -> "Born", "Left" -> virtualLeft[[g]],
  "Right" -> bornSumRight, "Spin" -> spinCase, "References" -> bornReferences|>, {g, Length[virtualLeft]}];
virtual = If[amplitudes["Virtual"] === {}, <|"UU" -> 0, "LL" -> 0|>,
  Association @ Table[spinCase -> Module[{results = runJobs[virtualJobs[spinCase]]},
    bornKinematics;
    Collect[I average[spinCase] Total[MapThread[toFamilies[#1["Value"], #2[[1]]] &, {results, routings}]] /. u -> -s - t,
      {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, colourFactors]], {spinCase, {"UU", "LL"}}]];
Require[FreeQ[virtual, Complex], "The virtual interference is not real", virtual];
{Length[virtual["UU"]], Union[Cases[virtual, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]] // Length}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*For every real-emission process of the class, the square of the sum of its diagrams, with the spin sums of the two cases, in terms of the invariants s, t, u, A, B, C and the hatted invariant spe3; the processes are added with their weights (classes.wl).  Reference vectors of the physical polarization sums: pa for the gluons p and k1, pb for k2 (and the other incoming momentum for an incoming gluon).*)


(* ::Input:: *)
realReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa, k2 -> pb|>;
realKinematics;
realSides = Association @ KeyValueMap[#1 -> sidesOf[#2, realReferences, realConservation] &, amplitudes["Real"]];
realSquares = Association @ Table[spinCase -> Association @ KeyValueMap[
    #1 -> square["real-" <> StringReplace[#1, " " -> "_"], #2, "Real", realReferences, spinCase] &, realSides], {spinCase, {"UU", "LL"}}];
realKinematics;
(* the sum is not brought to a common denominator (step 3 expands it into partial fractions anyway); the colour factors are polynomials in Nc = CA *)
real = Association @ Table[spinCase -> (average[spinCase] Total[KeyValueMap[weights[#1] #2 &, realSquares[spinCase]]] /. Nc -> CA),
  {spinCase, {"UU", "LL"}}];
Require[FreeQ[real, Complex], "The real-emission square is not real"];
Require[FreeQ[real, Pair | Momentum | LorentzIndex], "A scalar product was left unevaluated in the real-emission square", Union[Cases[real, _Pair, Infinity]]];
Require[Max[Exponent[#, spe3] & /@ termList[Expand[real["LL"]]]] <= 1, "The real-emission square has spe3^2 terms"];
LeafCount /@ real


(* ::Text:: *)
(*Check of the gauge invariance of the real-emission square: for the first process of the class (and the spin average), the physical polarization sums with other reference vectors (pb for p and k1, pa for k2) give the same result.  The amplitudes are prepared again, since the conditions eps(l).r = 0 depend on the reference vectors.*)


(* ::Input:: *)
alternativeReferences = <|pa -> pb, pb -> pa, p -> pb, k1 -> pb, k2 -> pa|>;
firstProcess = First[Keys[amplitudes["Real"]]];
realKinematics;
alternative = square["gauge-" <> StringReplace[firstProcess, " " -> "_"],
  sidesOf[amplitudes["Real"][firstProcess], alternativeReferences, realConservation], "Real", alternativeReferences, "UU"];
(* the difference is evaluated exactly at random rational points (all invariants, D and Nc), with C = -s4 - A - B *)
gaugePoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], Nc -> 2 + 10 v[[6]]}], {3}];
referenceDifference = Table[(alternative - realSquares["UU"][firstProcess]) /. C -> -(s + t + u) - A - B /. point, {point, gaugePoints}];
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
(*Save the Born square, the virtual interference and the real-emission square for step 3, and close the parallel kernels.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Virtual" -> virtual, "Real" -> real, "Class" -> class, "Spec" -> spec,
  "Families" -> families, "Variables" -> "s, t, u = -s-t; A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 = -s4-A-B; spe3; loopHat2 via HatPower[n]"|>;
SaveResult["Algebra", algebra, "FeynCalc/class-" <> class];
CloseKernels[];
Keys[algebra]
