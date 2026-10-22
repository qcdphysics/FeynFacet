(* Independent finite assembly of the seven transverse-spin-transfer channels.
   No FeynFacet package, process step 4, or distribution helper is loaded.
   Algebra/Reduced are shared mathematical inputs; Transfer is a target only.
   See DERIVATION.txt for the distribution pullback and dimensional factors. *)
$HistoryLength = 0;
$MaxExtraPrecision = 200;
workDirectory = DirectoryName[ExpandFileName[$InputFileName]];
sourceDirectory = Environment["INDEPENDENT_SOURCE"];
If[!StringQ[sourceDirectory], Print["Set INDEPENDENT_SOURCE to the repository root."]; Quit[2]];
selectedRoute = Environment["INDEPENDENT_ROUTE"];
If[!StringQ[selectedRoute], selectedRoute = "FeynCalc"];
classNumber = Environment["INDEPENDENT_CLASS"];
If[!StringQ[classNumber], classNumber = "10"];
physicalClassNumber = StringReplace[classNumber, RegularExpression["x$"] -> ""];
countertermsOnly = Environment["INDEPENDENT_CT_ONLY"] === "1";
outputRoot = Environment["INDEPENDENT_OUTPUT"];
If[!StringQ[outputRoot], Print["Set INDEPENDENT_OUTPUT to a separate output directory."]; Quit[2]];
outputDirectory = FileNameJoin[{outputRoot, "class-" <> classNumber}];
If[!DirectoryQ[outputDirectory], CreateDirectory[outputDirectory, CreateIntermediateDirectories -> True]];
started = AbsoluteTime[];
inputFiles = FileNameJoin[{sourceDirectory, "Processes", "HadronicTransfer", "Results",
   selectedRoute, "class-" <> classNumber, # <> ".wl"}] & /@ {"Algebra", "Reduced", "Transfer"};
inputHashes = AssociationThread[inputFiles, IntegerString[FileHash[#, "SHA256"], 16, 64] & /@ inputFiles];
{amplitudes, integrals, target} = Get /@ inputFiles;
assumptions = s > 0 && 0 < v < 1 && 0 < w < 1 && muR > 0 && muF > 0 && muFF > 0 && CA > 1 && Nf > 0;
checks = {}; failures = {};
(* Expand logarithms only when every real rational factor has a certified
   nonzero sign and the complete argument is positive. Unknown and complex
   factors preserve the original logarithm, including its phase. *)
Clear[certifiedFactorSign, positiveRationalLog, certifiedLogOf, certifiedLogs];
certifiedFactorSign[factor_, assumptions_] := certifiedFactorSign[factor, assumptions] = Which[
  !TrueQ[FullSimplify[Element[factor, Reals], assumptions]], 0,
  TrueQ[FullSimplify[factor > 0, assumptions]], 1,
  TrueQ[FullSimplify[factor < 0, assumptions]], -1,
  True, 0];

positiveRationalLog[value_] := If[(IntegerQ[value] || Head[value] === Rational) && value > 0,
  Total[#[[2]] Log[#[[1]]] & /@ FactorInteger[Numerator[value]]] -
    Total[#[[2]] Log[#[[1]]] & /@ FactorInteger[Denominator[value]]], Log[value]];

certifiedLogOf[argument_, assumptions_] := certifiedLogOf[argument, assumptions] = Module[
  {piPlaceholder, factored, factor, exponent, sign, productSign = 1, terms = {}, valid = True},
  (* Treat Pi as a formal positive factor, so a numeric coefficient such as
     16 Pi^2/9 is decomposed rather than left as one indivisible number. *)
  factored = Quiet[Check[FactorList[Together[argument /. Pi -> piPlaceholder]], $Failed]];
  If[factored === $Failed || !ListQ[factored] ||
    !AllTrue[factored, MatchQ[#, {_, _Integer}] &], Return[Log[argument]]];
  Do[
    factor = pair[[1]] /. piPlaceholder -> Pi;
    exponent = pair[[2]];
    sign = certifiedFactorSign[factor, assumptions];
    If[sign === 0, valid = False; Break[]];
    If[sign === -1, factor = Expand[-factor]; productSign *= (-1)^exponent];
    AppendTo[terms, exponent positiveRationalLog[factor]], {pair, factored}];
  If[!valid, Log[argument], If[productSign === 1, Total[terms], Log[argument]]]];

certifiedLogs[expression_, assumptions_] :=
  expression /. Log[argument_] :> certifiedLogOf[argument, assumptions];

certifiedZeroQ[expression_, assumptions_] := Module[{rational},
  rational = Together[Expand[certifiedLogs[expression, assumptions]]];
  TrueQ[rational === 0 || FullSimplify[rational, assumptions] === 0]];

zeroQ[expression_] := certifiedZeroQ[expression /. CF -> (CA^2 - 1)/(2 CA), assumptions];
tupleEqualQ[left_, right_] := And @@ MapThread[zeroQ[#1 - #2] &, {left, right}];
finiteSlotEqualQ[left_, right_, slot_] := zeroQ[left[[slot]] - right[[slot]]];
record[label_, result_] := (AppendTo[checks, <|"name" -> label, "passed" -> TrueQ[result]|>];
  If[!TrueQ[result], AppendTo[failures, label]];
  Print["CHECK ", label, " ", TrueQ[result], " ", DateString[{"Hour", ":", "Minute", ":", "Second"}]]; TrueQ[result]);
endValue[expression_] := Together[Normal[Series[expression /. w -> 1 - endpointDistance,
  {endpointDistance, 0, 0}]]];
laurent[expression_, power_] := SeriesCoefficient[expression, {eps, 0, power}];
channelParticles = <|"1" -> {"q", "qp", "q"}, "3" -> {"q", "qpb", "q"},
  "6" -> {"q", "q", "q"}, "10" -> {"q", "g", "qb"}, "11" -> {"q", "qb", "q"},
  "13" -> {"q", "g", "q"}, "17" -> {"qb", "q", "q"}|>;
If[!KeyExistsQ[channelParticles, physicalClassNumber], Print["Unsupported transfer class: ", classNumber]; Quit[2]];
{incomingPolarized, incomingUnpolarized, outgoingPolarized} = channelParticles[physicalClassNumber];

(* P(z) = regular(z) + cusp [1/(1-z)]_+ + endpoint delta(1-z).
   These LO kernels are written here from their defining splitting functions. *)
transverseKernel = {-2 CF, 2 CF, 3 CF/2};
unpolarizedKernel[child_, parent_] := Which[
  child === "g" && parent === "g", {2 CA (1/z - 2 + z (1 - z)), 2 CA, (11 CA - 2 Nf)/6},
  child =!= "g" && parent =!= "g" && child === parent, {-CF (1 + z), 2 CF, 3 CF/2},
  child =!= "g" && parent =!= "g", {0, 0, 0},
  child === "g", {CF (1 + (1 - z)^2)/z, 0, 0},
  True, {(z^2 + (1 - z)^2)/2, 0, 0}];
conjugation = <|"q" -> "qb", "qb" -> "q", "qp" -> "qpb", "qpb" -> "qp", "g" -> "g"|>;
bornSquare[spin_, particles_List] := Module[{names, found},
  names = (StringRiffle[Take[#, 2], " "] <> " -> " <> Last[#]) & /@
    {particles, Lookup[conjugation, particles]};
  found = Select[names, KeyExistsQ[amplitudes["BornTable"][spin], #] &];
  If[found === {}, 0, amplitudes["BornTable"][spin][First[found]]]];
bornDensity[square_, fraction_] := alphaS^2/(s fraction) *
  (square /. D -> 4 - 2 eps /. {t -> -s (1 - fraction), u -> -s fraction});

(* One generic pullback, rather than one tuple formula for each leg.
   z = map(w), W(w,eps) is the density after resolving the Born delta.
   A = W (1-w)/(1-map(w)); A(1)=W(1)/map'(1).
   The endpoint log follows from analytic regularization in z, not w. *)
pullback[kernel_List, map_, weight_] := Module[{cusp = kernel[[2]], endpoint = kernel[[3]],
  analyticWeight, endpointWeight, derivative},
  derivative = D[map, w] /. w -> 1;
  analyticWeight = Cancel[weight (1 - w)/(1 - map)];
  endpointWeight = endValue[analyticWeight];
  {Together[(kernel[[1]] /. z -> map) weight + cusp (analyticWeight - endpointWeight)/(1 - w)],
   cusp endpointWeight, 0, endpointWeight (endpoint + cusp Log[derivative])}];
subtraction[tuple_, scale_] := Module[{full},
  full = alphaS/(2 Pi eps) Exp[eps (Log[4 Pi] - EulerGamma + Log[muR^2/scale^2])] tuple;
  <|"Pole" -> (laurent[#, -1] & /@ full), "Finite" -> (laurent[#, 0] & /@ full)|>];
zeroTuple = {0, 0, 0, 0};
independentCounterterms[spin_, fragmentationExponent_ : 2] := Module[{sameSquare, born, aTerm, bTerm, cTerm, renTerm,
  other, otherSquare, otherKernel, flavourWeight, map, density, terms},
  sameSquare = bornSquare[spin, channelParticles[physicalClassNumber]];
  born = bornDensity[sameSquare, v];
  aTerm = subtraction[pullback[transverseKernel, w, born/w], muF];
  terms = Table[
    otherSquare = bornSquare[spin, {incomingPolarized, other, outgoingPolarized}];
    otherKernel = unpolarizedKernel[other, incomingUnpolarized];
    flavourWeight = If[MemberQ[{"qp", "qpb"}, other] &&
      Intersection[channelParticles[physicalClassNumber], {"qp", "qpb"}] === {}, Nf - 1, 1];
    map = (1 - v)/(1 - v w);
    density = bornDensity[otherSquare, v w] v w/(1 - v);
    flavourWeight pullback[otherKernel, map, density], {other, {"q", "qb", "qp", "qpb", "g"}}];
  bTerm = subtraction[Total[terms], muF];
  map = 1 - v + v w;
  density = v w/map^2 (born /. v -> v w/map) map^(fragmentationExponent eps);
  cTerm = subtraction[pullback[transverseKernel, map, density], muFF];
  renTerm = subtraction[-(11 CA - 2 Nf)/3 {0, 0, 0, born}, muR];
  <|"renormalization" -> renTerm, "incoming a" -> aTerm, "incoming b" -> bTerm, "fragmentation" -> cTerm|>];

ctResults = Association@Table[spin -> independentCounterterms[spin], {spin, {"XX", "NN"}}];
Do[Do[Do[record[spin <> "/counterterm/" <> name <> "/" <> ToString[part],
  zeroQ[ctResults[spin][name]["Finite"][[part]] - target["Counterterms"][spin][name][[part]]]],
  {part, 4}], {name, Keys[ctResults[spin]]}], {spin, {"XX", "NN"}}];

(* Resolve endpoint powers using valuations, not factor signs at a sample point.
   If base(h)=h^m a(h), a(0) != 0, then base^q=h^(m q) a^q.
   The integer h-power remains ordinary; tau records only its eps slope. *)
endpointPower[base_, exponent_] := Module[{rational, order, rest},
  If[FreeQ[base, thresholdDistance], Return[base^exponent]];
  rational = Together[base];
  order = Exponent[Numerator[rational], thresholdDistance, Min] -
    Exponent[Denominator[rational], thresholdDistance, Min];
  rest = Cancel[rational/thresholdDistance^order];
  rest^exponent thresholdDistance^(order (exponent /. eps -> 0)) *
    thresholdTag^(-order Coefficient[Expand[exponent], eps])];
endpointSeries[expression_] := Together[Normal[Series[expression,
  {thresholdDistance, 0, 0}]]];
(* Expand each distinct epsilon-dependent product only once. Numerator terms
   frequently differ by rational functions of s,v,w but share the same short
   Gamma/phase-space factor. Expanding their full sum repeatedly is expensive. *)
groupedLaurent[expression_, orders_List, endpointOrders_List : {}] := Module[
  {terms, pairs, groups, factors, result, endpoints, factorSeries, allOrders,
   coefficient, freeCoefficient, minimumOrder = 0, phaseStarted,
   groupingSeconds, seriesSeconds = 0, endpointSeconds = 0},
  phaseStarted = AbsoluteTime[];
  terms = With[{expanded = Expand[expression]}, If[Head[expanded] === Plus, List @@ expanded, {expanded}]];
  pairs = Table[factors = If[Head[term] === Times, List @@ term, {term}];
    {Times @@ Select[factors, FreeQ[#, eps] &],
     Times @@ Select[factors, !FreeQ[#, eps] &]}, {term, terms}];
  groups = GatherBy[pairs, Last];
  groupingSeconds = AbsoluteTime[] - phaseStarted;
  Print["    epsilon factors: ", Length[groups], " groups from ", Length[terms], " terms ", DateString[]];
  result = ConstantArray[0, Length[orders]];
  endpoints = ConstantArray[0, Length[endpointOrders]];
  allOrders = Union[orders, endpointOrders];
  Do[
    phaseStarted = AbsoluteTime[];
    factorSeries = Series[group[[1, 2]], {eps, 0, Max[allOrders]}];
    If[Head[factorSeries] === SeriesData,
      minimumOrder = Min[minimumOrder, factorSeries[[4]]/factorSeries[[6]]]];
    freeCoefficient = Total[group[[All, 1]]];
    result += (SeriesCoefficient[factorSeries, {eps, 0, #}] & /@ orders) freeCoefficient;
    seriesSeconds += AbsoluteTime[] - phaseStarted;
    (* Series each complete group contribution, retaining every negative h
       power until groups are summed. A full b1(h) is never constructed. *)
    phaseStarted = AbsoluteTime[];
    endpoints += (endpointSeries[SeriesCoefficient[factorSeries, {eps, 0, #}]
      freeCoefficient] & /@ endpointOrders);
    endpointSeconds += AbsoluteTime[] - phaseStarted, {group, groups}];
  Print["    seconds: grouping=", N[groupingSeconds], " Laurent=", N[seriesSeconds],
    " endpoint=", N[endpointSeconds]];
  <|"Bulk" -> result, "Endpoint" -> (Together /@ endpoints), "MinimumGroupOrder" -> minimumOrder|>];
independentReal[expression_, endpointRules_, masterRules_] := Module[{tagged, maxTag, branch, rate,
  pole2 = 0, pole1 = zeroTuple, finite = zeroTuple, expansion,
  bm2, bm1, b0, bm2End, bm1End, b0End, b1End},
  tagged = expression /. w -> 1 - thresholdDistance;
  tagged = tagged /. Power[base_, exponent_] /; !FreeQ[exponent, eps] :>
    endpointPower[base, exponent];
  tagged = Expand[tagged];
  maxTag = Exponent[tagged, thresholdTag];
  Print["  endpoint valuation complete; regulator order ", maxTag, " ", DateString[]];
  record["real/valid-regulator-exponents", PolynomialQ[tagged, thresholdTag] && IntegerQ[maxTag] && maxTag >= 0];
  branch = Coefficient[tagged, thresholdTag, 0];
  expansion = groupedLaurent[branch, {-2, -1, 0}];
  {bm2, pole1[[1]], finite[[1]]} = expansion["Bulk"];
  record["real/regular-supported-pole-order", expansion["MinimumGroupOrder"] >= -2 && zeroQ[bm2 /. masterRules]];
  Do[
    branch = thresholdDistance Coefficient[tagged, thresholdTag, rate];
    If[branch === 0, Continue[]];
    Print["  endpoint branch ", rate, " Laurent expansion ", DateString[]];
    expansion = groupedLaurent[branch, {-2, -1, 0}, {-2, -1, 0, 1}];
    {bm2, bm1, b0} = expansion["Bulk"];
    Print["  endpoint branch ", rate, " Laurent complete ", DateString[]];
    {bm2End, bm1End, b0End, b1End} = expansion["Endpoint"] /. endpointRules;
    record["real/branch-" <> ToString[rate] <> "/supported-pole-order",
      expansion["MinimumGroupOrder"] >= -2 && zeroQ[bm2 /. masterRules] && zeroQ[bm2End /. masterRules]];
    Print["  endpoint branch ", rate, " endpoint complete ", DateString[]];
    record["real/branch-" <> ToString[rate] <> "/finite-endpoint",
      FreeQ[{bm1End, b0End, b1End}, thresholdDistance | ComplexInfinity | Indeterminate | DirectedInfinity]];
    (* Subtract the endpoint test function before Laurent expansion.
       Integral_0^1 dh h^(-1-rate eps) = -1/(rate eps).
       The endpoint-removed density can be expanded as an ordinary function. *)
    pole2 += -bm1End/rate;
    pole1 += {(bm1 - bm1End)/thresholdDistance, bm1End, 0, -b0End/rate};
    finite += {(b0 - b0End - rate (bm1 - bm1End) Log[thresholdDistance])/thresholdDistance,
      b0End, -rate bm1End, -b1End/rate}, {rate, 1, maxTag}];
  <|"Pole2" -> (pole2 /. thresholdDistance -> 1 - w),
    "Pole1" -> (pole1 /. thresholdDistance -> 1 - w),
    "Finite" -> (finite /. thresholdDistance -> 1 - w)|>];

assembledResults = <||>;
If[!countertermsOnly,
  Do[
    Print["RECONSTRUCT ", classNumber, " ", spin, " real ", DateString[]];
    realDensity = 2 alphaS^3 muR^(2 eps) integrals["Real"][spin]/gs^6 /.
      {t -> -s (1 - v), u -> -s v w, s4 -> s v (1 - w)};
    realTerms = independentReal[realDensity, integrals["ThresholdValues"], integrals["MasterFunctions"]];
    virtualDensity = If[amplitudes["Spec"]["Virtual"],
      (2/(4 Pi)^(2 - eps)) (4 Pi alphaS)^3 muR^(2 eps) *
        integrals["Virtual"][spin]/gs^6/(16 Pi^2 s v) /.
        {t -> -s (1 - v), u -> -s v w, s4 -> s v (1 - w)}, 0];
    virtualExpansion = Series[virtualDensity, {eps, 0, 0}];
    virtualLeadingOrder = If[Head[virtualExpansion] === SeriesData,
      virtualExpansion[[4]]/virtualExpansion[[6]], If[virtualExpansion === 0, Infinity, 0]];
    record[spin <> "/virtual-supported-pole-order", TrueQ[virtualLeadingOrder >= -2]];
    virtualTerms = SeriesCoefficient[virtualExpansion, {eps, 0, #}] & /@ {-2, -1, 0};
    Print["RECONSTRUCT ", classNumber, " ", spin, " real complete ", DateString[]];
    reconstructed = (realTerms["Finite"] + {0, 0, 0, virtualTerms[[3]]} +
      Total[#["Finite"] & /@ Values[ctResults[spin]]]) /. integrals["MasterFunctions"];
    poleOne = (realTerms["Pole1"] + {0, 0, 0, virtualTerms[[2]]} +
      Total[#["Pole"] & /@ Values[ctResults[spin]]]) /. integrals["MasterFunctions"];
    record[spin <> "/pole2", zeroQ[realTerms["Pole2"] + virtualTerms[[1]]]];
    Do[
      record[spin <> "/pole1/" <> ToString[part], zeroQ[poleOne[[part]]]];
      record[spin <> "/real/" <> ToString[part], zeroQ[(realTerms["Finite"][[part]] -
        target["Real"][spin][[part]]) /. integrals["MasterFunctions"]]];
      record[spin <> "/finite/" <> ToString[part], finiteSlotEqualQ[reconstructed, target["Finite"][spin], part]];
      mutatedTarget = ReplacePart[target["Finite"][spin], part -> target["Finite"][spin][[part]] +
        alphaS^3/(Pi s v) If[part === 1, v (1 - v) w (1 - w), v (1 - v)]];
      record[spin <> "/control/finite-slot-" <> ToString[part],
        !finiteSlotEqualQ[reconstructed, mutatedTarget, part]],
      {part, 4}];
    record[spin <> "/virtual", zeroQ[virtualTerms[[3]] - target["Virtual"][spin]]];
    assembledResults[spin] = <|"Real" -> realTerms, "Virtual" -> virtualTerms,
      "Finite" -> reconstructed, "Pole1" -> poleOne|>;
    If[physicalClassNumber === "1",
      omittedBornTerms = Block[{amplitudes = amplitudes /. D -> 4}, independentCounterterms[spin]];
      omittedMeasureTerms = independentCounterterms[spin, 0];
      Do[omission = {omittedBornTerms, omittedMeasureTerms}[[omissionIndex]];
        omittedCorrection = Total[#["Finite"] & /@ Values[omission]] -
        Total[#["Finite"] & /@ Values[ctResults[spin]]];
        record[spin <> "/control/actual-omission-" <> ToString[omissionIndex],
          !tupleEqualQ[reconstructed + omittedCorrection, target["Finite"][spin]]],
        {omissionIndex, 2}]];
    , {spin, {"XX", "NN"}}]
];

(* The finite transfer tensor is isotropic in the two physical transverse
   axes after colour identities, renormalization and collinear subtraction.
   Test the independently assembled components, not their displayed target. *)
If[!countertermsOnly, Do[
  record["tensor/XX-equals-NN/" <> ToString[part],
    finiteSlotEqualQ[assembledResults["XX"]["Finite"], assembledResults["NN"]["Finite"], part]],
  {part, 4}]];

(* Sensitivity controls pass the perturbed tuple through the same comparator.
   A finite regular perturbation vanishes at both endpoints and changes no scale
   derivative, so universal pole/scale constraints cannot detect it alone. *)
Do[injected = ReplacePart[zeroTuple, part -> If[part === 1, v (1 - v) w (1 - w), v (1 - v)]];
  record["control/finite-slot-" <> ToString[part], !tupleEqualQ[injected, zeroTuple]], {part, 4}];
record["control/delta-normalization", !zeroQ[alphaS^3/(Pi s v)]];
referenceCounterterm = ctResults["XX"]["incoming a"]["Finite"];
record["control/finite-only-counterterm-constant", !tupleEqualQ[referenceCounterterm,
  referenceCounterterm + {0, 0, 0, alphaS^3 CF/(Pi s v)}]];
evanescentEffect = subtraction[pullback[transverseKernel, w, (1 + eps CF)/w], muF]["Finite"] -
  subtraction[pullback[transverseKernel, w, 1/w], muF]["Finite"];
record["control/evanescent-Born", !tupleEqualQ[evanescentEffect, zeroTuple]];
fragmentationMap = 1 - v + v w;
fragmentationEffect = subtraction[pullback[transverseKernel, fragmentationMap,
  (1 + w) fragmentationMap^(2 eps)], muFF]["Finite"] -
  subtraction[pullback[transverseKernel, fragmentationMap, 1 + w], muFF]["Finite"];
record["control/fragmentation-dimension", !tupleEqualQ[fragmentationEffect, zeroTuple]];

(* Direct integration of test functions in the original z variable gives a
   second check of endpoint terms and Jacobians in the generic pullback. *)
If[physicalClassNumber === "1", Do[
  zMap = {w, (1 - v)/(1 - v w), 1 - v + v w}[[leg]] /. v -> 2/5;
  inverseMap = w /. First[Solve[zMap == z, w]];
  lowerZ = zMap /. w -> 0;
  unitWeight = 1 + w + w^2;
  testFunction = 1 + w + w^3;
  unitKernel = {-2, 2, 3/2};
  transformed = pullback[unitKernel, zMap, unitWeight];
  fOriginal = Cancel[(unitWeight testFunction/D[zMap, w]) /. w -> inverseMap];
  fOne = fOriginal /. z -> 1;
  directIntegral = Integrate[unitKernel[[1]] fOriginal +
      unitKernel[[2]] Cancel[(fOriginal - fOne)/(1 - z)], {z, lowerZ, 1}] +
    fOne (unitKernel[[3]] + unitKernel[[2]] Log[1 - lowerZ]);
  transformedIntegral = Integrate[transformed[[1]] testFunction +
    transformed[[2]] Cancel[(testFunction - (testFunction /. w -> 1))/(1 - w)], {w, 0, 1}] +
    transformed[[4]] (testFunction /. w -> 1);
  record["test-function/pullback-leg-" <> ToString[leg], zeroQ[directIntegral - transformedIntegral]], {leg, 3}]];

record["input-integrity", inputHashes === AssociationThread[inputFiles,
  IntegerString[FileHash[#, "SHA256"], 16, 64] & /@ inputFiles]];
If[!countertermsOnly, record["full-assembly-both-components", Sort[Keys[assembledResults]] === {"NN", "XX"}]];
checkerHash = IntegerString[FileHash[$InputFileName, "SHA256"], 16, 64];
Put[<|"Class" -> classNumber, "Counterterms" -> ctResults, "Assembly" -> assembledResults, "Checks" -> checks,
  "CheckerSHA256" -> checkerHash,
  "Failures" -> failures, "Inputs" -> inputHashes, "Seconds" -> AbsoluteTime[] - started|>,
  FileNameJoin[{outputDirectory, "independent.wl"}]];
Export[FileNameJoin[{outputDirectory, "summary.json"}],
  <|"class" -> classNumber, "mode" -> If[countertermsOnly, "counterterms-only", "full-assembly"],
    "assemblyComponents" -> Keys[assembledResults], "passed" -> (failures === {}), "checkCount" -> Length[checks],
    "failures" -> failures, "seconds" -> N[AbsoluteTime[] - started], "inputs" -> inputHashes,
    "checkerSHA256" -> checkerHash|>, "RawJSON"];
Print["INDEPENDENT ASSEMBLY class=", classNumber, " checks=", Length[checks], " failures=", Length[failures]];
Quit[If[failures === {}, 0, 1]];
