(* Exact regression tests for the NLO distribution expansion.
   The reference expands a whole branch directly; the production routine
   groups repeated epsilon-dependent structures before expanding them. *)
$HistoryLength = 0;
root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{root, "FeynFacet", "FeynFacet.wl"}]];

CheckAbort[
passed = 0;
check[name_, condition_] := If[TrueQ[condition], passed++; Print["PASS ", name],
  Print["FAIL ", name]; Exit[1]];
region = 0 < v < 1 && 0 < w < 1;
slots[result_] := Join[{result["Pole2"]}, result["Pole1"], result["Finite"]];
equal[left_, right_] := And @@ (Together[#] === 0 & /@ (slots[left] - slots[right]));

reference[A_, branches_Association] := Module[{endpoint, pieces, series, b, ends, x},
  endpoint[expression_] := Together[Normal[Series[expression /. w -> 1 - x, {x, 0, 0}]]];
  pieces = KeyValueMap[Function[{lambda, B},
    series = Series[B, {eps, 0, 1}];
    b = Table[SeriesCoefficient[series, {eps, 0, n}], {n, -1, 1}];
    ends = endpoint /@ b;
    <|"Pole2" -> -ends[[1]]/lambda,
      "Pole1" -> {(b[[1]] - ends[[1]])/(1 - w), ends[[1]], -ends[[2]]/lambda},
      "Finite" -> {(b[[2]] - ends[[2]])/(1 - w) - lambda (b[[1]] - ends[[1]]) Log[1 - w]/(1 - w),
        ends[[2]], -lambda ends[[1]], -ends[[3]]/lambda}|>], branches];
  <|"Pole2" -> Total[#["Pole2"] & /@ pieces],
    "Pole1" -> {SeriesCoefficient[A, {eps, 0, -1}], 0, 0} + Total[#["Pole1"] & /@ pieces],
    "Finite" -> {SeriesCoefficient[A, {eps, 0, 0}], 0, 0, 0} + Total[#["Finite"] & /@ pieces]|>];

compare[name_, A_, branches_] := Module[{result, expected},
  result = EndpointExpand[A, branches, eps, w, region];
  expected = reference[A, branches];
  If[!equal[result, expected], Print["Coefficient differences: ", InputForm[slots[result] - slots[expected]]]];
  check[name, equal[result, expected]]];

zero = EndpointExpand[0, <||>, eps, w, region];
check["zero and empty branches", slots[zero] === ConstantArray[0, 8]];
check["distribution basis is retained", zero["Basis"] === {"regular", "D0", "D1", "delta"}];
compare["regular finite expression", a + b w, <||>];
compare["regular single pole", a/eps + b + eps c, <||>];
compare["constant branch", 0, <|1 -> 1|>];
compare["single epsilon pole in a branch", 0, <|2 -> 1/eps|>];
compare["epsilon term supplies a finite delta coefficient", 0, <|2 -> eps (1 + w)|>];
check["finite delta coefficient has the distribution sign",
  EndpointExpand[0, <|2 -> eps (1 + w)|>, eps, w]["Finite"] === {0, 0, 0, -1}];
compare["w-dependent coefficients retain the regular remainder", (1 + v w)/eps,
  <|1 -> (a + b w)/eps + c w^2 + eps d (1 - w), 2 -> v w/eps + eps w|>];

gammaFactor = Gamma[1 - eps]^2/Gamma[1 - 2 eps];
compare["shared Gamma factors and dimensional powers", 0,
  <|1 -> gammaFactor (v^(1 - eps)/(eps (1 + eps)) +
    w v^(-1 - eps)/(eps (1 + eps))),
    2 -> gammaFactor (1 + w) v^(-eps)/eps|>];
compare["mixed polynomial denominators are factored exactly", 0,
  <|1 -> (1 + w)/(eps - eps v + eps^2 - eps^2 v)|>];
compare["negative denominator signs are retained", 0,
  <|1 -> (1 + w)/(eps (v - 1) (1 - 2 eps))|>];
compare["rational scalar cancellation", 0,
  <|1 -> (1/(s - t) - (s + t)/(s^2 - t^2)) gammaFactor/eps|>];
compare["intermediate deeper epsilon poles cancel", 0,
  <|1 -> ((1 + eps)^(-1) - 1 + eps)/eps^2|>];
compare["removable logarithmic endpoint", 0,
  <|1 -> Log[w]/((1 - w) eps) + eps w|>];
compare["endpoint poles cancel between distinct structures", 0,
  <|1 -> (Log[w]/(1 - w)^2 + 1/(1 - w))/eps|>];
compare["scalar master functions retain their endpoint values", 0,
  <|1 -> master[w]/eps + eps master[w]^2|>];
compare["trace and integration dimensions remain distinct", 0,
  <|1 -> (D + d w) gammaFactor/eps|>];
compare["colour and flavour coefficients remain symbolic", 0,
  <|1 -> (CF CA + Nf w) gammaFactor/eps, 2 -> CF^2 v^(-eps)/eps|>];
compare["nonlinear analytic dimensional exponent", 0,
  <|1 -> v^(1 - eps + eps^2)/eps|>];

first = EndpointExpand[0, <|1 -> v^(-eps)/eps|>, eps, w];
second = EndpointExpand[0, <|1 -> w^(-eps)/eps|>, eps, w];
check["caches do not leak between calls", equal[first, reference[0, <|1 -> v^(-eps)/eps|>]] &&
  equal[second, reference[0, <|1 -> w^(-eps)/eps|>]]];
check["single-branch interface is retained",
  EndpointExpand[0, (1 + w)/eps, 2, eps, w] ===
    EndpointExpand[0, <|2 -> (1 + w)/eps|>, eps, w]];
check["epsilon symbol is supplied by the caller",
  EndpointExpand[0, <|1 -> (1 + w)/eta + eta v|>, eta, w] ===
    EndpointExpand[0, <|1 -> (1 + w)/eps + eps v|>, eps, w]];

rejected[name_, expression_] := check[name, CheckAbort[expression, $Aborted] === $Aborted];
SetAttributes[rejected, HoldRest];
rejected["a genuine endpoint pole is rejected", EndpointExpand[0, <|1 -> 1/(1 - w)|>, eps, w]];
rejected["an uncancelled deeper epsilon pole is rejected", EndpointExpand[0, <|1 -> eps^-2|>, eps, w]];
rejected["a deeper pole in the regular part is rejected", EndpointExpand[eps^-2, <||>, eps, w]];
rejected["non-Laurent epsilon logarithms are rejected", EndpointExpand[0, <|1 -> Log[eps]|>, eps, w]];
rejected["fractional Laurent powers are rejected", EndpointExpand[0, <|1 -> Sqrt[eps]|>, eps, w]];
rejected["failed input aborts the expansion", EndpointExpand[0, <|1 -> Failure["test", <||>]|>, eps, w]];
rejected["nonscalar input is rejected", EndpointExpand[{a, b}, <||>, eps, w]];
rejected["zero branch exponent is rejected", EndpointExpand[0, <|0 -> 1|>, eps, w]];

Print[passed, " distribution tests passed."];
Exit[0], Exit[1]];
