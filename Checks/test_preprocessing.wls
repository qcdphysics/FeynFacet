(* Small exact tests of the two pre-Kira representation changes. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "FeynFacet", "Preprocessing.wl"}]];

passed = 0;
check[label_, condition_] := If[TrueQ[condition], passed++; Print["PASS ", label],
  Print["FAIL ", label]; Exit[1]];
variables = {fA, fB, fC, fAp, fBp, fCp, spe3};
check["zero has no contribution terms", CombineLaurentTerms[0, variables] === {}];
check["scalar coefficients cancel before row conversion",
  CombineLaurentTerms[(1/(s - t) - (s + t)/(s^2 - t^2)) spe3/fA^2, variables] === {}];
combined = CombineLaurentTerms[c1 spe3 fB/fA^2 + c2 spe3 fB/fA^2, variables];
check["equal negative-power monomials are combined", Length[combined] === 1 &&
  Together[First[combined] - (c1 + c2) spe3 fB/fA^2] === 0];
check["distinct hat powers remain distinct",
  Length[CombineLaurentTerms[c1/fA + c2 spe3/fA, variables]] === 2];
check["trace and integration dimensions remain scalar coefficients",
  CombineLaurentTerms[D spe3/fA + d spe3/fA, variables] === {(d + D) spe3/fA}];
check["fractional powers are rejected", FailureQ[CombineLaurentTerms[Sqrt[fA], variables]]];
check["sums in a Laurent denominator are rejected", FailureQ[CombineLaurentTerms[1/(fA + fB), variables]]];
check["functional dependence is rejected", FailureQ[CombineLaurentTerms[Sin[fA], variables]]];
check["Failure input never enters grouping", FailureQ[CombineLaurentTerms[Failure["test", <||>], variables]]];
check["nonscalar input is rejected", FailureQ[CombineLaurentTerms[{fA}, variables]]];

calls = 0;
formOf[factor_] := (calls++; Which[
  Expand[factor - (a - b)] === 0, fA,
  Expand[factor + (a - b)] === 0, -fA,
  True, Failure["UnknownForm", <||>]]);
expr = c1/((s + t) (a - b)^2) + c2/((s + t) (a - b)^2);
mapped = MapIntegrationDenominators[expr, {}, {a, b}, formOf];
check["each distinct denominator is mapped once", calls === 1];
check["independent denominator factors and the integrand are preserved",
  Together[(mapped /. fA -> a - b) - expr] === 0];
mapped = MapIntegrationDenominators[c1/(b - a) + c2/(a - b)^2, {}, {a, b}, formOf];
check["denominator signs are preserved",
  Together[(mapped /. fA -> a - b) - (c1/(b - a) + c2/(a - b)^2)] === 0];
check["invariant substitutions precede mapping",
  Together[(MapIntegrationDenominators[c1/u, {u -> a - b}, {a, b}, formOf] /. fA -> a - b) - c1/(a - b)] === 0];
check["zero denominator conversion", MapIntegrationDenominators[0, {}, {a, b}, formOf] === 0];
check["failed nested mapping exits the whole conversion",
  FailureQ[MapIntegrationDenominators[c1/(a + b) + c2/(a - b), {}, {a, b}, formOf]]];
check["failed denominator conversion is rejected by grouping",
  FailureQ[CombineLaurentTerms[MapIntegrationDenominators[1/(a + b), {}, {a, b}, formOf], variables]]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "FeynFacet", "Tools.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "FeynFacet", "Kira.wl"}]];
check["known masters are accepted", AssociationQ[CheckReductionCoverage[{k[1]}, {}, {k}, {k[1]}]]];
check["missing targets are rejected before relabelling",
  CheckAbort[CheckReductionCoverage[{k[2]}, {}, {k}, {k[1]}], $Aborted] === $Aborted];
check["unknown residuals are rejected before relabelling",
  CheckAbort[CheckReductionCoverage[{k[2]}, {k[2] -> k[3]}, {k}, {k[1]}], $Aborted] === $Aborted];
Print[passed, " preprocessing tests passed."];
Exit[0];
