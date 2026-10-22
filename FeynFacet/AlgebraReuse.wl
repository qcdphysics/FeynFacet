(* Small linear-algebra helpers for explicit spin projections.

   Scalar coefficients are combined only when the remaining indexed structure
   is identical. No gamma-matrix identity, tensor basis, or spin projector is
   introduced here. The caller supplies the original evaluation operation.

   Cache keys must include the kinematic conventions and every parameter of
   that operation. Each Wolfram worker owns its own bounded cache. *)
BeginPackage["FeynFacet`"];

GroupIdenticalAlgebraTerms::usage = "GroupIdenticalAlgebraTerms[expr, pattern] returns scalar coefficients and identical structures containing pattern. Each scalar coefficient is free of pattern.";
CreateBoundedAlgebraCache::usage = "CreateBoundedAlgebraCache[evaluate, entries, bytes] returns Apply[key,input], Statistics[], and Clear[] functions. Cache keys are compared exactly; nothing is persisted to disk.";

Begin["`AlgebraReusePrivate`"];

GroupIdenticalAlgebraTerms[expr_, pattern_] := Module[
  {tag = Unique["algebraGrouping"], terms, separated, groups},
  Catch[
    If[!FreeQ[expr, _Failure | $Failed | $Aborted],
      Throw[Failure["AlgebraInput", <|"Message" -> "The input contains a failed calculation."|>], tag]];
    terms = Expand[expr];
    terms = If[Head[terms] === Plus, List @@ terms, {terms}];
    separated = Table[Module[{factors, indexed, scalar},
      factors = If[Head[term] === Times, List @@ term, {term}];
      indexed = Times @@ Select[factors, !FreeQ[#, pattern] &];
      scalar = Times @@ Select[factors, FreeQ[#, pattern] &];
      If[!FreeQ[scalar, pattern] || scalar indexed =!= term,
        Throw[Failure["AlgebraSplit", <|"Term" -> HoldForm[term]|>], tag]];
      {indexed, scalar}], {term, terms}];
    groups = GatherBy[separated, First];
    DeleteCases[Table[Module[{coefficient},
      coefficient = Total[Last /@ group];
      If[Length[group] > 1, coefficient = Together[coefficient]];
      If[coefficient === 0, Nothing,
        <|"Structure" -> First[First[group]], "Coefficient" -> coefficient|>]],
      {group, groups}], Nothing],
    tag]
];

CreateBoundedAlgebraCache[evaluate_, entries_Integer : 256, bytes_Integer : 134217728] := Module[
  {cache = <||>, usedBytes = 0, calls = 0, hits = 0, resets = 0,
   tag = Unique["algebraCache"], apply, statistics, clear},
  If[entries < 1 || bytes < 1,
    Return[Failure["AlgebraCacheBounds", <|"Entries" -> entries, "Bytes" -> bytes|>]]];
  apply = Function[{key, input}, Catch[Module[{value, size},
    If[!FreeQ[{key, input}, _Failure | $Failed | $Aborted],
      Throw[Failure["AlgebraCacheInput", <|"Message" -> "A failed calculation cannot be cached."|>], tag]];
    calls++;
    If[KeyExistsQ[cache, key],
      hits++; cache[key],
      value = evaluate[input];
      If[!FreeQ[value, _Failure | $Failed | $Aborted],
        Throw[Failure["AlgebraEvaluation", <|"Value" -> value|>], tag]];
      size = ByteCount[{key, value}];
      If[size <= bytes,
        If[Length[cache] >= entries || usedBytes + size > bytes,
          cache = <||>; usedBytes = 0; resets++];
        cache[key] = value; usedBytes += size];
      value]], tag]];
  statistics = Function[Null, <|"Calls" -> calls, "Hits" -> hits,
    "Evaluations" -> calls - hits, "Entries" -> Length[cache],
    "Bytes" -> usedBytes, "Resets" -> resets|>];
  clear = Function[Null, cache = <||>; usedBytes = 0; Null];
  <|"Apply" -> apply, "Statistics" -> statistics, "Clear" -> clear|>
];

End[];
EndPackage[];
