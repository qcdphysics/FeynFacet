(* ::Package:: *)

(* Keep additive fractions separate until the integration variables have been
   replaced by inert forms. Combine equal Laurent monomials before assigning
   integral families. Neither operation changes a spin or dimension convention. *)

BeginPackage["FeynFacet`"];

MapIntegrationDenominators::usage = "MapIntegrationDenominators[expr, rules, variables, formOf] maps each distinct term denominator once, after applying rules. Factors independent of variables are retained. Unsupported input returns Failure.";

CombineLaurentTerms::usage = "CombineLaurentTerms[expr, variables] returns a list of terms with equal integer Laurent powers combined. Together acts only on the scalar coefficients; exact zeros are removed. Unsupported variable dependence returns Failure.";


CreatePartialFractionConverter::usage = "CreatePartialFractionConverter[forms, hat, transform] returns a cached converter. Its input is the mapped rational integrand; its output contains grouped Terms and structural counts.";
CreateIntegralRowConverter::usage = "CreateIntegralRowConverter[variables, mapUnit, reconstruct] returns a cached converter for grouped terms. mapUnit uses the existing family mapping; reconstruct checks each scaled row list exactly.";

Begin["`Private`"];

additiveTerms[expr_] := Which[expr === 0, {}, Head[expr] === Plus, List @@ expr, True, {expr}];

MapIntegrationDenominators[expr_, rules_List, variables : {__Symbol}, formOf_] :=
 Module[{tag = Unique["denominatorFailure"], cache = <||>, pattern,
   reject, mappedDenominator, expanded},
  reject[reason_, detail_] := Throw[Failure["UnsupportedDenominator", <|
    "Reason" -> reason, "Expression" -> HoldForm[detail]|>], tag];
  Catch[
   If[!DuplicateFreeQ[variables] || !FreeQ[expr, _Failure | _List | _Association | _Rule | _RuleDelayed],
    reject["Expected a scalar expression and distinct integration variables", expr]];
   pattern = Alternatives @@ variables;
   mappedDenominator[den_] := Module[{factors, mapped},
    If[KeyExistsQ[cache, den], Return[cache[den]]];
    factors = Quiet[Check[FactorList[den], $Failed]];
    If[!MatchQ[factors, {{_, _Integer} ...}], reject["Denominator could not be factored", den]];
    mapped = Times @@ Map[Function[item, Module[{factor = item[[1]], power = item[[2]], image},
      image = If[FreeQ[factor, pattern], factor, formOf[factor]];
      If[FailureQ[image] || image === $Failed || !FreeQ[image, pattern],
       reject["Integration-dependent factor was not mapped", factor]];
      image^power]], factors];
    cache[den] = mapped];
   expanded = Expand[expr /. rules];
   Total[Map[Function[term, Numerator[term]/mappedDenominator[Denominator[term]]],
     additiveTerms[expanded]]], tag]];

CombineLaurentTerms[expr_, variables : {__Symbol}] :=
 Module[{tag = Unique["laurentFailure"], pattern, reject, rows, groups},
  reject[reason_, detail_] := Throw[Failure["UnsupportedLaurentTerm", <|
    "Reason" -> reason, "Expression" -> HoldForm[detail]|>], tag];
  Catch[
   If[!DuplicateFreeQ[variables] || !FreeQ[expr, _Failure | _List | _Association | _Rule | _RuleDelayed],
    reject["Expected a scalar expression and distinct monomial variables", expr]];
   pattern = Alternatives @@ variables;
   rows = Map[Function[term, Module[{powers, monomial, coefficient},
      powers = Exponent[term, #] & /@ variables;
      If[!VectorQ[powers, IntegerQ], reject["Laurent powers must be integers", term]];
      monomial = Times @@ MapThread[Power, {variables, powers}];
      coefficient = term/monomial;
      If[!FreeQ[coefficient, pattern] || !FreeQ[coefficient, _Failure],
       reject["The scalar coefficient still contains a monomial variable", term]];
      {powers, coefficient}]], additiveTerms[Expand[expr]]];
   groups = GroupBy[rows, First -> Last];
   KeyValueMap[Function[{powers, coefficients}, With[{coefficient = Together[Total[coefficients]]},
      If[coefficient === 0, Nothing,
       coefficient Times @@ MapThread[Power, {variables, powers}]]]], groups], tag]];

MapIntegrationDenominators[___] := Failure["InvalidPreprocessingArguments", <||>];
CombineLaurentTerms[___] := Failure["InvalidPreprocessingArguments", <||>];



(* The input terms must already be separated. Expanding here would distribute
   their potentially large scalar coefficients and undo the grouping. *)
monomialParts[term_, variables_] := Module[{powers, monomial, coefficient},
 powers = Exponent[term, #] & /@ variables;
 If[!VectorQ[powers, IntegerQ], Return[Failure["NonintegerLaurentPower", <||>]]];
 monomial = Times @@ MapThread[Power, {variables, powers}];
 coefficient = term/monomial;
 If[!FreeQ[coefficient, Alternatives @@ variables] || !FreeQ[coefficient, _Failure],
  Return[Failure["NonScalarCoefficient", <||>]]];
 {powers, coefficient}];

groupPairs[pairs_, variables_] := KeyValueMap[
 Function[{powers, coefficients}, With[{coefficient = Together[Total[coefficients]]},
  If[coefficient === 0, Nothing, coefficient Times @@ MapThread[Power, {variables, powers}]]]],
 GroupBy[pairs, First -> Last]];

CreatePartialFractionConverter[forms : {__Symbol}, hat_Symbol, transform_] :=
 Module[{cache = <||>, tag = Unique["sparsePartialFraction"], variables = Append[forms, hat]},
 Function[integrand, Module[{reject, terms, pieces, pairs = {}, mapped, powers, coefficient, output, hits = 0, before},
  reject[reason_, detail_] := Throw[Failure[reason, <|"Detail" -> HoldForm[detail]|>], tag];
  Catch[
   If[!DuplicateFreeQ[variables], reject["RepeatedVariable", variables]];
   terms = FeynFacet`CombineLaurentTerms[integrand, variables];
   If[!ListQ[terms] || !FreeQ[terms, _Failure], reject["InputGroupingFailed", terms]];
   before = Length[cache];
   Do[
    pieces = monomialParts[term, variables];
    If[FailureQ[pieces], reject["InvalidInputMonomial", pieces]];
    {powers, coefficient} = pieces;
    If[KeyExistsQ[cache, Most[powers]], hits++,
     mapped = transform[Times @@ MapThread[Power, {forms, Most[powers]}]];
     If[!FreeQ[mapped, _Failure | $Failed | $Aborted] || !FreeQ[mapped, hat],
      reject["GeometricPartialFractionFailed", mapped]];
     mapped = FeynFacet`CombineLaurentTerms[mapped, forms];
     If[!ListQ[mapped] || !FreeQ[mapped, _Failure], reject["GeometricGroupingFailed", mapped]];
     mapped = monomialParts[#, forms] & /@ mapped;
     If[!FreeQ[mapped, _Failure], reject["InvalidGeometricMonomial", mapped]];
     cache[Most[powers]] = mapped];
    Scan[Function[row, AppendTo[pairs, {Append[row[[1]], Last[powers]], coefficient row[[2]]}]], cache[Most[powers]]],
    {term, terms}];
   output = groupPairs[pairs, variables];
   <|"Terms" -> output, "InputStructures" -> Length[terms], "OutputStructures" -> Length[output],
     "NewGeometricMaps" -> Length[cache] - before, "CacheHits" -> hits|>, tag]]]];

CreateIntegralRowConverter[variables : {__Symbol}, mapUnit_, reconstruct_] :=
 Module[{cache = <||>, tag = Unique["sparseIntegralRows"]},
 Function[terms, Module[{reject, pieces, powers, coefficient, unit, rows, output = {}, before, hits = 0},
  reject[reason_, detail_] := Throw[Failure[reason, <|"Detail" -> HoldForm[detail]|>], tag];
  Catch[
   If[!ListQ[terms] || !FreeQ[terms, _Failure], reject["ExpectedGroupedTerms", terms]];
   If[!DuplicateFreeQ[variables], reject["RepeatedVariable", variables]];
   before = Length[cache];
   Do[
    pieces = monomialParts[term, variables];
    If[FailureQ[pieces], reject["InvalidInputMonomial", pieces]];
    {powers, coefficient} = pieces;
    If[KeyExistsQ[cache, powers], hits++,
     unit = Times @@ MapThread[Power, {variables, powers}];
     rows = mapUnit[unit];
     If[!ListQ[rows] || !AllTrue[rows, AssociationQ] || !FreeQ[rows, _Failure],
      reject["UnitRowConversionFailed", rows]];
     cache[powers] = rows];
    rows = Map[Function[row, Join[row, <|"Coefficient" -> Factor[coefficient row["Coefficient"]]|>]], cache[powers]];
    If[!TrueQ[reconstruct[term, rows]], reject["ScaledRowsDoNotReconstruct", term]];
    AppendTo[output, rows], {term, terms}];
   <|"Rows" -> Flatten[output, 1], "NewUnitMaps" -> Length[cache] - before, "CacheHits" -> hits|>, tag]]]];

CreatePartialFractionConverter[___] := Failure["InvalidPartialFractionArguments", <||>];
CreateIntegralRowConverter[___] := Failure["InvalidRowConverterArguments", <||>];

End[];
EndPackage[];
