(* ::Package:: *)

(* Distributions.wl: expansion of the real-emission threshold branch
   (1-w)^(-1-lambda eps) into delta(1-w), [1/(1-w)]_+ and [log(1-w)/(1-w)]_+ . *)

BeginPackage["FeynFacet`"];

EndpointExpand::usage = "EndpointExpand[A, branches, eps, w, assumptions] expands \
an integrated real-emission result written as\n\
    A(w,eps) + Sum_lambda (1-w)^(-1-lambda eps) B_lambda(w,eps),\n\
branches being <|lambda -> B_lambda, ...|> (the output of ThresholdSplit), with \
A and B_lambda regular at w = 1 and B = b[-1]/eps + b[0] + eps b[1] + ..., using\n\
    (1-w)^(-1-lambda eps) = -delta(1-w)/(lambda eps) + D0 - lambda eps D1 + ...,\n\
D0 = [1/(1-w)]_+ and D1 = [log(1-w)/(1-w)]_+.  A w-dependent coefficient of D0 or \
D1 is split into its value at w=1 times the distribution plus an ordinary \
remainder.  The result is an Association with keys\n\
  \"Pole2\"  -> coefficient of delta(1-w)/eps^2,\n\
  \"Pole1\"  -> {ordinary, D0, delta} coefficients of 1/eps,\n\
  \"Finite\" -> {ordinary, D0, D1, delta}.\n\
EndpointExpand[A, B, lambda, eps, w, assumptions] is the single-branch form.";

ThresholdSplit::usage = "ThresholdSplit[expr, w, eps] writes an integrated real-emission \
result, in which the threshold variable 1-w appears through rational functions and \
through powers (1-w)^(a - lambda eps) coming from the master integrals, as\n\
    A(w,eps) + Sum_lambda (1-w)^(-1-lambda eps) B_lambda(w,eps),\n\
and returns <|\"Regular\" -> A, \"Branches\" -> <|lambda -> B_lambda, ...|>|>.  \
Every base of a non-integer power is factored into factors that are positive in \
the physical region (an interior point can be given as the fourth argument, \
default s = 3, v = 1/3, w = 1/2), the terms are sorted by lambda, and B_lambda is \
(1-w) times the coefficient of (1-w)^(-lambda eps).  The function aborts if A or \
a B_lambda is singular at w = 1, i.e. if the expression is more singular than \
1/(1-w) in any branch.  The branches must be separated beforehand: a master \
integral is to be written as a sum of terms each with a single power of its \
threshold variable (see the connection formulas in the process notebooks).";

CanonicalDistribution::usage = "CanonicalDistribution[prefactor, {R, B, C, A}, w] \
converts a hard part written as prefactor(w) [R(w) + B(w) D0 + C(w) D1 + A delta(1-w)], \
the form used in the literature, into the basis {regular, D0, D1, delta} of this \
package, in which the coefficients of D0, D1 and delta(1-w) are numbers at w = 1 \
and everything else is an ordinary function.  It uses the identities\n\
   f(w) D0 = f(1) D0 + (f(w) - f(1))/(1-w),   f(w) D1 = f(1) D1 + (f(w) - f(1)) log(1-w)/(1-w),\n\
valid for any function f that is regular at w = 1 (here f = prefactor B etc.).  f(1) is the limit w -> 1, \
so a removable singularity at w = 1 (e.g. B = log(w)/(w-1), f(1) = 1) is allowed; a coefficient that is \
singular at w = 1 stops the calculation.";

ThresholdBasis::usage = "ThresholdBasis is the ordered basis {\"regular\", \"D0\", \"D1\", \"delta\"} \
used for every finite coefficient in this package.";

Begin["`Private`"];

ThresholdBasis = {"regular", "D0", "D1", "delta"};

coefficient[expr_, eps_, n_] := SeriesCoefficient[expr, {eps, 0, n}];
(* value at w = 1, through a series expansion so that removable singularities
   such as log(z)/z at z -> 0 are handled *)
atOne[expr_, w_, assumptions_] := Module[{wb, v},
  v = Normal[Series[expr /. w -> 1 - wb, {wb, 0, 0}]];
  Require[FreeQ[v, wb | ComplexInfinity | Indeterminate | DirectedInfinity],
    "Endpoint coefficient is singular at w = 1", expr];
  Together[v]];

EndpointExpand[A_, branches_Association, eps_Symbol, w_Symbol, assumptions_: True] :=
  Module[{denominators = <||>, epsilonSeries = <||>, endpointSeries = <||>,
    factorDenominator, expandStructure, laurent, endpoint, branch, parts, regular, wb},

    Require[AllTrue[Keys[branches], IntegerQ[#] && # > 0 &],
      "Endpoint branches require positive integer lambda", Keys[branches]];

    (* Separate each branch into C(w) F(w,eps), with C free of eps.  Expanding
       a whole branch repeats the same Gamma functions and dimensional powers
       for many rational coefficients.  Factor polynomial denominators first
       to expose their eps-only factors, and keep powers such as x^(n+a eps)
       as x^n dimensionPower[x,a eps] until the distinct F's have been grouped.
       No product of power bases is split here: ThresholdSplit supplies the
       physical-region prescription for those bases.  All caches are local
       to this call, so another spin projection or kinematic convention starts
       with an empty cache. *)
    factorDenominator[base_] := If[KeyExistsQ[denominators, base], denominators[base],
      denominators[base] = Factor[base]];

    expandStructure[structure_, order_] := Module[{key = {structure, order}, series, terms, pairs},
      If[!KeyExistsQ[epsilonSeries, key],
        series = Expand[Normal[Series[structure /. dimensionPower -> Power, {eps, 0, order}]]];
        Require[FreeQ[series, _Series | _SeriesData | _Failure | $Failed | $Aborted |
          ComplexInfinity | Indeterminate | DirectedInfinity],
          "Could not obtain a Laurent expansion in epsilon", structure];
        terms = If[series === 0, {}, If[Head[series] === Plus, List @@ series, {series}]];
        pairs = Map[Function[term, Module[{n = Exponent[term, eps], scalar},
          Require[IntegerQ[n], "Noninteger Laurent power in epsilon", term];
          scalar = term/eps^n;
          Require[FreeQ[scalar, eps], "Non-Laurent epsilon dependence", term];
          {n, scalar}]], terms];
        epsilonSeries[key] = Map[Total, GroupBy[pairs, First -> Last]]];
      epsilonSeries[key]];

    laurent[expression_, order_] := Module[{normalized, rows, expansions, powers, value},
      Require[!ListQ[expression] && !AssociationQ[expression],
        "EndpointExpand requires scalar expressions", expression];
      normalized = expression /. Power[base_, n_Integer] /;
        n < 0 && !FreeQ[base, eps] && PolynomialQ[base, eps] :> factorDenominator[base]^n;
      normalized = normalized /. Power[base_, power_] /;
        FreeQ[base, eps] && !FreeQ[power, eps] && PolynomialQ[power, eps] :>
          base^(power /. eps -> 0) dimensionPower[base, Expand[power - (power /. eps -> 0)]];
      rows = GroupIdenticalAlgebraTerms[normalized, eps];
      Require[ListQ[rows], "Epsilon-structure grouping failed", rows];
      expansions = expandStructure[#["Structure"], order] & /@ rows;
      powers = Union[Flatten[Keys /@ expansions]];
      value[n_] := Total[MapThread[#1["Coefficient"] Lookup[#2, n, 0] &, {rows, expansions}]];
      (* Individual structures may have deeper poles which cancel in their
         sum.  Check that cancellation before extracting the NLO coefficients. *)
      Do[Require[Together[value[n]] === 0, "Unsupported epsilon pole beyond the NLO branch form", n],
        {n, Select[powers, # < -1 &]}];
      Table[value[n], {n, -1, order}]];

    endpoint[expression_] := Module[{rows, value, structure},
      rows = GroupIdenticalAlgebraTerms[expression, w];
      Require[ListQ[rows], "Endpoint-structure grouping failed", rows];
      value = Together[Total[Map[Function[row,
        structure = row["Structure"];
        If[!KeyExistsQ[endpointSeries, structure],
          endpointSeries[structure] = Normal[Series[structure /. w -> 1 - wb, {wb, 0, 0}]]];
        row["Coefficient"] endpointSeries[structure]], rows]]];
      (* Keep negative powers of wb until all structures have been added:
         removable endpoint poles need not cancel within a single term. *)
      Require[FreeQ[value, wb | _Series | _SeriesData | _Failure | $Failed | $Aborted |
        ComplexInfinity | Indeterminate | DirectedInfinity],
        "Endpoint coefficient is singular at w = 1", expression];
      value];

    branch[B_, lambda_] := Module[{b = laurent[B, 1], ends},
      ends = endpoint /@ b;
      <|"Pole2" -> -ends[[1]]/lambda,
        "Pole1" -> {(b[[1]] - ends[[1]])/(1 - w), ends[[1]], -ends[[2]]/lambda},
        "Finite" -> {(b[[2]] - ends[[2]])/(1 - w) - lambda (b[[1]] - ends[[1]]) Log[1 - w]/(1 - w),
          ends[[2]], -lambda ends[[1]], -ends[[3]]/lambda}|>];

    parts = KeyValueMap[branch[#2, #1] &, branches];
    regular = laurent[A, 0];
    <|"Pole2" -> Total[#["Pole2"] & /@ parts],
      "Pole1" -> {regular[[1]], 0, 0} + Total[#["Pole1"] & /@ parts],
      "Finite" -> {regular[[2]], 0, 0, 0} + Total[#["Finite"] & /@ parts],
      "Basis" -> ThresholdBasis|>];

EndpointExpand[A_, B_, lambda_, eps_Symbol, w_Symbol, assumptions_: True] :=
  EndpointExpand[A, <|lambda -> B|>, eps, w, assumptions];

(* write a quantity that is positive in the physical region as a product of
   factors that are each positive there, so that the powers can be split
   without producing complex logarithms; "point" is an interior point *)
positiveFactors[base_, point_] := Module[{factors, sign = 1, result = 1},
  factors = FactorList[base];
  Do[If[N[f[[1]] /. point] > 0, result *= f[[1]]^f[[2]],
      result *= (Expand[-f[[1]]])^f[[2]]; sign *= (-1)^f[[2]]], {f, factors}];
  Require[sign == 1, "A power base is not positive in the physical region", base];
  result];

ThresholdSplit[expr_, w_Symbol, eps_Symbol, point_: {}] := Module[{e, wb, tp, interior, lambdas, regular, branches},
  interior = Join[point, {wb -> 1/2, Symbol["s"] -> 3, Symbol["v"] -> 1/3, Symbol["muR"] -> 1}];
  e = expr /. w -> 1 - wb;
  e = e /. Power[base_, n_] /; !FreeQ[n, eps] :> Power[positiveFactors[base, interior], n];
  (* (f1 f2 ...)^n -> f1^n f2^n ..., applied only to these powers *)
  e = e //. {Power[b_Times, n_] /; !FreeQ[n, eps] :> Times @@ (Power[#, n] & /@ List @@ b),
    Power[Power[f_, m_], n_] /; !FreeQ[n, eps] :> Power[f, Expand[m n]]};
  (* tag the threshold exponent: wb^(a - lambda eps) -> wb^a tp^lambda *)
  e = e /. Power[wb, n_] /; !FreeQ[n, eps] :> wb^(n /. eps -> 0) tp^(-Coefficient[Expand[n], eps]);
  (* the terms are kept separate (no common denominator): combining hundreds of
     terms into one rational function makes the later series expansions very slow *)
  e = Collect[Expand[e], tp];
  lambdas = Union[Cases[e, tp^n_ :> n, Infinity], Cases[e, tp :> 1, {0, Infinity}]];
  Require[AllTrue[lambdas, IntegerQ[#] && # > 0 &], "Unexpected threshold exponent in ThresholdSplit", lambdas];
  regular = Coefficient[e, tp, 0];
  branches = Association @ Table[lambda -> wb Coefficient[e, tp, lambda], {lambda, lambdas}];
  Require[FreeQ[Normal[Series[regular, {wb, 0, 0}]], wb], "The part without threshold branch is singular at w = 1"];
  <|"Regular" -> (regular /. wb -> 1 - w), "Branches" -> (branches /. wb -> 1 - w)|>];

(* f(1): the plain substitution w -> 1 when it is finite; otherwise (a removable singularity such as
   log(w)/(w-1), for which the substitution gives 0/0) the limit w -> 1, taken by atOne as the leading term of
   the series in 1-w, which aborts if f is genuinely singular at w = 1 *)
valueAtOne[f_, w_] := Module[{direct = Quiet[Check[f /. w -> 1, $Failed]]},
  If[direct =!= $Failed && FreeQ[direct, ComplexInfinity | Indeterminate | DirectedInfinity], direct, atOne[f, w, True]]];

CanonicalDistribution[prefactor_, {R_, B_, C_, A_}, w_Symbol] := Module[{fB, fC, fA, fB1, fC1, fA1},
  fB = prefactor B; fC = prefactor C; fA = prefactor A;
  {fB1, fC1, fA1} = valueAtOne[#, w] & /@ {fB, fC, fA};
  {prefactor R + (fB - fB1)/(1 - w) + (fC - fC1) Log[1 - w]/(1 - w), fB1, fC1, fA1}];

End[];
EndPackage[];
