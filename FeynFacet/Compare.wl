(* ::Package:: *)

(* Compare.wl: exact and numerical comparison of a computed coefficient with
   a reference expression, on the physical branch of every logarithm. *)

BeginPackage["FeynFacet`"];

ExactDifference::usage = "ExactDifference[a, b, assumptions] returns the \
simplified difference a - b after expanding logarithms and powers on the \
physical branch defined by assumptions (e.g. s > 0 && 0 < v < 1 && 0 < w < 1).  \
Zero means exact agreement.  Lists are compared element by element.";

NumericDifference::usage = "NumericDifference[a, b, points, digits] evaluates a - b \
at each replacement list in points to the given number of digits; results below \
10^(-digits/2) are set to zero.";

PhysicalPoints::usage = "PhysicalPoints is a list of rational test points in the \
physical region, with distinct scales, used by the numerical checks.";

Begin["`Private`"];

ExactDifference[a_, b_, assumptions_: True] := Module[{d},
  d = Together[a - b];
  d = FunctionExpand[d];
  d = PowerExpand[d, Assumptions -> assumptions];
  d = FullSimplify[d, assumptions];
  Factor[Together[d]]];

NumericDifference[a_, b_, points_List, digits_: 30] :=
  Table[Quiet[Chop[N[(a - b) /. point, digits], 10^(-digits/2)], N::meprec], {point, points}];

(* the symbols are resolved when PhysicalPoints is used, so they are the notebook's own s, v, w, ... *)
PhysicalPoints := With[{s = Symbol["s"], v = Symbol["v"], w = Symbol["w"], nf = Symbol["Nf"],
    muR = Symbol["muR"], muF = Symbol["muF"], muFF = Symbol["muFF"]}, {
  {s -> 17, v -> 2/5, w -> 3/7, nf -> 4, muR -> Sqrt[29], muF -> Sqrt[19], muFF -> Sqrt[23]},
  {s -> 31, v -> 5/11, w -> 4/9, nf -> 5, muR -> Sqrt[43], muF -> Sqrt[37], muFF -> Sqrt[41]},
  {s -> 7, v -> 3/4, w -> 2/3, nf -> 3, muR -> Sqrt[13], muF -> Sqrt[5], muFF -> Sqrt[11]}}];

End[];
EndPackage[];
