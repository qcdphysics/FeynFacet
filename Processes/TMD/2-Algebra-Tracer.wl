(* ::Title:: *)
(*TMD beam and soft functions, step 2 (Tracer route): squared matrix elements*)


(* ::Text:: *)
(*Same task and same output format as the FeynCalc route (read its introduction for the definitions), now with the Tracer package.  Input: the amplitude records of step 1.  Output: Results/Tracer/Algebra.wl.  Run this file in a fresh kernel: Tracer and FeynCalc cannot be loaded together.*)
(**)
(*How an operator matrix element becomes a Tracer trace.  A record of a quark-operator amplitude has one open chain Chain[Ubar[q], W..., Uspinor[p]] (or Vspinor[k] for the g -> q qbar channel), where Ubar[q] is the placeholder for the operator's open Dirac index.  Squaring, summing over spins and closing the operator index with gamma^+ = n-slash gives*)
(**)
(*   Tr[ n-slash W_i rho Wbar_j ],   rho = p-slash (incoming quark) or k-slash (outgoing antiquark),*)
(**)
(*in Tracer G[1, n, W_i..., p, Wbar_j...] once Spur[1] is set.  For the gluon operator there is no chain: the operator leg is the polarization vector PolVector[q, mu, -1], and the two operator indices are contracted with -g^{mu nu}.  External gluons are summed with the physical polarization sum of reference vector n.  Tracer vocabulary: S[p,q] = p.q, S[p,{mu}] = p^mu, S[{mu},{nu}] = g^{mu nu}; VectorDimension[D] sets the dimension; OnShell installs the kinematics.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadTracer[]


(* ::Input:: *)
records = LoadAmplitudes[]["Records"];
VectorDimension[D];
Spur[1];
ShowRecords[Flatten[records["BeamQQ"]]]


(* ::Section:: *)
(*From records to traces*)


(* ::Text:: *)
(*conjugate[record]: i -> -i, reversed gamma strings, transposed generators, conjugated polarization vectors, fresh names for the Lorentz indices; f^{abc}, the deltas and the eikonal denominators are real.  sewTrace builds the Tracer string of the one fermion line: the operator spinor Ubar[q] becomes n, an ordinary outgoing spinor (the quark of q -> g) its momentum, the incoming spinor its momentum.*)


(* ::Input:: *)
conjugate[record_] := Module[{indices, fresh},
  indices = Union[Cases[record, i_Symbol /; StringStartsQ[SymbolName[i], "Lor"] :> i, Infinity]];   (* every Lorentz index Lor1, Lor2, ... *)
  fresh = Thread[indices -> (Symbol["c" <> SymbolName[#]] & /@ indices)];
  <|"Coefficient" -> (record["Coefficient"] /. fresh /. {Complex[a_, b_] :> Complex[a, -b],
      ColorT[l_, i_, j_] :> ColorT[Reverse[l], j, i], PolVector[l_, i_, sign_] :> PolVector[l, i, -sign]}),
    "Chains" -> (Reverse /@ record["Chains"] /. fresh)|>];

incoming[chain_Chain] := First[Cases[chain, (Uspinor | Vspinor)[l_] :> l]];
outgoing[chain_Chain] := First[Cases[chain, Ubar[l_] :> l]] /. q -> n;      (* the operator leg ubar(q) closes with n-slash *)
sewTrace[chainL_, chainR_] := Inactive[G] @@ Join[{1, outgoing[chainL]}, TracerString[chainL], {incoming[chainL]}, TracerString[chainR]];   (* chainR is already reversed *)

conjugate[records["BeamQQ"][[2, 1]]]


(* ::Text:: *)
(*Colour: Tr(T^a T^a) = CF Nc, delta_ii = Nc, delta^aa = Nc^2-1, f^{abc} f^{abc} = CA (Nc^2-1).  Polarization sums: physical sum with reference vector n for the gluons k and p, -g^{mu nu} for the operator leg q and for the soft gluon; they are written with the record vocabulary (Metric, MomentumComponent).*)
(**)
(*The metric tensors (gluon propagators, polarization sums) are contracted before anything is given to Tracer: in every term, g^{ij} renames the index j to i wherever it occurs, including inside the gamma strings, and g^{ii} = D.  (Tracer's automatic contraction is not reliable for chains of several metric tensors.)  Afterwards scalars translates the remaining objects into Tracer's S: Den[q] = 1/q^2, Eikonal[q, v] = 1/(v.q), l^mu = S[l, {mu}], a.b = S[a, b], all expanded in sums of momenta so that S always sees single vectors.*)


(* ::Input:: *)
colourSum = {ColorT[{a_}, i_, j_] ColorT[{a_}, j_, i_] :> CF Nc, ColorDelta[i_, j_]^2 :> Nc, GluonDelta[a_, b_]^2 :> Nc^2 - 1,
  ColorF[a_, b_, c_]^2 :> CA (Nc^2 - 1)};
physicalSum = {PolVector[l_, i_, s_] PolVector[l_, j_, t_] /; l =!= q && s =!= t :>
    -Metric[i, j] + (MomentumComponent[l, i] MomentumComponent[n, j] + MomentumComponent[n, i] MomentumComponent[l, j])/MomentumProduct[l, n],
  PolVector[q, i_, s_] PolVector[q, j_, t_] /; s =!= t :> -Metric[i, j]};
feynmanSum = PolVector[l_, i_, s_] PolVector[l_, j_, t_] /; s =!= t :> -Metric[i, j];

contractTerm[term_] := Module[{t = term /. Metric[i_, i_] -> D /. Metric[i_, j_]^2 -> D, m},
  While[(m = FirstCase[t, Metric[i_, j_] :> {i, j}, None]) =!= None,
    t = (t /. Metric[m[[1]], m[[2]]] -> 1) /. m[[2]] -> m[[1]] /. Metric[i_, i_] -> D /. Metric[i_, j_]^2 -> D];
  t];
contractMetrics[expr_] := With[{e = Expand[expr]}, If[Head[e] === Plus, contractTerm /@ e, contractTerm[e]]];

dot[a_Plus, b_] := dot[#, b] & /@ a;  dot[a_, b_Plus] := dot[a, #] & /@ b;
dot[c_?NumericQ a_, b_] := c dot[a, b];  dot[a_, c_?NumericQ b_] := c dot[a, b];  dot[a_Symbol, b_Symbol] := S[a, b];
component[l_Plus, i_] := component[#, i] & /@ l;  component[c_?NumericQ l_, i_] := c component[l, i];  component[l_Symbol, i_] := S[l, {i}];
scalars = {Den[l_] :> 1/dot[l, l], Eikonal[l_, v_] :> 1/dot[v, l], Metric[i_, j_] :> S[{i}, {j}],
  MomentumComponent[l_, i_] :> component[l, i], MomentumProduct[a_, b_] :> dot[a, b]};

interference[termL_, termR_, polarizationSum_] := Module[{conj = conjugate[termR], coefficient, traces},
  coefficient = termL["Coefficient"] conj["Coefficient"] //. colourSum //. polarizationSum;   (* //.: several polarization pairs per term *)
  traces = Times @@ MapThread[sewTrace, {termL["Chains"], conj["Chains"]}];
  contractMetrics[coefficient traces] /. scalars];
square[amplitudes_, polarizationSum_] := Module[{terms = Flatten[amplitudes]},
  Sum[interference[terms[[i]], terms[[j]], polarizationSum], {i, Length[terms]}, {j, Length[terms]}]];


(* ::Section:: *)
(*Kinematics, Born normalization, real emission*)


(* ::Text:: *)
(*The kinematics is installed inside Internal`InheritedBlock, which keeps it local while retaining Tracer's own definitions of S.  Identical inactive traces are grouped before Activate, so their scalar coefficients are combined and each distinct trace is evaluated once. Products of traces and polarization sums are contracted when the result is expanded.*)


(* ::Input:: *)
beamKinematics := OnShell[on, {p, 0}, {k, 0}, {n, 0}, {nb, 0}, {n, nb, 2}, {p, n, pn}, {k, n, (1 - x) pn}, {p, k, kT2/(2 (1 - x))}, {k, nb, kT2/((1 - x) pn)}, {p, nb, 0}];
evaluate[expr_] := Internal`InheritedBlock[{S}, Module[{terms},
  beamKinematics;
  terms = GroupIdenticalAlgebraTerms[expr, HoldPattern[Inactive[G][__]]];
  Require[ListQ[terms], "The TMD expression could not be separated into scalar coefficients and traces", terms];
  Factor[Expand[Total[(#["Coefficient"] Activate[#["Structure"]]) & /@ terms]]]]];
colourFactors[expr_] := Factor[expr /. CF -> (Nc^2 - 1)/(2 Nc) /. CA -> Nc] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;

born = <|"q" -> colourFactors[evaluate[square[records["BornQ"], physicalSum]] (1/2) (1/Nc)/(2 pn)/gOp^2 /. x -> 1],
  "g" -> colourFactors[evaluate[square[records["BornG"], physicalSum]]/(D - 2)/(Nc^2 - 1)/pn^2/gOp^2 /. x -> 1]|>;
Require[born === <|"q" -> 1, "g" -> 1|>, "The Born matrix elements are not normalized to 1", born];
born


(* ::Input:: *)
real = colourFactors /@ <|
  "qq" -> evaluate[square[records["BeamQQ"], physicalSum]] (1/2) (1/Nc)/(2 pn)/gOp^2,
  "qg" -> evaluate[square[records["BeamQG"], physicalSum]]/(D - 2)/(Nc^2 - 1)/(2 pn)/gOp^2,
  "gq" -> evaluate[square[records["BeamGQ"], physicalSum]] (1/2) (1/Nc)/(x pn^2)/gOp^2,
  "gg" -> evaluate[square[records["BeamGG"], physicalSum]]/(D - 2)/(Nc^2 - 1)/(x pn^2)/gOp^2|>;
Require[FreeQ[real, pn], "The real matrix elements depend on n.p", real];
kernels = Map[Factor[# kT2/gs^2] &, real];
expected = <|"qq" -> 2 CF ((1 + x^2)/(1 - x) - eps (1 - x)), "qg" -> 2 (1/2) (1 - 2 x (1 - x)/(1 - eps)),
  "gq" -> 2 CF ((1 + (1 - x)^2)/x - eps x), "gg" -> 4 CA (x/(1 - x) + (1 - x)/x + x (1 - x))|>;
Require[And @@ Table[Simplify[kernels[c]/(1 - x) /. D -> 4 - 2 eps] === Simplify[expected[c]], {c, Keys[kernels]}],
  "The real-emission kernels differ from the D-dimensional splitting functions", kernels];
TableForm[Normal[kernels], TableHeadings -> {None, {"channel", "kT2 M2/gs^2"}}]


(* ::Text:: *)
(*Gauge-invariance check: the qq channel with -g^{mu nu} for the emitted gluon.*)


(* ::Input:: *)
Require[Simplify[colourFactors[evaluate[square[records["BeamQQ"], feynmanSum]] (1/2) (1/Nc)/(2 pn)/gOp^2] - real["qq"]] === 0,
  "The qq matrix element depends on the gluon polarization sum"];


(* ::Section:: *)
(*The soft matrix element*)


(* ::Input:: *)
soft = colourFactors /@ Internal`InheritedBlock[{S},
  OnShell[on, {k, 0}, {n, 0}, {nb, 0}, {n, nb, 2}, {k, n, kn}, {k, nb, kb}];
  <|"q" -> Factor[Expand[Activate[square[records["SoftQ"], feynmanSum]]]]/Nc/gSoft^2,
    "g" -> Factor[Expand[Activate[square[records["SoftG"], feynmanSum]]]]/(Nc^2 - 1)/gSoft^2|>];
Require[Simplify[soft["q"] - 4 gs^2 CF/(kn kb)] === 0 && Simplify[soft["g"] - 4 gs^2 CA/(kn kb)] === 0, "Soft matrix elements differ from 4 gs^2 C_i/((n.k)(nb.k))", soft];
soft


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Soft" -> soft|>;
SaveResult["Algebra", algebra, "Tracer"];
Keys[algebra]
