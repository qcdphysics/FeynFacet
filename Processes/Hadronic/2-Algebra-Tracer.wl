(* ::Title:: *)
(*Proton-proton classes, step 2 (Tracer route): spin and colour algebra*)


(* ::Text:: *)
(*Same task, same input and same output as 2-Algebra-FeynCalc.wl; read its introduction first, which defines the three parts of the output (Born, Virtual, Real), the invariants s, t, u, A, B, C, spe3 and the scalar-integral families of the one-loop interference.  Here the Dirac and Lorentz algebra is done by the classic Tracer package (M. Jamin and M. E. Lautenbacher, Comput. Phys. Commun. 74 (1993) 265), and the colour algebra by an explicit SU(N) reduction written out in this file.  Input: the amplitude records of step 1 (Results/class-<n>/Amplitudes.wl, key "Records").  Output: Results/Tracer/class-<n>/Algebra.wl, which the last cell compares exactly with Results/FeynCalc/class-<n>/Algebra.wl.  Run this file in a fresh kernel: Tracer and FeynCalc cannot be loaded together.*)
(**)
(*Spin sums (HVBM scheme for LL), as in the FeynCalc route.  The physical momenta pa, pb, p are four-dimensional; the unobserved k1, k2 and the loop momentum have (D-4)-dimensional ("hatted") parts.  Incoming quark: p-slash/2 (UU) or +-p-slash gamma5/2 (LL; plus for a quark, minus for an antiquark, since u ubar = p-slash (1 + h gamma5)/2 and v vbar = p-slash (1 - h gamma5)/2 for a massless fermion of helicity h); the factor 1/2 is part of the average below.  Incoming gluon: the physical polarization sum averaged with 1/(D-2) (UU), or the circular density (LL), which already contains the factor 1/2 of the helicity difference.  Outgoing gluons: the physical sum -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r) with a reference vector r.*)
(**)
(*Tracer's vocabulary (summarized in FeynFacet/Tracer.wl and introduced in Processes/LeptonQuark/2-Algebra-Tracer.wl):*)
(*  G[line, p, {mu}, ...]   the string p-slash gamma^mu ... on a fermion line; it is a trace once Spur[line] has been declared; G5 inside it is gamma5;*)
(*  S[p, q] = p.q,  S[p, {mu}] = p^mu,  S[{mu}, {nu}] = g^{mu nu}  (D-dimensional); an index that occurs twice in an expanded product is contracted automatically, also when the other occurrence is inside any expression f[..., {mu}, ...]: S[{mu}, p] f[..., {mu}, ...] = f[..., p, ...];*)
(*  H[p, q]   the hatted, (D-4)-dimensional part of p.q;*)
(*  Eps[a, b, c, d]   the four-dimensional epsilon tensor (Tracer stores it internally as eps[{a, b, c, d}]);*)
(*  VectorDimension[D]   D-dimensional vectors and gamma matrices;  AntiCommute[off]   the HVBM gamma5, which anticommutes with the four-dimensional gamma matrices and commutes with the hatted ones.*)
(**)
(*Sign of the epsilon tensor.  Tracer defines gamma5 = (i/4!) eps_{mu nu rho sigma} gamma^mu gamma^nu gamma^rho gamma^sigma, so that Tr[gamma5 a-slash b-slash c-slash d-slash] = +4 i Eps[a, b, c, d]: the opposite sign to FeynCalc with $LeviCivitaSign = -1.  The circular density of a polarized incoming gluon of momentum l, i eps^{mu nu l r}/(2 l.r) in the FeynCalc route, therefore reads -i Eps[mu, nu, l, r]/(2 l.r) here.  This sign is tested by the Born squares, which contain a gamma5 and a circular density together in the q g channels and two circular densities in the g g channels.*)


(* ::Section:: *)
(*Load*)


(* ::Text:: *)
(*The FeynFacet package (file locations, the record format, the checks Require[...]) and then Tracer.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*Tracer itself.  It redefines Dot and NonCommutativeMultiply on startup; LoadTracer restores them, so that ordinary Mathematica expressions keep working.*)


(* ::Input:: *)
LoadTracer[]


(* ::Text:: *)
(*Choices: the class; the number of parallel kernels for the squares ("auto" = one per processor core, at most 8; "1" = none); and whether the pair files of an earlier run of this same, unchanged file may be reused ("no" recomputes everything; "yes" only continues an interrupted run: a product file is reused only if its companion file <name>.input holds the fingerprint of the current inputs, the amplitudes of step 1 and this step file; a file written for other inputs, or without that companion (for example, after an interrupted write), is computed again).*)
(**)
(*The records were written in step 1 by a kernel with FeynCalc loaded, in which the number of light flavours is FeynCalc's symbol Nf (context FeynCalc`); it is renamed to the plain symbol Nf of the FeynCalc route's output.  Tracer is set to D dimensions with the HVBM gamma5, and the lines 1, 2, 3 are declared traces: lines 1 and 2 for the (at most two) traces formed by closing the open quark lines of an interference, line 3 for a closed quark loop of a one-loop graph.*)


(* ::Input:: *)
class = "13";
workers = "auto";
resume = "no";
loaded = LoadAmplitudes["class-" <> class];
plainNf = Thread[Select[Union[Cases[loaded, _Symbol, Infinity, Heads -> True]], SymbolName[#] === "Nf" && Context[#] =!= "Global`" &] -> Nf];
records = loaded["Records"] /. plainNf; spec = loaded["Spec"] /. plainNf; weights = loaded["Weights"] /. plainNf;
VectorDimension[D];
AntiCommute[off];
Spur[1, 2, 3];
spec


(* ::Text:: *)
(*The first record of the first diagram of each part, as an example of the input: a coefficient (couplings, colour factors ColorT, ColorF, ColorDelta, GluonDelta, propagators Den[q] = 1/q^2, metric tensors, momentum components, polarization vectors PolVector[l, mu, +1 incoming / -1 outgoing]) and the open fermion chains Chain[Ubar or Vbar, gamma matrices, Uspinor or Vspinor].  A one-loop record with a closed quark loop also has "Traces": the loop's Dirac trace, a string of gamma matrices, which is traced.  Legacy records also list there the other factors of the loop (numbers, couplings, generators) as separate entries; these multiply.  Everything FeynArts attaches to the closed loop (its sign, the factor Nf of the flavours in it) is already part of the record.*)


(* ::Input:: *)
ShowRecords[Join[If[records["Born"] === {}, {}, {records["Born"][[1, 1]]}], {First[First[records["Real"]]][[1]]}]]


(* ::Text:: *)
(*The QCD amplitudes contain only gamma^mu and slashed momenta, no gamma5 (the gamma5 of the LL spin sums is inserted below).  This matters for the complex conjugation below, which reverses a gamma string; that is the Dirac adjoint only for gamma^mu and slashed momenta (gamma0 gamma5^dagger gamma0 = -gamma5 would need a sign).*)


(* ::Input:: *)
Require[FreeQ[records, GammaFive], "A record contains gamma5; the conjugation below does not handle it"];
"no gamma5 in the amplitudes"


(* ::Section:: *)
(*Conventions*)


(* ::Text:: *)
(*Incoming partons and averages, as in the FeynCalc route: a quark is averaged with 1/(2 Nc) (spin and colour), a gluon with 1/((D-2)(Nc^2-1)) for UU and 1/(Nc^2-1) for LL (the 1/2 of the helicity difference is in the circular density); polarizedQuarkSign is +1 for a quark and -1 for an antiquark.  colourFactors writes the result with CA = Nc and CF = (Nc^2-1)/(2 Nc).*)
(**)
(*Kinematics are installed by functions, because every parallel kernel installs them for each job; they are always called inside Internal`InheritedBlock[{S, H}, ...], which keeps the declared scalar products local to one job while retaining Tracer's own rules for S and H.  The physical momenta pa, pb, p have no hatted components (H = 0; H is also made linear, since Tracer leaves H of a momentum sum unexpanded).  In the real emission the unobserved k1 and k2 have the hatted products H[k1, k1] = H[k2, k2] = -H[k1, k2] = spe3 (k1hat = -k2hat, since the physical momenta have no hatted parts); in the Born and one-loop kinematics k1 is four-dimensional too, and the loop momentum keeps its hatted part, H[loop, loop] = loopHat2 (see the one-loop section).  Inside the amplitudes k1 = pa + pb - p (Born and one loop) and k2 = pa + pb - p - k1 (real) are eliminated, except in the spinors, whose spin sums give the slashed momentum itself.*)


(* ::Input:: *)
incoming = (spec["Born"] /. None -> spec["Real"][[1, 1]])[[1 ;; 2]];       (* the two incoming partons *)
polarizedQuarkSign = <|pa -> If[StringEndsQ[incoming[[1]], "b"], -1, 1], pb -> If[StringEndsQ[incoming[[2]], "b"], -1, 1]|>;
incomingGluons = Select[{pa, pb}, incoming[[Position[{pa, pb}, #][[1, 1]]]] === "g" &];
incomingQuarks = Complement[{pa, pb}, incomingGluons];
average[spinCase_] := Times @@ (Switch[incoming[[#]], "g", If[spinCase === "LL", 1, 1/(D - 2)]/(Nc^2 - 1), _, 1/(2 Nc)] & /@ {1, 2});
colourFactors[expr_] := Factor[expr] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;

declare[products_List] := (Unprotect[S]; Scan[(S[#[[1]], #[[2]]] = #[[3]]) &, products]; Protect[S]);
fourDimensional[momenta_List] := (Unprotect[H];
  H[x_, y_Plus] := H[x, #] & /@ y; H[x_Plus, y_] := H[#, y] & /@ x;
  H[x_, c_?NumericQ y_] := c H[x, y]; H[c_?NumericQ x_, y_] := c H[x, y];
  Scan[(H[#] = 0; H[#, _] = 0) &, momenta]; Protect[H]);
bornKinematics := (fourDimensional[{pa, pb, p, k1}];
  declare[{{pa, pa, 0}, {pb, pb, 0}, {p, p, 0}, {k1, k1, 0}, {pa, pb, s/2}, {pa, p, -t/2}, {pb, p, -u/2},
    {pa, k1, -u/2}, {pb, k1, -t/2}, {p, k1, s/2}}]);
realKinematics := (fourDimensional[{pa, pb, p}];
  Unprotect[H]; H[k1, k1] = spe3; H[k2, k2] = spe3; H[k1, k2] = -spe3; Protect[H];
  declare[{{pa, pa, 0}, {pb, pb, 0}, {p, p, 0}, {k1, k1, 0}, {k2, k2, 0}, {pa, pb, s/2}, {pa, p, -t/2}, {pb, p, -u/2},
    {pa, k1, -A/2}, {pb, k1, -B/2}, {p, k1, C/2}, {pa, k2, (s + t + A)/2}, {pb, k2, (s + u + B)/2}, {p, k2, (-t - u - C)/2},
    {k1, k2, (s + t + u)/2}}]);
installKinematics[name_] := If[name === "Born", bornKinematics, realKinematics];
bornConservation = k1 -> pa + pb - p;
realConservation = k2 -> pa + pb - p - k1;
conservedMomentum[rule_][q_] := Expand[q /. rule];      (* the map of momenta inside an amplitude: eliminate k1 or k2 *)
epsHead = Head[Eps[a1, a2, a3, a4]];          (* Tracer's internal head of the epsilon tensor *)
{incomingQuarks, incomingGluons, average["UU"]}


(* ::Section:: *)
(*Colour: an explicit SU(N) reduction*)


(* ::Text:: *)
(*The colour factor of an interference is a closed colour network: the colour structure of a term of the amplitude times the complex conjugate of a colour structure of the other amplitude, all indices summed.  It is reduced with three identities, in the normalization Tr(T^a T^b) = delta^{ab}/2 of the FeynCalc route:*)
(**)
(*  (1) structure constants as traces:  f^{abc} = -2 i ( Tr[T^a T^b T^c] - Tr[T^b T^a T^c] )   (from [T^a, T^b] = i f^{abc} T^c);*)
(*  (2) the Fierz identity, which removes a summed gluon index:  (T^a)_{ij} (T^a)_{kl} = 1/2 ( delta_{il} delta_{kj} - delta_{ij} delta_{kl}/Nc );*)
(*  (3) closed loops:  delta_{ii} = Nc,  Tr[T^a] = 0,  and delta^{aa} = Nc^2 - 1 for a gluon delta.*)
(**)
(*Notation: generators[{a, b, ...}, i, j] = (T^a T^b ...)_{ij} is a string of generators between the quark colours i and j (generators[{}, i, j] = delta_{ij}; a closed trace Tr[T^a T^b ...] is a string whose two ends carry the same summed index); adjointDelta[a, b] = delta^{ab}.  A network is first written in these objects: each ColorT is a string, each ColorDelta an empty string, each GluonDelta an adjointDelta, and each ColorF becomes, by (1), the difference of two closed traces with new summed quark indices (fresh names, one pair for every f, also when the same f occurs twice as f^{abc} f^{abc}).  Then the following rules are applied until nothing changes:*)
(**)
(*  adjoint deltas: delta^{aa} = Nc^2 - 1, delta^{ab} delta^{ab} = Nc^2 - 1, delta^{ab} delta^{bc} = delta^{ac}, and delta^{ab} inside a string renames b to a;*)
(*  closed loops (3): generators[{}, i, i] = Nc, generators[{a}, i, i] = 0, delta_{ij} delta_{ij} = Nc;*)
(*  joining: a delta_{ij} next to a string renames the index; two strings that share a quark index j, (X)_{ij} (Y)_{jk} = (X Y)_{ik}, are joined into one;*)
(*  the Fierz identity (2), for a gluon index a that occurs twice in one string or once in each of two strings:*)
(*     generators[{X, a, Y, a, Z}, i, j] = 1/2 ( Tr[Y] generators[{X, Z}, i, j] - generators[{X, Y, Z}, i, j]/Nc ),*)
(*     generators[{X, a, Y}, i, j] generators[{U, a, V}, k, l] = 1/2 ( generators[{X, V}, i, l] generators[{U, Y}, k, j] - generators[{X, Y}, i, j] generators[{U, V}, k, l]/Nc ).*)
(**)
(*Every application of (2) removes one gluon index, so the reduction ends with a polynomial in Nc and 1/Nc; a network that still contains a colour object afterwards is reported as a failure.  The complex conjugate of a colour structure transposes the generators, (T^a T^b)_{ij}^* = (T^b T^a)_{ji} (f and the deltas are real), and renames the indices summed inside it so that they cannot coincide with those of the amplitude.  Each distinct network is evaluated once and remembered (colourCache); networks that differ only by the names of their summed indices are recognized by a canonical renaming: the factors are sorted by their shape with the indices blanked out, and the indices are then numbered colour1, colour2, ... in order of first appearance.*)


(* ::Input:: *)
colourHeads = ColorT | ColorF | ColorDelta | GluonDelta;
factorList[term_] := Flatten[Replace[If[Head[term] === Times, List @@ term, {term}], Power[x_, n_Integer?Positive] :> ConstantArray[x, n], {1}]];
termList[expr_] := If[Head[expr] === Plus, List @@ expr, {expr}];
colourIndices[colour_] := Flatten[factorList[colour] /. {ColorT[l_List, i_, j_] :> {l, i, j}, (ColorF | ColorDelta | GluonDelta)[x__] :> {x}, _?NumericQ -> {}}];
internalColour[colour_] := Select[Tally[colourIndices[colour]], #[[2]] > 1 &][[All, 1]];      (* the indices summed inside one colour structure *)
conjugateColour[colour_] := colour /. Thread[internalColour[colour] -> (Symbol["c" <> SymbolName[#]] & /@ internalColour[colour])] /.
  ColorT[l_List, i_, j_] :> ColorT[Reverse[l], j, i];

(* canonical names: factors sorted by shape (indices blanked), then indices numbered by first appearance *)
canonicalColour[network_] := Module[{factors = factorList[network], sorted, order},
  sorted = SortBy[factors, {ToString[# /. {ColorT[l_, _, _] :> ColorT[Length[l], 0, 0], (h : ColorF | ColorDelta | GluonDelta)[x__] :> h[]}] &, ToString}];
  order = DeleteDuplicates[Flatten[colourIndices /@ sorted]];
  Times @@ sorted /. Thread[order -> Table[Symbol["colour" <> ToString[n]], {n, Length[order]}]]];

(* a network in generator strings; identity (1) for every f, factor by factor, so that each copy of a repeated f gets its own fresh indices *)
freshIndex[] := Unique["x"];
SetAttributes[adjointDelta, Orderless];
toGeneratorStrings[network_] := Times @@ (Replace[#, {ColorT[l_List, i_, j_] :> generators[l, i, j], ColorDelta[i_, j_] :> generators[{}, i, j], GluonDelta[a_, b_] :> adjointDelta[a, b],
  ColorF[a_, b_, c_] :> With[{x = freshIndex[], y = freshIndex[]}, -2 I (generators[{a, b, c}, x, x] - generators[{b, a, c}, y, y])]}] & /@ factorList[network]);
(* the reduction rules, in the order of the text above *)
colourRules = {
  adjointDelta[a_, a_] :> Nc^2 - 1, adjointDelta[a_, b_]^2 :> Nc^2 - 1,
  adjointDelta[a_, b_] adjointDelta[b_, c_] :> adjointDelta[a, c],
  adjointDelta[a_, b_] generators[{x___, b_, y___}, i_, j_] :> generators[{x, a, y}, i, j],
  generators[{}, i_, i_] :> Nc, generators[{_}, i_, i_] :> 0, generators[{}, i_, j_]^2 :> Nc,
  generators[{}, i_, j_] generators[l_, i_, k_] :> generators[l, j, k], generators[{}, i_, j_] generators[l_, k_, j_] :> generators[l, k, i],
  generators[l1_, i_, j_] generators[l2_, j_, k_] :> generators[Join[l1, l2], i, k],
  generators[{x1___, a_, x2___, a_, x3___}, i_, j_] :> With[{z = freshIndex[]}, 1/2 (generators[{x2}, z, z] generators[{x1, x3}, i, j] - generators[{x1, x2, x3}, i, j]/Nc)],
  generators[{x1___, a_, x2___}, i_, j_] generators[{y1___, a_, y2___}, k_, l_] :>
    1/2 (generators[{x1, y2}, i, l] generators[{y1, x2}, k, j] - generators[{x1, x2}, i, j] generators[{y1, y2}, k, l]/Nc)};
colourCache = <||>;
colourValue[network_] := Module[{key = canonicalColour[network], value},
  If[KeyExistsQ[colourCache, key], Return[colourCache[key]]];
  value = Factor[FixedPoint[Expand[# /. colourRules] &, Expand[toGeneratorStrings[key]]]];
  If[!FreeQ[value, generators | adjointDelta | colourHeads], Return[Failure["Colour", <|"Network" -> key, "Value" -> value|>]]];
  colourCache[key] = value];
(* the colour sum of an amplitude structure c1 times the complex conjugate of c2 *)
colourSum[c1_, c2_] := colourValue[c1 conjugateColour[c2]];
"colour reduction defined"


(* ::Text:: *)
(*Checks of the reduction on known colour sums (CF = (Nc^2-1)/(2 Nc)): Tr(T^a T^a) = CF Nc; f^{abc} f^{abc} = Nc (Nc^2-1); (T^a T^b)_{ij} ((T^a T^b)_{ij})^* = CF^2 Nc; (T^a T^b)_{ij} ((T^b T^a)_{ij})^* = -CF/2; f^{abc} (T^c)_{ij} ((T^a T^b)_{ij})^* = f^{abc} Tr[T^c T^b T^a] = -(i/4) f^{abc} f^{abc} = -i Nc (Nc^2-1)/4, which tests the sign of identity (1); the closed loop f^{a1 b1 c1} f^{a1 b2 c2} f^{b1 b2 d} f^{c1 c2 d} = Nc^2 (Nc^2-1)/2; and f^{z1 z2 z3} f^{z1 z2 z3}, which differs from the second network only by the names of its indices and must be found in the cache (the number of cached networks must not grow).  The full test of the colour algebra is the comparison of the Born squares with Processes/PartonicBorn below.*)


(* ::Input:: *)
colourTests = {
  colourSum[ColorT[{a1}, i1, i2], ColorT[{a1}, i1, i2]] - (Nc^2 - 1)/2,
  colourSum[ColorF[a1, a2, a3], ColorF[a1, a2, a3]] - Nc (Nc^2 - 1),
  colourSum[ColorT[{a1, a2}, i1, i2], ColorT[{a1, a2}, i1, i2]] - Nc ((Nc^2 - 1)/(2 Nc))^2,
  colourSum[ColorT[{a1, a2}, i1, i2], ColorT[{a2, a1}, i1, i2]] + (Nc^2 - 1)/(4 Nc),
  colourSum[ColorF[a1, a2, b] ColorT[{b}, i1, i2], ColorT[{a1, a2}, i1, i2]] + I Nc (Nc^2 - 1)/4,
  colourValue[ColorF[a1, b1, c1] ColorF[a1, b2, c2] ColorF[b1, b2, d] ColorF[c1, c2, d]] - Nc^2 (Nc^2 - 1)/2};
cacheSize = Length[colourCache];
colourValue[ColorF[z1, z2, z3] ColorF[z1, z2, z3]];
Require[Simplify[colourTests] === {0, 0, 0, 0, 0, 0} && Length[colourCache] === cacheSize, "The colour reduction fails a known colour sum", colourTests];
colourTests


(* ::Section:: *)
(*From records to Tracer expressions*)


(* ::Text:: *)
(*Each record becomes a Tracer expression: Metric[mu, nu] -> S[{mu}, {nu}], MomentumComponent[q, mu] -> S[q, {mu}], PolVector[l, mu, _] -> S[e_l, {mu}] with a vector e_l for the polarization of the gluon l (named epa, epb, ep, ek1, ek2; in the complex-conjugate amplitude epaStar, ...), Den[q] -> 1/S[q, q] for a tree propagator (evaluated with the installed kinematics) or 1/DjOf[q] for a propagator that carries the loop momentum (see the one-loop section).  An open fermion chain ubar(a) W u(b) becomes the inert object word[a, W..., b], with the spinor momenta at its ends and W the Tracer arguments ({mu} for gamma^mu, the momentum q for q-slash).  A closed quark loop becomes G[3, ...], a trace that Tracer evaluates at once, times its other factors.  The momenta inside the amplitude go through the map mom (momentum conservation, and for a one-loop graph also the shift of the loop momentum).  Because the record's factors are multiplied out, Tracer contracts every index summed inside the amplitude: indices of the metric tensors and momentum components disappear, and those attached to a fermion chain move into its word.*)


(* ::Input:: *)
polarization[l_] := Symbol["e" <> ToString[l]];
polarizationStar[l_] := Symbol["e" <> ToString[l] <> "Star"];
spinorMomentum[(Ubar | Vbar | Uspinor | Vspinor)[q_]] := q;
tracerArgument[GammaMatrix[i_], mom_] := {i};
tracerArgument[Slash[q_], mom_] := mom[q];
propagator[q_] := If[FreeQ[q, loop], 1/S[q, q], 1/DjOf[q]];
closedLoop[pieces_List, mom_] := Module[{gammas = Select[pieces, MatchQ[#, _GammaMatrix | _Slash] &]},
  (Times @@ Select[pieces, !MatchQ[#, _GammaMatrix | _Slash] &]) If[gammas === {}, 1, G @@ Join[{3}, tracerArgument[#, mom] & /@ gammas]]];
convert[record_, mom_] := Module[{coefficient, words, loops},
  coefficient = record["Coefficient"] /. {Metric[i_, j_] :> S[{i}, {j}], MomentumComponent[q_, i_] :> S[mom[q], {i}],
    PolVector[l_, i_, _] :> S[polarization[l], {i}], Den[q_] :> propagator[mom[q]]};
  words = Times @@ (word @@ Join[{spinorMomentum[First[#]]}, tracerArgument[#, mom] & /@ (List @@ #)[[2 ;; -2]], {spinorMomentum[Last[#]]}] & /@ record["Chains"]);
  loops = If[KeyExistsQ[record, "Traces"], Times @@ (closedLoop[List @@ #, mom] & /@ record["Traces"]), 1];
  Expand[coefficient words loops]];
"conversion defined"


(* ::Section:: *)
(*The squaring method*)


(* ::Text:: *)
(*As in the FeynCalc route, a square sum_{i,j} M_i M_j^* is not formed as one product:*)
(**)
(*1. Every amplitude (one diagram) is converted and simplified once, before any product, and the conditions e_l.l = 0 and e_l.r = 0 are imposed for every gluon l, where r is the reference vector of the same gluon's polarization sum (both the physical sum and the circular density are transverse to l and to r).  For the gluon whose momentum was eliminated by momentum conservation, e_x.x = 0 is imposed in its eliminated form: x = pa + pb - p (- k1) gives e_x.q = -(the other terms) for the first momentum q of x that is not its reference vector.*)
(*2. Colour is separated from the Lorentz algebra: each amplitude is written as sum_x c_x L_x with distinct colour structures c_x; the colour factor of c_x c_y^* is computed by colourValue above.  The complex conjugate of a Lorentz part L_y changes i -> -i, e_l -> e_lStar, renames its summed Lorentz indices, and replaces each word ubar(a) W u(b) by its complex conjugate ubar(b) W^R u(a), with W^R the gamma string in reversed order (wordBar[b, W^R..., a]; reversing is the Dirac adjoint gamma0 W^dagger gamma0 of a string of gamma^mu and slashed momenta).*)
(*3. Only the products with i <= j are computed: for i < j, M_i M_j^* + M_j M_i^* = 2 Re(M_i M_j^*) (all invariants are real).*)
(*4. The fermion lines are closed: in the product of the amplitude's words and the conjugate amplitude's words, the word that ends on the spinor of momentum q is followed by the conjugate word that starts there, with the spin sum in between: q-slash for u(q) ubar(q) and v(q) vbar(q) alike, and, for an incoming quark in the LL case, q-slash gamma5 times polarizedQuarkSign (+1 quark, -1 antiquark; gamma5 is placed directly after q-slash, and since q is four-dimensional this is the same as -gamma5 q-slash).  One continues with the next word until the first word comes back: each closed walk is one trace, trace[line, ...], kept inert for now.  The interference has one or two traces, depending on how the quark lines of the two diagrams connect.*)
(*5. The gluon polarization sums are applied one gluon at a time, outgoing gluons first, before the traces are taken: the terms are grouped by the factors that carry that gluon's polarization vectors e_l and e_lStar; in each distinct group e_l -> {mu_l}, e_lStar -> {nu_l}, and the product with the density -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r) (or, for a polarized incoming gluon, the circular density -i Eps[mu, nu, l, r]/(2 l.r)) is expanded once, Tracer contracting the indices, also into the inert traces.  Contracting before tracing keeps the traces short: a contracted pair gamma^mu ... gamma_mu inside a trace is reduced by Tracer's D-dimensional contraction formulas before the trace is expanded.  (Tracer does not move an index of Eps into a trace, so the two indices of a circular density stay free in the trace and are contracted with Eps after it.)*)
(*6. The traces are then evaluated (each distinct one once) and the products expanded; every remaining scalar product is a kinematic invariant.  An epsilon tensor of the eliminated momentum is rewritten with momentum conservation.  The products (i, j) are independent: they run on parallel kernels, and each is written to its own file with the time it took, so that the progress is visible and an interrupted run can be continued.*)


(* ::Input:: *)
(* ---- 1. one amplitude, with the transversality conditions ---- *)
transversality[references_, conservation_] := Module[{x = conservation[[1]], momenta, solved, coefficient},
  momenta = Cases[termList[Expand[conservation[[2]]]], c_. q_Symbol :> q];
  solved = First[DeleteCases[momenta, references[x]]];
  coefficient[q_] := Coefficient[conservation[[2]], q];
  {{S[polarization[x], solved] -> -Sum[coefficient[q] S[polarization[x], q], {q, DeleteCases[momenta, solved]}]/coefficient[solved]},
   Flatten[Table[{S[polarization[l], l] -> 0, S[polarization[l], references[l]] -> 0}, {l, Keys[references]}]]}];
prepare[diagram_, references_, conservation_, mom_] := Module[{a = Total[convert[#, mom] & /@ diagram], rules = transversality[references, conservation]},
  Expand[Expand[a /. rules[[1]]] /. rules[[2]]]];

(* ---- 2. colour separation and complex conjugation ---- *)
colourPart[term_] := Times @@ Select[factorList[term], !FreeQ[#, colourHeads] &];
(* {{c_x, L_x}, ...}: the distinct colour structures of an amplitude with their Lorentz parts *)
colourGroups[amp_] := Module[{groups = GatherBy[termList[Expand[amp]], colourPart]},
  Table[{colourPart[First[g]], Total[(#/colourPart[#]) & /@ g]}, {g, groups}]];
lorentzIndices[expr_] := Union[Cases[expr, {i_Symbol} /; StringStartsQ[SymbolName[i], "Lor"] :> i, Infinity]];    (* FeynArts names them Lor1, Lor2, ... *)
gluonsOf[expr_] := Union[Cases[expr, s_Symbol /; StringMatchQ[SymbolName[s], "e" ~~ ("pa" | "pb" | "p" | "k1" | "k2")] :> ToExpression[StringDrop[SymbolName[s], 1]], Infinity]];
conjugateLorentz[lorentz_] := Module[{gluons = gluonsOf[lorentz], indices = lorentzIndices[lorentz]},
  lorentz /. Thread[polarization /@ gluons -> polarizationStar /@ gluons] /. Thread[indices -> (Symbol["c" <> SymbolName[#]] & /@ indices)] /.
    {Complex[x_, y_] :> Complex[x, -y], word[a_, w___, b_] :> wordBar[b, Sequence @@ Reverse[{w}], a]}];
conjugateGroups[amp_] := {conjugateColour[#[[1]]], conjugateLorentz[#[[2]]]} & /@ colourGroups[amp];

(* ---- 4. closing the fermion lines ---- *)
(* the spin sum at the spinor of momentum q: {Tracer arguments, sign} *)
junction[q_, spinCase_] := If[spinCase === "LL" && MemberQ[incomingQuarks, q], {{q, G5}, polarizedQuarkSign[q]}, {{q}, 1}];
sew[term_, spinCase_] := Module[{factors = factorList[term], words, bars, rest, traces = {}, sign = 1, current, start, chain, q, j},
  words = Association[#[[1]] -> # & /@ Cases[factors, _word]];          (* each word, keyed by the momentum it starts with *)
  bars = Association[#[[1]] -> # & /@ Cases[factors, _wordBar]];        (* each conjugate word, keyed likewise *)
  rest = Times @@ DeleteCases[factors, _word | _wordBar];
  While[Length[words] > 0,
    current = First[words]; start = current[[1]]; chain = {};
    While[True,
      KeyDropFrom[words, current[[1]]];
      q = Last[current]; chain = Join[chain, List @@ current[[2 ;; -2]]];
      j = junction[q, spinCase]; chain = Join[chain, j[[1]]]; sign *= j[[2]];
      Require[KeyExistsQ[bars, q], "A fermion line does not close", term];
      current = bars[q]; q = Last[current]; chain = Join[chain, List @@ current[[2 ;; -2]]];
      j = junction[q, spinCase]; chain = Join[chain, j[[1]]]; sign *= j[[2]];
      If[q === start, Break[]];
      Require[KeyExistsQ[words, q], "A fermion line does not close", term];
      current = words[q]];
    AppendTo[traces, trace @@ Join[{Length[traces] + 1}, chain]]];
  sign rest Times @@ traces];

(* ---- 5. polarization sums ---- *)
physicalSum[l_, r_] := With[{mu = {Symbol["mu" <> ToString[l]]}, nu = {Symbol["nu" <> ToString[l]]}},
  -S[mu, nu] + (S[l, mu] S[r, nu] + S[r, mu] S[l, nu])/S[l, r]];
circular[l_, r_] := With[{mu = {Symbol["mu" <> ToString[l]]}, nu = {Symbol["nu" <> ToString[l]]}},
  -I Eps[mu, nu, l, r]/(2 S[l, r])];                          (* Tracer's sign of Eps, see the introduction *)
tensorPart[term_, l_] := Times @@ Select[factorList[term], !FreeQ[#, polarization[l] | polarizationStar[l]] &];
(* the contraction uses the installed scalar products, which differ between the Born and the real kinematics: the kinematics are part of the key *)
projectionCache = <||>; currentKinematics = None;
projectOnce[tensor_, l_, density_] := With[{key = {currentKinematics, tensor, density}},
  If[KeyExistsQ[projectionCache, key], projectionCache[key],
    With[{value = Expand[(tensor /. {polarization[l] -> {Symbol["mu" <> ToString[l]]}, polarizationStar[l] -> {Symbol["nu" <> ToString[l]]}}) density]},
      If[FreeQ[tensor, trace], projectionCache[key] = value, value]]]];     (* tensors with traces are rarely repeated: not stored *)
project[expr_, l_, density_] := If[expr === 0, 0, Module[{groups = GatherBy[termList[expr], tensorPart[#, l] &]},
  Require[And @@ (!FreeQ[#, polarization[l]] && !FreeQ[#, polarizationStar[l]] & /@ groups[[All, 1]]), "A term lacks a polarization vector", l];
  Expand[Total[Table[With[{tensor = tensorPart[First[g], l]}, Total[(#/tensor) & /@ g] projectOnce[tensor, l, density]], {g, groups}]]]]];

(* ---- 6. traces ---- *)
evaluateTraces[expr_] := Module[{cache = <||>},
  Expand[expr /. trace[n_, a___] :> Lookup[cache, Key[{n, a}], cache[{n, a}] = Expand[G[n, a]]]]];

(* one product (i, j): left = {{c_x, L_x}, ...} of the amplitude, right = {{c_y^*, L_y^*}, ...} of the conjugate amplitude *)
interference[left_, right_, spinCase_, references_, conservation_] := Module[{product, gluons, expr},
  product = Expand[Total[Flatten[Table[colourValue[x[[1]] y[[1]]] x[[2]] y[[2]], {x, left}, {y, right}]]]];
  If[!FreeQ[product, Failure], Return[FirstCase[product, _Failure, Null, Infinity]]];
  If[product === 0, Return[0]];          (* a graph can vanish once e_l.l = e_l.r = 0 are imposed *)
  gluons = gluonsOf[{left[[All, 2]]}];
  expr = Expand[Total[sew[#, spinCase] & /@ termList[product]]];
  Do[expr = project[expr, l, physicalSum[l, references[l]]], {l, Complement[gluons, {pa, pb}]}];       (* outgoing gluons first *)
  Do[expr = project[expr, l, If[spinCase === "LL", circular[l, references[l]], physicalSum[l, references[l]]]], {l, Intersection[gluons, {pa, pb}]}];
  expr = evaluateTraces[expr];
  (* the eliminated momentum inside a four-dimensional epsilon tensor (in the Born kinematics Eps[pa, pb, p, k1] = 0) *)
  Expand[expr /. epsHead[x_List] :> (Eps @@ (x /. conservation))]];
(* the real part: after expanding, every term is a number times real symbols (the denominators are real propagators), so the real part of the number is kept *)
realPartOf[expr_] := Expand[expr] /. Complex[a_, b_] :> a;
imaginaryPartOf[expr_] := (expr - (expr /. Complex[a_, b_] :> Complex[a, -b]))/2;
"squaring functions defined"


(* ::Text:: *)
(*Parallel kernels (choice workers above): each loads FeynFacet and Tracer (Tracer works in parallel kernels).*)


(* ::Input:: *)
feynFacetFile = FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}];
nWorkers = If[workers === "auto", Min[$ProcessorCount, 8], ToExpression[workers]];   (* "auto": one per processor core, at most 8 *)
If[nWorkers > 1,
  LaunchKernels[nWorkers];
  With[{file = feynFacetFile}, ParallelEvaluate[Get[file]; LoadTracer[]]]];
"parallel kernels: " <> ToString[Max[1, Length[Kernels[]]]]


(* ::Text:: *)
(*Each parallel kernel sets the dimension, the HVBM gamma5 and the traced lines, and receives the functions above with empty caches.  A job is one product (i, j) of the squaring method: it installs its kinematics (Born or real) inside Internal`InheritedBlock, computes the interference, and is written to Results/Tracer/class-<n>/pairs/<name>.wl with its time.  square adds the products of one square: a diagonal product as it is, a product with i < j twice its real part (step 3).  A diagonal product |M_i|^2 is real, but its epsilon-tensor terms may leave an imaginary part that vanishes only after simplification: that imaginary part is checked to be exactly zero at a rational point (checkPoint) and dropped.  sidesOf prepares every diagram once in the given kinematics, separates its colour, and forms its conjugate.*)


(* ::Input:: *)
If[nWorkers > 1,
  ParallelEvaluate[VectorDimension[D]; AntiCommute[off]; Spur[1, 2, 3]];
  DistributeDefinitions[declare, fourDimensional, bornKinematics, realKinematics, installKinematics, epsHead,
    colourHeads, factorList, termList, colourIndices, internalColour, conjugateColour, canonicalColour, freshIndex, adjointDelta, toGeneratorStrings, colourRules, colourValue,
    polarization, polarizationStar, junction, sew, physicalSum, circular, tensorPart, projectOnce, project, evaluateTraces, gluonsOf, interference,
    incomingQuarks, incomingGluons, polarizedQuarkSign];
  ParallelEvaluate[colourCache = <||>; projectionCache = <||>; currentKinematics = None]];
pairDirectory = FileNameJoin[{ResultsDirectory["Tracer/class-" <> class], "pairs"}];
If[!DirectoryQ[pairDirectory], CreateDirectory[pairDirectory]];
(* the fingerprint of the inputs of the product files: the amplitudes of step 1 and this step file (see the choice resume above) *)
inputFingerprint = InputFingerprint[loaded, {FileNameJoin[{NotebookDirectory[], "2-Algebra-Tracer.wl"}]}];
runJob[job_] := Module[{file = FileNameJoin[{pairDirectory, job["Name"] <> ".wl"}], value, seconds},
  If[resume === "yes" && ReusableQ[file, inputFingerprint], file,
    {seconds, value} = AbsoluteTiming[Internal`InheritedBlock[{S, H}, installKinematics[job["Kinematics"]]; currentKinematics = job["Kinematics"];
      interference[job["Left"], job["Right"], job["Spin"], job["References"], job["Conservation"]]]];
    If[!FreeQ[value, Failure], value,
      Put[<|"Name" -> job["Name"], "Value" -> value, "Seconds" -> seconds|>, file]; MarkReusable[file, inputFingerprint]]]];
runJobs[jobs_] := Module[{files},
  If[nWorkers > 1, DistributeDefinitions[runJob, pairDirectory, resume, inputFingerprint]];
  files = If[nWorkers > 1, ParallelMap[runJob, jobs, Method -> "FinestGrained"], Map[runJob, jobs]];
  Require[FreeQ[files, Failure], "A job failed", FirstCase[files, _Failure, Null, Infinity]];
  Get /@ files];
(* ---- 3. the pairs i <= j of a square, and their sum ---- *)
squareJobs[name_, sides_, kinematics_, references_, conservation_, spinCase_] := Flatten[Table[
  <|"Name" -> name <> "-" <> spinCase <> "-" <> ToString[i] <> "-" <> ToString[j], "Kinematics" -> kinematics,
    "Left" -> sides[[1, i]], "Right" -> sides[[2, j]], "Spin" -> spinCase, "References" -> references, "Conservation" -> conservation, "Pair" -> {i, j}|>,
  {i, Length[sides[[1]]]}, {j, i, Length[sides[[1]]]}]];
checkPoint = {s -> 1, t -> -3/10, u -> -2/5, A -> -1/7, B -> -2/9, C -> -3/10 + 1/7 + 2/9, D -> 37/10, Nc -> 3, spe3 -> 1/11, gs -> 1, Nf -> 5};
square[name_, sides_, kinematics_, references_, conservation_, spinCase_] := Module[{jobs = squareJobs[name, sides, kinematics, references, conservation, spinCase], results},
  results = runJobs[jobs];
  Require[And @@ MapThread[#2["Pair"][[1]] != #2["Pair"][[2]] || (imaginaryPartOf[#1["Value"]] /. checkPoint) === 0 &, {results, jobs}],
    "A diagonal product has a nonzero imaginary part", name];
  Total[MapThread[If[#2["Pair"][[1]] == #2["Pair"][[2]], realPartOf[#1["Value"]], 2 realPartOf[#1["Value"]]] &, {results, jobs}]]];
sidesOf[diagrams_, kinematics_, references_, conservation_] := Module[{prepared},
  prepared = Internal`InheritedBlock[{S, H}, installKinematics[kinematics];
    prepare[#, references, conservation, conservedMomentum[conservation]] & /@ diagrams];
  {colourGroups /@ prepared, conjugateGroups /@ prepared}];
"jobs defined"


(* ::Section:: *)
(*Born*)


(* ::Text:: *)
(*The Born square with the reference vectors of the FeynCalc route: pa -> pb, pb -> pa (incoming gluons), p -> pa, k1 -> pa; averaged over the incoming spins and colours, and written with CA, CF and u = -s - t.  No identical-particle factor: the tagged parton p is distinguished from the unobserved one.*)


(* ::Input:: *)
bornReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa|>;
bornSides = If[records["Born"] === {}, None, sidesOf[records["Born"], "Born", bornReferences, bornConservation]];
born = If[bornSides === None, <|"UU" -> 0, "LL" -> 0|>,
  Association @ Table[spinCase -> With[{sq = square["born", bornSides, "Born", bornReferences, bornConservation, spinCase]},
    Factor[colourFactors[Together[average[spinCase] sq]] /. u -> -s - t]], {spinCase, {"UU", "LL"}}]];
Require[FreeQ[born, Complex | S | H | epsHead | trace | word | wordBar], "The Born square is not a real scalar", born];
born


(* ::Text:: *)
(*The Born squares must agree with the table of Processes/PartonicBorn (checked against the textbook); this tests the colour reduction, the closing of the fermion lines, the gluon sums and the signs of the polarized densities.*)


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
(*The families, the routing of the loop momentum and the reduction of the loop momentum's scalar products are those of the FeynCalc route:*)
(**)
(*   VirtualST[n1,n2,n3,n4]: 1/[(loop^2)^n1 ((loop+pa)^2)^n2 ((loop+pa+pb)^2)^n3 ((loop+p)^2)^n4]   (external legs pa, pb, -pd, -p in cyclic order: invariants s and t),*)
(*   VirtualSU: ... ((loop+pa+pb)^2)^n3 ((loop+pa+pb-p)^2)^n4   (legs pa, pb, -p, -pd: s and u),      VirtualTU: ((loop+pa)^2)^n2 ((loop+pa-p)^2)^n3 ((loop+pa-p+pb)^2)^n4   (legs pa, -p, pb, -pd: t and u),*)
(**)
(*with negative indices for numerators.  FeynArts routes the loop momentum differently in every graph; for every graph the routing loop -> sigma loop + delta (sigma = +-1) that maps its loop propagators onto one family is found by trying the shifts, the same search in the same order as in the FeynCalc route.*)


(* ::Input:: *)
families = <|
  "VirtualST" -> {loop, loop + pa, loop + pa + pb, loop + p},                 (* legs pa, pb, -pd, -p in cyclic order: invariants s and t *)
  "VirtualSU" -> {loop, loop + pa, loop + pa + pb, loop + pa + pb - p},       (* legs pa, pb, -p, -pd: s and u *)
  "VirtualTU" -> {loop, loop + pa, loop + pa - p, loop + pa - p + pb}|>;      (* legs pa, -p, pb, -pd: t and u *)
familyHeads = <|"VirtualST" -> KVirtualST, "VirtualSU" -> KVirtualSU, "VirtualTU" -> KVirtualTU|>;
shifts = Flatten[Table[c1 pa + c2 pb + c3 p, {c1, -1, 1}, {c2, -1, 1}, {c3, -1, 1}]];
propagatorMomenta[diagram_] := Union[Cases[diagram, Den[q_] /; !FreeQ[q, loop] :> Expand[q /. bornConservation], Infinity]];   (* the loop propagators *)
sameUpToSign[q1_, q2_] := Expand[q1 - q2] === 0 || Expand[q1 + q2] === 0;
(* find sigma, delta and the family such that every loop propagator of the graph, with loop -> sigma loop + delta, is one of the family's *)
routing[diagram_] := Module[{props = propagatorMomenta[diagram], found},
  found = Catch[Do[If[And @@ Table[Or @@ Table[sameUpToSign[q /. loop -> sigma loop + delta, f], {f, families[name]}], {q, props}],
      Throw[{name, sigma, delta}]], {name, Keys[families]}, {sigma, {1, -1}}, {delta, shifts}]; None];
  Require[found =!= None, "No family fits the propagators of a one-loop graph", props]; found];
routings = routing /@ records["Virtual"];
Counts[routings[[All, 1]]]


(* ::Text:: *)
(*For each graph: convert it with the map mom = (shift of the loop momentum, then k1 -> pa + pb - p), so that each loop propagator becomes 1/DjOf[q] with q one of the family's momenta up to sign (a tree propagator next to the loop is evaluated); prepare it once and separate its colour; its products with the conjugate of the summed Born amplitude are the jobs below (one per graph and spin case).*)
(**)
(*After the spin sums and the traces the loop momentum appears only through S[loop, loop], S[loop, pa], S[loop, pb], S[loop, p] (S[loop, k1], from the spin sum of an outgoing quark k1, is rewritten with k1 = pa + pb - p) and the hatted H[loop, loop] = loopHat2.  Tracer's S is D-dimensional and H its (D-4)-dimensional part, so no four-dimensional product needs to be rewritten (in the FeynCalc route loop4.loop4 = loop.loop - loopHat2 is substituted).  The four D-dimensional products are expressed through the inverse propagators Dj of the family (the four equations loop^2 = D1, 2 loop.q_j = D_j - D_1 - q_j^2 are solved for loop^2, loop.pa, loop.pb, loop.p), and each term is read off as coefficient x family[indices] x HatPower[n].  loopHat2 shifts the dimension of the integral in step 3, not the dimension D of the coefficients.*)


(* ::Input:: *)
loopProducts = {L2 -> S[loop, loop], Lpa -> S[loop, pa], Lpb -> S[loop, pb], Lp -> S[loop, p]};
inverseProps[name_] := Module[{q = families[name], sol},
  sol = Internal`InheritedBlock[{S, H}, bornKinematics;
    Solve[Table[Dj[j] == Expand[S[q[[j]], q[[j]]]] /. Reverse /@ loopProducts, {j, 4}] /. u -> -s - t, {L2, Lpa, Lpb, Lp}]];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators"]; First[sol]];
loopRules[name_] := loopRules[name] = With[{sol = inverseProps[name]}, Table[(x /. loopProducts) -> (x /. sol), {x, {L2, Lpa, Lpb, Lp}}]];
shiftedMomentum[{name_, sigma_, delta_}][q_] := Expand[q /. loop -> sigma loop + delta /. bornConservation];
(* the contracted interference of one graph -> sum of coefficient x family[indices] x HatPower[n] *)
toFamilies[expr0_, name_] := Module[{expr},
  expr = Expand[expr0 /. S[loop, k1] :> S[loop, pa + pb - p] /. H[loop, loop] -> loopHat2 /. loopRules[name]];
  Require[FreeQ[expr /. _DjOf -> 1, loop], "A product of the loop momentum was not expressed through the propagators",
    Short[Union[Cases[expr, h_[a___] /; !FreeQ[{a}, loop] && FreeQ[h, Times | Plus | Power | DjOf], Infinity]], 8]];
  expr = Expand[expr /. DjOf[q_] :> Dj @@ (Position[families[name], f_ /; sameUpToSign[f, q]][[1]])];
  (* every term is a coefficient times powers of the Dj (negative for propagators, positive for numerators) and of loopHat2 *)
  Total[Table[With[{powers = Table[Exponent[term, Dj[j]], {j, 4}], hat = Exponent[term, loopHat2]},
    (term/(Times @@ Table[Dj[j]^powers[[j]], {j, 4}] loopHat2^hat)) familyHeads[name] @@ (-powers) HatPower[hat]], {term, DeleteCases[termList[expr], 0]}]]];
virtualLeft = Internal`InheritedBlock[{S, H}, bornKinematics;
  MapThread[colourGroups[prepare[#1, bornReferences, bornConservation, shiftedMomentum[#2]]] &, {records["Virtual"], routings}]];
bornSumRight = If[records["Born"] === {}, {}, Internal`InheritedBlock[{S, H}, bornKinematics;
  conjugateGroups[Total[prepare[#, bornReferences, bornConservation, conservedMomentum[bornConservation]] & /@ records["Born"]]]]];
Length[virtualLeft]


(* ::Text:: *)
(*The interference 2 Re(M1 M0*), normalized as in the FeynCalc route.  FeynArts writes the one-loop amplitude with the measure d^D loop/(2 Pi)^D; the scalar integrals are defined with the measure d^D loop/(i Pi^(D/2)), so that M1 = i/(4 Pi)^(D/2) x (sum of coefficients c_j times integrals I_j).  The products (coefficient x M0*) formed by the spin sums are real, and*)
(**)
(*   2 Re(M1 M0*) = (2/(4 Pi)^(D/2)) Sum_j [i c_j M0*] Re(I_j),*)
(**)
(*with the real parts of the integrals in the physical region (s > 0) taken in step 3.  The association virtual below holds Sum_j [i c_j M0*] I_j, averaged, for both spin cases, with the loop-momentum hat powers as HatPower[n]; the factor 2/(4 Pi)^(D/2) is included in step 4.*)


(* ::Input:: *)
virtualJobs[spinCase_] := Table[<|"Name" -> "virtual-" <> spinCase <> "-" <> ToString[g], "Kinematics" -> "Born", "Left" -> virtualLeft[[g]],
  "Right" -> bornSumRight, "Spin" -> spinCase, "References" -> bornReferences, "Conservation" -> bornConservation|>, {g, Length[virtualLeft]}];
virtual = If[records["Virtual"] === {}, <|"UU" -> 0, "LL" -> 0|>,
  Association @ Table[spinCase -> Module[{results = runJobs[virtualJobs[spinCase]]},
    Collect[I average[spinCase] Total[MapThread[toFamilies[#1["Value"], #2[[1]]] &, {results, routings}]] /. u -> -s - t,
      {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, colourFactors]], {spinCase, {"UU", "LL"}}]];
Require[FreeQ[virtual, Complex | S | H | epsHead | trace], "The virtual interference is not a real combination of scalar integrals", virtual];
{Length[virtual["UU"]], Union[Cases[virtual, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]] // Length}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*For every real-emission process of the class, the square of the sum of its diagrams with the spin sums of the two cases, in terms of s, t, u, A, B, C and spe3; the processes are added with their weights from classes.wl, which also carry the identical-particle factors (1/2 for two identical unobserved gluons; 1 when the tagged parton has an identical unobserved partner, since the two ways of tagging cancel the 1/2; Nf - 1 or Nf for a quark pair of other flavours).  Reference vectors: pa for the gluons p and k1, pb for k2, the other incoming momentum for an incoming gluon.  In the real kinematics a four-dimensional epsilon tensor Eps[pa, pb, p, k1] can appear (k1 has a four-dimensional part out of the plane of pa, pb, p); it always comes with an imaginary coefficient and drops out of the real part, which the checks at the end of the cell below confirm.  The sum over processes is not brought to a common denominator (step 3 expands it into partial fractions anyway).*)


(* ::Input:: *)
realReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa, k2 -> pb|>;
realSides = Association @ KeyValueMap[#1 -> sidesOf[#2, "Real", realReferences, realConservation] &, records["Real"]];
realSquares = Association @ Table[spinCase -> Association @ KeyValueMap[
    #1 -> square["real-" <> StringReplace[#1, " " -> "_"], #2, "Real", realReferences, realConservation, spinCase] &, realSides], {spinCase, {"UU", "LL"}}];
real = Association @ Table[spinCase -> (average[spinCase] Total[KeyValueMap[weights[#1] #2 &, realSquares[spinCase]]] /. Nc -> CA),
  {spinCase, {"UU", "LL"}}];
Require[FreeQ[real, Complex], "The real-emission square is not real"];
Require[FreeQ[real, S | H | epsHead | trace | word | wordBar | List], "A scalar product or tensor was left in the real-emission square", Union[Cases[real, _S | _H | _epsHead, Infinity]]];
Require[Max[Exponent[#, spe3] & /@ termList[Expand[real["LL"]]]] <= 1, "The real-emission square has spe3^2 terms"];
LeafCount /@ real


(* ::Text:: *)
(*Check of the gauge invariance of the real-emission square, as in the FeynCalc route: for the first process of the class (spin average), the physical polarization sums with other reference vectors (pb for p and k1, pa for k2) give the same result.  The amplitudes are prepared again, since the conditions e_l.r = 0 depend on the reference vectors.  The difference is evaluated exactly at random rational points (all invariants, D, Nc and spe3), with C = -s4 - A - B.*)


(* ::Input:: *)
alternativeReferences = <|pa -> pb, pb -> pa, p -> pb, k1 -> pb, k2 -> pa|>;
firstProcess = First[Keys[records["Real"]]];
alternative = square["gauge-" <> StringReplace[firstProcess, " " -> "_"],
  sidesOf[records["Real"][firstProcess], "Real", alternativeReferences, realConservation], "Real", alternativeReferences, realConservation, "UU"];
gaugePoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], Nc -> 2 + 10 v[[6]], spe3 -> v[[7]]}], {3}];
referenceDifference = Table[(alternative - realSquares["UU"][firstProcess]) /. C -> -(s + t + u) - A - B /. point, {point, gaugePoints}];
Require[referenceDifference === {0, 0, 0}, "The real-emission square depends on the gluon reference vectors", referenceDifference];
"gauge invariance of the real emission checked"


(* ::Section:: *)
(*Timing and save*)


(* ::Text:: *)
(*The number of products computed, and their total and largest time in seconds (from the pair files).*)


(* ::Input:: *)
pairTimes = Association @ Table[FileBaseName[f] -> Get[f]["Seconds"], {f, FileNames["*.wl", pairDirectory]}];
{Length[pairTimes], Round[Total[pairTimes]], Round[Max[pairTimes]]}


(* ::Text:: *)
(*The result, in the format of the FeynCalc route, is saved to Results/Tracer/class-<n>/Algebra.wl, and the parallel kernels are closed.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Virtual" -> virtual, "Real" -> real, "Class" -> class, "Spec" -> spec,
  "Families" -> families, "Variables" -> "s, t, u = -s-t; A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 = -s4-A-B; spe3; loopHat2 via HatPower[n]"|>;
SaveResult["Algebra", algebra, "Tracer/class-" <> class];
CloseKernels[];
Keys[algebra]


(* ::Section:: *)
(*Comparison with the FeynCalc route*)


(* ::Text:: *)
(*The two routes must give the same Algebra.wl.  Born, virtual and real parts are compared, in exact arithmetic, at two random rational points (a sampled comparison, not a symbolic identity) of s, t, u, A, B, D, CA, spe3 and Nf, with C = -(s+t+u) - A - B and CF = (CA^2-1)/(2 CA); in the virtual part every scalar integral K...[n1,n2,n3,n4] is replaced by its own random rational number and HatPower[n] by hp^n with a random rational hp, so that the coefficient of every integral and every hat power is compared.  (If the FeynCalc route has not been run for this class, the comparison is skipped and says so.)*)


(* ::Input:: *)
feynCalcFile = FileNameJoin[{ResultsDirectory["FeynCalc/class-" <> class], "Algebra.wl"}];
comparison = If[!FileExistsQ[feynCalcFile], "FeynCalc result not available: comparison skipped",
  Module[{other = Get[feynCalcFile], integrals, points, valueAt},
    integrals = Union[Cases[{virtual, other["Virtual"]}, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]];
    points = Table[Join[With[{v = RandomInteger[{1, 997}, 9]/RandomInteger[{998, 1997}, 9]},
        {s -> 1 + v[[1]], t -> -v[[2]], u -> -v[[3]], A -> -v[[4]], B -> -v[[5]], D -> 4 - v[[6]], CA -> 2 + 10 v[[7]], spe3 -> v[[8]], Nf -> 1 + 5 v[[9]], hp -> 1 + v[[1]] v[[2]]}],
      Thread[integrals -> RandomInteger[{1, 10^6}, Length[integrals]]/RandomInteger[{1, 10^6}, Length[integrals]]]], {2}];
    valueAt[expr_, point_] := expr /. HatPower[n_] :> hp^n /. point /. C -> -(s + t + u) - A - B /. CF -> (CA^2 - 1)/(2 CA) /. point;
    Association @ Flatten[Table[part <> " " <> sc -> Union[Table[Together[valueAt[algebra[part][sc], point] - valueAt[other[part][sc], point]], {point, points}]],
      {part, {"Born", "Virtual", "Real"}}, {sc, {"UU", "LL"}}]]]];
Require[!AssociationQ[comparison] || Union[Flatten[Values[comparison]]] === {0}, "The Tracer route differs from the FeynCalc route", comparison];
comparison
