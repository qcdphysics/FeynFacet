(* ::Title:: *)
(*Transversely polarized proton-proton classes (A_TT), step 2 (Tracer route): spin and colour algebra*)


(* ::Text:: *)
(*Same task, same input and same output as 2-Algebra-FeynCalc.wl; read its introduction first: it defines the four parts of the output (Born, Virtual, Real and the transversely polarized Born table BornTable), the spin densities of the two incoming partons, the in-plane axis eX and the normal n, the four components XX, NN, XN, NX, the projections TT = (XX - NN)/2 (the cos(2 Phi) coefficient, the hard part of A_TT) and Symmetric = (XX + NN)/2 (evanescent in the HVBM scheme), the reduction of the normal components of the unobserved and loop momenta, and the invariants s, t, u, A, B, C, spe3.  Here the Dirac and Lorentz algebra is done by the classic Tracer package (M. Jamin and M. E. Lautenbacher, Comput. Phys. Commun. 74 (1993) 265) and the colour algebra by the explicit SU(N) reduction of Processes/Hadronic/2-Algebra-Tracer.wl, whose squaring machinery (records -> Tracer words, colour separation, products i <= j, closing of the fermion lines, polarization sums before the traces, parallel kernels) is taken over unchanged; read that file for Tracer's vocabulary and for the method.  Input: the amplitude records of step 1 (Results/class-<n>/Amplitudes.wl, key "Records").  Output: Results/Tracer/class-<n>/Algebra.wl, which the last cell compares exactly with Results/FeynCalc/class-<n>/Algebra.wl.  Run this file in a fresh kernel: Tracer and FeynCalc cannot be loaded together.*)
(**)
(*What differs from the unpolarized and longitudinal classes.  (1) The spin sums of the two incoming partons are the chiral-odd densities pa-slash gamma5 sa-slash/2 and pb-slash gamma5 sb-slash/2 (the same expression for a quark and for an antiquark, docs/TRANSVERSE-PP.md), inserted in Tracer as the string q, G5, s when a fermion line is closed at the spinor of momentum q = pa, pb; the factor 1/2 of each density is part of the average.  (2) The spin vectors are Tracer vectors of their own: the unnormalized in-plane axis aX = (u/s) pa + (t/s) pb + p (aX^2 = -t u/s; the normalization s/(t u) that turns it into the unit vector eX is applied to the XX component afterwards, as in the FeynCalc route) and the normal nn, with nn.pa = nn.pb = nn.p = 0, nn.nn = -1, nn.aX = 0.  Their scalar products with every momentum are declared with the kinematics: aX.x follows from the definition of aX, and the normal components nn.k1 = nr = -nn.k2 (real emission; nn.k1 = 0 in the Born kinematics, where k1 = pa + pb - p lies in the plane) and nn.loop = nLoop (one loop) stay symbolic until the normal reduction, (nn.k1)^2 = k1par^2 + spe3 and (nn.loop)^2 = looppar^2 - loop^2 + loopHat2, with the odd powers dropped.  Both axes are four-dimensional (H = 0), as the HVBM scheme requires of spin vectors.  (3) The four components are computed as four spin cases of the same jobs, and the checks of the FeynCalc route are repeated: XN = NX = 0, XX + NN = 0 in four dimensions, the Born table against the class's own Born, under a <-> b relabelling, charge conjugation, vanishing control channels and the literature, and the gauge invariance of the real emission.*)
(**)
(*Tracer's gamma5 and epsilon tensor have the opposite sign to FeynCalc's (gamma5 = (i/4!) eps gamma gamma gamma gamma, Tr[gamma5 a b c d] = +4 i Eps[a, b, c, d]).  Every trace here contains either two gamma5's (both polarized partons on one fermion line) or one gamma5 paired with a second trace that contains the other one; in the first case the two gamma5's meet and no epsilon tensor is produced, in the second the two epsilon tensors are contracted into scalar products by Tracer as soon as the traces are multiplied.  The sign of Eps therefore drops out of every result, which is tested where it would show: the Born table against the lowest-order cross sections of the literature and the exact comparison with the FeynCalc route at the end.*)


(* ::Section:: *)
(*Load*)


(* ::Text:: *)
(*The FeynFacet package (file locations, the record format, the checks Require[...]) and then Tracer, which redefines Dot and NonCommutativeMultiply on startup; LoadTracer restores them.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadTracer[]


(* ::Text:: *)
(*Choices: the class (5, 6, 7, 11, 12); the number of parallel kernels for the squares ("auto" = one per processor core, at most 8; "1" = none); and whether the pair files of an earlier run of this same, unchanged file may be reused ("no" recomputes everything; "yes" only continues an interrupted run: a product file is reused only if its companion file <name>.input holds the fingerprint of the current inputs, the amplitudes of step 1 and this step file; a file written for other inputs, or without that companion (for example, after an interrupted write), is computed again).*)
(**)
(*The records were written in step 1 by a kernel with FeynCalc loaded, in which the number of light flavours is FeynCalc's symbol Nf (context FeynCalc`); it is renamed to the plain symbol Nf of the FeynCalc route's output.  Tracer is set to D dimensions with the HVBM gamma5 (AntiCommute[off]: gamma5 anticommutes with the four-dimensional gamma matrices and commutes with the hatted ones), and the lines 1, 2, 3 are declared traces: lines 1 and 2 for the (at most two) traces formed by closing the open quark lines of an interference, line 3 for a closed quark loop of a one-loop graph.*)


(* ::Input:: *)
class = "5";
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
(*The first record of the first diagram of each part, as an example of the input: a coefficient (couplings, colour factors ColorT, ColorF, ColorDelta, GluonDelta, propagators Den[q] = 1/q^2, metric tensors, momentum components, polarization vectors PolVector[l, mu, -1] of the outgoing gluons) and the open fermion chains Chain[Ubar or Vbar, gamma matrices, Uspinor or Vspinor].  A one-loop record with a closed quark loop also has "Traces".  The amplitudes contain no gamma5 (the gamma5's of the densities are inserted below), which the complex conjugation of a gamma string by reversal relies on.*)


(* ::Input:: *)
Require[FreeQ[records, GammaFive], "A record contains gamma5; the conjugation below does not handle it"];
ShowRecords[Join[If[records["Born"] === {}, {}, {records["Born"][[1, 1]]}], {First[First[records["Real"]]][[1]]}]]


(* ::Section:: *)
(*Conventions*)


(* ::Text:: *)
(*Both incoming partons must be quarks or antiquarks (a gluon has no transversity).  Each is averaged with 1/(2 Nc): the 1/2 of its density (1 + gamma5 s-slash)/2 and the colour average; colourFactors writes the result with CA = Nc and CF = (Nc^2-1)/(2 Nc).  The spin axes of pa and pb for each component, and the normalization that turns the axis aX into the unit vector eX: XX is divided by -aX^2 = t u/s; XN and NX are computed up to the positive factor Sqrt[s/(t u)], which does not matter since they must vanish.*)


(* ::Input:: *)
incoming = (spec["Born"] /. None -> spec["Real"][[1, 1]])[[1 ;; 2]];       (* the two incoming partons *)
Require[FreeQ[incoming, "g"], "The transversely polarized process needs two incoming quarks or antiquarks", incoming];
components = {"XX", "NN", "XN", "NX"};
spinAxes = <|"XX" -> <|pa -> aX, pb -> aX|>, "NN" -> <|pa -> nn, pb -> nn|>, "XN" -> <|pa -> aX, pb -> nn|>, "NX" -> <|pa -> nn, pb -> aX|>|>;
normalization[spinCase_] := If[spinCase === "XX", s/(t u), 1];
average[spinCase_] := 1/(2 Nc)^2;                                   (* 1/(2 Nc) per incoming quark: the 1/2 of the density and the colour average *)
colourFactors[expr_] := Factor[expr] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
{incoming, average["XX"]}


(* ::Text:: *)
(*Kinematics are installed by functions, because every parallel kernel installs them for each job; they are always called inside Internal`InheritedBlock[{S, H}, ...], which keeps the declared scalar products local to one job while retaining Tracer's own rules for S and H.  The physical momenta pa, pb, p and the two spin axes have no hatted components (H = 0; H is also made linear, since Tracer leaves H of a momentum sum unexpanded).  In the real emission the unobserved k1 and k2 have the hatted products H[k1, k1] = H[k2, k2] = -H[k1, k2] = spe3 (k1hat = -k2hat, since the physical momenta have no hatted parts); in the Born and one-loop kinematics k1 is four-dimensional too, and the loop momentum keeps its hatted part, H[loop, loop] = loopHat2 (see the one-loop section).  Inside the amplitudes k1 = pa + pb - p (Born and one loop) and k2 = pa + pb - p - k1 (real) are eliminated, except in the spinors, whose spin sums give the slashed momentum itself.*)
(**)
(*The scalar products are kept as lists {x, y, x.y} (productOf reads one in either order), from which the products of the axes are built: aX.x = (u/s) pa.x + (t/s) pb.x + p.x for every momentum x and aX.aX from the products of pa, pb, p; nn.pa = nn.pb = nn.p = 0, nn.nn = -1, nn.aX = 0, and the normal components of the unobserved momenta are given explicitly (nn.k1 = 0 in the Born kinematics, where k1 = pa + pb - p lies in the plane; nn.k1 = nr, nn.k2 = -nr in the real kinematics).  The loop momentum of the one-loop graphs has aX.loop = (u/s) loop.pa + (t/s) loop.pb + loop.p and the symbolic normal component nn.loop = nLoop.  A check that the axis is what the introduction says: aX.pa = aX.pb = 0, aX.aX = aX.p = -t u/s, aX.k1 = t u/s (Born kinematics), nn.aX = 0, nn.nn = -1.*)


(* ::Input:: *)
declare[products_List] := (Unprotect[S]; Scan[(S[#[[1]], #[[2]]] = #[[3]]) &, products]; Protect[S]);
fourDimensional[momenta_List] := (Unprotect[H];
  H[x_, y_Plus] := H[x, #] & /@ y; H[x_Plus, y_] := H[#, y] & /@ x;
  H[x_, c_?NumericQ y_] := c H[x, y]; H[c_?NumericQ x_, y_] := c H[x, y];
  Scan[(H[#] = 0; H[#, _] = 0) &, momenta]; Protect[H]);
productOf[products_, x_, y_] := FirstCase[products, {x, y, v_} | {y, x, v_} :> v, 0];
bornProducts = {{pa, pa, 0}, {pb, pb, 0}, {p, p, 0}, {k1, k1, 0}, {pa, pb, s/2}, {pa, p, -t/2}, {pb, p, -u/2}, {pa, k1, -u/2}, {pb, k1, -t/2}, {p, k1, s/2}};
realProducts = {{pa, pa, 0}, {pb, pb, 0}, {p, p, 0}, {k1, k1, 0}, {k2, k2, 0}, {pa, pb, s/2}, {pa, p, -t/2}, {pb, p, -u/2},
  {pa, k1, -A/2}, {pb, k1, -B/2}, {p, k1, C/2}, {pa, k2, (s + t + A)/2}, {pb, k2, (s + u + B)/2}, {p, k2, (-t - u - C)/2}, {k1, k2, (s + t + u)/2}};
axisCoefficients = {u/s, t/s, 1};                      (* aX = (u/s) pa + (t/s) pb + p *)
axisProducts[products_, momenta_List, normalComponents_List] := Module[{base = {pa, pb, p}, ca = axisCoefficients},
  Join[Table[{aX, x, Factor[Sum[ca[[i]] productOf[products, base[[i]], x], {i, 3}]]}, {x, momenta}],
    {{aX, aX, Factor[Sum[ca[[i]] ca[[j]] productOf[products, base[[i]], base[[j]]], {i, 3}, {j, 3}]]}, {aX, nn, 0}, {nn, nn, -1}},
    Table[{nn, x, 0}, {x, base}], normalComponents]];
bornKinematics := (fourDimensional[{pa, pb, p, k1, aX, nn}];
  declare[Join[bornProducts, axisProducts[bornProducts, {pa, pb, p, k1}, {{nn, k1, 0}}],
    {{aX, loop, axisCoefficients . {S[loop, pa], S[loop, pb], S[loop, p]}}, {nn, loop, nLoop}}]]);
realKinematics := (fourDimensional[{pa, pb, p, aX, nn}];
  Unprotect[H]; H[k1, k1] = spe3; H[k2, k2] = spe3; H[k1, k2] = -spe3; Protect[H];
  declare[Join[realProducts, axisProducts[realProducts, {pa, pb, p, k1, k2}, {{nn, k1, nr}, {nn, k2, -nr}}]]]);
installKinematics[name_] := If[name === "Born", bornKinematics, realKinematics];
bornConservation = k1 -> pa + pb - p;
realConservation = k2 -> pa + pb - p - k1;
conservedMomentum[rule_][q_] := Expand[q /. rule];      (* the map of momenta inside an amplitude: eliminate k1 or k2 *)
epsHead = Head[Eps[a1, a2, a3, a4]];          (* Tracer's internal head of the epsilon tensor *)
axisCheck = Internal`InheritedBlock[{S, H}, bornKinematics;
  Simplify[{S[aX, pa], S[aX, pb], S[aX, aX] + t u/s, S[aX, p] + t u/s, S[aX, k1] - t u/s, S[nn, aX], S[nn, nn] + 1, H[aX, k1], H[nn, loop]} /. u -> -s - t]];
Require[axisCheck === {0, 0, 0, 0, 0, 0, 0, 0, 0}, "The in-plane axis or the normal is not what the text says", axisCheck];
"axes checked"


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
(*Checks of the reduction on known colour sums (CF = (Nc^2-1)/(2 Nc)): Tr(T^a T^a) = CF Nc; f^{abc} f^{abc} = Nc (Nc^2-1); (T^a T^b)_{ij} ((T^a T^b)_{ij})^* = CF^2 Nc; (T^a T^b)_{ij} ((T^b T^a)_{ij})^* = -CF/2; f^{abc} (T^c)_{ij} ((T^a T^b)_{ij})^* = f^{abc} Tr[T^c T^b T^a] = -(i/4) f^{abc} f^{abc} = -i Nc (Nc^2-1)/4, which tests the sign of identity (1); the closed loop f^{a1 b1 c1} f^{a1 b2 c2} f^{b1 b2 d} f^{c1 c2 d} = Nc^2 (Nc^2-1)/2; and f^{z1 z2 z3} f^{z1 z2 z3}, which differs from the second network only by the names of its indices and must be found in the cache (the number of cached networks must not grow).  The full test of the colour algebra is the Born table below, against the literature and the FeynCalc route.*)


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
(*Each record becomes a Tracer expression: Metric[mu, nu] -> S[{mu}, {nu}], MomentumComponent[q, mu] -> S[q, {mu}], PolVector[l, mu, _] -> S[e_l, {mu}] with a vector e_l for the polarization of the gluon l (named ep, ek1, ek2; in the complex-conjugate amplitude epStar, ...), Den[q] -> 1/S[q, q] for a tree propagator (evaluated with the installed kinematics) or 1/DjOf[q] for a propagator that carries the loop momentum (see the one-loop section).  An open fermion chain ubar(a) W u(b) becomes the inert object word[a, W..., b], with the spinor momenta at its ends and W the Tracer arguments ({mu} for gamma^mu, the momentum q for q-slash).  A closed quark loop becomes G[3, ...], a trace that Tracer evaluates at once, times its other factors.  The momenta inside the amplitude go through the map mom (momentum conservation, and for a one-loop graph also the shift of the loop momentum).  Because the record's factors are multiplied out, Tracer contracts every index summed inside the amplitude: indices of the metric tensors and momentum components disappear, and those attached to a fermion chain move into its word.*)


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
(*As in Processes/Hadronic (both routes), a square sum_{i,j} M_i M_j^* is not formed as one product:*)
(**)
(*1. Every amplitude (one diagram) is converted and simplified once, before any product, and the conditions e_l.l = 0 and e_l.r = 0 are imposed for every gluon l, where r is the reference vector of the same gluon's physical polarization sum.  For the gluon whose momentum was eliminated by momentum conservation, e_x.x = 0 is imposed in its eliminated form: x = pa + pb - p (- k1) gives e_x.q = -(the other terms) for the first momentum q of x that is not its reference vector.*)
(*2. Colour is separated from the Lorentz algebra: each amplitude is written as sum_x c_x L_x with distinct colour structures c_x; the colour factor of c_x c_y^* is computed by colourValue above.  The complex conjugate of a Lorentz part L_y changes i -> -i, e_l -> e_lStar, renames its summed Lorentz indices, and replaces each word ubar(a) W u(b) by its complex conjugate ubar(b) W^R u(a), with W^R the gamma string in reversed order (wordBar[b, W^R..., a]; reversing is the Dirac adjoint gamma0 W^dagger gamma0 of a string of gamma^mu and slashed momenta).*)
(*3. Only the products with i <= j are computed: for i < j, M_i M_j^* + M_j M_i^* = 2 Re(M_i M_j^*) (all invariants are real).*)
(*4. The fermion lines are closed: in the product of the amplitude's words and the conjugate amplitude's words, the word that ends on the spinor of momentum q is followed by the conjugate word that starts there, with the spin sum in between.  This is where the transverse polarization enters: at an incoming parton q = pa, pb the spin sum is the density q-slash gamma5 s-slash (the string q, G5, s in Tracer, with s the axis of that parton in the component being computed: aX or nn), the same for u(q) ubar(q) and v(q) vbar(q); at an outgoing quark it is q-slash, for u(q) ubar(q) and v(q) vbar(q) alike.  One continues with the next word until the first word comes back: each closed walk is one trace, trace[line, ...], kept inert for now.  The interference has one trace (both polarized partons on one line: two gamma5's) or two (one gamma5 each), depending on how the quark lines of the two diagrams connect.*)
(*5. The gluon polarization sums are applied one gluon at a time, before the traces are taken: the terms are grouped by the factors that carry that gluon's polarization vectors e_l and e_lStar; in each distinct group e_l -> {mu_l}, e_lStar -> {nu_l}, and the product with the physical density -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r) is expanded once, Tracer contracting the indices, also into the inert traces.  Contracting before tracing keeps the traces short: a contracted pair gamma^mu ... gamma_mu inside a trace is reduced by Tracer's D-dimensional contraction formulas before the trace is expanded.  (Every external gluon of these classes is outgoing.)*)
(*6. The traces are then evaluated (each distinct one once) and the products expanded; every remaining scalar product is a kinematic invariant, a normal component nr or nLoop, or a product of the loop momentum.  An epsilon tensor of the eliminated momentum is rewritten with momentum conservation (in the Born kinematics Eps[pa, pb, p, k1] = 0).  The products (i, j) are independent: they run on parallel kernels, and each is written to its own file with the time it took, so that the progress is visible and an interrupted run can be continued.*)


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
(* the spin sum at the spinor of momentum q, as Tracer arguments: the transverse density q-slash gamma5 s-slash at an incoming parton, q-slash at an outgoing one *)
junction[q_, spinCase_] := If[MemberQ[{pa, pb}, q], {q, G5, spinAxes[spinCase][q]}, {q}];
sew[term_, spinCase_] := Module[{factors = factorList[term], words, bars, rest, traces = {}, current, start, chain, q},
  words = Association[#[[1]] -> # & /@ Cases[factors, _word]];          (* each word, keyed by the momentum it starts with *)
  bars = Association[#[[1]] -> # & /@ Cases[factors, _wordBar]];        (* each conjugate word, keyed likewise *)
  rest = Times @@ DeleteCases[factors, _word | _wordBar];
  While[Length[words] > 0,
    current = First[words]; start = current[[1]]; chain = {};
    While[True,
      KeyDropFrom[words, current[[1]]];
      q = Last[current]; chain = Join[chain, List @@ current[[2 ;; -2]], junction[q, spinCase]];
      Require[KeyExistsQ[bars, q], "A fermion line does not close", term];
      current = bars[q]; q = Last[current]; chain = Join[chain, List @@ current[[2 ;; -2]], junction[q, spinCase]];
      If[q === start, Break[]];
      Require[KeyExistsQ[words, q], "A fermion line does not close", term];
      current = words[q]];
    AppendTo[traces, trace @@ Join[{Length[traces] + 1}, chain]]];
  rest Times @@ traces];

(* ---- 5. polarization sums ---- *)
physicalSum[l_, r_] := With[{mu = {Symbol["mu" <> ToString[l]]}, nu = {Symbol["nu" <> ToString[l]]}},
  -S[mu, nu] + (S[l, mu] S[r, nu] + S[r, mu] S[l, nu])/S[l, r]];
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
  Do[expr = project[expr, l, physicalSum[l, references[l]]], {l, gluons}];
  expr = evaluateTraces[expr];
  (* the eliminated momentum inside a four-dimensional epsilon tensor (in the Born kinematics Eps[pa, pb, p, k1] = 0) *)
  Expand[expr /. epsHead[x_List] :> (Eps @@ (x /. conservation))]];
(* the real part: after expanding, every term is a number times real symbols (the denominators are real propagators), so the real part of the number is kept *)
realPartOf[expr_] := Expand[expr] /. Complex[a_, b_] :> a;
imaginaryPartOf[expr_] := (expr - (expr /. Complex[a_, b_] :> Complex[a, -b]))/2;
"squaring functions defined"


(* ::Text:: *)
(*The normal component of an unobserved or loop momentum x, as in the FeynCalc route: the interference is a polynomial in nx = nn.x; its odd powers integrate to zero, and nx^2 = xpar^2 - x4.x4, with xpar^2 = d.G^-1.d (d the products of x with pa, pb, p, G their Gram matrix) and x4.x4 the four-dimensional square of x (x = xpar + c nn + xhat, nx = -c, x4.x4 = xpar^2 - c^2).  For k1 in the real emission x4.x4 = -spe3, so nr^2 = k1par^2 + spe3; for the loop momentum x4.x4 = loop.loop - loopHat2.  The part odd in nx is returned separately (the virtual and real sections report its size and drop it).*)


(* ::Input:: *)
gramMatrix[products_] := Table[productOf[products, x, y], {x, {pa, pb, p}}, {y, {pa, pb, p}}];
parallelSquare[d_List, gram_] := Factor[d . Inverse[gram] . d];
normalReduction[expr_, nx_, parallel_, fourSquare_] := Module[{coefficients = CoefficientList[Expand[expr], nx]},
  If[coefficients === {}, coefficients = {0}];                                   (* a vanishing expression (a graph can vanish exactly) *)
  Require[Length[coefficients] <= 3, "A power of the normal component above two appeared", Length[coefficients] - 1];
  <|"Even" -> coefficients[[1]] + If[Length[coefficients] == 3, coefficients[[3]] (parallel - fourSquare), 0],
    "Odd" -> If[Length[coefficients] >= 2, coefficients[[2]], 0]|>];
realParallelSquare = parallelSquare[{productOf[realProducts, pa, k1], productOf[realProducts, pb, k1], productOf[realProducts, p, k1]}, gramMatrix[realProducts]];     (* k1par^2 *)
loopParallelSquare = parallelSquare[{S[loop, pa], S[loop, pb], S[loop, p]}, gramMatrix[bornProducts]];                                                             (* looppar^2 *)
realParallelSquare


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
(*Each parallel kernel sets the dimension, the HVBM gamma5 and the traced lines, and receives the functions above with empty caches.  A job is one product (i, j) of the squaring method for one spin component: it installs its kinematics (Born or real) inside Internal`InheritedBlock, computes the interference, and is written to Results/Tracer/class-<n>/pairs/<name>.wl with its time.  square adds the products of one square: a diagonal product as it is, a product with i < j twice its real part (step 3).  A diagonal product |M_i|^2 is real, but its epsilon-tensor terms may leave an imaginary part that vanishes only after simplification: that imaginary part is checked to be exactly zero at a rational point (checkPoint) and dropped.  sidesOf prepares every diagram once in the given kinematics, separates its colour, and forms its conjugate.*)


(* ::Input:: *)
If[nWorkers > 1,
  ParallelEvaluate[VectorDimension[D]; AntiCommute[off]; Spur[1, 2, 3]];
  DistributeDefinitions[declare, fourDimensional, productOf, bornProducts, realProducts, axisCoefficients, axisProducts, bornKinematics, realKinematics, installKinematics, epsHead,
    colourHeads, factorList, termList, colourIndices, internalColour, conjugateColour, canonicalColour, freshIndex, adjointDelta, toGeneratorStrings, colourRules, colourValue,
    polarization, polarizationStar, spinAxes, junction, sew, physicalSum, tensorPart, projectOnce, project, evaluateTraces, gluonsOf, interference];
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
checkPoint = {s -> 1, t -> -3/10, u -> -2/5, A -> -1/7, B -> -2/9, C -> -3/10 + 1/7 + 2/9, D -> 37/10, Nc -> 3, spe3 -> 1/11, nr -> 3/13, gs -> 1, Nf -> 5};
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
(*The Born square with the reference vectors of the FeynCalc route (pa for the gluons p and k1), for the four components, averaged over the incoming spins and colours, normalized (XX divided by -aX^2) and written with CA, CF and u = -s - t.  Then the checks XN = NX = 0 (exactly) and XX + NN = 0 in four dimensions, and the two projections TT = (XX - NN)/2 and Symmetric = (XX + NN)/2.  The evanescent symmetric part, proportional to (D - 4), is displayed.  No identical-particle factor: the tagged parton p is distinguished from the unobserved one.*)


(* ::Input:: *)
bornReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa|>;
bornSides = If[records["Born"] === {}, None, sidesOf[records["Born"], "Born", bornReferences, bornConservation]];
bornComponents = If[bornSides === None, Association[# -> 0 & /@ components],
  Association @ Table[sc -> With[{sq = square["born", bornSides, "Born", bornReferences, bornConservation, sc]},
    Factor[colourFactors[Together[average[sc] normalization[sc] sq]] /. u -> -s - t]], {sc, components}]];
Require[FreeQ[bornComponents, Complex | S | H | epsHead | trace | word | wordBar | nr | nLoop], "The Born square is not a real scalar", bornComponents];
Require[bornComponents["XN"] === 0 && bornComponents["NX"] === 0, "The mixed components of the Born square do not vanish", bornComponents];
Require[Factor[(bornComponents["XX"] + bornComponents["NN"]) /. D -> 4] === 0, "XX + NN of the Born square does not vanish in four dimensions", bornComponents];
born = <|"TT" -> Factor[(bornComponents["XX"] - bornComponents["NN"])/2], "Symmetric" -> Factor[(bornComponents["XX"] + bornComponents["NN"])/2],
  "XX" -> bornComponents["XX"], "NN" -> bornComponents["NN"]|>;
born


(* ::Section:: *)
(*The transversely polarized Born table*)


(* ::Text:: *)
(*The D-dimensional Born squares |M|^2/gs^4 (colours and spins of the unobserved parton summed, colours of a and b averaged, the two transverse densities inserted) of every 2 -> 2 channel a b -> c that the collinear counterterms of the five classes need: q qb -> q' (and qb'), q q -> q, q qb -> q (and qb), q qb -> g, their a <-> b relabellings (qb q -> ...), and the charge conjugate qb qb -> qb; plus the control channels q q' -> q, q q' -> q', q qb' -> q, q qb' -> qb', which must vanish (classes 1-4: the two polarized partons on different fermion lines).  The FeynCalc route generates these amplitudes with FeynArts; here, where FeynArts' output cannot be converted (that needs FeynCalc), the tree amplitudes are written by hand in the record format of step 1 and squared with the machinery above, which makes the table an independent check of both.  The records follow the conventions of the step-1 records (the quark-gluon vertex -i gs gamma^mu T^a, Den[q] = 1/q^2, with the i's of the propagators absorbed in the coefficient; the overall phase of an amplitude is irrelevant for its square): the exchange of one gluon between two quark lines, with the colour generators T^a_{c3 c1} T^a_{c4 c2} (c1, c2 the colours of pa, pb; c3, c4 those of p and the unobserved parton k1); the annihilation of q qb into a gluon, T^a_{c2 c1} T^a_{c3 c4}; and for q qb -> g g the three-gluon vertex, gs f^{abc} [g^{mu nu} (k1 - p)^rho + g^{nu rho} (2 k1 + p)^mu - g^{rho mu} (k1 + 2 p)^nu] with the quark line's T^c, and the two quark-exchange graphs with T^a T^b in both orders.  The relative sign between the two graphs of q q -> q q (t- and u-channel exchange) and of q qb -> q qb (annihilation and exchange) is -1 (Fermi statistics).  The derived channels: observing the other outgoing parton is the renaming p <-> k1; a <-> b is pa <-> pb; charge conjugation reverses every quark line, ubar(p) gamma^mu u(pa) -> vbar(pa) gamma^mu v(p), transposes its colour generators and gives a sign -1 per quark-gluon vertex (an even number here).  Keys as in Processes/PartonicBorn ("q qb -> qp"); each entry holds XX, NN, TT and Symmetric.*)
(**)
(*Checks: XN = NX = 0 and XX + NN = 0 in four dimensions for every channel; the class's own Born square (from the FeynArts records of step 1) equals its table entry (from the hand-written records); the a <-> b relabelled channel equals the original with t <-> u (the axis aX is symmetric, the normal changes sign for both densities); the charge conjugate equals its representative (the densities are the same for quarks and antiquarks); the control channels vanish; and at D = 4, Nc = 3 the TT entries reproduce Table I of Soffer, Stratmann and Vogelsang: with sigma = s E dsigma/d^3p = 2 s^2 dsigma/(dt du dPhi) and dsigma/(dt du dPhi) = dsigma/(dt dPhi) delta(s4), the cos(2 Phi) coefficient of their Eq. (6) corresponds to |M|^2_TT/gs^4 = (4/9) (t u/s^2) delta|M|^2.*)


(* ::Input:: *)
(* the hand-written tree amplitudes: a list of diagrams, each a list of records *)
tableRecords = <|
  "q qb -> qp" -> {
    {<|"Coefficient" -> -I gs^2 ColorT[{Glu5}, Col2, Col1] ColorT[{Glu5}, Col3, Col4] Den[k1 + p] Metric[Lor1, Lor2],
       "Chains" -> {Chain[Ubar[p], GammaMatrix[Lor2], Vspinor[k1]], Chain[Vbar[pb], GammaMatrix[Lor1], Uspinor[pa]]}|>}},
  "q q -> q" -> {
    {<|"Coefficient" -> -I gs^2 ColorT[{Glu5}, Col3, Col1] ColorT[{Glu5}, Col4, Col2] Den[k1 - pb] Metric[Lor1, Lor2],
       "Chains" -> {Chain[Ubar[k1], GammaMatrix[Lor2], Uspinor[pb]], Chain[Ubar[p], GammaMatrix[Lor1], Uspinor[pa]]}|>},
    {<|"Coefficient" -> I gs^2 ColorT[{Glu5}, Col3, Col2] ColorT[{Glu5}, Col4, Col1] Den[p - pb] Metric[Lor1, Lor2],
       "Chains" -> {Chain[Ubar[k1], GammaMatrix[Lor1], Uspinor[pa]], Chain[Ubar[p], GammaMatrix[Lor2], Uspinor[pb]]}|>}},
  "q qb -> q" -> {
    {<|"Coefficient" -> -I gs^2 ColorT[{Glu5}, Col2, Col1] ColorT[{Glu5}, Col3, Col4] Den[k1 + p] Metric[Lor1, Lor2],
       "Chains" -> {Chain[Ubar[p], GammaMatrix[Lor2], Vspinor[k1]], Chain[Vbar[pb], GammaMatrix[Lor1], Uspinor[pa]]}|>},
    {<|"Coefficient" -> I gs^2 ColorT[{Glu5}, Col2, Col4] ColorT[{Glu5}, Col3, Col1] Den[k1 - pb] Metric[Lor1, Lor2],
       "Chains" -> {Chain[Ubar[p], GammaMatrix[Lor1], Uspinor[pa]], Chain[Vbar[pb], GammaMatrix[Lor2], Vspinor[k1]]}|>}},
  "q qb -> g" -> {
    {<|"Coefficient" -> -gs^2 ColorF[Glu3, Glu4, Glu5] ColorT[{Glu5}, Col2, Col1] Den[k1 + p] Metric[Lor2, Lor4] Metric[Lor3, Lor4] MomentumComponent[2 k1 + p, Lor1] PolVector[k1, Lor2, -1] PolVector[p, Lor1, -1],
       "Chains" -> {Chain[Vbar[pb], GammaMatrix[Lor3], Uspinor[pa]]}|>,
     <|"Coefficient" -> -gs^2 ColorF[Glu3, Glu4, Glu5] ColorT[{Glu5}, Col2, Col1] Den[k1 + p] Metric[Lor1, Lor4] Metric[Lor3, Lor4] MomentumComponent[-k1 - 2 p, Lor2] PolVector[k1, Lor2, -1] PolVector[p, Lor1, -1],
       "Chains" -> {Chain[Vbar[pb], GammaMatrix[Lor3], Uspinor[pa]]}|>,
     <|"Coefficient" -> -gs^2 ColorF[Glu3, Glu4, Glu5] ColorT[{Glu5}, Col2, Col1] Den[k1 + p] Metric[Lor1, Lor2] Metric[Lor3, Lor4] MomentumComponent[-k1 + p, Lor4] PolVector[k1, Lor2, -1] PolVector[p, Lor1, -1],
       "Chains" -> {Chain[Vbar[pb], GammaMatrix[Lor3], Uspinor[pa]]}|>},
    {<|"Coefficient" -> -I gs^2 ColorT[{Glu3}, Col5, Col1] ColorT[{Glu4}, Col2, Col5] Den[k1 - pb] PolVector[k1, Lor2, -1] PolVector[p, Lor1, -1],
       "Chains" -> {Chain[Vbar[pb], GammaMatrix[Lor2], Slash[k1 - pb], GammaMatrix[Lor1], Uspinor[pa]]}|>},
    {<|"Coefficient" -> -I gs^2 ColorT[{Glu3}, Col2, Col5] ColorT[{Glu4}, Col5, Col1] Den[p - pb] PolVector[k1, Lor2, -1] PolVector[p, Lor1, -1],
       "Chains" -> {Chain[Vbar[pb], GammaMatrix[Lor1], Slash[p - pb], GammaMatrix[Lor2], Uspinor[pa]]}|>}},
  "q qp -> q" -> {
    {<|"Coefficient" -> -I gs^2 ColorT[{Glu5}, Col3, Col1] ColorT[{Glu5}, Col4, Col2] Den[k1 - pb] Metric[Lor1, Lor2],
       "Chains" -> {Chain[Ubar[k1], GammaMatrix[Lor2], Uspinor[pb]], Chain[Ubar[p], GammaMatrix[Lor1], Uspinor[pa]]}|>}},
  "q qpb -> q" -> {
    {<|"Coefficient" -> I gs^2 ColorT[{Glu5}, Col2, Col4] ColorT[{Glu5}, Col3, Col1] Den[k1 - pb] Metric[Lor1, Lor2],
       "Chains" -> {Chain[Ubar[p], GammaMatrix[Lor1], Uspinor[pa]], Chain[Vbar[pb], GammaMatrix[Lor2], Vspinor[k1]]}|>}}|>;
(* the derived channels *)
chargeConjugate[diagrams_] := diagrams /. {Chain[Ubar[a_], w___, Uspinor[b_]] :> Chain[Vbar[b], Sequence @@ Reverse[{w}], Vspinor[a]], ColorT[l_List, i_, j_] :> ColorT[Reverse[l], j, i]};
Do[tableRecords[pair[[2]]] = tableRecords[pair[[1]]] /. {p -> k1, k1 -> p}, {pair, {{"q qb -> qp", "q qb -> qpb"}, {"q qb -> q", "q qb -> qb"}, {"q qp -> q", "q qp -> qp"}, {"q qpb -> q", "q qpb -> qpb"}}}];
Do[tableRecords["qb q -> " <> c] = tableRecords["q qb -> " <> c] /. {pa -> pb, pb -> pa}, {c, {"qp", "qpb", "q", "qb", "g"}}];
tableRecords["qb qb -> qb"] = chargeConjugate[tableRecords["q q -> q"]];
tableChannels = {"q qb -> qp", "q qb -> qpb", "q q -> q", "q qb -> q", "q qb -> qb", "q qb -> g", "qb q -> qp", "qb q -> qpb", "qb q -> q", "qb q -> qb", "qb q -> g", "qb qb -> qb"};
controlChannels = {"q qp -> q", "q qp -> qp", "q qpb -> q", "q qpb -> qpb"};
jobName[name_] := StringReplace[name, {" -> " -> "_to_", " " -> "_"}];
tableEntry[name_] := Module[{sides = sidesOf[tableRecords[name], "Born", bornReferences, bornConservation], xx, nnv, xn, nx},
  {xx, nnv, xn, nx} = Table[With[{sq = square["table-" <> jobName[name], sides, "Born", bornReferences, bornConservation, sc]},
    Factor[colourFactors[Together[average[sc] normalization[sc] sq/gs^4]] /. u -> -s - t]], {sc, components}];
  Require[FreeQ[{xx, nnv, xn, nx}, Complex | S | H | epsHead | trace | word | wordBar | nr | gs], "A table entry is not a real scalar", name];
  Require[xn === 0 && nx === 0, "The mixed components of a table entry do not vanish", {name, xn, nx}];
  Require[Factor[(xx + nnv) /. D -> 4] === 0, "XX + NN of a table entry does not vanish in four dimensions", {name, xx, nnv}];
  <|"XX" -> xx, "NN" -> nnv, "TT" -> Factor[(xx - nnv)/2], "Symmetric" -> Factor[(xx + nnv)/2]|>];
tableEntries = Association @ Table[name -> tableEntry[name], {name, Join[tableChannels, controlChannels]}];
bornTable = Association @ Table[key -> Map[#[key] &, tableEntries], {key, {"TT", "Symmetric", "XX", "NN"}}];
TableForm[Normal[Map[Simplify[# /. D -> 4] &, bornTable["TT"]]]]


(* ::Input:: *)
tableCheck[expr_] := Factor[expr /. u -> -s - t];
Require[And @@ Table[tableCheck[bornTable["TT"][name]] === 0 && tableCheck[bornTable["Symmetric"][name]] === 0, {name, controlChannels}],
  "A control channel of the Born table does not vanish", Map[bornTable["TT"], controlChannels]];
(* a <-> b: pa <-> pb exchanges t and u, i.e. t -> -s - t in the entries (u is eliminated) *)
swapped[name_] := Module[{parts = StringSplit[name, " "]}, parts[[2]] <> " " <> parts[[1]] <> " -> " <> parts[[4]]];
relabelDifferences = Association @ Flatten[Table[{key, name} -> tableCheck[bornTable[key][swapped[name]] - (bornTable[key][name] /. t -> -s - t)],
    {key, {"TT", "Symmetric"}}, {name, Select[tableChannels, StringStartsQ[#, "q qb"] &]}]];
Require[And @@ (# === 0 & /@ Values[relabelDifferences]), "An a <-> b relabelled entry is not the t <-> u image of the original", Select[relabelDifferences, # =!= 0 &]];
Require[tableCheck[bornTable["TT"]["qb qb -> qb"] - bornTable["TT"]["q q -> q"]] === 0 && tableCheck[bornTable["Symmetric"]["qb qb -> qb"] - bornTable["Symmetric"]["q q -> q"]] === 0,
  "The charge-conjugate channel differs from its representative"];
(* the class's own Born square (the FeynArts records of step 1) against the table (the hand-written records) *)
If[spec["Born"] =!= None,
  Require[And @@ Table[tableCheck[born[key]/gs^4 - bornTable[key][StringRiffle[spec["Born"][[1 ;; 2]], " "] <> " -> " <> spec["Born"][[3]]]] === 0, {key, {"TT", "Symmetric", "XX", "NN"}}],
    "The Born square of the class differs from its table entry", born]];
(* the literature: Table I of Soffer, Stratmann and Vogelsang at D = 4, Nc = 3 *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "SSV-born.wl"}]];
literatureDifferences = Association @ Table[name -> Factor[(bornTable["TT"][name] /. {D -> 4, CA -> 3, CF -> 4/3} /. u -> -s - t) - ((4/9) (t u/s^2) SSVMatrixElement[name] /. u -> -s - t)],
  {name, Keys[SSVMatrixElement]}];
Require[And @@ (# === 0 & /@ Values[literatureDifferences]), "A transversely polarized Born square differs from Soffer, Stratmann and Vogelsang", literatureDifferences];
"Born table checked: control channels vanish, a <-> b relabelling, charge conjugation, the class Born, and the literature at D = 4, Nc = 3"


(* ::Section:: *)
(*Virtual: the one-loop interference in terms of scalar integrals*)


(* ::Text:: *)
(*The families, the routing of the loop momentum and the reduction of the loop momentum's scalar products are those of Processes/Hadronic (both routes):*)
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
(*For each graph: convert it with the map mom = (shift of the loop momentum, then k1 -> pa + pb - p), so that each loop propagator becomes 1/DjOf[q] with q one of the family's momenta up to sign (a tree propagator next to the loop is evaluated); prepare it once and separate its colour; its products with the conjugate of the summed Born amplitude are the jobs below (one per graph and spin component).*)
(**)
(*After the spin sums and the traces the loop momentum appears through S[loop, loop], S[loop, pa], S[loop, pb], S[loop, p] (S[loop, k1], from the spin sum of an outgoing quark k1, is rewritten with k1 = pa + pb - p), the hatted H[loop, loop] = loopHat2, and the normal component nLoop = nn.loop (aX.loop was declared as its combination of loop.pa, loop.pb, loop.p).  Tracer's S is D-dimensional and H its (D-4)-dimensional part, so no four-dimensional product needs to be rewritten.  The normal component is reduced first: its odd powers are dropped (the tensor loop integrals are built from pa, pb, p only, which have no normal component), and nLoop^2 = looppar^2 - loop.loop + loopHat2 with looppar^2 = d.G^-1.d, d = (loop.pa, loop.pb, loop.p).  Then the four D-dimensional products are expressed through the inverse propagators Dj of the family (the four equations loop^2 = D1, 2 loop.q_j = D_j - D_1 - q_j^2 are solved for loop^2, loop.pa, loop.pb, loop.p), and each term is read off as coefficient x family[indices] x HatPower[n].  loopHat2 shifts the dimension of the integral in step 3, not the dimension D of the coefficients.*)


(* ::Input:: *)
loopProducts = {L2 -> S[loop, loop], Lpa -> S[loop, pa], Lpb -> S[loop, pb], Lp -> S[loop, p]};
inverseProps[name_] := Module[{q = families[name], sol},
  sol = Internal`InheritedBlock[{S, H}, bornKinematics;
    Solve[Table[Dj[j] == Expand[S[q[[j]], q[[j]]]] /. Reverse /@ loopProducts, {j, 4}] /. u -> -s - t, {L2, Lpa, Lpb, Lp}]];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators"]; First[sol]];
loopRules[name_] := loopRules[name] = With[{sol = inverseProps[name]}, Table[(x /. loopProducts) -> (x /. sol), {x, {L2, Lpa, Lpb, Lp}}]];
shiftedMomentum[{name_, sigma_, delta_}][q_] := Expand[q /. loop -> sigma loop + delta /. bornConservation];
(* the contracted interference of one graph -> sum of coefficient x family[indices] x HatPower[n]; the part odd in nn.loop is dropped and its size reported *)
toFamilies[expr0_, name_] := Module[{expr, reduced},
  expr = Expand[expr0 /. S[loop, k1] :> S[loop, pa + pb - p] /. H[loop, loop] -> loopHat2];
  reduced = normalReduction[expr, nLoop, loopParallelSquare, S[loop, loop] - loopHat2];
  expr = Expand[reduced["Even"] /. loopRules[name]];
  Require[FreeQ[expr /. _DjOf -> 1, loop], "A product of the loop momentum was not expressed through the propagators",
    Short[Union[Cases[expr, h_[a___] /; !FreeQ[{a}, loop] && FreeQ[h, Times | Plus | Power | DjOf], Infinity]], 8]];
  expr = Expand[expr /. DjOf[q_] :> Dj @@ (Position[families[name], f_ /; sameUpToSign[f, q]][[1]])];
  (* every term is a coefficient times powers of the Dj (negative for propagators, positive for numerators) and of loopHat2 *)
  <|"Even" -> Total[Table[With[{powers = Table[Exponent[term, Dj[j]], {j, 4}], hat = Exponent[term, loopHat2]},
      (term/(Times @@ Table[Dj[j]^powers[[j]], {j, 4}] loopHat2^hat)) familyHeads[name] @@ (-powers) HatPower[hat]], {term, DeleteCases[termList[expr], 0]}]],
    "OddLeafCount" -> LeafCount[reduced["Odd"]]|>];
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
(*with the real parts of the integrals in the physical region (s > 0) taken in step 3.  The association virtual below holds Sum_j [i c_j M0*] I_j per component, averaged and normalized, with the loop-momentum hat powers as HatPower[n]; the factor 2/(4 Pi)^(D/2) is included in step 4.  Checks: XN and NX vanish (their even parts in nn.loop; the odd parts, which integrate to zero, are dropped), and XX + NN vanishes in four dimensions once the hatted terms (which integrate to order eps) are dropped.*)


(* ::Input:: *)
virtualJobs[spinCase_] := Table[<|"Name" -> "virtual-" <> spinCase <> "-" <> ToString[g], "Kinematics" -> "Born", "Left" -> virtualLeft[[g]],
  "Right" -> bornSumRight, "Spin" -> spinCase, "References" -> bornReferences, "Conservation" -> bornConservation|>, {g, Length[virtualLeft]}];
virtualComponents = If[records["Virtual"] === {}, Association[# -> 0 & /@ components],
  Association @ Table[spinCase -> Module[{results = runJobs[virtualJobs[spinCase]], graphs},
    graphs = MapThread[toFamilies[#1["Value"], #2[[1]]] &, {results, routings}];
    oddLeafCounts[spinCase] = Total[#["OddLeafCount"] & /@ graphs];
    Collect[I average[spinCase] normalization[spinCase] Total[#["Even"] & /@ graphs] /. u -> -s - t,
      {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, colourFactors]], {spinCase, components}]];
Require[FreeQ[virtualComponents, Complex | S | H | epsHead | trace | nLoop], "The virtual interference is not a real combination of scalar integrals or still contains nn.loop", Short[virtualComponents, 5]];
Require[virtualComponents["XN"] === 0 && virtualComponents["NX"] === 0, "The mixed components of the virtual interference do not vanish", {virtualComponents["XN"], virtualComponents["NX"]}];
Require[Collect[(virtualComponents["XX"] + virtualComponents["NN"]) /. HatPower[n_ /; n > 0] -> 0 /. HatPower[0] -> 1 /. D -> 4, {_KVirtualST, _KVirtualSU, _KVirtualTU}, Factor] === 0,
  "XX + NN of the virtual interference does not vanish in four dimensions"];
virtual = <|"TT" -> Collect[(virtualComponents["XX"] - virtualComponents["NN"])/2, {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, Factor],
  "Symmetric" -> Collect[(virtualComponents["XX"] + virtualComponents["NN"])/2, {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, Factor]|>;
{Length[virtual["TT"]], Union[Cases[virtual, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]] // Length, Map[Union[Cases[#, _HatPower, Infinity]] &, virtual],
  "leaf count of the dropped odd parts in nn.loop" -> If[records["Virtual"] === {}, None, oddLeafCounts /@ components]}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*For every real-emission process of the class, the square of the sum of its diagrams for the four components, in terms of the invariants s, t, u, A, B, C, the hatted invariant spe3 and the normal component nr; the processes are added with their weights from classes.wl, which also carry the identical-particle factors (1/2 for two identical unobserved gluons; 1 when the tagged parton has an identical unobserved partner, since the two ways of tagging cancel the 1/2; Nf - 1 for a quark pair of other flavours).  Reference vectors of the physical polarization sums: pa for the gluons p and k1, pb for k2.  Then the normal reduction: the parts odd in nr are dropped (they integrate to zero), nr^2 is replaced by k1par^2 + spe3.  Checks: the even parts of XN and NX vanish, and XX + NN vanishes in four dimensions without the hatted terms; both are evaluated exactly at random rational points of all invariants (the expressions have tens of thousands of terms).  The sum over processes is not brought to a common denominator (step 3 expands it into partial fractions anyway); the colour factors are polynomials in Nc = CA.*)


(* ::Input:: *)
realReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa, k2 -> pb|>;
realSides = Association @ KeyValueMap[#1 -> sidesOf[#2, "Real", realReferences, realConservation] &, records["Real"]];
realSquares = Association @ Table[spinCase -> Association @ KeyValueMap[
    #1 -> square["real-" <> StringReplace[#1, " " -> "_"], #2, "Real", realReferences, realConservation, spinCase] &, realSides], {spinCase, components}];
Require[FreeQ[realSquares, S | H | epsHead | trace | word | wordBar | List], "A scalar product or tensor was left in a real-emission square", Union[Cases[realSquares, _S | _H | _epsHead, Infinity]]];
realReduced = Association @ Table[spinCase -> Association @ KeyValueMap[
    #1 -> normalReduction[#2, nr, realParallelSquare, -spe3] &, realSquares[spinCase]], {spinCase, components}];       (* k1_4.k1_4 = -spe3: nr^2 = k1par^2 + spe3 *)
realComponents = Association @ Table[spinCase -> (average[spinCase] normalization[spinCase] Total[KeyValueMap[weights[#1] #2["Even"] &, realReduced[spinCase]]] /. Nc -> CA), {spinCase, components}];
Require[FreeQ[realComponents, Complex | nr | nn | aX], "The real-emission square is not real or still contains the normal component"];
{realParallelSquare, Map[LeafCount, realComponents]}


(* ::Input:: *)
(* exact rational points of all invariants (C = -s4 - A - B), D, Nc, spe3 *)
realPoints = Table[With[{r = RandomInteger[{1, 97}, 8]/RandomInteger[{98, 197}, 8]},
  {s -> 1, t -> -r[[1]], u -> -r[[2]] (1 - r[[1]]), A -> -r[[3]], B -> -r[[4]], D -> 4 - r[[5]], CA -> 2 + 10 r[[6]], CF -> 1 + r[[7]], spe3 -> r[[8]], Nf -> 5}], {3}];
atPoint[expr_, point_] := expr /. C -> -(s + t + u) - A - B /. point;
Require[And @@ Table[atPoint[realComponents["XN"], point] === 0 && atPoint[realComponents["NX"], point] === 0, {point, realPoints}],
  "The even parts of the mixed components of the real emission do not vanish"];
Require[And @@ Table[atPoint[(realComponents["XX"] + realComponents["NN"]) /. {D -> 4, spe3 -> 0}, point] === 0, {point, realPoints}],
  "XX + NN of the real emission does not vanish in four dimensions"];
real = <|"TT" -> (realComponents["XX"] - realComponents["NN"])/2, "Symmetric" -> (realComponents["XX"] + realComponents["NN"])/2|>;
Require[Max[Exponent[#, spe3] & /@ termList[Expand[real["TT"] + real["Symmetric"]]]] <= 1, "The real-emission square has spe3^2 terms (step 3 has the masters in D + 2 only)"];
{"odd parts dropped (leaf counts)" -> Map[Total[KeyValueMap[LeafCount[#2["Odd"]] &, #]] &, realReduced], LeafCount /@ real}


(* ::Text:: *)
(*Check of the gauge invariance of the real-emission square: for the first process of the class and the XX component (the one without normal components), the physical polarization sums with other reference vectors (pb for p and k1, pa for k2) give the same result.  The amplitudes are prepared again, since the conditions e_l.r = 0 depend on the reference vectors.  The difference is evaluated exactly at random rational points (all invariants, D and Nc), with C = -s4 - A - B.*)


(* ::Input:: *)
alternativeReferences = <|pa -> pb, pb -> pa, p -> pb, k1 -> pb, k2 -> pa|>;
firstProcess = First[Keys[records["Real"]]];
alternative = square["gauge-" <> StringReplace[firstProcess, " " -> "_"],
  sidesOf[records["Real"][firstProcess], "Real", alternativeReferences, realConservation], "Real", alternativeReferences, realConservation, "XX"];
gaugePoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], Nc -> 2 + 10 v[[6]], spe3 -> v[[7]]}], {3}];
referenceDifference = Table[(alternative - realSquares["XX"][firstProcess]) /. C -> -(s + t + u) - A - B /. point, {point, gaugePoints}];
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
(*The result, in the format of the FeynCalc route (the Born square, the virtual interference and the real-emission square for the projections TT and Symmetric, and the Born table), is saved to Results/Tracer/class-<n>/Algebra.wl for steps 3 and 4, and the parallel kernels are closed.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Virtual" -> virtual, "Real" -> real, "BornTable" -> bornTable, "Class" -> class, "Spec" -> spec,
  "Families" -> families, "Projections" -> "TT = (XX - NN)/2 (the cos(2 Phi) coefficient, the hard part of A_TT); Symmetric = (XX + NN)/2 (evanescent in the HVBM scheme)",
  "Variables" -> "s, t, u = -s-t; A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 = -s4-A-B; spe3; loopHat2 via HatPower[n]; Born table: |M|^2/gs^4 in s, t, D, CF, CA"|>;
SaveResult["Algebra", algebra, "Tracer/class-" <> class];
CloseKernels[];
Keys[algebra]


(* ::Section:: *)
(*Comparison with the FeynCalc route*)


(* ::Text:: *)
(*The two routes must give the same Algebra.wl.  Born, virtual and real parts (projections TT and Symmetric) are compared, in exact arithmetic, at two random rational points (a sampled comparison, not a symbolic identity) of s, t, u, A, B, D, CA, spe3 and Nf, with C = -(s+t+u) - A - B and CF = (CA^2-1)/(2 CA); in the virtual part every scalar integral K...[n1,n2,n3,n4] is replaced by its own random rational number and HatPower[n] by hp^n with a random rational hp, so that the coefficient of every integral and every hat power is compared.  The Born table is compared entry by entry for its four keys (TT, Symmetric, XX, NN) and all channels, at the same points.  (If the FeynCalc route has not been run for this class, the comparison is skipped and says so.)*)


(* ::Input:: *)
feynCalcFile = FileNameJoin[{ResultsDirectory["FeynCalc/class-" <> class], "Algebra.wl"}];
comparison = If[!FileExistsQ[feynCalcFile], "FeynCalc result not available: comparison skipped",
  Module[{other = Get[feynCalcFile], integrals, points, valueAt},
    integrals = Union[Cases[{virtual, other["Virtual"]}, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]];
    points = Table[Join[With[{v = RandomInteger[{1, 997}, 9]/RandomInteger[{998, 1997}, 9]},
        {s -> 1 + v[[1]], t -> -v[[2]], u -> -v[[3]], A -> -v[[4]], B -> -v[[5]], D -> 4 - v[[6]], CA -> 2 + 10 v[[7]], spe3 -> v[[8]], Nf -> 1 + 5 v[[9]], hp -> 1 + v[[1]] v[[2]]}],
      Thread[integrals -> RandomInteger[{1, 10^6}, Length[integrals]]/RandomInteger[{1, 10^6}, Length[integrals]]]], {2}];
    valueAt[expr_, point_] := expr /. HatPower[n_] :> hp^n /. point /. C -> -(s + t + u) - A - B /. CF -> (CA^2 - 1)/(2 CA) /. point;
    Join[
      Association @ Flatten[Table[part <> " " <> sc -> Union[Table[Together[valueAt[algebra[part][sc], point] - valueAt[other[part][sc], point]], {point, points}]],
        {part, {"Born", "Virtual", "Real"}}, {sc, {"TT", "Symmetric"}}]],
      Association @ Table["BornTable " <> key -> Union[Flatten[Table[Together[valueAt[bornTable[key][name], point] - valueAt[other["BornTable"][key][name], point]],
        {name, Keys[other["BornTable"][key]]}, {point, points}]]], {key, {"TT", "Symmetric", "XX", "NN"}}]]]];
Require[!AssociationQ[comparison] || Union[Flatten[Values[comparison]]] === {0}, "The Tracer route differs from the FeynCalc route", comparison];
comparison
