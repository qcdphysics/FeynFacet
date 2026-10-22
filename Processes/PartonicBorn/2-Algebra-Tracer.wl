(* ::Title:: *)
(*Partonic Born cross sections, Tracer route: the same table with Tracer and an explicit colour reduction*)


(* ::Text:: *)
(*Same task and same output format as 1-Born.wl (the FeynCalc route), with the classic Tracer package for the Dirac and Lorentz algebra and the colour sums written out below.  Input: the amplitude records of all 29 channels written by 1-Born.wl (Results/Amplitudes.wl), and its table Results/Born.wl for the comparison.  Output: Results/Tracer/Born.wl.  Run this file in a fresh kernel: Tracer and FeynCalc cannot be loaded together.*)
(**)
(*Tracer's vocabulary is summarized in FeynFacet/Tracer.wl and explained in Processes/LeptonQuark/2-Algebra-Tracer.wl: G[line, p, {mu}, ...] is the string p-slash gamma^mu ... on a fermion line, traced when Spur[line] is set; S[p, q] = p.q, S[p, {mu}] = p^mu, S[{mu}, {nu}] = g^{mu nu}, and repeated indices in an expanded product are contracted; H is the (D-4)-dimensional part of a scalar product; Eps is the four-dimensional epsilon tensor; VectorDimension[D] and AntiCommute[off] select D dimensions and the HVBM gamma5.  Tracer's gamma5 has the opposite sign to FeynCalc's with $LeviCivitaSign = -1 (Tr[gamma5 gamma^a gamma^b gamma^c gamma^d] = +4 i Eps[a,b,c,d] in Tracer), so the circular gluon density below carries the opposite sign to the FeynCalc route: every LL term contains exactly two such objects (two gamma5's, a gamma5 and a circular density, or two circular densities), so the table is the same.*)
(**)
(*Kinematics as in 1-Born.wl: a(pa) + b(pb) -> c(p) + d(pd), s = 2 pa.pb, t = -2 pa.p, u = -2 pb.p; all four momenta are four-dimensional.*)


(* ::Section:: *)
(*Load*)


(* ::Text:: *)
(*Load the FeynFacet package (amplitude records, the conversion of fermion chains to Tracer strings, and the checks Require).*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*Load Tracer (its location is set in Paths.wl).*)


(* ::Input:: *)
LoadTracer[]


(* ::Text:: *)
(*The amplitude records of 1-Born.wl, and Tracer's settings: D-dimensional gamma matrices and vectors, the HVBM gamma5 (AntiCommute[off]), and fermion lines 1 and 2 traced (an interference has at most two closed quark lines).  The output shows the records of q qb -> q: every record is one term of the amplitude, a coefficient (colour factor times Lorentz factor) and its open fermion chains.*)


(* ::Input:: *)
amplitudes = LoadAmplitudes[];
channels = amplitudes["Channels"];
records = amplitudes["Records"];
VectorDimension[D];
AntiCommute[off];
Spur[1]; Spur[2];
ShowRecords[records["q qb -> q"]]


(* ::Text:: *)
(*Tracer's gamma5 convention, checked on a separate traced line 3: Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = +4 i Eps[a,b,c,d], with the trace written as Tracer's G[3, {a}, {b}, {c}, {d}, G5].  This is minus the FeynCalc convention ($LeviCivitaSign = -1, -4 i eps^{abcd}) used in 1-Born.wl, the reason for the sign of the circular gluon density below.*)


(* ::Input:: *)
Spur[3];
tracerG5 = Expand[G[3, {mu1}, {mu2}, {mu3}, {mu4}, G5]];
Require[Expand[tracerG5 - 4 I Eps[{mu1}, {mu2}, {mu3}, {mu4}]] === 0, "Tracer's gamma5 trace is not +4 i Eps", tracerG5];
tracerG5


(* ::Section:: *)
(*Colour: an explicit SU(N) reduction*)


(* ::Text:: *)
(*Every record coefficient is a colour factor (a product of ColorT, ColorF, ColorDelta, GluonDelta) times a Lorentz factor.  The colour sum of an interference, amplitude colour C1 times the complex conjugate of C2 summed over the colours of all four partons, is a number, computed here from three identities with the normalization Tr(T^a T^b) = delta^{ab}/2 of the FeynCalc route:*)
(**)
(*  (1) structure constants as traces:  f^{abc} = -2 i ( Tr[T^a T^b T^c] - Tr[T^b T^a T^c] )   (from [T^a, T^b] = i f^{abc} T^c);*)
(*  (2) the Fierz identity, which removes every summed gluon index:  (T^a)_{ij} (T^a)_{kl} = 1/2 ( delta_{il} delta_{kj} - delta_{ij} delta_{kl}/Nc );*)
(*  (3) closed quark-colour loops:  delta_{ii} = Nc, and delta^{aa} = Nc^2 - 1 for a gluon delta.*)
(**)
(*In this notation a generator is gen[a, i, j] = (T^a)_{ij}, a quark delta deltaF[i, j], a gluon delta deltaA[a, b]; products of generators and traces are written with explicit fresh quark-colour indices.  The complex conjugate transposes the generators, (T^a)_{ij}^* = (T^a)_{ji} (f and the deltas are real), and renames the colour indices summed inside the amplitude, so that they are not confused with those of the other factor.*)


(* ::Input:: *)
(* elementary objects with explicit indices; fresh[] makes a new summed quark-colour index *)
fresh[] := Unique["x"];
generatorString[{a_}, i_, j_] := gen[a, i, j];
generatorString[{a_, rest__}, i_, j_] := With[{x = fresh[]}, gen[a, i, x] generatorString[{rest}, x, j]];
colourTrace[list_List] := With[{x = fresh[]}, generatorString[list, x, x]];
elementaryFactor[ColorT[list_List, i_, j_]] := generatorString[list, i, j];
elementaryFactor[ColorDelta[i_, j_]] := deltaF[i, j];
elementaryFactor[GluonDelta[a_, b_]] := deltaA[a, b];
elementaryFactor[ColorF[a_, b_, c_]] := -2 I (colourTrace[{a, b, c}] - colourTrace[{b, a, c}]);                (* identity (1) *)
elementaryFactor[x_?NumericQ] := x;
(* factor by factor, a power such as f^{abc} f^{abc} = ColorF[a,b,c]^2 counted as separate copies, so that every copy gets its own fresh indices *)
factorList[colour_] := Flatten[Replace[If[Head[colour] === Times, List @@ colour, {colour}], Power[x_, n_Integer?Positive] :> ConstantArray[x, n], {1}]];
elementary[colour_] := Times @@ (elementaryFactor /@ factorList[colour]);

(* reduction of one monomial: gluon deltas, then Fierz, then quark deltas *)
reduceRules = {
  deltaA[a_, a_] :> Nc^2 - 1,                                                                                   (* (3) *)
  deltaA[a_, b_] gen[b_, i_, j_] :> gen[a, i, j], deltaA[b_, a_] gen[b_, i_, j_] :> gen[a, i, j],
  deltaA[a_, b_] deltaA[b_, c_] :> deltaA[a, c], deltaA[a_, b_] deltaA[c_, b_] :> deltaA[a, c], deltaA[b_, a_] deltaA[b_, c_] :> deltaA[a, c],
  gen[a_, i_, j_] gen[a_, k_, l_] :> 1/2 (deltaF[i, l] deltaF[k, j] - deltaF[i, j] deltaF[k, l]/Nc),                               (* (2) *)
  gen[a_, i_, j_]^2 :> 1/2 (deltaF[i, j]^2 - deltaF[i, j]^2/Nc),                                          (* (2) with k = i, l = j *)
  deltaF[i_, i_] :> Nc,                                                                                         (* (3) *)
  deltaF[i_, j_]^2 :> Nc,
  deltaF[i_, j_] deltaF[j_, k_] :> deltaF[i, k], deltaF[i_, j_] deltaF[k_, j_] :> deltaF[i, k], deltaF[j_, i_] deltaF[j_, k_] :> deltaF[i, k],
  deltaF[i_, j_] gen[a_, j_, k_] :> gen[a, i, k], deltaF[j_, i_] gen[a_, j_, k_] :> gen[a, i, k],
  deltaF[i_, j_] gen[a_, k_, j_] :> gen[a, k, i], deltaF[j_, i_] gen[a_, k_, j_] :> gen[a, k, i]};
reduce[expr_] := FixedPoint[Expand[# //. reduceRules] &, Expand[expr]];

(* internal (summed) colour indices of a colour factor: those that occur twice *)
colourIndices[colour_] := Join[Cases[colour, ColorT[l_List, i_, j_] :> Sequence @@ Join[l, {i, j}], {0, Infinity}],
  Cases[colour, (ColorDelta | GluonDelta | ColorF)[x__] :> Sequence[x], {0, Infinity}]];
internalColour[colour_] := Select[Tally[colourIndices[colour]], #[[2]] > 1 &][[All, 1]];
conjugateColour[colour_] := colour /. Thread[internalColour[colour] -> (Symbol["c" <> SymbolName[#]] & /@ internalColour[colour])] /.
  ColorT[list_List, i_, j_] :> ColorT[Reverse[list], j, i];

(* canonical form: internal indices renamed in a fixed order, so that equal structures share one evaluation *)
canonical[colour_] := colour /. Thread[internalColour[colour] -> Table[Symbol["k" <> ToString[n]], {n, Length[internalColour[colour]]}]];
colourSum[c1_, c2_] := colourSumCanonical[canonical[c1], canonical[c2]];
colourSumCanonical[c1_, c2_] := colourSumCanonical[c1, c2] = Module[{value = reduce[elementary[c1 conjugateColour[c2]]]},
  Require[FreeQ[value, gen | deltaF | deltaA | elementaryFactor], "Colour reduction left unknown or uncontracted objects", value];
  Factor[value]];


(* ::Text:: *)
(*Checks of the reduction on known colour sums (CF = (Nc^2-1)/(2 Nc)): Tr(T^a T^a) = CF Nc; f^{abc} f^{abc} = Nc (Nc^2-1); Tr(T^a T^b T^b T^a) = CF^2 Nc; Tr(T^a T^b T^a T^b) = (CF - Nc/2) CF Nc = -CF/2; and f^{abc} Tr(T^c T^b T^a) = -(i/4) f^{abc} f^{abc} = -i Nc (Nc^2-1)/4, which tests the sign of identity (1).  In each case the second argument of colourSum is the factor that is complex conjugated; the output is the list of differences (computed minus known), all zero.*)


(* ::Input:: *)
colourTests = {
  colourSum[ColorT[{a1}, i1, i2], ColorT[{a1}, i1, i2]] - (Nc^2 - 1)/2,
  colourSum[ColorF[a1, a2, a3], ColorF[a1, a2, a3]] - Nc (Nc^2 - 1),
  colourSum[ColorT[{a1, a2}, i1, i2], ColorT[{a1, a2}, i1, i2]] - Nc ((Nc^2 - 1)/(2 Nc))^2,
  colourSum[ColorT[{a1, a2}, i1, i2], ColorT[{a2, a1}, i1, i2]] - (-(Nc^2 - 1)/(4 Nc)),
  colourSum[ColorF[a1, a2, b] ColorT[{b}, i1, i2], ColorT[{a1, a2}, i1, i2]] - (-I Nc (Nc^2 - 1)/4)};
Require[Simplify[colourTests] === {0, 0, 0, 0, 0}, "The colour reduction fails a known colour sum", colourTests];
Simplify[colourTests]


(* ::Section:: *)
(*From records to traces*)


(* ::Text:: *)
(*Records with the same colour factor and the same fermion chains are summed first (their Lorentz factors add), so that each interference is computed once per pair of such groups; the output is the number of groups per channel.  The complex conjugate of a Lorentz factor: i -> -i, conjugated polarization vectors (the third argument of PolVector, +1 for eps and -1 for eps*, flipped), and primed names (prefix c) for the Lorentz indices, so that the summed indices of the amplitude and of its conjugate stay independent.*)


(* ::Input:: *)
colourPart[coefficient_] := Times @@ Cases[If[Head[coefficient] === Times, List @@ coefficient, {coefficient}], _ColorT | _ColorF | _ColorDelta | _GluonDelta];
groups[recordList_] := KeyValueMap[<|"Colour" -> #1[[1]], "Chains" -> #1[[2]], "Lorentz" -> Total[#2]|> &,
  GroupBy[recordList, {colourPart[#["Coefficient"]], #["Chains"]} & -> (#["Coefficient"]/colourPart[#["Coefficient"]] &)]];
lorentzIndices[expr_] := Union[Cases[expr, i_Symbol /; StringStartsQ[SymbolName[i], "Lor"], Infinity]];
primed[expr_] := expr /. Thread[lorentzIndices[expr] -> (Symbol["c" <> SymbolName[#]] & /@ lorentzIndices[expr])];
conjugateLorentz[lorentz_] := primed[lorentz] /. {Complex[x_, y_] :> Complex[x, -y], PolVector[q_, i_, sign_] :> PolVector[q, i, -sign]};
Length /@ Association @ Table[channel -> groups[records[channel]], {channel, channels}]


(* ::Text:: *)
(*Closing the fermion lines.  Each chain is a segment with a left end (the barred spinor ubar(q) or vbar(q)), a gamma string, and a right end (the spinor u(q') or v(q')); both ends are labelled {"U" or "V", momentum}.  The complex conjugate of a chain swaps its two ends and reverses its gamma string (the Dirac adjoint of ubar(a) W u(b) is ubar(b) Wbar u(a)).  In the product of amplitude and conjugate, the segment that ends on u(q) is followed by the one that starts with ubar(q): the spin sum u(q) ubar(q) = q-slash (v(q) vbar(q) = q-slash) joins them.  Following the joins until the start is reached again closes a trace; the interference has one or two traces, depending on how the quark lines of the two diagrams connect.  At every junction a density marker rho[{type, q}] is left for the spin case to fill in.  The output is an example: the first record of q qb -> q times the conjugate of its second record (the t-channel scattering graph times the conjugate of the s-channel annihilation graph), which gives a single trace running through all four quark spin sums.*)


(* ::Input:: *)
spinorEnd[(Ubar | Uspinor)[q_]] := {"U", q};
spinorEnd[(Vbar | Vspinor)[q_]] := {"V", q};
segment[chain_Chain] := <|"Left" -> spinorEnd[First[chain]], "Word" -> TracerString[chain], "Right" -> spinorEnd[Last[chain]]|>;
conjugateSegment[chain_Chain] := <|"Left" -> spinorEnd[Last[chain]], "Word" -> TracerStringReversed[chain], "Right" -> spinorEnd[First[chain]]|>;
closeTraces[segments_List] := Module[{rest = segments, traces = {}, start, current, word, next},
  While[rest =!= {},
    current = First[rest]; rest = Rest[rest]; start = current["Left"]; word = current["Word"];
    While[current["Right"] =!= start,
      next = FirstPosition[rest, s_ /; s["Left"] === current["Right"], Missing[], {1}, Heads -> False];
      Require[!MissingQ[next], "A fermion line does not close", current];
      word = Join[word, {rho[current["Right"]]}, rest[[First[next]]]["Word"]];
      current = rest[[First[next]]]; rest = Delete[rest, next]];
    AppendTo[traces, Append[word, rho[current["Right"]]]]];
  traces];
closeTraces[Join[segment /@ records["q qb -> q"][[1]]["Chains"], conjugateSegment /@ records["q qb -> q"][[2]]["Chains"]]]


(* ::Text:: *)
(*Spin densities.  Outgoing quarks, and incoming quarks for UU: the spin sum q-slash.  For LL an incoming quark gets q-slash gamma5 (the factor 1/2 is in the average) and an incoming antiquark -q-slash gamma5, the densities of the FeynCalc route.  Gluons: the physical polarization sum -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r) with the reference vectors of 1-Born.wl, and for LL an incoming gluon gets the circular density -i Eps^{mu nu l r}/(2 l.r) (minus the FeynCalc form, see the introduction), mu being the index of the amplitude.*)


(* ::Input:: *)
reference = <|pa -> pb, pb -> pa, p -> pa, pd -> pa|>;
density[{type_, q_}, spinCase_] := If[spinCase === "LL" && MemberQ[{pa, pb}, q],
  {If[type === "V", -1, 1], {q, G5}},
  {1, {q}}];
traceProduct[traces_List, spinCase_] := Times @@ MapIndexed[Function[{word, line}, Module[{pieces},
    pieces = Replace[word, {rho[end_] :> density[end, spinCase], x_ :> {1, {x}}}, {1}];
    (Times @@ pieces[[All, 1]]) G @@ Join[line, Flatten[pieces[[All, 2]], 1]]]], traces];
polarizationSums[spinCase_] := Join[
  If[spinCase === "LL", Table[With[{l = l, r = reference[l]},
    PolVector[l, i_, 1] PolVector[l, j_, -1] :> -I Eps[{i}, {j}, l, r]/(2 S[l, r])], {l, {pa, pb}}], {}],
  Table[With[{l = l, r = reference[l]},
    PolVector[l, i_, _] PolVector[l, j_, _] :> -S[{i}, {j}] + (S[l, {i}] S[r, {j}] + S[r, {i}] S[l, {j}])/S[l, r]], {l, {pa, pb, p, pd}}]];
toTracer = {Den[q_] :> 1/S[q, q], Metric[i_, j_] :> S[{i}, {j}], MomentumComponent[q_, i_] :> S[q, {i}]};


(* ::Section:: *)
(*Squared matrix elements*)


(* ::Text:: *)
(*interference[g1, g2, spinCase] is the Lorentz part of group g1 times the conjugate of group g2, with the traces of the closed fermion lines and the polarization sums; the colour sum multiplies it.  The four momenta are four-dimensional and on shell: their hatted ((D-4)-dimensional) parts H vanish, and H is made linear in its arguments, because Tracer leaves H of a momentum sum such as H[p - pb, p - pb] (from the propagator momenta inside the HVBM traces) unexpanded; the products are expanded so that Tracer contracts the repeated indices.  The averages are those of 1-Born.wl: 1/(2 Nc) for an incoming quark, 1/((D-2)(Nc^2-1)) for an unpolarized incoming gluon and 1/(Nc^2-1) for a polarized one.  square[channel, spinCase] sums the interferences of all pairs of groups with their colour sums, applies the averages, divides by gs^4 and sets u = -s-t; Tracer's OnShell declares the scalar products (massless momenta, pa.pb = p.pd = s/2, pa.p = pb.pd = -t/2, pa.pd = pb.p = -u/2), and the definitions for H and S stay local to square (InheritedBlock).  The output is the simplest channel, q qp -> q (one t-channel gluon exchange), UU; at D = 4 it is CF/CA (s^2+u^2)/t^2 = (Nc^2-1)/(2 Nc^2) (s^2+u^2)/t^2.*)


(* ::Input:: *)
interference[g1_, g2_, spinCase_] := Module[{traces, expr},
  traces = closeTraces[Join[segment /@ g1["Chains"], conjugateSegment /@ primed[g2["Chains"]]]];
  expr = Expand[g1["Lorentz"] conjugateLorentz[g2["Lorentz"]]] //. polarizationSums[spinCase];      (* expanded: one pair of polarization vectors per gluon in every term *)
  Require[FreeQ[expr, PolVector], "A polarization vector was left", expr];
  Expand[(expr /. toTracer) traceProduct[traces, spinCase]]];
incoming[channel_String] := Take[StringSplit[channel], 2];
average[channel_, spinCase_] := Times @@ (Switch[#, "g", If[spinCase === "LL", 1, 1/(D - 2)]/(Nc^2 - 1), _, 1/(2 Nc)] & /@ incoming[channel]);
fourDimensional[momenta_List] := (Unprotect[H];
  H[x_Plus] := H /@ x; H[c_?NumericQ x_] := c H[x];
  H[x_, y_Plus] := H[x, #] & /@ y; H[x_Plus, y_] := H[#, y] & /@ x;
  H[x_, c_?NumericQ y_] := c H[x, y]; H[c_?NumericQ x_, y_] := c H[x, y];
  Scan[(H[#] = 0; H[#, _] = 0; H[_, #] = 0) &, momenta]; Protect[H]);
square[channel_, spinCase_] := Internal`InheritedBlock[{S, H},
  fourDimensional[{pa, pb, p, pd}];
  OnShell[on, {pa, 0}, {pb, 0}, {p, 0}, {pd, 0}, {pa, pb, s/2}, {p, pd, s/2}, {pa, p, -t/2}, {pb, pd, -t/2}, {pa, pd, -u/2}, {pb, p, -u/2}];
  Module[{channelGroups = groups[records[channel]]},
    Factor[Together[average[channel, spinCase] Sum[colourSum[g1["Colour"], g2["Colour"]] interference[g1, g2, spinCase],
      {g1, channelGroups}, {g2, channelGroups}]/gs^4 /. u -> -s - t]]]];
square["q qp -> q", "UU"]


(* ::Section:: *)
(*All channels, both spin cases*)


(* ::Text:: *)
(*The table of all 29 channels, UU and LL, as rational functions of s, t, D and Nc.  The output is the number of distinct colour sums that were evaluated (each pair of canonical colour structures once, remembered by colourSumCanonical).*)


(* ::Input:: *)
bornNc = Association @ Table[spinCase -> Association @ Table[channel -> square[channel, spinCase], {channel, channels}], {spinCase, {"UU", "LL"}}];
Length[DownValues[colourSumCanonical]] - 1


(* ::Text:: *)
(*The comparison: the FeynCalc table is written with CF and CA; with CF = (Nc^2-1)/(2 Nc) and CA = Nc both tables are rational functions of s, t, D and Nc, and their difference must vanish identically, channel by channel.*)


(* ::Input:: *)
feynCalc = LoadResult["Born"]["Born"];
toNc = {CF -> (Nc^2 - 1)/(2 Nc), CA -> Nc};
difference = Association @ Table[spinCase -> Association @ Table[channel ->
    Together[bornNc[spinCase][channel] - (feynCalc[spinCase][channel] /. toNc)], {channel, channels}], {spinCase, {"UU", "LL"}}];
Require[Union[Flatten[Values /@ Values[difference]]] === {0}, "The Tracer table differs from the FeynCalc table", difference];
TableForm[Table[{channel, difference["UU"][channel], difference["LL"][channel]}, {channel, channels}],
  TableHeadings -> {None, {"channel", "UU: Tracer - FeynCalc", "LL: Tracer - FeynCalc"}}]


(* ::Section:: *)
(*Save*)


(* ::Text:: *)
(*Written in the format of the FeynCalc route: CF and CA restored as in 1-Born.wl ((Nc-1)(Nc+1) = 2 Nc CF, then Nc = CA).  The saved table is checked once more against the FeynCalc one.*)


(* ::Input:: *)
colourFactors[expr_] := Factor[expr] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
born = Map[colourFactors, bornNc, {2}];
Require[Union[Flatten[Table[Together[(born[sc][ch] - feynCalc[sc][ch]) /. toNc], {sc, {"UU", "LL"}}, {ch, channels}]]] === {0},
  "The saved Tracer table differs from the FeynCalc table"];
result = <|"Born" -> born, "Channels" -> channels, "Variables" -> "s, t (u = -s-t), D, CF, CA; divided by gs^4"|>;
SaveResult["Born", result, "Tracer"];
Keys[born["UU"]]
