(* ::Title:: *)
(*Transverse spin transfer at Born level in parton-parton scattering: the resolved-photon channels*)


(* ::Text:: *)
(*The lowest-order transverse-spin-transfer squares of the massless 2 -> 2 QCD channels a(pa) b(pb) -> c(p) d(pd) with the first incoming parton transversely polarized and the transverse spin of the tagged parton c measured, needed as the resolved-photon Born terms of the Collins asymmetry (Eq. (resolved-transverse-born) of the paper).  As in Processes/PartonicBorn: s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2, s + t + u = 0; the result is divided by gs^4; the spins and colours of b are averaged, those of d summed; the same FeynArts generation and polarization sums.  The transverse projections are those of 2-Algebra-FeynCalc.wl with pa in place of the polarized quark k and pb in place of the lepton l: T = sigma(a up -> c up) - sigma(a up -> c down), density pa-slash gamma5 a-slash/2 for a (the same expression for an incoming antiquark, v(s) vbar(s) = p-slash (1 + gamma5 s-slash)/2) and the projection p-slash gamma5 b-slash for c (quark or antiquark), with the normal axis (T_NN) and the in-plane axes (T_XX), and the Collins projection (T_NN + T_XX)/2.  A tagged gluon carries no transversity; transfer to a quark of another flavour must vanish (control channels: q qp -> qp, q qb -> qp, q qb -> qpb and their charge conjugates, so that the table holds every channel together with its charge conjugate).*)
(**)
(*Expected (the paper, with K = (Nc^2-1)/Nc^2, in four dimensions): q q' -> q and q qb' -> q: -K s u/t^2;  q q -> q: -K s u/t^2 + K s/(Nc t);  q qb -> q: -K s u/t^2 + K u/(Nc t);  q qb -> qb: -K/Nc;  q g -> q: -2 s u/t^2 + K; the charge-conjugate channels have the same values.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[]


(* ::Section:: *)
(*Channels and diagrams*)


(* ::Input:: *)
field = <|"q" -> F[3, {1}], "qb" -> -F[3, {1}], "qp" -> F[3, {2}], "qpb" -> -F[3, {2}], "g" -> V[5]|>;
excluded = {S[_], U[_], V[1 | 2 | 3], F[1 | 2 | 4], F[3, {3}]};
representatives = {{"q", "qp", "q"}, {"q", "qpb", "q"}, {"q", "q", "q"}, {"q", "qb", "q"}, {"q", "qb", "qb"}, {"q", "g", "q"}};
controls = {{"q", "qp", "qp"}, {"q", "qb", "qp"}, {"q", "qb", "qpb"}};
conjugate[x_] := x /. {"q" -> "qb", "qb" -> "q", "qp" -> "qpb", "qpb" -> "qp"};
channels = Join[representatives, controls, conjugate /@ representatives, conjugate /@ controls];
name[{a_, b_, c_}] := a <> " " <> b <> " -> " <> c;
flavourNumber[x_] := Switch[x, "q", {1, 0}, "qb", {-1, 0}, "qp", {0, 1}, "qpb", {0, -1}, "g", {0, 0}];
partner[{a_, b_, c_}] := Module[{ok = Select[Keys[field], Total[flavourNumber /@ {a, b}] === Total[flavourNumber /@ {c, #}] &]},
  Require[Length[ok] == 1, "The unobserved parton of " <> name[{a, b, c}] <> " is not unique", ok]; First[ok]];
diagrams = Association @ Table[name[ch] -> InsertFields[CreateTopologies[0, 2 -> 2],
  {field[ch[[1]]], field[ch[[2]]]} -> {field[ch[[3]]], field[partner[ch]]},
  InsertionLevel -> {Classes}, Model -> "SMQCD", ExcludeParticles -> excluded], {ch, channels}];
Map[Length, diagrams]


(* ::Section:: *)
(*Squared matrix elements with the transverse projections*)


(* ::Text:: *)
(*As in Processes/PartonicBorn: the spinors of the polarized legs are renamed so that the slashes FermionSpinSum produces can be replaced by the densities; external gluons get the physical polarization sum with a reference vector; colour sums by SUNSimplify.*)


(* ::Input:: *)
FCClearScalarProducts[];
SetMandelstam[s, t, u, pa, pb, -p, -pd, 0, 0, 0, 0];
Scan[(Momentum[#, D - 4] = 0) &, {pa, pb, p, pd, nn}];
SPD[nn] = -1; SP[nn] = -1; Scan[(SPD[nn, #] = 0; SP[nn, #] = 0) &, {pa, pb, p, pd}];
FCSetDiracGammaScheme["BMHV"];
$LeviCivitaSign = -1;
(* the spin axes as vectors of their own with explicit scalar products, as in 2-Algebra-FeynCalc.wl: aX = (t/s) pb + (u/s) pa + p, bX = t/(t+u) pb - u/(t+u) pa - s(u-t)/(t+u)^2 p *)
Scan[(Momentum[#, D - 4] = 0) &, {aX, bX}];
axisCoefficients = <|aX -> {t/s, u/s, 1}, bX -> {t/(t + u), -u/(t + u), -s (u - t)/(t + u)^2}|>;
Do[With[{c = axisCoefficients[axis], base = {pb, pa, p}},
  Do[SP[axis, x] = SPD[axis, x] = Factor[Sum[c[[i]] SP[base[[i]], x], {i, 3}]], {x, {pa, pb, p, pd}}];
  SP[axis, nn] = SPD[axis, nn] = 0], {axis, {aX, bX}}];
Do[With[{c1 = axisCoefficients[pair[[1]]], c2 = axisCoefficients[pair[[2]]], base = {pb, pa, p}},
  SP[pair[[1]], pair[[2]]] = SPD[pair[[1]], pair[[2]]] = Factor[Sum[c1[[i]] c2[[j]] SP[base[[i]], base[[j]]], {i, 3}, {j, 3}]]], {pair, {{aX, aX}, {bX, bX}, {aX, bX}}}];
axes = <|"aX" -> aX, "bX" -> bX, "Lengths" -> -t u/(t + u)|>;
Require[Simplify[{SP[aX, pb], SP[aX, pa], SP[bX, p], SP[bX, pa] + SP[bX, pb], SP[aX] + t u/s, SP[bX] + s t u/(t + u)^2} /. u -> -s - t] === {0, 0, 0, 0, 0, 0},
  "The in-plane axes are not what the text says"];
amplitude[ch_] := Total[FCFAConvert[CreateFeynAmp[diagrams[name[ch]], Truncated -> False, PreFactor -> 1],
  IncomingMomenta -> {pa, pb}, OutgoingMomenta -> {p, pd}, UndoChiralSplittings -> True, ChangeDimension -> D,
  List -> True, SMP -> False, Contract -> False, DropSumOver -> True,
  FinalSubstitutions -> {FCGV["MU"] -> 0, FCGV["MC"] -> 0}] /. SMP["g_s"] -> gs /. GaugeXi[_] -> 1];
reference = <|pa -> pb, pb -> pa, p -> pa, pd -> pa|>;
polarizationSum[l_, mu_, nu_] := With[{r = reference[l]}, -MTD[mu, nu] + (FVD[l, mu] FVD[r, nu] + FVD[r, mu] FVD[l, nu])/SPD[l, r]];
average[ch_] := 1/2 Switch[ch[[2]], "g", 1/(D - 2)/(Nc^2 - 1), _, 1/(2 Nc)] Switch[ch[[1]], "g", 1/(Nc^2 - 1), _, 1/Nc];   (* 1/2: the fixed spin of a *)
colourFactors[expr_] := Factor[expr /. SUNN -> Nc];
tag[l_] := Symbol[ToString[l] <> "Spin"];
spinVectors = <|"NN" -> {nn, nn}, "XX" -> {axes["aX"], axes["bX"]}|>;
square[ch_, spinCase_] := Module[{amp = amplitude[ch], conj, expr, gluons, indices, sa, sc},
  {sa, sc} = spinVectors[spinCase];
  gluons = Union[Cases[amp, Polarization[l_, __] :> l, Infinity]];
  indices = Association @ Table[gluons[[i]] -> {Symbol["mu" <> ToString[i]], Symbol["nu" <> ToString[i]]}, {i, Length[gluons]}];
  conj = ComplexConjugate[amp];
  expr = amp conj;
  expr = expr /. (With[{l = #, lTag = tag[#]}, Spinor[c_. Momentum[l, D], m_, o_] :> Spinor[c Momentum[lTag, D], m, o]] & /@ {pa, p});
  expr = expr /. Flatten[Table[{Momentum[Polarization[l, I, ___], D] -> LorentzIndex[indices[l][[1]], D],
    Momentum[Polarization[l, -I, ___], D] -> LorentzIndex[indices[l][[2]], D]}, {l, gluons}]];
  expr = expr Times @@ Table[polarizationSum[l, indices[l][[1]], indices[l][[2]]], {l, gluons}];
  expr = FermionSpinSum[FeynAmpDenominatorExplicit[expr]];
  expr = expr /. {DiracGamma[Momentum[tag[pa], D], D] :> GS[pa].GA[5].GS[sa], DiracGamma[Momentum[tag[p], D], D] :> GS[p].GA[5].GS[sc]};
  expr = DiracSimplify[expr, DiracTraceEvaluate -> True] // Contract;
  expr = SUNSimplify[expr, Explicit -> True, SUNNToCACF -> False] // Contract // ExpandScalarProduct;
  colourFactors[Together[average[ch] expr/gs^4 /. u -> -s - t]]/If[spinCase === "XX", axes["Lengths"], 1]];


(* ::Input:: *)
transfer = Association @ Table[name[ch] -> Module[{nn4 = square[ch, "NN"], xx4 = square[ch, "XX"]},
  <|"NN" -> nn4, "XX" -> xx4, "Collins" -> Factor[(nn4 + xx4)/2]|>], {ch, channels}];
TableForm[Normal[Map[Simplify[#["Collins"] /. D -> 4] &, transfer]]]


(* ::Section:: *)
(*Checks against the paper*)


(* ::Text:: *)
(*In four dimensions the two components agree in every channel; the Collins projection is compared with the table of the paper, the control channels and their charge conjugates must vanish, and the charge-conjugate channels must equal their representatives.*)


(* ::Input:: *)
K = (Nc^2 - 1)/Nc^2;
expected = <|"q qp -> q" -> -K s u/t^2, "q qpb -> q" -> -K s u/t^2, "q q -> q" -> -K s u/t^2 + K s/(Nc t),
  "q qb -> q" -> -K s u/t^2 + K u/(Nc t), "q qb -> qb" -> -K/Nc, "q g -> q" -> -2 s u/t^2 + K|>;
check[expr_] := Simplify[expr /. D -> 4 /. u -> -s - t];
Require[And @@ Table[check[transfer[ch]["NN"] - transfer[ch]["XX"]] === 0, {ch, Keys[transfer]}], "The two transverse components differ in four dimensions"];
differences = Association @ KeyValueMap[#1 -> check[transfer[#1]["Collins"] - #2] &, expected];
Require[And @@ (# === 0 & /@ Values[differences]), "A resolved transverse Born term differs from the paper", differences];
Require[And @@ Table[check[transfer[name[ch]]["Collins"]] === 0, {ch, Join[controls, conjugate /@ controls]}], "A control channel does not vanish"];
Require[And @@ Table[check[transfer[name[ch]]["Collins"] - transfer[name[conjugate[ch]]]["Collins"]] === 0, {ch, representatives}],
  "A charge-conjugate channel differs from its representative"];
differences


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
result = <|"Transfer" -> transfer, "Channels" -> Map[name, channels], "Variables" -> "s, t (u = -s-t), D, Nc; divided by gs^4; first incoming parton polarized, tagged parton c"|>;
SaveResult["ResolvedBorn", result];
Keys[transfer]
