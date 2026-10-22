(* ::Title:: *)
(*Partonic Born cross sections: every massless 2 -> 2 QCD channel in D dimensions*)


(* ::Text:: *)
(*The lowest-order hard cross sections of all massless two-parton scattering channels with a tagged (observed) final-state parton, in D = 4 - 2 eps dimensions, for the spin average (UU) and for the double-longitudinal helicity difference of the incoming partons (LL).  They are used as the Born terms and as the lower-order ingredients of the collinear counterterms of the NLO proton-proton classes (Processes/Hadronic), and as the resolved-photon Born terms of the lepton-proton cross sections.  (The transverse spin transfers needed for the Collins effect are computed in Processes/Collins.)*)
(**)
(*Kinematics: a(pa) + b(pb) -> c(p) + d(pd), s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2, s + t + u = 0.  Output: Results/Born.wl with, for every channel "a b -> c", the squared matrix element averaged over the spins and colours of a and b and summed over those of c and d, divided by gs^4, as a function of s, t, D, CF, CA (TF = 1/2 is FeynCalc's normalization; u = -s-t).  The hard cross section s E dsigma/d^3p of the tagged parton is |M|^2/(16 Pi^2 s v) delta(1-w) with v = 1 + t/s (docs/CONVENTIONS.md).*)
(**)
(*Flavours: q is a quark of flavour i, qp a quark of a different flavour k, qb and qpb their antiquarks, g the gluon; FeynArts' SMQCD model represents q and qp by the up quark (generation 1) and the charm quark (generation 2); the photon, weak bosons, leptons and down-type quarks are excluded.  Identical particles: when the tagged parton has an identical unobserved partner (q q -> q q) the identical-particle factor 1/2 is cancelled by the two ways of tagging, so no factor is applied; this is the convention of the reference codes.*)
(**)
(*Spin densities (HVBM scheme; all four external momenta are four-dimensional): UU averages with 1/2 per incoming quark and 1/(D-2) per incoming gluon; LL inserts p-slash gamma5/2 for an incoming quark (-p-slash gamma5/2 for an incoming antiquark, whose positive-helicity spinor sum is (1 - gamma5) p-slash/2) and the circular density i eps^{mu nu l r}/(2 l.r) for an incoming gluon l (helicity + minus helicity -, including the factor 1/2; r is the other incoming momentum), so that the LL result is (sigma(++) - sigma(+-))/2 for the helicities of a and b, the other spins summed.  $LeviCivitaSign = -1 (Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = -4 i eps^{abcd}), as in the lepton processes.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[]


(* ::Section:: *)
(*Channels and diagrams*)


(* ::Text:: *)
(*The particle names map to FeynArts fields; a channel is {a, b, c} with the tagged parton c, the unobserved parton d follows from flavour conservation.  All channels needed by the collinear counterterms of the 16 proton-proton classes and by the resolved-photon terms are included; incoming-leg permutations are separate channels because t is defined with the first incoming parton.*)


(* ::Input:: *)
field = <|"q" -> F[3, {1}], "qb" -> -F[3, {1}], "qp" -> F[3, {2}], "qpb" -> -F[3, {2}], "g" -> V[5]|>;
excluded = {S[_], U[_], V[1 | 2 | 3], F[1 | 2 | 4], F[3, {3}]};
channels = {
  {"q", "qp", "q"}, {"q", "qp", "qp"}, {"qp", "q", "q"}, {"q", "qpb", "q"}, {"q", "qpb", "qpb"}, {"qpb", "q", "q"},
  {"q", "q", "q"}, {"q", "qb", "q"}, {"q", "qb", "qb"}, {"qb", "q", "q"}, {"q", "qb", "qp"}, {"q", "qb", "qpb"}, {"q", "qb", "g"},
  {"q", "g", "q"}, {"g", "q", "q"}, {"q", "g", "g"}, {"g", "q", "g"}, {"qb", "g", "qb"}, {"g", "qb", "qb"}, {"qb", "g", "g"},
  {"g", "g", "g"}, {"g", "g", "q"}, {"g", "g", "qb"},
  {"qb", "qpb", "qb"}, {"qb", "qp", "qb"}, {"qb", "qb", "qb"}, {"qb", "q", "qb"}, {"qb", "q", "g"}, {"qb", "q", "qp"}};
name[{a_, b_, c_}] := a <> " " <> b <> " -> " <> c;
flavourNumber[x_] := Switch[x, "q", {1, 0}, "qb", {-1, 0}, "qp", {0, 1}, "qpb", {0, -1}, "g", {0, 0}];
partner[{a_, b_, c_}] := Module[{ok = Select[Keys[field], Total[flavourNumber /@ {a, b}] === Total[flavourNumber /@ {c, #}] &]},
  Require[Length[ok] == 1, "The unobserved parton of " <> name[{a, b, c}] <> " is not unique", ok]; First[ok]];
diagrams = Association @ Table[name[ch] -> InsertFields[CreateTopologies[0, 2 -> 2],
  {field[ch[[1]]], field[ch[[2]]]} -> {field[ch[[3]]], field[partner[ch]]},
  InsertionLevel -> {Classes}, Model -> "SMQCD", ExcludeParticles -> excluded], {ch, channels}];
Map[Length, diagrams]


(* ::Input:: *)
Paint[diagrams["q g -> q"], ColumnsXRows -> {3, 1}, Numbering -> Simple, SheetHeader -> None, ImageSize -> {600, 200}]


(* ::Section:: *)
(*Squared matrix elements*)


(* ::Text:: *)
(*FCFAConvert gives the amplitude of each channel; FeynCalc forms the square: the quark spinors are closed by FermionSpinSum (with the densities of the spin case inserted for the incoming quarks), external gluons are summed with the physical polarization sum -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r) with a reference vector r (another external momentum; gauge invariance makes the result independent of it), or with the circular density for a polarized incoming gluon; the colour sums are done by SUNSimplify.*)


(* ::Input:: *)
FCClearScalarProducts[];
SetMandelstam[s, t, u, pa, pb, -p, -pd, 0, 0, 0, 0];
Scan[(Momentum[#, D - 4] = 0) &, {pa, pb, p, pd}];
FCSetDiracGammaScheme["BMHV"];
$LeviCivitaSign = -1;
amplitude[ch_] := Total[FCFAConvert[CreateFeynAmp[diagrams[name[ch]], Truncated -> False, PreFactor -> 1],
  IncomingMomenta -> {pa, pb}, OutgoingMomenta -> {p, pd}, UndoChiralSplittings -> True, ChangeDimension -> D,
  List -> True, SMP -> False, Contract -> False, DropSumOver -> True,
  FinalSubstitutions -> {FCGV["MU"] -> 0, FCGV["MC"] -> 0}] /. SMP["g_s"] -> gs /. GaugeXi[_] -> 1];     (* massless quarks, Feynman gauge *)
reference = <|pa -> pb, pb -> pa, p -> pa, pd -> pa|>;         (* reference vectors of the physical polarization sums *)
polarizationSum[l_, mu_, nu_] := With[{r = reference[l]}, -MTD[mu, nu] + (FVD[l, mu] FVD[r, nu] + FVD[r, mu] FVD[l, nu])/SPD[l, r]];
circular[l_, mu_, nu_] := Contract[I LC[mu, nu, al, be] FV[l, al] FV[reference[l], be]/(2 SP[l, reference[l]])];    (* four-dimensional; contracted here so that the dummy indices do not repeat *)
(* spin and colour averages of the incoming partons; a polarized gluon (LL) carries its 1/2 in the circular density *)
average[ch_, spinCase_] := Times @@ (Switch[#, "g", If[spinCase === "LL", 1, 1/(D - 2)]/(Nc^2 - 1), _, 1/(2 Nc)] & /@ ch[[1 ;; 2]]);
colourFactors[expr_] := Factor[expr /. SUNN -> Nc] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;

(* density matrices: which incoming (and tagged) momenta get gamma5 insertions *)
tag[l_] := Symbol[ToString[l] <> "Spin"];
square[ch_, spinCase_] := Module[{amp = amplitude[ch], conj, expr, gluons, indices, polarizedQuarks, polarizedGluons, taggedQuark},
  polarizedQuarks = If[spinCase === "LL", {pa, pb}[[Position[ch[[1 ;; 2]], Except["g"], {1}, Heads -> False][[All, 1]]]], {}];
  polarizedGluons = If[spinCase === "LL", {pa, pb}[[Position[ch[[1 ;; 2]], "g", {1}, Heads -> False][[All, 1]]]], {}];
  gluons = Union[Cases[amp, Polarization[l_, __] :> l, Infinity]];
  indices = Association @ Table[gluons[[i]] -> {Symbol["mu" <> ToString[i]], Symbol["nu" <> ToString[i]]}, {i, Length[gluons]}];
  conj = ComplexConjugate[amp];
  expr = amp conj;
  (* rename the spinors of polarized legs so that the slash FermionSpinSum produces can be replaced by the density *)
  expr = expr /. (With[{l = #, lTag = tag[#]}, Spinor[c_. Momentum[l, D], m_, o_] :> Spinor[c Momentum[lTag, D], m, o]] & /@ polarizedQuarks);
  expr = expr /. Flatten[Table[{Momentum[Polarization[l, I, ___], D] -> LorentzIndex[indices[l][[1]], D],
    Momentum[Polarization[l, -I, ___], D] -> LorentzIndex[indices[l][[2]], D]}, {l, gluons}]];
  expr = expr Times @@ Table[If[MemberQ[polarizedGluons, l], circular[l, indices[l][[1]], indices[l][[2]]], polarizationSum[l, indices[l][[1]], indices[l][[2]]]], {l, gluons}];
  expr = FermionSpinSum[FeynAmpDenominatorExplicit[expr]];
  (* p-slash gamma5 for a quark, -p-slash gamma5 for an antiquark (v vbar = (1 - lambda gamma5) p-slash/2); the 1/2 is in the average *)
  expr = expr /. (With[{l = #, lTag = tag[#], sign = If[StringEndsQ[ch[[Position[{pa, pb}, #][[1, 1]]]], "b"], -1, 1]},
    DiracGamma[Momentum[lTag, D], D] :> sign GSD[l].GA[5]] & /@ polarizedQuarks);
  expr = DiracSimplify[expr, DiracTraceEvaluate -> True] // Contract;
  expr = SUNSimplify[expr, Explicit -> True, SUNNToCACF -> False] // Contract // ExpandScalarProduct;
  colourFactors[Together[average[ch, spinCase] expr/gs^4 /. u -> -s - t]]];


(* ::Section:: *)
(*Amplitude records for the FORM and Tracer routes*)


(* ::Text:: *)
(*The same FeynArts amplitudes, written as plain records (one coefficient and the open fermion chains per term; vocabulary in FeynFacet/Amplitudes.wl), so that 2-Algebra-FORM.wl and 2-Algebra-Tracer.wl square exactly the amplitudes squared here, with their own Dirac, Lorentz and colour algebra.  The four-gluon vertex carries a FeynCalc dummy colour index ($AL$nnn, a different number in every session); it is renamed Glu6 (FeynArts numbers the external gluon colours Glu1..Glu4 and the propagator Glu5), so that the file is reproducible.*)


(* ::Input:: *)
fourGluonDummy[x_Symbol] := StringStartsQ[SymbolName[x], "$AL"];
records = Association @ Table[name[ch] -> (AmplitudeRecords[amplitude[ch]] /. x_Symbol /; fourGluonDummy[x] :> Glu6), {ch, channels}];
Require[FreeQ[records, x_Symbol /; StringStartsQ[SymbolName[x], "$"]], "A FeynCalc dummy symbol is left in the records"];
SaveAmplitudes[<|"Channels" -> Map[name, channels], "Records" -> records|>];
ShowRecords[records["q g -> q"]]


(* ::Section:: *)
(*All channels, both spin cases*)


(* ::Input:: *)
born = Association @ Table[spinCase -> Association @ Table[name[ch] -> square[ch, spinCase], {ch, channels}], {spinCase, {"UU", "LL"}}];
TableForm[Normal[Map[Simplify[# /. D -> 4] &, born["UU"]]]]


(* ::Text:: *)
(*Checks in four dimensions against the textbook table (Ellis, Stirling and Webber, QCD and Collider Physics, Table 7.1, at Nc = 3; the polarized quark channels follow from helicity conservation, e.g. q qb -> qp qpb only through opposite helicities), written with s + t + u = 0.  The identical-quark channel q q -> q keeps no factor 1/2 (two ways of tagging), so it is twice the table entry.  The Born terms of all ten Born classes, including the polarized gluon channels, are also compared with the reference codes in Processes/Hadronic.*)


(* ::Input:: *)
textbook = <|
  "q qp -> q" -> CF/CA (s^2 + u^2)/t^2,
  "q qpb -> q" -> CF/CA (s^2 + u^2)/t^2,
  "q q -> q" -> 2 (CF/CA ((s^2 + u^2)/t^2 + (s^2 + t^2)/u^2) - 2 CF/CA^2 s^2/(t u))/2,
  "q qb -> q" -> CF/CA ((s^2 + u^2)/t^2 + (u^2 + t^2)/s^2) - 2 CF/CA^2 u^2/(s t),
  "q qb -> qp" -> CF/CA (t^2 + u^2)/s^2,
  "q qb -> g" -> 2 CF^2/CA (t^2 + u^2)/(t u) - 2 CF (t^2 + u^2)/s^2,
  "q g -> q" -> -CF/CA (s^2 + u^2)/(s u) + (s^2 + u^2)/t^2,
  "g g -> g" -> 9/2 (3 - t u/s^2 - s u/t^2 - s t/u^2),
  "g g -> q" -> 1/(2 CA) (t^2 + u^2)/(t u) - CA/(CA^2 - 1) (t^2 + u^2)/s^2|>;
textbookLL = <|
  "q qp -> q" -> CF/CA (s^2 - u^2)/t^2,
  "q qb -> qp" -> -CF/CA (t^2 + u^2)/s^2,
  "q g -> q" -> -CF/CA (s^2 - u^2)/(s u) + (s^2 - u^2)/t^2,
  "g g -> q" -> -(1/(2 CA) (t^2 + u^2)/(t u) - CA/(CA^2 - 1) (t^2 + u^2)/s^2)|>;
check[spinCase_, table_] := Association @ KeyValueMap[
  #1 -> Simplify[(born[spinCase][#1] /. D -> 4) - #2 /. u -> -s - t /. {CF -> 4/3, CA -> 3}] &, table];
Require[And @@ (# === 0 & /@ Values[check["UU", textbook]]) && And @@ (# === 0 & /@ Values[check["LL", textbookLL]]),
  "A Born square differs from the textbook table", {check["UU", textbook], check["LL", textbookLL]}];
{check["UU", textbook], check["LL", textbookLL]}


(* ::Text:: *)
(*Crossing: the incoming-leg permutation b a -> c equals a b -> c with t <-> u; the LL difference is symmetric under the exchange as well.*)


(* ::Input:: *)
swap[expr_] := expr /. t -> -s - t;         (* t <-> u with u = -s-t *)
crossingPairs = {{"q qp -> q", "qp q -> q"}, {"q g -> q", "g q -> q"}, {"q g -> g", "g q -> g"}, {"q qb -> q", "qb q -> q"}};
Require[And @@ Flatten[Table[Simplify[born[sc][pair[[2]]] - swap[born[sc][pair[[1]]]]] === 0, {sc, {"UU", "LL"}}, {pair, crossingPairs}]],
  "Crossing relation between incoming-leg permutations fails"];
"crossing relations hold"


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
result = <|"Born" -> born, "Channels" -> Map[name, channels], "Variables" -> "s, t (u = -s-t), D, CF, CA; divided by gs^4"|>;
SaveResult["Born", result];
Keys[born["UU"]]
