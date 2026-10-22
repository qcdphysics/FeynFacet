(* ::Title:: *)
(*Lepton-quark scattering with transverse spin transfer (Collins), step 1: diagrams and amplitudes*)


(* ::Text:: *)
(*Process: l(l) + q(k) -> q(p) + X with a transversely polarized incoming quark and the transverse spin of the observed quark measured (the transverse spin transfer that feeds the Collins asymmetry of Section 6 of the paper).  The diagrams are those of Processes/LeptonQuark, so this step is the same as there: it generates every Feynman diagram with FeynArts and converts the amplitudes into a plain, program-independent form that the three algebra routes read in step 2.  The spin projections are made in step 2.*)
(**)
(*Four sets of diagrams are needed:*)
(*  Born:     l q -> l q by photon exchange (one diagram);*)
(*  Real:     l q -> l q g, gluon emitted from the incoming or outgoing quark (two diagrams);*)
(*  Virtual:  the one-loop QCD corrections to l q -> l q (vertex and two external-leg graphs);*)
(*  Photon:   gamma q -> q g, the Born process needed to subtract the collinear photon emitted by the lepton (two diagrams).*)
(**)
(*Momenta: l, k incoming; lp (scattered lepton), p (observed quark), r (gluon) outgoing; loop is the loop momentum; the photon in the last process carries q.  Couplings: ee is the electromagnetic coupling, eq the quark charge in units of ee, gs the strong coupling.  FeynArts' SMQCD model has an up quark of charge 2/3; the prefactor 3 eq/2 replaces it by the general charge eq.*)


(* ::Section:: *)
(*Load the package*)


(* ::Text:: *)
(*The package must be loaded in a cell of its own, before any cell that uses its functions is parsed.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[]


(* ::Section:: *)
(*Generate the diagrams*)


(* ::Text:: *)
(*F[2,{1}] is the electron, F[3,{1}] the up quark, V[1] the photon and V[5] the gluon in FeynArts' SMQCD model.  We exclude Higgs and Goldstone bosons S, ghosts U, the weak bosons V[2|3], neutrinos F[1] and down-type quarks F[4]: only QED and QCD interactions remain.*)


(* ::Input:: *)
process = {F[2, {1}], F[3, {1}]} -> {F[2, {1}], F[3, {1}]};
excluded = {S[_], U[_], V[2 | 3], F[1 | 4]};
generate[topologies_, fields_, extra_: {}] := InsertFields[topologies, fields,
  InsertionLevel -> {Classes}, Model -> "SMQCD", ExcludeParticles -> Join[excluded, extra]];

bornDiagrams = generate[CreateTopologies[0, 2 -> 2], process, {V[5]}];
realDiagrams = DiagramSelect[
  generate[CreateTopologies[0, 2 -> 3], {F[2, {1}], F[3, {1}]} -> {F[2, {1}], F[3, {1}], V[5]}],
  !FreeQ[TreeFields[##], V[1, ___]] &];                      (* keep photon exchange, drop Z and gluon exchange *)
virtualDiagrams = DiagramSelect[
  generate[CreateTopologies[1, 2 -> 2, ExcludeTopologies -> {Tadpoles}], process],
  !FreeQ[LoopFields[##], V[5, ___]] && !FreeQ[TreeFields[##], V[1, ___]] &];   (* gluon in the loop, photon exchanged *)
vertexDiagrams = DiagramSelect[virtualDiagrams, WFCorrectionFields[##] === {} &];   (* the vertex graph alone *)
photonDiagrams = generate[CreateTopologies[0, 2 -> 2], {V[1], F[3, {1}]} -> {F[3, {1}], V[5]}, {F[2]}];

<|"Born" -> Length[bornDiagrams], "Real" -> Length[realDiagrams], "Virtual" -> Length[virtualDiagrams],
  "VirtualVertex" -> Length[vertexDiagrams], "Photon" -> Length[photonDiagrams]|>


(* ::Text:: *)
(*The two external-leg (self-energy) graphs among the three virtual diagrams are scaleless on massless external lines and vanish in dimensional regularization; only the vertex graph is used in step 2.  Expected counts: Born 1, Real 2, Virtual 3, VirtualVertex 1, Photon 2.*)


(* ::Input:: *)
Require[Length /@ {bornDiagrams, realDiagrams, virtualDiagrams, vertexDiagrams, photonDiagrams} === {1, 2, 3, 1, 2},
  "Unexpected number of diagrams"];
Paint[realDiagrams, ColumnsXRows -> {2, 1}, Numbering -> Simple, SheetHeader -> None, ImageSize -> {500, 220}]


(* ::Input:: *)
Paint[vertexDiagrams, ColumnsXRows -> {1, 1}, Numbering -> Simple, SheetHeader -> None, ImageSize -> {250, 220}]


(* ::Section:: *)
(*Amplitudes*)


(* ::Text:: *)
(*CreateFeynAmp gives the amplitude of each diagram with external spinors and polarization vectors in place (Truncated -> False).  FCFAConvert rewrites it in FeynCalc notation; with Contract -> False nothing is contracted or simplified, so the open spinor chains are exactly those of the diagram.  All momenta and Dirac matrices are D-dimensional (ChangeDimension -> D); masses are set to zero.*)


(* ::Input:: *)
amplitude[diagrams_] := CreateFeynAmp[diagrams, Truncated -> False, PreFactor -> 1, GaugeRules -> {}];
convert[diagrams_, incoming_, outgoing_, loops_: {}] := FCFAConvert[amplitude[diagrams],
  IncomingMomenta -> incoming, OutgoingMomenta -> outgoing, LoopMomenta -> loops,
  UndoChiralSplittings -> True, ChangeDimension -> D, List -> True, SMP -> True, Contract -> False,
  DropSumOver -> True, Prefactor -> 3 eq/2,
  FinalSubstitutions -> {SMP["m_e"] -> 0, SMP["m_u"] -> 0}] /. GaugeXi[_] -> 1;   (* Feynman gauge *)

feyncalcAmplitudes = <|
  "Born" -> convert[bornDiagrams, {l, k}, {lp, p}],
  "Real" -> convert[realDiagrams, {l, k}, {lp, p, r}],
  "VirtualVertex" -> convert[vertexDiagrams, {l, k}, {lp, p}, {loop}],
  "Photon" -> convert[photonDiagrams, {q, k}, {p, r}]|>;
feyncalcAmplitudes["Born"]


(* ::Section:: *)
(*The same amplitudes as plain records*)


(* ::Text:: *)
(*AmplitudeRecords writes each term of an amplitude as a coefficient (couplings, propagators Den[q] = 1/q^2, colour tensors, metric tensors, polarization vectors) times open fermion chains Chain[Ubar[p], GammaMatrix[mu], Slash[k - r], ..., Uspinor[k]].  This is the form FORM and Tracer start from; it is also the easiest form to compare with a hand-written amplitude.*)


(* ::Input:: *)
records = Map[AmplitudeRecords, feyncalcAmplitudes, {2}];
ShowRecords[Flatten[records["Real"]]]


(* ::Text:: *)
(*The two real-emission diagrams: emission from the incoming quark line (propagator Den[k - r]) and from the outgoing one (Den[p + r]), each with the photon propagator Den[k - p - r].*)


(* ::Input:: *)
ShowRecords[Flatten[records["VirtualVertex"]]]


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
SaveAmplitudes[<|"FeynCalc" -> feyncalcAmplitudes, "Records" -> records|>]
