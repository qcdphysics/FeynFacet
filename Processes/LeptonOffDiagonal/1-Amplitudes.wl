(* ::Title:: *)
(*Lepton-quark and lepton-gluon scattering with a different tagged parton, step 1: diagrams and amplitudes*)


(* ::Text:: *)
(*The two off-diagonal channels of single-inclusive parton production in lepton-proton collisions,*)
(*   l(l) + q(k) -> g(p) + X   (the gluon is observed)   and   l(l) + g(k) -> q(p) + X   (the quark is observed),*)
(*which begin at order alphaEM^2 alphaS: there is no Born term, so the NLO coefficient consists of the real emission l q -> l q g (two diagrams, the gluon tagged) or l g -> l q qbar (two diagrams, the quark tagged), the collinear counterterms (fragmentation q -> g for the first channel, the parton distribution g -> q for the second), and the photon-in-lepton subtraction with the Born processes gamma q -> g q and gamma g -> q qbar.  Both the spin average (UU) and the double-longitudinal helicity difference (LL) are computed; the results are compared with Hinderer, Schlegel and Vogelsang, arXiv:1505.06415 (Eqs. 27, 28) and arXiv:1703.10872 (Eqs. 13, 14).*)
(**)
(*Everything follows the lepton-quark process (Processes/LeptonQuark): the same momenta (l, k incoming; lp the scattered lepton, p the observed parton, r the other unobserved parton), invariants, scalar-integral families and master integrals; only the tagged particle changes, and with it the matrix elements and the counterterms.  The photon Born terms use the photon momentum q and the same observed parton p.  The Born process l q -> l q (one diagram) is generated as well: the collinear counterterms of both channels are built from it.*)


(* ::Section:: *)
(*Load the package*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[]


(* ::Section:: *)
(*Generate the diagrams*)


(* ::Text:: *)
(*FeynArts' SMQCD model: F[2,{1}] the electron, F[3,{1}] the up quark, V[1] the photon, V[5] the gluon; the observed parton is the first outgoing particle.  Only photon exchange between the lepton and the quark line is kept (Z exchange is excluded by the particle list).*)


(* ::Input:: *)
excluded = {S[_], U[_], V[2 | 3], F[1 | 4]};
generate[topologies_, fields_, extra_: {}] := InsertFields[topologies, fields,
  InsertionLevel -> {Classes}, Model -> "SMQCD", ExcludeParticles -> Join[excluded, extra]];
bornDiagrams = DiagramSelect[generate[CreateTopologies[0, 2 -> 2], {F[2, {1}], F[3, {1}]} -> {F[3, {1}], F[2, {1}]}], !FreeQ[TreeFields[##], V[1, ___]] &];
realQG = DiagramSelect[generate[CreateTopologies[0, 2 -> 3], {F[2, {1}], F[3, {1}]} -> {V[5], F[2, {1}], F[3, {1}]}], !FreeQ[TreeFields[##], V[1, ___]] &];
realGQ = DiagramSelect[generate[CreateTopologies[0, 2 -> 3], {F[2, {1}], V[5]} -> {F[3, {1}], F[2, {1}], -F[3, {1}]}], !FreeQ[TreeFields[##], V[1, ___]] &];
photonQG = generate[CreateTopologies[0, 2 -> 2], {V[1], F[3, {1}]} -> {V[5], F[3, {1}]}, {F[2]}];
photonGQ = generate[CreateTopologies[0, 2 -> 2], {V[1], V[5]} -> {F[3, {1}], -F[3, {1}]}, {F[2]}];
Require[Length /@ {bornDiagrams, realQG, realGQ, photonQG, photonGQ} === {1, 2, 2, 2, 2}, "Unexpected number of diagrams"];
Paint[realGQ, ColumnsXRows -> {2, 1}, Numbering -> Simple, SheetHeader -> None, ImageSize -> {500, 220}]


(* ::Section:: *)
(*Amplitudes*)


(* ::Text:: *)
(*Converted exactly as in the lepton-quark process: FCFAConvert with Contract -> False, D-dimensional, massless, Feynman gauge, and the general quark charge eq in units of ee through the prefactor 3 eq/2 (the model's up quark has charge 2/3).  The outgoing momenta are {p, lp, r}: the observed parton first, then the scattered lepton, then the unobserved parton (the quark of l q -> g or the antiquark of l g -> q); for the photon processes {p, r} with the incoming photon momentum q; for the Born process {p, lp}.*)


(* ::Input:: *)
convert[diagrams_, incoming_, outgoing_] := FCFAConvert[CreateFeynAmp[diagrams, Truncated -> False, PreFactor -> 1, GaugeRules -> {}],
  IncomingMomenta -> incoming, OutgoingMomenta -> outgoing, UndoChiralSplittings -> True, ChangeDimension -> D,
  List -> True, SMP -> True, Contract -> False, DropSumOver -> True, Prefactor -> 3 eq/2,
  FinalSubstitutions -> {SMP["m_e"] -> 0, SMP["m_u"] -> 0}] /. GaugeXi[_] -> 1;      (* Feynman gauge *)
amplitudes = <|
  "Born" -> convert[bornDiagrams, {l, k}, {p, lp}],
  "RealQG" -> convert[realQG, {l, k}, {p, lp, r}],
  "RealGQ" -> convert[realGQ, {l, k}, {p, lp, r}],
  "PhotonQG" -> convert[photonQG, {q, k}, {p, r}],
  "PhotonGQ" -> convert[photonGQ, {q, k}, {p, r}]|>;
amplitudes["RealQG"]


(* ::Input:: *)
records = Map[AmplitudeRecords, amplitudes, {2}];
Require[FreeQ[records, $Failed], "A record could not be formed"];
SaveAmplitudes[<|"FeynCalc" -> amplitudes, "Records" -> records|>];
ShowRecords[Flatten[records["RealGQ"]]]
