(* ::Title:: *)
(*TMD beam and soft functions, step 1: the Wilson-line model, diagrams and amplitudes*)


(* ::Text:: *)
(*This process computes the NLO matching coefficients of the quark and gluon TMD beam functions (I_qq, I_qg, I_gq, I_gg), the NLO vacuum soft function (quark and gluon), and from them the TMD parton distributions f_ij = I_ij Sqrt[S_i].  All of them are matrix elements of light-cone operators with Wilson lines, not scattering amplitudes.  This step shows how such operator matrix elements are generated with FeynArts, checks that the generated currents are gauge invariant, draws the diagrams, and converts the amplitudes into the plain records that the three algebra routes read in step 2.*)
(**)
(*The operators.  With two light-like vectors n and nb (n^2 = nb^2 = 0, n.nb = 2) and an incoming parton of momentum p along nb (so that n.p is its large light-cone component),*)
(**)
(*   quark:  O_q(b) = qbar(b) W_n(b)^dagger n-slash W_n(0) q(0),      W_n(x) = P exp(-i gs Integral_0^Infinity ds n.A(x + s n)),*)
(*   gluon:  O_g(b) = n_nu F^{nu mu, a}(b) W_n^A(b)^{ab dagger} W_n^A(0)^{bc} n_rho F^{rho}_mu^{c}(0),*)
(*   soft:   O_S(b) = Tr[ (W_n^dagger W_nb)(b)^dagger (W_n^dagger W_nb)(0) ] / Nc,*)
(**)
(*and the beam function is B(x, b_T) = Integral over b^- of e^{-i x n.p b^-} <p| O(b) |p>, normalized so that B = delta(1-x) at lowest order; the soft function is the vacuum expectation value of O_S, equal to 1 at lowest order (see docs/CONVENTIONS.md).  Inserting a complete set of states, <p| O(b) |p> = Sum over X of <p| O^dagger-part |X> <X| O-part(0) |p>, so what must be generated is the amplitude <X| W_n(0) q(0) |p>, i.e. the process p -> X + "operator".  The Wilson line runs from the operator point to infinity.*)
(**)
(*The model.  models/WilsonLineQCD.mod (classes) and models/WilsonLine.gen (generic) implement QCD plus:*)
(*  - Wilson-line fields S[1] = W_n, S[2] = W_nb (fundamental), S[3] = W_n^A, S[4] = W_nb^A (adjoint), drawn as dashed lines with an arrow.  Their Feynman rules are those of a heavy quark with velocity n (HQET): propagator i/(n.q), vertex -i gs T^a n^mu.  The far end of a Wilson line at infinity is an external leg carrying no momentum.*)
(*  - "sink" fields that represent the operator vertex: O_q = F[2] (carries the open Dirac index of the quark operator), O_g = V[2] (the open Lorentz index of F^{nu mu}), O_S = S[5].  They are external legs with a single coupling each: O_q q W_n, O_g g W_n^A (plus the two-gluon term of the field strength), O_S W_n W_nb.  The couplings gOp and gSoft only mark the operator insertion and are divided out.*)
(*In the diagrams the parton enters from the left, the operator leg points up, the Wilson line leaves to the right.  The precise Feynman rules and sign conventions are explained in the two model files.*)
(**)
(*Momenta: p incoming parton; k the emitted parton (gluon, or quark for the g -> q channel); q = p - k the momentum flowing into the operator; w = 0 the momentum of the Wilson-line end at infinity.  For the soft function the operator O_S carries the momentum k of the emitted gluon and both Wilson-line ends carry zero momentum.*)


(* ::Section:: *)
(*Load the package*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[]


(* ::Input:: *)
$ModelPath = Prepend[$ModelPath, FileNameJoin[{NotebookDirectory[], "models"}]];
$ShapeDataDir = FileNameJoin[{NotebookDirectory[], "models", "ShapeData"}];   (* layouts of the diagrams, see models/ShapeData *)
FileNames["*.mod" | "*.gen", FileNameJoin[{NotebookDirectory[], "models"}]]


(* ::Section:: *)
(*Generate the diagrams*)


(* ::Text:: *)
(*Every operator matrix element is generated as the process "parton -> operator leg + Wilson-line end + emitted parton".  In FeynArts notation F[1] is the quark, V[1] the gluon, F[2], V[2], S[5] the operator legs and S[1], S[3] the Wilson-line ends; a minus sign denotes the antiparticle (the outgoing antiquark in g -> q qbar, the Wilson line W_nb of the soft operator whose colour charge flows in from infinity).  The one-loop (virtual) corrections are generated too, only to be drawn: with p^2 = n^2 = 0 the loop integrals have no scale and vanish in dimensional regularization, so the NLO result comes from the real emission alone.*)


(* ::Input:: *)
generate[loops_, process_] := InsertFields[CreateTopologies[loops, Length[process[[1]]] -> Length[process[[2]]],
    ExcludeTopologies -> If[loops == 0, {}, {Tadpoles, WFCorrections}]], process,
  InsertionLevel -> {Classes}, Model -> "WilsonLineQCD", GenericModel -> "WilsonLine"];

diagrams = <|
  "BornQ"  -> generate[0, {F[1]} -> {F[2], S[1]}],                (* q -> O_q W_n *)
  "BeamQQ" -> generate[0, {F[1]} -> {F[2], S[1], V[1]}],          (* q -> O_q W_n g *)
  "BeamQG" -> generate[0, {V[1]} -> {F[2], S[1], -F[1]}],         (* g -> O_q W_n qbar *)
  "BornG"  -> generate[0, {V[1]} -> {V[2], S[3]}],                (* g -> O_g W_n^A *)
  "BeamGQ" -> generate[0, {F[1]} -> {V[2], S[3], F[1]}],          (* q -> O_g W_n^A q *)
  "BeamGG" -> generate[0, {V[1]} -> {V[2], S[3], V[1]}],          (* g -> O_g W_n^A g *)
  "SoftQ"  -> generate[0, {S[5]} -> {S[1], -S[2], V[1]}],         (* O_S -> W_n W_nb g, fundamental *)
  "SoftG"  -> generate[0, {S[5]} -> {S[3], -S[4], V[1]}]|>;       (* the same with adjoint Wilson lines *)
virtualDiagrams = DiagramSelect[generate[1, {F[1]} -> {F[2], S[1]}], FreeQ[LoopFields[##], F[2] | V[2] | S[5]] &];
Length /@ Append[diagrams, "VirtualQ" -> virtualDiagrams]


(* ::Text:: *)
(*Expected counts.  Quark operator: one Born graph; two real graphs (gluon from the quark, gluon from the Wilson line); two g -> q graphs (the gluon splits into q qbar, or the Wilson line absorbs the gluon and the operator creates the antiquark).  Gluon operator: one Born graph; one q -> g graph; four g -> g graphs (three-gluon vertex, two-gluon term of the field strength, emission from and absorption by the Wilson line).  Soft: two graphs each (emission from either line).  Virtual: the one vertex graph (self-energy insertions are excluded; the selection removes graphs in which an operator leg would propagate inside the loop).*)


(* ::Input:: *)
Require[Length /@ Values[diagrams] === {1, 2, 2, 1, 1, 4, 2, 2} && Length[virtualDiagrams] === 1, "Unexpected number of diagrams"];


(* ::Section:: *)
(*Draw the diagrams*)


(* ::Text:: *)
(*The generated figures are written to Results/figures/ as PDF files.  The figures/ directory contains the shipped previews.  The layout of each topology is prescribed in models/ShapeData (FeynArts' own shape files, editable with Shape[]).*)


(* ::Input:: *)
figures = FileNameJoin[{ResultsDirectory[], "figures"}];
If[!DirectoryQ[figures], CreateDirectory[figures, CreateIntermediateDirectories -> True]];
draw[name_, tops_] := Module[{n = Length[tops], columns, picture},
  columns = Min[n, 4];
  picture = Paint[tops, ColumnsXRows -> {columns, Ceiling[n/columns]}, Numbering -> Simple, SheetHeader -> None,
    AutoEdit -> False, ImageSize -> {170 columns, 170 Ceiling[n/columns]}, DisplayFunction -> Identity];
  Export[FileNameJoin[{figures, name <> ".pdf"}], picture];
  picture];
draw["VirtualQ", virtualDiagrams];
Column[KeyValueMap[Column[{Style[#1, Bold], draw[#1, #2]}] &, diagrams]]


(* ::Section:: *)
(*Amplitudes*)


(* ::Text:: *)
(*CreateFeynAmp gives the amplitudes; FCFAConvert rewrites them in FeynCalc notation without contracting anything.  The operator leg keeps the symbolic momentum q (its spinor or polarization vector is the placeholder for the open index); everywhere else q is replaced by p - k and the Wilson-line end momentum w by 0.*)


(* ::Input:: *)
convert[tops_, incoming_, outgoing_] := Module[{amp},
  amp = FCFAConvert[CreateFeynAmp[tops, Truncated -> False, PreFactor -> 1],
    IncomingMomenta -> incoming, OutgoingMomenta -> outgoing, UndoChiralSplittings -> True,
    ChangeDimension -> D, List -> True, SMP -> False, Contract -> False, DropSumOver -> True];
  amp /. {Spinor[Momentum[q, d___], r__] :> Spinor[Momentum[qOperator, d], r], Polarization[q, r__] :> Polarization[qOperator, r]} /.
    {q -> p - k, w -> 0} /. qOperator -> q];
amplitudes = <|
  "BornQ"  -> convert[diagrams["BornQ"], {p}, {q, w}] /. k -> 0,       (* the operator carries q = p *)
  "BeamQQ" -> convert[diagrams["BeamQQ"], {p}, {q, w, k}],
  "BeamQG" -> convert[diagrams["BeamQG"], {p}, {q, w, k}],
  "BornG"  -> convert[diagrams["BornG"], {p}, {q, w}] /. k -> 0,
  "BeamGQ" -> convert[diagrams["BeamGQ"], {p}, {q, w, k}],
  "BeamGG" -> convert[diagrams["BeamGG"], {p}, {q, w, k}],
  "SoftQ"  -> convert[diagrams["SoftQ"], {k}, {w, w, k}],
  "SoftG"  -> convert[diagrams["SoftG"], {k}, {w, w, k}]|>;
amplitudes["BeamQQ"]


(* ::Text:: *)
(*Each amplitude is linear in the operator marker gOp or gSoft; FAEikonal[q, n] stands for the eikonal denominator 1/(n.q).*)


(* ::Input:: *)
Require[And @@ Flatten[Map[Exponent[# /. Dot -> Times, gOp] + Exponent[# /. Dot -> Times, gSoft] === 1 &, Values[amplitudes], {2}]],
  "An amplitude is not linear in the operator marker"];


(* ::Section:: *)
(*Check: gauge invariance of the generated currents*)


(* ::Text:: *)
(*The operators are gauge invariant, so the sum of the graphs of each process must vanish when the polarization vector of a gluon is replaced by its momentum, the other external particles being physical: k_mu J^mu = 0.  This is a non-trivial check of the eikonal Feynman rules, of the orientation of the Wilson lines and of the relative signs of the field-strength couplings: the quark-emission and Wilson-line-emission graphs cancel only if the eikonal vertex, the eikonal propagator and the operator vertex are consistent.  The operator index is left open, so the identity holds for every component of the operator.*)
(**)
(*Kinematics: p^2 = k^2 = n^2 = nb^2 = 0, n.nb = 2; e is the polarization vector of the other gluon, transverse to its momentum.*)


(* ::Input:: *)
FCClearScalarProducts[];
SPD[p] = 0; SPD[k] = 0; SPD[n] = 0; SPD[nb] = 0; SPD[n, nb] = 2;
SPD[p, n] = pn; SPD[k, n] = kn; SPD[p, k] = pk; SPD[k, nb] = kb; SPD[p, nb] = pb;
SPD[e, n] = en; SPD[e, nb] = eb; SPD[e, p] = ep; SPD[e, k] = ek; SPD[e] = ee2;
eikonal = FAEikonal[mom_, v_] :> 1/ExpandScalarProduct[SPD[v, mom]];
current[name_, replacements_] := Factor[DiracSimplify[Contract[ExpandScalarProduct[FeynAmpDenominatorExplicit[
  Total[amplitudes[name]] /. eikonal /. replacements]]]]];
wardIdentities = <|
  "BeamQQ: gluon k" -> current["BeamQQ", Momentum[Polarization[k, -I], D] -> Momentum[k, D]],
  "BeamQG: gluon p" -> current["BeamQG", Momentum[Polarization[p, I], D] -> Momentum[p, D]],
  "BeamGG: gluon k, gluon p physical" -> current["BeamGG", {Momentum[Polarization[q, -I], D] -> LorentzIndex[mu, D],
    Momentum[Polarization[k, -I], D] -> Momentum[k, D], Momentum[Polarization[p, I], D] -> Momentum[e, D]}] /. ep -> 0,
  "BeamGG: gluon p, gluon k physical" -> current["BeamGG", {Momentum[Polarization[q, -I], D] -> LorentzIndex[mu, D],
    Momentum[Polarization[p, I], D] -> Momentum[p, D], Momentum[Polarization[k, -I], D] -> Momentum[e, D]}] /. ek -> 0,
  "SoftQ: gluon k" -> current["SoftQ", Momentum[Polarization[k, -I], D] -> Momentum[k, D]],
  "SoftG: gluon k" -> current["SoftG", Momentum[Polarization[k, -I], D] -> Momentum[k, D]]|>;
Require[And @@ (# === 0 & /@ Values[wardIdentities]), "A Ward identity fails", wardIdentities];
wardIdentities


(* ::Text:: *)
(*(The q -> g channel has no external gluon; its single graph needs no check.)  A deliberately wrong model, e.g. with the opposite sign of the eikonal vertex or of the two-gluon operator coupling, fails these identities; they were used to fix the conventions of the model files.*)


(* ::Section:: *)
(*Records*)


(* ::Text:: *)
(*AmplitudeRecords turns each amplitude into a list of records, <|"Coefficient" -> ..., "Chains" -> {Chain[...]}|> (vocabulary in FeynFacet/Amplitudes.wl).  The open Dirac index of the quark operator appears as the spinor Ubar[q] on the left of the chain; the open Lorentz index of the gluon operator as the polarization vector PolVector[q, mu, -1]; the eikonal denominators as Eikonal[k, n] = 1/(n.k).*)


(* ::Input:: *)
FCClearScalarProducts[];           (* the kinematics of the Ward check must not enter the records *)
records = Map[AmplitudeRecords, amplitudes, {2}];
ShowRecords[Flatten[records["BeamQQ"]]]


(* ::Input:: *)
ShowRecords[Flatten[records["BeamGG"]]]


(* ::Input:: *)
ShowRecords[Flatten[records["SoftQ"]]]


(* ::Input:: *)
SaveAmplitudes[<|"FeynCalc" -> amplitudes, "Records" -> records|>];
FileNameJoin[{ResultsDirectory[], "Amplitudes.wl"}]
