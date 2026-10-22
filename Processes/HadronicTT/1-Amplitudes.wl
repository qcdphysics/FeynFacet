(* ::Title:: *)
(*Transversely polarized proton-proton classes (A_TT), step 1: diagrams and amplitudes of one class*)


(* ::Text:: *)
(*Single-inclusive parton production in collisions of two transversely polarized protons, a(pa) + b(pb) -> c(p) + X, at NLO in QCD for the five channel classes of classes.wl (5, 6, 7, 11, 12 in the numbering of Processes/Hadronic): both incoming partons are quarks or antiquarks carrying a transverse spin, the observed parton is unpolarized.  The variable class below selects the class; every step of this process starts with it, and the results are stored under Results/<route>/class-<n>/.*)
(**)
(*The diagrams are those of Processes/Hadronic (the same FeynArts generation, the same amplitude records): for the chosen class this step generates the Born process a b -> c d (when it exists), its one-loop corrections, and every real-emission process a b -> c d e of the class with its weight (identical unobserved partons, flavour multiplicities).  The transverse spin enters only in step 2, through the spin densities of the two incoming partons.*)
(**)
(*Momenta: pa, pb incoming; p the observed parton; k1, k2 the unobserved partons (k1 + k2 = pa + pb - p); loop the loop momentum.  Flavours: q = quark of flavour i (up quark of FeynArts' SMQCD model), qp = quark of flavour k (charm), g = gluon; closed quark loops run over all three generations of the model and are counted with the symbol Nf in step 2.*)


(* ::Section:: *)
(*Load the package and the class*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[]


(* ::Input:: *)
class = "5";
Get[FileNameJoin[{NotebookDirectory[], "classes.wl"}]];
Require[KeyExistsQ[HadronicTTClasses, ToExpression[class]], "Not a class of the transversely polarized process (5, 6, 7, 11, 12)", class];
spec = HadronicTTClasses[ToExpression[class]];
Print["Class ", class, ": ", spec["Name"]];
spec


(* ::Section:: *)
(*Generate the diagrams*)


(* ::Text:: *)
(*FeynArts fields: F[3,{1}] the up quark (flavour i), F[3,{2}] the charm quark (flavour k), V[5] the gluon, U[5] the gluon ghost (only in loops).  Everything else of the SMQCD model (photon, weak bosons, Higgs, leptons, down-type quarks, the third generation as an external particle) is excluded.  The observed parton is the first outgoing particle of each process; the tadpole topologies and the self-energy corrections on external legs are excluded from the one-loop generation (they are scaleless and vanish for massless partons).*)


(* ::Input:: *)
field = <|"q" -> F[3, {1}], "qb" -> -F[3, {1}], "qp" -> F[3, {2}], "qpb" -> -F[3, {2}], "g" -> V[5]|>;
excluded = {S[_], V[1 | 2 | 3], F[1 | 2 | 4]};
generate[loops_, in_, out_] := InsertFields[CreateTopologies[loops, Length[in] -> Length[out],
    ExcludeTopologies -> If[loops == 0, {}, {Tadpoles, WFCorrections}]],
  (field /@ in) -> (field /@ out), InsertionLevel -> {Classes}, Model -> "SMQCD",
  ExcludeParticles -> Join[excluded, If[loops == 0, {U[_]}, {}]]];
bornDiagrams = If[spec["Born"] === None, None, generate[0, spec["Born"][[1 ;; 2]], spec["Born"][[3 ;; 4]]]];
virtualDiagrams = If[spec["Virtual"], generate[1, spec["Born"][[1 ;; 2]], spec["Born"][[3 ;; 4]]], None];
realDiagrams = Association @ Table[StringRiffle[proc[[1]], " "] -> generate[0, proc[[1]][[1 ;; 2]], proc[[1]][[3 ;; 5]]], {proc, spec["Real"]}];
<|"Born" -> If[bornDiagrams === None, 0, Length[bornDiagrams]], "Virtual" -> If[virtualDiagrams === None, 0, Length[virtualDiagrams]],
  "Real" -> Map[Length, realDiagrams]|>


(* ::Input:: *)
If[bornDiagrams =!= None, Paint[bornDiagrams, ColumnsXRows -> {Min[Length[bornDiagrams], 4], 1}, Numbering -> Simple, SheetHeader -> None, ImageSize -> {170 Min[Length[bornDiagrams], 4], 170}]]


(* ::Section:: *)
(*Amplitudes*)


(* ::Text:: *)
(*FCFAConvert rewrites the FeynArts amplitudes in FeynCalc notation without contracting anything; quark masses are set to zero, the gauge parameter to one (Feynman gauge).  The one-loop amplitudes keep the loop momentum loop.  The third-generation quark running in closed loops appears with the symbolic FeynArts generation sum, which step 2 replaces by the number of flavours Nf.*)


(* ::Input:: *)
convert[diagrams_, in_, out_, loops_: {}] := FCFAConvert[CreateFeynAmp[diagrams, Truncated -> False, PreFactor -> 1],
  IncomingMomenta -> in, OutgoingMomenta -> out, LoopMomenta -> loops, UndoChiralSplittings -> True, ChangeDimension -> D,
  List -> True, SMP -> False, Contract -> False, DropSumOver -> False,
  FinalSubstitutions -> {FCGV["MU"] -> 0, FCGV["MC"] -> 0, FCGV["MT"] -> 0, SMP["g_s"] -> gs}] /. GaugeXi[_] -> 1 /. MQU[_] -> 0;   (* MQU: mass of the loop quark of generic generation *)
amplitudes = <|
  "Born" -> If[bornDiagrams === None, {}, convert[bornDiagrams, {pa, pb}, {p, k1}]],
  "Virtual" -> If[virtualDiagrams === None, {}, convert[virtualDiagrams, {pa, pb}, {p, k1}, {loop}]],
  "Real" -> Map[convert[#, {pa, pb}, {p, k1, k2}] &, realDiagrams]|>;
Length /@ Append[Map[Length, amplitudes["Real"]], "Virtual" -> amplitudes["Virtual"]]


(* ::Text:: *)
(*Closed quark loops: FeynArts sums the generation index of a closed fermion loop (SumOver[Index[Generation, n], 3]); the loop amplitudes are the same for every massless flavour, so each such amplitude counts Nf times.  Amplitudes with an external-generation sum do not occur because every external flavour is fixed.*)


(* ::Input:: *)
flavourSum = SumOver[Index[Generation, _], 3] -> Nf;
Require[FreeQ[amplitudes, SumOver[Index[Generation, _], _, External]], "An external generation sum appeared"];
amplitudes = amplitudes /. flavourSum /. SumOver[__] -> 1;      (* the colour sums are implicit in the colour algebra of step 2 *)
Count[amplitudes["Virtual"], Nf, Infinity]


(* ::Section:: *)
(*Records*)


(* ::Input:: *)
records = <|"Born" -> AmplitudeRecords /@ amplitudes["Born"], "Virtual" -> AmplitudeRecords /@ amplitudes["Virtual"],
  "Real" -> Map[AmplitudeRecords /@ # &, amplitudes["Real"]]|>;
Require[FreeQ[records, $Failed], "A record could not be formed"];
weights = Association @ Table[StringRiffle[proc[[1]], " "] -> proc[[2]], {proc, spec["Real"]}];
SaveAmplitudes[<|"FeynCalc" -> amplitudes, "Records" -> records, "Class" -> class, "Spec" -> spec, "Weights" -> weights|>, "class-" <> class];
ShowRecords[Flatten[records["Born"]]]
