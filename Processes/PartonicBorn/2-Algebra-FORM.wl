(* ::Title:: *)
(*Partonic Born cross sections, FORM route: the same table with FORM and its colour library*)


(* ::Text:: *)
(*Same task and same output format as 1-Born.wl (the FeynCalc route), with FORM doing all the spin, Lorentz and colour algebra.  Input: the amplitude records of all 29 channels written by 1-Born.wl (Results/Amplitudes.wl), and its table Results/Born.wl for the comparison.  Output: Results/FORM/Born.wl.*)
(**)
(*Mathematica only writes the amplitudes to a file (amplitudes.inc, one local expression Amp1 ... Amp29 per channel), runs the FORM program form/Born.frm twice (SPIN = U, spin sum; SPIN = P, helicity difference of the incoming partons) and reads the squares back.  Open form/common.h and form/Born.frm to follow the algebra; for every channel the program*)
(**)
(*  1. multiplies the amplitude by its complex conjugate (reversed gamma strings, transposed generators, primed summed indices);*)
(*  2. sums over all colours with FORM's colour library color.h (#call docolor), before any Dirac algebra, and specializes the group invariants to SU(Nc) with Tr(T^a T^b) = delta^{ab}/2, the normalization of the FeynCalc route;*)
(*  3. closes the fermion lines into traces with the spin sums of the quarks: one trace or two, depending on how the quark lines of the two interfering diagrams connect;*)
(*  4. inserts the spin densities: p-slash gamma5 (-p-slash gamma5 for an antiquark) for a polarized incoming quark, with the HVBM gamma5 written through a four-dimensional epsilon tensor, gamma5 = -(i/24) eps_{mu nu rho sigma} gamma^mu gamma^nu gamma^rho gamma^sigma, normalized so that Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = -4 i eps^{abcd} (FeynCalc with $LeviCivitaSign = -1; checked below); the circular density i eps^{mu nu l r}/(2 l.r) for a polarized incoming gluon l, r being the other incoming momentum; for every other gluon the physical polarization sum -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r) with the reference vectors of 1-Born.wl, r = pb for l = pa and r = pa for l = pb, p, pd;*)
(*  5. evaluates the D-dimensional traces (tracen), contracts the indices (indices inside epsilon tensors are four-dimensional) and inserts the kinematics s = 2 pa.pb, t = -2 pa.p, u = -2 pb.p.  FORM contracts two epsilon tensors with the Euclidean sign (eps_{abcd} eps_{abcd} = +24 instead of the Minkowski -24); every LL term contains exactly two of them, one per polarized incoming parton, so the LL squares are multiplied by -1.*)
(**)
(*The spin and colour averages of the incoming partons are applied here afterwards, exactly as in 1-Born.wl.  The FORM squares must equal the FeynCalc table exactly, for every channel, both spin cases, arbitrary D and arbitrary number of colours Nc.*)


(* ::Section:: *)
(*Load and export*)


(* ::Text:: *)
(*Load the FeynFacet package (amplitude records, the FORM interface WriteFORMAmplitudes, RunFORM, ReadFORM, and the checks Require).*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*The amplitude records of 1-Born.wl, channel by channel, are written to amplitudes.inc; the FORM expression Amp<i> is the amplitude of the i-th channel of this list, and the table shows how many records (terms with one colour factor and one set of fermion chains) each amplitude has.  With "ColorLibrary" -> True the colour objects are written in the conventions of color.h: generators T(i,j,a) = (T^a)_{ij} with the adjoint index last, structure constants f(a,b,c), d_ for both colour deltas, quark colours with dimension NR and gluon colours with dimension NA.  The declarations of all symbols, vectors and indices go to declarations.inc, in the same directory.*)


(* ::Input:: *)
amplitudes = LoadAmplitudes[];
channels = amplitudes["Channels"];
work = ResultsDirectory["FORM"];
WriteFORMAmplitudes[Association @ Table["Amp" <> ToString[i] -> amplitudes["Records"][channels[[i]]], {i, Length[channels]}], work, "ColorLibrary" -> True];
TableForm[Table[{"Amp" <> ToString[i], channels[[i]], Length[amplitudes["Records"][channels[[i]]]]}, {i, Length[channels]}],
  TableHeadings -> {None, {"FORM expression", "channel", "records"}}]


(* ::Text:: *)
(*As an example, the amplitude of q g -> q in FORM notation (three diagrams: s- and u-channel quark exchange, t-channel gluon exchange with the three-gluon vertex):*)


(* ::Input:: *)
First[Select[StringSplit[Import[FileNameJoin[{work, "amplitudes.inc"}], "Text"], "L "],
  StringStartsQ[#, "Amp" <> ToString[First[FirstPosition[channels, "q g -> q"]]] <> " "] &]]


(* ::Section:: *)
(*Run the FORM program*)


(* ::Text:: *)
(*Born.frm loops over the channels and writes the squares, in channel order, to BornU.m (SPIN = U) and BornP.m (SPIN = P).  It runs with TFORM, the multi-threaded version of FORM, on three worker threads; a FORM error stops the notebook.  The output shows the two result files.*)


(* ::Input:: *)
program = FileNameJoin[{NotebookDirectory[], "form", "Born.frm"}];
Do[RunFORM[program, work, "Workers" -> 3, "Define" -> {"SPIN=" <> spin, "CHANNELS=" <> ToString[Length[channels]]}], {spin, {"U", "P"}}];
FileNames["Born*.m", work]


(* ::Text:: *)
(*Born.frm (SPIN = P) prints G5Check, the difference between the trace Tr[gamma^a gamma^b gamma^c gamma^d gamma5] computed with the explicit HVBM gamma5 and -4 i eps^{abcd} (FeynCalc's convention with $LeviCivitaSign = -1); it must vanish.*)


(* ::Input:: *)
formLog = Import[FileNameJoin[{work, "Born-SPIN=P-CHANNELS=" <> ToString[Length[channels]] <> ".log"}], "Text"];
g5Check = StringDelete[First[StringCases[formLog, "G5Check =" ~~ Shortest[x__] ~~ ";" :> x]], Whitespace | "(" | ")"];
Require[g5Check === "0", "gamma5 normalization check failed in Born.frm", g5Check];
"G5Check = " <> g5Check


(* ::Section:: *)
(*Averages and comparison with the FeynCalc table*)


(* ::Text:: *)
(*FORM returns |M|^2 summed over the colours of all four partons and over the spins (UU), or with the spin densities of the two incoming partons (LL), in terms of s, t, u, D, Nc and gs.  The averages are those of 1-Born.wl: 1/(2 Nc) for an incoming quark (2 spin states times Nc colours; for LL the 1/2 turns p-slash gamma5 into the helicity density p-slash gamma5/2); 1/((D-2)(Nc^2-1)) for an unpolarized incoming gluon (D-2 polarization states in D dimensions), 1/(Nc^2-1) for a polarized one (its 1/2 is in the circular density).  The result is divided by gs^4 and written with u = -s-t; the output is the channel q g -> q as an example.*)


(* ::Input:: *)
incoming[channel_String] := Take[StringSplit[channel], 2];
average[channel_, spinCase_] := Times @@ (Switch[#, "g", If[spinCase === "LL", 1, 1/(D - 2)]/(Nc^2 - 1), _, 1/(2 Nc)] & /@ incoming[channel]);
formSquares = <|"UU" -> ReadFORM[work, "BornU.m"], "LL" -> ReadFORM[work, "BornP.m"]|>;
Require[Length /@ Values[formSquares] === {Length[channels], Length[channels]}, "FORM did not return one square per channel"];
bornNc = Association @ Table[spinCase -> Association @ Table[channels[[i]] ->
    Factor[Together[average[channels[[i]], spinCase] formSquares[spinCase][[i]]/gs^4 /. u -> -s - t]], {i, Length[channels]}], {spinCase, {"UU", "LL"}}];
bornNc["UU"]["q g -> q"]


(* ::Text:: *)
(*The comparison: the FeynCalc table is written with CF and CA; with CF = (Nc^2-1)/(2 Nc) and CA = Nc both tables are rational functions of s, t, D and Nc, and their difference must vanish identically, channel by channel.*)


(* ::Input:: *)
feynCalc = LoadResult["Born"]["Born"];
toNc = {CF -> (Nc^2 - 1)/(2 Nc), CA -> Nc};
difference = Association @ Table[spinCase -> Association @ Table[channel ->
    Together[bornNc[spinCase][channel] - (feynCalc[spinCase][channel] /. toNc)], {channel, channels}], {spinCase, {"UU", "LL"}}];
Require[Union[Flatten[Values /@ Values[difference]]] === {0}, "The FORM table differs from the FeynCalc table", difference];
TableForm[Table[{channel, difference["UU"][channel], difference["LL"][channel]}, {channel, channels}],
  TableHeadings -> {None, {"channel", "UU: FORM - FeynCalc", "LL: FORM - FeynCalc"}}]


(* ::Section:: *)
(*Save*)


(* ::Text:: *)
(*Written in the format of the FeynCalc route: CF and CA restored as in 1-Born.wl ((Nc-1)(Nc+1) = 2 Nc CF, then Nc = CA).  The saved table is checked once more against the FeynCalc one; the output lists the saved channels.*)


(* ::Input:: *)
colourFactors[expr_] := Factor[expr] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
born = Map[colourFactors, bornNc, {2}];
Require[Union[Flatten[Table[Together[(born[sc][ch] - feynCalc[sc][ch]) /. toNc], {sc, {"UU", "LL"}}, {ch, channels}]]] === {0},
  "The saved FORM table differs from the FeynCalc table"];
result = <|"Born" -> born, "Channels" -> channels, "Variables" -> "s, t (u = -s-t), D, CF, CA; divided by gs^4"|>;
SaveResult["Born", result, "FORM"];
Keys[born["UU"]]
