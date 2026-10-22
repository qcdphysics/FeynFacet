(* ::Title:: *)
(*Lepton-quark scattering, step 2 (FORM route): spin and colour algebra*)


(* ::Text:: *)
(*Same task and same output format as the FeynCalc route, now with FORM doing the Dirac algebra.  Input: the amplitude records of step 1.  Output: Results/FORM/Algebra.wl.*)
(**)
(*Mathematica only writes the amplitudes to a file (amplitudes.inc, in the vocabulary described in form/common.h), runs the four FORM programs in the form/ directory, and reads their output back.  Open those programs to follow the algebra: Born.frm, Real.frm and Photon.frm square the amplitudes (conjugation, spin sums, density matrices, gamma traces, polarization and colour sums, kinematics), Virtual.frm projects the vertex graph onto the Born current.  Each program is run twice, with the preprocessor variable SPIN = U (spin average) and SPIN = P (helicity difference).*)


(* ::Section:: *)
(*Load and export*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
records = LoadAmplitudes[]["Records"];
work = ResultsDirectory["FORM"];
WriteFORMAmplitudes[<|"Born1" -> records["Born"][[1]], "Real1" -> records["Real"][[1]], "Real2" -> records["Real"][[2]],
  "Vertex1" -> records["VirtualVertex"][[1]], "Photon1" -> records["Photon"][[1]], "Photon2" -> records["Photon"][[2]]|>, work]


(* ::Input:: *)
Import[FileNameJoin[{work, "amplitudes.inc"}], "Text"]


(* ::Section:: *)
(*Run the FORM programs*)


(* ::Input:: *)
program[name_] := FileNameJoin[{NotebookDirectory[], "form", name}];
Do[RunFORM[program[name], work, "Define" -> {"SPIN=" <> spin}], {name, {"Born.frm", "Real.frm", "Photon.frm"}}, {spin, {"U", "P"}}];
RunFORM[program["Virtual.frm"], work, "Define" -> {"SPIN=U"}];
FileNames["*.m", work]


(* ::Text:: *)
(*Born.frm (SPIN=P) prints G5Check, the difference between the trace Tr[gamma^a gamma^b gamma^c gamma^d gamma5] computed with the explicit gamma5 and -4 i eps^{abcd}; it must vanish.*)


(* ::Input:: *)
bornLog = Import[FileNameJoin[{work, "Born-SPIN=P.log"}], "Text"];
g5Check = StringDelete[First[StringCases[bornLog, "G5Check =" ~~ Shortest[x__] ~~ ";" :> x]], Whitespace | "(" | ")"];
Require[g5Check === "0", "gamma5 normalization check failed in Born.frm", g5Check];
"G5Check = " <> g5Check


(* ::Section:: *)
(*Read the results*)


(* ::Text:: *)
(*FORM writes its results in Mathematica syntax.  Born and real squares are read as functions of D (through the symbol D of the declarations), the photon squares get the photon average 1/(D-2) and the Mandelstam names sh, th, uh; the virtual projection is the ratio of the numerator trace and the Born-current trace, in the family VirtualIntegral[n0,n1,n2], times the colour factor CF.*)


(* ::Input:: *)
born = <|"Unpolarized" -> ReadFORM[work, "BornU.m"], "Polarized" -> ReadFORM[work, "BornP.m"]|>


(* ::Input:: *)
Require[Simplify[(born["Unpolarized"] - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];
Require[Simplify[(born["Polarized"] - 2 ee^4 eq^2 (s^2 - u^2)/t^2) /. u -> -s - t] === 0, "Polarized Born square differs from the expected result"];
real = <|"Unpolarized" -> ReadFORM[work, "RealU.m"], "Polarized" -> ReadFORM[work, "RealP.m"]|>;
Factor[real["Polarized"] /. spe3 -> 0]


(* ::Input:: *)
photon = <|"Unpolarized" -> ReadFORM[work, "PhotonU.m"]/(D - 2), "Polarized" -> ReadFORM[work, "PhotonP.m"]|>


(* ::Input:: *)
{numerator, bornCurrent} = ReadFORM[work, "Virtual.m"];
virtualFamily = Factor[numerator/bornCurrent]


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Photon" -> photon, "VirtualFamily" -> virtualFamily|>;
SaveResult["Algebra", algebra, "FORM"];
Keys[algebra]
