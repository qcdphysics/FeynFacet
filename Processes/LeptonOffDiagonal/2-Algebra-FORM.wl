(* ::Title:: *)
(*Lepton-quark and lepton-gluon scattering with a different tagged parton, step 2 (FORM route): spin and colour algebra*)


(* ::Text:: *)
(*Same task and same output format as the FeynCalc route, now with FORM doing the Dirac algebra.  Input: the amplitude records of step 1.  Output: Results/FORM/Algebra.wl.*)
(**)
(*Mathematica only writes the amplitudes to a file (amplitudes.inc, in the vocabulary described in form/common.h), runs the FORM programs in the form/ directory, and reads their output back.  Open those programs to follow the algebra: Born.frm, RealQG.frm, RealGQ.frm, PhotonQG.frm and PhotonGQ.frm square the amplitudes (conjugation, spin sums, density matrices, gamma traces, polarization and colour sums, kinematics).  Each program is run twice, with the preprocessor variable SPIN = U (spin average) and SPIN = P (helicity difference).*)


(* ::Section:: *)
(*Load and export*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
records = LoadAmplitudes[]["Records"];
work = ResultsDirectory["FORM"];
WriteFORMAmplitudes[<|"Born1" -> records["Born"][[1]],
  "RealQG1" -> records["RealQG"][[1]], "RealQG2" -> records["RealQG"][[2]], "RealGQ1" -> records["RealGQ"][[1]], "RealGQ2" -> records["RealGQ"][[2]],
  "PhotonQG1" -> records["PhotonQG"][[1]], "PhotonQG2" -> records["PhotonQG"][[2]], "PhotonGQ1" -> records["PhotonGQ"][[1]], "PhotonGQ2" -> records["PhotonGQ"][[2]]|>, work]


(* ::Input:: *)
Import[FileNameJoin[{work, "amplitudes.inc"}], "Text"]


(* ::Section:: *)
(*Run the FORM programs*)


(* ::Input:: *)
program[name_] := FileNameJoin[{NotebookDirectory[], "form", name}];
Do[RunFORM[program[name], work, "Define" -> {"SPIN=" <> spin}], {name, {"Born.frm", "RealQG.frm", "RealGQ.frm", "PhotonQG.frm", "PhotonGQ.frm"}}, {spin, {"U", "P"}}];
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
(*FORM writes its results in Mathematica syntax.  The squares are functions of D (through the symbol D of the declarations); the unpolarized results of the channels with an incoming gluon or photon still need the averages 1/(D-2) per vector boson (FORM does not divide by polynomials).*)


(* ::Input:: *)
born = <|"Unpolarized" -> ReadFORM[work, "BornU.m"], "Polarized" -> ReadFORM[work, "BornP.m"]|>;
Require[Simplify[(born["Unpolarized"] - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];
Require[Simplify[(born["Polarized"] - 2 ee^4 eq^2 (s^2 - u^2)/t^2) /. u -> -s - t] === 0, "Polarized Born square differs from the expected result"];
born


(* ::Input:: *)
real = <|"QG" -> <|"Unpolarized" -> ReadFORM[work, "RealQGU.m"], "Polarized" -> ReadFORM[work, "RealQGP.m"]|>,
  "GQ" -> <|"Unpolarized" -> ReadFORM[work, "RealGQU.m"]/(D - 2), "Polarized" -> ReadFORM[work, "RealGQP.m"]|>|>;
photon = <|"QG" -> <|"Unpolarized" -> ReadFORM[work, "PhotonQGU.m"]/(D - 2), "Polarized" -> ReadFORM[work, "PhotonQGP.m"]|>,
  "GQ" -> <|"Unpolarized" -> ReadFORM[work, "PhotonGQU.m"]/(D - 2)^2, "Polarized" -> ReadFORM[work, "PhotonGQP.m"]|>|>;
Require[FreeQ[photon, den], "A photon propagator was not recognized by the FORM programs", Cases[photon, _den, Infinity]];
Factor[real["GQ"]["Polarized"] /. spe3 -> 0]


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Photon" -> photon|>;
SaveResult["Algebra", algebra, "FORM"];
Keys[algebra]
