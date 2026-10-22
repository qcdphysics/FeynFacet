(* ::Title:: *)
(*Lepton-quark scattering with transverse spin transfer, step 2 (FORM route): spin and colour algebra*)


(* ::Text:: *)
(*Same task and same output format as the FeynCalc route, now with FORM doing the Dirac algebra (the programs are in form/; read Processes/LeptonQuark/2-Algebra-FORM.wl and its form/common.h first).  Born.frm, Real.frm and Photon.frm are each run three times, with SPIN = U (spin average), NN and XX (transverse spin transfer with the normal and the in-plane axes); Virtual.frm projects the vertex graph onto the Born current.  The in-plane axes are used with polynomial coefficients, aXs = s aX and bXs = (t+u)^2 bX (FORM has no polynomial denominators), so the XX results are divided here by s (t+u)^2 and by the product of the lengths -t u/(t+u), i.e. by -s t u (t+u); the normal component nr of the real emission is reduced here as in the FeynCalc route.*)


(* ::Section:: *)
(*Load and export*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
records = LoadAmplitudes[]["Records"];
work = ResultsDirectory["FORM"];
WriteFORMAmplitudes[<|"Born1" -> records["Born"][[1]], "Real1" -> records["Real"][[1]], "Real2" -> records["Real"][[2]],
  "Vertex1" -> records["VirtualVertex"][[1]], "Photon1" -> records["Photon"][[1]], "Photon2" -> records["Photon"][[2]]|>, work]


(* ::Section:: *)
(*Run the FORM programs*)


(* ::Input:: *)
program[name_] := FileNameJoin[{NotebookDirectory[], "form", name}];
Do[RunFORM[program[name], work, "Define" -> {"SPIN=" <> spin}], {name, {"Born.frm", "Real.frm", "Photon.frm"}}, {spin, {"U", "NN", "XX"}}];
RunFORM[program["Virtual.frm"], work, "Define" -> {"SPIN=U"}];
FileNames["*.m", work]


(* ::Section:: *)
(*Read the results*)


(* ::Input:: *)
spinCases = {"Unpolarized", "NN", "XX"};
formName = <|"Unpolarized" -> "U", "NN" -> "NN", "XX" -> "XX"|>;
xxScale = <|"Unpolarized" -> 1, "NN" -> 1, "XX" -> -s t u (t + u)|>;
xxScalePhoton = <|"Unpolarized" -> 1, "NN" -> 1, "XX" -> -sh th uh (th + uh)|>;
born = Association @ Table[spin -> Factor[ReadFORM[work, "Born" <> formName[spin] <> ".m"]/(xxScale[spin] /. u -> -s - t)], {spin, spinCases}];
born["Transverse"] = Factor[(born["NN"] + born["XX"])/2];
Require[Simplify[(born["Unpolarized"] - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];
Require[Simplify[{born["NN"] - ee^4 eq^2 (-4 s u/t^2 - (D - 4)), born["XX"] - ee^4 eq^2 (-4 s u/t^2 + (D - 4))} /. u -> -s - t] === {0, 0},
  "Transverse Born tensors differ from Eq. (transverse-born-tensor) of the paper"];
born


(* ::Input:: *)
gram = {{0, s/2, -u/2}, {s/2, 0, -t/2}, {-u/2, -t/2, 0}};                     (* l.l, l.k, l.p; ... *)
parallelSquare = Factor[{(s + t + u + a + b)/2, -a/2, b/2}.Inverse[gram].{(s + t + u + a + b)/2, -a/2, b/2}];
normalReduction[expr_] := Module[{coefficients = CoefficientList[Expand[expr], nr]},
  Require[Length[coefficients] <= 3, "The real-emission tensor has a power of nr above two", Length[coefficients] - 1];
  Factor[coefficients[[1]] + If[Length[coefficients] == 3, coefficients[[3]] (parallelSquare + spe3), 0]]];
real = Association @ Table[spin -> normalReduction[ReadFORM[work, "Real" <> formName[spin] <> ".m"]/xxScale[spin]], {spin, spinCases}];
real["Transverse"] = Factor[(real["NN"] + real["XX"])/2];
Require[Factor[(real["NN"] - real["XX"]) /. {D -> 4, spe3 -> 0}] === 0, "The two transverse components differ in four dimensions"];
Factor[real["Transverse"] /. {D -> 4, spe3 -> 0}]


(* ::Input:: *)
photon = Association @ Table[spin -> Factor[ReadFORM[work, "Photon" <> formName[spin] <> ".m"]/(D - 2)/(xxScalePhoton[spin] /. uh -> -sh - th)], {spin, spinCases}];
photon["Transverse"] = Factor[(photon["NN"] + photon["XX"])/2];
Require[Simplify[(photon["Transverse"]/(CF ee^2 eq^2 gs^2) - (4 - 2 eps^2 th^2/((1 - eps) sh uh))) /. D -> 4 - 2 eps /. uh -> -sh - th] === 0,
  "The photon transverse tensor differs from Eq. (transverse-photon-born) of the paper"];
photon


(* ::Input:: *)
{numerator, bornCurrent} = ReadFORM[work, "Virtual.m"];
virtualFamily = Factor[numerator/bornCurrent]


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Photon" -> photon, "VirtualFamily" -> virtualFamily|>;
SaveResult["Algebra", algebra, "FORM"];
Keys[algebra]
