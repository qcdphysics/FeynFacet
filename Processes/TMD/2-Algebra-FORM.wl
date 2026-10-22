(* ::Title:: *)
(*TMD beam and soft functions, step 2 (FORM route): squared matrix elements*)


(* ::Text:: *)
(*Same task and same output format as the FeynCalc route (read its introduction for the definitions), now with FORM.  Input: the amplitude records of step 1.  Output: Results/FORM/Algebra.wl.  The amplitudes are written to a FORM include file; the program form/Square.frm squares the amplitudes of one channel with the normalization given on its command line, and its header form/common.h contains the conjugation, the closing of the operator index, the polarization and colour sums.  The notebook only collects the results and performs the same checks as the FeynCalc route.*)


(* ::Section:: *)
(*Load and export the amplitudes*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
records = LoadAmplitudes[]["Records"];
workDir = FileNameJoin[{ResultsDirectory["FORM"], "form"}];
(* one FORM expression per diagram: BeamQQ1, BeamQQ2, ... *)
labelled = Association @ Flatten[KeyValueMap[Function[{name, list}, Table[name <> ToString[i] -> list[[i]], {i, Length[list]}]], records]];
WriteFORMAmplitudes[labelled, workDir];
Keys[labelled]


(* ::Input:: *)
Import[FileNameJoin[{workDir, "amplitudes.inc"}], "Text"]


(* ::Section:: *)
(*Square every channel*)


(* ::Text:: *)
(*The normalizations are those of the FeynCalc route: quark operator 1/(2 n.p) with the averages 1/2 and 1/Nc (incoming quark) or 1/(D-2) and 1/(Nc^2-1) (incoming gluon); gluon operator 1/(x (n.p)^2) with the same averages; soft function 1/Nc or 1/(Nc^2-1) for the colour trace.  The operator markers gOp, gSoft are divided out.*)


(* ::Input:: *)
program = FileNameJoin[{NotebookDirectory[], "form", "Square.frm"}];
channels = <|
  "BornQ"  -> {"BornQ1", "1/(2*2*Nc*pn*gOp^2)", "beam"},
  "BornG"  -> {"BornG1", "1/((D-2)*(Nc^2-1)*pn^2*gOp^2)", "beam"},
  "BeamQQ" -> {"BeamQQ1+BeamQQ2", "1/(2*2*Nc*pn*gOp^2)", "beam"},
  "BeamQG" -> {"BeamQG1+BeamQG2", "1/((D-2)*(Nc^2-1)*2*pn*gOp^2)", "beam"},
  "BeamGQ" -> {"BeamGQ1", "1/(2*Nc*x*pn^2*gOp^2)", "beam"},
  "BeamGG" -> {"BeamGG1+BeamGG2+BeamGG3+BeamGG4", "1/((D-2)*(Nc^2-1)*x*pn^2*gOp^2)", "beam"},
  "SoftQ"  -> {"SoftQ1+SoftQ2", "1/(Nc*gSoft^2)", "soft"},
  "SoftG"  -> {"SoftG1+SoftG2", "1/((Nc^2-1)*gSoft^2)", "soft"}|>;
Do[RunFORM[program, workDir, "Define" -> {"AMP=" <> channels[name][[1]], "NORM=" <> channels[name][[2]], "KIN=" <> channels[name][[3]], "OUT=" <> name <> ".m"}], {name, Keys[channels]}];
squares = AssociationMap[Factor[ReadFORM[workDir, # <> ".m"]] &, Keys[channels]]


(* ::Text:: *)
(*FORM leaves the invariants as symbols: pn = n.p, kn = n.k, pk = p.k, kb = nb.k.  For the beam channels n.k = (1-x) n.p, p.k = kT2/(2(1-x)) and nb.k = kT2/((1-x) n.p); the soft matrix element keeps kn and kb.*)


(* ::Input:: *)
beamKinematics = {kn -> (1 - x) pn, pk -> kT2/(2 (1 - x)), kb -> kT2/((1 - x) pn)};
squares = Association @ KeyValueMap[#1 -> If[channels[#1][[3]] === "beam", Factor[#2 /. beamKinematics], #2] &, squares]


(* ::Section:: *)
(*Checks and output*)


(* ::Text:: *)
(*Colour factors through CF and CA (FORM's results contain Nc, CF and CA); then the same checks as in the FeynCalc route: Born = 1, the real-emission kernels are the D-dimensional splitting functions, the soft matrix element is 4 gs^2 C_i/((n.k)(nb.k)).*)


(* ::Input:: *)
colourFactors[expr_] := Factor[expr /. CF -> (Nc^2 - 1)/(2 Nc) /. CA -> Nc] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
born = <|"q" -> colourFactors[squares["BornQ"]], "g" -> colourFactors[squares["BornG"]]|>;
Require[born === <|"q" -> 1, "g" -> 1|>, "The Born matrix elements are not normalized to 1", born];
real = colourFactors /@ <|"qq" -> squares["BeamQQ"], "qg" -> squares["BeamQG"], "gq" -> squares["BeamGQ"], "gg" -> squares["BeamGG"]|>;
soft = colourFactors /@ <|"q" -> squares["SoftQ"], "g" -> squares["SoftG"]|>;
Require[FreeQ[real, pn], "The real matrix elements depend on n.p", real];
kernels = Map[Factor[# kT2/gs^2] &, real];
expected = <|"qq" -> 2 CF ((1 + x^2)/(1 - x) - eps (1 - x)), "qg" -> 2 (1/2) (1 - 2 x (1 - x)/(1 - eps)),
  "gq" -> 2 CF ((1 + (1 - x)^2)/x - eps x), "gg" -> 4 CA (x/(1 - x) + (1 - x)/x + x (1 - x))|>;
Require[And @@ Table[Simplify[kernels[c]/(1 - x) /. D -> 4 - 2 eps] === Simplify[expected[c]], {c, Keys[kernels]}],
  "The real-emission kernels differ from the D-dimensional splitting functions", kernels];
Require[Simplify[soft["q"] - 4 gs^2 CF/(kn kb)] === 0 && Simplify[soft["g"] - 4 gs^2 CA/(kn kb)] === 0, "Soft matrix elements differ from 4 gs^2 C_i/((n.k)(nb.k))", soft];
TableForm[Normal[kernels], TableHeadings -> {None, {"channel", "kT2 M2/gs^2"}}]


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Soft" -> soft|>;
SaveResult["Algebra", algebra, "FORM"];
Keys[algebra]
