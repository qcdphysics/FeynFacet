(* Compare the transverse spin-transfer q g -> q class through FORM
   and FeynCalc. This is the final comparison cell of
   Processes/HadronicTransfer/2-Algebra-FORM.wl, with the same mathematics.

   Set FEYNFACET_SOURCE to the repository containing computed Results/ files.
   Set FEYNFACET_BENCHMARK_OUTPUT to a NEW local output directory.
   Run with wolframscript -file Checks/native/TransferRouteFORM.wls.
   Two exact rational points are saved; these are sampled algebra checks,
   rather than a claim of a symbolic identity or independent integration. *)
$HistoryLength = 0;
Require[condition_,message_String,details_:None] := If[TrueQ[condition],True,
  Print[message,If[details===None,"",": "<>ToString[details,InputForm]]];Quit[1]];
sourceRoot=Environment["FEYNFACET_SOURCE"];
outputRoot=Environment["FEYNFACET_BENCHMARK_OUTPUT"];
Require[And@@(StringQ[#]&&StringLength[StringTrim[#]]>0&/@{sourceRoot,outputRoot}),
  "Set FEYNFACET_SOURCE and FEYNFACET_BENCHMARK_OUTPUT."];
sourceRoot=ExpandFileName[sourceRoot]; outputRoot=ExpandFileName[outputRoot];
Require[DirectoryQ[sourceRoot],"The source directory does not exist."];
Require[!FileExistsQ[outputRoot],"Choose a new benchmark output directory."];
processDirectory=FileNameJoin[{sourceRoot,"Processes","HadronicTransfer"}];
sourceFile=FileNameJoin[{processDirectory,"2-Algebra-FORM.wl"}];
routeFile=FileNameJoin[{processDirectory,"Results","FORM","class-13","Algebra.wl"}];
feynCalcFile=FileNameJoin[{processDirectory,"Results","FeynCalc","class-13","Algebra.wl"}];
Require[And@@(FileExistsQ/@{sourceFile,routeFile,feynCalcFile}),
  "Run the FeynCalc and FORM class-13 algebra notebooks first."];
sha[file_]:=IntegerString[FileHash[file,"SHA256"],16,64];
scriptFile=ExpandFileName[$InputFileName];
inputHashes=Association[#->sha[#]&/@{sourceFile,routeFile,feynCalcFile,scriptFile}];
CreateDirectory[outputRoot,CreateIntermediateDirectories->True];
pointsFile=FileNameJoin[{outputRoot,"exact-points.wl"}];
expectedChannels={"q qp -> q","q qpb -> q","q q -> q","q qb -> q",
  "q qb -> qb","q g -> q","qb q -> q"};
algebra=Get[routeFile]; canonical=Get[feynCalcFile];
Require[And@@(AssociationQ/@{algebra,canonical}),"An algebra file is invalid."];
Do[
 Require[item["Class"]==="13" && item["Spins"]==={"XX","NN"},
  "Expected the class-13 transverse spin-transfer tensor."];
 Require[Sort[Keys[item["BornTable"]]]===Sort[{"XX","NN"}] &&
  And@@Table[Sort[Keys[item["BornTable"][component]]]===Sort[expectedChannels],{component,{"XX","NN"}}],
  "The seven transfer Born channels are incomplete."],{item,{algebra,canonical}}];
Clear[canonical];
bornTable=algebra["BornTable"]; virtual=algebra["Virtual"];
families=algebra["Families"]; outputSpins=algebra["Spins"]; class="13";
ResultsDirectory[route_String]:=FileNameJoin[{processDirectory,"Results",route}];

(* Check that the copied native mathematics still matches this source version.
   Only the point-saving statement is removed before the held AST comparison. *)
inputCells[text_]:=Module[{pieces,style=None,cells={}},
 pieces=StringSplit[text,RegularExpression["(?m)^\\(\\* ::(\\w+):: \\*\\)\\s*$"]->"$1"];
 Do[If[StringMatchQ[piece,RegularExpression["^\\w+$"]],style=piece,
   If[style==="Input" && StringTrim[piece]=!="",AppendTo[cells,StringTrim[piece]]]],{piece,pieces}];cells];
sourceCells=inputCells[Import[sourceFile,"Text"]];
Require[sourceCells=!={},"The source has no executable Input cells."];
scriptParts=StringSplit[Import[scriptFile,"Text"],{
 RegularExpression["(?m)^\\(\\* NATIVE-BEGIN \\*\\)\\s*$"],
 RegularExpression["(?m)^\\(\\* NATIVE-END \\*\\)\\s*$"]}];
Require[Length[scriptParts]===3,"The benchmark native-code markers are invalid."];
captureLine="    Put[points, pointsFile];\n";
Require[StringCount[scriptParts[[2]],captureLine]===1,"The point recorder differs from its specified form."];
sourceHeld=ToExpression[Last[sourceCells],InputForm,HoldComplete];
benchmarkNative=StringTrim[StringReplace[scriptParts[[2]],captureLine->""]];
Require[SameQ[Last[sourceCells],benchmarkNative],
 "The benchmark comparison text differs from the source notebook."];
benchmarkHeld=ToExpression[benchmarkNative,InputForm,HoldComplete];
Require[MatchQ[sourceHeld,_HoldComplete] && SameQ[sourceHeld,benchmarkHeld],
 "The benchmark comparison block differs from the source notebook."];
nativeCellASTSHA256=IntegerString[Hash[sourceHeld,"SHA256"],16,64];

bornProducts = <|{pa, pa} -> 0, {pb, pb} -> 0, {p, p} -> 0, {k1, k1} -> 0, {pa, pb} -> s/2, {pa, p} -> -t/2, {pb, p} -> -u/2,
  {pa, k1} -> -u/2, {pb, k1} -> -t/2, {p, k1} -> s/2|>;
dot[q1_, q2_, table_] := Module[{v = {pa, pb, p, k1, k2}},
  Expand[Sum[Coefficient[q1, a] Coefficient[q2, b] Lookup[table, Key[Sort[{a, b}]], Lookup[table, Key[{a, b}], Lookup[table, Key[{b, a}], 0]]], {a, v}, {b, v}]]];
gramMatrix[table_] := Table[dot[x, y, table], {x, {pa, pb, p}}, {y, {pa, pb, p}}];
parallelSquare[d_List, gram_] := Factor[d . Inverse[gram] . d];
inverseProps[name_] := Module[{q = families[name], sol},
  sol = Solve[Table[Dj[j] == L2 + 2 (Coefficient[q[[j]], pa] Lpa + Coefficient[q[[j]], pb] Lpb + Coefficient[q[[j]], p] Lp) + dot[q[[j]] - loop, q[[j]] - loop, bornProducts] /. u -> -s - t, {j, 4}], {L2, Lpa, Lpb, Lp}];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators"]; First[sol]];
loopRules = Association @ Table[name -> inverseProps[name], {name, Keys[families]}];
loopPerp2 = L2 - parallelSquare[{Lpa, Lpb, Lp}, gramMatrix[bornProducts] /. u -> -s - t];
perpMoment[a_, b_] := (-1)^a Pochhammer[1/2, a] Pochhammer[(D - 4)/2, b]/Pochhammer[(D - 3)/2, a + b];

SeedRandom[2026101701];
nativeResult=MemoryConstrained[TimeConstrained[CheckAbort[(
(* NATIVE-BEGIN *)
feynCalcFile = FileNameJoin[{ResultsDirectory["FeynCalc/class-" <> class], "Algebra.wl"}];
routeComparison = If[!FileExistsQ[feynCalcFile], "the FeynCalc route has not been run for this class",
  Module[{other = Get[feynCalcFile], integrals, integralValues, points, evaluate, parts},
    integrals = Union[Cases[{algebra["Virtual"], other["Virtual"]}, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]];
    integralValues = Table[With[{v = RandomInteger[{1, 10^9}, Length[integrals]]/RandomInteger[{1, 10^9}, Length[integrals]]},
      Require[DuplicateFreeQ[v], "The random values of the integrals are not distinct"]; v], {2}];
    points = Table[Join[Thread[{s, t, u, A, B, D, CA, spe3, Nf, gs, hp} -> RandomInteger[{1, 10^6}, 11]/RandomInteger[{1, 10^6}, 11] {1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1}],
      Thread[integrals -> integralValues[[k]]]], {k, 2}];
    Put[points, pointsFile];
    evaluate[x_, point_] := Together[x /. HatPower[n_] :> hp^n /. point /. C -> -(s + t + u) - A - B /. CF -> (CA^2 - 1)/(2 CA) /. point];
    (* the hat projection of a virtual interference, at the kinematic part of a point (the integrals stay symbols until the shifts are made) *)
    projectHat[x_, point_] := Module[{kin = Select[point, FreeQ[#[[1]], KVirtualST | KVirtualSU | KVirtualTU] &], expr},
      expr = Expand[x /. CF -> (CA^2 - 1)/(2 CA) /. kin];
      Expand[expr /. HatPower[b_] h_[n__] /; b > 0 && MemberQ[{KVirtualST, KVirtualSU, KVirtualTU}, h] :>
        With[{name = StringDrop[SymbolName[h], 1]}, (perpMoment[0, b] /. kin) Total[Map[(#[[2]] HatPower[0] h @@ ({n} - #[[1]])) &,
          CoefficientRules[Expand[(loopPerp2 /. loopRules[name] /. kin)^b], {Dj[1], Dj[2], Dj[3], Dj[4]}]]]]]];
    evaluateProjected[x_, point_] := Together[projectHat[x, point] /. HatPower[0] -> 1 /. point];
    Require[Keys[other["BornTable"]["XX"]] === Keys[algebra["BornTable"]["XX"]], "The two routes have different Born tables"];
    parts = Join[Flatten[Table[{part, spinCase} -> {algebra[part][spinCase], other[part][spinCase]}, {part, {"Born", "Virtual", "Real"}}, {spinCase, outputSpins}]],
      Flatten[Table[{"BornTable", spinCase, name} -> {algebra["BornTable"][spinCase][name], other["BornTable"][spinCase][name]}, {spinCase, outputSpins}, {name, Keys[algebra["BornTable"]["XX"]]}]]];
    Join[Association @ Table[item[[1]] -> (And @@ Table[If[item[[1, 1]] === "Virtual", evaluateProjected, evaluate][item[[2, 1]] - item[[2, 2]], point] === 0, {point, points}]), {item, parts}],
      <|"VirtualLiteral" -> Association @ Table[spinCase -> (And @@ Table[evaluate[algebra["Virtual"][spinCase] - other["Virtual"][spinCase], point] === 0, {point, points}]), {spinCase, outputSpins}]|>]]];
Require[StringQ[routeComparison] || And @@ Values[KeyDrop[routeComparison, "VirtualLiteral"]], "The FORM route differs from the FeynCalc route", Select[KeyDrop[routeComparison, "VirtualLiteral"], ! # &]];
If[StringQ[routeComparison], routeComparison, {Counts[Values[KeyDrop[routeComparison, "VirtualLiteral"]]], "literal virtual comparison" -> routeComparison["VirtualLiteral"]}]
(* NATIVE-END *)
; routeComparison),$Aborted],180,$Aborted],4*1024^3,$Aborted];
Require[AssociationQ[nativeResult],"The native comparison was skipped, aborted or failed."];

requiredKeys=Join[
 Flatten[Table[{part,component},{part,{"Born","Virtual","Real"}},{component,outputSpins}],1],
 Flatten[Table[{"BornTable",component,name},{component,outputSpins},{name,expectedChannels}],1]];
Require[Sort[Keys[KeyDrop[nativeResult,"VirtualLiteral"]]]===Sort[requiredKeys] &&
 And@@(TrueQ/@Values[KeyDrop[nativeResult,"VirtualLiteral"]]),
 "The 20 required FORM tensor/table verdicts are incomplete or nonzero."];
Require[KeyExistsQ[nativeResult,"VirtualLiteral"] && Sort[Keys[nativeResult["VirtualLiteral"]]]===Sort[outputSpins],
 "The literal virtual verdicts were not retained."];
(* The NN literal verdict can be False. The required virtual comparison uses
   the angular projection in the unchanged native block; both are reported. *)

Require[FileExistsQ[pointsFile],"The exact sample points were not saved."];
exactPoints=Get[pointsFile];
Require[ListQ[exactPoints] && Length[exactPoints]===2,"Expected exactly two rational sample points."];
Require[And@@KeyValueMap[FileExistsQ[#1]&&sha[#1]===#2&,inputHashes],
 "An input changed during the benchmark."];
jsonResult=KeyMap[If[StringQ[#],#,ToString[#,InputForm]]&,nativeResult];
result=<|"Passed"->True,"Route"->"FORM","ReferenceRoute"->"FeynCalc","Class"->"13",
 "Seed"->2026101701,"ExactRationalPoints"->2,"RequiredTensorComponents"->{"XX","NN"},
 "RequiredBornChannels"->expectedChannels,"RequiredScalarComparisons"->40,
 "NativeVerdicts"->jsonResult,"NativeCellHeldASTMatchesSource"->True,
 "NativeCellASTSHA256"->nativeCellASTSHA256,"InputSHA256"->inputHashes,
 "ExactPointsSHA256"->sha[pointsFile],"WolframVersion"->$Version,"MaximumMemoryBytes"->MaxMemoryUsed[]|>;
Export[FileNameJoin[{outputRoot,"comparison.json"}],result,"RawJSON"];
Print["FORM/FeynCalc: 40 exact sampled scalar comparisons passed; exact points saved."];
Quit[0];
