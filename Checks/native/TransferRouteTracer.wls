(* Compare the transverse spin-transfer q g -> q class through Tracer
   and FeynCalc. This is the final comparison cell of
   Processes/HadronicTransfer/2-Algebra-Tracer.wl, with the same mathematics.

   Set FEYNFACET_SOURCE to the repository containing computed Results/ files.
   Set FEYNFACET_BENCHMARK_OUTPUT to a NEW local output directory.
   Run with wolframscript -file Checks/native/TransferRouteTracer.wls.
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
sourceFile=FileNameJoin[{processDirectory,"2-Algebra-Tracer.wl"}];
routeFile=FileNameJoin[{processDirectory,"Results","Tracer","class-13","Algebra.wl"}];
feynCalcFile=FileNameJoin[{processDirectory,"Results","FeynCalc","class-13","Algebra.wl"}];
Require[And@@(FileExistsQ/@{sourceFile,routeFile,feynCalcFile}),
  "Run the FeynCalc and Tracer class-13 algebra notebooks first."];
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



SeedRandom[2026101702];
nativeResult=MemoryConstrained[TimeConstrained[CheckAbort[(
(* NATIVE-BEGIN *)
feynCalcFile = FileNameJoin[{ResultsDirectory["FeynCalc/class-" <> class], "Algebra.wl"}];
comparison = If[!FileExistsQ[feynCalcFile], "FeynCalc result not available: comparison skipped",
  Module[{other = Get[feynCalcFile], integrals, points, valueAt},
    integrals = Union[Cases[{virtual, other["Virtual"]}, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]];
    points = Table[Join[With[{v = RandomInteger[{1, 997}, 9]/RandomInteger[{998, 1997}, 9]},
        {s -> 1 + v[[1]], t -> -v[[2]], u -> -v[[3]], A -> -v[[4]], B -> -v[[5]], D -> 4 - v[[6]], CA -> 2 + 10 v[[7]], spe3 -> v[[8]], Nf -> 1 + 5 v[[9]], hp -> 1 + v[[1]] v[[2]]}],
      Thread[integrals -> RandomInteger[{1, 10^6}, Length[integrals]]/RandomInteger[{1, 10^6}, Length[integrals]]]], {2}];
    Put[points, pointsFile];
    valueAt[expr_, point_] := expr /. HatPower[n_] :> hp^n /. point /. C -> -(s + t + u) - A - B /. CF -> (CA^2 - 1)/(2 CA) /. point;
    Join[
      Association @ Flatten[Table[part <> " " <> sc -> Union[Table[Together[valueAt[algebra[part][sc], point] - valueAt[other[part][sc], point]], {point, points}]],
        {part, {"Born", "Virtual", "Real"}}, {sc, outputSpins}]],
      Association @ Table["BornTable " <> sc -> Union[Flatten[Table[Together[valueAt[bornTable[sc][name], point] - valueAt[other["BornTable"][sc][name], point]],
        {name, Keys[other["BornTable"][sc]]}, {point, points}]]], {sc, outputSpins}]]]];
Require[!AssociationQ[comparison] || Union[Flatten[Values[comparison]]] === {0}, "The Tracer route differs from the FeynCalc route", comparison];
comparison
(* NATIVE-END *)
; comparison),$Aborted],180,$Aborted],4*1024^3,$Aborted];
Require[AssociationQ[nativeResult],"The native comparison was skipped, aborted or failed."];

requiredKeys=Join[Flatten[Table[part<>" "<>component,
 {part,{"Born","Virtual","Real"}},{component,outputSpins}]],("BornTable "<>#&/@outputSpins)];
Require[Sort[Keys[nativeResult]]===Sort[requiredKeys] && And@@(#==={0}&/@Values[nativeResult]),
 "The six tensor and two seven-channel Born-table comparisons are incomplete or nonzero."];

Require[FileExistsQ[pointsFile],"The exact sample points were not saved."];
exactPoints=Get[pointsFile];
Require[ListQ[exactPoints] && Length[exactPoints]===2,"Expected exactly two rational sample points."];
Require[And@@KeyValueMap[FileExistsQ[#1]&&sha[#1]===#2&,inputHashes],
 "An input changed during the benchmark."];
jsonResult=KeyMap[If[StringQ[#],#,ToString[#,InputForm]]&,nativeResult];
result=<|"Passed"->True,"Route"->"Tracer","ReferenceRoute"->"FeynCalc","Class"->"13",
 "Seed"->2026101702,"ExactRationalPoints"->2,"RequiredTensorComponents"->{"XX","NN"},
 "RequiredBornChannels"->expectedChannels,"RequiredScalarComparisons"->40,
 "NativeVerdicts"->jsonResult,"NativeCellHeldASTMatchesSource"->True,
 "NativeCellASTSHA256"->nativeCellASTSHA256,"InputSHA256"->inputHashes,
 "ExactPointsSHA256"->sha[pointsFile],"WolframVersion"->$Version,"MaximumMemoryBytes"->MaxMemoryUsed[]|>;
Export[FileNameJoin[{outputRoot,"comparison.json"}],result,"RawJSON"];
Print["Tracer/FeynCalc: 40 exact sampled scalar comparisons passed; exact points saved."];
Quit[0];
