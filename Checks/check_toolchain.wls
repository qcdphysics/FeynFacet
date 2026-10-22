(* Functional checks of the toolchain, one test per fresh kernel.
   Usage: wolframscript -file Checks/check_toolchain.wls <test> <scratch directory>
   Tests: kernel, tools, feyncalc, form, tracer, kira, polymake, subtropica.
   Each run prints one line "FEYNFACET_CHECK {json}" with the outcome; feynfacet.py check runs them all.
   Scratch files are written only below the given directory. *)
$HistoryLength = 0;
Get[FileNameJoin[{DirectoryName[$InputFileName], "ScriptArguments.wl"}]];
{checkName, checkScratch} = Take[Join[ScriptArguments[], {"", ""}], 2];
emit[status_String, detail_] := (
  Print["FEYNFACET_CHECK " <> ExportString[<|"test" -> checkName, "status" -> status,
    "detail" -> detail|>, "RawJSON", "Compact" -> True]];
  Exit[If[status === "pass", 0, 1]]);
(* Require prints a message and aborts; turn any abort or message-producing failure into a report *)
SetAttributes[guarded, HoldFirst];
guarded[expr_] := Module[{out, printed = {}},
  out = Block[{Print = (AppendTo[printed, StringJoin[ToString /@ {##}]] &)}, CheckAbort[Quiet[expr], $Aborted]];
  If[out === $Aborted, emit["fail", StringRiffle[Take[printed, -Min[5, Length[printed]]], " | "]], out]];

root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
If[checkName =!= "kernel", guarded[Get[FileNameJoin[{root, "FeynFacet", "FeynFacet.wl"}]]]];
If[StringQ[checkScratch] && checkScratch =!= "",
  Quiet[CreateDirectory[checkScratch, CreateIntermediateDirectories -> True]];
  SetEnvironment["FEYNFACET_PROCESS_DIRECTORY" -> checkScratch]];

Switch[checkName,

 "kernel",
 emit["pass", <|"version" -> $Version, "temporaryDirectory" -> $TemporaryDirectory,
   "userBaseDirectory" -> $UserBaseDirectory|>],

 "tools",
 Module[{locate, rows},
  locate[expr_] := Block[{Print}, Quiet[CheckAbort[Check[expr, Missing[]], Missing[]]]];
  SetAttributes[locate, HoldFirst];
  rows = <|
    "FeynCalc" -> locate[FindFile["FeynCalc`"]],
    "SubTropica" -> locate[FindFile["SubTropica`"]],
    "Tracer" -> locate[ToolPath["Tracer"]],
    "FORM" -> locate[ToolPath["FORM"]], "TFORM" -> locate[ToolPath["TFORM"]],
    "FORMColor" -> locate[ToolPath["FORMColor"]],
    "Kira" -> locate[ToolPath["Kira"]], "Fermat" -> locate[ToolPath["Fermat"]],
    "Polymake" -> locate[ToolPath["Polymake"]]|>;
  emit["pass", Map[If[StringQ[#], #, Null] &, rows]]],

 "feyncalc",
 guarded[LoadFeynCalc[]];
 Module[{trace, topologies},
  trace = guarded[FeynCalc`DiracSimplify[FeynCalc`DiracTrace[FeynCalc`GAD[al, be, al, be]],
    FeynCalc`DiracTraceEvaluate -> True]];
  topologies = guarded[Length[FeynArts`CreateTopologies[0, 2 -> 2]]];
  If[Simplify[trace - 4 D (2 - D)] === 0 && topologies == 4,
   emit["pass", <|"trace" -> ToString[trace, InputForm], "treeTopologies" -> topologies|>],
   emit["fail", <|"trace" -> ToString[trace, InputForm], "treeTopologies" -> topologies,
     "expected" -> "trace 4 D (2 - D), 4 tree topologies"|>]]],

 "form",
 Module[{dir = FileNameJoin[{checkScratch, "form"}], program, value},
  Quiet[CreateDirectory[dir, CreateIntermediateDirectories -> True]];
  program = FileNameJoin[{dir, "check.frm"}];
  Export[program, StringRiffle[{
    "#include- color.h", "Dimension D;", "Indices al, be;", "Local F = g_(1, al, be, al, be);",
    "tracen, 1;", ".sort", "Format Mathematica;", "#write <check.m> \"%E\", F", ".end"}, "\n"], "Text"];
  guarded[RunFORM[program, dir]];
  value = guarded[ToExpression[StringReplace[Import[FileNameJoin[{dir, "check.m"}], "Text"], WhitespaceCharacter -> ""]]];
  If[Simplify[value - 4 D (2 - D)] === 0,
   emit["pass", <|"trace" -> ToString[value, InputForm]|>],
   emit["fail", <|"trace" -> ToString[value, InputForm], "expected" -> "4 D (2 - D)", "log" -> dir|>]]],

 "tracer",
 guarded[LoadTracer[]];
 Module[{value},
  (* Tracer's symbols exist only after loading: look them up at run time *)
  ToExpression["VectorDimension[D]; Spur[1]"];
  value = guarded[Expand[ToExpression["G[1, {al}, {be}, {al}, {be}]"]]];
  If[Simplify[value - 4 D (2 - D)] === 0,
   emit["pass", <|"trace" -> ToString[value, InputForm]|>],
   emit["fail", <|"trace" -> ToString[value, InputForm], "expected" -> "4 D (2 - D)"|>]]],

 "kira",
 Module[{fixtures = FileNameJoin[{root, "Checks", "kira"}], specs, work, targets,
   rules, expected, differences, results = {}, executable},
  executable = guarded[ToolPath["Kira"]];
  (* Printed before execution so the selected binary is visible even if Kira aborts. *)
  Print["Checking Kira executable: " <> executable];
  specs = Import[FileNameJoin[{fixtures, "expectations.json"}], "RawJSON"];
  Do[
   work = FileNameJoin[{checkScratch, spec["directory"]}];
   targets = ToExpression /@ Lookup[spec["cases"], "target"];
   expected = (ToExpression /@ Lookup[spec["cases"], "coefficient"]) ToExpression[spec["master"]];
   guarded[PrepareKira[FileNameJoin[{fixtures, spec["directory"]}], work]];
   guarded[WriteKiraTargets[work, targets, {}]];
   guarded[RunKira[work]];
   rules = guarded[KiraRules[work, spec["family"], {}]];
   (* Exact rational identities, including explicit zeros; a right master with
      a wrong coefficient, missing rule, or unknown master must fail. *)
   differences = Together /@ ((targets /. rules) - expected);
   If[!ListQ[rules] || !And @@ (SameQ[#, 0] & /@ differences),
    emit["fail", <|"executable" -> executable, "family" -> spec["family"],
      "differences" -> ToString[differences, InputForm],
      "log" -> FileNameJoin[{work, "kira.log"}]|>]];
   AppendTo[results, <|"family" -> spec["family"], "exactIdentities" -> Length[targets]|>],
   {spec, specs}];
  emit["pass", <|"executable" -> executable, "reductions" -> results|>]],

 "polymake",
 (* a native calculation in the two polymake applications that SubTropica uses: a version
    banner alone does not show that polymake's Perl modules load *)
 Module[{dir = FileNameJoin[{checkScratch, "polymake"}], polymake, script, run},
  Quiet[CreateDirectory[dir, CreateIntermediateDirectories -> True]];
  polymake = guarded[ToolPath["Polymake"]];
  script = FileNameJoin[{dir, "check.pl"}];
  Export[script, StringRiffle[{"use application \"polytope\";", "my $c = cube(3);",
    "print $c->N_FACETS, \"\\n\";", "use application \"fan\";",
    "print normal_fan($c)->N_RAYS, \"\\n\";"}, "\n"] <> "\n", "Text"];
  run = RunProcess[{polymake, "--script", script}];
  If[run["ExitCode"] === 0 && StringSplit[run["StandardOutput"]] === {"6", "6"} &&
     !StringContainsQ[run["StandardError"], "ERROR"],
   emit["pass", <|"executable" -> polymake, "cube" -> "6 facets, normal fan with 6 rays"|>],
   emit["fail", <|"executable" -> polymake, "exitCode" -> run["ExitCode"],
     "output" -> StringTake[run["StandardOutput"], UpTo[400]],
     "errors" -> StringTake[run["StandardError"], UpTo[800]]|>]]],

 "subtropica",
 guarded[LoadSubTropica[]];
 Module[{value, difference, box, poles, calls},
  value = guarded[EulerIntegral[1/Gamma[1 - eps], x^(-eps) (1 + x)^(2 eps - 2), {x}, {}, 2, eps]];
  difference = Simplify[Normal[Series[value - Gamma[1 - eps]/Gamma[2 - 2 eps], {eps, 0, 2}]]];
  (* the one-loop massless box of the hadronic processes: a three-variable integral whose
     evaluation needs polymake; a failed polymake call used to give 0 *)
  box = guarded[Quiet[EulerIntegral[Gamma[2 + eps], (1 + x1 + x2 + x3)^(2 eps) (x1 x3 + r x2)^(-2 - eps),
    {x1, x2, x3}, {r}, 0, eps], Series::ztest1]];
  calls = Lookup[Last[$SubTropicaTimings], "PolymakeCalls", 0];
  poles = Simplify[{Coefficient[Expand[box], eps, -2] - 4/r,
    Coefficient[Expand[box], eps, -1] + 2 (2 EulerGamma + Log[r])/r}];
  If[difference === 0 && poles === {0, 0} && calls > 0,
   emit["pass", <|"integrals" -> "Gamma(1-eps)/Gamma(2-2eps) through eps^2; box poles",
     "polymakeCalls" -> calls|>],
   emit["fail", <|"difference" -> ToString[difference, InputForm],
     "boxPoles" -> ToString[poles, InputForm], "polymakeCalls" -> calls|>]]],

 _,
 emit["fail", "unknown test " <> checkName]];
