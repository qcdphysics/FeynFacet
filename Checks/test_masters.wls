(* Regression tests for external failures that SubTropica can turn into zero.
   Usage: wolframscript -file Checks/test_masters.wls /absolute/scratch/path
   The scratch directory must be new. Failed-integral workspaces are retained. *)
$HistoryLength = 0;
root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{root, "Checks", "ScriptArguments.wl"}]];
Get[FileNameJoin[{root, "FeynFacet", "FeynFacet.wl"}]];

CheckAbort[
arguments = ScriptArguments[];
Require[Length[arguments] === 1, "Give one new absolute scratch directory"];
scratch = ExpandFileName[First[arguments]];
Require[!FileExistsQ[scratch], "The test scratch directory must be new", scratch];
CreateDirectory[scratch, CreateIntermediateDirectories -> True];
SetEnvironment["FEYNFACET_PROCESS_DIRECTORY" -> scratch];
passed = 0;
check[name_, condition_] := If[TrueQ[condition], passed++; Print["PASS ", name],
  Print["FAIL ", name]; Exit[1]];
script[name_, body_] := Module[{path = FileNameJoin[{scratch, name}]},
  Export[path, "#!/bin/sh\n" <> body <> "\n", "Text"];
  Require[RunProcess[{"chmod", "755", path}, "ExitCode"] === 0, "chmod failed", path]; path];
LoadSubTropica[];
realPolymake = SubTropica`$PolymakeCommand;
bad = script["fail-polymake", "echo 'backend failed' >&2\nexit 42"];
late = script["late-polymake", "printf '%0800d\\n' 0 >&2\necho 'ERROR: backend failed' >&2\nexit 0"];
interrupted = script["interrupted-polymake", "kill -TERM \"$PPID\"\nexit 0"];
box[] := EulerIntegral[Gamma[2 + eps],
  (1 + x1 + x2 + x3)^(2 eps) (x1 x3 + r x2)^(-2 - eps), {x1, x2, x3}, {r}, 0, eps];
workspaces[] := FileNames["SubTropica-*", ResultsDirectory[]];
SetAttributes[abortsAndKeeps, HoldRest];
abortsAndKeeps[name_, expression_] := Module[{before = workspaces[], result, added},
  result = CheckAbort[expression, $Aborted];
  added = Complement[workspaces[], before];
  check[name, result === $Aborted && Length[added] === 1]];

(* A healthy box immediately before breaking the backend is the acceptance
   test: a successful prerequisite check must not mask later failures. *)
value = box[];
check["healthy box double pole", Together[Coefficient[value, eps, -2] - 4/r] === 0];
check["healthy box single pole", Together[Coefficient[value, eps, -1] + 2 (2 EulerGamma + Log[r])/r] === 0];
check["healthy workspace removed", workspaces[] === {}];
SubTropica`$PolymakeCommand = bad;
abortsAndKeeps["nonzero exit aborts after healthy calculation", box[]];
abortsAndKeeps["failed data is not reused on same-kernel retry", box[]];
SubTropica`$PolymakeCommand = late;
abortsAndKeeps["ERROR after character 600 aborts", box[]];
SubTropica`$PolymakeCommand = interrupted;
abortsAndKeeps["unfinished wrapper aborts", box[]];

SubTropica`$PolymakeCommand = realPolymake;
oldPath = Environment["PATH"];
stubDirectory = CreateDirectory[FileNameJoin[{scratch, "stub-bin"}]];
stub = script["stub-bin/mktemp", "echo 'simulated record creation failure' >&2\nexit 1"];
Internal`WithLocalSettings[SetEnvironment["PATH" -> stubDirectory <> ":" <> oldPath],
  abortsAndKeeps["missing call record aborts", box[]], SetEnvironment["PATH" -> oldPath]];

unexecutable = FileNameJoin[{scratch, "unexecutable-wrapper"}];
Export[unexecutable, "#!/bin/sh\nexit 0\n", "Text"];
abortsAndKeeps["wrapper that cannot start aborts",
  Block[{FeynFacet`Private`polymakeWrapper = Function[{a, b, c}, unexecutable]}, box[]]];

before = workspaces[];
value = box[];
check["healthy recovery after failures", Together[Coefficient[value, eps, -2] - 4/r] === 0 && workspaces[] === before];
value = EulerIntegral[0, x^(-eps) (1 + x)^(2 eps - 2), {x}, {}, 0, eps];
check["legitimate zero is accepted", value === 0];
check["zero-call integral is accepted", Last[$SubTropicaTimings]["PolymakeCalls"] === 0];
abortsAndKeeps["backend Failure value aborts",
  Block[{SubTropica`STIntegrate = Function[Null, Failure["test", <||>]]}, box[]]];
abortsAndKeeps["backend aborted value aborts",
  Block[{SubTropica`STIntegrate = Function[Null, $Aborted]}, box[]]];
abortsAndKeeps["backend infinite value aborts",
  Block[{SubTropica`STIntegrate = Function[Null, Infinity]}, box[]]];

(* Test launcher coverage and restoration directly, including background
   batches, without relying on SubTropica's current choice of launch API. *)
RunProcess[{"/bin/true"}];
process = StartProcess[{"/bin/true"}];
While[ProcessStatus[process] === "Running", Pause[.01]]; KillProcess[process];
definitions[] := {Attributes[Run], DownValues[Run], Attributes[RunProcess], DownValues[RunProcess],
  Attributes[StartProcess], DownValues[StartProcess]};
originalDefinitions = definitions[];
recordDirectory = CreateDirectory[FileNameJoin[{scratch, "launcher-records"}]];
wrapper = FeynFacet`Private`polymakeWrapper[bad, $TemporaryDirectory, recordDirectory];
Block[{FeynFacet`Private`$PolymakeLaunches = {}},
  FeynFacet`Private`withPolymakeAudit[
    Run[ShellQuote[wrapper] <> " 'exit 0'"], wrapper];
  FeynFacet`Private`withPolymakeAudit[RunProcess[{wrapper, "exit 0"}], wrapper];
  FeynFacet`Private`withPolymakeAudit[
    process = StartProcess[{"/bin/sh", "-c", ShellQuote[wrapper] <> " & " <> ShellQuote[wrapper] <> " & wait"}];
    While[ProcessStatus[process] === "Running", Pause[.01]]; KillProcess[process], wrapper];
  check["Run, RunProcess and two-call StartProcess batch counted",
    Total[Lookup[FeynFacet`Private`$PolymakeLaunches, "Calls"]] === 4]];
{count, failures} = FeynFacet`Private`polymakeCalls[recordDirectory];
check["parallel records and exit-like arguments preserve failures", count === 4 && Length[failures] === 4 &&
  And @@ (#["Exit"] === "42" & /@ failures)];
check["process definitions restored after success", definitions[] === originalDefinitions];
CheckAbort[FeynFacet`Private`withPolymakeAudit[Abort[], wrapper], Null];
check["process definitions restored after abort", definitions[] === originalDefinitions];
DeleteFile[wrapper];
Print[passed, " master-integral failure tests passed. Diagnostics: ", scratch];
Exit[0], Exit[1]];
