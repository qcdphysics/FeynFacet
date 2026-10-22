(* ::Package:: *)

(* Masters.wl: evaluate Euler-type master integrals with SubTropica
   (HyperIntica backend) and return their epsilon expansions in the
   notebook's own symbol eps.  The integral representations are written in
   the notebooks; this module only supplies the call and the bookkeeping. *)

BeginPackage["FeynFacet`"];

LoadSubTropica::usage = "LoadSubTropica[] loads SubTropica and HyperIntica, \
tells SubTropica where polymake is, and keeps its symbols out of the Global \
context.  Call once per kernel before EulerIntegral.";

EulerIntegral::usage = "EulerIntegral[prefactor, integrand, variables, parameters, order, eps] \
integrates integrand over every variable from 0 to infinity, as an expansion in \
eps around 0 through order, using SubTropica`STIntegrate with the HyperIntica \
backend.  Hyperlogarithms are converted to PolyLog and Zeta.  Example: \
EulerIntegral[1/Gamma[1-eps], x^(-eps) (1+x)^(2eps-2), {x}, {}, 2, eps].";

$SubTropicaTimings::usage = "$SubTropicaTimings records elapsed seconds, expansion order and the number of polymake calls for each STIntegrate call in the current kernel. It does not change the integral or its evaluation options.";
$SubTropicaTimings = {};

Begin["`Private`"];

LoadSubTropica[] := Module[{},
  Needs["SubTropica`"];
  SubTropica`$PolymakeCommand = ToolPath["Polymake"];
  SetOptions[HyperIntica`HyperInt, "EvaluatePeriodsQ" -> True];
  (* the notebooks use their own eps, x, z ...; do not let SubTropica's
     exported names shadow them *)
  $ContextPath = DeleteCases[$ContextPath, "SubTropica`" | "HyperIntica`"];
  "SubTropica " <> ToString[SubTropica`$SubTropicaVersion]];

(* SubTropica can turn a failed external call into an ordinary zero. Audit
   both sides of the boundary: count launches in the calling kernel, then
   require one complete record per invocation. This also detects a wrapper
   that cannot start or create its record. KernelsAvailable -> 1 below keeps
   every launch in this dynamically scoped audit, including shell batches. *)
shellQuote[s_String] := "'" <> StringReplace[s, "'" -> "'\\''"] <> "'";

polymakeWrapper[polymake_String, wrapperDirectory_String, recordDirectory_String] := Module[{path},
  path = FileNameJoin[{wrapperDirectory, "feynfacet-polymake-" <> CreateUUID[]}];
  Require[StringFreeQ[path, WhitespaceCharacter],
    "The polymake wrapper path contains a space, which SubTropica cannot call", path];
  Require[StringQ[Export[path, StringRiffle[{
    "#!/bin/sh",
    "# One record per invocation; status is separate from arbitrary arguments.",
    "record=$(mktemp " <> shellQuote[FileNameJoin[{recordDirectory, "call.XXXXXX"}]] <> ") || exit 70",
    "printf '%s\\n' \"$@\" > \"$record\" || exit 70",
    shellQuote[polymake] <> " \"$@\" > \"$record.stdout\" 2> \"$record.stderr\"",
    "status=$?",
    "cat \"$record.stdout\" || exit 70",
    "printf '%s\\n' \"$status\" > \"$record.status\" || exit 70",
    "exit \"$status\""}, "\n"] <> "\n", "Text"]],
    "Cannot write the polymake wrapper", path];
  Require[RunProcess[{"chmod", "755", path}, "ExitCode"] === 0,
    "Cannot make the polymake wrapper executable", path];
  path];

recordPolymakeLaunch[arguments_List, wrapper_String] := Module[{count},
  count = Total[StringCount[#, wrapper] & /@
    Cases[{First[arguments, {}]}, _String, Infinity]];
  If[count > 0, AppendTo[$PolymakeLaunches,
    <|"Calls" -> count, "Arguments" -> ToString[arguments, InputForm]|>]]];

SetAttributes[withPolymakeAudit, HoldFirst];
withPolymakeAudit[expression_, wrapper_String] :=
  Internal`InheritedBlock[{Run, RunProcess, StartProcess},
    Unprotect[Run, RunProcess, StartProcess];
    Run[args___] /; TrueQ[$AuditPolymakeLaunches] :=
      (recordPolymakeLaunch[{args}, wrapper]; Block[{$AuditPolymakeLaunches = False}, Run[args]]);
    RunProcess[args___] /; TrueQ[$AuditPolymakeLaunches] :=
      (recordPolymakeLaunch[{args}, wrapper]; Block[{$AuditPolymakeLaunches = False}, RunProcess[args]]);
    StartProcess[args___] /; TrueQ[$AuditPolymakeLaunches] :=
      (recordPolymakeLaunch[{args}, wrapper]; Block[{$AuditPolymakeLaunches = False}, StartProcess[args]]);
    (* ProcessLink has more specific built-in definitions. Put our conditional
       rule first, then delegate with the flag off to the untouched originals. *)
    DownValues[RunProcess] = Join[
      Select[DownValues[RunProcess], !FreeQ[#, $AuditPolymakeLaunches] &],
      Select[DownValues[RunProcess], FreeQ[#, $AuditPolymakeLaunches] &]];
    DownValues[StartProcess] = Join[
      Select[DownValues[StartProcess], !FreeQ[#, $AuditPolymakeLaunches] &],
      Select[DownValues[StartProcess], FreeQ[#, $AuditPolymakeLaunches] &]];
    Protect[Run, RunProcess, StartProcess];
    Block[{$AuditPolymakeLaunches = True}, expression]];

(* ReadString returns EndOfFile for an empty file; an unreadable/missing file
   must remain distinct from an empty, successful diagnostic stream. *)
polymakeRecordText[path_String] := If[!FileExistsQ[path], $Failed,
  Replace[Quiet[Check[ReadString[path], $Failed]], EndOfFile -> ""]];

polymakeCalls[recordDirectory_String] := Module[{records, rows},
  records = Select[FileNames["call.*", recordDirectory],
    Function[path, !Or @@ (StringEndsQ[path, #] & /@ {".status", ".stderr", ".stdout"})]];
  rows = Map[Function[record, Module[{arguments, errors, output, status, failed},
      arguments = polymakeRecordText[record];
      errors = polymakeRecordText[record <> ".stderr"];
      output = polymakeRecordText[record <> ".stdout"];
      status = polymakeRecordText[record <> ".status"];
      failed = !And @@ (StringQ /@ {arguments, errors, output, status});
      If[!failed, failed = StringTrim[status] =!= "0" ||
        StringContainsQ[errors <> output, "ERROR"]];
      <|"Record" -> record, "Call" -> arguments,
        "Exit" -> If[StringQ[status], StringTrim[status], "did not finish"],
        "Errors" -> If[StringQ[errors], StringTake[errors, UpTo[1200]], "missing or unreadable stderr"],
        "Output" -> If[StringQ[output], StringTake[output, UpTo[1200]], "missing or unreadable stdout"],
        "Failed" -> failed|>]], records];
  {Length[records], Select[rows, TrueQ[#["Failed"]] &]}];

EulerIntegral[prefactor_, integrand_, variables_List, parameters_List, order_Integer, eps_Symbol] :=
  Module[{raw, value, elapsed, workDirectory = None, temporaryDirectory, recordDirectory,
      wrapper = None, calls, failures, launches, expectedCalls, keep = True, polymake = SubTropica`$PolymakeCommand,
      systemTemporary = $TemporaryDirectory,
      originalDirectory = Directory[], originalDepth = Length[DirectoryStack[]]},
    Require[StringQ[polymake] && FileExistsQ[polymake],
      "polymake is not configured; call LoadSubTropica[] first", polymake];
    (* SubTropica writes relative integration files and clears its tracking
       directory. Give each call its own workspace, including the tracking
       paths that SubTropica fixes when the package is first loaded.  The
       workspace is deleted after a successful integration and kept otherwise. *)
    Internal`WithLocalSettings[Null,
      workDirectory = CreateDirectory[FileNameJoin[
        {ResultsDirectory[], "SubTropica-" <> CreateUUID[]}]];
      Require[StringQ[workDirectory], "Cannot create a SubTropica workspace"];
      temporaryDirectory = CreateDirectory[FileNameJoin[{workDirectory, "tmp"}]];
      Require[StringQ[temporaryDirectory], "Cannot create a SubTropica temporary directory"];
      recordDirectory = CreateDirectory[FileNameJoin[{workDirectory, "polymake-calls"}]];
      Require[StringQ[recordDirectory], "Cannot create the polymake record directory"];
      wrapper = polymakeWrapper[polymake, systemTemporary, recordDirectory];
      SetDirectory[workDirectory];
      Block[{$TemporaryDirectory = temporaryDirectory, SubTropica`$PolymakeCommand = wrapper,
          (* SubTropica also caches $Failed tropical data. A previous failure
             must not bypass this invocation's backend audit on retry. *)
          SubTropica`$STTropicalDataCache = <||>,
          $PolymakeLaunches = {},
          SubTropica`$STJobTrackingDir = FileNameJoin[{temporaryDirectory, "STJobTracking"}],
          SubTropica`$STCompletedJobsLog = FileNameJoin[
            {temporaryDirectory, "STJobTracking", "completed_jobs.m"}]},
        {elapsed, raw} = AbsoluteTiming[withPolymakeAudit[SubTropica`STIntegrate[
          {prefactor, integrand, variables, parameters} /. eps -> SubTropica`eps,
          "Order" -> order,
          "Integrator" -> "HyperIntica", "LROrderBackend" -> "HyperIntica",
          "KernelsAvailable" -> 1, "Rationalize" -> False,
          "SimplifyOutput" -> Identity, "Verbose" -> False, "ShowTimings" -> False,
          "ReuseExistingResults" -> False], wrapper]];
        launches = $PolymakeLaunches];
      expectedCalls = Total[Lookup[launches, "Calls", 0]];
      Export[FileNameJoin[{workDirectory, "polymake-launches.json"}], launches, "RawJSON"];
      {calls, failures} = polymakeCalls[recordDirectory];
      AppendTo[$SubTropicaTimings, <|"Seconds" -> elapsed, "ExpansionOrder" -> order,
        "PolymakeCalls" -> calls, "PolymakeLaunches" -> expectedCalls|>];
      Require[calls === expectedCalls,
        "polymake was launched without a complete audit record; the workspace is kept for diagnosis",
        <|"Workspace" -> workDirectory, "ExpectedCalls" -> expectedCalls, "RecordedCalls" -> calls,
          "Launches" -> launches|>];
      Require[failures === {}, "polymake failed during a SubTropica integral, so its result " <>
        "cannot be trusted; the workspace is kept for diagnosis",
        <|"Workspace" -> workDirectory, "FailedCalls" -> failures|>];
      Require[FreeQ[raw, $Failed | $Aborted | _Failure | Indeterminate | _DirectedInfinity] &&
        Head[raw] =!= SubTropica`STIntegrate,
        "SubTropica did not return a result; the workspace is kept for diagnosis",
        <|"Workspace" -> workDirectory, "Result" -> raw|>];
      value = raw /. HyperIntica`Hlog[arg_, word_List] :> HyperIntica`HlogAsMpl[arg, word];
      value = FixedPoint[(# /. {HyperIntica`mzv[n_Integer] :> Zeta[n],
        HyperIntica`Mpl[{n_Integer}, {arg_}] :> PolyLog[n, arg]}) &, value];
      value = Collect[Expand[Normal[Series[value, {SubTropica`eps, 0, order}]]],
        SubTropica`eps, Together] /. SubTropica`eps -> eps;
      Require[FreeQ[value, $Failed | $Aborted | _Failure | _Series | _SeriesData |
        Indeterminate | _DirectedInfinity],
        "SubTropica expansion did not finish; the workspace is kept for diagnosis", workDirectory];
      keep = False,
      AbortProtect[
        (* Also unwind any directory changes left by an interrupted backend. *)
        While[Length[DirectoryStack[]] > originalDepth, ResetDirectory[]];
        If[Directory[] =!= originalDirectory, SetDirectory[originalDirectory]];
        If[StringQ[wrapper] && FileExistsQ[wrapper], DeleteFile[wrapper]];
        If[!keep && StringQ[workDirectory] && DirectoryQ[workDirectory],
          DeleteDirectory[workDirectory, DeleteContents -> True]]]];
    value];

End[];
EndPackage[];
