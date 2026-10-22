(* ::Package:: *)

(* Tools.wl: locate installed programs, name the output directories, and
   stop a calculation when a required check fails.  No physics. *)

BeginPackage["FeynFacet`"];

ToolPath::usage = "ToolPath[\"Kira\"], ToolPath[\"FORM\"], ToolPath[\"TFORM\"], \
ToolPath[\"Polymake\"], ToolPath[\"Tracer\"] or ToolPath[\"FORMColor\"] returns the \
full path of that installed program or file.  Paths.wl in the repository root \
can override the automatic search through $ToolPaths.";

$ToolPaths::usage = "$ToolPaths is an Association of manual tool locations, \
normally set in Paths.wl, e.g. $ToolPaths[\"Kira\"] = \"/opt/kira/bin/kira\".";

Require::usage = "Require[condition, message, details] prints a Failure and \
aborts the evaluation unless condition is True.  Every physics check in the \
notebooks is guarded by Require, so a run that reaches the end has passed them all.";

ProcessDirectory::usage = "ProcessDirectory[] is the directory of the notebook \
being evaluated (the process folder).  In a headless run it is supplied by the \
runner through the environment variable FEYNFACET_PROCESS_DIRECTORY.";

ResultsDirectory::usage = "ResultsDirectory[] and ResultsDirectory[route] return \
<process>/Results and <process>/Results/<route>, creating them if needed.";

SaveResult::usage = "SaveResult[\"Name\", expr, route] writes expr to \
Results/<route>/Name.wl (route omitted: Results/Name.wl) and returns expr.";

LoadResult::usage = "LoadResult[\"Name\", route] reads Results/<route>/Name.wl.";

CheckInstallation::usage = "CheckInstallation[] runs every installed tool once \
and returns a table with the location and version of each, and the install \
command for anything that is missing.";

RunSteps::usage = "RunSteps[\"LeptonQuark\"] or RunSteps[\"LeptonQuark\", \"FORM\"] evaluates \
every numbered step file of that process in order, each in a fresh kernel, \
exactly as if the notebooks were evaluated by hand (the default route is \
FeynCalc).  It prints one line per step and stops at the first failed check; \
the kernel output is in Processes/<process>/Results/<step>-<route>[-<override values>].log.  \
RunSteps[process, route, {3, 4}] runs only the given step numbers; a fourth \
argument such as \"class=2\" or \"spin=LL\" overrides the top-level choices of the step files.";

StepReceiptValidQ::usage = "StepReceiptValidQ[stepFile, choices, started] checks a completed \
step receipt against its exact choices and the current source, input and output \
hashes.  Omit started when explicitly reusing a previously completed step.";

InputFingerprint::usage = "InputFingerprint[expr, files] is a SHA-256 hex string \
of the expression expr (e.g. the amplitudes of step 1) and of the contents of the \
given files (e.g. the step file itself): the inputs of the intermediate files of \
a step.  It changes when any of them changes.";

ReusableQ::usage = "ReusableQ[file, fingerprint] is True if file exists and its \
companion file <file>.input, written by MarkReusable after file was complete, \
holds the same fingerprint.  A file without companion (written by an earlier \
version, or interrupted while being written) or with another fingerprint is \
not reusable.";

MarkReusable::usage = "MarkReusable[file, fingerprint] writes <file>.input with \
the fingerprint of the inputs file was computed from; call it after file is complete.";

ShellQuote::usage = "ShellQuote[\"string\"] quotes a string for bash (single \
quotes, with every ' written as '\\''), so that a path with spaces or quotes \
stays one argument.";

Begin["`Private`"];

(* ---- optional manual configuration ------------------------------------- *)
If[!AssociationQ[$ToolPaths], $ToolPaths = <||>];
With[{paths = FileNameJoin[{$FeynFacetDirectory, "Paths.wl"}]},
  If[FileExistsQ[paths], Get[paths]]];

(* ---- failure guard ------------------------------------------------------ *)
Require[condition_, message_String, details_: None] :=
  If[TrueQ[condition], True,
    Print["CHECK FAILED: ", message, If[details === None, "", "\n  details: " <> ToString[details, InputForm]]];
    Abort[]];

(* ---- directories -------------------------------------------------------- *)
ProcessDirectory[] := Module[{env = Environment["FEYNFACET_PROCESS_DIRECTORY"], nb},
  If[StringQ[env] && DirectoryQ[env], Return[ExpandFileName[env]]];
  nb = Quiet[NotebookDirectory[]];
  If[StringQ[nb], Return[nb]];
  Require[False, "Cannot determine the process directory. Open the notebook from \
its process folder, or set FEYNFACET_PROCESS_DIRECTORY."]];

ResultsDirectory[] := ensureDirectory[FileNameJoin[{ProcessDirectory[], "Results"}]];
ResultsDirectory[route_String] :=
  ensureDirectory[FileNameJoin[{ProcessDirectory[], "Results", route}]];
ensureDirectory[dir_] := (If[!DirectoryQ[dir],
  CreateDirectory[dir, CreateIntermediateDirectories -> True]]; dir);

SaveResult[name_String, expr_] :=
  (Put[expr, FileNameJoin[{ResultsDirectory[], name <> ".wl"}]]; expr);
SaveResult[name_String, expr_, route_String] :=
  (Put[expr, FileNameJoin[{ResultsDirectory[route], name <> ".wl"}]]; expr);
LoadResult[name_String] := loadFile[FileNameJoin[{ResultsDirectory[], name <> ".wl"}]];
LoadResult[name_String, route_String] :=
  loadFile[FileNameJoin[{ResultsDirectory[route], name <> ".wl"}]];
loadFile[file_] := (Require[FileExistsQ[file],
  "Missing input file; run the preceding step first", file]; Get[file]);

(* ---- program search ----------------------------------------------------- *)
searchDirectories[] := DeleteDuplicates[Select[Join[
  StringSplit[Replace[Environment["PATH"], Except[_String] -> ""],
    If[$OperatingSystem === "Windows", ";", ":"]],
  {FileNameJoin[{$HomeDirectory, ".local", "bin"}], "/usr/local/bin", "/usr/bin",
   "/opt/homebrew/bin", "/opt/local/bin"}], DirectoryQ]];

findProgram[names_List] := SelectFirst[
  Flatten[Table[FileNameJoin[{dir, name}], {dir, searchDirectories[]}, {name, names}]],
  FileExistsQ[#] && !DirectoryQ[#] &, $Failed];

(* Fermat: the environment variable FERMATPATH, or a file fer64 (fer on
   32-bit systems) in the usual places where its archive is extracted *)
findFermat[] := Module[{env = Environment["FERMATPATH"], roots, found},
  If[StringQ[env] && FileExistsQ[env], Return[env]];
  roots = Select[{FileNameJoin[{$HomeDirectory, "Applications"}], FileNameJoin[{$HomeDirectory, "Fermat"}],
    FileNameJoin[{$HomeDirectory, ".local", "opt"}], FileNameJoin[{$HomeDirectory, "Documents"}],
    FileNameJoin[{$HomeDirectory, "Downloads"}], "/usr/local", "/opt", "/Applications"}, DirectoryQ];
  found = Quiet[FileNames[{"fer64", "fer"}, roots, 4]];
  SelectFirst[found, !DirectoryQ[#] &, $Failed]];

colorCandidates[] := Join[
  Flatten[Table[FileNameJoin[{dir, rel}], {dir, searchDirectories[]},
    {rel, {"../share/form/color/color.h", "../share/form-packages/color/color.h"}}]],
  {"/opt/form-packages/color/color.h", "/usr/share/form/color/color.h",
   "/usr/local/share/form/color/color.h",
   FileNameJoin[{$HomeDirectory, ".local", "share", "form", "color", "color.h"}]}];

ToolPath[name_String] := Module[{manual = Lookup[$ToolPaths, name, None], found},
  If[StringQ[manual] && manual =!= "",
    found = ExpandFileName[manual];
    Require[FileExistsQ[found], "Configured tool not found: " <> name, manual];
    Return[found]];
  found = Switch[name,
    "FORM", findProgram[{"form"}],
    "TFORM", findProgram[{"tform"}],
    "Kira", findProgram[{"kira"}],
    "Polymake", findProgram[{"polymake"}],
    "Tracer", SelectFirst[Quiet[FindFile /@ {"Tracer`", "Tracer.m", "TRACER.m"}], StringQ, $Failed],
    "FORMColor", SelectFirst[colorCandidates[], FileExistsQ, $Failed],
    "Fermat", findFermat[],
    _, $Failed];
  Require[StringQ[found], "Cannot find " <> name <>
    ". Install it (see INSTALL.md) or set $ToolPaths[\"" <> name <> "\"] in Paths.wl."];
  ExpandFileName[found]];

(* ---- inputs of reusable intermediate files ------------------------------ *)
InputFingerprint[expr_, files_List] := (
  Require[AllTrue[files, FileExistsQ], "InputFingerprint: a file does not exist", Select[files, !FileExistsQ[#] &]];
  Hash[{Hash[expr, "SHA256"], FileHash[#, "SHA256"] & /@ files}, "SHA256", "HexString"]);
ReusableQ[file_String, fingerprint_String] := FileExistsQ[file] && FileExistsQ[file <> ".input"] &&
  StringTrim[Quiet[Import[file <> ".input", "Text"]] /. Except[_String] -> ""] === fingerprint;
MarkReusable[file_String, fingerprint_String] := (Export[file <> ".input", fingerprint, "Text"]; file);

ShellQuote[x_String] := "'" <> StringReplace[x, "'" -> "'\\''"] <> "'";

(* ---- one-shot installation check ---------------------------------------- *)
version[cmd_List] := Module[{r = Quiet[Check[RunProcess[cmd], $Failed]]},
  If[!AssociationQ[r], "not runnable",
    First[Select[StringSplit[r["StandardOutput"] <> r["StandardError"], "\n"],
      StringLength[#] > 0 &], "exit " <> ToString[r["ExitCode"]]]]];

installHint = <|
  "FeynCalc" -> "In Mathematica: Import[\"https://raw.githubusercontent.com/FeynCalc/feyncalc/master/install.m\"]; InstallFeynCalc[AutoInstallFeynArts -> True]",
  "FeynArts" -> "installed together with FeynCalc (AutoInstallFeynArts -> True)",
  "Tracer" -> "download Tracer.m from library.wolfram.com/infocenter/MathSource/2987 into " <> FileNameJoin[{$UserBaseDirectory, "Applications"}],
  "FORM" -> "brew install form (macOS) or the binary from github.com/form-dev/form/releases into ~/.local/bin",
  "TFORM" -> "comes with FORM",
  "FORM color.h" -> "curl -fL https://www.nikhef.nl/~form/maindir/packages/color/color.h -o ~/.local/share/form/color/color.h",
  "Kira" -> "kira.hepforge.org: Linux binary into ~/.local/bin, or build from source (docs/INSTALL-macOS.md)",
  "Fermat" -> "download from home.bway.net/lewis and extract under ~/Applications/Fermat (or set $ToolPaths[\"Fermat\"] in Paths.wl)",
  "SubTropica" -> "In Mathematica: PacletInstall[\"https://github.com/SubTropica/SubTropica/releases/download/v1.2.10/SubTropica-1.2.10.paclet\"]",
  "polymake" -> "brew tap apaffenholz/polymake && brew install polymake (macOS) or apt-get install polymake (Debian/Ubuntu)",
  "python3" -> "optional, only for Checks/run_process.py"|>;

(* probe[expr]: the value of expr, or "" when it fails.  ToolPath stops a calculation with Require (Print and
   Abort) when a program is missing; Check alone does not catch an Abort, so CheckAbort turns it into "" here
   and the table lists every tool, found or not, with what to do for the missing ones *)
SetAttributes[probe, HoldFirst];
probe[expr_] := Block[{Print}, Quiet[CheckAbort[Check[expr, ""], ""]]];

CheckInstallation[] := Module[{rows, status, hint},
  status[name_, location_, version_] := If[StringQ[location] && (FileExistsQ[location] || StringStartsQ[version, "Python"]),
    {name, location, version}, {name, Style["NOT FOUND", Red], hint[name]}];
  (* a path set in Paths.wl that does not exist is reported as such *)
  hint[name_] := With[{configured = Lookup[$ToolPaths, Replace[name, {"FORM color.h" -> "FORMColor", "polymake" -> "Polymake"}], None]},
    If[StringQ[configured] && configured =!= "", "$ToolPaths[\"" <> Replace[name, {"FORM color.h" -> "FORMColor", "polymake" -> "Polymake"}] <>
      "\"] in Paths.wl is " <> configured <> ", which does not exist; correct it or remove it. " <> installHint[name], installHint[name]]];
  rows = {
    {"Mathematica", $Version, ""},
    status["FeynCalc", probe[FindFile["FeynCalc`"]], ""],
    status["FeynArts", probe[FileNameJoin[{DirectoryName[FindFile["FeynCalc`"], 2], "FeynArts", "FeynArts.m"}]], ""],
    status["Tracer", probe[ToolPath["Tracer"]], ""],
    status["FORM", probe[ToolPath["FORM"]], probe[version[{ToolPath["FORM"], "-v"}]]],
    status["TFORM", probe[ToolPath["TFORM"]], ""],
    status["FORM color.h", probe[ToolPath["FORMColor"]], ""],
    status["Kira", probe[ToolPath["Kira"]], probe[version[{ToolPath["Kira"], "--version"}]]],
    status["Fermat", probe[ToolPath["Fermat"]], ""],
    status["SubTropica", probe[FindFile["SubTropica`"]], ""],
    status["polymake", probe[ToolPath["Polymake"]], probe[version[{ToolPath["Polymake"], "--version"}]]],
    {"python3 (optional)", Replace[probe[version[{"python3", "--version"}]], "" | "not runnable" -> "not found; only needed for Checks/run_process.py"], ""}};
  Grid[Prepend[rows, Style[#, Bold] & /@ {"Tool", "Location", "Version / what to do"}], Alignment -> Left,
    Dividers -> {False, {2 -> True}}, Spacings -> {2, 0.5}]];

(* ---- run the steps of a process from inside Mathematica ----------------- *)
stepFiles[process_String, route_String, numbers_] := Module[{dir, files},
  dir = FileNameJoin[{$FeynFacetDirectory, "Processes", process}];
  Require[DirectoryQ[dir], "No such process folder", dir];
  files = Select[FileNames["*.wl", dir], StringMatchQ[FileNameTake[#], DigitCharacter ~~ ___] &];
  files = Select[files, !StringContainsQ[FileNameTake[#], "-Algebra-"] || StringContainsQ[FileNameTake[#], "-Algebra-" <> route] &];
  If[numbers =!= All, files = Select[files, MemberQ[numbers, ToExpression[First[StringSplit[FileNameTake[#], "-"]]]] &]];
  Sort[files]];

stepReceiptFile[file_String, choices_List] := FileNameJoin[{DirectoryName[file], "Results",
  FileBaseName[file] <> If[choices === {}, "", "-" <> StringRiffle[choices[[All, 2]], "-"]] <> ".receipt.json"}];

StepReceiptValidQ[file_String, choices_List, started_: None, previousHash_: None] := Quiet[Check[Module[
  {receipt = stepReceiptFile[file, choices], data, cells, provenance, root, validFiles, expectedCells},
  If[!FileExistsQ[file] || !FileExistsQ[receipt], Return[False]];
  If[NumberQ[started] && AbsoluteTime[FileDate[receipt]] < Floor[started], Return[False]];
  If[StringQ[previousHash] && IntegerString[FileHash[receipt, "SHA256"], 16, 64] === previousHash, Return[False]];
  data = Import[receipt, "RawJSON"];
  If[!AssociationQ[data] || !TrueQ[Lookup[data, "completed", False]] ||
    ExpandFileName[Lookup[data, "file", ""]] =!= ExpandFileName[file] ||
    Lookup[data, "overrides", {}] =!= choices, Return[False]];
  cells = Lookup[data, "cells", {}];
  expectedCells = StringCount[Import[file, "Text"], "(* ::Input:: *)"];
  If[!ListQ[cells] || Length[cells] =!= expectedCells || cells === {} ||
    !AllTrue[cells, AssociationQ] || Lookup[cells, "cell", None] =!= Range[expectedCells] ||
    !AllTrue[Lookup[cells, "failed", True], # === False &], Return[False]];
  provenance = Lookup[data, "provenance", <||>];
  If[!AssociationQ[provenance] || !AllTrue[Lookup[provenance, {"sources", "inputs", "outputs"}, None], AssociationQ] ||
    provenance["sources"] === <||>, Return[False]];
  root = StringReplace[ExpandFileName[$FeynFacetDirectory], $PathnameSeparator ~~ EndOfString -> ""];
  validFiles[bindings_] := And @@ KeyValueMap[Function[{name, digest}, Module[{path},
    If[!StringQ[name] || !StringQ[digest], Return[False]];
    path = ExpandFileName[FileNameJoin[{root, name}]];
    StringStartsQ[path, root <> $PathnameSeparator] && FileExistsQ[path] && !DirectoryQ[path] &&
      IntegerString[FileHash[path, "SHA256"], 16, 64] === digest]], bindings];
  KeyExistsQ[provenance["sources"], StringDrop[ExpandFileName[file], StringLength[root] + 1]] &&
    And @@ (validFiles /@ Lookup[provenance, {"sources", "inputs", "outputs"}])
], False]];

RunSteps[process_String, route_String: "FeynCalc", numbers_: All, overrides_String: ""] := Module[
  {kernel = First[$CommandLine], script, runner, files, choices, suffix, env, log,
   command, run, output, passed, receipt, previousHash, start, summary = {}},
  runner = FileNameJoin[{$FeynFacetDirectory, "Checks", "run_step.wls"}];
  files = stepFiles[process, route, numbers];
  Require[files =!= {}, "No step files found for " <> process];
  choices = Prepend[If[overrides === "", {}, StringSplit[#, "=", 2] & /@ StringSplit[overrides, ";"]], {"route", route}];
  Require[AllTrue[choices, MatchQ[#, {_String, _String}] &] &&
    AllTrue[choices, StringMatchQ[StringRiffle[#, "="], RegularExpression["[A-Za-z][A-Za-z0-9]*=[A-Za-z0-9_,.+-]+"]] &] &&
    DuplicateFreeQ[choices[[All, 1]]], "Step choices must be distinct NAME=VALUE pairs", choices];
  suffix = "-" <> StringRiffle[choices[[All, 2]], "-"];
  script = findProgram[{"wolframscript"}];
  Do[
    start = AbsoluteTime[];
    receipt = stepReceiptFile[file, choices];
    previousHash = If[FileExistsQ[receipt], IntegerString[FileHash[receipt, "SHA256"], 16, 64], None];
    Print[FileBaseName[file], " ... "];
    (* Explicit env arguments carry the same choices to both supported launchers. *)
    env = {"FEYNFACET_SET=" <> StringRiffle[StringRiffle[#, "="] & /@ choices, ";"], "FEYNFACET_STEP=" <> file};
    command = Join[{"/usr/bin/env"}, env, If[StringQ[script],
      {script, "-file", runner, file}, {kernel, "-script", runner}]];
    run = Quiet[Check[RunProcess[command, ProcessDirectory -> DirectoryName[file]], $Failed]];
    output = If[AssociationQ[run], Lookup[run, "StandardOutput", ""] <> "\n" <> Lookup[run, "StandardError", ""],
      "The step kernel could not be started."];
    log = FileNameJoin[{DirectoryName[file], "Results", FileBaseName[file] <> suffix <> ".log"}];
    If[!DirectoryQ[DirectoryName[log]], CreateDirectory[DirectoryName[log]]];
    Export[log, output, "Text"];
    passed = AssociationQ[run] && Lookup[run, "ExitCode", 1] === 0 && StepReceiptValidQ[file, choices, start, previousHash];
    AppendTo[summary, {FileBaseName[file], If[passed, "passed", "FAILED"], Round[AbsoluteTime[] - start]}];
    Print["   ", If[passed, "passed", "FAILED"], ", ", Round[AbsoluteTime[] - start], " s"];
    If[!passed,
      Print[StringTake[output, -Min[2000, StringLength[output]]]];
      Print["Stopped: the step needs a successful exit and a fresh, hash-bound receipt. Full output: ", log];
      Break[]],
    {file, files}];
  Grid[Prepend[summary, {"step", "result", "seconds"}], Alignment -> Left, Dividers -> {False, {2 -> True}}]];

End[];
EndPackage[];
