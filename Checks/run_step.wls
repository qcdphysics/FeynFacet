(* Evaluate one notebook file (a .wl file with (* ::Input:: *) cells) in this
   kernel, cell by cell, the way a user would with Shift+Enter.  Any message,
   Failure, $Failed, Abort or False result stops the run.  Used by
   run_process.py; can also be called directly:
     wolframscript -file Checks/run_step.wls Processes/LeptonQuark/1-Amplitudes.wl
   Overrides such as spin or route are passed through the environment variable
   FEYNFACET_SET, e.g. FEYNFACET_SET="route=FORM;spin=Polarized".  The receipt
   Results/<step>[-<override values>].receipt.json records every cell and the
   provenance of the run (SHA-256 of sources, inputs and outputs; Git commit). *)
$HistoryLength = 0;
(* The runner's own variables live in the context FeynFacetRunner`, so that none
   of them (cells, rows, current, ...) is visible to the evaluated cells; the
   cells themselves are parsed in Global` (see the Block below). *)
Begin["FeynFacetRunner`"];
argumentHelper = FileNameJoin[{DirectoryName[$InputFileName], "ScriptArguments.wl"}];
Get[argumentHelper];
(* the step file: last command-line argument (wolframscript -file ... step.wl), or the
   environment variable FEYNFACET_STEP when launched from RunSteps[] *)
stepArguments = ScriptArguments[];
stepFile = If[stepArguments =!= {}, Last[stepArguments], Environment["FEYNFACET_STEP"]];
If[!StringQ[stepFile], Print["No step file given"]; Quit[2]];
stepFile = ExpandFileName[stepFile];
If[!FileExistsQ[stepFile], Print["No such file: ", stepFile]; Quit[2]];
SetEnvironment["FEYNFACET_PROCESS_DIRECTORY" -> DirectoryName[stepFile]];
(* the notebooks locate their files with NotebookDirectory[]; give it the same meaning here *)
Unprotect[NotebookDirectory]; NotebookDirectory[] = DirectoryName[stepFile]; Protect[NotebookDirectory];
overrides = Select[StringSplit[Replace[Environment["FEYNFACET_SET"], Except[_String] -> ""], ";"],
  StringContainsQ[#, "="] &];
overrides = (StringSplit[#, "=", 2] &) /@ overrides;

(* provenance: SHA-256 of the sources (the step file, the package, the other files of the
   process folder except notebooks and Markdown, this runner) at the start, of every repository file the cells read
   (Get, Import; hashed when read) and of every file they write (Put, Export, SaveResult;
   hashed at the end).  Files under kira/ working directories are left out (they are
   intermediate).  The receipt also records the Git commit and the sources that differ from it. *)
repoRoot = ExpandFileName[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], ".."}]];
processDir = ExpandFileName[DirectoryName[stepFile]];
sha[f_] := IntegerString[FileHash[f, "SHA256"], 16, 64];
relative[f_] := StringReplace[ExpandFileName[f], StartOfString ~~ repoRoot ~~ "/" -> ""];
loggedQ[f_] := StringQ[f] && StringStartsQ[f, repoRoot <> "/"] && !StringContainsQ[f, "/kira/"] &&
  !StringStartsQ[f, FileNameJoin[{repoRoot, "FeynFacet"}] <> "/"];
sourceFiles = Union[Join[{stepFile, ExpandFileName[$InputFileName], argumentHelper},
  FileNames["*.wl", FileNameJoin[{repoRoot, "FeynFacet"}], Infinity],
  Select[FileNames[All, processDir, Infinity],
    !DirectoryQ[#] && !StringStartsQ[#, FileNameJoin[{processDir, "Results"}]] &&
      !StringEndsQ[#, ".nb" | ".md"] &&
      (* FeynArts creates and removes these FORM scratch streams at runtime.
         They are generated data, not scientific source dependencies. *)
      !StringMatchQ[FileNameTake[#], RegularExpression["xform[0-9]+\\.str"]] &]]];
sources = Association[relative[#] -> sha[#] & /@ sourceFiles];
git[args__] := Quiet[Check[StringTrim[RunProcess[{"git", "-C", repoRoot, args}, "StandardOutput"]], ""]];
gitHead = git["rev-parse", "HEAD"];
gitDirty = Select[StringTrim[StringDrop[#, 2]] & /@ Select[StringSplit[git["status", "--porcelain"], "\n"], StringLength[#] > 3 &],
  KeyExistsQ[sources, #] &];
inputs = <||>; written = {};
logRead[f_String] := With[{path = Quiet[FindFile[f]]}, If[loggedQ[path] && !KeyExistsQ[inputs, relative[path]] && !MemberQ[written, path],
  inputs[relative[path]] = sha[path]]];
logWrite[f_String] := With[{path = ExpandFileName[f]}, If[loggedQ[path], written = Union[Append[written, path]]]];
Unprotect[Get, Import, Put, Export];
Get[f_String, rest___] /; !TrueQ[$inLog] := Block[{$inLog = True}, logRead[f]; Get[f, rest]];
Import[f_String, rest___] /; !TrueQ[$inLog] := Block[{$inLog = True}, logRead[f]; Import[f, rest]];
Put[args__, f_String] /; !TrueQ[$inLog] := Block[{$inLog = True}, With[{r = Put[args, f]}, logWrite[f]; r]];
Export[f_String, rest___] /; !TrueQ[$inLog] := Block[{$inLog = True}, With[{r = Export[f, rest]}, logWrite[f]; r]];
Protect[Get, Import, Put, Export];
provenance[final_] := <|"gitHead" -> gitHead, "sourcesDifferingFromGit" -> gitDirty, "sources" -> sources, "inputs" -> inputs,
  "outputs" -> If[final, Association[relative[#] -> sha[#] & /@ Select[written, FileExistsQ]], Association[relative[#] -> Null & /@ written]]|>;

text = Import[stepFile, "Text"];
(* split into cells at the package-notebook markers *)
pieces = StringSplit[text, RegularExpression["(?m)^\\(\\* ::(\\w+):: \\*\\)\\s*$"] -> "$1"];
cells = {}; current = None;
Do[If[StringMatchQ[piece, RegularExpression["^\\w+$"]], current = piece,
  If[current === "Input" && StringTrim[piece] =!= "", AppendTo[cells, piece]]], {piece, pieces}];
If[cells === {}, Print["No input cells found in ", stepFile]; Quit[2]];

allowedMessages = {"HoldForm[ParallelMap::nopar]", "HoldForm[ParallelTable::nopar]",
  "HoldForm[ParallelDo::nopar]", "HoldForm[General::stop]"};
rows = {};
receiptFile = FileNameJoin[{DirectoryName[stepFile], "Results",
  FileBaseName[stepFile] <> If[overrides === {}, "", "-" <> StringRiffle[overrides[[All, 2]], "-"]] <> ".receipt.json"}];
If[!DirectoryQ[DirectoryName[receiptFile]], CreateDirectory[DirectoryName[receiptFile], CreateIntermediateDirectories -> True]];

(* the loop variable has an unusual name on purpose: a plain i, j, n ... would be
   visible to the evaluated cells and to the model files that FeynArts reads *)
Do[
  code = StringTrim[cells[[runStepCellNumber]]];
  Do[code = StringReplace[code,
      RegularExpression["(?m)^(\\s*)" <> ov[[1]] <> "\\s*=\\s*\"[^\"]*\""] -> ("$1" <> ov[[1]] <> " = \"" <> ov[[2]] <> "\"")],
    {ov, overrides}];
  Print["CELL ", runStepCellNumber, "/", Length[cells], "  ", DateString["Time"], "  ", StringTake[StringTrim[code], UpTo[70]]];
  Block[{$MessageList = {}, $Context = "Global`"},
    held = ToExpression[code, InputForm, HoldComplete];
    (* several statements in one cell: evaluate them as one compound expression *)
    held = Replace[held, HoldComplete[args__] :> HoldComplete[CompoundExpression[args]]];
    parseMessages = ToString[#, InputForm] & /@ $MessageList];
  If[!MatchQ[held, _HoldComplete], Print["FAIL: cell did not parse"]; Quit[1]];
  start = AbsoluteTime[];
  Block[{$MessageList = {}, $Context = "Global`"},
    value = CheckAbort[Catch[ReleaseHold[held], _,
      Function[{v, tag}, Failure["UncaughtThrow", <|"Value" -> v, "Tag" -> tag|>]]], $Aborted];
    messages = ToString[#, InputForm] & /@ $MessageList];
  unexpected = Union[Complement[messages, allowedMessages], parseMessages];
  failed = !FreeQ[value, $Failed | $Aborted | _Failure] || value === False || unexpected =!= {};
  AppendTo[rows, <|"cell" -> runStepCellNumber, "seconds" -> N[AbsoluteTime[] - start],
    "messages" -> messages, "result" -> ToString[Head[value], InputForm], "failed" -> failed|>];
  Block[{$inLog = True}, Export[receiptFile, <|"file" -> stepFile, "overrides" -> overrides, "cells" -> rows, "completed" -> False,
    "provenance" -> provenance[False]|>, "RawJSON"]];
  shown = StringTake[ToString[value, InputForm], UpTo[600]];
  If[failed, Print["FAIL in cell ", runStepCellNumber, ": ", unexpected, "  result: ", shown]; Quit[1]];
  Print["  -> ", shown];
  Print["  done in ", ToString[NumberForm[Last[rows]["seconds"], {6, 1}]], " s"],
  {runStepCellNumber, Length[cells]}];
(* these lines are read after the cells have run, when the cells may have created Global` symbols of the same names:
   the runner's variables are written with their full context *)
Block[{FeynFacetRunner`$inLog = True}, Export[FeynFacetRunner`receiptFile, <|"file" -> FeynFacetRunner`stepFile, "overrides" -> FeynFacetRunner`overrides,
  "cells" -> FeynFacetRunner`rows, "completed" -> True, "provenance" -> FeynFacetRunner`provenance[True]|>, "RawJSON"]];
Print["COMPLETED ", FileNameTake[FeynFacetRunner`stepFile], " (", Length[FeynFacetRunner`cells], " cells, ",
  ToString[NumberForm[Total[FeynFacetRunner`rows[[All, "seconds"]]], {8, 1}]], " s)"];
End[];
Quit[0];
