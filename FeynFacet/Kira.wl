(* ::Package:: *)

(* Kira.wl: write target lists, run Kira, read its reduction rules and
   restore the physical scale.  The integral families themselves are the
   YAML files shipped with each process; this module never invents one. *)

BeginPackage["FeynFacet`"];

PrepareKira::usage = "PrepareKira[sourceDir, workDir] copies the Kira \
configuration tree sourceDir (config/, jobs.yaml, preferred) to a fresh workDir, \
deleting any earlier reduction there, and returns workDir.";

WriteKiraTargets::usage = "WriteKiraTargets[workDir, targets, nameRules] writes the \
file workDir/targets, one integral per line, translating Wolfram heads to Kira \
family names with nameRules, e.g. {KRealAB -> \"real_ab\"}.";

RunKira::usage = "RunKira[workDir] runs 'kira jobs.yaml' in workDir, saves \
kira.log there, and aborts on a non-zero exit code.  Kira needs Fermat, which \
is found automatically (ToolPath[\"Fermat\"]) or given in Paths.wl.";

KiraRules::usage = "KiraRules[workDir, family, nameRules] reads \
workDir/results/<family>/kira_targets.m and returns Wolfram replacement rules \
target -> combination of masters, translating family names back to heads.";

RestoreScale::usage = "RestoreScale[rules, scale, invariantRules, heads] undoes Kira's \
symbol_to_replace_by_one: every master integral (an expression with a head from \
the list heads and integer indices) is multiplied by scale^(sum of master \
indices - sum of target indices), and the dimensionless invariants are \
restored with invariantRules, e.g. {t -> t/s, s4 -> s4/s}.";

CheckReductionCoverage::usage = "CheckReductionCoverage[targets, rules, heads, masters] rejects a missing target rule or a residual integral outside the explicitly supported master basis, before any master relabelling.";

Begin["`Private`"];

CheckReductionCoverage[targets_List, rules_List, heads_List, masters_List] :=
 Module[{missing, residual, unknown},
  missing = Complement[targets, Union[First /@ rules, masters]];
  residual = Union[Cases[targets /. rules, (Alternatives @@ heads)[__Integer], Infinity]];
  unknown = Complement[residual, masters];
  Require[missing === {} && unknown === {},
   "Kira did not reduce every target to the supported master basis",
   <|"MissingTargetRules" -> missing, "UnsupportedResidualIntegrals" -> unknown|>];
  <|"Targets" -> Length[targets], "ResidualMasters" -> Length[residual]|>];

PrepareKira[sourceDir_String, workDir_String] := (
  Require[DirectoryQ[sourceDir], "Kira source directory not found", sourceDir];
  If[DirectoryQ[workDir], DeleteDirectory[workDir, DeleteContents -> True]];
  If[!DirectoryQ[DirectoryName[workDir]], CreateDirectory[DirectoryName[workDir], CreateIntermediateDirectories -> True]];
  Require[CopyDirectory[sourceDir, workDir] =!= $Failed, "Could not copy the Kira configuration", {sourceDir, workDir}];
  workDir);

WriteKiraTargets[workDir_String, targets_List, nameRules_List] := Module[{text},
  text = StringRiffle[
    StringReplace[ToString[# /. nameRules, InputForm], " " -> ""] & /@ targets, "\n"] <> "\n";
  text = StringReplace[text, "\"" -> ""];
  Export[FileNameJoin[{workDir, "targets"}], text, "Text"];
  targets];

RunKira[workDir_String] := Module[{run, log, fermat},
  fermat = ToolPath["Fermat"];
  (* FERMATPATH through /usr/bin/env rather than ProcessEnvironment, which can drop a variable already set (see RunSteps) *)
  run = RunProcess[{"/usr/bin/env", "FERMATPATH=" <> fermat, ToolPath["Kira"], "jobs.yaml"}, ProcessDirectory -> workDir];
  log = FileNameJoin[{workDir, "kira.log"}];
  Export[log, run["StandardOutput"] <> "\n" <> run["StandardError"], "Text"];
  Require[run["ExitCode"] === 0, "Kira failed; see " <> log, run["ExitCode"]];
  run["ExitCode"]];

KiraRules[workDir_String, family_String, nameRules_List] := Module[{file, text},
  file = FileNameJoin[{workDir, "results", family, "kira_targets.m"}];
  Require[FileExistsQ[file], "Kira did not write its result file", file];
  text = Import[file, "Text"];
  (* translate Kira family names (strings) into Wolfram heads *)
  text = StringReplace[text,
    Reverse /@ (nameRules /. Rule[head_Symbol, name_String] :> Rule[SymbolName[head], name]) /.
      Rule[name_String, headName_String] :> (name -> headName)];
  ToExpression[text, InputForm]];

RestoreScale[rules_List, scale_Symbol, invariantRules_List, heads_List] :=
  rules /. Rule[lhs_, rhs_] :> Rule[lhs,
    (rhs /. invariantRules) /.
      (master : (Alternatives @@ heads)[__Integer]) :>
        scale^(Total[List @@ master] - Total[List @@ lhs]) master];

End[];
EndPackage[];
