(* Exact comparison of analytical coefficient tables in two directories. Metadata describing
   which routes were run is excluded; every mathematical entry is compared.
   Usage: wolframscript -file Checks/compare_coefficient_tables.wls OLD NEW REPORT *)
$HistoryLength = 0;
Get[FileNameJoin[{DirectoryName[$InputFileName], "ScriptArguments.wl"}]];
comparisonArguments = ScriptArguments[];
If[Length[comparisonArguments] =!= 3, Print["Give OLD NEW REPORT directories."]; Exit[2]];
{comparisonOld, comparisonNew, comparisonReport} = ExpandFileName /@ comparisonArguments;
comparisonFiles = {"PartonicBorn", "ProtonProton", "ProtonProtonTT", "ProtonProtonTransfer",
 "LeptonQuark", "LeptonOffDiagonal", "Collins", "TMD"};
comparisonRows = {}; comparisonFailures = {};
compareEntry[before_Association, after_Association, path_] := Module[{a, b},
 a = KeyDrop[before, {"Routes", "StepsThreeFour"}];
 b = KeyDrop[after, {"Routes", "StepsThreeFour"}];
 If[Sort[Keys[a]] =!= Sort[Keys[b]], AppendTo[comparisonFailures, {path, "Different keys"}]; Return[]];
 Do[compareEntry[a[key], b[key], Append[path, key]], {key, Keys[a]}]];
compareEntry[before_List, after_List, path_] := If[Length[before] =!= Length[after],
 AppendTo[comparisonFailures, {path, "Different lengths"}],
 Do[compareEntry[before[[entry]], after[[entry]], Append[path, entry]], {entry, Length[before]}]];
compareEntry[before_String, after_String, path_] := If[before =!= after,
 AppendTo[comparisonFailures, {path, "Changed descriptive text; inspect separately"}]];
compareEntry[before_, after_, path_] := Module[{difference},
 difference = TimeConstrained[Together[before - after], 600, $Aborted];
 AppendTo[comparisonRows, <|"Path" -> ToString[path, InputForm], "ExactAgreement" -> (difference === 0)|>];
 If[difference =!= 0, AppendTo[comparisonFailures, {path, Short[difference]}]]];
Do[
 Print["Exact coefficient comparison: ", comparisonFile];
 compareEntry[Get[FileNameJoin[{comparisonOld, comparisonFile <> ".wl"}]],
  Get[FileNameJoin[{comparisonNew, comparisonFile <> ".wl"}]], {comparisonFile}],
 {comparisonFile, comparisonFiles}];
Export[comparisonReport, <|"passed" -> (comparisonFailures === {}),
 "MathematicalEntries" -> Length[comparisonRows], "Checks" -> comparisonRows,
 "Failures" -> (ToString[#, InputForm] & /@ comparisonFailures)|>, "RawJSON"];
Exit[If[comparisonFailures === {}, 0, 1]];
