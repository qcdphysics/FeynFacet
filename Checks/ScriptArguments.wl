(* Arguments supplied after the script filename. wolframscript fills
   $ScriptCommandLine = {script, args...}; a direct WolframKernel -script
   leaves it empty and puts {kernel, ..., -script, script, args...} in
   $CommandLine instead. Locate source files with $InputFileName. *)
BeginPackage["FeynFacet`ScriptArguments`"];
ScriptArguments::usage = "ScriptArguments[] returns the user arguments for \
wolframscript -file or WolframKernel -script. ScriptArguments[script, command] \
accepts explicit command-line lists for testing.";
Begin["`Private`"];

ScriptArguments[] := ScriptArguments[$ScriptCommandLine, $CommandLine];
ScriptArguments[script_List, command_List] := Module[{position},
  If[script =!= {}, Return[Rest[script]]];
  position = FirstPosition[command, "-script" | "-file", Missing["NotFound"]];
  If[MissingQ[position], {}, Drop[command, Min[First[position] + 1, Length[command]]]]];

End[];
EndPackage[];
