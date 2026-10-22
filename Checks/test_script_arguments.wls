(* Command-line formats used by the two supported Wolfram launchers. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "ScriptArguments.wl"}]];
passed = 0;
check[label_, condition_] := If[TrueQ[condition], passed++; Print["PASS ", label],
  Print["FAIL ", label]; Exit[1]];

check["wolframscript arguments", ScriptArguments[{"check.wls", "kernel", "/tmp/check"}, {}] ===
  {"kernel", "/tmp/check"}];
check["wolframscript with no user arguments", ScriptArguments[{"check.wls"}, {}] === {}];
check["prefer the script argument list", ScriptArguments[{"check.wls", "manifest-only"},
  {"WolframKernel", "-script", "other.wls", "other"}] === {"manifest-only"}];
check["direct kernel arguments", ScriptArguments[{},
  {"WolframKernel", "-script", "check.wls", "kernel", "/tmp/check"}] === {"kernel", "/tmp/check"}];
check["direct kernel with no user arguments", ScriptArguments[{},
  {"WolframKernel", "-script", "check.wls"}] === {}];
check["kernel options before the script", ScriptArguments[{},
  {"WolframKernel", "-noprompt", "-script", "check.wls", "FORM"}] === {"FORM"}];
check["file-mode fallback", ScriptArguments[{},
  {"wolframscript", "-file", "check.wls", "FORM"}] === {"FORM"}];
check["paths containing spaces stay single arguments", ScriptArguments[{},
  {"/Wolfram App/WolframKernel", "-script", "/a b/check.wls", "/a b/step.wl"}] === {"/a b/step.wl"}];
check["user options and launcher-like arguments are preserved", ScriptArguments[{},
  {"WolframKernel", "-script", "check.wls", "--flag", "-script", "value"}] ===
  {"--flag", "-script", "value"}];
check["a script path among user arguments is preserved", ScriptArguments[{},
  {"WolframKernel", "-script", "check.wls", "check.wls"}] === {"check.wls"}];
check["an interactive kernel has no script arguments", ScriptArguments[{}, {"WolframKernel"}] === {}];
check["empty command lines", ScriptArguments[{}, {}] === {}];
Print[passed, " script-argument tests passed."];
Exit[0];
