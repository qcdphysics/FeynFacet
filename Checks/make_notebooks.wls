(* Convert every numbered step file Processes/<P>/N-*.wl into a Mathematica
   notebook N-*.nb next to it.  The .wl files are the source; the notebooks
   are a convenience for readers who prefer the notebook interface.
     wolframscript -file Checks/make_notebooks.wls [Processes/LeptonQuark]
   Each marker (* ::Section:: *) etc. becomes a cell of that style; the text
   cells are the (* ... *) comments, the input cells the code. *)
root = DirectoryName[ExpandFileName[$InputFileName], 2];
Get[FileNameJoin[{root, "Checks", "ScriptArguments.wl"}]];
arguments = ScriptArguments[];
folders = If[arguments =!= {}, {FileNameJoin[{root, First[arguments]}]},
  Select[FileNames["*", FileNameJoin[{root, "Processes"}]], DirectoryQ]];
cellsOf[text_] := Module[{pieces, cells = {}, style = None},
  pieces = StringSplit[text, RegularExpression["(?m)^\\(\\* ::(\\w+):: \\*\\)\\s*$"] -> "$1"];
  Do[If[StringMatchQ[piece, RegularExpression["^\\w+$"]], style = piece,
     If[style =!= None && StringTrim[piece] =!= "",
       AppendTo[cells, If[style === "Input",
         Cell[StringTrim[piece], "Input"],
         Cell[StringTrim[StringReplace[StringTrim[piece], {"(*" -> "", "*)" -> ""}]], style]]]]],
    {piece, pieces}];
  cells];
(* the Start notebook in the repository root *)
startFile = FileNameJoin[{root, "Start.wl"}];
If[FileExistsQ[startFile],
  Put[Notebook[cellsOf[Import[startFile, "Text"]], WindowSize -> {1000, 800}], FileNameJoin[{root, "Start.nb"}]];
  Print["wrote ", FileNameJoin[{root, "Start.nb"}]]];
Do[
  Do[
    nb = Notebook[cellsOf[Import[file, "Text"]], WindowSize -> {1000, 800}];
    Put[nb, StringReplace[file, ".wl" ~~ EndOfString -> ".nb"]];
    Print["wrote ", StringReplace[file, ".wl" ~~ EndOfString -> ".nb"]],
    {file, Select[FileNames["*.wl", folder], StringMatchQ[FileNameTake[#], DigitCharacter ~~ ___] &]}],
  {folder, folders}];
