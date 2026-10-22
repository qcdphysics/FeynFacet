(* ::Package:: *)

(* FORM.wl: write amplitude records as FORM expressions, run a FORM program,
   and read back what it writes.  The FORM programs next to each process
   contain the algebra and are meant to be read.

   The FORM vocabulary mirrors the record vocabulary of Amplitudes.wl:
     Chain(Ubar(p),G(mu),Sl(k-r),G(nu),U(k))   open fermion line
     G5                                        gamma_5 inside a chain
     Den(q)        1/q^2                       Pol(p,mu)   polarization vector
     T(a,i,j)      (T^a)_{ij}                  dF(i,j)     fundamental delta
     f(a,b,c)      structure constant          dA(a,b)     adjoint delta
     Eik(q,n)      1/(n.q), eikonal (Wilson-line) denominator
     d_(mu,nu), p(mu), e_(..)                  FORM's own metric, vector, epsilon
   Couplings ee, gs, eq are symbols; the imaginary unit is i_. *)

BeginPackage["FeynFacet`"];

WriteFORMAmplitudes::usage = "WriteFORMAmplitudes[<|\"Born1\" -> records, ...|>, workDir] \
writes workDir/amplitudes.inc, in which each labelled list of records becomes \
one local expression (L Born1 = ...;), and workDir/declarations.inc with the \
symbols, vectors, indices and functions these expressions use.  With the option \
\"ColorLibrary\" -> True the colour objects are written for FORM's colour library \
color.h instead: declarations.inc includes color.h, generators become its tensor \
T(i,j,a) (adjoint index last), both deltas become d_, fundamental indices get \
dimension NR and adjoint ones NA, and the preprocessor variable COLOUR lists every \
internal (summed) colour index with its primed copy, for the conjugation procedure.";

RunFORM::usage = "RunFORM[program, workDir] runs the FORM program (a *.frm file) \
with workDir as working directory, so that #include amplitudes.inc and \
declarations.inc find the exported files; the program's own directory and the \
FORM colour library are also on the include path.  Output goes to \
workDir/<program>.log; a non-zero exit aborts.  RunFORM[program, workDir, \
\"Workers\" -> n] uses TFORM with n threads; \"Define\" -> {\"SPIN=P\"} passes preprocessor \
definitions (FORM's -d option).";

ReadFORM::usage = "ReadFORM[workDir, \"File.m\"] reads a Mathematica-format \
expression written by FORM's #write, turning rat(x,y) into x/y.";

Begin["`Private`"];

(* ---- record -> FORM text --------------------------------------------------- *)
(* colorLibrary is True while WriteFORMAmplitudes runs with "ColorLibrary" -> True:
   generators are then written in the index order of color.h, T(i,j,a), and the
   deltas as FORM's d_ (whose dimension follows from the declared indices) *)
colorLibrary = False;
generator[a_, i_, j_] := If[colorLibrary, FORMTc[i, j, a], FORMT[a, i, j]];
formString[x_] := StringReplace[ToString[x, InputForm], {
  " " -> "", "I" -> "i_", "*^" -> "e", "\\[" -> "", "]" -> ")", "[" -> "("}];

toFORM[Global`Chain[args__]] := "Chain(" <> StringRiffle[chainItem /@ {args}, ","] <> ")";
chainItem[Global`Ubar[p_]] := "Ubar(" <> formString[p] <> ")";
chainItem[Global`Uspinor[p_]] := "U(" <> formString[p] <> ")";
chainItem[Global`Vbar[p_]] := "Vbar(" <> formString[p] <> ")";
chainItem[Global`Vspinor[p_]] := "V(" <> formString[p] <> ")";
chainItem[Global`GammaMatrix[a_]] := "G(" <> formString[a] <> ")";
chainItem[Global`Slash[p_]] := "Sl(" <> formString[p] <> ")";
chainItem[Global`GammaFive] := "G5";

(* momentum components and scalar products are linear: (p-k)(mu) = p(mu) - k(mu) *)
linear[Global`MomentumComponent[p_Plus, a_]] := linear[Global`MomentumComponent[#, a]] & /@ p;
linear[Global`MomentumComponent[c_?NumericQ p_, a_]] := c linear[Global`MomentumComponent[p, a]];
linear[Global`MomentumComponent[p_Symbol, a_]] := FORMvec[p][a];
linear[Global`MomentumProduct[p_Plus, q_]] := linear[Global`MomentumProduct[#, q]] & /@ p;
linear[Global`MomentumProduct[p_, q_Plus]] := linear[Global`MomentumProduct[p, #]] & /@ q;
linear[Global`MomentumProduct[c_?NumericQ p_, q_]] := c linear[Global`MomentumProduct[p, q]];
linear[Global`MomentumProduct[p_, c_?NumericQ q_]] := c linear[Global`MomentumProduct[p, q]];
linear[Global`MomentumProduct[p_Symbol, q_Symbol]] := FORMsp[p, q];
scalarToFORM[x_] := formString[x /. {
  Global`Metric[a_, b_] :> FORMd[a, b],
  m_Global`MomentumComponent :> linear[m],
  m_Global`MomentumProduct :> linear[m],
  Global`PolVector[p_, a_, _] :> FORMPol[p, a],
  Global`ColorT[{a_}, i_, j_] :> generator[a, i, j],
  Global`ColorT[list_List, i_, j_] :> FORMtprod[list, i, j],
  Global`ColorDelta[i_, j_] :> If[colorLibrary, FORMd[i, j], FORMdF[i, j]],
  Global`GluonDelta[a_, b_] :> If[colorLibrary, FORMd[a, b], FORMdA[a, b]],
  Global`ColorF[a_, b_, c_] :> FORMf[a, b, c],
  Global`Epsilon4[a__] :> FORMe[a],
  Global`Eikonal[q_, v_] :> FORMEik[q, v]}];

(* products of generators: T(a,i,x1)*T(b,x1,j) with fresh dummies *)
tprod[list_, i_, j_] := Module[{n = Length[list], dummies},
  dummies = Table[Symbol["Global`Cx" <> ToString[k]], {k, n - 1}];
  Times @@ Table[generator[list[[k]], If[k == 1, i, dummies[[k - 1]]], If[k == n, j, dummies[[k]]]], {k, n}]];

(* FORM gamma strings take one vector per slot: Slash[k - r] is split into Slash[k] - Slash[r] *)
expandSlashes[chain_Global`Chain] := Module[{terms = {{1, {}}}},
  Do[terms = Flatten[Table[
      If[MatchQ[item, Global`Slash[_Plus]],
        Table[{term[[1]] Coefficient[item[[1]], m], Append[term[[2]], Global`Slash[m]]}, {m, Variables[item[[1]]]}],
        {{term[[1]] If[MatchQ[item, Global`Slash[Times[_?NumericQ, _]]], item[[1, 1]], 1],
          Append[term[[2]], item /. Global`Slash[Times[_?NumericQ, m_]] :> Global`Slash[m]]}}],
      {term, terms}], 1], {item, List @@ chain}];
  {#[[1]], Global`Chain @@ #[[2]]} & /@ terms];

recordToFORM[r_Association] := Module[{chainTerms},
  chainTerms = Tuples[expandSlashes /@ r["Chains"]];
  StringRiffle[Table["(" <> scalarToFORM[r["Coefficient"] Times @@ ct[[All, 1]] /.
      FORMtprod[l_, i_, j_] :> tprod[l, i, j]] <> ")" <>
    StringJoin["*" <> toFORM[#] & /@ ct[[All, 2]]], {ct, chainTerms}], "\n +"]];

finish[text_] := StringReplace[StringReplace[text, RegularExpression["FORMsp\\((\\w+),(\\w+)\\)"] -> "$1.$2"],
  {"FeynFacet`Private`" -> "", "Global`" -> "",
  "FORMd(" -> "d_(", "FORMe(" -> "e_(", "FORMvec(" -> "", ")(" -> "(",
  "FORMPol(" -> "Pol(", "FORMT(" -> "T(", "FORMTc(" -> "T(", "FORMdF(" -> "dF(", "FORMdA(" -> "dA(", "FORMf(" -> "f(",
  "FORMEik(" -> "Eik("}];

Options[WriteFORMAmplitudes] = {"ColorLibrary" -> False};
WriteFORMAmplitudes[amplitudes_Association, workDir_String, OptionsPattern[]] := Module[
  {lines, all, symbols, vectors, indices, lorentz, text, library = TrueQ[OptionValue["ColorLibrary"]],
   colourSlots, fundamental, adjoint, internal, dummies},
  If[!DirectoryQ[workDir], CreateDirectory[workDir, CreateIntermediateDirectories -> True]];
  lines = Block[{colorLibrary = library}, KeyValueMap[Function[{label, records},
    "L " <> label <> " =\n  " <> finish[StringRiffle[recordToFORM /@ records, "\n +"]] <> ";"],
    amplitudes]];
  text = "* Amplitudes written by FeynFacet (WriteFORMAmplitudes); see Amplitudes.wl for the vocabulary.\n" <>
    StringRiffle[lines, "\n"] <> "\n";
  Export[FileNameJoin[{workDir, "amplitudes.inc"}], text, "Text"];
  all = Flatten[Values[amplitudes]];
  vectors = Union[Cases[all, (Global`Ubar | Global`Uspinor | Global`Vbar | Global`Vspinor | Global`Slash |
      Global`Den | Global`MomentumComponent)[p_, ___] :> Variables[p], Infinity] // Flatten,
    Cases[all, Global`PolVector[p_, __] :> p, Infinity],
    Cases[all, Global`Eikonal[q_, v_] :> Join[Variables[q], Variables[v]], Infinity] // Flatten];
  indices = Union[Cases[all, (Global`GammaMatrix | Global`PolVector | Global`MomentumComponent)[_, a_, ___] :> a, Infinity],
    Cases[all, Global`Metric[a_, b_] :> {a, b}, Infinity] // Flatten,
    Cases[all, Global`GammaMatrix[a_] :> a, Infinity],
    Cases[all, (Global`ColorT[_, i_, j_] | Global`ColorDelta[i_, j_] | Global`GluonDelta[i_, j_]) :> {i, j}, Infinity] // Flatten,
    Cases[all, Global`ColorT[l_List, __] :> l, Infinity] // Flatten,
    Cases[all, Global`ColorF[a__] :> {a}, Infinity] // Flatten];
  symbols = Complement[Union[Cases[all[[All, "Coefficient"]], _Symbol, Infinity]],
    vectors, indices, {I, D}];
  lorentz = Select[indices, StringStartsQ[SymbolName[#], "Lor"] &];
  (* colour indices, for the colour library: fundamental (quark colours) and adjoint
     (gluon colours); an index that occurs twice in one coefficient is summed inside
     the amplitude (internal) and gets a primed copy for the complex conjugate *)
  colourSlots[c_] := Join[Cases[c, Global`ColorT[l_List, i_, j_] :> Join[l, {i, j}], {0, Infinity}] // Flatten,
    Cases[c, (Global`ColorDelta | Global`GluonDelta | Global`ColorF)[x__] :> {x}, {0, Infinity}] // Flatten];
  fundamental = Union[Flatten[Cases[all, (Global`ColorT[_, i_, j_] | Global`ColorDelta[i_, j_]) :> {i, j}, Infinity]]];
  adjoint = Complement[indices, lorentz, fundamental];
  internal = Union @@ Table[Select[Tally[colourSlots[r["Coefficient"]]], #[[2]] > 1 &][[All, 1]], {r, all}];
  dummies = {"Cx1", "Cx2", "Cx3", "Cx4"};
  Export[FileNameJoin[{workDir, "declarations.inc"}], StringJoin[
    "* Declarations written by FeynFacet (WriteFORMAmplitudes).\n",
    If[library, "* FORM's colour library: tensors T(i,j,a) and f(a,b,c), symbols NA, NR, cA, cR, I2R;\n" <>
      "* it sets Dimension NA, so Dimension D is restored right after it\n#include color.h\n", ""],
    "Symbols D", If[symbols === {}, "", "," <> StringRiffle[ToString /@ symbols, ","]], ";\n",
    "Dimension D;\n",
    If[vectors === {}, "", "Vectors " <> StringRiffle[ToString /@ vectors, ","] <> ";\n"],
    If[library,
      StringJoin[
        If[lorentz === {}, "", "Indices " <> StringRiffle[ToString /@ lorentz, ","] <> ";\n"],
        "* colour indices: fundamental ones have dimension NR, adjoint ones NA\n",
        "Indices " <> StringRiffle[# <> "=NR" & /@ Join[ToString /@ fundamental, dummies, "c" <> # & /@ dummies,
          "c" <> ToString[#] & /@ Intersection[internal, fundamental]], ","] <> ";\n",
        If[adjoint === {}, "", "Indices " <> StringRiffle[# <> "=NA" & /@ Join[ToString /@ adjoint,
          "c" <> ToString[#] & /@ Intersection[internal, adjoint]], ","] <> ";\n"],
        "* COLOUR lists every internal colour index with its primed copy, for the conjugation procedure\n",
        "#define COLOUR \"" <> StringRiffle[Flatten[{#, "c" <> #} & /@ Join[ToString /@ internal, dummies]], ","] <> "\"\n"],
      StringJoin[
        If[indices === {}, "", "Indices " <> StringRiffle[ToString /@ indices, ","] <> ";\n"],
        "Indices Cx1,Cx2,Cx3,Cx4;\n"]],
    "Indices " <> StringRiffle["c" <> ToString[#] & /@ lorentz, ","] <> ";\n",
    "* LORENTZ lists every Lorentz index with its primed copy, for the conjugation procedure\n",
    "#define LORENTZ \"" <> StringRiffle[Flatten[{ToString[#], "c" <> ToString[#]} & /@ lorentz], ","] <> "\"\n",
    If[library, "CFunctions Chain,Ubar,U,Vbar,V,G,Sl,Den,Pol,Eik;\n", "CFunctions Chain,Ubar,U,Vbar,V,G,Sl,Den,Pol,T,dF,dA,f,Eik;\n"],
    "Symbol G5;\n",
    "* four-dimensional twins of the Lorentz indices, and the procedure that moves\n",
    "* every index inside an epsilon tensor to its twin: in the HVBM scheme the\n",
    "* epsilon tensor lives in four dimensions, so indices contracted through it\n",
    "* must be summed over four dimensions only (call after tracen, before contract)\n",
    "Indices " <> StringRiffle[ToString[#] <> "4=4" & /@ Join[lorentz, "c" <> ToString[#] & /@ lorentz], ","] <> ";\n",
    "Indices y1,y2,y3;\n",
    "#procedure ProjectEpsilon\n  repeat;\n",
    StringJoin[Table[With[{x = ToString[index]},
      "    id e_(" <> x <> ",y1?,y2?,y3?) = e_(" <> x <> "4,y1,y2,y3);\n" <>
      "    id e_(y1?," <> x <> ",y2?,y3?) = e_(y1," <> x <> "4,y2,y3);\n" <>
      "    id e_(y1?,y2?," <> x <> ",y3?) = e_(y1,y2," <> x <> "4,y3);\n" <>
      "    id e_(y1?,y2?,y3?," <> x <> ") = e_(y1,y2,y3," <> x <> "4);\n"],
      {index, Join[lorentz, "c" <> ToString[#] & /@ lorentz]}]],
    "  endrepeat;\n#endprocedure\n"], "Text"];
  <|"amplitudes" -> Keys[amplitudes], "symbols" -> symbols, "vectors" -> vectors, "indices" -> indices|>];

(* ---- running FORM and reading its output ---------------------------------- *)
Options[RunFORM] = {"Workers" -> 1, "Define" -> {}};
RunFORM[program_String, workDir_String, OptionsPattern[]] := Module[
  {exe, path, run, log, workers = OptionValue["Workers"], defines = OptionValue["Define"]},
  Require[FileExistsQ[program], "FORM program not found", program];
  exe = If[workers > 1, ToolPath["TFORM"], ToolPath["FORM"]];
  path = StringRiffle[{workDir, DirectoryName[program], DirectoryName[ToolPath["FORMColor"]]}, ":"];
  run = RunProcess[Join[{exe}, If[workers > 1, {"-w" <> ToString[workers]}, {}],
      Flatten[{"-d", #} & /@ defines], {"-p", path, program}], ProcessDirectory -> workDir];
  log = FileNameJoin[{workDir, FileBaseName[program] <>
    StringReplace[StringJoin["-" <> # & /@ defines], Except[WordCharacter | "-" | "="] -> "_"] <> ".log"}];
  Export[log, run["StandardOutput"] <> "\n" <> run["StandardError"], "Text"];
  Require[run["ExitCode"] === 0, "FORM failed; see " <> log, run["ExitCode"]];
  run["ExitCode"]];

ReadFORM[workDir_String, file_String] := Module[{path = FileNameJoin[{workDir, file}]},
  Require[FileExistsQ[path], "FORM did not write the expected file", path];
  ReleaseHold[Get[path]] /. Global`rat[x_, y_] :> x/y];

End[];
EndPackage[];
