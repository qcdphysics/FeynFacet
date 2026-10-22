#!/usr/bin/env wolframscript
(* ::Package:: *)

(* make_figures.wls: regenerate every figure of the results website.

   Usage (from anywhere):
       wolframscript -file site/make_figures.wls              all figures
       wolframscript -file site/make_figures.wls diagrams     Feynman-diagram tiles only
       wolframscript -file site/make_figures.wls plots        coefficient plots only

   Output, all under site/figures/:
       diagrams/<process>-<name>.svg and .png (PNG at 2x), one diagram per file,
       diagrams/index.json and diagrams/index.wl (tile list with captions),
       plots/<name>.svg and .png (PNG at 2x).
   The captions of every figure are in site/figures/captions.md.

   Diagrams.  They are generated with FeynArts exactly as in
   Processes/*/1-Amplitudes.wl (same CreateTopologies / InsertFields calls,
   model options and diagram selections); nothing is read from Results/.
   Paint is called with DisplayFunction -> Identity, which returns a
   FeynArtsGraphics object without needing a notebook front end; each
   DiagramGraphics in it is wrapped into its own FeynArtsGraphics, turned into
   an ordinary Graphics with Show, and exported.  Export to SVG/PNG starts a
   headless front end by itself, so the script runs in a plain terminal.
   FeynArts' particle labels of the SMQCD model (e, u, c) are replaced by
   the generic names of the paper (l, q, q').

   Plots.  The coefficients are read with Get from Coefficients/*.wl; only
   plain Mathematica is needed for this part. *)

$HistoryLength = 0;
Off[General::shdw];   (* FeynCalc and the coefficient files both use names such as Nf and CF *)
siteDir = DirectoryName[ExpandFileName[$InputFileName]];
rootDir = ParentDirectory[siteDir];
figDir = FileNameJoin[{siteDir, "figures"}];
diagDir = FileNameJoin[{figDir, "diagrams"}];
plotDir = FileNameJoin[{figDir, "plots"}];
Scan[If[!DirectoryQ[#], CreateDirectory[#, CreateIntermediateDirectories -> True]] &, {diagDir, plotDir}];

parts = Replace[Rest[$ScriptCommandLine], {} -> {"diagrams", "plots"}];
check[cond_, msg_] := If[!TrueQ[cond], Print["FAILED: ", msg]; Exit[1]];

(* PNG fallback at twice the nominal size (ImageSize is in printer's points, 72 dpi) *)
exportBoth[base_String, g_] := (
  Export[base <> ".svg", g];
  Export[base <> ".png", g, ImageResolution -> 144];
  Print["  wrote ", FileNameTake[base, -2], ".svg/.png"]);


(* ::Section:: *)
(*Part 1: Feynman-diagram tiles*)


(* FeynCalc and FeynArts are loaded before the definition below is read, so that
   its symbols (InsertFields, F, V, Paint, ...) are parsed in FeynArts' contexts *)
If[MemberQ[parts, "diagrams"],
  Get[FileNameJoin[{rootDir, "FeynFacet", "FeynFacet.wl"}]];
  Print[ToExpression["LoadFeynCalc[]"]]];


makeDiagrams[] := Module[{tiles = {}, tile, paintAll, smqcd, generate, excluded, born, real, virtual, vertex,
    realQG, realGQ, field, hexcl, hgen, classes, ggggg, tmdShapes, tmdGen, tmdSets, tmdVirtual, oldShapeDir, n},


  (* all DiagramGraphics of a TopologyList, in FeynArts' order *)
  paintAll[diagrams_] := Cases[
    Paint[diagrams, ColumnsXRows -> {1, 1}, Numbering -> None,
      SheetHeader -> None, AutoEdit -> False, DisplayFunction -> Identity],
    FeynArts`DiagramGraphics[___][___], Infinity];

  (* one tile: relabel the propagators, white background, fixed size *)
  tile[id_String, process_String, group_String, dg_, relabel_List, caption_String] := Module[{g},
    g = FeynArts`FeynArtsGraphics[][{{dg /. FeynArts`PropagatorGraphics[a_, l_, b_] :>
          FeynArts`PropagatorGraphics[a, Replace[l, relabel], b]}}, ImageSize -> {240, 200}];
    g = Show[First[Show[g]], Background -> White, ImageSize -> {240, 200}, ImagePadding -> {{6, 6}, {10, 6}}];
    exportBoth[FileNameJoin[{diagDir, id}], g];
    AppendTo[tiles, <|"id" -> id, "process" -> process, "group" -> group,
      "svg" -> "diagrams/" <> id <> ".svg", "png" -> "diagrams/" <> id <> ".png", "caption" -> caption|>]];

  (* ---- lepton-quark, as in Processes/LeptonQuark/1-Amplitudes.wl ---- *)
  excluded = {S[_], U[_], V[2 | 3], F[1 | 4]};
  generate[topologies_, fields_, extra_: {}] := InsertFields[topologies, fields,
    InsertionLevel -> {Classes}, Model -> "SMQCD", ExcludeParticles -> Join[excluded, extra]];
  smqcd = {"e" -> "\[ScriptL]", "u" -> "q", "c" -> "q'"};

  born = generate[CreateTopologies[0, 2 -> 2], {F[2, {1}], F[3, {1}]} -> {F[2, {1}], F[3, {1}]}, {V[5]}];
  real = DiagramSelect[
    generate[CreateTopologies[0, 2 -> 3], {F[2, {1}], F[3, {1}]} -> {F[2, {1}], F[3, {1}], V[5]}],
    !FreeQ[TreeFields[##], V[1, ___]] &];
  virtual = DiagramSelect[
    generate[CreateTopologies[1, 2 -> 2, ExcludeTopologies -> {Tadpoles}], {F[2, {1}], F[3, {1}]} -> {F[2, {1}], F[3, {1}]}],
    !FreeQ[LoopFields[##], V[5, ___]] && !FreeQ[TreeFields[##], V[1, ___]] &];
  vertex = DiagramSelect[virtual, WFCorrectionFields[##] === {} &];
  born = paintAll[born]; real = paintAll[real]; vertex = paintAll[vertex];
  check[Length /@ {born, real, vertex} === {1, 2, 1}, "lepton-quark diagram counts"];
  tile["leptonquark-born", "LeptonQuark", "Born", born[[1]], smqcd,
    "Born diagram of \[ScriptL] q \[RightArrow] \[ScriptL] q: one-photon exchange."];
  tile["leptonquark-vertex", "LeptonQuark", "Virtual", vertex[[1]], smqcd,
    "One-loop QCD vertex correction to \[ScriptL] q \[RightArrow] \[ScriptL] q (the external-leg self energies are scaleless and vanish)."];
  tile["leptonquark-real-1", "LeptonQuark", "Real", real[[1]], smqcd,
    "Real emission \[ScriptL] q \[RightArrow] \[ScriptL] q g: gluon radiated from the incoming quark."];
  tile["leptonquark-real-2", "LeptonQuark", "Real", real[[2]], smqcd,
    "Real emission \[ScriptL] q \[RightArrow] \[ScriptL] q g: gluon radiated from the outgoing quark."];

  (* ---- Collins: same diagrams as lepton-quark, transversely polarized quark ---- *)
  tile["collins-born", "Collins", "Born", born[[1]], {"e" -> "\[ScriptL]", "u" -> "q\[UpArrow]"},
    "Born diagram of \[ScriptL] q\[UpArrow] \[RightArrow] \[ScriptL] q\[UpArrow]: transverse spin is transferred from the incoming to the observed quark. The NLO diagrams are those of lepton-quark scattering."];

  (* ---- lepton off-diagonal, as in Processes/LeptonOffDiagonal/1-Amplitudes.wl ---- *)
  realQG = DiagramSelect[generate[CreateTopologies[0, 2 -> 3], {F[2, {1}], F[3, {1}]} -> {V[5], F[2, {1}], F[3, {1}]}],
    !FreeQ[TreeFields[##], V[1, ___]] &];
  realGQ = DiagramSelect[generate[CreateTopologies[0, 2 -> 3], {F[2, {1}], V[5]} -> {F[3, {1}], F[2, {1}], -F[3, {1}]}],
    !FreeQ[TreeFields[##], V[1, ___]] &];
  realQG = paintAll[realQG]; realGQ = paintAll[realGQ];
  check[Length /@ {realQG, realGQ} === {2, 2}, "lepton off-diagonal diagram counts"];
  Do[tile["leptonoffdiagonal-qg-" <> ToString[i], "LeptonOffDiagonal", "Real q \[RightArrow] g", realQG[[i]], smqcd,
    "\[ScriptL] q \[RightArrow] g + X (observed gluon): real-emission diagram " <> ToString[i] <> " of 2."], {i, 2}];
  Do[tile["leptonoffdiagonal-gq-" <> ToString[i], "LeptonOffDiagonal", "Real g \[RightArrow] q", realGQ[[i]], smqcd,
    "\[ScriptL] g \[RightArrow] q + X (observed quark): photon-gluon fusion, diagram " <> ToString[i] <> " of 2."], {i, 2}];

  (* ---- proton-proton, as in Processes/Hadronic/1-Amplitudes.wl ---- *)
  field = <|"q" -> F[3, {1}], "qb" -> -F[3, {1}], "qp" -> F[3, {2}], "qpb" -> -F[3, {2}], "g" -> V[5]|>;
  hexcl = {S[_], V[1 | 2 | 3], F[1 | 2 | 4]};
  hgen[loops_, in_, out_] := InsertFields[CreateTopologies[loops, Length[in] -> Length[out],
      ExcludeTopologies -> If[loops == 0, {}, {Tadpoles, WFCorrections}]],
    (field /@ in) -> (field /@ out), InsertionLevel -> {Classes}, Model -> "SMQCD",
    ExcludeParticles -> Join[hexcl, If[loops == 0, {U[_]}, {}]]];
  Get[FileNameJoin[{rootDir, "Processes", "Hadronic", "classes.wl"}]];
  classes = <|1 -> "q q' \[RightArrow] q q'", 13 -> "q g \[RightArrow] q g", 15 -> "g g \[RightArrow] g g"|>;
  KeyValueMap[Function[{c, label}, Module[{spec = HadronicClasses[c], dgs},
      dgs = paintAll[hgen[0, spec["Born"][[1 ;; 2]], spec["Born"][[3 ;; 4]]]];
      Print["class ", c, ": ", Length[dgs], " Born diagrams"];
      Do[tile["pp-class" <> ToString[c] <> "-born-" <> ToString[i], "Hadronic", "Born (class " <> ToString[c] <> ")", dgs[[i]], smqcd,
        "Born diagram " <> ToString[i] <> " of " <> ToString[Length[dgs]] <> " of " <> label <>
          " (channel class " <> ToString[c] <> ", " <> spec["Name"] <> ")."], {i, Length[dgs]}]]], classes];
  (* g g -> g g g, the first real-emission process of class 15 *)
  ggggg = paintAll[hgen[0, {"g", "g"}, {"g", "g", "g"}]];
  n = Length[ggggg]; Print["g g -> g g g: ", n, " diagrams"];
  check[n == 25, "g g -> g g g diagram count"];
  Do[tile["pp-ggggg-real-" <> ToString[i], "Hadronic", "Real (class 15)", ggggg[[i]], smqcd,
      "Real emission g g \[RightArrow] g g g (channel class 15): diagram " <> ToString[i] <> " of " <> ToString[n] <> "."],
    {i, {1, 13, 15, 16}}];

  (* ---- TMD, as in Processes/TMD/1-Amplitudes.wl, with the project's diagram layouts ---- *)
  (* the shape files are copied so that FeynArts cannot write into the repository *)
  tmdShapes = FileNameJoin[{$TemporaryDirectory, "feynfacet-site-shapes-" <> ToString[$ProcessID]}];
  CopyDirectory[FileNameJoin[{rootDir, "Processes", "TMD", "models", "ShapeData"}], tmdShapes];
  oldShapeDir = $ShapeDataDir; $ShapeDataDir = tmdShapes;
  $ModelPath = Prepend[$ModelPath, FileNameJoin[{rootDir, "Processes", "TMD", "models"}]];
  tmdGen[loops_, process_] := InsertFields[CreateTopologies[loops, Length[process[[1]]] -> Length[process[[2]]],
      ExcludeTopologies -> If[loops == 0, {}, {Tadpoles, WFCorrections}]], process,
    InsertionLevel -> {Classes}, Model -> "WilsonLineQCD", GenericModel -> "WilsonLine"];
  tmdSets = {
    {"bornq", tmdGen[0, {F[1]} -> {F[2], S[1]}], "Quark beam function, lowest order: q \[RightArrow] O_q W_n"},
    {"beamqq", tmdGen[0, {F[1]} -> {F[2], S[1], V[1]}], "Quark beam function I_qq at NLO: q \[RightArrow] O_q W_n g"},
    {"beamqg", tmdGen[0, {V[1]} -> {F[2], S[1], -F[1]}], "Quark beam function I_qg at NLO: g \[RightArrow] O_q W_n q\:0304"},
    {"borng", tmdGen[0, {V[1]} -> {V[2], S[3]}], "Gluon beam function, lowest order: g \[RightArrow] O_g W_n^A"},
    {"beamgq", tmdGen[0, {F[1]} -> {V[2], S[3], F[1]}], "Gluon beam function I_gq at NLO: q \[RightArrow] O_g W_n^A q"},
    {"beamgg", tmdGen[0, {V[1]} -> {V[2], S[3], V[1]}], "Gluon beam function I_gg at NLO: g \[RightArrow] O_g W_n^A g"},
    {"softq", tmdGen[0, {S[5]} -> {S[1], -S[2], V[1]}], "Quark soft function at NLO: O_S \[RightArrow] W_n W_n\:0304 g"},
    {"softg", tmdGen[0, {S[5]} -> {S[3], -S[4], V[1]}], "Gluon soft function at NLO: O_S \[RightArrow] W_n^A W_n\:0304^A g"},
    {"virtualq", DiagramSelect[tmdGen[1, {F[1]} -> {F[2], S[1]}], FreeQ[LoopFields[##], F[2] | V[2] | S[5]] &],
      "Quark beam function, one-loop virtual graph (scaleless, vanishes in dimensional regularization)"}};
  Scan[Function[set, Module[{dgs = paintAll[set[[2]]]},
      Do[tile["tmd-" <> set[[1]] <> If[Length[dgs] > 1, "-" <> ToString[i], ""], "TMD", set[[1]], dgs[[i]], {},
        set[[3]] <> If[Length[dgs] > 1, ", diagram " <> ToString[i] <> " of " <> ToString[Length[dgs]], ""] <> "."],
        {i, Length[dgs]}]]], tmdSets];
  $ShapeDataDir = oldShapeDir;
  Quiet[DeleteDirectory[tmdShapes, DeleteContents -> True]];

  Export[FileNameJoin[{diagDir, "index.json"}], tiles, "RawJSON", "Compact" -> False];
  Put[tiles, FileNameJoin[{diagDir, "index.wl"}]];
  Print[Length[tiles], " diagram tiles"]];


(* ::Section:: *)
(*Part 2: coefficient plots*)


makePlots[] := Module[{readCoefficients, pp, lq, col, ppSub, colours, dashes, panel, legend, sel, classLabel, regular, threshold,
    frame, lqSub, lqRegular, ratio, fRange, born, thr, vs},

  (* read with only System` on the context path, so that every symbol (Nf, D, LR, v, w, ...) is Global`
     even when FeynCalc has been loaded for the diagrams *)
  readCoefficients[name_] := Block[{$ContextPath = {"System`"}, $Context = "Global`"},
    Get[FileNameJoin[{rootDir, "Coefficients", name <> ".wl"}]]];
  pp = readCoefficients["ProtonProton"];
  lq = readCoefficients["LeptonQuark"];
  col = readCoefficients["Collins"];
  check[Keys[pp["Classes"][1]] === {"Process", "Born", "UU", "LL"}, "ProtonProton.wl structure"];
  check[Sort[Keys[col["Results"]]] === {"Transverse", "Unpolarized"}, "Collins.wl structure"];

  ppSub = {Global`Nc -> 3, Global`Nf -> 5, Global`LR -> 0, Global`LF -> 0, Global`LFF -> 0};
  lqSub = {Global`LR -> 0, Global`LI -> 0, Global`LO -> 0, Global`LL -> 0};
  (* categorical slots 1-4 of the reference palette, with a second encoding by dash *)
  colours = RGBColor /@ {"#2a78d6", "#eb6834", "#1baf7a", "#eda100"};
  dashes = {{}, {6, 3}, {1.5, 3}, {8, 3, 1.5, 3}};
  frame = {Frame -> True, FrameTicks -> Automatic, FrameStyle -> GrayLevel[0.25], FrameTicksStyle -> GrayLevel[0.15],
    LabelStyle -> Directive[FontFamily -> "Helvetica", FontSize -> 11, GrayLevel[0.1]], GridLines -> None,
    Axes -> {False, True}, AxesStyle -> GrayLevel[0.75], AspectRatio -> 0.72, ImageSize -> 270,
    ImagePadding -> {{58, 12}, {40, 26}}, Background -> White};
  sel = {1, 13, 15, 16};
  classLabel[c_] := "class " <> ToString[c] <> ": " <> StringReplace[pp["Classes"][c]["Process"],
    {"qpb" -> "q\:0304'", "qp" -> "q'", "qb" -> "q\:0304", "->" -> "\[RightArrow]"}];
  fRange[fs_, x_, {a_, b_}] := Module[{vals = Flatten[Table[Evaluate[fs], {x, Subdivide[a, b, 200]}]], lo, hi},
    {lo, hi} = MinMax[Select[vals, Im[#] == 0 &]]; {lo, hi} + {-1, 1} 0.06 (hi - lo)];
  legend[labels_, n_] := LineLegend[Table[Directive[colours[[i]], AbsoluteThickness[2], AbsoluteDashing[dashes[[i]]]], {i, n}],
    labels, LegendLayout -> "Row", LabelStyle -> Directive[FontFamily -> "Helvetica", FontSize -> 12, GrayLevel[0.1]],
    LegendMarkerSize -> {28, 12}];

  (* (a) proton-proton: regular coefficient vs w at v = 1/2; rows = classes, columns = UU, LL *)
  regular = Table[Module[{fs = Table[pp["Classes"][c][s][[1]] /. ppSub /. Global`v -> 1/2, {s, {"UU", "LL"}}], range},
      range = fRange[fs /. Global`w -> x, x, {1/10, 99/100}];
      Table[Plot[Evaluate[fs[[k]] /. Global`w -> x], {x, 1/10, 99/100}, PlotRange -> {{0.1, 1}, range},
        PlotStyle -> Directive[colours[[k]], AbsoluteThickness[2], AbsoluteDashing[dashes[[k]]]],
        PlotLabel -> classLabel[c] <> ", " <> {"UU", "LL"}[[k]],
        FrameLabel -> {Style["w", Italic], Row[{Subscript[Style["c", Italic], "regular"], "(", Style["v", Italic], " = 0.5, ", Style["w", Italic], ")"}]}, Evaluate[Sequence @@ frame]], {k, 2}]],
    {c, sel}];
  exportBoth[FileNameJoin[{plotDir, "pp-regular-v05"}], Legended[GraphicsGrid[regular, Spacings -> {10, 4}, ImageSize -> 2 * 290],
    Placed[legend[{"UU (spin average)", "LL (helicity difference)"}, 2], Top]]];

  (* (b) proton-proton: D0, D1, delta coefficients vs v; rows = classes, columns = D0, D1, delta, curves UU and LL *)
  threshold = Table[Table[Module[{fs = Table[pp["Classes"][c][s][[k]] /. ppSub /. Global`v -> x, {s, {"UU", "LL"}}]},
      Plot[Evaluate[fs], {x, 1/10, 9/10}, PlotRange -> {{0.1, 0.9}, fRange[fs, x, {1/10, 9/10}]},
        PlotStyle -> Table[Directive[colours[[i]], AbsoluteThickness[2], AbsoluteDashing[dashes[[i]]]], {i, 2}],
        PlotLabel -> classLabel[c],
        FrameLabel -> {Style["v", Italic], Subscript[Style["c", Italic], {"", "D0", "D1", "\[Delta]"}[[k]]]}, Evaluate[Sequence @@ frame]]],
      {k, 2, 4}], {c, sel}];
  exportBoth[FileNameJoin[{plotDir, "pp-threshold"}],
    Legended[GraphicsGrid[threshold, Spacings -> {10, 4}, ImageSize -> 3 * 290],
      Placed[legend[{"UU (spin average)", "LL (helicity difference)"}, 2], Top]]];

  (* (c1) lepton-quark: regular coefficient vs w at v = 1/2: unpolarized, helicity, transverse (Collins) *)
  lqRegular = {lq["Results"]["Unpolarized"]["Coefficients"][[1]], lq["Results"]["Polarized"]["Coefficients"][[1]],
    col["Results"]["Transverse"]["Coefficients"][[1]]} /. lqSub /. Global`v -> 1/2;
  check[Simplify[(col["Results"]["Unpolarized"]["Coefficients"] - lq["Results"]["Unpolarized"]["Coefficients"]) /. lqSub] ===
    {0, 0, 0, 0}, "Collins.wl and LeptonQuark.wl unpolarized coefficients agree"];
  exportBoth[FileNameJoin[{plotDir, "leptonquark-regular-v05"}],
    Legended[Plot[Evaluate[lqRegular /. Global`w -> x], {x, 1/100, 99/100},
      PlotRange -> {{0, 1}, fRange[lqRegular /. Global`w -> x, x, {1/100, 99/100}]},
      PlotStyle -> Table[Directive[colours[[i]], AbsoluteThickness[2], AbsoluteDashing[dashes[[i]]]], {i, 3}],
      FrameLabel -> {Style["w", Italic], Row[{Subscript[Style["c", Italic], "regular"], "(", Style["v", Italic], " = 0.5, ", Style["w", Italic], ")"}]},
      Evaluate[Sequence @@ ReplaceAll[frame, {(ImageSize -> _) -> (ImageSize -> 420), (AspectRatio -> _) -> (AspectRatio -> 0.62)}]]],
      Placed[legend[{"unpolarized", "helicity (LL)", "transverse transfer"}, 3], Top]]];

  (* (c2) spin-transfer ratios at threshold: Born and the coefficients of D0, D1, delta(1-w), polarized / unpolarized *)
  born = {lq["Results"]["Polarized"]["Born"]/lq["Results"]["Unpolarized"]["Born"],
    col["Results"]["Transverse"]["Born"]/col["Results"]["Unpolarized"]["Born"]} // Simplify;
  thr = {lq["Results"]["Polarized"]["Coefficients"][[2 ;; 4]]/lq["Results"]["Unpolarized"]["Coefficients"][[2 ;; 4]],
    col["Results"]["Transverse"]["Coefficients"][[2 ;; 4]]/col["Results"]["Unpolarized"]["Coefficients"][[2 ;; 4]]} /. lqSub;
  Print["Born ratios (helicity, transverse): ", born];
  Print["threshold ratios minus Born ratio (all scale logs 0): ", Simplify[thr - born]];
  check[AllTrue[Flatten[Simplify[thr - born]], # === 0 &], "threshold coefficient ratios equal the Born ratios"];
  vs = Subdivide[1/20, 19/20, 18];
  exportBoth[FileNameJoin[{plotDir, "leptonquark-spin-transfer"}],
    Legended[Show[
      Plot[Evaluate[born /. Global`v -> x], {x, 0, 1}, PlotRange -> {{0, 1}, {0, 1.05}},
        PlotStyle -> Table[Directive[colours[[{2, 3}[[i]]]], AbsoluteThickness[2], AbsoluteDashing[dashes[[{2, 3}[[i]]]]]], {i, 2}],
        FrameLabel -> {Style["v", Italic], "polarized / unpolarized"},
        Evaluate[Sequence @@ ReplaceAll[frame, {(ImageSize -> _) -> (ImageSize -> 420), (AspectRatio -> _) -> (AspectRatio -> 0.62)}]]],
      ListPlot[Flatten[Table[Table[{x, N[thr[[r, k]] /. Global`v -> x]}, {x, vs}], {r, 2}, {k, 3}], 1],
        PlotMarkers -> Flatten[Table[{{"\[EmptyCircle]", 9}, {"\[EmptySquare]", 9}, {"\[EmptyUpTriangle]", 9}}, 2], 1],
        PlotStyle -> Flatten[Table[GrayLevel[0.15], 6]]]],
      Placed[Column[{
        LineLegend[Table[Directive[colours[[i]], AbsoluteThickness[2], AbsoluteDashing[dashes[[i]]]], {i, {2, 3}}],
          {"helicity (LL), Born", "transverse transfer, Born"}, LegendLayout -> "Row", LegendMarkerSize -> {28, 12},
          LabelStyle -> Directive[FontFamily -> "Helvetica", FontSize -> 12, GrayLevel[0.1]]],
        PointLegend[{GrayLevel[0.15], GrayLevel[0.15], GrayLevel[0.15]}, {"D0 coefficient", "D1 coefficient", "\[Delta](1\[Minus]w) coefficient"},
          LegendMarkers -> {{"\[EmptyCircle]", 9}, {"\[EmptySquare]", 9}, {"\[EmptyUpTriangle]", 9}}, LegendLayout -> "Row",
          LabelStyle -> Directive[FontFamily -> "Helvetica", FontSize -> 12, GrayLevel[0.1]]]}], Top]]];
  ];


(* ::Section:: *)
(*Run*)


If[MemberQ[parts, "diagrams"], Print["Diagrams"]; makeDiagrams[]];
If[MemberQ[parts, "plots"], Print["Plots"]; makePlots[]];
Print["done"];
