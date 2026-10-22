(* make_evaluator.wls: write Coefficients/evaluate.py, a stand-alone Python
   (mpmath) evaluator of the coefficients in

     Coefficients/ProtonProton.wl, ProtonProtonTT.wl, ProtonProtonTransfer.wl,
     LeptonQuark.wl, LeptonOffDiagonal.wl, Collins.wl

   with reference values for its built-in self-check.

     wolframscript -file Checks/make_evaluator.wls

   How it works:
     1. every coefficient, Born term and prefactor is translated into a Python
        expression by the recursive translator toPy below (Plus, Times, Power,
        Rational, Integer, Pi, Log of the six arguments v, 1-v, w, 1-w, 1-v w,
        1-v+v w, and the symbols v, w, Nc, Nf, LR, LF, LFF, LI, LO, LL); any other
        head or symbol stops the script, so nothing is translated silently wrong;
     2. the reference values are computed here, from the .wl expressions
        themselves, at exact rational points with N[..., 50], and embedded in
        evaluate.py as strings;
     3. finally `python3 Coefficients/evaluate.py --check` is run (if python3 with
        mpmath is available), which recomputes every reference value from the
        translated expressions at mp.dps = 60 and requires 40-digit agreement.

   Normalization in evaluate.py: the coefficient functions return the bare
   entries {c_regular, c_D0, c_D1, c_delta} of the .wl files (no Prefactor);
   Born terms and prefactors are returned with all couplings = 1 and s = 1,
   with CF = (Nc^2-1)/(2 Nc) and TR = 1/2 in the lepton prefactors.
   Nothing outside Coefficients/evaluate.py is written. *)

$HistoryLength = 0;
root = DirectoryName[ExpandFileName[$InputFileName], 2];
coeffDir = FileNameJoin[{root, "Coefficients"}];
outFile = FileNameJoin[{coeffDir, "evaluate.py"}];
fail[msg__] := (Print["FAILED: ", msg]; Quit[1]);

pp = Get[FileNameJoin[{coeffDir, "ProtonProton.wl"}]];
lq = Get[FileNameJoin[{coeffDir, "LeptonQuark.wl"}]];
lod = Get[FileNameJoin[{coeffDir, "LeptonOffDiagonal.wl"}]];
col = Get[FileNameJoin[{coeffDir, "Collins.wl"}]];
(* the transversely polarized proton-proton files hold the classes collected so far ("NotYetCollected" lists the others) *)
ppTT = Get[FileNameJoin[{coeffDir, "ProtonProtonTT.wl"}]];
ppTr = Get[FileNameJoin[{coeffDir, "ProtonProtonTransfer.wl"}]];

(* ------------------------------------------------------------------------- *)
(* 1. the translator                                                          *)
(* ------------------------------------------------------------------------- *)

symbolNames = <|v -> "v", w -> "w", Nc -> "Nc", Nf -> "Nf", LR -> "LR", LF -> "LF",
   LFF -> "LFF", LI -> "LI", LO -> "LO", LL -> "LL"|>;
(* the logarithms are computed once per function, as local variables *)
logTable = {
   {Log[v], "L_v", "log(v)"},
   {Log[1 - v], "L_1mv", "log(1 - v)"},
   {Log[w], "L_w", "log(w)"},
   {Log[1 - w], "L_1mw", "log(1 - w)"},
   {Log[1 - v w], "L_1mvw", "log(1 - v*w)"},
   {Log[1 - v + v w], "L_1mvpvw", "log(1 - v + v*w)"}};
logNames = Association[Rule @@@ logTable[[All, {1, 2}]]];

negQ[x_?NumberQ] := x < 0;
negQ[Times[c_?NumberQ, ___]] := c < 0;
negQ[_] := False;

(* sums with more terms than this are split into parenthesized groups, which
   keeps the depth of Python's syntax tree small *)
maxTerms = 64;

ClearAll[toPy, factorPy, basePy];
toPy[n_Integer] := ToString[n];
toPy[Rational[a_, b_]] := "R(" <> ToString[a] <> ", " <> ToString[b] <> ")";
toPy[Pi] := "pi";
toPy[s_Symbol] := Lookup[symbolNames, s, fail["unknown symbol ", s]];
toPy[l_Log] := Lookup[logNames, l, fail["unexpected logarithm ", l]];
toPy[e_Plus] := toPySum[List @@ e];
toPySum[terms_List] := If[Length[terms] > maxTerms,
   "(" <> StringRiffle[toPySum /@ Partition[terms, UpTo[maxTerms]], " + "] <> ")",
   "(" <> MapIndexed[
     Function[{t, i}, If[negQ[t], If[i == {1}, "-", " - "] <> toPy[-t],
       If[i == {1}, "", " + "] <> toPy[t]]], terms] <> ")"];
toPy[e_Times] := Module[{c, rest, num = {}, den = {}},
   If[negQ[e], Return["-" <> toPy[-e]]];
   {c, rest} = If[NumberQ[First[e]], {First[e], Rest[List @@ e]}, {1, List @@ e}];
   Which[
    Head[c] === Rational, If[Numerator[c] != 1, AppendTo[num, ToString[Numerator[c]]]];
     AppendTo[den, ToString[Denominator[c]]],
    c =!= 1, AppendTo[num, toPy[c]]];
   Do[If[MatchQ[f, Power[_, _Integer?Negative]],
     AppendTo[den, factorPy[Power[f[[1]], -f[[2]]]]],
     AppendTo[num, factorPy[f]]], {f, rest}];
   If[num === {}, num = {"1"}];
   StringRiffle[num, "*"] <> Which[
     den === {}, "",
     Length[den] == 1, "/" <> First[den],
     True, "/(" <> StringRiffle[den, "*"] <> ")"]];
toPy[Power[b_, n_Integer]] /; n > 0 := basePy[b] <> "**" <> ToString[n];
toPy[Power[b_, n_Integer]] /; n < 0 := "1/" <> toPy[Power[b, -n]];
toPy[Power[b_, r_Rational]] := basePy[b] <> "**" <> toPy[r];
toPy[e_] := fail["cannot translate ", Short[e]];
(* a factor of a product, or the base of a power *)
factorPy[e_Times] := "(" <> toPy[e] <> ")";
factorPy[e_] := toPy[e];
basePy[b_Integer] /; b > 0 := ToString[b];
basePy[b : (_Symbol | _Log | _Plus)] := toPy[b];
basePy[b_] := "(" <> toPy[b] <> ")";

(* wrap a long expression into lines (breaks only at the spaces around + and -,
   all inside the enclosing parentheses, so Python continues the lines) *)
wrap[str_String, indent_String, width_Integer] := Module[{tokens, lines = {}, line = ""},
   (* attach every operator to the operand that follows it *)
   tokens = StringReplace[str, {" + " -> " +\[RawEscape]", " - " -> " -\[RawEscape]"}];
   tokens = StringReplace[StringSplit[tokens, " "], "\[RawEscape]" -> " "];
   Do[If[line === "", line = t,
     If[StringLength[line] + 1 + StringLength[t] > width,
      AppendTo[lines, line]; line = t,
      line = line <> " " <> t]], {t, tokens}];
   AppendTo[lines, line];
   StringRiffle[indent <> # & /@ lines, "\n"]];

pyExpr[e_] := Module[{s = toPy[e]},
   (* a sum is printed with its own parentheses, drop the outer pair *)
   If[Head[e] === Plus && StringStartsQ[s, "("] && StringEndsQ[s, ")"], StringTake[s, {2, -2}], s]];

(* ------------------------------------------------------------------------- *)
(* 2. the results                                                             *)
(* ------------------------------------------------------------------------- *)

couplingsToOne = {alphaS -> 1, alphaEM -> 1, eq -> 1, s -> 1};
colourRules = {CF -> (Nc^2 - 1)/(2 Nc), TR -> 1/2};
basisNames = {"c_regular", "c_D0", "c_D1", "c_delta"};
allowed = Keys[symbolNames];

(* a result: name of the Python function, its arguments, the expressions it
   returns (with their names), and a docstring *)
result[name_, args_, exprs_, names_, doc_] := Module[{free},
   free = Complement[Cases[exprs, _Symbol, {0, Infinity}, Heads -> False], {Pi}];
   If[! SubsetQ[args, free], fail[name, ": free symbols ", Complement[free, args], " not in ", args]];
   <|"name" -> name, "args" -> args, "exprs" -> exprs, "names" -> names, "doc" -> doc|>];

checkFour[name_, c_] := Module[{},
   If[! (ListQ[c] && Length[c] == 4), fail[name, ": not a list of four coefficients"]];
   If[! FreeQ[c[[2 ;; 4]], w], fail[name, ": a distribution coefficient depends on w"]];
   c];

ppSpins = {"UU", "LL"};
ppSpinText = <|"UU" -> "spin average (UU)", "LL" -> "helicity difference (LL)"|>;
ppArgs = {v, w, Nc, Nf, LR, LF, LFF};
lepSpins = {"Unpolarized", "Polarized"};
classes = Keys[pp["Classes"]];
If[classes =!= Range[16], fail["ProtonProton.wl: classes ", classes]];

If[Simplify[pp["Prefactor"] - alphaS^3/(Pi s v w)] =!= 0, fail["unexpected pp prefactor"]];
If[Simplify[ppTT["Prefactor"] - alphaS^3/(Pi s v w)] =!= 0 || Simplify[ppTr["Prefactor"] - alphaS^3/(Pi s v w)] =!= 0,
  fail["unexpected TT or transfer prefactor"]];
ttClasses = Keys[ppTT["Classes"]];
trClasses = Keys[ppTr["Classes"]];
trComps = {"XX", "NN"};
Print["TT classes: ", ttClasses, " (not yet collected: ", ppTT["NotYetCollected"], "); transfer classes: ", trClasses,
  " (not yet collected: ", ppTr["NotYetCollected"], ")"];

results = Join[
   (* proton-proton *)
   Flatten@Table[{
      result["_pp_" <> IntegerString[k, 10, 2] <> "_" <> sp, ppArgs,
       checkFour[{k, sp}, pp["Classes"][k][sp]], basisNames,
       "Class " <> ToString[k] <> ", " <> pp["Classes"][k]["Process"] <> ", " <> ppSpinText[sp] <>
        ": (c_regular, c_D0, c_D1, c_delta)."],
      result["_pp_born_" <> IntegerString[k, 10, 2] <> "_" <> sp, {v, Nc},
       {pp["Classes"][k]["Born"][sp] /. couplingsToOne}, {"born"},
       "Class " <> ToString[k] <> ", " <> pp["Classes"][k]["Process"] <> ", " <> ppSpinText[sp] <>
        ": Born term with s = 1, alphaS = 1."]},
     {k, classes}, {sp, ppSpins}],
   {result["_pp_prefactor", {v, w}, {pp["Prefactor"] /. couplingsToOne}, {"prefactor"},
     "Proton-proton prefactor alphaS^3/(pi s v w) at s = 1, alphaS = 1."]},
   (* p^ p^ -> c + X, cos(2 Phi) coefficient *)
   Flatten@Table[{
      result["_tt_" <> IntegerString[k, 10, 2], ppArgs,
       checkFour[{"TT", k}, ppTT["Classes"][k]["TT"]], basisNames,
       "TT class " <> ToString[k] <> ", " <> ppTT["Classes"][k]["Process"] <>
        ", cos(2 Phi) coefficient of the double-transverse cross section: (c_regular, c_D0, c_D1, c_delta)."],
      result["_tt_born_" <> IntegerString[k, 10, 2], {v, Nc},
       {ppTT["Classes"][k]["Born"] /. couplingsToOne}, {"born"},
       "TT class " <> ToString[k] <> ", " <> ppTT["Classes"][k]["Process"] <> ": Born term with s = 1, alphaS = 1."]},
     {k, ttClasses}],
   (* a^ b -> c^ + X, components XX and NN *)
   Flatten@Table[{
      result["_tr_" <> IntegerString[k, 10, 2] <> "_" <> cp, ppArgs,
       checkFour[{"Transfer", k, cp}, ppTr["Classes"][k][cp]], basisNames,
       "Transfer class " <> ToString[k] <> ", " <> ppTr["Classes"][k]["Process"] <> ", component " <> cp <>
        ": (c_regular, c_D0, c_D1, c_delta)."],
      result["_tr_born_" <> IntegerString[k, 10, 2] <> "_" <> cp, {v, Nc},
       {ppTr["Classes"][k]["Born"][cp] /. couplingsToOne}, {"born"},
       "Transfer class " <> ToString[k] <> ", " <> ppTr["Classes"][k]["Process"] <> ", component " <> cp <>
        ": Born term with s = 1, alphaS = 1."]},
     {k, trClasses}, {cp, trComps}],
   (* l q -> q *)
   Flatten@Table[{
      result["_lq_" <> sp, {v, w, LR},
       checkFour[{"LeptonQuark", sp}, lq["Results"][sp]["Coefficients"]], basisNames,
       "l q -> q + X, " <> sp <> ": (c_regular, c_D0, c_D1, c_delta)."],
      result["_lq_born_" <> sp, {v}, {lq["Results"][sp]["Born"] /. couplingsToOne}, {"born"},
       "l q -> q + X, " <> sp <> ": Born term with s = 1, alphaEM = eq = 1."],
      result["_lq_prefactor_" <> sp, {v, w, Nc},
       {lq["Results"][sp]["Prefactor"] /. couplingsToOne /. colourRules}, {"prefactor"},
       "l q -> q + X, " <> sp <> ": prefactor at s = 1, alphaEM = alphaS = eq = 1, CF = (Nc^2-1)/(2 Nc)."]},
     {sp, lepSpins}],
   (* off-diagonal lepton channels *)
   Flatten@Table[{
      result["_lod_" <> ch <> "_" <> sp, {v, w, LR},
       checkFour[{"LeptonOffDiagonal", ch, sp}, lod["Results"][sp][ch]["Coefficients"]], basisNames,
       Lookup[<|"QG" -> "l q -> g + X", "GQ" -> "l g -> q + X"|>, ch] <> " (" <> ch <> "), " <> sp <>
        ": (c_regular, c_D0, c_D1, c_delta)."],
      result["_lod_prefactor_" <> ch <> "_" <> sp, {v, w, Nc},
       {lod["Results"][sp][ch]["Prefactor"] /. couplingsToOne /. colourRules}, {"prefactor"},
       ch <> ", " <> sp <> ": prefactor at s = 1, alphaEM = alphaS = eq = 1, CF = (Nc^2-1)/(2 Nc), TR = 1/2."]},
     {sp, lepSpins}, {ch, {"QG", "GQ"}}],
   (* Collins *)
   Flatten@Table[{
      result["_collins_" <> pt, {v, w, LR, LI, LO, LL},
       checkFour[{"Collins", pt}, col["Results"][pt]["Coefficients"]], basisNames,
       "l q -> q + X, Collins.wl part " <> pt <> ": (c_regular, c_D0, c_D1, c_delta)."],
      result["_collins_born_" <> pt, {v}, {col["Results"][pt]["Born"] /. couplingsToOne}, {"born"},
       "l q -> q + X, Collins.wl part " <> pt <> ": Born term with s = 1, alphaEM = eq = 1."],
      result["_collins_prefactor_" <> pt, {v, w, Nc},
       {col["Results"][pt]["Prefactor"] /. couplingsToOne /. colourRules}, {"prefactor"},
       pt <> ": prefactor at s = 1, alphaEM = alphaS = eq = 1, CF = (Nc^2-1)/(2 Nc)."]},
     {pt, {"Unpolarized", "Transverse"}}]];

If[Keys[lod["Results"]] =!= lepSpins || Keys[lq["Results"]] =!= lepSpins, fail["unexpected lepton keys"]];
If[Sort[Keys[col["Results"]]] =!= Sort[{"Unpolarized", "Transverse"}], fail["unexpected Collins keys"]];
Print[Length[results], " Python functions to write"];

pyFunction[r_] := Module[{exprs = r["exprs"], logs, body},
   logs = Select[logTable, ! FreeQ[exprs, #[[1]]] &];
   body = Join[
     {"    '''" <> r["doc"] <> "'''"},
     ("    " <> #[[2]] <> " = " <> #[[3]]) & /@ logs,
     MapThread[
      Function[{nm, e}, Module[{s = pyExpr[e]},
        If[StringLength[s] <= 80, "    " <> nm <> " = " <> s,
         "    " <> nm <> " = (\n" <> wrap[s, "        ", 96] <> "\n    )"]]],
      {r["names"], exprs}],
     {"    return " <> If[Length[exprs] == 1, First[r["names"]], StringRiffle[r["names"], ", "]]}];
   "def " <> r["name"] <> "(" <> StringRiffle[symbolNames /@ r["args"], ", "] <> "):\n" <>
    StringRiffle[body, "\n"] <> "\n"];

t0 = AbsoluteTime[];
generated = StringRiffle[pyFunction /@ results, "\n\n"];
Print["translated in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* ------------------------------------------------------------------------- *)
(* 3. reference values                                                        *)
(* ------------------------------------------------------------------------- *)

points = {
   <|v -> 3/10, w -> 1/5, Nc -> 3, Nf -> 5, LR -> 0, LF -> 0, LFF -> 0, LI -> 0, LO -> 0, LL -> 0|>,
   <|v -> 7/11, w -> 5/9, Nc -> 3, Nf -> 5, LR -> 1/3, LF -> -2/7, LFF -> 3/5, LI -> -2/7, LO -> 3/5, LL -> -1/4|>,
   <|v -> 1/7, w -> 13/14, Nc -> 4, Nf -> 3, LR -> -1/2, LF -> 1/4, LFF -> -1/5, LI -> 1/4, LO -> -1/5, LL -> 2/3|>};

digitString[x_] := Module[{d, e},
   If[x == 0, Return["0"]];
   {d, e} = RealDigits[x, 10, 50];
   If[x < 0, "-", ""] <> "0." <> StringJoin[ToString /@ d] <> "e" <> ToString[e]];
refValue[e_, p_] := Module[{x = e /. Normal[p]},
   If[! FreeQ[x, Alternatives @@ allowed], fail["point does not fix all symbols"]];
   If[x === 0, "0", digitString[N[x, 50]]]];

pyRational[q_] := "'" <> ToString[q, InputForm] <> "'";
pointsPy = "_POINTS = [\n" <> StringRiffle[Function[p,
     "    {" <> StringRiffle[(("'" <> symbolNames[#] <> "': " <> pyRational[p[#]]) & /@ Keys[p]), ", "] <> "},"] /@
    points, "\n"] <> "\n]\n";

t0 = AbsoluteTime[];
referencePy = "_REFERENCE = [\n" <> StringRiffle[Flatten@Table[
      "    ('" <> r["name"] <> "', " <> ToString[ip - 1] <> ", (" <>
       If[Length[r["exprs"]] == 1, "'" <> refValue[First[r["exprs"]], points[[ip]]] <> "',",
        "\n" <> StringRiffle[("        '" <> refValue[#, points[[ip]]] <> "'") & /@ r["exprs"], ",\n"]] <> ")),",
      {r, results}, {ip, Length[points]}], "\n"] <> "\n]\n";
nChecks = Length[points] Total[Length[#["exprs"]] & /@ results];
Print[nChecks, " reference values computed in ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* ------------------------------------------------------------------------- *)
(* 4. the module                                                              *)
(* ------------------------------------------------------------------------- *)

header = "#!/usr/bin/env python3
'''
evaluate.py -- stand-alone numerical evaluator of the NLO coefficients in Coefficients/

GENERATED by Checks/make_evaluator.wls from ProtonProton.wl, ProtonProtonTT.wl,
ProtonProtonTransfer.wl, LeptonQuark.wl, LeptonOffDiagonal.wl and Collins.wl.  Do not edit by hand: rerun

    wolframscript -file Checks/make_evaluator.wls

Requires only mpmath (Python 3.8 or later).  All numbers are mpmath mpf at the
current precision mp.dps (set mpmath.mp.dps before calling, default 15).

WHAT THE NUMBERS ARE (see Coefficients/README.md for the full formula)

The partonic cross section sigma = s E dsigma/d^3p is written in v = 1 + t/s,
w = -u/(s + t), 0 < v, w < 1, as

    sigma = Born delta(1-w)
          + P(v, w) c_regular(v, w)
          + P(v, 1) [ c_D0(v) D0 + c_D1(v) D1 + c_delta(v) delta(1-w) ],
    D0 = [1/(1-w)]_+ ,  D1 = [log(1-w)/(1-w)]_+ .

The coefficient functions (pp, tt, transfer, lepton_quark, lepton_offdiagonal,
collins) return
a dict with the keys 'regular', 'D0', 'D1', 'delta': exactly the four entries
{c_regular, c_D0, c_D1, c_delta} of the .wl files.  They do NOT include the
prefactor P.  The prefactor and the Born term are separate functions, evaluated
at s = 1 with all couplings (alphaS, alphaEM, eq) set to 1:

    pp_prefactor(v, w)        = 1/(pi v w)       [P = alphaS^3/(pi s v w)]
    born(cls, spin, v, Nc)    = Born * s/alphaS^2   (also tt_born, transfer_born)
    lepton_*_prefactor        = CF/(pi v w) or TR/(pi v w), CF = (Nc^2-1)/(2 Nc), TR = 1/2
                                [P = alphaEM^2 alphaS eq^2 CF/(pi s v w), CF -> TR for GQ]
    lepton_quark_born, collins_born = Born * s/(alphaEM^2 eq^2)

so that, e.g., the physical proton-proton Born term is born(...) * alphaS^2 / s
and the physical regular part is pp_prefactor(v, w) * alphaS^3 / s * c_regular.
c_D0, c_D1, c_delta do not depend on w (checked when this file was generated);
the w passed to the coefficient functions matters only for c_regular, which needs
0 < w < 1.

Scales: LR = log(muR^2/s), LF = log(muF^2/s), LFF = log(muFF^2/s) in pp, tt and
transfer (LF for both incoming partons, LFF for the fragmenting parton); the lepton
results use one common scale LR; collins() takes LR and LI, LO, LL = log(mu^2/s)
for the incoming quark, the outgoing quark and the photon-in-lepton subtraction.
Colour SU(Nc), CA = Nc, CF = (Nc^2-1)/(2 Nc); Nf light flavours.  MSbar; helicity
differences in the HVBM scheme with the helicity-restoring finite renormalization;
tt and transfer in the HVBM scheme with four-dimensional spin vectors (no finite
renormalization).  tt is the coefficient of cos(2 Phi) of the double-transverse
cross section (1/2)[sigma(up,up) - sigma(up,down)], Processes/HadronicTT;
transfer(cls, 'XX' or 'NN', ...) are the in-plane and normal components of the
transverse spin transfer a^ b -> c^ + X, Processes/HadronicTransfer (equal for
the finite coefficients).  Only the classes collected so far are present
(tt_classes(), transfer_classes()).

Arguments may be numbers, mpf, fractions.Fraction or strings such as '0.3' or
'3/10' (a string is converted at the current precision; a Python float such as
0.3 is taken with its binary value).

FUNCTIONS
    pp(cls, spin, v, w, Nc=3, Nf=5, LR=0, LF=0, LFF=0)   cls = 1..16, spin = 'UU' or 'LL'
    born(cls, spin, v, Nc=3),  pp_prefactor(v, w),  pp_process(cls)
    tt(cls, v, w, Nc=3, Nf=5, LR=0, LF=0, LFF=0)          cls in tt_classes() (subset of 5, 6, 7, 11, 12)
    tt_born(cls, v, Nc=3),  tt_process(cls),  tt_classes()
    transfer(cls, comp, v, w, Nc=3, Nf=5, LR=0, LF=0, LFF=0)   cls in transfer_classes() (subset of
                                                          1, 3, 6, 10, 11, 13, 17), comp = 'XX' or 'NN'
    transfer_born(cls, comp, v, Nc=3),  transfer_process(cls),  transfer_classes()
    lepton_quark(spin, v, w, LR=0)                        spin = 'Unpolarized'/'UU' or 'Polarized'/'LL'
    lepton_quark_born(spin, v),  lepton_quark_prefactor(spin, v, w, Nc=3)
    lepton_offdiagonal(channel, spin, v, w, LR=0)         channel = 'QG' (l q -> g) or 'GQ' (l g -> q)
    lepton_offdiagonal_prefactor(channel, spin, v, w, Nc=3)
    collins(part, v, w, LR=0, LI=0, LO=0, LL=0)           part = 'Transverse' or 'Unpolarized'
    collins_born(part, v),  collins_prefactor(part, v, w, Nc=3)
    run_checks(verbose=False)

COMMAND LINE
    python3 evaluate.py pp 13 LL 0.3 0.2 --digits 30 [--Nc 3 --Nf 5 --LR 0 --LF 0 --LFF 0]
    python3 evaluate.py tt 5 0.3 0.2 [same options as pp]
    python3 evaluate.py transfer 1 XX 0.3 0.2 [same options as pp]
    python3 evaluate.py lepton_quark Polarized 0.3 0.2 [--LR 0]
    python3 evaluate.py lepton_offdiagonal QG Unpolarized 0.3 0.2 [--LR 0]
    python3 evaluate.py collins Transverse 0.3 0.2 [--LR 0 --LI 0 --LO 0 --LL 0]
    python3 evaluate.py --check [--verbose]

SELF-CHECK
    --check recomputes, at mp.dps = 60, every coefficient, Born term and prefactor
    at three exact rational points (with nonzero scale logarithms and Nc = 3, 4)
    and compares with 50-digit reference values that make_evaluator.wls computed
    directly from the .wl files; agreement to 40 significant digits is required
    (absolute 1e-40 where the reference is exactly zero).  The reference values
    come from the same saved .wl expressions, so this check tests the translation
    of those expressions into Python (and the arithmetic of this module), not the
    physics: the physics checks are those of the process notebooks
    (Processes/<process>/) and of Checks/collect_coefficients.wls.
'''

import sys
import argparse
from fractions import Fraction
from mpmath import mp, mpf, log, pi, nstr

__all__ = ['pp', 'born', 'pp_prefactor', 'pp_process',
           'tt', 'tt_born', 'tt_process', 'tt_classes',
           'transfer', 'transfer_born', 'transfer_process', 'transfer_classes',
           'lepton_quark', 'lepton_quark_born', 'lepton_quark_prefactor',
           'lepton_offdiagonal', 'lepton_offdiagonal_prefactor',
           'collins', 'collins_born', 'collins_prefactor', 'run_checks', 'BASIS']

BASIS = ('regular', 'D0', 'D1', 'delta')


def R(a, b):
    '''The rational number a/b at the current precision.'''
    return mpf(a) / b


def _num(x):
    '''Convert an argument (number, mpf, Fraction or string, e.g. 3/10) to mpf.'''
    if isinstance(x, Fraction):
        return mpf(x.numerator) / x.denominator
    if isinstance(x, str):
        s = x.strip()
        if '/' in s:
            a, b = s.split('/')
            return mpf(a.strip()) / mpf(b.strip())
        return mpf(s)
    return mpf(x)


# ---------------------------------------------------------------------------
# generated functions: one per result, the coefficients as separate expressions
# ---------------------------------------------------------------------------

";

api = "

# ---------------------------------------------------------------------------
# public interface
# ---------------------------------------------------------------------------

_PP_PROCESS = @@PROCESSES@@

_TT_PROCESS = @@TTPROCESSES@@

_TR_PROCESS = @@TRPROCESSES@@

_TR_COMP = {'XX': 'XX', 'NN': 'NN'}

_PP_SPIN = {'UU': 'UU', 'LL': 'LL', 'UNPOLARIZED': 'UU', 'POLARIZED': 'LL'}
_LEP_SPIN = {'UNPOLARIZED': 'Unpolarized', 'UU': 'Unpolarized',
             'POLARIZED': 'Polarized', 'LL': 'Polarized'}
_CHANNEL = {'QG': 'QG', 'GQ': 'GQ'}
_PART = {'TRANSVERSE': 'Transverse', 'UNPOLARIZED': 'Unpolarized'}


def _key(table, x, what):
    try:
        return table[str(x).upper()]
    except KeyError:
        raise ValueError('unknown %s %r (allowed: %s)' % (what, x, ', '.join(sorted(set(table.values())))))


def _cls(cls):
    k = int(cls)
    if k not in _PP_PROCESS:
        raise ValueError('class must be 1..16, got %r' % (cls,))
    return '%02d' % k


def _four(vals):
    return dict(zip(BASIS, vals))


def pp_process(cls):
    '''The partonic channel of proton-proton class cls (1..16), e.g. q g -> q.'''
    return _PP_PROCESS[int(cls)]


def pp(cls, spin, v, w, Nc=3, Nf=5, LR=0, LF=0, LFF=0):
    '''Proton-proton class cls (1..16), spin UU or LL: dict of c_regular, c_D0,
    c_D1, c_delta (without the prefactor alphaS^3/(pi s v w)).'''
    f = globals()['_pp_' + _cls(cls) + '_' + _key(_PP_SPIN, spin, 'spin')]
    return _four(f(_num(v), _num(w), _num(Nc), _num(Nf), _num(LR), _num(LF), _num(LFF)))


def born(cls, spin, v, Nc=3):
    '''Born term of proton-proton class cls at s = 1, alphaS = 1.'''
    f = globals()['_pp_born_' + _cls(cls) + '_' + _key(_PP_SPIN, spin, 'spin')]
    return f(_num(v), _num(Nc))


def pp_prefactor(v, w):
    '''Proton-proton prefactor at s = 1, alphaS = 1: 1/(pi v w).'''
    return _pp_prefactor(_num(v), _num(w))


def tt_classes():
    '''The TT classes present in this file (collected from all three routes).'''
    return sorted(_TT_PROCESS)


def tt_process(cls):
    '''The partonic channel of TT class cls, e.g. q qb -> qp.'''
    return _TT_PROCESS[int(cls)]


def _tt_cls(cls):
    k = int(cls)
    if k not in _TT_PROCESS:
        raise ValueError('TT class must be one of %s, got %r' % (tt_classes(), cls))
    return '%02d' % k


def tt(cls, v, w, Nc=3, Nf=5, LR=0, LF=0, LFF=0):
    '''p^ p^ -> c + X, TT class cls: dict of c_regular, c_D0, c_D1, c_delta of the
    cos(2 Phi) coefficient (without the prefactor alphaS^3/(pi s v w)).'''
    f = globals()['_tt_' + _tt_cls(cls)]
    return _four(f(_num(v), _num(w), _num(Nc), _num(Nf), _num(LR), _num(LF), _num(LFF)))


def tt_born(cls, v, Nc=3):
    '''Born term of TT class cls at s = 1, alphaS = 1 (cos(2 Phi) coefficient).'''
    return globals()['_tt_born_' + _tt_cls(cls)](_num(v), _num(Nc))


def transfer_classes():
    '''The transfer classes present in this file (collected from all three routes).'''
    return sorted(_TR_PROCESS)


def transfer_process(cls):
    '''The partonic channel of transfer class cls, e.g. q qp -> q (first parton
    polarized, observed parton after the arrow).'''
    return _TR_PROCESS[int(cls)]


def _tr_cls(cls):
    k = int(cls)
    if k not in _TR_PROCESS:
        raise ValueError('transfer class must be one of %s, got %r' % (transfer_classes(), cls))
    return '%02d' % k


def transfer(cls, comp, v, w, Nc=3, Nf=5, LR=0, LF=0, LFF=0):
    '''a^ b -> c^ + X, transfer class cls, component XX (in plane) or NN (normal):
    dict of c_regular, c_D0, c_D1, c_delta (without the prefactor alphaS^3/(pi s v w)).'''
    f = globals()['_tr_' + _tr_cls(cls) + '_' + _key(_TR_COMP, comp, 'component')]
    return _four(f(_num(v), _num(w), _num(Nc), _num(Nf), _num(LR), _num(LF), _num(LFF)))


def transfer_born(cls, comp, v, Nc=3):
    '''Born term of transfer class cls, component XX or NN, at s = 1, alphaS = 1.'''
    return globals()['_tr_born_' + _tr_cls(cls) + '_' + _key(_TR_COMP, comp, 'component')](_num(v), _num(Nc))


def lepton_quark(spin, v, w, LR=0):
    '''l q -> q + X, Unpolarized or Polarized: dict of the four coefficients
    (without the prefactor alphaEM^2 alphaS CF eq^2/(pi s v w)).'''
    f = globals()['_lq_' + _key(_LEP_SPIN, spin, 'spin')]
    return _four(f(_num(v), _num(w), _num(LR)))


def lepton_quark_born(spin, v):
    '''Born term of l q -> q + X at s = 1, alphaEM = eq = 1.'''
    return globals()['_lq_born_' + _key(_LEP_SPIN, spin, 'spin')](_num(v))


def lepton_quark_prefactor(spin, v, w, Nc=3):
    '''Prefactor of l q -> q + X at s = 1 and unit couplings: CF/(pi v w).'''
    return globals()['_lq_prefactor_' + _key(_LEP_SPIN, spin, 'spin')](_num(v), _num(w), _num(Nc))


def lepton_offdiagonal(channel, spin, v, w, LR=0):
    '''l q -> g + X (QG) or l g -> q + X (GQ), Unpolarized or Polarized: dict of
    the four coefficients (without the prefactor; no Born term in these channels).'''
    f = globals()['_lod_' + _key(_CHANNEL, channel, 'channel') + '_' + _key(_LEP_SPIN, spin, 'spin')]
    return _four(f(_num(v), _num(w), _num(LR)))


def lepton_offdiagonal_prefactor(channel, spin, v, w, Nc=3):
    '''Prefactor of the off-diagonal lepton channels at s = 1 and unit couplings:
    CF/(pi v w) for QG, TR/(pi v w) with TR = 1/2 for GQ.'''
    f = globals()['_lod_prefactor_' + _key(_CHANNEL, channel, 'channel') + '_' + _key(_LEP_SPIN, spin, 'spin')]
    return f(_num(v), _num(w), _num(Nc))


def collins(part, v, w, LR=0, LI=0, LO=0, LL=0):
    '''Collins.wl, part Transverse or Unpolarized: dict of the four coefficients
    (without the prefactor alphaEM^2 alphaS CF eq^2/(pi s v w)).'''
    f = globals()['_collins_' + _key(_PART, part, 'part')]
    return _four(f(_num(v), _num(w), _num(LR), _num(LI), _num(LO), _num(LL)))


def collins_born(part, v):
    '''Born term of Collins.wl part Transverse or Unpolarized at s = 1, alphaEM = eq = 1.'''
    return globals()['_collins_born_' + _key(_PART, part, 'part')](_num(v))


def collins_prefactor(part, v, w, Nc=3):
    '''Prefactor of Collins.wl part at s = 1 and unit couplings: CF/(pi v w).'''
    return globals()['_collins_prefactor_' + _key(_PART, part, 'part')](_num(v), _num(w), _num(Nc))


# ---------------------------------------------------------------------------
# self-check against reference values computed by make_evaluator.wls from the
# .wl files at exact rational points, with 50 significant digits
# ---------------------------------------------------------------------------

";

checks = "

def run_checks(verbose=False):
    '''Compare every generated function with the embedded reference values at
    mp.dps = 60, requiring 40 significant digits.  Returns True if all pass.'''
    saved = mp.dps
    mp.dps = 60
    try:
        tol = mpf(10) ** -40
        npass, failures, worst = 0, [], mpf(0)
        for name, ip, refs in _REFERENCE:
            f = globals()[name]
            argnames = f.__code__.co_varnames[:f.__code__.co_argcount]
            vals = f(*[_num(_POINTS[ip][a]) for a in argnames])
            if not isinstance(vals, tuple):
                vals = (vals,)
            if len(vals) != len(refs):
                failures.append((name, ip, 'wrong number of values'))
                continue
            for i, (val, ref) in enumerate(zip(vals, refs)):
                r = mpf(ref)
                dev = abs(val - r) / abs(r) if r != 0 else abs(val)
                worst = max(worst, dev)
                if dev <= tol:
                    npass += 1
                else:
                    failures.append((name, ip, i, nstr(val, 45), ref, nstr(dev, 3)))
                if verbose:
                    print('%-30s point %d value %d: deviation %s' % (name, ip, i, nstr(dev, 3)))
        for fl in failures:
            print('FAILED:', *fl)
        print('%d checks passed, %d failed (%d functions x %d points; worst relative deviation %s)'
              % (npass, len(failures), len(set(r[0] for r in _REFERENCE)), len(_POINTS), nstr(worst, 3)))
        return not failures
    finally:
        mp.dps = saved


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------

def _print(title, values, digits):
    print(title)
    for k, x in values:
        print('  %-10s = %s' % (k, nstr(x, digits)))


def main(argv=None):
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--digits', type=int, default=16, help='significant digits printed (default 16)')
    p = argparse.ArgumentParser(description='Evaluate the NLO coefficients of Coefficients/ (see the module docstring).')
    p.add_argument('--check', action='store_true', help='run the built-in checks')
    p.add_argument('--verbose', action='store_true', help='with --check: print every comparison')
    sub = p.add_subparsers(dest='cmd')
    a = sub.add_parser('pp', help='proton-proton class: pp CLS SPIN v w', parents=[common])
    a.add_argument('cls', type=int)
    a.add_argument('spin')
    a.add_argument('v')
    a.add_argument('w')
    for o, d in (('--Nc', '3'), ('--Nf', '5'), ('--LR', '0'), ('--LF', '0'), ('--LFF', '0')):
        a.add_argument(o, default=d)
    a = sub.add_parser('tt', help='p^ p^ -> c, cos(2 Phi) coefficient: tt CLS v w', parents=[common])
    a.add_argument('cls', type=int)
    a.add_argument('v')
    a.add_argument('w')
    for o, d in (('--Nc', '3'), ('--Nf', '5'), ('--LR', '0'), ('--LF', '0'), ('--LFF', '0')):
        a.add_argument(o, default=d)
    a = sub.add_parser('transfer', help='a^ b -> c^ spin transfer: transfer CLS XX|NN v w', parents=[common])
    a.add_argument('cls', type=int)
    a.add_argument('comp')
    a.add_argument('v')
    a.add_argument('w')
    for o, d in (('--Nc', '3'), ('--Nf', '5'), ('--LR', '0'), ('--LF', '0'), ('--LFF', '0')):
        a.add_argument(o, default=d)
    a = sub.add_parser('lepton_quark', aliases=['lq'], help='l q -> q: lepton_quark SPIN v w', parents=[common])
    a.add_argument('spin')
    a.add_argument('v')
    a.add_argument('w')
    a.add_argument('--LR', default='0')
    a.add_argument('--Nc', default='3', help='only for the printed prefactor')
    a = sub.add_parser('lepton_offdiagonal', aliases=['lod'], help='QG / GQ: lepton_offdiagonal CHANNEL SPIN v w',
                       parents=[common])
    a.add_argument('channel')
    a.add_argument('spin')
    a.add_argument('v')
    a.add_argument('w')
    a.add_argument('--LR', default='0')
    a.add_argument('--Nc', default='3', help='only for the printed prefactor')
    a = sub.add_parser('collins', help='Collins.wl: collins PART v w', parents=[common])
    a.add_argument('part')
    a.add_argument('v')
    a.add_argument('w')
    for o in ('--LR', '--LI', '--LO', '--LL'):
        a.add_argument(o, default='0')
    a.add_argument('--Nc', default='3', help='only for the printed prefactor')
    args = p.parse_args(argv)

    if args.check:
        return 0 if run_checks(args.verbose) else 1
    if args.cmd is None:
        p.print_help()
        return 2
    mp.dps = max(args.digits, 15) + 10
    if args.cmd == 'pp':
        c = pp(args.cls, args.spin, args.v, args.w, args.Nc, args.Nf, args.LR, args.LF, args.LFF)
        title = ('pp class %d (%s), %s, v = %s, w = %s, Nc = %s, Nf = %s, LR = %s, LF = %s, LFF = %s'
                 % (args.cls, pp_process(args.cls), _key(_PP_SPIN, args.spin, 'spin'), args.v, args.w,
                    args.Nc, args.Nf, args.LR, args.LF, args.LFF))
        extra = [('Born', born(args.cls, args.spin, args.v, args.Nc)),
                 ('prefactor', pp_prefactor(args.v, args.w))]
        note = '  (Born and prefactor at s = 1, alphaS = 1)'
    elif args.cmd == 'tt':
        c = tt(args.cls, args.v, args.w, args.Nc, args.Nf, args.LR, args.LF, args.LFF)
        title = ('TT class %d (%s), cos(2 Phi) coefficient, v = %s, w = %s, Nc = %s, Nf = %s, LR = %s, LF = %s, LFF = %s'
                 % (args.cls, tt_process(args.cls), args.v, args.w, args.Nc, args.Nf, args.LR, args.LF, args.LFF))
        extra = [('Born', tt_born(args.cls, args.v, args.Nc)), ('prefactor', pp_prefactor(args.v, args.w))]
        note = '  (Born and prefactor at s = 1, alphaS = 1)'
    elif args.cmd == 'transfer':
        c = transfer(args.cls, args.comp, args.v, args.w, args.Nc, args.Nf, args.LR, args.LF, args.LFF)
        title = ('transfer class %d (%s), component %s, v = %s, w = %s, Nc = %s, Nf = %s, LR = %s, LF = %s, LFF = %s'
                 % (args.cls, transfer_process(args.cls), _key(_TR_COMP, args.comp, 'component'), args.v, args.w,
                    args.Nc, args.Nf, args.LR, args.LF, args.LFF))
        extra = [('Born', transfer_born(args.cls, args.comp, args.v, args.Nc)), ('prefactor', pp_prefactor(args.v, args.w))]
        note = '  (Born and prefactor at s = 1, alphaS = 1)'
    elif args.cmd in ('lepton_quark', 'lq'):
        c = lepton_quark(args.spin, args.v, args.w, args.LR)
        title = ('l q -> q + X, %s, v = %s, w = %s, LR = %s'
                 % (_key(_LEP_SPIN, args.spin, 'spin'), args.v, args.w, args.LR))
        extra = [('Born', lepton_quark_born(args.spin, args.v)),
                 ('prefactor', lepton_quark_prefactor(args.spin, args.v, args.w, args.Nc))]
        note = '  (Born and prefactor at s = 1, alphaEM = alphaS = eq = 1)'
    elif args.cmd in ('lepton_offdiagonal', 'lod'):
        c = lepton_offdiagonal(args.channel, args.spin, args.v, args.w, args.LR)
        title = ('lepton off-diagonal %s, %s, v = %s, w = %s, LR = %s'
                 % (_key(_CHANNEL, args.channel, 'channel'), _key(_LEP_SPIN, args.spin, 'spin'),
                    args.v, args.w, args.LR))
        extra = [('prefactor', lepton_offdiagonal_prefactor(args.channel, args.spin, args.v, args.w, args.Nc))]
        note = '  (prefactor at s = 1, alphaEM = alphaS = eq = 1; no Born term)'
    else:
        c = collins(args.part, args.v, args.w, args.LR, args.LI, args.LO, args.LL)
        title = ('Collins %s, v = %s, w = %s, LR = %s, LI = %s, LO = %s, LL = %s'
                 % (_key(_PART, args.part, 'part'), args.v, args.w, args.LR, args.LI, args.LO, args.LL))
        extra = [('Born', collins_born(args.part, args.v)),
                 ('prefactor', collins_prefactor(args.part, args.v, args.w, args.Nc))]
        note = '  (Born and prefactor at s = 1, alphaEM = alphaS = eq = 1)'
    _print(title, [('c_' + k, c[k]) for k in BASIS] + extra, args.digits)
    print(note)
    return 0


if __name__ == '__main__':
    sys.exit(main())
";

processDict[assoc_, keys_] := If[keys === {}, "{}", "{\n" <> StringRiffle[("    " <> ToString[#] <> ": '" <> assoc[#]["Process"] <> "',") & /@ keys, "\n"] <> "\n}"];
processesPy = processDict[pp["Classes"], classes];
module = header <> generated <> StringReplace[api, {"@@PROCESSES@@" -> processesPy,
     "@@TTPROCESSES@@" -> processDict[ppTT["Classes"], ttClasses], "@@TRPROCESSES@@" -> processDict[ppTr["Classes"], trClasses]}] <>
   pointsPy <> "\n" <> referencePy <> checks;

stream = OpenWrite[outFile, CharacterEncoding -> "UTF-8"];
WriteString[stream, module];
Close[stream];
Print["wrote ", outFile, " (", FileByteCount[outFile], " bytes, ", Length[StringSplit[module, "\n"]], " lines)"];

(* ------------------------------------------------------------------------- *)
(* 5. run the self-check                                                      *)
(* ------------------------------------------------------------------------- *)

python = Environment["PYTHON"]; If[! StringQ[python], python = "python3"];
run = RunProcess[{python, outFile, "--check"}];
If[run === $Failed || run["ExitCode"] =!= 0 && StringContainsQ[run["StandardError"], "No module named 'mpmath'"],
  Print["(could not run the self-check: ", python, " with mpmath not found; set PYTHON=...)"],
  Print[StringTrim[run["StandardOutput"]], StringTrim[run["StandardError"]]];
  If[run["ExitCode"] =!= 0, fail["evaluate.py --check"]]];
