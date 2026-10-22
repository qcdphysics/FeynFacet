(* ::Title:: *)
(*Proton-proton classes, step 2 (FORM route): spin and colour algebra*)


(* ::Text:: *)
(*Input: the amplitude records of step 1 for the class chosen below.  Output: Results/FORM/class-<n>/Algebra.wl, with the same content, conventions and format as the file of the FeynCalc route (2-Algebra-FeynCalc.wl):*)
(*  Born:     the Born square |M|^2 (UU and LL), summed over final spins and colours, averaged over the initial ones, in D dimensions;*)
(*  Virtual:  the interference 2 Re(M_1loop M_Born^*) (UU and LL), averaged and summed in the same way, written as Sum_j [i c_j M0*] I_j in terms of the scalar one-loop integrals of the families VirtualST, VirtualSU, VirtualTU, with the powers of the (D-4)-dimensional part of the loop momentum as HatPower[n];*)
(*  Real:     the sum over the real-emission processes of the class, each weighted as in classes.wl, of |M|^2 (UU and LL) in the invariants s, t, u, A, B, C and the hatted invariant spe3.*)
(*The last section checks that the file is identical to the one of the FeynCalc route.*)
(**)
(*Spin sums (HVBM scheme for LL), as in the FeynCalc route: the physical momenta pa, pb, p are four-dimensional; the unobserved momenta k1, k2 and the loop momentum have (D-4)-dimensional parts.  Incoming quark: q-slash (UU) or sign(q) q-slash gamma5 (LL, sign -1 for an antiquark), with the average 1/2 applied separately; incoming gluon: physical polarization sum with the average 1/(D-2) (UU), or the circular density i eps^{mu nu l r}/(2 l.r) (LL), which contains the average 1/2 itself; outgoing gluons: physical polarization sum -g^{mu nu} + (l^mu r^nu + r^mu l^nu)/(l.r) with a reference vector r.  The LL result is (sigma(++) - sigma(+-))/2 for the helicities of the incoming partons.*)
(**)
(*Division of the work.  Mathematica reads the amplitude records of step 1, applies momentum conservation, splits every diagram into its colour structures and Lorentz parts and writes it as a FORM expression; it starts the FORM jobs, reads their results and puts in the kinematic values of the propagator denominators.  FORM does all the spin and colour algebra.  The FORM programs are in the directory form/ next to this file and are meant to be read:*)
(*  form/common.h:   declarations; the spin sums and densities of the quarks (with the HVBM gamma5); the gluon polarization sums and circular densities; the four-dimensional contraction of the epsilon tensors; the kinematics;*)
(*  form/square.frm: one product M_I M_J^* of two diagrams: transversality of the gluon polarizations, colour with the package color.h, traces, polarization sums, kinematics;*)
(*  form/checks.frm: the gamma5 normalization and the short form of q-slash gamma5 used in common.h.*)


(* ::Section:: *)
(*Load*)


(* ::Text:: *)
(*The package FeynFacet: the checks Require, the reading of the amplitudes of step 1, the location of the FORM program, and the saving of the results.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*Choices: the class; the number of FORM jobs that run at the same time ("auto" = one per processor core, at most 8; "1" = one job at a time); and whether the product files of an earlier run of this same, unchanged file may be reused ("no" recomputes everything; "yes" only continues an interrupted run: a product file is reused only if its companion file <name>.input holds the fingerprint of the current inputs, the amplitudes of step 1 and this step file and of the FORM programs in form/; a file written for other inputs, or without that companion (for example, after an interrupted write), is computed again).  FeynCalc is not needed: the records of step 1 are plain expressions (FeynFacet/Amplitudes.wl explains their vocabulary).  The flavour number Nf of the closed quark loops is written by step 1 with FeynCalc's symbol; it is renamed to the plain Nf of the other routes.*)


(* ::Input:: *)
class = "13";
workers = "auto";
resume = "no";
loaded = LoadAmplitudes["class-" <> class] /. Symbol["FeynCalc`Nf"] -> Nf;
records = loaded["Records"];
spec = loaded["Spec"]; weights = loaded["Weights"];
spec


(* ::Section:: *)
(*Conventions*)


(* ::Text:: *)
(*As in the FeynCalc route: the two incoming partons; the sign of the density of a polarized incoming quark (-1 for an antiquark); the averages over the initial spins and colours (1/(2 Nc) for a quark; 1/((D-2)(Nc^2-1)) for a gluon in the spin average, 1/(Nc^2-1) for the circular density, which carries its 1/2 itself); and the reference vectors r of the gluon polarization sums: Born pa -> pb, pb -> pa, p -> pa, k1 -> pa; real emission in addition k2 -> pb.  Colour factors are written with CF = (Nc^2-1)/(2 Nc) and CA = Nc.  One unobserved momentum is not independent: k1 = pa + pb - p (Born and one loop) and k2 = pa + pb - p - k1 (real emission) inside propagators and vertices; the external spinors keep their own momentum, so the spin sums give k1-slash or k2-slash, whose scalar products are part of the kinematics.*)


(* ::Input:: *)
incoming = (spec["Born"] /. None -> spec["Real"][[1, 1]])[[1 ;; 2]];       (* the two incoming partons *)
polarizedQuarkSign = <|pa -> If[StringEndsQ[incoming[[1]], "b"], -1, 1], pb -> If[StringEndsQ[incoming[[2]], "b"], -1, 1]|>;
incomingGluons = Select[{pa, pb}, incoming[[Position[{pa, pb}, #][[1, 1]]]] === "g" &];
incomingQuarks = Complement[{pa, pb}, incomingGluons];
average[spinCase_] := Times @@ (Switch[incoming[[#]], "g", If[spinCase === "LL", 1, 1/(D - 2)]/(Nc^2 - 1), _, 1/(2 Nc)] & /@ {1, 2});
colourFactors[expr_] := Factor[expr] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
bornReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa|>;
realReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa, k2 -> pb|>;
bornConservation = k1 -> pa + pb - p;
realConservation = k2 -> pa + pb - p - k1;
{incomingQuarks, incomingGluons, average["UU"]}


(* ::Text:: *)
(*The D-dimensional scalar products of the two kinematics (the FORM procedures bornKinematics and realKinematics of common.h contain the same table): s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2; in the real emission A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 and s4 = (k1+k2)^2 = s + t + u = -A - B - C.  Mathematica needs them for the values of the propagator denominators and of the factors 1/(l.r) of the polarization sums, which FORM keeps as symbols.  The check: the eliminated momenta are massless, k2^2 = 0 with C = -(s+t+u) - A - B, and k1^2 = 0 at the Born point u = -s - t.*)


(* ::Input:: *)
bornProducts = <|{pa, pa} -> 0, {pb, pb} -> 0, {p, p} -> 0, {k1, k1} -> 0, {pa, pb} -> s/2, {pa, p} -> -t/2, {pb, p} -> -u/2,
  {pa, k1} -> -u/2, {pb, k1} -> -t/2, {p, k1} -> s/2|>;
realProducts = <|{pa, pa} -> 0, {pb, pb} -> 0, {p, p} -> 0, {k1, k1} -> 0, {k2, k2} -> 0, {pa, pb} -> s/2, {pa, p} -> -t/2, {pb, p} -> -u/2,
  {pa, k1} -> -A/2, {pb, k1} -> -B/2, {p, k1} -> C/2, {pa, k2} -> (s + t + A)/2, {pb, k2} -> (s + u + B)/2, {p, k2} -> (-t - u - C)/2, {k1, k2} -> (s + t + u)/2|>;
(* the scalar product of two momenta (linear combinations of pa, pb, p, k1, k2) in a table of products *)
dot[q1_, q2_, table_] := Module[{v = {pa, pb, p, k1, k2}},
  Expand[Sum[Coefficient[q1, a] Coefficient[q2, b] Lookup[table, Key[Sort[{a, b}]], Lookup[table, Key[{a, b}], Lookup[table, Key[{b, a}], 0]]], {a, v}, {b, v}]]];
Require[Expand[dot[pa + pb - p - k1, pa + pb - p - k1, realProducts] /. C -> -(s + t + u) - A - B] === 0 &&
  Expand[dot[pa + pb - p, pa + pb - p, bornProducts] /. u -> -s - t] === 0, "The kinematic tables are not consistent"];
"kinematics defined"


(* ::Section:: *)
(*The squaring method*)


(* ::Text:: *)
(*The method is that of the FeynCalc route (read its five steps); here it is carried out by FORM:*)
(**)
(*1. Every diagram is written once as a FORM expression, before any product, with momentum conservation applied; the conditions eps(l).l = 0 and eps(l).r = 0 for every gluon l, where r is the reference vector of the same gluon's later polarization sum, are imposed on it at the start of every job (FORM procedure transversality).*)
(*2. Colour is separated from the Lorentz algebra: each diagram is written as sum_x c_x L_x with distinct colour structures c_x, and in each job the colour factors of the pairs c_x c_y^* are evaluated by color.h before the Dirac algebra.*)
(*3. Only the products with i <= j are computed: for i < j, M_i M_j^* + M_j M_i^* = 2 Re(M_i M_j^*).*)
(*4. After the spin sums and the traces, the gluon polarization sums (or circular densities) are multiplied in and contracted by FORM.*)
(*5. The products (i, j) are independent FORM jobs, run in parallel, each with its own result file (Results/FORM/class-<n>/pairs/).*)
(**)
(*This section writes the amplitudes as FORM expressions (first three cells), prepares the process directories, and defines the job runner and the sum over the pairs.*)


(* ::Text:: *)
(*A record is one term of a diagram (FeynFacet/Amplitudes.wl): a coefficient (couplings, colour objects, propagators 1/q^2 as Den[q], metric tensors, momentum components, polarization vectors) and open fermion chains (and, for a graph with a closed quark loop, the loop's Dirac trace, listed under "Traces"; legacy records also list there the scalar factors of the loop's vertices, as separate entries, which join the coefficient).  Each record is turned into*)
(*  - its colour part: the product of the colour objects, in the notation of color.h (T(i,j,a,b,...) = (T^a T^b ...)_{ij}, f(a,b,c), d_);*)
(*  - its Lorentz part: everything else, with each polarization vector of a gluon l removed and its index renamed to the open index mu<l> (nu<l> in the conjugate), each chain ubar(q1) ... u(q2) written as Str(q1, ..., q2) (a gamma matrix as its index, a slashed momentum as a vector; a momentum sum gets an auxiliary vector Q<n>, expanded by FORM after the traces), the quark loop as the trace g_(3, ...), and each propagator 1/q^2 as a symbol den<n> (dj<n> for a loop propagator) whose value is put in later.*)
(*Momentum conservation (and, for a one-loop graph, the routing of the loop momentum chosen in the Virtual section) is applied inside the propagators, slashes and vertices, not to the spinors.*)
(*The complex conjugate of a record (for the right-hand side M_J^* of a product): i -> -i; the generators transposed, (T^a T^b)_{ij}^* = (T^b T^a)_{ji}; the chains reversed, (ubar(q1) X1 ... Xn u(q2))^* = ubar(q2) Xn ... X1 u(q1) (the QCD vertices contain no gamma5); and every Lorentz index and every summed colour index renamed with a prefix c, the open gluon index mu<l> to nu<l>, so that the indices of the two sides stay apart.  The colour indices of the external partons (those that appear once in a term) are shared by the two sides and summed in the product.*)


(* ::Input:: *)
termsOf[x_] := If[Head[x] === Plus, List @@ x, {x}];
factorsOf[x_] := If[Head[x] === Times, List @@ x, {x}];
colourQ[x_] := MatchQ[x, _ColorT | _ColorF | _ColorDelta | _GluonDelta | Power[_ColorT | _ColorF | _ColorDelta | _GluonDelta, _]];
gammaQ[x_] := MatchQ[x, _GammaMatrix | _Slash | GammaFive];
(* a record as {coefficient, chains, loop trace items}: the scalar factors listed among the traces join the coefficient.
   The four-gluon vertex brings summed adjoint indices named by FeynCalc ($AL$n); they are renamed Glu901, Glu902, ... *)
recordParts[record_] := Module[{r, loops, gammaLoops, coefficient, dummies},
  dummies = Union[Cases[record, x_Symbol /; StringStartsQ[SymbolName[x], "$AL"], Infinity]];
  r = record /. Table[dummies[[k]] -> Symbol["Glu" <> ToString[900 + k]], {k, Length[dummies]}];
  loops = Lookup[r, "Traces", {}]; coefficient = r["Coefficient"];
  gammaLoops = Select[loops, !FreeQ[#, GammaMatrix | Slash] &];
  Require[Length[gammaLoops] <= 1, "A record with more than one quark loop", r];
  coefficient = coefficient Times @@ Flatten[Select[List @@ #, !gammaQ[#] &] & /@ loops];
  {coefficient, r["Chains"], If[gammaLoops === {}, None, Select[List @@ First[gammaLoops], gammaQ]]}];
(* momentum conservation (rule) and, for a one-loop graph, the routing of the loop momentum (shift), inside propagators, slashes and vertices (not in the spinors) *)
applyMomenta[{coefficient_, chains_, loop_}, rule_, shift_] := Module[{f = Expand[(# /. shift) /. rule] &},
  {coefficient /. {Den[q_] :> Den[f[q]], MomentumComponent[q_, mu_] :> MomentumComponent[f[q], mu]},
   chains /. Slash[q_] :> Slash[f[q]], If[loop === None, None, loop /. Slash[q_] :> Slash[f[q]]]}];
(* the polarization vectors are removed: the index of the gluon l becomes the open index mu<l> *)
stripPolarizations[{coefficient_, chains_, loop_}] := Module[{pols = Cases[factorsOf[coefficient], PolVector[l_, mu_, _] :> {l, mu}], renaming},
  renaming = (#[[2]] -> Symbol["mu" <> ToString[#[[1]]]]) & /@ pols;
  {(coefficient /. PolVector[__] -> 1) /. renaming, chains /. renaming, loop /. renaming}];
colourIndicesOf[c_] := Cases[c, x_Symbol /; StringMatchQ[SymbolName[x], ("Col" | "Glu") ~~ DigitCharacter ..], Infinity];
externalColour[c_] := Select[Union[colourIndicesOf[c]], Count[colourIndicesOf[c], #] == 1 &];
(* the complex conjugate, with the renamed indices of the right-hand side *)
conjugateRecord[{coefficient_, chains_, loop_}] := Module[{colour = Times @@ Select[factorsOf[coefficient], colourQ], external, renaming},
  external = externalColour[colour];
  renaming = Join[(# -> Symbol["c" <> SymbolName[#]]) & /@ Complement[Union[colourIndicesOf[colour]], external],
    (# -> Symbol["c" <> SymbolName[#]]) & /@ Union[Cases[{coefficient, chains}, (x_Symbol /; StringStartsQ[SymbolName[x], "Lor"]), Infinity]],
    (# -> Symbol["nu" <> StringDrop[SymbolName[#], 2]]) & /@ Union[Cases[{coefficient, chains}, (x_Symbol /; StringStartsQ[SymbolName[x], "mu"]), Infinity]]];
  Require[loop === None, "A quark loop on the conjugate side"];
  {(coefficient /. Complex[a_, b_] :> Complex[a, -b] /. ColorT[l_List, i_, j_] :> ColorT[Reverse[l], j, i]) /. renaming,
   (Reverse /@ chains) /. renaming, None}];
"record functions defined"


(* ::Text:: *)
(*The FORM text of the colour and Lorentz parts.  A small printer turns them into FORM syntax.  While one process is written, three registries collect the auxiliary vectors Q<n> (with the momentum sums they stand for), the denominator symbols den<n> (with the momentum of the propagator; q and -q share one symbol) and every index, for the declarations of process.h.  A diagram file amp<i>.inc (bar<j>.inc for a conjugate) defines Amp = sum_x c_x LorentzPartA(x), its number NAMP of colour structures, and the Lorentz parts LorA<x> (LorentzPartB, NBAR and LorB<y> for the conjugate), as square.frm expects.*)


(* ::Input:: *)
newRegistry[] := (auxiliaryVectors = <||>; denominatorSymbols = <||>; usedIndices = {});
canonicalMomentum[q_] := With[{e = Expand[q]}, First[Sort[{e, Expand[-e]}]]];
denominatorSymbol[q_] := With[{key = canonicalMomentum[q]},
  If[!KeyExistsQ[denominatorSymbols, key], denominatorSymbols[key] = Symbol["den" <> ToString[Length[denominatorSymbols] + 1]]];
  denominatorSymbols[key]];
(* a slashed momentum: {numerical factor, vector}; a sum of momenta becomes an auxiliary vector *)
slashVector[q_] := Which[MatchQ[q, _Symbol], {1, q}, MatchQ[q, c_?NumericQ v_Symbol], {q[[1]], q[[2]]},
  True, If[!KeyExistsQ[auxiliaryVectors, q], auxiliaryVectors[q] = Symbol["Q" <> ToString[Length[auxiliaryVectors] + 1]]]; {1, auxiliaryVectors[q]}];
formNumber[x_] := Which[IntegerQ[x], ToString[x], Head[x] === Rational, "(" <> ToString[Numerator[x]] <> "/" <> ToString[Denominator[x]] <> ")",
  Head[x] === Complex && Re[x] == 0, "(" <> formNumber[Im[x]] <> "*i_)",
  Head[x] === Complex, "(" <> formNumber[Re[x]] <> "+" <> formNumber[Im[x]] <> "*i_)", True, Require[False, "Unknown number", x]];
(* formText: an expression as FORM text; every index met on the way is registered *)
formText[x_?NumericQ] := formNumber[x];
formText[x_Symbol] := (If[StringMatchQ[SymbolName[x], ("Lor" | "cLor" | "mu" | "nu" | "Col" | "cCol" | "Glu" | "cGlu") ~~ ___], AppendTo[usedIndices, x]]; SymbolName[x]);
formText[x_Plus] := "(" <> StringRiffle[formText /@ List @@ x, "+"] <> ")";
formText[x_Times] := StringRiffle[formText /@ List @@ x, "*"];
formText[Power[x_, n_Integer]] := "(" <> formText[x] <> ")^" <> ToString[n];
formText[Den[q_]] := formText[denominatorSymbol[q]];
formText[Metric[a_, b_]] := "d_(" <> formText[a] <> "," <> formText[b] <> ")";
formText[MomentumComponent[q_, mu_]] := formText[Sum[Coefficient[q, v] FORMComponent[v, mu], {v, Variables[q]}]];
formText[FORMComponent[v_, mu_]] := SymbolName[v] <> "(" <> formText[mu] <> ")";
formText[ColorT[l_List, i_, j_]] := "T(" <> StringRiffle[formText /@ Join[{i, j}, l], ","] <> ")";
formText[ColorF[a_, b_, c_]] := "f(" <> StringRiffle[formText /@ {a, b, c}, ","] <> ")";
formText[ColorDelta[i_, j_] | GluonDelta[i_, j_]] := "d_(" <> formText[i] <> "," <> formText[j] <> ")";
formText[x_] := Require[False, "No FORM form for", x];
(* a linear combination of momenta, c1 q1 + c2 q2 + ..., as FORM text *)
linearForm[q_] := StringJoin[Table[With[{c = Coefficient[q, v]}, Which[c === 0, "", c === 1, "+" <> SymbolName[v], c === -1, "-" <> SymbolName[v],
    True, If[c > 0, "+", "-"] <> formNumber[Abs[c]] <> "*" <> SymbolName[v]]], {v, Variables[q]}]];
(* the gamma matrices of a chain or loop: {product of the numerical factors of the slashes, list of FORM arguments} *)
gammaArguments[items_List] := Module[{parts = Map[Switch[Head[#], GammaMatrix, {1, formText[#[[1]]]},
    Slash, With[{vector = slashVector[#[[1]]]}, {vector[[1]], SymbolName[vector[[2]]]}]] &, items]},
  {Times @@ parts[[All, 1]], parts[[All, 2]]}];
(* a chain as Str(q1, items, q2), the quark loop as g_(3, items) *)
chainForm[Chain[left_, items___, right_]] := With[{g = gammaArguments[{items}]},
  {g[[1]], "Str(" <> StringRiffle[Join[{SymbolName[left[[1]]]}, g[[2]], {SymbolName[right[[1]]]}], ","] <> ")"}];
loopForm[items_] := With[{g = gammaArguments[items]}, {g[[1]], "g_(3," <> StringRiffle[g[[2]], ","] <> ")"}];
(* one prepared record -> {colour part, FORM text of the Lorentz part} *)
recordForm[{coefficient_, chains_, loop_}] := Module[{colour = Times @@ Select[factorsOf[coefficient], colourQ], rest, chainParts, loopPart, factor},
  rest = coefficient/colour;
  chainParts = chainForm /@ chains; loopPart = If[loop === None, {1, None}, loopForm[loop]];
  factor = rest (Times @@ chainParts[[All, 1]]) loopPart[[1]];
  {colour, StringRiffle[Join[{formText[factor]}, chainParts[[All, 2]], If[loop === None, {}, {loopPart[[2]]}]], "*"]}];
(* a diagram (list of prepared records) -> the text of amp<i>.inc (side "A") or bar<j>.inc (side "B") *)
diagramText[prepared_List, side_, comment_] := Module[{forms = recordForm /@ prepared, groups, name = If[side === "A", "Amp", "Bar"]},
  groups = GatherBy[forms, First];
  StringJoin["* ", comment, "\n",
    "L ", name, " = ", StringRiffle[Table[formText[groups[[x, 1, 1]]] <> "*LorentzPart" <> side <> "(" <> ToString[x] <> ")", {x, Length[groups]}], "\n  + "], ";\n",
    "#define N", If[side === "A", "AMP", "BAR"], " \"", ToString[Length[groups]], "\"\n",
    Table["L Lor" <> side <> ToString[x] <> " =\n    " <> StringRiffle[groups[[x, All, 2]], "\n  + "] <> ";\n", {x, Length[groups]}]]];
"FORM printer defined"


(* ::Text:: *)
(*process.h: the declarations of one process and what the FORM procedures need to know about it: the kinematics (Born, Virtual or Real), the momenta with (D-4)-dimensional parts, the gluons with their reference vectors, the incoming quarks with the signs of their densities, the indices (Lorentz indices D-dimensional, colour indices in the fundamental (NR) or adjoint (NA) representation), the auxiliary vectors and their definitions, and the transversality conditions eps(l).l = eps(l).r = 0 on both sides of the product.  For the gluon x whose momentum was eliminated by momentum conservation, eps(x).x = 0 is imposed in its eliminated form, solved for the last momentum that is not x's reference vector (as in the FeynCalc route).*)


(* ::Input:: *)
transversalityLines[references_, rule_, gluons_] := Module[{x = rule[[1]], momenta, solved, coefficient, lines = {}},
  Do[With[{index = Symbol[side <> ToString[l]]},
      If[l === x,
        momenta = Cases[termsOf[Expand[rule[[2]]]], c_. q_Symbol :> q];
        solved = Last[DeleteCases[momenta, references[l]]];
        coefficient[q_] := Coefficient[rule[[2]], q];
        AppendTo[lines, "  id " <> SymbolName[solved] <> "(" <> SymbolName[index] <> ") = " <>
          StringReplace[linearForm[Expand[-Sum[coefficient[q] q, {q, DeleteCases[momenta, solved]}]/coefficient[solved]]], v : (LetterCharacter ~~ WordCharacter ...) :> v <> "(" <> SymbolName[index] <> ")"] <> ";"]];
      AppendTo[lines, "  id " <> SymbolName[l] <> "(" <> SymbolName[index] <> ") = 0;"];
      AppendTo[lines, "  id " <> SymbolName[references[l]] <> "(" <> SymbolName[index] <> ") = 0;"]], {l, gluons}, {side, {"mu", "nu"}}];
  lines];
processHeader[name_, kinematics_, references_, rule_, gluons_, hatted_] := Module[{lorentz, colour},
  usedIndices = Union[usedIndices, Flatten[{Symbol["mu" <> ToString[#]], Symbol["nu" <> ToString[#]]} & /@ gluons]];
  lorentz = Select[usedIndices, StringMatchQ[SymbolName[#], ("Lor" | "cLor" | "mu" | "nu") ~~ ___] &];
  colour = Complement[usedIndices, lorentz];
  StringJoin["* process.h: ", name, " (written by 2-Algebra-FORM.wl)\n",
    "#define KINEMATICS \"", kinematics, "\"\n",
    "#define HATTED \"", StringRiffle[SymbolName /@ hatted, ","], "\"\n",
    "* gluons: incoming, outgoing, and the reference vectors of their polarization sums\n",
    "#define INGLUONS \"", StringRiffle[SymbolName /@ Intersection[gluons, incomingGluons], ","], "\"\n",
    "#define OUTGLUONS \"", StringRiffle[SymbolName /@ Complement[gluons, {pa, pb}], ","], "\"\n",
    Table["#define REF" <> SymbolName[l] <> " \"" <> SymbolName[references[l]] <> "\"\n", {l, gluons}],
    "* incoming quarks and the signs of their densities (helicity difference)\n",
    If[MemberQ[incomingQuarks, pa], "#define QUARKa \"pa\"\n#define SIGNa \"" <> ToString[polarizedQuarkSign[pa]] <> "\"\n", ""],
    If[MemberQ[incomingQuarks, pb], "#define QUARKb \"pb\"\n#define SIGNb \"" <> ToString[polarizedQuarkSign[pb]] <> "\"\n", ""],
    If[lorentz === {}, "", "Indices " <> StringRiffle[SymbolName /@ lorentz, ","] <> ";\n"],
    If[colour === {}, "", "Indices " <> StringRiffle[SymbolName[#] <> If[StringContainsQ[SymbolName[#], "Col"], "=NR", "=NA"] & /@ colour, ","] <> ";\n"],
    If[auxiliaryVectors === <||>, "", "Vectors " <> StringRiffle[SymbolName /@ Values[auxiliaryVectors], ","] <> ";\n"],
    "Symbols ", StringRiffle[Join[SymbolName /@ Values[denominatorSymbols], {"dj1", "dj2", "dj3", "dj4"}, ("invlr" <> SymbolName[#]) & /@ gluons], ","], ";\n",
    "* the auxiliary vectors: momentum sums inside the gamma strings\n#procedure auxiliary\n",
    KeyValueMap["  id " <> SymbolName[#2] <> " = " <> linearForm[#1] <> ";\n" &, auxiliaryVectors],
    "#endprocedure\n",
    "* transversality of the gluon polarizations: eps(l).l = 0 and eps(l).r = 0 (r the reference vector)\n#procedure transversality\n",
    StringRiffle[transversalityLines[references, rule, gluons], "\n"], "\n#endprocedure\n"]];
"process header defined"


(* ::Text:: *)
(*Each process (the Born process, the one-loop interference, every real-emission process, and the process of the gauge check) gets a directory Results/FORM/class-<n>/work/<process>/ with process.h, the diagrams amp<i>.inc and the conjugate diagrams bar<j>.inc.  The values of the denominator symbols den<n> = 1/q^2 and of the factors invlr<l> = 1/(l.r) of the polarization sums are kept (symbolValues) for the reading of the results.*)


(* ::Input:: *)
classDirectory = ResultsDirectory["FORM/class-" <> class];
workDirectory[name_] := With[{dir = FileNameJoin[{classDirectory, "work", name}]},
  If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]]; dir];
momentaOf[particles_] := {pa, pb, p, k1, k2}[[1 ;; Length[particles]]];
gluonsOf[particles_] := Pick[momentaOf[particles], particles, "g"];
prepareDiagram[recs_, rule_, shift_] := stripPolarizations[applyMomenta[recordParts[#], rule, shift]] & /@ recs;
symbolValues = <||>;
(* write one process: amp<i>.inc for the left-hand diagrams, bar<j>.inc for the conjugates of the right-hand ones *)
writeProcess[name_, kinematics_, references_, rule_, particles_, left_List, right_List] := Module[{dir = workDirectory[name], gluons = gluonsOf[particles], table},
  newRegistry[];
  Do[Export[FileNameJoin[{dir, "amp" <> ToString[i] <> ".inc"}], diagramText[left[[i]], "A", name <> ": diagram " <> ToString[i]], "Text"], {i, Length[left]}];
  Do[Export[FileNameJoin[{dir, "bar" <> ToString[j] <> ".inc"}], diagramText[conjugateRecord /@ right[[j]], "B", name <> ": conjugate of diagram " <> ToString[j]], "Text"], {j, Length[right]}];
  Export[FileNameJoin[{dir, "process.h"}], processHeader[name, kinematics, references, rule, gluons, If[kinematics === "Real", {k1, k2}, If[kinematics === "Virtual", {loop}, {}]]], "Text"];
  table = If[kinematics === "Real", realProducts, bornProducts];
  symbolValues[name] = Join[KeyValueMap[#2 -> 1/dot[#1, #1, table] &, KeySelect[denominatorSymbols, FreeQ[#, loop] &]],
    Table[Symbol["invlr" <> SymbolName[l]] -> 1/dot[l, references[l], table], {l, gluons}]];
  dir];
"process writer defined"


(* ::Text:: *)
(*The jobs.  A job is one product: diagram I of a process times the conjugate of diagram J, for one spin case.  At most `workers` FORM processes run at the same time; each writes its result (the list {real part, imaginary part}) under its name in the working directory of its process, from where it is moved to Results/FORM/class-<n>/pairs/<name>.m, next to the fingerprint of its inputs, <name>.m.input; its output to <name>.log (kept only if the job fails, and shown in the error message), and its running time in seconds to <name>.time.  The FORM program square.frm finds common.h in form/, color.h in the FORM colour package, and the process files in its working directory.*)


(* ::Input:: *)
nWorkers = If[workers === "auto", Min[$ProcessorCount, 8], ToExpression[workers]];   (* "auto": one per processor core, at most 8 *)
formDirectory = FileNameJoin[{NotebookDirectory[], "form"}];
pairDirectory = FileNameJoin[{classDirectory, "pairs"}];
If[!DirectoryQ[pairDirectory], CreateDirectory[pairDirectory]];
If[resume === "no", DeleteFile /@ FileNames["*", pairDirectory]];
resultFile[job_] := FileNameJoin[{pairDirectory, job["Name"] <> ".m"}];
logFile[job_] := FileNameJoin[{pairDirectory, job["Name"] <> ".log"}];
(* the shell command of a job, run in the working directory of its process.  Every argument is quoted for bash (ShellQuote), so that a path with
   spaces or quotes stays one argument.  FORM writes its result under the plain job name in that working directory (a full path as the value of the
   macro OUT could contain characters that FORM's preprocessor does not accept); the finished file is then moved to the pairs directory *)
command[job_] := StringRiffle[Join[ShellQuote /@ {ToolPath["FORM"],
    "-d", "SPIN=" <> If[job["Spin"] === "UU", "U", "P"], "-d", "I=" <> ToString[job["I"]], "-d", "J=" <> ToString[job["J"]], "-d", "OUT=" <> job["Name"],
    "-p", StringRiffle[{workDirectory[job["Process"]], formDirectory, DirectoryName[ToolPath["FORMColor"]]}, ":"],
    FileNameJoin[{formDirectory, "square.frm"}]}, {">", ShellQuote[logFile[job]], "2>&1"}], " "];
localResult[job_] := FileNameJoin[{workDirectory[job["Process"]], job["Name"] <> ".m"}];
(* the fingerprint of the inputs of the product files: the amplitudes of step 1, this step file and the FORM programs (see the choice resume above) *)
inputFingerprint = InputFingerprint[loaded, Join[{FileNameJoin[{NotebookDirectory[], "2-Algebra-FORM.wl"}]}, Select[FileNames["*", formDirectory], !DirectoryQ[#] &]]];
(* a queue of jobs: start new ones while fewer than nWorkers run, collect the finished ones; stop everything at the first failure *)
runJobs[jobs_] := Module[{queue, running = {}, finished, code, job, process},
  queue = Select[jobs, !(resume === "yes" && ReusableQ[resultFile[#], inputFingerprint]) &];
  While[queue =!= {} || running =!= {},
    While[Length[running] < nWorkers && queue =!= {},
      job = First[queue]; queue = Rest[queue];
      If[FileExistsQ[localResult[job]], DeleteFile[localResult[job]]];
      (* exec makes the tracked process the FORM worker, so failure cleanup
         also stops the worker rather than leaving it behind a shell. *)
      process = StartProcess[{"bash", "-c", "exec " <> command[job]}, ProcessDirectory -> workDirectory[job["Process"]]];
      If[!MatchQ[process, _ProcessObject], KillProcess /@ running[[All, 2]];
        Require[False, "Could not start a FORM job; check the shell and operating-system process limits", job["Name"]]];
      AppendTo[running, {job, process, AbsoluteTime[]}]];
    Pause[0.05];
    finished = Select[running, ProcessStatus[#[[2]]] === "Finished" &];
    Do[code = ProcessInformation[f[[2]], "ExitCode"];
      running = DeleteCases[running, f];
      (* A finished ProcessObject still owns native streams. KillProcess
         releases those resources as well as closing the Wolfram streams. *)
      KillProcess[f[[2]]];
      If[code =!= 0, KillProcess /@ running[[All, 2]];
        Require[False, "A FORM job failed; see its log", {f[[1]]["Name"], Import[logFile[f[[1]]], "Text"]}]];
      If[!FileExistsQ[localResult[f[[1]]]], KillProcess /@ running[[All, 2]]];
      Require[FileExistsQ[localResult[f[[1]]]], "A FORM job wrote no result; see its log", {f[[1]]["Name"], Import[logFile[f[[1]]], "Text"]}];
      If[FileExistsQ[resultFile[f[[1]]]], DeleteFile[resultFile[f[[1]]]]];
      RenameFile[localResult[f[[1]]], resultFile[f[[1]]]];
      Put[AbsoluteTime[] - f[[3]], FileNameJoin[{pairDirectory, f[[1]]["Name"] <> ".time"}]];
      MarkReusable[resultFile[f[[1]]], inputFingerprint];
      DeleteFile[logFile[f[[1]]]], {f, finished}]];
  (* the values {real part, imaginary part} of all jobs, in their order *)
  Table[Get[resultFile[j]], {j, jobs}]];
"job runner defined"


(* ::Text:: *)
(*Before the squares: checks.frm verifies the two forms of the HVBM gamma5 used in common.h, Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = -4 i e_(a,b,c,d) (the normalization, the same as in the FeynCalc and Tracer routes) and q-slash gamma5 = -i/6 e_(q,a,b,c) gamma^a gamma^b gamma^c inside a trace with D-dimensional momenta; FORM must print Gamma5Checks = 0.*)


(* ::Input:: *)
checkRun = RunProcess[{ToolPath["FORM"], "-p", StringRiffle[{formDirectory, DirectoryName[ToolPath["FORMColor"]]}, ":"], FileNameJoin[{formDirectory, "checks.frm"}]}];
Require[checkRun["ExitCode"] === 0 && StringContainsQ[checkRun["StandardOutput"], "Gamma5Checks = 0;"], "The gamma5 checks of checks.frm failed", checkRun["StandardOutput"]];
"gamma5 checks passed"


(* ::Text:: *)
(*The square of a process: the jobs of all pairs i <= j, and their sum (diagonal products as they are, i < j twice the real part), with the values of den<n> and invlr<l> put in.  A diagonal product |M_i|^2 is real; its imaginary part must vanish exactly at a rational point (as in the FeynCalc route).*)


(* ::Input:: *)
checkPoint = {s -> 1, t -> -3/10, u -> -2/5, A -> -1/7, B -> -2/9, C -> -3/10 + 1/7 + 2/9, D -> 37/10, Nc -> 3, spe3 -> 1/11, gs -> 1, Nf -> 5};
squareJobs[name_, n_, spinCase_] := Flatten[Table[<|"Name" -> name <> "-" <> spinCase <> "-" <> ToString[i] <> "-" <> ToString[j],
  "Process" -> name, "I" -> i, "J" -> j, "Spin" -> spinCase|>, {i, n}, {j, i, n}]];
square[name_, n_, spinCase_] := Module[{jobs = squareJobs[name, n, spinCase], results, values = symbolValues[name]},
  results = runJobs[jobs];
  Require[And @@ MapThread[#2["I"] != #2["J"] || (#1[[2]] /. values /. checkPoint) === 0 &, {results, jobs}], "A diagonal product has a nonzero imaginary part", name];
  Total[MapThread[If[#2["I"] == #2["J"], 1, 2] #1[[1]] &, {results, jobs}]] /. values];
"parallel FORM jobs: " <> ToString[nWorkers]


(* ::Section:: *)
(*Born*)


(* ::Text:: *)
(*The Born square with the reference vectors pa -> pb, pb -> pa (incoming gluons), p -> pa, k1 -> pa, averaged over the initial spins and colours, with u = -s - t.*)


(* ::Input:: *)
bornDiagrams = prepareDiagram[#, bornConservation, {}] & /@ records["Born"];
If[bornDiagrams =!= {}, writeProcess["born", "Born", bornReferences, bornConservation, spec["Born"], bornDiagrams, bornDiagrams]];
born = If[bornDiagrams === {}, <|"UU" -> 0, "LL" -> 0|>,
  Association @ Table[spinCase -> Factor[colourFactors[Together[average[spinCase] square["born", Length[bornDiagrams], spinCase]]] /. u -> -s - t],
    {spinCase, {"UU", "LL"}}]];
Require[FreeQ[born, Complex], "The Born square is not real", born];
born


(* ::Text:: *)
(*The Born squares must agree with the table of Processes/PartonicBorn (which was checked against the textbook).*)


(* ::Input:: *)
(* the table is read only for a class with a Born term (the classes without one do not need Processes/PartonicBorn here) *)
If[spec["Born"] =!= None,
  bornTable = Get[FileNameJoin[{NotebookDirectory[], "..", "PartonicBorn", "Results", "Born.wl"}]];
  Require[And @@ Table[Simplify[born[spinCase]/gs^4 - bornTable["Born"][spinCase][StringRiffle[spec["Born"][[1 ;; 2]], " "] <> " -> " <> spec["Born"][[3]]]] === 0, {spinCase, {"UU", "LL"}}],
    "The Born square differs from Processes/PartonicBorn", born]];
"Born checked"


(* ::Section:: *)
(*Virtual: the one-loop interference in terms of scalar integrals*)


(* ::Text:: *)
(*The families are those of the FeynCalc route (read its explanation):*)
(*   VirtualST[n1,n2,n3,n4]: 1/[(loop^2)^n1 ((loop+pa)^2)^n2 ((loop+pa+pb)^2)^n3 ((loop+p)^2)^n4]   (external legs pa, pb, -pd, -p in cyclic order: invariants s and t),*)
(*   VirtualSU: ... ((loop+pa+pb)^2)^n3 ((loop+pa+pb-p)^2)^n4   (legs pa, pb, -p, -pd: s and u),      VirtualTU: ((loop+pa)^2)^n2 ((loop+pa-p)^2)^n3 ((loop+pa-p+pb)^2)^n4   (legs pa, -p, pb, -pd: t and u),*)
(*with negative indices for numerators.  For every graph the routing loop -> sigma loop + delta that maps its loop propagators onto one family is found by trying the shifts in the same order as in the FeynCalc route, so that both routes label the integrals identically.  Each loop propagator then becomes the symbol dj<j> = 1/D_j of the family's j-th propagator; the graphs, with these symbols, are written as the left-hand diagrams of the process "virtual", and the conjugate of the summed Born amplitude as its only right-hand diagram.*)


(* ::Input:: *)
families = <|
  "VirtualST" -> {loop, loop + pa, loop + pa + pb, loop + p},                 (* legs pa, pb, -pd, -p in cyclic order: invariants s and t *)
  "VirtualSU" -> {loop, loop + pa, loop + pa + pb, loop + pa + pb - p},       (* legs pa, pb, -p, -pd: s and u *)
  "VirtualTU" -> {loop, loop + pa, loop + pa - p, loop + pa - p + pb}|>;      (* legs pa, -p, pb, -pd: t and u *)
familyHeads = <|"VirtualST" -> KVirtualST, "VirtualSU" -> KVirtualSU, "VirtualTU" -> KVirtualTU|>;
shifts = Flatten[Table[c1 pa + c2 pb + c3 p, {c1, -1, 1}, {c2, -1, 1}, {c3, -1, 1}]];
sameUpToSign[q1_, q2_] := Expand[q1 - q2] === 0 || Expand[q1 + q2] === 0;
propagatorMomenta[graph_] := Union[Cases[graph, Den[q_] /; !FreeQ[q, loop] :> Expand[q /. bornConservation], Infinity]];   (* the loop propagators; tree propagators stay as they are *)
(* find sigma, delta and the family such that every propagator of the graph, with loop -> sigma loop + delta, is one of the family's *)
routing[graph_] := Module[{props = propagatorMomenta[graph], found},
  found = Catch[Do[If[And @@ Table[Or @@ Table[sameUpToSign[q /. loop -> sigma loop + delta, f], {f, families[name]}], {q, props}],
      Throw[{name, sigma, delta}]], {name, Keys[families]}, {sigma, {1, -1}}, {delta, shifts}]; None];
  Require[found =!= None, "No family fits the propagators of a one-loop graph", props]; found];
routings = routing /@ records["Virtual"];
(* the graph with the routing applied and its loop propagators as the symbols dj<j> *)
loopGraph[recs_, {name_, sigma_, delta_}] := Module[{prepared = prepareDiagram[recs, bornConservation, loop -> sigma loop + delta]},
  prepared /. Den[q_] /; !FreeQ[q, loop] :> Symbol["dj" <> ToString[Position[families[name], f_ /; sameUpToSign[f, q], {1}, Heads -> False][[1, 1]]]]];
virtualGraphs = MapThread[loopGraph, {records["Virtual"], routings}];
Require[FreeQ[virtualGraphs, Den[q_] /; !FreeQ[q, loop]], "A loop propagator was not mapped onto its family"];
If[virtualGraphs =!= {}, writeProcess["virtual", "Virtual", bornReferences, bornConservation, spec["Born"], virtualGraphs, {Flatten[bornDiagrams, 1]}]];
Counts[routings[[All, 1]]]


(* ::Text:: *)
(*One job per graph and spin case: the graph times the conjugate of the summed Born amplitude, M1_graph M0*, with its real and imaginary parts.  FORM returns the loop momentum through the symbols L2 = loop.loop, Lpa = loop.pa, Lpb = loop.pb, Lp = loop.p and Lhat2 = loophat.loophat (the four-dimensional loop4.loop4 = L2 - Lhat2 in the HVBM scheme); the four products are written through the inverse propagators D_j of the graph's family (the four equations loop^2 = D1, 2 loop.q_j = D_j - D_1 - q_j^2 are solved for L2, Lpa, Lpb, Lp), and every term is read off as coefficient x family[indices] x HatPower[n], with n the power of Lhat2 (it shifts the dimension of the integral in step 3, not the dimension D of the coefficients).*)
(**)
(*Normalization, as in the FeynCalc route: FeynArts writes the one-loop amplitude with the measure d^D loop/(2 Pi)^D; the scalar integrals are defined with the measure d^D loop/(i Pi^(D/2)), so that M1 = i/(4 Pi)^(D/2) x (sum of coefficients times integrals).  The products (coefficient x M0*) are real, and*)
(**)
(*   2 Re(M1 M0*) = (2/(4 Pi)^(D/2)) Sum_j [i c_j M0*] Re(I_j),*)
(**)
(*with the real parts of the integrals in the physical region (s > 0) taken in step 3.  The association virtual below holds Sum_j [i c_j M0*] I_j for both spin cases (i times the complex product of FORM, which must come out real); the factor 2/(4 Pi)^(D/2) is included in step 4.*)


(* ::Input:: *)
inverseProps[name_] := Module[{q = families[name], sol},
  sol = Solve[Table[Dj[j] == L2 + 2 (Coefficient[q[[j]], pa] Lpa + Coefficient[q[[j]], pb] Lpb + Coefficient[q[[j]], p] Lp) + dot[q[[j]] - loop, q[[j]] - loop, bornProducts] /. u -> -s - t, {j, 4}], {L2, Lpa, Lpb, Lp}];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators"]; First[sol]];
loopRules = Association @ Table[name -> inverseProps[name], {name, Keys[families]}];
(* the interference of one graph -> sum of coefficient x family[indices] x HatPower[n] *)
toFamilies[value_, name_] := Module[{expr},
  expr = Expand[value /. symbolValues["virtual"] /. loopRules[name] /. Table[Symbol["dj" <> ToString[j]] -> 1/Dj[j], {j, 4}]];
  Require[FreeQ[expr, L2 | Lpa | Lpb | Lp], "A scalar product of the loop momentum was left"];
  Total[Table[With[{powers = Table[Exponent[term, Dj[j]], {j, 4}], hat = Exponent[term, Lhat2]},
    (term/(Times @@ Table[Dj[j]^powers[[j]], {j, 4}] Lhat2^hat)) familyHeads[name] @@ (-powers) HatPower[hat]], {term, DeleteCases[termsOf[expr], 0]}]]];
virtualJobs[spinCase_] := Table[<|"Name" -> "virtual-" <> spinCase <> "-" <> ToString[g], "Process" -> "virtual", "I" -> g, "J" -> 1, "Spin" -> spinCase|>, {g, Length[virtualGraphs]}];
virtual = If[virtualGraphs === {}, <|"UU" -> 0, "LL" -> 0|>,
  Association @ Table[spinCase -> Module[{results = runJobs[virtualJobs[spinCase]]},
    Collect[I average[spinCase] Total[MapThread[toFamilies[#1[[1]] + I #1[[2]], #2[[1]]] &, {results, routings}]] /. u -> -s - t,
      {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, colourFactors]], {spinCase, {"UU", "LL"}}]];
Require[FreeQ[virtual, Complex], "The virtual interference is not real", virtual];
{Length[virtual["UU"]], Union[Cases[virtual, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]] // Length}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*For every real-emission process of the class, the square of the sum of its diagrams in both spin cases, in terms of the invariants s, t, u, A, B, C and the hatted invariant spe3; the processes are added with their weights (classes.wl: 1/2 for two identical unobserved gluons, Nf - 1 or Nf for the unobserved quark pairs of other flavours).  Reference vectors of the physical polarization sums: pa for the gluons p and k1, pb for k2 (and the other incoming momentum for an incoming gluon).  As in the FeynCalc route the sum is not brought to a common denominator (step 3 expands it into partial fractions anyway); the colour factors are polynomials in Nc = CA.  In the HVBM scheme the LL square is at most linear in spe3.*)


(* ::Input:: *)
processName[key_] := "real-" <> StringReplace[key, " " -> "_"];
realDiagrams = Association @ KeyValueMap[#1 -> (prepareDiagram[#, realConservation, {}] & /@ #2) &, records["Real"]];
KeyValueMap[writeProcess[processName[#1], "Real", realReferences, realConservation, StringSplit[#1], #2, #2] &, realDiagrams];
realSquares = Association @ Table[spinCase -> Association @ KeyValueMap[#1 -> square[processName[#1], Length[#2], spinCase] &, realDiagrams], {spinCase, {"UU", "LL"}}];
real = Association @ Table[spinCase -> (average[spinCase] Total[KeyValueMap[weights[#1] #2 &, realSquares[spinCase]]] /. Nc -> CA), {spinCase, {"UU", "LL"}}];
Require[FreeQ[real, Complex], "The real-emission square is not real"];
Require[Max[Exponent[#, spe3] & /@ termsOf[Expand[real["LL"]]]] <= 1, "The real-emission square has spe3^2 terms"];
LeafCount /@ real


(* ::Text:: *)
(*Check of the gauge invariance of the real-emission square: for the first process of the class (and the spin average), the physical polarization sums with other reference vectors (pb for p and k1, pa for k2) give the same result, exactly, at three random rational points (all invariants, D and Nc, with C = -s4 - A - B).  The diagrams are written again (process "gauge"), because the transversality conditions eps(l).r = 0 depend on the reference vectors.*)


(* ::Input:: *)
alternativeReferences = <|pa -> pb, pb -> pa, p -> pb, k1 -> pb, k2 -> pa|>;
firstProcess = First[Keys[realDiagrams]];
writeProcess["gauge", "Real", alternativeReferences, realConservation, StringSplit[firstProcess], realDiagrams[firstProcess], realDiagrams[firstProcess]];
alternative = square["gauge", Length[realDiagrams[firstProcess]], "UU"];
gaugePoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], Nc -> 2 + 10 v[[6]]}], {3}];
referenceDifference = Table[(alternative - realSquares["UU"][firstProcess]) /. C -> -(s + t + u) - A - B /. point, {point, gaugePoints}];
Require[referenceDifference === {0, 0, 0}, "The real-emission square depends on the gluon reference vectors", referenceDifference];
"gauge invariance of the real emission checked"


(* ::Section:: *)
(*Timing and save*)


(* ::Text:: *)
(*The running times of the FORM jobs of this class: their number, the total and the longest, in seconds.*)


(* ::Input:: *)
pairTimes = Association @ Table[FileBaseName[f] -> Get[f], {f, FileNames["*.time", pairDirectory]}];
{Length[pairTimes], Round[Total[pairTimes]], Round[Max[pairTimes]]}


(* ::Text:: *)
(*The result file Results/FORM/class-<n>/Algebra.wl, in the format of the FeynCalc route.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Virtual" -> virtual, "Real" -> real, "Class" -> class, "Spec" -> spec,
  "Families" -> families, "Variables" -> "s, t, u = -s-t; A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 = -s4-A-B; spe3; loopHat2 via HatPower[n]"|>;
SaveResult["Algebra", algebra, "FORM/class-" <> class];
Keys[algebra]


(* ::Section:: *)
(*Comparison with the FeynCalc route*)


(* ::Text:: *)
(*If the FeynCalc route has been run for this class, the two files must agree (otherwise the comparison is skipped and the cell says so): Born, Virtual and Real of both spin cases are compared at two exact random rational points of s, t, u, A, B, D, CA, spe3, Nf, gs (with C = -(s+t+u) - A - B and CF = (CA^2 - 1)/(2 CA)); in the virtual every scalar integral is replaced by its own random rational number and HatPower[n] by hp^n, so that the coefficient of every integral is compared.  The arithmetic is exact, but the comparison is sampled at two points, not a symbolic identity.*)


(* ::Input:: *)
feynCalcFile = FileNameJoin[{ResultsDirectory["FeynCalc/class-" <> class], "Algebra.wl"}];
routeComparison = If[!FileExistsQ[feynCalcFile], "the FeynCalc route has not been run for this class",
  Module[{other = Get[feynCalcFile], integrals, points, evaluate},
    integrals = Union[Cases[{algebra["Virtual"], other["Virtual"]}, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]];
    points = Table[Join[Thread[{s, t, u, A, B, D, CA, spe3, Nf, gs, hp} -> RandomInteger[{1, 10^6}, 11]/RandomInteger[{1, 10^6}, 11] {1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1}],
      Thread[integrals -> RandomInteger[{1, 10^9}, Length[integrals]]/RandomInteger[{1, 10^9}, Length[integrals]]]], {2}];
    evaluate[x_, point_] := Together[x /. HatPower[n_] :> hp^n /. point /. C -> -(s + t + u) - A - B /. CF -> (CA^2 - 1)/(2 CA) /. point];
    Association @ Flatten[Table[{part, spinCase} -> (And @@ Table[evaluate[algebra[part][spinCase] - other[part][spinCase], point] === 0, {point, points}]),
      {part, {"Born", "Virtual", "Real"}}, {spinCase, {"UU", "LL"}}]]]];
Require[StringQ[routeComparison] || And @@ Values[routeComparison], "The FORM route differs from the FeynCalc route", routeComparison];
routeComparison
