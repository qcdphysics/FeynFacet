(* ::Title:: *)
(*Transversely polarized proton-proton classes (A_TT), step 2 (FORM route): spin and colour algebra*)


(* ::Text:: *)
(*Input: the amplitude records of step 1 for the class chosen below.  Output: Results/FORM/class-<n>/Algebra.wl, with the same content, conventions and format as the file of the FeynCalc route (2-Algebra-FeynCalc.wl, read it first: it defines the observable, the spin densities, the axes and the four components):*)
(*  Born:       the Born square for the two transverse spin projections TT = (XX - NN)/2 and Symmetric = (XX + NN)/2 (and the components XX, NN themselves), summed over final spins and colours, averaged over the initial colours, in D dimensions;*)
(*  Virtual:    the interference 2 Re(M_1loop M_Born^*) for the same projections, written as Sum_j [i c_j M0*] I_j in terms of the scalar one-loop integrals of the families VirtualST, VirtualSU, VirtualTU, with the powers of the (D-4)-dimensional part of the loop momentum as HatPower[n];*)
(*  Real:       the weighted sum over the real-emission processes of the class of |M|^2 for the same projections, in the invariants s, t, u, A, B, C and the hatted invariant spe3;*)
(*  BornTable:  the D-dimensional transversely polarized Born squares |M|^2/gs^4 of every 2 -> 2 channel needed by the collinear counterterms, in the format of Processes/PartonicBorn.*)
(*The last section checks that the file agrees exactly with the one of the FeynCalc route.*)
(**)
(*Spin densities (docs/TRANSVERSE-PP.md).  Both incoming partons are quarks or antiquarks with a transverse spin, inserted as pa-slash gamma5 sa-slash/2 and pb-slash gamma5 sb-slash/2 (the same expression for a quark and for an antiquark; the HVBM gamma5, four-dimensional spin vectors).  The spin vectors are written in the basis of the scattering plane of pa, pb, p: the in-plane axis eX = Sqrt[s/(t u)] aX with aX = (u/s) pa + (t/s) pb + p (aX.pa = aX.pb = 0, aX.aX = aX.p = -t u/s), and the unit normal n (n.pa = n.pb = n.p = 0, n.n = -1); sa = cos(phia) eX + sin(phia) n, sb likewise, and the squared matrix element is XX cos cos + NN sin sin + XN cos sin + NX sin cos.  The four components are computed separately (XX with the unnormalized axis aX on both legs, divided by -aX.aX = t u/s; XN and NX, which must vanish by parity, up to the positive factor Sqrt[s/(t u)]).  The hard part of A_TT is the coefficient of cos(2 Phi), TT = (XX - NN)/2 (the covariant projector F = (s/(t u)) [2 (p.sa)(p.sb) + (t u/s) sa.sb] of Mukherjee, Stratmann and Vogelsang); the symmetric part (XX + NN)/2 vanishes in four dimensions but not in the HVBM scheme, and is carried along as a check.*)
(**)
(*Normal components.  n.k1 = -n.k2 = nk1 appears in NN, XN and NX of the real emission; odd powers integrate to zero, and nk1^2 = k1par^2 + spe3 with k1par^2 = d.G^-1.d, G the Gram matrix of pa, pb, p and d = (pa.k1, pb.k1, p.k1) (Processes/Collins, normalReduction).  The loop momentum has a normal component nloop = n.loop too: (n.loop)^2 = looppar^2 - loop^2 + loopHat2, odd powers dropped.  Both reductions are done here, in Mathematica, on the results of FORM, exactly as in the FeynCalc route.*)
(**)
(*Division of the work.  Mathematica reads the amplitude records of step 1, applies momentum conservation, splits every diagram into its colour structures and Lorentz parts and writes it as a FORM expression; it starts the FORM jobs, reads their results, puts in the kinematic values of the propagator denominators and reduces the normal components.  FORM does all the spin and colour algebra, including the transverse densities.  FeynCalc and FeynArts are loaded only to generate the diagrams of the 2 -> 2 channels of the Born table (which step 1 does not provide) and to turn them into records; their squares are again computed by FORM.  The FORM programs are in the directory form/ next to this file and are meant to be read:*)
(*  form/common.h:   declarations; the spin sums and the transverse densities of the two incoming partons (with the HVBM gamma5 and the symbolic spin vectors sa, sb); the gluon polarization sums; the four-dimensional contraction of the two epsilon tensors; the substitution of the axes aX, n for the spin vectors; the kinematics;*)
(*  form/square.frm: one product M_I M_J^* of two diagrams for one component: transversality of the gluon polarizations, colour with the package color.h, traces, polarization sums, epsilon contraction, axes, kinematics;*)
(*  form/checks.frm: the gamma5 normalization, the short form of q-slash gamma5 and the short form of the transverse density used in common.h.*)
(**)
(*Invariants as in Processes/Hadronic: s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2; real emission: s4 = (k1+k2)^2 = s + t + u, A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 with A + B + C = -s4; spe3 = k1hat.k1hat.*)


(* ::Section:: *)
(*Load*)


(* ::Text:: *)
(*The package FeynFacet: the checks Require, the reading of the amplitudes of step 1, the location of the FORM program, and the saving of the results.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*FeynCalc with FeynArts is loaded first (as in the FeynCalc route), so that every symbol of this file (Nf, CA, CF) is the same symbol as in the records of step 1 and in the result file of the FeynCalc route; it is used only to generate the diagrams of the Born table.  Choices: the class; the number of FORM jobs that run at the same time ("auto" = one per processor core, at most 8; "1" = one job at a time); and whether the product files of an earlier run of this same, unchanged file may be reused ("no" recomputes everything; "yes" only continues an interrupted run: a product file is reused only if its companion file <name>.input holds the fingerprint of the current inputs, the amplitudes of step 1 and this step file and of the FORM programs in form/; a file written for other inputs, or without that companion (for example, after an interrupted write), is computed again).  The records of step 1 are plain expressions (FeynFacet/Amplitudes.wl explains their vocabulary).*)


(* ::Input:: *)
LoadFeynCalc[];
class = "5";
workers = "auto";
resume = "no";
loaded = LoadAmplitudes["class-" <> class];
records = loaded["Records"];
spec = loaded["Spec"]; weights = loaded["Weights"];
spec


(* ::Section:: *)
(*Conventions*)


(* ::Text:: *)
(*As in the FeynCalc route: the two incoming partons must be quarks or antiquarks (a gluon has no transversity); the average 1/(2 Nc) per incoming parton (the 1/2 of the density (1 + gamma5 s-slash)/2 and the colour average); the normalization s/(t u) of the XX component (the axis aX is not a unit vector); the reference vectors r of the gluon polarization sums: Born p -> pa, k1 -> pa; real emission in addition k2 -> pb (pa -> pb, pb -> pa are never used: no incoming gluons).  Colour factors are written with CF = (Nc^2-1)/(2 Nc) and CA = Nc.  One unobserved momentum is not independent: k1 = pa + pb - p (Born and one loop) and k2 = pa + pb - p - k1 (real emission) inside propagators and vertices; the external spinors keep their own momentum, so the spin sums give k1-slash or k2-slash, whose scalar products are part of the kinematics.  The spin axes of the four components are set inside FORM (procedure spinVectors of common.h, selected by SPIN = XX, NN, XN, NX); the association spinAxes only documents them.*)


(* ::Input:: *)
incoming = (spec["Born"] /. None -> spec["Real"][[1, 1]])[[1 ;; 2]];       (* the two incoming partons *)
Require[FreeQ[incoming, "g"], "The transversely polarized process needs two incoming quarks or antiquarks", incoming];
components = {"XX", "NN", "XN", "NX"};
spinAxes = <|"XX" -> <|pa -> aX, pb -> aX|>, "NN" -> <|pa -> nn, pb -> nn|>, "XN" -> <|pa -> aX, pb -> nn|>, "NX" -> <|pa -> nn, pb -> aX|>|>;
normalization[spinCase_] := If[spinCase === "XX", s/(t u), 1];
average[spinCase_] := 1/(2 Nc)^2;
colourFactors[expr_] := Factor[expr] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
bornReferences = <|p -> pa, k1 -> pa|>;
realReferences = <|p -> pa, k1 -> pa, k2 -> pb|>;
bornConservation = k1 -> pa + pb - p;
realConservation = k2 -> pa + pb - p - k1;
axisCoefficients = {u/s, t/s, 1};                      (* aX = (u/s) pa + (t/s) pb + p, as in common.h *)
{incoming, average["XX"], normalization["XX"]}


(* ::Text:: *)
(*The D-dimensional scalar products of the two kinematics (the FORM procedures bornKinematics and realKinematics of common.h contain the same table): s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2; in the real emission A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 and s4 = (k1+k2)^2 = s + t + u = -A - B - C.  Mathematica needs them for the values of the propagator denominators and of the factors 1/(l.r) of the polarization sums, which FORM keeps as symbols, and for the Gram matrix of the normal reduction.  Checks: the eliminated momenta are massless (k2^2 = 0 with C = -(s+t+u) - A - B, k1^2 = 0 at the Born point u = -s - t); and the in-plane axis is what the text says: aX.pa = aX.pb = 0, aX.aX = aX.p = -t u/s, aX.k1 = t u/s at the Born point (the same products that common.h derives from the definition of aX).*)


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
axisVector = axisCoefficients . {pa, pb, p};
Require[Simplify[{dot[axisVector, pa, bornProducts], dot[axisVector, pb, bornProducts], dot[axisVector, axisVector, bornProducts] + t u/s,
    dot[axisVector, p, bornProducts] + t u/s, dot[axisVector, k1, bornProducts] - t u/s} /. u -> -s - t] === {0, 0, 0, 0, 0}, "The in-plane axis is not what the text says"];
(* the Gram matrix of pa, pb, p and the square of the projection of a momentum x onto their plane, xpar^2 = d.G^-1.d with d = (pa.x, pb.x, p.x) *)
gramMatrix[table_] := Table[dot[x, y, table], {x, {pa, pb, p}}, {y, {pa, pb, p}}];
parallelSquare[d_List, gram_] := Factor[d . Inverse[gram] . d];
"kinematics defined"


(* ::Section:: *)
(*The squaring method*)


(* ::Text:: *)
(*The method is that of Processes/Hadronic (both routes; read the five steps there), carried out by FORM:*)
(**)
(*1. Every diagram is written once as a FORM expression, before any product, with momentum conservation applied; the conditions eps(l).l = 0 and eps(l).r = 0 for every gluon l, where r is the reference vector of the same gluon's later polarization sum, are imposed on it at the start of every job (FORM procedure transversality).*)
(*2. Colour is separated from the Lorentz algebra: each diagram is written as sum_x c_x L_x with distinct colour structures c_x, and in each job the colour factors of the pairs c_x c_y^* are evaluated by color.h before the Dirac algebra.*)
(*3. Only the products with i <= j are computed: for i < j, M_i M_j^* + M_j M_i^* = 2 Re(M_i M_j^*).*)
(*4. After the spin sums with the two transverse densities and the traces, the gluon polarization sums are multiplied in and contracted by FORM; the two epsilon tensors of the two gamma5's are contracted in four dimensions, and the spin vectors are set to the axes of the component.*)
(*5. The products (i, j) of every component are independent FORM jobs, run in parallel, each with its own result file (Results/FORM/class-<n>/pairs/).*)
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
(*process.h: the declarations of one process and what the FORM procedures need to know about it: the kinematics (Born, Virtual or Real), the momenta with (D-4)-dimensional parts, the gluons with their reference vectors (all outgoing), the indices (Lorentz indices D-dimensional, colour indices in the fundamental (NR) or adjoint (NA) representation), the auxiliary vectors and their definitions, and the transversality conditions eps(l).l = eps(l).r = 0 on both sides of the product.  For the gluon x whose momentum was eliminated by momentum conservation, eps(x).x = 0 is imposed in its eliminated form, solved for the last momentum that is not x's reference vector (as in the FeynCalc route).*)


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
  Require[FreeQ[gluons, pa | pb], "An incoming gluon in a transversely polarized process", gluons];
  usedIndices = Union[usedIndices, Flatten[{Symbol["mu" <> ToString[#]], Symbol["nu" <> ToString[#]]} & /@ gluons]];
  lorentz = Select[usedIndices, StringMatchQ[SymbolName[#], ("Lor" | "cLor" | "mu" | "nu") ~~ ___] &];
  colour = Complement[usedIndices, lorentz];
  StringJoin["* process.h: ", name, " (written by 2-Algebra-FORM.wl)\n",
    "#define KINEMATICS \"", kinematics, "\"\n",
    "#define HATTED \"", StringRiffle[SymbolName /@ hatted, ","], "\"\n",
    "* gluons (all outgoing) and the reference vectors of their polarization sums\n",
    "#define INGLUONS \"\"\n",
    "#define OUTGLUONS \"", StringRiffle[SymbolName /@ gluons, ","], "\"\n",
    Table["#define REF" <> SymbolName[l] <> " \"" <> SymbolName[references[l]] <> "\"\n", {l, gluons}],
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
(*Each process (the Born process, the one-loop interference, every real-emission process, every channel of the Born table, and the process of the gauge check) gets a directory Results/FORM/class-<n>/work/<process>/ with process.h, the diagrams amp<i>.inc and the conjugate diagrams bar<j>.inc.  The values of the denominator symbols den<n> = 1/q^2 and of the factors invlr<l> = 1/(l.r) of the polarization sums are kept (symbolValues) for the reading of the results.*)


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
(*The jobs.  A job is one product: diagram I of a process times the conjugate of diagram J, for one component (SPIN = XX, NN, XN or NX).  At most `workers` FORM processes run at the same time; each writes its result (the list {real part, imaginary part}) under its name in the working directory of its process, from where it is moved to Results/FORM/class-<n>/pairs/<name>.m, next to the fingerprint of its inputs, <name>.m.input; its output to <name>.log (kept only if the job fails, and shown in the error message), and its running time in seconds to <name>.time.  The FORM program square.frm finds common.h in form/, color.h in the FORM colour package, and the process files in its working directory.*)


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
    "-d", "SPIN=" <> job["Spin"], "-d", "I=" <> ToString[job["I"]], "-d", "J=" <> ToString[job["J"]], "-d", "OUT=" <> job["Name"],
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
(*Before the squares: checks.frm verifies the forms of the HVBM gamma5 and of the transverse density used in common.h: Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = -4 i e_(a,b,c,d) (the normalization, the same as in the FeynCalc route), q-slash gamma5 = -i/6 e_(q,a,b,c) gamma^a gamma^b gamma^c inside a trace with D-dimensional momenta, and the density pa-slash gamma5 sa-slash in this short form against the four-matrix gamma5 inside a trace with both densities; FORM must print Gamma5Checks = 0.*)


(* ::Input:: *)
checkRun = RunProcess[{ToolPath["FORM"], "-p", StringRiffle[{formDirectory, DirectoryName[ToolPath["FORMColor"]]}, ":"], FileNameJoin[{formDirectory, "checks.frm"}]}];
Require[checkRun["ExitCode"] === 0 && StringContainsQ[checkRun["StandardOutput"], "Gamma5Checks = 0;"], "The gamma5 checks of checks.frm failed", checkRun["StandardOutput"]];
"gamma5 checks passed"


(* ::Text:: *)
(*The square of a process for one component: the jobs of all pairs i <= j, and their sum (diagonal products as they are, i < j twice the real part), with the values of den<n> and invlr<l> put in.  A diagonal product |M_i|^2 is real; its imaginary part must vanish exactly at a rational point (as in the FeynCalc route; the point includes the normal component nk1).  The normal reduction of a square (odd powers of nk1 dropped, nk1^2 = k1par^2 + spe3) is the normalReduction of the FeynCalc route, applied below to the real-emission squares (with nk1) and to the one-loop interference (with nloop).*)


(* ::Input:: *)
checkPoint = {s -> 1, t -> -3/10, u -> -2/5, A -> -1/7, B -> -2/9, C -> -3/10 + 1/7 + 2/9, D -> 37/10, Nc -> 3, spe3 -> 1/11, nk1 -> 3/13, gs -> 1, Nf -> 5};
squareJobs[name_, n_, spinCase_] := Flatten[Table[<|"Name" -> name <> "-" <> spinCase <> "-" <> ToString[i] <> "-" <> ToString[j],
  "Process" -> name, "I" -> i, "J" -> j, "Spin" -> spinCase|>, {i, n}, {j, i, n}]];
square[name_, n_, spinCase_] := Module[{jobs = squareJobs[name, n, spinCase], results, values = symbolValues[name]},
  results = runJobs[jobs];
  Require[And @@ MapThread[#2["I"] != #2["J"] || (#1[[2]] /. values /. checkPoint) === 0 &, {results, jobs}], "A diagonal product has a nonzero imaginary part", name];
  Total[MapThread[If[#2["I"] == #2["J"], 1, 2] #1[[1]] &, {results, jobs}]] /. values];
(* the normal component nx of an unobserved or loop momentum x: polynomial in nx; odd powers integrate to zero, nx^2 = xpar^2 - x4.x4
   with xpar^2 = d.G^-1.d (d the products of x with pa, pb, p; G their Gram matrix) and x4.x4 the four-dimensional square of x *)
normalReduction[expr_, nx_, parallel_, fourSquare_] := Module[{coefficients = CoefficientList[Expand[expr], nx]},
  If[coefficients === {}, coefficients = {0}];                                   (* a vanishing expression (a graph can vanish exactly) *)
  Require[Length[coefficients] <= 3, "A power of the normal component above two appeared", Length[coefficients] - 1];
  <|"Even" -> coefficients[[1]] + If[Length[coefficients] == 3, coefficients[[3]] (parallel - fourSquare), 0],
    "Odd" -> If[Length[coefficients] >= 2, coefficients[[2]], 0]|>];
"parallel FORM jobs: " <> ToString[nWorkers]


(* ::Section:: *)
(*Born*)


(* ::Text:: *)
(*The Born square with the reference vectors p -> pa, k1 -> pa for the gluons, averaged over the initial colours (with the 1/2 of each density), normalized, with u = -s - t, for the four components; then the checks XN = NX = 0 (exactly) and XX + NN = 0 in four dimensions, and the two projections TT = (XX - NN)/2 and Symmetric = (XX + NN)/2.  The evanescent symmetric part is displayed.*)


(* ::Input:: *)
bornDiagrams = prepareDiagram[#, bornConservation, {}] & /@ records["Born"];
If[bornDiagrams =!= {}, writeProcess["born", "Born", bornReferences, bornConservation, spec["Born"], bornDiagrams, bornDiagrams]];
bornComponents = If[bornDiagrams === {}, Association[# -> 0 & /@ components],
  Association @ Table[sc -> Factor[colourFactors[Together[average[sc] normalization[sc] square["born", Length[bornDiagrams], sc]]] /. u -> -s - t], {sc, components}]];
Require[FreeQ[bornComponents, Complex | nk1], "The Born square is not a real scalar", bornComponents];
Require[bornComponents["XN"] === 0 && bornComponents["NX"] === 0, "The mixed components of the Born square do not vanish", bornComponents];
Require[Factor[(bornComponents["XX"] + bornComponents["NN"]) /. D -> 4] === 0, "XX + NN of the Born square does not vanish in four dimensions", bornComponents];
born = <|"TT" -> Factor[(bornComponents["XX"] - bornComponents["NN"])/2], "Symmetric" -> Factor[(bornComponents["XX"] + bornComponents["NN"])/2],
  "XX" -> bornComponents["XX"], "NN" -> bornComponents["NN"]|>;
born


(* ::Section:: *)
(*The transversely polarized Born table*)


(* ::Text:: *)
(*The D-dimensional Born squares |M|^2/gs^4 (colours and spins of the unobserved parton summed, colours of a and b averaged, the two transverse densities inserted) of every 2 -> 2 channel a b -> c that the collinear counterterms of the five classes need: q qb -> q' (and qb'), q q -> q, q qb -> q (and qb), q qb -> g, their a <-> b relabellings (qb q -> ...), and the charge conjugate qb qb -> qb; plus the control channels q q' -> q, q q' -> q', q qb' -> q, q qb' -> qb', which must vanish (classes 1-4: the two polarized partons on different fermion lines).  The diagrams are generated here with FeynArts (the same fields, model and exclusions as in the FeynCalc route), converted into records with AmplitudeRecords (as step 1 does for the class), and squared by FORM with the machinery above, each channel as a Born process with the unobserved parton k1 and the reference vector pa for the gluons.  Keys as in Processes/PartonicBorn ("q qb -> qp"); each entry holds XX, NN, TT and Symmetric.*)
(**)
(*Checks (those of the FeynCalc route): XN = NX = 0 and XX + NN = 0 in four dimensions for every channel; the a <-> b relabelled channel equals the original with t <-> u (the axis aX is symmetric, the normal changes sign for both densities); the charge conjugate equals its representative; the control channels vanish; the class's own Born square (the diagrams of step 1) agrees with its table entry; and at D = 4, Nc = 3 the TT entries reproduce Table I of Soffer, Stratmann and Vogelsang: with sigma = s E dsigma/d^3p = 2 s^2 dsigma/(dt du dPhi) and dsigma/(dt du dPhi) = dsigma/(dt dPhi) delta(s4), the cos(2 Phi) coefficient of their Eq. (6) corresponds to |M|^2_TT/gs^4 = (4/9) (t u/s^2) delta|M|^2.*)


(* ::Input:: *)
tableChannels = {{"q", "qb", "qp"}, {"q", "qb", "qpb"}, {"q", "q", "q"}, {"q", "qb", "q"}, {"q", "qb", "qb"}, {"q", "qb", "g"},
  {"qb", "q", "qp"}, {"qb", "q", "qpb"}, {"qb", "q", "q"}, {"qb", "q", "qb"}, {"qb", "q", "g"}, {"qb", "qb", "qb"}};
controlChannels = {{"q", "qp", "q"}, {"q", "qp", "qp"}, {"q", "qpb", "q"}, {"q", "qpb", "qpb"}};
field = <|"q" -> F[3, {1}], "qb" -> -F[3, {1}], "qp" -> F[3, {2}], "qpb" -> -F[3, {2}], "g" -> V[5]|>;
channelName[{a_, b_, c_}] := a <> " " <> b <> " -> " <> c;
flavourNumber[x_] := Switch[x, "q", {1, 0}, "qb", {-1, 0}, "qp", {0, 1}, "qpb", {0, -1}, "g", {0, 0}];
partner[{a_, b_, c_}] := Module[{ok = Select[Keys[field], Total[flavourNumber /@ {a, b}] === Total[flavourNumber /@ {c, #}] &]},
  Require[Length[ok] == 1, "The unobserved parton of " <> channelName[{a, b, c}] <> " is not unique", ok]; First[ok]];
tableDiagrams[ch_] := InsertFields[CreateTopologies[0, 2 -> 2], {field[ch[[1]]], field[ch[[2]]]} -> {field[ch[[3]]], field[partner[ch]]},
  InsertionLevel -> {Classes}, Model -> "SMQCD", ExcludeParticles -> {S[_], U[_], V[1 | 2 | 3], F[1 | 2 | 4], F[3, {3}]}];
(* the amplitudes of a channel as records, one list per diagram, with the unobserved parton k1 (as in step 1) *)
tableRecords[ch_] := Module[{amps},
  amps = FCFAConvert[CreateFeynAmp[tableDiagrams[ch], Truncated -> False, PreFactor -> 1],
    IncomingMomenta -> {pa, pb}, OutgoingMomenta -> {p, k1}, UndoChiralSplittings -> True, ChangeDimension -> D,
    List -> True, SMP -> False, Contract -> False, DropSumOver -> True,
    FinalSubstitutions -> {FCGV["MU"] -> 0, FCGV["MC"] -> 0}] /. SMP["g_s"] -> gs /. GaugeXi[_] -> 1;
  amps = AmplitudeRecords /@ amps;
  Require[FreeQ[amps, $Failed], "A record of the Born table could not be formed", channelName[ch]]; amps];
tableName[ch_] := "table-" <> StringRiffle[ch, "_"];
tableEntry[ch_] := Module[{diagrams = prepareDiagram[#, bornConservation, {}] & /@ tableRecords[ch], xx, nnv, xn, nx},
  writeProcess[tableName[ch], "Born", bornReferences, bornConservation, Append[ch, partner[ch]], diagrams, diagrams];
  {xx, nnv, xn, nx} = Table[Factor[colourFactors[Together[average[sc] normalization[sc] square[tableName[ch], Length[diagrams], sc]/gs^4]] /. u -> -s - t], {sc, components}];
  Require[FreeQ[{xx, nnv, xn, nx}, Complex | nk1], "A table entry is not a real scalar", channelName[ch]];
  Require[xn === 0 && nx === 0, "The mixed components of a table entry do not vanish", {channelName[ch], xn, nx}];
  Require[Factor[(xx + nnv) /. D -> 4] === 0, "XX + NN of a table entry does not vanish in four dimensions", {channelName[ch], xx, nnv}];
  <|"XX" -> xx, "NN" -> nnv, "TT" -> Factor[(xx - nnv)/2], "Symmetric" -> Factor[(xx + nnv)/2]|>];
tableEntries = Association @ Table[channelName[ch] -> tableEntry[ch], {ch, Join[tableChannels, controlChannels]}];
bornTable = Association @ Table[key -> Map[#[key] &, tableEntries], {key, {"TT", "Symmetric", "XX", "NN"}}];
TableForm[Normal[Map[Simplify[# /. D -> 4] &, bornTable["TT"]]]]


(* ::Input:: *)
tableCheck[expr_] := Factor[expr /. u -> -s - t];
Require[And @@ Table[tableCheck[bornTable["TT"][channelName[ch]]] === 0 && tableCheck[bornTable["Symmetric"][channelName[ch]]] === 0, {ch, controlChannels}],
  "A control channel of the Born table does not vanish", Map[bornTable["TT"], channelName /@ controlChannels]];
(* a <-> b: pa <-> pb exchanges t and u, i.e. t -> -s - t in the entries (u is eliminated) *)
swapped[name_] := Module[{parts = StringSplit[name, " "]}, parts[[2]] <> " " <> parts[[1]] <> " -> " <> parts[[4]]];
relabelDifferences = Association @ Flatten[Table[{key, name} -> tableCheck[bornTable[key][swapped[name]] - (bornTable[key][name] /. t -> -s - t)],
    {key, {"TT", "Symmetric"}}, {name, channelName /@ Select[tableChannels, #[[1]] === "q" && #[[2]] === "qb" &]}]];
Require[And @@ (# === 0 & /@ Values[relabelDifferences]), "An a <-> b relabelled entry is not the t <-> u image of the original", Select[relabelDifferences, # =!= 0 &]];
Require[tableCheck[bornTable["TT"]["qb qb -> qb"] - bornTable["TT"]["q q -> q"]] === 0 && tableCheck[bornTable["Symmetric"]["qb qb -> qb"] - bornTable["Symmetric"]["q q -> q"]] === 0,
  "The charge-conjugate channel differs from its representative"];
(* the class's own Born square (the diagrams of step 1) against the table (the diagrams generated here) *)
If[spec["Born"] =!= None,
  Require[And @@ Table[tableCheck[born[key]/gs^4 - bornTable[key][StringRiffle[spec["Born"][[1 ;; 2]], " "] <> " -> " <> spec["Born"][[3]]]] === 0, {key, {"TT", "Symmetric", "XX", "NN"}}],
    "The Born square of the class differs from its table entry", born]];
(* the literature: Table I of Soffer, Stratmann and Vogelsang at D = 4, Nc = 3 *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "SSV-born.wl"}]];
literatureDifferences = Association @ Table[name -> Factor[(bornTable["TT"][name] /. {D -> 4, CA -> 3, CF -> 4/3} /. u -> -s - t) - ((4/9) (t u/s^2) SSVMatrixElement[name] /. u -> -s - t)],
  {name, Keys[SSVMatrixElement]}];
Require[And @@ (# === 0 & /@ Values[literatureDifferences]), "A transversely polarized Born square differs from Soffer, Stratmann and Vogelsang", literatureDifferences];
"Born table checked: control channels vanish, a <-> b relabelling, charge conjugation, the class Born, and the literature at D = 4, Nc = 3"


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
(*One job per graph and component: the graph times the conjugate of the summed Born amplitude, M1_graph M0*, with its real and imaginary parts.  FORM returns the loop momentum through the symbols L2 = loop.loop, Lpa = loop.pa, Lpb = loop.pb, Lp = loop.p and Lhat2 = loophat.loophat (the four-dimensional loop4.loop4 = L2 - Lhat2 in the HVBM scheme), and its normal component through nloop = n.loop (the product aX.loop was already expanded into Lpa, Lpb, Lp by FORM).  First the normal reduction: the odd powers of nloop are dropped (the tensor loop integrals are built from pa, pb, p only; the dropped parts must vanish for XN and NX, whose even parts are checked below) and nloop^2 = looppar^2 - (L2 - Lhat2) with looppar^2 = d.G^-1.d, d = (Lpa, Lpb, Lp).  Then the four products are written through the inverse propagators D_j of the graph's family (the four equations loop^2 = D1, 2 loop.q_j = D_j - D_1 - q_j^2 are solved for L2, Lpa, Lpb, Lp), and every term is read off as coefficient x family[indices] x HatPower[n], with n the power of Lhat2 (it shifts the dimension of the integral in step 3, not the dimension D of the coefficients).*)
(**)
(*Normalization, as in the FeynCalc route: FeynArts writes the one-loop amplitude with the measure d^D loop/(2 Pi)^D; the scalar integrals are defined with the measure d^D loop/(i Pi^(D/2)), so that M1 = i/(4 Pi)^(D/2) x (sum of coefficients times integrals).  The products (coefficient x M0*) are real, and 2 Re(M1 M0*) = (2/(4 Pi)^(D/2)) Sum_j [i c_j M0*] Re(I_j), with the real parts of the integrals in the physical region (s > 0) taken in step 3.  The association virtual below holds Sum_j [i c_j M0*] I_j for the two projections (i times the complex product of FORM, which must come out real); the factor 2/(4 Pi)^(D/2) is included in step 4.  Checks: XN and NX vanish (their even parts in nloop), and XX + NN vanishes in four dimensions once the hatted terms (which integrate to order eps) are dropped.*)


(* ::Input:: *)
inverseProps[name_] := Module[{q = families[name], sol},
  sol = Solve[Table[Dj[j] == L2 + 2 (Coefficient[q[[j]], pa] Lpa + Coefficient[q[[j]], pb] Lpb + Coefficient[q[[j]], p] Lp) + dot[q[[j]] - loop, q[[j]] - loop, bornProducts] /. u -> -s - t, {j, 4}], {L2, Lpa, Lpb, Lp}];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators"]; First[sol]];
loopRules = Association @ Table[name -> inverseProps[name], {name, Keys[families]}];
loopParallelSquare = parallelSquare[{Lpa, Lpb, Lp}, gramMatrix[bornProducts]];                     (* looppar^2 = d.G^-1.d *)
(* the interference of one graph -> sum of coefficient x family[indices] x HatPower[n]; the odd part in nloop is returned as its leaf count *)
toFamilies[value_, name_] := Module[{expr, reduced},
  expr = Expand[value /. symbolValues["virtual"]];
  reduced = normalReduction[expr, nloop, loopParallelSquare, L2 - Lhat2];
  expr = Expand[reduced["Even"] /. loopRules[name] /. Table[Symbol["dj" <> ToString[j]] -> 1/Dj[j], {j, 4}]];
  Require[FreeQ[expr, L2 | Lpa | Lpb | Lp | nloop], "A scalar product of the loop momentum was left"];
  <|"Even" -> Total[Table[With[{powers = Table[Exponent[term, Dj[j]], {j, 4}], hat = Exponent[term, Lhat2]},
      (term/(Times @@ Table[Dj[j]^powers[[j]], {j, 4}] Lhat2^hat)) familyHeads[name] @@ (-powers) HatPower[hat]], {term, DeleteCases[termsOf[expr], 0]}]],
    "OddLeafCount" -> LeafCount[reduced["Odd"]]|>];
virtualJobs[spinCase_] := Table[<|"Name" -> "virtual-" <> spinCase <> "-" <> ToString[g], "Process" -> "virtual", "I" -> g, "J" -> 1, "Spin" -> spinCase|>, {g, Length[virtualGraphs]}];
oddLeafCounts = <||>;
virtualComponents = If[virtualGraphs === {}, Association[# -> 0 & /@ components],
  Association @ Table[spinCase -> Module[{results = runJobs[virtualJobs[spinCase]], graphs},
    graphs = MapThread[toFamilies[#1[[1]] + I #1[[2]], #2[[1]]] &, {results, routings}];
    oddLeafCounts[spinCase] = Total[#["OddLeafCount"] & /@ graphs];
    Collect[I average[spinCase] normalization[spinCase] Total[#["Even"] & /@ graphs] /. u -> -s - t,
      {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, colourFactors]], {spinCase, components}]];
Require[FreeQ[virtualComponents, Complex | nloop], "The virtual interference is not real or still contains n.loop", Short[virtualComponents, 5]];
Require[virtualComponents["XN"] === 0 && virtualComponents["NX"] === 0, "The mixed components of the virtual interference do not vanish", {virtualComponents["XN"], virtualComponents["NX"]}];
Require[Collect[(virtualComponents["XX"] + virtualComponents["NN"]) /. HatPower[n_ /; n > 0] -> 0 /. HatPower[0] -> 1 /. D -> 4, {_KVirtualST, _KVirtualSU, _KVirtualTU}, Factor] === 0,
  "XX + NN of the virtual interference does not vanish in four dimensions"];
virtual = <|"TT" -> Collect[(virtualComponents["XX"] - virtualComponents["NN"])/2, {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, Factor],
  "Symmetric" -> Collect[(virtualComponents["XX"] + virtualComponents["NN"])/2, {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, Factor]|>;
{Length[virtual["TT"]], Union[Cases[virtual, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]] // Length, Map[Union[Cases[#, _HatPower, Infinity]] &, virtual],
  "leaf count of the dropped odd parts in n.loop" -> oddLeafCounts}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*For every real-emission process of the class, the square of the sum of its diagrams for the four components, in terms of the invariants s, t, u, A, B, C, the hatted invariant spe3 and the normal component nk1; the processes are added with their weights (classes.wl: 1/2 for two identical unobserved partons, Nf - 1 for the unobserved quark pairs of other flavours).  Reference vectors of the physical polarization sums: pa for the gluons p and k1, pb for k2.  Then the normal reduction: the parts odd in nk1 are dropped (they integrate to zero), nk1^2 is replaced by k1par^2 + spe3.  As in the FeynCalc route the sum is not brought to a common denominator (step 3 expands it into partial fractions anyway); the colour factors are polynomials in Nc = CA.  Checks: the even parts of XN and NX vanish, and XX + NN vanishes in four dimensions without the hatted terms; both are evaluated exactly at random rational points of all invariants (the expressions have tens of thousands of terms).  In the HVBM scheme the squares are at most linear in spe3.*)


(* ::Input:: *)
processName[key_] := "real-" <> StringReplace[key, " " -> "_"];
realDiagrams = Association @ KeyValueMap[#1 -> (prepareDiagram[#, realConservation, {}] & /@ #2) &, records["Real"]];
KeyValueMap[writeProcess[processName[#1], "Real", realReferences, realConservation, StringSplit[#1], #2, #2] &, realDiagrams];
realSquares = Association @ Table[spinCase -> Association @ KeyValueMap[#1 -> square[processName[#1], Length[#2], spinCase] &, realDiagrams], {spinCase, components}];
realParallelSquare = parallelSquare[{-A/2, -B/2, C/2}, gramMatrix[realProducts]];      (* k1par^2 = d.G^-1.d with d = (pa.k1, pb.k1, p.k1) *)
realReduced = Association @ Table[spinCase -> Association @ KeyValueMap[
    #1 -> normalReduction[#2, nk1, realParallelSquare, -spe3] &, realSquares[spinCase]], {spinCase, components}];       (* k1.k1 (four-dimensional) = -spe3: nk1^2 = k1par^2 + spe3 *)
realComponents = Association @ Table[spinCase -> (average[spinCase] normalization[spinCase] Total[KeyValueMap[weights[#1] #2["Even"] &, realReduced[spinCase]]] /. Nc -> CA), {spinCase, components}];
Require[FreeQ[realComponents, Complex | nk1], "The real-emission square is not real or still contains the normal component"];
{realParallelSquare, Map[LeafCount, realComponents]}


(* ::Input:: *)
(* exact rational points of all invariants (C = -s4 - A - B), D, Nc, spe3 *)
realPoints = Table[With[{r = RandomInteger[{1, 97}, 8]/RandomInteger[{98, 197}, 8]},
  {s -> 1, t -> -r[[1]], u -> -r[[2]] (1 - r[[1]]), A -> -r[[3]], B -> -r[[4]], D -> 4 - r[[5]], CA -> 2 + 10 r[[6]], CF -> 1 + r[[7]], spe3 -> r[[8]], Nf -> 5}], {3}];
atPoint[expr_, point_] := expr /. C -> -(s + t + u) - A - B /. point;
Require[And @@ Table[atPoint[realComponents["XN"], point] === 0 && atPoint[realComponents["NX"], point] === 0, {point, realPoints}],
  "The even parts of the mixed components of the real emission do not vanish"];
Require[And @@ Table[atPoint[(realComponents["XX"] + realComponents["NN"]) /. {D -> 4, spe3 -> 0}, point] === 0, {point, realPoints}],
  "XX + NN of the real emission does not vanish in four dimensions"];
real = <|"TT" -> (realComponents["XX"] - realComponents["NN"])/2, "Symmetric" -> (realComponents["XX"] + realComponents["NN"])/2|>;
Require[Max[Exponent[#, spe3] & /@ termsOf[Expand[real["TT"] + real["Symmetric"]]]] <= 1, "The real-emission square has spe3^2 terms (step 3 has the masters in D + 2 only)"];
{"odd parts dropped (leaf counts)" -> Map[Total[KeyValueMap[LeafCount[#2["Odd"]] &, #]] &, realReduced], LeafCount /@ real}


(* ::Text:: *)
(*Check of the gauge invariance of the real-emission square: for the first process of the class and the XX component (the one without normal components), the physical polarization sums with other reference vectors (pb for p and k1, pa for k2) give the same result, exactly, at three random rational points (all invariants, D and Nc, with C = -s4 - A - B).  The diagrams are written again (process "gauge"), because the transversality conditions eps(l).r = 0 depend on the reference vectors.*)


(* ::Input:: *)
alternativeReferences = <|p -> pb, k1 -> pb, k2 -> pa|>;
firstProcess = First[Keys[realDiagrams]];
writeProcess["gauge", "Real", alternativeReferences, realConservation, StringSplit[firstProcess], realDiagrams[firstProcess], realDiagrams[firstProcess]];
alternative = square["gauge", Length[realDiagrams[firstProcess]], "XX"];
gaugePoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], Nc -> 2 + 10 v[[6]], spe3 -> v[[7]]}], {3}];
referenceDifference = Table[(alternative - realSquares["XX"][firstProcess]) /. C -> -(s + t + u) - A - B /. point, {point, gaugePoints}];
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
(*The result file Results/FORM/class-<n>/Algebra.wl, in the format of the FeynCalc route: the Born square, the virtual interference and the real-emission square (projections TT and Symmetric), and the Born table for steps 3 and 4.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Virtual" -> virtual, "Real" -> real, "BornTable" -> bornTable, "Class" -> class, "Spec" -> spec,
  "Families" -> families, "Projections" -> "TT = (XX - NN)/2 (the cos(2 Phi) coefficient, the hard part of A_TT); Symmetric = (XX + NN)/2 (evanescent in the HVBM scheme)",
  "Variables" -> "s, t, u = -s-t; A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 = -s4-A-B; spe3; loopHat2 via HatPower[n]; Born table: |M|^2/gs^4 in s, t, D, CF, CA"|>;
SaveResult["Algebra", algebra, "FORM/class-" <> class];
Keys[algebra]


(* ::Section:: *)
(*Comparison with the FeynCalc route*)


(* ::Text:: *)
(*If the FeynCalc route has been run for this class, the two files must agree (otherwise the comparison is skipped and the cell says so): Born (TT, Symmetric, XX, NN), Virtual and Real (TT, Symmetric) and every entry of the Born table (all four keys) are compared at two exact random rational points of s, t, u, A, B, D, CA, spe3, Nf, gs (with C = -(s+t+u) - A - B and CF = (CA^2 - 1)/(2 CA)); in the virtual every scalar integral is replaced by its own, distinct random rational number and HatPower[n] by hp^n, so that the coefficient of every integral is compared.  The arithmetic is exact, but the comparison is sampled at two points, not a symbolic identity.*)


(* ::Input:: *)
feynCalcFile = FileNameJoin[{ResultsDirectory["FeynCalc/class-" <> class], "Algebra.wl"}];
routeComparison = If[!FileExistsQ[feynCalcFile], "the FeynCalc route has not been run for this class",
  Module[{other = Get[feynCalcFile], integrals, integralValues, points, evaluate, parts},
    integrals = Union[Cases[{algebra["Virtual"], other["Virtual"]}, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]];
    integralValues = Table[With[{v = RandomInteger[{1, 10^9}, Length[integrals]]/RandomInteger[{1, 10^9}, Length[integrals]]},
      Require[DuplicateFreeQ[v], "The random values of the integrals are not distinct"]; v], {2}];
    points = Table[Join[Thread[{s, t, u, A, B, D, CA, spe3, Nf, gs, hp} -> RandomInteger[{1, 10^6}, 11]/RandomInteger[{1, 10^6}, 11] {1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1}],
      Thread[integrals -> integralValues[[k]]]], {k, 2}];
    evaluate[x_, point_] := Together[x /. HatPower[n_] :> hp^n /. point /. C -> -(s + t + u) - A - B /. CF -> (CA^2 - 1)/(2 CA) /. point];
    parts = Join[Flatten[Table[{part, key} -> {algebra[part][key], other[part][key]}, {part, {"Born", "Virtual", "Real"}}, {key, Keys[algebra[part]]}]],
      Flatten[Table[{"BornTable", key, name} -> {algebra["BornTable"][key][name], other["BornTable"][key][name]}, {key, Keys[algebra["BornTable"]]}, {name, Keys[algebra["BornTable"]["TT"]]}]]];
    Require[Keys[other["BornTable"]["TT"]] === Keys[algebra["BornTable"]["TT"]], "The two routes have different Born tables"];
    Association @ Table[item[[1]] -> (And @@ Table[evaluate[item[[2, 1]] - item[[2, 2]], point] === 0, {point, points}]), {item, parts}]]];
Require[StringQ[routeComparison] || And @@ Values[routeComparison], "The FORM route differs from the FeynCalc route", Select[routeComparison, ! # &]];
If[StringQ[routeComparison], routeComparison, Counts[Values[routeComparison]]]
