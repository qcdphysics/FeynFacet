(* ::Title:: *)
(*Transverse spin transfer in proton-proton collisions, step 2 (FORM route): spin and colour algebra*)


(* ::Text:: *)
(*Input: the amplitude records of step 1 for the class chosen below.  Output: Results/FORM/class-<n>/Algebra.wl, with the same content, conventions and format as the file of the FeynCalc route (2-Algebra-FeynCalc.wl, read it first: it defines the observable, the spin density of the polarized incoming parton a, the spin projection of the observed parton c, the axes and the four components):*)
(*  Born:      the Born transfer square in the two transverse components XX and NN, summed over the final spins and colours of the unobserved partons, averaged over the colours of a and over the spin and colour of b, in D dimensions;*)
(*  Virtual:   the interference 2 Re(M_1loop M_Born^*) with the same projections, written as Sum_j [i c_j M0*] I_j in terms of the scalar one-loop integrals of the families VirtualST, VirtualSU, VirtualTU, with the powers of the (D-4)-dimensional part of the loop momentum as HatPower[n];*)
(*  Real:      the sum over the real-emission processes of the class, each weighted as in classes.wl, of the projected |M|^2 in the invariants s, t, u, A, B, C and the hatted invariant spe3;*)
(*  BornTable: the D-dimensional transfer squares |M|^2/gs^4 of the 2 -> 2 channels of classes.wl (TransferBornChannels), the lower Born terms of the collinear counterterms of step 4.*)
(*The last section checks that the file agrees exactly with the one of the FeynCalc route.*)
(**)
(*Spin projections (docs/TRANSVERSE-PP.md).  The incoming parton a (momentum pa, a quark or an antiquark) carries a transverse spin s_a, inserted as the density pa-slash gamma5 sa-slash/2 (the same expression for a quark and for an antiquark); the transverse spin s_c of the observed parton c (momentum p) is measured, inserted as the projection p-slash gamma5 sc-slash (no 1/2: the partonic quantity is sigma(s_c) - sigma(-s_c)); parton b is unpolarized.  HVBM gamma5, four-dimensional spin vectors.  The spin vectors are written in the basis of the scattering plane of pa, pb, p: s_a = cos(phi_a) eX + sin(phi_a) n and s_c = cos(phi_c) eC + sin(phi_c) n, with the in-plane axis eX = (u/s) pa + (t/s) pb + p of the incoming partons (eX.pa = eX.pb = 0, eX^2 = -t u/s), the in-plane axis eC = -u/(t+u) pa + t/(t+u) pb - s (u-t)/(t+u)^2 p of the observed parton (eC.p = 0, eC.(pa+pb) = 0, eC^2 = -s t u/(t+u)^2) and the unit normal n (n.pa = n.pb = n.p = 0, n^2 = -1).  The squared matrix element is XX cos cos + NN sin sin + XN cos sin + NX sin cos; the four components are computed separately (XX with the unnormalized axes, divided by the product of their lengths -t u/(t+u); XN and NX must vanish by parity, up to positive factors), and XX and NN are both kept: they coincide at Born level in four dimensions and differ at NLO.*)
(**)
(*Normal components.  n.k1 = -n.k2 = nk1 appears in the components with a spin vector along the normal of the real emission; odd powers integrate to zero, and nk1^2 = k1par^2 + spe3 with k1par^2 = d.G^-1.d, G the Gram matrix of pa, pb, p and d = (pa.k1, pb.k1, p.k1) (Processes/Collins, normalReduction).  The loop momentum has a normal component nloop = n.loop too, reduced by the symmetric integration over the (D-3)-dimensional space orthogonal to the scattering plane (section "Virtual", exactly as in the FeynCalc route).  Both reductions are done here, in Mathematica, on the results of FORM.*)
(**)
(*Division of the work.  Mathematica reads the amplitude records of step 1 (the class's diagrams and the Born diagrams of the transfer channels of the table), applies momentum conservation, splits every diagram into its colour structures and Lorentz parts and writes it as a FORM expression; it starts the FORM jobs, reads their results, puts in the kinematic values of the propagator denominators and reduces the normal components.  FORM does all the spin and colour algebra, including the transverse density and projection.  The FORM programs are in the directory form/ next to this file and are meant to be read:*)
(*  form/common.h:   declarations; the spin sums, the density of a and the projection of c (with the HVBM gamma5 and the symbolic spin vectors sa, sc); the gluon polarization sums; the four-dimensional contraction of the two epsilon tensors; the substitution of the axes eXs, eCs, n for the spin vectors; the kinematics;*)
(*  form/square.frm: one product M_I M_J^* of two diagrams for one component: transversality of the gluon polarizations, colour with the package color.h, traces, polarization sums, epsilon contraction, axes, kinematics;*)
(*  form/checks.frm: the gamma5 normalization, the short form of q-slash gamma5, and the short forms of the density and the projection used in common.h.*)
(**)
(*Invariants as in Processes/Hadronic: s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2; real emission: s4 = (k1+k2)^2 = s + t + u, A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 with A + B + C = -s4; spe3 = k1hat.k1hat.*)


(* ::Section:: *)
(*Load*)


(* ::Text:: *)
(*The package FeynFacet: the checks Require, the reading of the amplitudes of step 1, the location of the FORM program, and the saving of the results.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*Choices: the class; the number of FORM jobs that run at the same time ("auto" = one per processor core, at most 8; "1" = one job at a time); and whether the product files of an earlier run of this same, unchanged file may be reused ("no" recomputes everything; "yes" only continues an interrupted run: a product file is reused only if its companion file <name>.input holds the fingerprint of the current inputs, the amplitudes of step 1 and this step file and of the FORM programs in form/; a file written for other inputs, or without that companion (for example, after an interrupted write), is computed again).  FeynCalc is not needed: the records of step 1 are plain expressions (FeynFacet/Amplitudes.wl explains their vocabulary), and the Born diagrams of the transfer channels of the table are among them.  The flavour number Nf of the closed quark loops is written by step 1 with FeynCalc's symbol; it is renamed to the plain Nf of the result files.*)


(* ::Input:: *)
class = "1";
workers = "auto";
resume = "no";
loaded = LoadAmplitudes["class-" <> class] /. Symbol["FeynCalc`Nf"] -> Nf;
records = loaded["Records"];
spec = loaded["Spec"]; weights = loaded["Weights"];
spec


(* ::Section:: *)
(*Conventions*)


(* ::Text:: *)
(*As in the FeynCalc route: the two incoming partons, of which a must be a quark or an antiquark (a gluon has no transversity); the averages 1/(2 Nc) for a (the 1/2 of its density and its colour average) and for a quark b, 1/((D-2)(Nc^2-1)) for a gluon b; the normalization of the XX component (FORM uses the scaled axes eXs = s eX and eCs = (t+u)^2 eC with polynomial coefficients: XX is divided by s (t+u)^2 and by the product of the lengths of eX and eC, -t u/(t+u), in all by -s t u (t+u)); the reference vectors r of the gluon polarization sums: Born pa -> pb, pb -> pa, p -> pa, k1 -> pa; real emission in addition k2 -> pb.  Colour factors are written with CF = (Nc^2-1)/(2 Nc) and CA = Nc.  One unobserved momentum is not independent: k1 = pa + pb - p (Born and one loop) and k2 = pa + pb - p - k1 (real emission) inside propagators and vertices; the external spinors keep their own momentum, so the spin sums give k1-slash or k2-slash, whose scalar products are part of the kinematics.  The spin axes of the four components are set inside FORM (procedure spinVectors of common.h, selected by SPIN = XX, NN, XN, NX); the association spinAxes only documents them.*)


(* ::Input:: *)
incoming = (spec["Born"] /. None -> spec["Real"][[1, 1]])[[1 ;; 2]];       (* the two incoming partons *)
Require[incoming[[1]] =!= "g", "The polarized parton a must be a quark or an antiquark", incoming];
spinCases = {"XX", "NN", "XN", "NX"};
outputSpins = {"XX", "NN"};
spinAxes = <|"XX" -> {eX, eC}, "NN" -> {nn, nn}, "XN" -> {eX, nn}, "NX" -> {nn, eC}|>;     (* {s_a, s_c} of the four components *)
axisLengths = -t u/(t + u);                                                                   (* |eX| |eC| = sqrt(eX^2 eC^2), t, u < 0 *)
normalization[spinCase_] := If[spinCase === "XX", 1/(s (t + u)^2 axisLengths), 1];
average[in_List] := Times @@ (Switch[in[[#]], "g", 1/(D - 2)/(Nc^2 - 1), _, 1/(2 Nc)] & /@ {1, 2});
colourFactors[expr_] := Factor[expr] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;
bornReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa|>;
realReferences = <|pa -> pb, pb -> pa, p -> pa, k1 -> pa, k2 -> pb|>;
bornConservation = k1 -> pa + pb - p;
realConservation = k2 -> pa + pb - p - k1;
axisCoefficients = <|eX -> {u/s, t/s, 1}, eC -> {-u/(t + u), t/(t + u), -s (u - t)/(t + u)^2}|>;    (* on the basis {pa, pb, p}; eXs = s eX, eCs = (t+u)^2 eC in common.h *)
{incoming, average[incoming], normalization["XX"]}


(* ::Text:: *)
(*The D-dimensional scalar products of the two kinematics (the FORM procedures bornKinematics and realKinematics of common.h contain the same table): s = (pa+pb)^2, t = (pa-p)^2, u = (pb-p)^2; in the real emission A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 and s4 = (k1+k2)^2 = s + t + u = -A - B - C.  Mathematica needs them for the values of the propagator denominators and of the factors 1/(l.r) of the polarization sums, which FORM keeps as symbols, and for the Gram matrix of the normal reductions.  Checks: the eliminated momenta are massless (k2^2 = 0 with C = -(s+t+u) - A - B, k1^2 = 0 at the Born point u = -s - t); and the axes are what the text says: eX.pa = eX.pb = 0, eC.p = 0, eC.(pa+pb) = 0, eX^2 = -t u/s, eC^2 = -s t u/(t+u)^2 (the same products that common.h derives from the definitions of eXs and eCs).*)


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
axisVector[axis_] := axisCoefficients[axis] . {pa, pb, p};
Require[Simplify[{dot[axisVector[eX], pa, realProducts], dot[axisVector[eX], pb, realProducts], dot[axisVector[eC], p, realProducts], dot[axisVector[eC], pa + pb, realProducts],
    dot[axisVector[eX], axisVector[eX], realProducts] + t u/s, dot[axisVector[eC], axisVector[eC], realProducts] + s t u/(t + u)^2}] === {0, 0, 0, 0, 0, 0},
  "The axes are not what the text says"];
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
(*4. After the spin sums with the density of a and the projection of c and the traces, the gluon polarization sums are multiplied in and contracted by FORM; the two epsilon tensors of the two gamma5's are contracted in four dimensions, and the spin vectors are set to the axes of the component.*)
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
(*process.h: the declarations of one process and what the FORM procedures need to know about it: the kinematics (Born, Virtual or Real), the momenta with (D-4)-dimensional parts, the gluons (incoming and outgoing) with their reference vectors, the indices (Lorentz indices D-dimensional, colour indices in the fundamental (NR) or adjoint (NA) representation), the auxiliary vectors and their definitions, and the transversality conditions eps(l).l = eps(l).r = 0 on both sides of the product.  For the gluon x whose momentum was eliminated by momentum conservation, eps(x).x = 0 is imposed in its eliminated form, solved for the last momentum that is not x's reference vector (as in the FeynCalc route).*)


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
  Require[FreeQ[gluons, pa], "The polarized parton a is a gluon in " <> name, gluons];
  usedIndices = Union[usedIndices, Flatten[{Symbol["mu" <> ToString[#]], Symbol["nu" <> ToString[#]]} & /@ gluons]];
  lorentz = Select[usedIndices, StringMatchQ[SymbolName[#], ("Lor" | "cLor" | "mu" | "nu") ~~ ___] &];
  colour = Complement[usedIndices, lorentz];
  StringJoin["* process.h: ", name, " (written by 2-Algebra-FORM.wl)\n",
    "#define KINEMATICS \"", kinematics, "\"\n",
    "#define HATTED \"", StringRiffle[SymbolName /@ hatted, ","], "\"\n",
    "* gluons: incoming, outgoing, and the reference vectors of their polarization sums\n",
    "#define INGLUONS \"", StringRiffle[SymbolName /@ Intersection[gluons, {pa, pb}], ","], "\"\n",
    "#define OUTGLUONS \"", StringRiffle[SymbolName /@ Complement[gluons, {pa, pb}], ","], "\"\n",
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
(*Each process (the Born process, the one-loop interference, every real-emission process, every channel of the Born table, and the process of the gauge check) gets a directory Results/FORM/class-<n>/work/<process>/ with process.h, the diagrams amp<i>.inc and the conjugate diagrams bar<j>.inc.  The external gluons of a process are read off its records (the polarization vectors), so that the Born diagrams of the table channels need no particle list; for the processes of the class the list is checked against classes.wl.  The values of the denominator symbols den<n> = 1/q^2 and of the factors invlr<l> = 1/(l.r) of the polarization sums are kept (symbolValues) for the reading of the results, and the kinematics of every process (processKinematics) for the normal reduction of its squares.*)


(* ::Input:: *)
classDirectory = ResultsDirectory["FORM/class-" <> class];
workDirectory[name_] := With[{dir = FileNameJoin[{classDirectory, "work", name}]},
  If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]]; dir];
momentaOf[particles_] := {pa, pb, p, k1, k2}[[1 ;; Length[particles]]];
gluonsOf[particles_] := Pick[momentaOf[particles], particles, "g"];
gluonsIn[recs_] := Union[Cases[recs, PolVector[l_, __] :> l, Infinity]];                        (* the external gluons of a list of diagrams (records) *)
prepareDiagram[recs_, rule_, shift_] := stripPolarizations[applyMomenta[recordParts[#], rule, shift]] & /@ recs;
symbolValues = <||>; processKinematics = <||>;
(* write one process: amp<i>.inc for the left-hand diagrams, bar<j>.inc for the conjugates of the right-hand ones *)
writeProcess[name_, kinematics_, references_, rule_, gluons_List, left_List, right_List] := Module[{dir = workDirectory[name], table},
  newRegistry[];
  Do[Export[FileNameJoin[{dir, "amp" <> ToString[i] <> ".inc"}], diagramText[left[[i]], "A", name <> ": diagram " <> ToString[i]], "Text"], {i, Length[left]}];
  Do[Export[FileNameJoin[{dir, "bar" <> ToString[j] <> ".inc"}], diagramText[conjugateRecord /@ right[[j]], "B", name <> ": conjugate of diagram " <> ToString[j]], "Text"], {j, Length[right]}];
  Export[FileNameJoin[{dir, "process.h"}], processHeader[name, kinematics, references, rule, gluons, If[kinematics === "Real", {k1, k2}, If[kinematics === "Virtual", {loop}, {}]]], "Text"];
  table = If[kinematics === "Real", realProducts, bornProducts];
  symbolValues[name] = Join[KeyValueMap[#2 -> 1/dot[#1, #1, table] &, KeySelect[denominatorSymbols, FreeQ[#, loop] &]],
    Table[Symbol["invlr" <> SymbolName[l]] -> 1/dot[l, references[l], table], {l, gluons}]];
  processKinematics[name] = kinematics;
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
(*Before the squares: checks.frm verifies the forms of the HVBM gamma5, of the density of a and of the projection of c used in common.h: Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = -4 i e_(a,b,c,d) (the normalization, the same as in the FeynCalc route), q-slash gamma5 = -i/6 e_(q,a,b,c) gamma^a gamma^b gamma^c inside a trace with D-dimensional momenta, and the density pa-slash gamma5 sa-slash together with the projection p-slash gamma5 sc-slash in this short form against the four-matrix gamma5, in one trace and in two separate traces; FORM must print Gamma5Checks = 0.*)


(* ::Input:: *)
checkRun = RunProcess[{ToolPath["FORM"], "-p", StringRiffle[{formDirectory, DirectoryName[ToolPath["FORMColor"]]}, ":"], FileNameJoin[{formDirectory, "checks.frm"}]}];
Require[checkRun["ExitCode"] === 0 && StringContainsQ[checkRun["StandardOutput"], "Gamma5Checks = 0;"], "The gamma5 checks of checks.frm failed", checkRun["StandardOutput"]];
"gamma5 checks passed"


(* ::Text:: *)
(*The square of a process for one component: the jobs of all pairs i <= j, and their sum (diagonal products as they are, i < j twice the real part), with the values of den<n> and invlr<l> put in.  A diagonal product |M_i|^2 is real; its imaginary part must vanish exactly at a rational point (as in the FeynCalc route; the point includes the normal components nk1 and nloop).  In the real kinematics the normal component is then reduced as in the FeynCalc route (normalReduction, applied there to every product): the square is a polynomial in nk1 of degree at most two, the odd part integrates to zero (it is odd under the reflection of the scattering plane), and nk1^2 = k1par^2 + spe3 with k1par^2 = d.G^-1.d, d = (pa.k1, pb.k1, p.k1) = (-A/2, -B/2, C/2), from k1 = k1par + k1n + k1hat and k1^2 = 0.*)


(* ::Input:: *)
checkPoint = {s -> 1, t -> -3/10, u -> -2/5, A -> -1/7, B -> -2/9, C -> -3/10 + 1/7 + 2/9, D -> 37/10, Nc -> 3, spe3 -> 1/11, nk1 -> 3/13, nloop -> 2/13, gs -> 1, Nf -> 5};
k1Parallel2 = parallelSquare[{-A/2, -B/2, C/2}, gramMatrix[realProducts]];                          (* k1par^2 = d.G^-1.d *)
normalReduction[expr_] := Module[{coefficients = CoefficientList[Expand[expr], nk1]},
  Require[Length[coefficients] <= 3, "The real-emission square has a power of n.k1 above two", Length[coefficients] - 1];
  If[coefficients === {}, 0,     (* a product that vanishes identically (CoefficientList[0, nk1] is the empty list) *)
    coefficients[[1]] + If[Length[coefficients] == 3, coefficients[[3]] (k1Parallel2 + spe3), 0]]];
squareJobs[name_, n_, spinCase_] := Flatten[Table[<|"Name" -> name <> "-" <> spinCase <> "-" <> ToString[i] <> "-" <> ToString[j],
  "Process" -> name, "I" -> i, "J" -> j, "Spin" -> spinCase|>, {i, n}, {j, i, n}]];
square[name_, n_, spinCase_] := Module[{jobs = squareJobs[name, n, spinCase], results, values = symbolValues[name], sum},
  results = runJobs[jobs];
  Require[And @@ MapThread[#2["I"] != #2["J"] || (#1[[2]] /. values /. checkPoint) === 0 &, {results, jobs}], "A diagonal product has a nonzero imaginary part", name];
  sum = Total[MapThread[If[#2["I"] == #2["J"], 1, 2] #1[[1]] &, {results, jobs}]] /. values;
  If[processKinematics[name] === "Real", normalReduction[sum], sum]];
"parallel FORM jobs: " <> ToString[nWorkers]


(* ::Section:: *)
(*The transfer Born table*)


(* ::Text:: *)
(*The D-dimensional transfer squares |M|^2/gs^4 of the 2 -> 2 channels a^ b -> c^ of classes.wl (first parton polarized, observed parton c, unobserved parton k1, u = -s-t), from the Born records that step 1 wrote for them, with the reference vectors pa -> pb, pb -> pa, p -> pa, k1 -> pa for the gluons.  All four components are formed; XN and NX must vanish exactly.  Keys as in the FeynCalc route ("q qp -> q"); the table holds the components XX and NN.*)


(* ::Input:: *)
channelPartons[name_] := StringSplit[StringReplace[name, "->" -> ""]];                                (* {a, b, c} *)
tableName[name_] := "table-" <> StringReplace[name, {" " -> "_", "->" -> "to"}];
bornTableRaw = Association @ KeyValueMap[Function[{name, recs}, name -> Module[{diagrams = prepareDiagram[#, bornConservation, {}] & /@ recs, in = channelPartons[name][[1 ;; 2]]},
    writeProcess[tableName[name], "Born", bornReferences, bornConservation, gluonsIn[recs], diagrams, diagrams];
    Association @ Table[spinCase -> Factor[colourFactors[Together[average[in] normalization[spinCase] square[tableName[name], Length[diagrams], spinCase]/gs^4]] /. u -> -s - t], {spinCase, spinCases}]]],
  records["BornTable"]];
Require[FreeQ[bornTableRaw, Complex | nk1 | nloop], "A transfer Born square is not a real scalar", bornTableRaw];
Require[And @@ Flatten[Table[bornTableRaw[name][mixed] === 0, {name, Keys[bornTableRaw]}, {mixed, {"XN", "NX"}}]],
  "A mixed component XN or NX of a transfer Born square does not vanish", Map[KeyTake[#, {"XN", "NX"}] &, bornTableRaw]];
bornTable = Association @ Table[spinCase -> Map[#[spinCase] &, bornTableRaw], {spinCase, outputSpins}];
Dataset[Map[Simplify[# /. D -> 4] &, bornTable, {2}]]


(* ::Text:: *)
(*Checks of the table (those of the FeynCalc route).  (1) In four dimensions the two components coincide in every channel (the evanescent terms differ).  (2) The channels whose leading-order transverse spin transfer is in the literature: the Collins hard factors of F. Yuan, Phys. Rev. Lett. 100 (2008) 032003 [arXiv:0709.3272], Eq. (6), transcribed verbatim in reference/Yuan-2008-Collins-hard-factors.wl, written as H = alphaS^2 Pi/s^2 x h(s,t,u); the unpolarized hard factor in the same normalization is dsigma/dt = alphaS^2 Pi/s^2 x |M|^2/gs^4 (averaged), so h is compared with the table directly (no free factor).  (3) The q qb channels are not listed by Yuan; they are checked against charge conjugation (qb q -> q is the conjugate of q qb -> qb) and, when Processes/Collins has been run, against its independently written resolved-photon Born table (6-ResolvedBorn.wl, Results/ResolvedBorn.wl).*)


(* ::Input:: *)
fourDimensional[expr_] := Simplify[expr /. D -> 4 /. CA -> Nc /. CF -> (Nc^2 - 1)/(2 Nc) /. u -> -s - t];
Require[And @@ Table[fourDimensional[bornTable["XX"][name] - bornTable["NN"][name]] === 0, {name, Keys[bornTableRaw]}],
  "The two transverse components of a transfer Born square differ in four dimensions"];
Get[FileNameJoin[{NotebookDirectory[], "reference", "Yuan-2008-Collins-hard-factors.wl"}]];
literature = YuanCollinsHardFactors[s, t, u, Nc];                         (* h(s,t,u) for the channels of Eq. (6) *)
literatureDifferences = Association @ KeyValueMap[#1 -> fourDimensional[bornTable["XX"][#1] - #2] &, literature];
Require[And @@ (# === 0 & /@ Values[literatureDifferences]), "A transfer Born square differs from the literature", literatureDifferences];
Require[fourDimensional[bornTable["XX"]["qb q -> q"] - bornTable["XX"]["q qb -> qb"]] === 0 && Simplify[bornTable["NN"]["qb q -> q"] - bornTable["NN"]["q qb -> qb"]] === 0,
  "The charge-conjugate channels qb q -> q and q qb -> qb differ"];
collinsFile = FileNameJoin[{NotebookDirectory[], "..", "Collins", "Results", "ResolvedBorn.wl"}];
collinsDifferences = If[FileExistsQ[collinsFile], With[{resolved = Get[collinsFile]["Transfer"]},
    Association @ Table[name -> fourDimensional[bornTable["XX"][name] - resolved[name]["Collins"]], {name, Intersection[Keys[bornTableRaw], Keys[resolved]]}]],
  "Processes/Collins/Results/ResolvedBorn.wl not found: comparison skipped"];
Require[!AssociationQ[collinsDifferences] || And @@ (# === 0 & /@ Values[collinsDifferences]), "A transfer Born square differs from Processes/Collins' resolved Born table", collinsDifferences];
{literatureDifferences, collinsDifferences}


(* ::Section:: *)
(*Born*)


(* ::Text:: *)
(*The Born square of the class with the reference vectors pa -> pb, pb -> pa, p -> pa, k1 -> pa, averaged over the initial colours (with the 1/2 of the density) and over the spin of b, normalized, with u = -s - t, for the four components; XN = NX = 0 exactly; it must equal the table's entry for the class.*)


(* ::Input:: *)
bornDiagrams = prepareDiagram[#, bornConservation, {}] & /@ records["Born"];
If[bornDiagrams =!= {},
  Require[Sort[gluonsIn[records["Born"]]] === Sort[gluonsOf[spec["Born"]]], "The gluons of the Born records differ from classes.wl", {gluonsIn[records["Born"]], spec["Born"]}];
  writeProcess["born", "Born", bornReferences, bornConservation, gluonsOf[spec["Born"]], bornDiagrams, bornDiagrams]];
bornAll = If[bornDiagrams === {}, Association @ Table[spinCase -> 0, {spinCase, spinCases}],
  Association @ Table[spinCase -> Factor[colourFactors[Together[average[incoming] normalization[spinCase] square["born", Length[bornDiagrams], spinCase]]] /. u -> -s - t], {spinCase, spinCases}]];
Require[FreeQ[bornAll, Complex | nk1 | nloop], "The Born square is not a real scalar", bornAll];
Require[bornAll["XN"] === 0 && bornAll["NX"] === 0, "A mixed component XN or NX of the Born square does not vanish", bornAll];
born = KeyTake[bornAll, outputSpins];
If[spec["Born"] =!= None,
  Require[And @@ Table[Simplify[born[sc]/gs^4 - bornTable[sc][StringRiffle[spec["Born"][[1 ;; 2]], " "] <> " -> " <> spec["Born"][[3]]]] === 0, {sc, outputSpins}],
    "The Born square differs from the transfer Born table", born]];
born


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
If[virtualGraphs =!= {},
  Require[Sort[gluonsIn[records["Virtual"]]] === Sort[gluonsOf[spec["Born"]]], "The gluons of the one-loop records differ from classes.wl"];
  writeProcess["virtual", "Virtual", bornReferences, bornConservation, gluonsOf[spec["Born"]], virtualGraphs, {Flatten[bornDiagrams, 1]}]];
Counts[routings[[All, 1]]]


(* ::Text:: *)
(*One job per graph and component: the graph times the conjugate of the summed Born amplitude, M1_graph M0*, with its real and imaginary parts.  FORM returns the loop momentum through the symbols L2 = loop.loop, Lpa = loop.pa, Lpb = loop.pb, Lp = loop.p and Lhat2 = loophat.loophat (the four-dimensional loop4.loop4 = L2 - Lhat2 in the HVBM scheme), and its normal component through nloop = n.loop (the products of the axes with the loop momentum were already expanded into Lpa, Lpb, Lp by FORM).*)
(**)
(*The normal component of the loop momentum, exactly as in the FeynCalc route.  The integrand depends on the loop momentum through loop^2, loop.pa, loop.pb, loop.p and through its components orthogonal to the scattering plane: n.loop and the hatted part.  The orthogonal space has D - 3 dimensions (the direction n and the D - 4 extra dimensions, all spacelike), and the integration over its orientations is symmetric, so that with m = D - 3 and the orthogonal square loopPerp^2 = loop^2 - loopPar^2, loopPar^2 = d.G^-1.d, d = (Lpa, Lpb, Lp),*)
(**)
(*   (n.loop)^(2a) Lhat2^b  ->  (-1)^a (loopPerp^2)^(a+b) (1/2)_a ((m-1)/2)_b / (m/2)_(a+b)     (odd powers of n.loop vanish),*)
(**)
(*the moments of a direction on the unit sphere of m dimensions: for a = 1, b = 0 this is -loopPerp^2/(D-3), for a = 1, b = 1 it is -(D-4) (loopPerp^2)^2/((D-3)(D-1)).  Terms without n.loop keep their Lhat2 power as HatPower[b] for the dimension shift of step 3, as in Processes/Hadronic.  Then the four products are written through the inverse propagators D_j of the graph's family (the four equations loop^2 = D1, 2 loop.q_j = D_j - D_1 - q_j^2 are solved for L2, Lpa, Lpb, Lp), and every term is read off as coefficient x family[indices] x HatPower[n].*)
(**)
(*Normalization, as in the FeynCalc route: FeynArts writes the one-loop amplitude with the measure d^D loop/(2 Pi)^D; the scalar integrals are defined with the measure d^D loop/(i Pi^(D/2)), so that M1 = i/(4 Pi)^(D/2) x (sum of coefficients times integrals).  The products (coefficient x M0*) are real, and 2 Re(M1 M0*) = (2/(4 Pi)^(D/2)) Sum_j [i c_j M0*] Re(I_j), with the real parts of the integrals in the physical region (s > 0) taken in step 3.  The association virtual below holds Sum_j [i c_j M0*] I_j for the two components XX and NN (i times the complex product of FORM, which must come out real); the factor 2/(4 Pi)^(D/2) is included in step 4.  The mixed components are formed as well and must vanish.*)


(* ::Input:: *)
inverseProps[name_] := Module[{q = families[name], sol},
  sol = Solve[Table[Dj[j] == L2 + 2 (Coefficient[q[[j]], pa] Lpa + Coefficient[q[[j]], pb] Lpb + Coefficient[q[[j]], p] Lp) + dot[q[[j]] - loop, q[[j]] - loop, bornProducts] /. u -> -s - t, {j, 4}], {L2, Lpa, Lpb, Lp}];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators"]; First[sol]];
loopRules = Association @ Table[name -> inverseProps[name], {name, Keys[families]}];
loopPerp2 = L2 - parallelSquare[{Lpa, Lpb, Lp}, gramMatrix[bornProducts] /. u -> -s - t];           (* the orthogonal square of the loop momentum *)
perpMoment[a_, b_] := (-1)^a Pochhammer[1/2, a] Pochhammer[(D - 4)/2, b]/Pochhammer[(D - 3)/2, a + b];      (* m = D - 3 *)
normalLoopReduction[expr_] := Module[{terms = DeleteCases[termsOf[Expand[expr]], 0]},
  Total[Table[With[{a = Exponent[term, nloop], b = Exponent[term, Lhat2]},
    Which[OddQ[a], 0, a == 0, term, True, (term/(nloop^a Lhat2^b)) perpMoment[a/2, b] loopPerp2^(a/2 + b)]], {term, terms}]]];
(* the interference of one graph -> sum of coefficient x family[indices] x HatPower[n] *)
toFamilies[value_, name_] := Module[{expr},
  expr = normalLoopReduction[Expand[value /. symbolValues["virtual"]]] /. loopRules[name];
  Require[FreeQ[expr, nloop], "The normal component of the loop momentum was not reduced"];
  expr = Expand[expr /. Table[Symbol["dj" <> ToString[j]] -> 1/Dj[j], {j, 4}]];
  Require[FreeQ[expr, L2 | Lpa | Lpb | Lp], "A scalar product of the loop momentum was left"];
  Total[Table[With[{powers = Table[Exponent[term, Dj[j]], {j, 4}], hat = Exponent[term, Lhat2]},
    (term/(Times @@ Table[Dj[j]^powers[[j]], {j, 4}] Lhat2^hat)) familyHeads[name] @@ (-powers) HatPower[hat]], {term, DeleteCases[termsOf[expr], 0]}]]];
virtualJobs[spinCase_] := Table[<|"Name" -> "virtual-" <> spinCase <> "-" <> ToString[g], "Process" -> "virtual", "I" -> g, "J" -> 1, "Spin" -> spinCase|>, {g, Length[virtualGraphs]}];
virtualAll = If[virtualGraphs === {}, Association @ Table[spinCase -> 0, {spinCase, spinCases}],
  Association @ Table[spinCase -> Module[{results = runJobs[virtualJobs[spinCase]]},
    Collect[I average[incoming] normalization[spinCase] Total[MapThread[toFamilies[#1[[1]] + I #1[[2]], #2[[1]]] &, {results, routings}]] /. u -> -s - t,
      {_KVirtualST, _KVirtualSU, _KVirtualTU, _HatPower}, colourFactors]], {spinCase, spinCases}]];
Require[FreeQ[virtualAll, Complex | nloop], "The virtual interference is not real or still contains n.loop", Short[virtualAll, 5]];
Require[virtualAll["XN"] === 0 && virtualAll["NX"] === 0, "A mixed component XN or NX of the virtual interference does not vanish", KeyTake[virtualAll, {"XN", "NX"}]];
virtual = KeyTake[virtualAll, outputSpins];
{Length[virtual["XX"]], Length[virtual["NN"]], Union[Cases[virtual, (KVirtualST | KVirtualSU | KVirtualTU)[__], Infinity]] // Length, Union[Cases[virtual, _HatPower, Infinity]]}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*For every real-emission process of the class, the square of the sum of its diagrams for the four components, in terms of the invariants s, t, u, A, B, C and the hatted invariant spe3 (the normal component nk1 reduced in every square); the processes are added with their weights (classes.wl: 1/2 for two identical unobserved partons, Nf - 1 for the unobserved quark pairs of other flavours).  Reference vectors of the physical polarization sums: pa for the gluons p and k1, pb for k2 (and the other incoming momentum for an incoming gluon).  As in the FeynCalc route the sum is not brought to a common denominator (step 3 expands it into partial fractions anyway); the colour factors are polynomials in Nc = CA.  The mixed components must vanish; they are checked at exact rational points (the symbolic sum of thousands of terms is not brought to a common denominator).  In the HVBM scheme the squares are at most linear in spe3.*)


(* ::Input:: *)
processName[key_] := "real-" <> StringReplace[key, " " -> "_"];
realDiagrams = Association @ KeyValueMap[#1 -> (prepareDiagram[#, realConservation, {}] & /@ #2) &, records["Real"]];
KeyValueMap[(Require[Sort[gluonsIn[records["Real"][#1]]] === Sort[gluonsOf[StringSplit[#1]]], "The gluons of the real-emission records differ from classes.wl", #1];
  writeProcess[processName[#1], "Real", realReferences, realConservation, gluonsOf[StringSplit[#1]], #2, #2]) &, realDiagrams];
realSquares = Association @ Table[spinCase -> Association @ KeyValueMap[#1 -> normalization[spinCase] square[processName[#1], Length[#2], spinCase] &, realDiagrams], {spinCase, spinCases}];
realAll = Association @ Table[spinCase -> (average[incoming] Total[KeyValueMap[weights[#1] #2 &, realSquares[spinCase]]] /. Nc -> CA), {spinCase, spinCases}];
Require[FreeQ[realAll, Complex], "The real-emission square is not real"];
Require[FreeQ[realAll, nk1 | nloop], "The normal component was left in the real-emission square"];
mixedPoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], CA -> 2 + 10 v[[6]], spe3 -> v[[7]]/5, Nf -> 5}], {3}];
mixedValues = Table[Together[realAll[mixed] /. C -> -(s + t + u) - A - B /. point], {mixed, {"XN", "NX"}}, {point, mixedPoints}];
Require[mixedValues === {{0, 0, 0}, {0, 0, 0}}, "A mixed component XN or NX of the real-emission square does not vanish", mixedValues];
real = KeyTake[realAll, outputSpins];
Require[Max[Exponent[#, spe3] & /@ Join[termsOf[Expand[real["XX"]]], termsOf[Expand[real["NN"]]]]] <= 1, "The real-emission square has spe3^2 terms"];
LeafCount /@ real


(* ::Text:: *)
(*Check of the gauge invariance of the real-emission square: for the first process of the class and the XX component, the physical polarization sums with other reference vectors (pb for p and k1, pa for k2) give the same result, exactly, at three random rational points (all invariants, D and Nc, with C = -s4 - A - B).  The diagrams are written again (process "gauge"), because the transversality conditions eps(l).r = 0 depend on the reference vectors.*)


(* ::Input:: *)
alternativeReferences = <|pa -> pb, pb -> pa, p -> pb, k1 -> pb, k2 -> pa|>;
firstProcess = First[Keys[realDiagrams]];
writeProcess["gauge", "Real", alternativeReferences, realConservation, gluonsOf[StringSplit[firstProcess]], realDiagrams[firstProcess], realDiagrams[firstProcess]];
alternative = normalization["XX"] square["gauge", Length[realDiagrams[firstProcess]], "XX"];
gaugePoints = Table[With[{v = RandomInteger[{1, 97}, 7]/RandomInteger[{98, 197}, 7]},
  {s -> 1, t -> -v[[1]], u -> -v[[2]] (1 - v[[1]]), A -> -v[[3]], B -> -v[[4]], D -> 4 - v[[5]], Nc -> 2 + 10 v[[6]], spe3 -> v[[7]]/5}], {3}];
referenceDifference = Table[Together[(alternative - realSquares["XX"][firstProcess]) /. C -> -(s + t + u) - A - B /. point], {point, gaugePoints}];
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
(*The result file Results/FORM/class-<n>/Algebra.wl, in the format of the FeynCalc route: the Born square, the virtual interference and the real-emission square in the components XX and NN, and the transfer Born table for steps 3 and 4.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Virtual" -> virtual, "Real" -> real, "BornTable" -> bornTable, "Class" -> class, "Spec" -> spec, "Spins" -> outputSpins,
  "Families" -> families, "Variables" -> "s, t, u = -s-t; A = (pa-k1)^2, B = (pb-k1)^2, C = (p+k1)^2 = -s4-A-B; spe3; loopHat2 via HatPower[n]; XX divided by |eX||eC| = -t u/(t+u)"|>;
SaveResult["Algebra", algebra, "FORM/class-" <> class];
Keys[algebra]


(* ::Section:: *)
(*Comparison with the FeynCalc route*)


(* ::Text:: *)
(*If the FeynCalc route has been run for this class, the two files must agree (otherwise the comparison is skipped and the cell says so): Born, Virtual and Real of both components and every entry of the Born table (both components) are compared at two exact random rational points of s, t, u, A, B, D, CA, spe3, Nf, gs (with C = -(s+t+u) - A - B and CF = (CA^2 - 1)/(2 CA)); in the virtual every scalar integral is replaced by its own, distinct random rational number and HatPower[n] by hp^n, so that the coefficient of every integral is compared.  The arithmetic is exact, but the comparison is sampled at two points, not a symbolic identity.*)
(**)
(*The one-loop interference is compared in two ways.  Literally, coefficient by coefficient of every integral and hat power.  And after a projection that removes the one freedom of representation of the loop integrands: the normal component and the hatted square of the loop momentum are not independent of its other products in four dimensions, (n.loop)^2 = loopPar^2 - loop4.loop4 = loopPar^2 - L2 + Lhat2 (the Gram determinant of the five four-dimensional vectors pa, pb, p, n, loop4 vanishes), and the two routes may contract the two epsilon tensors into determinants that differ by a multiple of this identity (they do in the NN component of class 13, the only transfer class whose loop graphs contain an incoming gluon; the difference is exactly such a multiple graph by graph).  Such a multiple integrates to zero: (n.loop)^2 -> -loopPerp^2/(D-3) and Lhat2 -> (D-4)/(D-3) loopPerp^2 (the exact averages over the (D-3)-dimensional orthogonal space used in the Virtual section, perpMoment[1, 0] and perpMoment[0, 1]; step 3 evaluates HatPower[1] by the equivalent dimension shift), so that n^2 + perp^2 - Lhat2 -> loopPerp^2 (-1/(D-3) + 1 - (D-4)/(D-3)) = 0.  The projection replaces every HatPower[b] by its average perpMoment[0, b] (loopPerp^2)^b HatPower[0], with loopPerp^2 = L2 - loopPar^2 written through the inverse propagators of the family (a power of D_j in the numerator lowers the j-th index); both files are projected the same way.  Equality after the projection is required; the literal comparison is shown.*)


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
    (* the hat projection of a virtual interference, at the kinematic part of a point (the integrals stay symbols until the shifts are made) *)
    projectHat[x_, point_] := Module[{kin = Select[point, FreeQ[#[[1]], KVirtualST | KVirtualSU | KVirtualTU] &], expr},
      expr = Expand[x /. CF -> (CA^2 - 1)/(2 CA) /. kin];
      Expand[expr /. HatPower[b_] h_[n__] /; b > 0 && MemberQ[{KVirtualST, KVirtualSU, KVirtualTU}, h] :>
        With[{name = StringDrop[SymbolName[h], 1]}, (perpMoment[0, b] /. kin) Total[Map[(#[[2]] HatPower[0] h @@ ({n} - #[[1]])) &,
          CoefficientRules[Expand[(loopPerp2 /. loopRules[name] /. kin)^b], {Dj[1], Dj[2], Dj[3], Dj[4]}]]]]]];
    evaluateProjected[x_, point_] := Together[projectHat[x, point] /. HatPower[0] -> 1 /. point];
    Require[Keys[other["BornTable"]["XX"]] === Keys[algebra["BornTable"]["XX"]], "The two routes have different Born tables"];
    parts = Join[Flatten[Table[{part, spinCase} -> {algebra[part][spinCase], other[part][spinCase]}, {part, {"Born", "Virtual", "Real"}}, {spinCase, outputSpins}]],
      Flatten[Table[{"BornTable", spinCase, name} -> {algebra["BornTable"][spinCase][name], other["BornTable"][spinCase][name]}, {spinCase, outputSpins}, {name, Keys[algebra["BornTable"]["XX"]]}]]];
    Join[Association @ Table[item[[1]] -> (And @@ Table[If[item[[1, 1]] === "Virtual", evaluateProjected, evaluate][item[[2, 1]] - item[[2, 2]], point] === 0, {point, points}]), {item, parts}],
      <|"VirtualLiteral" -> Association @ Table[spinCase -> (And @@ Table[evaluate[algebra["Virtual"][spinCase] - other["Virtual"][spinCase], point] === 0, {point, points}]), {spinCase, outputSpins}]|>]]];
Require[StringQ[routeComparison] || And @@ Values[KeyDrop[routeComparison, "VirtualLiteral"]], "The FORM route differs from the FeynCalc route", Select[KeyDrop[routeComparison, "VirtualLiteral"], ! # &]];
If[StringQ[routeComparison], routeComparison, {Counts[Values[KeyDrop[routeComparison, "VirtualLiteral"]]], "literal virtual comparison" -> routeComparison["VirtualLiteral"]}]
