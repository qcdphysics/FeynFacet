(* ::Title:: *)
(*Lepton-quark scattering with transverse spin transfer, step 2 (Tracer route): spin and colour algebra*)


(* ::Text:: *)
(*Same task and same output format as the FeynCalc route, with the classic Tracer package for the Dirac algebra; the construction is that of Processes/LeptonQuark/2-Algebra-Tracer.wl (read its introduction for Tracer's vocabulary and conventions: Tracer's gamma5 and Eps have the opposite sign to FeynCalc's, which does not matter here because every transverse tensor contains two gamma5's).  Run this file in a fresh kernel: Tracer and FeynCalc cannot be loaded together.*)
(**)
(*The transverse spin transfer T replaces the incoming quark density k-slash/2 by k-slash gamma5 a-slash/2 and the outgoing spin sum p-slash by p-slash gamma5 b-slash, with the spin vectors of the FeynCalc route: the normal n (n.l = n.k = n.p = 0, n^2 = -1, n.r = nr) for T_NN, and the in-plane axes aX = (t/s) l + (u/s) k + p, bX = t/(t+u) l - u/(t+u) k - s(u-t)/(t+u)^2 p for T_XX, divided by the product of their lengths -t u/(t+u).  In Tracer the spin vectors are vectors of their own whose scalar products are declared with OnShell; nr^2 is reduced as in the FeynCalc route (nr^2 = rpar^2 + spe3) and odd powers of nr are dropped.  All vectors except lp and r are four-dimensional (H = 0).*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadTracer[]


(* ::Input:: *)
records = LoadAmplitudes[]["Records"];
VectorDimension[D];
AntiCommute[off];
Spur[1]; Spur[2];
ShowRecords[Flatten[records["Real"]]]


(* ::Section:: *)
(*From records to traces*)


(* ::Text:: *)
(*conjugate and the colour rules as in the lepton-quark process.  sewTrace builds the string of one fermion line for a spin case: "Unpolarized" inserts pOut, pIn; "NN" and "XX" insert pOut, G5, sOut and pIn, G5, sIn on the quark line (line 2; the lepton line is always unpolarized).*)


(* ::Input:: *)
conjugate[record_] := Module[{indices, fresh},
  indices = Union[Cases[record, i_Symbol /; StringStartsQ[SymbolName[i], "Lor"] :> i, Infinity]];
  fresh = Thread[indices -> (Symbol["c" <> SymbolName[#]] & /@ indices)];
  <|"Coefficient" -> (record["Coefficient"] /. fresh /. {Complex[x_, y_] :> Complex[x, -y],
      ColorT[a_, i_, j_] :> ColorT[Reverse[a], j, i], PolVector[p_, i_, sign_] :> PolVector[p, i, -sign]}),
    "Chains" -> (Reverse /@ record["Chains"] /. fresh)|>];
outgoing[chain_Chain] := First[Cases[chain, Ubar[p_] :> p]];
incoming[chain_Chain] := First[Cases[chain, Uspinor[p_] :> p]];
spinVectors = <|"NN" -> {nn, nn}, "XX" -> {aX, bX}|>;
quarkLineQ[chain_Chain] := incoming[chain] === k;
sewTrace[line_, chainL_, chainR_, spin_] := Inactive[G] @@ Join[{line, outgoing[chainL]},
  If[spin =!= "Unpolarized" && quarkLineQ[chainL], {G5, spinVectors[spin][[2]]}, {}], TracerString[chainL],
  {incoming[chainL]}, If[spin =!= "Unpolarized" && quarkLineQ[chainL], {G5, spinVectors[spin][[1]]}, {}], TracerString[chainR]];
colourSum = {ColorDelta[i_, i_] :> Nc, ColorDelta[i_, j_] /; !OrderedQ[{i, j}] :> ColorDelta[j, i],
  ColorT[{a_}, i_, x_] ColorT[{a_}, x_, j_] :> CF ColorDelta[i, j],
  ColorDelta[i_, j_]^2 :> Nc, ColorDelta[i_, j_] ColorDelta[j_, i_] :> Nc};
scalars = {Den[q_] :> 1/S[q, q], Metric[i_, j_] :> S[{i}, {j}]};
interference[termL_, termR_, spin_, bosonSums_] := Module[{conj = conjugate[termR], coefficient, traces},
  coefficient = termL["Coefficient"] conj["Coefficient"] //. colourSum //. bosonSums /. scalars;
  traces = Times @@ MapThread[sewTrace[#1, #2, #3, spin] &, {Range[Length[termL["Chains"]]], termL["Chains"], conj["Chains"]}];
  coefficient traces];
square[amplitudes_, spin_, bosonSums_, average_] := Module[{terms = Flatten[amplitudes]},
  average Sum[interference[terms[[i]], terms[[j]], spin, bosonSums], {i, Length[terms]}, {j, Length[terms]}]];
(* four-dimensional vectors have no hatted components; H is also made linear in its arguments (Tracer leaves H of a momentum sum unexpanded) *)
fourDimensional[momenta_List] := (Unprotect[H];
  H[x_, y_Plus] := Plus @@ (H[x, #] & /@ List @@ y); H[x_Plus, y_] := Plus @@ (H[#, y] & /@ List @@ x);
  H[x_, c_?NumericQ y_] := c H[x, y]; H[c_?NumericQ x_, y_] := c H[x, y];
  Scan[(H[#] = 0; H[#, _] = 0; H[_, #] = 0) &, momenta]; Protect[H]);
lengths[spin_, tIn_, uIn_] := If[spin === "XX", -tIn uIn/(tIn + uIn), 1];
(* the scalar products of the spin vectors with the momenta of a process, from the definitions of the axes *)
axisProducts[lIn_, kIn_, pOut_, sIn_, tIn_, uIn_, momenta_List, dots_Association] := Module[{base = {lIn, kIn, pOut}, ca, cb, dot},
  ca = {tIn/sIn, uIn/sIn, 1}; cb = {tIn/(tIn + uIn), -uIn/(tIn + uIn), -sIn (uIn - tIn)/(tIn + uIn)^2};
  dot[x_, y_] := Lookup[dots, Key[Sort[{x, y}]], 0];
  Join[Table[{aX, x, Factor[Sum[ca[[i]] dot[base[[i]], x], {i, 3}]]}, {x, momenta}], Table[{bX, x, Factor[Sum[cb[[i]] dot[base[[i]], x], {i, 3}]]}, {x, momenta}],
    {{aX, aX, Factor[Sum[ca[[i]] ca[[j]] dot[base[[i]], base[[j]]], {i, 3}, {j, 3}]]}, {bX, bX, Factor[Sum[cb[[i]] cb[[j]] dot[base[[i]], base[[j]]], {i, 3}, {j, 3}]]},
     {aX, bX, Factor[Sum[ca[[i]] cb[[j]] dot[base[[i]], base[[j]]], {i, 3}, {j, 3}]]}, {aX, nn, 0}, {bX, nn, 0}, {nn, nn, -1}}, Table[{nn, x, 0}, {x, base}]]];


(* ::Section:: *)
(*Born: l q -> l q*)


(* ::Input:: *)
bornDots = <|{k, l} -> s/2, {lp, p} -> s/2, {l, lp} -> -t/2, {k, p} -> -t/2, {l, p} -> -u/2, {k, lp} -> -u/2|>;
bornProducts = axisProducts[l, k, p, s, t, u, {l, k, p, lp}, bornDots];
born = Association @ Table[spin -> Internal`InheritedBlock[{S, H},
  fourDimensional[{l, k, lp, p, nn, aX, bX}];
  OnShell @@ Join[{on, {l, 0}, {k, 0}, {lp, 0}, {p, 0}}, KeyValueMap[{#1[[1]], #1[[2]], #2} &, bornDots], bornProducts, {{nn, lp, 0}}];
  Factor[Expand[Activate[square[records["Born"], spin, {}, 1/(4 Nc)]]]/lengths[spin, t, u]]], {spin, {"Unpolarized", "NN", "XX"}}];
born["Transverse"] = Factor[(born["NN"] + born["XX"])/2];
Require[Simplify[(born["Unpolarized"] - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];
Require[Simplify[{born["NN"] - ee^4 eq^2 (-4 s u/t^2 - (D - 4)), born["XX"] - ee^4 eq^2 (-4 s u/t^2 + (D - 4))} /. u -> -s - t] === {0, 0},
  "Transverse Born tensors differ from Eq. (transverse-born-tensor) of the paper"];
born


(* ::Section:: *)
(*Real emission: l q -> l q g*)


(* ::Input:: *)
realDots = <|{k, l} -> s/2, {k, p} -> -t/2, {l, p} -> -u/2, {k, r} -> -a/2, {p, r} -> b/2, {l, lp} -> -(a + b + t)/2, {l, r} -> (s + t + u + a + b)/2,
  {k, lp} -> (s + t + a)/2, {lp, p} -> -(t + u + b)/2, {lp, r} -> (s + t + u)/2|>;
realProducts = axisProducts[l, k, p, s, t, u, {l, k, p, lp, r}, realDots];
gluonSum = {PolVector[r, i_, -1] PolVector[r, j_, 1] :> -S[{i}, {j}]};
gram = Table[Lookup[realDots, Key[Sort[{x, y}]], 0], {x, {l, k, p}}, {y, {l, k, p}}];
parallelSquare = Factor[{realDots[{l, r}], realDots[{k, r}], realDots[{p, r}]}.Inverse[gram].{realDots[{l, r}], realDots[{k, r}], realDots[{p, r}]}];
normalReduction[expr_] := Module[{coefficients = CoefficientList[Expand[expr], nr]},
  Require[Length[coefficients] <= 3, "The real-emission tensor has a power of nr above two", Length[coefficients] - 1];
  Factor[coefficients[[1]] + If[Length[coefficients] == 3, coefficients[[3]] (parallelSquare + spe3), 0]]];
real = Association @ Table[spin -> Internal`InheritedBlock[{S, H},
  fourDimensional[{l, k, p, nn, aX, bX}];
  Unprotect[H]; H[r, r] = spe3; H[lp, lp] = spe3; H[lp, r] = -spe3; Protect[H];
  OnShell @@ Join[{on, {l, 0}, {k, 0}, {lp, 0}, {p, 0}, {r, 0}}, KeyValueMap[{#1[[1]], #1[[2]], #2} &, realDots], realProducts, {{nn, r, nr}, {nn, lp, -nr}}];
  normalReduction[Expand[Activate[square[records["Real"], spin, gluonSum, 1/(4 Nc)]]]/lengths[spin, t, u]]], {spin, {"Unpolarized", "NN", "XX"}}];
real["Transverse"] = Factor[(real["NN"] + real["XX"])/2];
Require[Factor[(real["NN"] - real["XX"]) /. {D -> 4, spe3 -> 0}] === 0, "The two transverse components differ in four dimensions"];
Factor[real["Transverse"] /. {D -> 4, spe3 -> 0}]


(* ::Section:: *)
(*Photon-quark Born: gamma q -> q g*)


(* ::Input:: *)
photonDots = <|{k, q} -> sh/2, {k, p} -> -th/2, {p, q} -> -uh/2, {p, r} -> sh/2, {k, r} -> -uh/2, {q, r} -> -th/2|>;
photonProducts = axisProducts[q, k, p, sh, th, uh, {q, k, p, r}, photonDots];
photonSums = {PolVector[q, i_, 1] PolVector[q, j_, -1] :> -S[{i}, {j}]/(D - 2), PolVector[r, i_, -1] PolVector[r, j_, 1] :> -S[{i}, {j}]};
photon = Association @ Table[spin -> Internal`InheritedBlock[{S, H},
  fourDimensional[{q, k, p, r, nn, aX, bX}];
  OnShell @@ Join[{on, {q, 0}, {k, 0}, {p, 0}, {r, 0}}, KeyValueMap[{#1[[1]], #1[[2]], #2} &, photonDots], photonProducts, {{nn, r, 0}}];
  Factor[Expand[Activate[square[records["Photon"], spin, photonSums, 1/(2 Nc)]]]/lengths[spin, th, uh]]], {spin, {"Unpolarized", "NN", "XX"}}];
photon["Transverse"] = Factor[(photon["NN"] + photon["XX"])/2];
Require[Simplify[(photon["Transverse"]/(CF ee^2 eq^2 gs^2) - (4 - 2 eps^2 th^2/((1 - eps) sh uh))) /. D -> 4 - 2 eps /. uh -> -sh - th] === 0,
  "The photon transverse tensor differs from Eq. (transverse-photon-born) of the paper"];
photon


(* ::Section:: *)
(*Virtual: the one-loop vertex*)


(* ::Text:: *)
(*Exactly as in the lepton-quark process: the vertex projected onto the Born current gives the scalar integrand of the form factor in the family VirtualIntegral[n0,n1,n2].*)


(* ::Input:: *)
vertex = records["VirtualVertex"][[1, 1]];
bornRecord = records["Born"][[1, 1]];
quarkChain = First[Select[vertex["Chains"], outgoing[#] === p &]];
word = TracerString[quarkChain] /. Metric -> Sequence;
gluonIndices = Cases[vertex["Coefficient"], Metric[i_, j_] /; FreeQ[{i, j}, Lor1 | Lor2] :> {i, j}][[1]];
word = word /. {gluonIndices[[2]]} -> {gluonIndices[[1]]};
currentIndex = First[Complement[Cases[word, {i_} :> i], gluonIndices]];
projector = Internal`InheritedBlock[{S, H},
  OnShell[on, {k, 0}, {p, 0}, {k, p, -t/2}];
  Factor[(G @@ Join[{1, p}, word, {k, {currentIndex}}])/G[1, p, {currentIndex}, k, {currentIndex}]]];
numerator = Expand[projector /. {S[loop, loop] -> d0, S[loop, k] -> (d0 - d1)/2, S[loop, p] -> (d0 - d2)/2}];
virtualFamily = Total[(#[[2]] VirtualIntegral @@ ({1, 1, 1} - #[[1]])) & /@ CoefficientRules[numerator, {d0, d1, d2}]];
relativeCoefficient = Simplify[(vertex["Coefficient"]/bornRecord["Coefficient"]) //. colourSum /.
  {Den[loop] -> 1, Den[k - loop] -> 1, Den[p - loop] -> 1, Den[loop - k] -> 1, Den[loop - p] -> 1, Metric[gluonIndices[[1]], gluonIndices[[2]]] -> 1}];
Require[Simplify[relativeCoefficient - (-I) CF gs^2] === 0, "Unexpected vertex coefficient", relativeCoefficient];
virtualFamily


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Photon" -> photon, "VirtualFamily" -> virtualFamily|>;
SaveResult["Algebra", algebra, "Tracer"];
Keys[algebra]
