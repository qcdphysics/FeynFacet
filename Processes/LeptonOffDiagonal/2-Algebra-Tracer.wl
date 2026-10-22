(* ::Title:: *)
(*Lepton-quark and lepton-gluon scattering with a different tagged parton, step 2 (Tracer route): spin and colour algebra*)


(* ::Text:: *)
(*Same task and same output format as the FeynCalc route, now with the classic Tracer package for the Dirac algebra; the construction is that of Processes/LeptonQuark/2-Algebra-Tracer.wl (read its introduction for the vocabulary and the conventions of Tracer, in particular that Tracer's Eps and gamma5 have the opposite sign to FeynCalc's with $LeviCivitaSign = -1, so that every circular density below carries a minus sign).  Input: the amplitude records of step 1.  Output: Results/Tracer/Algebra.wl.  Run this file in a fresh kernel: Tracer and FeynCalc cannot be loaded together.*)
(**)
(*New here: a fermion line ending on an outgoing antiquark, ubar(p) W v(r).  Its spin sum is Tr[p-slash W r-slash Wbar] with no density matrix (both fermions are outgoing); sewTrace recognizes the Vspinor of the record.  The incoming gluon of l g -> l q qbar is summed with -g^{mu nu} (average 1/(D-2)) or replaced by the circular density; its colour sum gives Tr[T^a T^a]/(Nc^2-1) = TR.*)


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
ShowRecords[Flatten[records["RealGQ"]]]


(* ::Section:: *)
(*From records to traces*)


(* ::Input:: *)
conjugate[record_] := Module[{indices, fresh},
  indices = Union[Cases[record, i_Symbol /; StringStartsQ[SymbolName[i], "Lor"] :> i, Infinity]];   (* every Lorentz index Lor1, Lor2, ... *)
  fresh = Thread[indices -> (Symbol["c" <> SymbolName[#]] & /@ indices)];
  <|"Coefficient" -> (record["Coefficient"] /. fresh /. {Complex[x_, y_] :> Complex[x, -y],
      ColorT[a_, i_, j_] :> ColorT[Reverse[a], j, i], PolVector[p_, i_, sign_] :> PolVector[p, i, -sign]}),
    "Chains" -> (Reverse /@ record["Chains"] /. fresh)|>];

outgoing[chain_Chain] := First[Cases[chain, Ubar[p_] :> p]];
incoming[chain_Chain] := First[Cases[chain, (Uspinor | Vspinor)[p_] :> p]];
antiquarkQ[chain_Chain] := !FreeQ[chain, _Vspinor];
(* Tr[pOut-slash W rho(pIn) Wbar]: rho = pIn-slash/2 or pIn-slash gamma5/2 for an incoming fermion, the plain spin sum r-slash for an outgoing antiquark *)
sewTrace[line_, chainL_, chainR_, spin_] := Inactive[G] @@ Join[{line, outgoing[chainL]}, TracerString[chainL],
  {incoming[chainL]}, If[spin === "Polarized" && !antiquarkQ[chainL], {G5}, {}], TracerString[chainR]];   (* chainR is already reversed *)

colourSum = {ColorDelta[i_, i_] :> Nc, ColorDelta[i_, j_] /; !OrderedQ[{i, j}] :> ColorDelta[j, i],
  ColorT[{a_}, i_, x_] ColorT[{a_}, x_, j_] :> CF ColorDelta[i, j],
  ColorDelta[i_, j_]^2 :> Nc, ColorDelta[i_, j_] ColorDelta[j_, i_] :> Nc};
scalars = {Den[q_] :> 1/S[q, q], Metric[i_, j_] :> S[{i}, {j}]};

interference[termL_, termR_, spin_, bosonSums_] := Module[{conj = conjugate[termR], coefficient, traces},
  coefficient = termL["Coefficient"] conj["Coefficient"] //. colourSum //. bosonSums /. scalars;    (* //. : a product may contain two polarization pairs *)
  traces = Times @@ MapThread[sewTrace[#1, #2, #3, spin] &, {Range[Length[termL["Chains"]]], termL["Chains"], conj["Chains"]}];
  coefficient traces];

(* square of a list of amplitudes (lists of records) with the given polarization sums and the overall average *)
square[amplitudes_, spin_, bosonSums_, average_] := Module[{terms = Flatten[amplitudes]},
  average Sum[interference[terms[[i]], terms[[j]], spin, bosonSums], {i, Length[terms]}, {j, Length[terms]}]];
fourDimensional[momenta_List] := (Unprotect[H]; Scan[(H[#] = 0; H[#, _] = 0) &, momenta]; Protect[H]);
(* the colour sum of the gluon-initiated channels, CF Nc, divided by the gluon average Nc^2 - 1, is TR *)
gluonColour = {CF Nc -> TR (Nc^2 - 1)};
conjugate[records["RealGQ"][[1, 1]]]


(* ::Section:: *)
(*Born: l q -> l q*)


(* ::Input:: *)
born = Association @ Table[spin -> Internal`InheritedBlock[{S, H},
  fourDimensional[{l, k, lp, p}];
  OnShell[on, {l, 0}, {k, 0}, {lp, 0}, {p, 0}, {l, k, s/2}, {lp, p, s/2}, {l, lp, -t/2}, {k, p, -t/2}, {l, p, -u/2}, {k, lp, -u/2}];
  Factor[Expand[Activate[square[records["Born"], spin, {}, 1/(4 Nc)]]]]], {spin, {"Unpolarized", "Polarized"}}];
Require[Simplify[(born["Unpolarized"] - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];
Require[Simplify[(born["Polarized"] - 2 ee^4 eq^2 (s^2 - u^2)/t^2) /. u -> -s - t] === 0, "Polarized Born square differs from the expected result"];
born


(* ::Section:: *)
(*Real emission: l q -> l g q and l g -> l q qbar*)


(* ::Text:: *)
(*l, k, p are four-dimensional, lp and r have hatted components H[r,r] = H[lp,lp] = -H[lp,r] = spe3.  Gluon sums: the observed gluon p of the first channel with -g^{mu nu}; the incoming gluon k of the second with -g^{mu nu}/(D-2) (unpolarized) or the circular density -i Eps[mu, nu, k, l]/(2 k.l) (Tracer's sign).  Averages: 1/4 and 1/Nc for l q; 1/2 and 1/(Nc^2-1) for l g.*)


(* ::Input:: *)
realKinematics := (fourDimensional[{l, k, p}];
  Unprotect[H]; H[r, r] = spe3; H[lp, lp] = spe3; H[lp, r] = -spe3; Protect[H];
  OnShell[on, {l, 0}, {k, 0}, {lp, 0}, {p, 0}, {r, 0}, {l, k, s/2}, {k, p, -t/2}, {l, p, -u/2}, {k, r, -a/2}, {p, r, b/2},
    {l, lp, -(a + b + t)/2}, {l, r, (s + t + u + a + b)/2}, {k, lp, (s + t + a)/2}, {lp, p, -(t + u + b)/2}, {lp, r, (s + t + u)/2}]);
observedGluonSum = {PolVector[p, i_, -1] PolVector[p, j_, 1] :> -S[{i}, {j}]};
incomingGluonSum = <|"Unpolarized" -> {PolVector[k, i_, 1] PolVector[k, j_, -1] :> -S[{i}, {j}]/(D - 2)},
  "Polarized" -> {PolVector[k, i_, 1] PolVector[k, j_, -1] :> -I Eps[{i}, {j}, k, l]/(2 S[k, l])}|>;
real = <|
  "QG" -> Association @ Table[spin -> Internal`InheritedBlock[{S, H}, realKinematics;
    Factor[Expand[Activate[square[records["RealQG"], spin, observedGluonSum, 1/(4 Nc)]]]]], {spin, {"Unpolarized", "Polarized"}}],
  "GQ" -> Association @ Table[spin -> Internal`InheritedBlock[{S, H}, realKinematics;
    Factor[Expand[Activate[square[records["RealGQ"], spin, incomingGluonSum[spin], 1/(2 (Nc^2 - 1))]]] /. gluonColour]], {spin, {"Unpolarized", "Polarized"}}]|>;
Require[FreeQ[real["GQ"], CF | Nc], "The colour factor of l g -> l q qbar is not TR"];
Factor[real["GQ"]["Polarized"] /. spe3 -> 0]


(* ::Section:: *)
(*Photon Born processes: gamma q -> g q and gamma g -> q qbar*)


(* ::Text:: *)
(*All momenta four-dimensional.  The photon is averaged with -g^{mu nu}/(D-2) or circularly polarized with the density -i Eps[mu, nu, q, k]/(2 q.k); the gluon of gamma g -> q qbar likewise with the reference momentum q.*)


(* ::Input:: *)
photonKinematics := (fourDimensional[{q, k, p, r}];
  OnShell[on, {q, 0}, {k, 0}, {p, 0}, {r, 0}, {q, k, sh/2}, {k, p, -th/2}, {q, p, -uh/2}, {p, r, sh/2}, {k, r, -uh/2}, {q, r, -th/2}]);
photonSums = <|"Unpolarized" -> {PolVector[q, i_, 1] PolVector[q, j_, -1] :> -S[{i}, {j}]/(D - 2)},
  "Polarized" -> {PolVector[q, i_, 1] PolVector[q, j_, -1] :> -I Eps[{i}, {j}, q, k]/(2 S[q, k])}|>;
gluonSums = <|"Unpolarized" -> {PolVector[k, i_, 1] PolVector[k, j_, -1] :> -S[{i}, {j}]/(D - 2)},
  "Polarized" -> {PolVector[k, i_, 1] PolVector[k, j_, -1] :> -I Eps[{i}, {j}, k, q]/(2 S[k, q])}|>;
photon = <|
  "QG" -> Association @ Table[spin -> Internal`InheritedBlock[{S, H}, photonKinematics;
    Factor[Expand[Activate[square[records["PhotonQG"], spin, Join[photonSums[spin], observedGluonSum], 1/(2 Nc)]]]]], {spin, {"Unpolarized", "Polarized"}}],
  "GQ" -> Association @ Table[spin -> Internal`InheritedBlock[{S, H}, photonKinematics;
    Factor[Expand[Activate[square[records["PhotonGQ"], spin, Join[photonSums[spin], gluonSums[spin]], 1/(Nc^2 - 1)]]] /. gluonColour]], {spin, {"Unpolarized", "Polarized"}}]|>;
Require[FreeQ[photon["GQ"], CF | Nc], "The colour factor of gamma g -> q qbar is not TR"];
photon


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Photon" -> photon|>;
SaveResult["Algebra", algebra, "Tracer"];
Keys[algebra]
