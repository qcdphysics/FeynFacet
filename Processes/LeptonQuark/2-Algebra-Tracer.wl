(* ::Title:: *)
(*Lepton-quark scattering, step 2 (Tracer route): spin and colour algebra*)


(* ::Text:: *)
(*Same task and same output format as the FeynCalc route, now with the classic Tracer package for the Dirac algebra.  Input: the amplitude records of step 1 (Results/Amplitudes.wl).  Output: Results/Tracer/Algebra.wl.  Run this file in a fresh kernel: Tracer and FeynCalc cannot be loaded together.*)
(**)
(*How a squared amplitude becomes a Tracer trace.  Every term of an amplitude is a coefficient times one open spinor chain per fermion line, ubar(pOut) W u(pIn), where W is a string of gamma matrices.  Summing over the spins of the outgoing fermion and inserting the density matrix of the incoming one gives, for the interference of term i with term j,*)
(**)
(*   Tr[ pOut-slash W_i rho(pIn) Wbar_j ],   rho = pIn-slash/2 (spin average) or pIn-slash gamma5/2 (helicity difference),*)
(**)
(*with Wbar_j the reversed string of the conjugate term.  In Tracer this trace is G[line, pOut, W_i..., pIn, (G5,) Wbar_j...] once Spur[line] has been declared.  The coefficients are multiplied by the complex conjugate of the other coefficient; gluon polarizations are summed with -g^{mu nu}; the colour sums are the few SU(N) identities written below.*)
(**)
(*Tracer vocabulary (see Tracer.pdf, the manual distributed with the package): G[line, p, {mu}, ...] is the string p-slash gamma^mu ... on the given fermion line, traced when Spur[line] has been set; S[p,q] = p.q, S[p,{mu}] = p^mu, S[{mu},{nu}] = g^{mu nu}, and repeated indices in a product are contracted automatically; H[p,q] is the (D-4)-dimensional part of p.q; Eps is the four-dimensional epsilon tensor; VectorDimension[D] sets the dimension and AntiCommute[off] selects the HVBM gamma5, defined as gamma5 = (i/4!) eps_{mu nu rho sigma} gamma^mu gamma^nu gamma^rho gamma^sigma so that Tr[gamma5 gamma^a gamma^b gamma^c gamma^d] = +4 i Eps[a,b,c,d]; this is minus the FeynCalc convention used in the other routes, which is why the circular-photon density below carries the opposite sign.  Conventions as in the FeynCalc route: s = 2 l.k, t = -2 k.p, u = -2 l.p, D = 4 - 2 eps; real emission adds a = (k-r)^2, b = (p+r)^2, s4 = (lp+r)^2, and spe3 = rhat.rhat.*)


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
ShowRecords[Flatten[records["Born"]]]


(* ::Section:: *)
(*From records to traces*)


(* ::Text:: *)
(*conjugate[record] forms the complex conjugate: i -> -i, reversed gamma strings, transposed colour generators, conjugated polarization vectors, and fresh names for the Lorentz indices so that they are not confused with those of the amplitude.  sewTrace[line, chainL, chainR, spin] builds the Tracer string of one fermion line.  The spinors of the records tell which momentum is incoming (Uspinor) and outgoing (Ubar).*)


(* ::Input:: *)
conjugate[record_] := Module[{indices, fresh},
  indices = Union[Cases[record, i_Symbol /; StringStartsQ[SymbolName[i], "Lor"] :> i, Infinity]];   (* every Lorentz index Lor1, Lor2, ... *)
  fresh = Thread[indices -> (Symbol["c" <> SymbolName[#]] & /@ indices)];
  <|"Coefficient" -> (record["Coefficient"] /. fresh /. {Complex[x_, y_] :> Complex[x, -y],
      ColorT[a_, i_, j_] :> ColorT[Reverse[a], j, i], PolVector[p_, i_, sign_] :> PolVector[p, i, -sign]}),
    "Chains" -> (Reverse /@ record["Chains"] /. fresh)|>];

outgoing[chain_Chain] := First[Cases[chain, Ubar[p_] :> p]];
incoming[chain_Chain] := First[Cases[chain, Uspinor[p_] :> p]];
sewTrace[line_, chainL_, chainR_, spin_] := Inactive[G] @@ Join[{line, outgoing[chainL]}, TracerString[chainL],
  {incoming[chainL]}, If[spin === "Polarized", {G5}, {}], TracerString[chainR]];   (* chainR is already reversed *)

conjugate[records["Born"][[1, 1]]]


(* ::Text:: *)
(*Colour.  The quark line carries either delta_{ij} (Born) or a generator (T^a)_{ij}.  Summed over external colours and averaged with 1/Nc: (T^a T^a)_{ij} = CF delta_{ij}, delta_{ii} = Nc, so that T^a_{ij} T^a_{ji} = CF Nc and delta_{ij} delta_{ij} = Nc.  Everything else in a coefficient is a scalar that Tracer understands: Den[q] = 1/q^2 -> 1/S[q,q], Metric -> S[{mu},{nu}]; polarization vectors are left for the explicit sums below.*)


(* ::Input:: *)
colourSum = {ColorDelta[i_, i_] :> Nc, ColorDelta[i_, j_] /; !OrderedQ[{i, j}] :> ColorDelta[j, i],
  ColorT[{a_}, i_, x_] ColorT[{a_}, x_, j_] :> CF ColorDelta[i, j],
  ColorDelta[i_, j_]^2 :> Nc, ColorDelta[i_, j_] ColorDelta[j_, i_] :> Nc};
scalars = {Den[q_] :> 1/S[q, q], Metric[i_, j_] :> S[{i}, {j}]};
gluonSum = PolVector[r, i_, -1] PolVector[r, j_, 1] :> -S[{i}, {j}];

(* interference of amplitude term i with conjugate term j, for the given spin choice *)
interference[termL_, termR_, spin_] := Module[{conj = conjugate[termR], coefficient, traces},
  coefficient = termL["Coefficient"] conj["Coefficient"] //. colourSum /. gluonSum /. scalars;
  traces = Times @@ MapThread[sewTrace[#1, #2, #3, spin] &, {Range[Length[termL["Chains"]]], termL["Chains"], conj["Chains"]}];
  coefficient traces];

(* square of a list of amplitudes (lists of records), averaged over the two incoming spins (1/4) and the quark colour (1/Nc) *)
square[amplitudes_, spin_] := Module[{terms = Flatten[amplitudes]},
  Sum[interference[terms[[i]], terms[[j]], spin], {i, Length[terms]}, {j, Length[terms]}]/(4 Nc)];


(* ::Section:: *)
(*Born: l q -> l q*)


(* ::Text:: *)
(*All Born momenta are four-dimensional: their hatted components H vanish.  Activate evaluates the traces; Tracer contracts the Lorentz indices shared by the two traces when the product is expanded.  The kinematics (OnShell) and the hatted components are declared inside Internal`InheritedBlock, which keeps them local to one process while retaining Tracer's own definitions of S and H.  The trace is shown before it is evaluated.*)


(* ::Input:: *)
fourDimensional[momenta_List] := (Unprotect[H]; Scan[(H[#] = 0; H[#, _] = 0) &, momenta]; Protect[H]);
bornTrace = square[records["Born"], "Polarized"]


(* ::Input:: *)
born = Association @ Table[spin -> Internal`InheritedBlock[{S, H},
  fourDimensional[{l, k, lp, p}];
  OnShell[on, {l, 0}, {k, 0}, {lp, 0}, {p, 0}, {l, k, s/2}, {lp, p, s/2}, {l, lp, -t/2}, {k, p, -t/2}, {l, p, -u/2}, {k, lp, -u/2}];
  Factor[Expand[Activate[square[records["Born"], spin]]]]], {spin, {"Unpolarized", "Polarized"}}]


(* ::Input:: *)
Require[Simplify[(born["Unpolarized"] - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];
Require[Simplify[(born["Polarized"] - 2 ee^4 eq^2 (s^2 - u^2)/t^2) /. u -> -s - t] === 0, "Polarized Born square differs from the expected result"];


(* ::Section:: *)
(*Real emission: l q -> l q g*)


(* ::Text:: *)
(*Now l, k, p are four-dimensional while the unobserved lp and r have hatted components, H[r,r] = H[lp,lp] = -H[lp,r] = spe3.  The two amplitude terms give a 2 x 2 matrix of interferences; the gluon sum -g^{mu nu} and the colour factor CF are included by interference[].*)


(* ::Input:: *)
real = Association @ Table[spin -> Internal`InheritedBlock[{S, H},
  fourDimensional[{l, k, p}];
  Unprotect[H]; H[r, r] = spe3; H[lp, lp] = spe3; H[lp, r] = -spe3; Protect[H];
  OnShell[on, {l, 0}, {k, 0}, {lp, 0}, {p, 0}, {r, 0}, {l, k, s/2}, {k, p, -t/2}, {l, p, -u/2}, {k, r, -a/2}, {p, r, b/2},
    {l, lp, -(a + b + t)/2}, {l, r, (s + t + u + a + b)/2}, {k, lp, (s + t + a)/2}, {lp, p, -(t + u + b)/2}, {lp, r, (s + t + u)/2}];
  Factor[Expand[Activate[square[records["Real"], spin]]]]], {spin, {"Unpolarized", "Polarized"}}];
Factor[real["Polarized"] /. spe3 -> 0]


(* ::Section:: *)
(*Photon-quark Born: gamma q -> q g*)


(* ::Text:: *)
(*Here only one fermion line exists.  The unpolarized photon is averaged with -g^{mu nu}/(D-2) and the quark spin with 1/2; square[] supplies 1/4, so the unpolarized result is multiplied by 2/(D-2).  For the helicity difference the photon density is -i eps^{mu nu q k}/(2 q.k) (Tracer's Eps is minus FeynCalc's with $LeviCivitaSign = -1) and the quark density k-slash gamma5/2: the factor 1/4 of square[] is replaced by 1/2.*)


(* ::Input:: *)
photonSums = {PolVector[q, i_, 1] PolVector[q, j_, -1] :> -S[{i}, {j}], PolVector[r, i_, -1] PolVector[r, j_, 1] :> -S[{i}, {j}]};
circularPhoton = {PolVector[q, i_, 1] PolVector[q, j_, -1] :> -I Eps[{i}, {j}, q, k]/(2 S[q, k]), PolVector[r, i_, -1] PolVector[r, j_, 1] :> -S[{i}, {j}]};
photon = Internal`InheritedBlock[{S, H},
  fourDimensional[{q, k, p, r}];
  OnShell[on, {q, 0}, {k, 0}, {p, 0}, {r, 0}, {q, k, sh/2}, {k, p, -th/2}, {q, p, -uh/2}, {p, r, sh/2}, {k, r, -uh/2}, {q, r, -th/2}];
  <|"Unpolarized" -> Factor[Expand[Activate[2/(D - 2) square[records["Photon"], "Unpolarized"] /. photonSums]]],
    "Polarized" -> Factor[Expand[Activate[2 square[records["Photon"], "Polarized"] /. circularPhoton]]]|>]


(* ::Section:: *)
(*Virtual: the one-loop vertex*)


(* ::Text:: *)
(*The vertex graph has the quark chain ubar(p) gamma^a (p-slash - loop-slash) gamma^mu (k-slash - loop-slash) gamma_a u(k) (indices named Lor3, Lor4 are contracted by the gluon propagator's metric).  After loop integration this current is a form factor times the Born current gamma^mu, so the ratio*)
(**)
(*   Tr[p-slash W_vertex k-slash gamma_mu] / Tr[p-slash gamma^mu k-slash gamma_mu]*)
(**)
(*is the scalar integrand of the form factor.  Written in the inverse propagators d0 = loop^2, d1 = (loop-k)^2, d2 = (loop-p)^2 (2 loop.k = d0 - d1, 2 loop.p = d0 - d2), each monomial lowers the power of a propagator: the result is a combination of VirtualIntegral[n0,n1,n2] = integral of 1/(d0^n0 d1^n1 d2^n2), the same family as in the FeynCalc route.  The coefficient of the vertex record relative to the Born record gives the normalization: gs^2 CF times the three propagators (and the factor that converts FeynArts' loop measure, fixed in step 3).*)


(* ::Input:: *)
vertex = records["VirtualVertex"][[1, 1]];
bornRecord = records["Born"][[1, 1]];
quarkChain = First[Select[vertex["Chains"], outgoing[#] === p &]];
quarkChain


(* ::Input:: *)
word = TracerString[quarkChain] /. Metric -> Sequence;
gluonIndices = Cases[vertex["Coefficient"], Metric[i_, j_] /; FreeQ[{i, j}, Lor1 | Lor2] :> {i, j}][[1]];
word = word /. {gluonIndices[[2]]} -> {gluonIndices[[1]]};         (* contract the gluon propagator metric *)
currentIndex = First[Complement[Cases[word, {i_} :> i], gluonIndices]];
projector = Internal`InheritedBlock[{S, H},
  OnShell[on, {k, 0}, {p, 0}, {k, p, -t/2}];
  Factor[(G @@ Join[{1, p}, word, {k, {currentIndex}}])/G[1, p, {currentIndex}, k, {currentIndex}]]]


(* ::Input:: *)
numerator = Expand[projector /. {S[loop, loop] -> d0, S[loop, k] -> (d0 - d1)/2, S[loop, p] -> (d0 - d2)/2}];
virtualFamily = Total[(#[[2]] VirtualIntegral @@ ({1, 1, 1} - #[[1]])) & /@ CoefficientRules[numerator, {d0, d1, d2}]]


(* ::Text:: *)
(*Normalization relative to the FeynCalc route.  There the family multiplied 2 CF Pi^2 gs^2 after dividing the interference 2 Re(M1 M0*) by |M0|^2.  Here the ratio of the two records' coefficients is CF gs^2 (with the colour sum) times the loop propagators and a factor -I from FeynArts' loop measure; the convention of step 3, which inserts the master integral with measure d^D loop/(i Pi^(D/2)), is fixed by comparing the virtual double pole with the real one.  We keep the same family normalization as the FeynCalc route so that step 3 is route independent.*)


(* ::Input:: *)
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
