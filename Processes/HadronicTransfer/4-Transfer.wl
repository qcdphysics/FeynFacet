(* ::Title:: *)
(*Transverse spin transfer in proton-proton collisions, step 4: the NLO hard coefficient of one class*)


(* ::Text:: *)
(*Same structure as Processes/Hadronic/4-Unpolarized.wl (read that first), for the transverse spin transfer a^ b -> c^: the incoming parton a transversely polarized, the transverse spin of the observed parton c measured, in the two components XX (spins in the scattering plane, along e_X for a and e_C for c) and NN (both along the normal) of docs/TRANSVERSE-PP.md.  Input: Results/<route>/class-<n>/Algebra.wl and Reduced.wl (steps 2 and 3), which hold the squares of both components and the D-dimensional transfer Born table.  Output: Results/<route>/class-<n>/Transfer.wl with the finite hard coefficients of both components in the basis {regular, D0, D1, delta}.*)
(**)
(*The hard coefficient is sigma(s,v,w) = s E dsigma/d^3p of the tagged parton (docs/CONVENTIONS.md), t = -s(1-v), u = -s v w, s4 = s v (1-w):*)
(**)
(*   sigma = alphaS^2 b(v) delta(1-w) + (alphaS^3/Pi) [ R(v,w)/w + B(v) D0(w) + C(v) D1(w) + A(v) delta(1-w) ] + O(alphaS^4),*)
(**)
(*with the Born b = |M|^2/(gs^4 s v) [at D = 4], the regular coefficient quoted with the explicit factor 1/w of the literature, the distributions D0 = [1/(1-w)]_+ and D1 = [log(1-w)/(1-w)]_+ defined at w = 1.  Scales: muR (renormalization), muF (parton distributions of both incoming hadrons), muFF (fragmentation).*)
(**)
(*What differs from the unpolarized case.  (1) The HVBM scheme: the hatted invariants became integrals in D + 2 dimensions in step 3, and the normal components of the unobserved and loop momenta were reduced in step 2.  (2) The collinear counterterms: the polarized incoming parton a and the observed parton c are subtracted with the transversity kernel Delta_T P_qq(z) = CF [2 z/(1-z)_+ + 3/2 delta(1-z)] (SplittingKernel["Transverse qq"]) and no other species (there is no gluon transversity: the transverse spin is transferred along the quark line only), the unpolarized parton b with the unpolarized kernels of all species; every counterterm multiplies the lower transfer Born of the same component, from the D-dimensional table of step 2 (a splitting does not change the spin axes: e_X is invariant under pa -> z pa and pb -> z pb, e_C and e_X change by multiples of p and pa under p -> p/z and pa -> z pa, which the projections p-slash gamma5 s-slash do not see, and the normal is the same).  (3) No helicity-restoring term (it belongs to the longitudinal case).  (4) Both components are carried through every step and checked separately; the scale derivatives are checked against the splitting convolutions.*)
(**)
(*Steps: (1) Born through order eps; (2) real emission: threshold branches (1-w)^(-1-lambda eps) and plus distributions; (3) virtual correction with its normalization; (4) counterterms: coupling renormalization, collinear subtractions of the two incoming legs and of the fragmenting parton, with the transfer Born table; (5) pole cancellation; (6) finite coefficients, scale derivatives, comparisons.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
class = "1";
algebra = LoadResult["Algebra", route <> "/class-" <> class];
reduced = LoadResult["Reduced", route <> "/class-" <> class];
spec = algebra["Spec"];
spins = algebra["Spins"];                                   (* {"XX", "NN"} *)
bornTable = algebra["BornTable"];                           (* the D-dimensional transfer squares of the 2 -> 2 channels, divided by gs^4 *)
Get[FileNameJoin[{NotebookDirectory[], "classes.wl"}]];
region = s > 0 && muR > 0 && muF > 0 && muFF > 0 && 0 < v < 1 && 0 < w < 1 && Nf > 0;
toVW = {t -> -s (1 - v), s4 -> s v (1 - w), u -> -s v w};
couplings = {gs -> Sqrt[4 Pi alphaS]};
(* the transfer Born square of a^ b -> c^ from the table (zero when there is no such 2 -> 2 process, e.g. q g -> qb).  A channel absent from the table is
   looked up as its charge conjugate (every quark <-> antiquark), which has the same transfer square (checked in step 2 on qb q -> q = q qb -> qb) *)
conjugateName[name_] := StringReplace[name, {"qpb" -> "qp", "qp" -> "qpb", "qb" -> "q", "q" -> "qb"}];
bornOf[spin_][a_, b_, c_] := Module[{name = a <> " " <> b <> " -> " <> c, candidates},
  candidates = Select[{name, conjugateName[name]}, KeyExistsQ[bornTable[spin], #] &];
  If[candidates === {}, 0, bornTable[spin][First[candidates]]]];
hard[square_] := square (4 Pi alphaS)^2/(16 Pi^2 s v);        (* the hard cross section of a Born square |M|^2/gs^4 *)
{spec["Name"], spins}


(* ::Section:: *)
(*Born*)


(* ::Input:: *)
bornHard = Association @ Table[spin -> If[spec["Born"] === None, 0, Series[hard[algebra["Born"][spin]/gs^4] /. D -> 4 - 2 eps /. toVW, {eps, 0, 1}]], {spin, spins}];
born0 = Map[If[# === 0, 0, Simplify[SeriesCoefficient[#, {eps, 0, 0}]]] &, bornHard];
born1 = Map[If[# === 0, 0, Simplify[SeriesCoefficient[#, {eps, 0, 1}]]] &, bornHard];
Require[Simplify[born0["XX"] - born0["NN"]] === 0, "The two components of the Born term differ in four dimensions", born0];
{born0, born1}


(* ::Section:: *)
(*Real emission and its threshold expansion*)


(* ::Text:: *)
(*As in Processes/Hadronic: the integrated real contribution contains s4^(-eps) and, from the unmixed families, z^(-1-eps) and z^(-eps) with z proportional to s4 = s v (1-w); ThresholdSplit makes every power of (1-w) explicit and EndpointExpand expands the branches (1-w)^(-1-lambda eps) into the distributions.  The measure of the two unobserved partons gives 1/(32 Pi^3); gs^6 -> (4 Pi alphaS)^3 with one factor muR^(2 eps) kept for the extra coupling.*)


(* ::Input:: *)
realHard = Association @ Table[spin -> (muR^(2 eps) (4 Pi alphaS)^3 reduced["Real"][spin]/gs^6/(32 Pi^3) /. toVW), {spin, spins}];
split = Map[ThresholdSplit[#, w, eps] &, realHard];
Map[Keys[#["Branches"]] &, split]


(* ::Input:: *)
realExpanded = Map[EndpointExpand[#["Regular"], #["Branches"], eps, w, region] /. reduced["ThresholdValues"] &, split];   (* the inert master functions at w = 1 are numbers *)
Map[Simplify[#["Pole2"], region] &, realExpanded]


(* ::Section:: *)
(*Virtual correction*)


(* ::Text:: *)
(*The interference of step 2 with the loop masters of step 3 inserted.  Its normalization follows from the conventions, not from the poles.  FeynArts (through FCFAConvert) writes the one-loop amplitude as M1 = integral d^D loop/(2 Pi)^D of a sum of terms c_j x (propagators), with every factor i of the vertices and propagators and the sign of a closed quark loop already in c_j; the overall phase convention of FeynArts is the same for the tree and the loop amplitude and cancels in M1 M0*.  The scalar integrals I_j of steps 2 and 3 are defined with the measure d^D loop/(i Pi^(D/2)), and*)
(**)
(*   d^D loop/(2 Pi)^D = (i Pi^(D/2)/(2 Pi)^D) d^D loop/(i Pi^(D/2)) = (i/(4 Pi)^(D/2)) d^D loop/(i Pi^(D/2)),   so   M1 = (i/(4 Pi)^(D/2)) Sum_j c_j I_j.*)
(**)
(*Step 2 stores virtual = Sum_j [i c_j M0*] I_j (the spin sums make each i c_j M0* real; step 3 takes the real parts of the I_j in the physical region), hence*)
(**)
(*   2 Re(M1 M0*) = (2/(4 Pi)^(D/2)) x virtual = (2/(4 Pi)^2) (4 Pi)^eps x virtual,*)
(**)
(*with phase +1, set explicitly below.  Two independent checks follow: here, the double pole of the virtual correction must be the universal -(alphaS/(2 Pi)) (C_a + C_b + C_c + C_d) times the Born term, C = CF for a quark and CA for a gluon (Catani, Phys. Lett. B 427 (1998) 161); in the section Pole cancellation, it must cancel the double pole of the real emission.  The ratio of the two double poles is still computed and must be 2/(4 Pi)^2 up to a phase.*)
(**)
(*The same normalization multiplies both components XX and NN; the double-pole ratio must come out the same for both.*)


(* ::Input:: *)
virtualRaw = Association @ Table[spin -> If[spec["Virtual"], (4 Pi alphaS)^3 muR^(2 eps) reduced["Virtual"][spin]/gs^6/(16 Pi^2 s v) /. toVW, 0], {spin, spins}];
virtualSeries = Map[If[# === 0, 0, Series[#, {eps, 0, 0}]] &, virtualRaw];
virtualPole2Raw = Map[If[# === 0, 0, Simplify[SeriesCoefficient[#, {eps, 0, -2}], region]] &, virtualSeries];
loopNormalization = Association @ Table[spin -> If[virtualRaw[spin] === 0, 0, 2/(4 Pi)^2 (4 Pi)^eps], {spin, spins}];   (* 2/(4 Pi)^(D/2), D/2 = 2 - eps, phase +1 (derivation above) *)
(* the ratio of the two double poles ("Pole2" is the coefficient of delta(1-w)/eps^2) must be the constant 2/(4 Pi)^2 of the loop measure, up to a phase;
   the SU(N) relation CF = (CA^2-1)/(2 CA) is needed to see it.  The phase itself is tested by the cancellation of the double poles below *)
poleRatio = Association @ Table[spin -> If[virtualRaw[spin] === 0, 0, Simplify[-realExpanded[spin]["Pole2"]/virtualPole2Raw[spin] /. CF -> (CA^2 - 1)/(2 CA), region]], {spin, spins}];
Do[Require[virtualRaw[spin] === 0 || FreeQ[poleRatio[spin], v | w | s] && MemberQ[{1, -1, I, -I}, Simplify[poleRatio[spin]/(2 (4 Pi)^(-2))]],
  "The virtual normalization is not 2/(4 Pi)^2 up to a phase", {spin, poleRatio[spin]}], {spin, spins}];
Require[poleRatio["XX"] === poleRatio["NN"], "The virtual normalization differs between the two components", poleRatio];
virtualHard = Association @ Table[spin -> loopNormalization[spin] virtualSeries[spin], {spin, spins}];
virtualPole2 = Map[If[# === 0, 0, SeriesCoefficient[#, {eps, 0, -2}]] &, virtualHard];
(* the universal double pole of a one-loop amplitude: -(alphaS/(2 Pi)) (sum of the Casimirs of the four partons) x Born, for each component *)
casimirSum = If[spec["Born"] === None, 0, Total[If[# === "g", CA, CF] & /@ spec["Born"]]];
Do[Require[virtualRaw[spin] === 0 || Simplify[(virtualPole2[spin] + alphaS/(2 Pi) casimirSum born0[spin]) /. CF -> (CA^2 - 1)/(2 CA), region] === 0,
  "The virtual double pole is not -(alphaS/(2 Pi)) (C_a + C_b + C_c + C_d) times the Born term", {spin, virtualPole2[spin], casimirSum born0[spin]}], {spin, spins}];
virtualPole1 = Map[If[# === 0, 0, Simplify[SeriesCoefficient[#, {eps, 0, -1}], region]] &, virtualHard];
virtualFinite = Map[If[# === 0, 0, Simplify[SeriesCoefficient[#, {eps, 0, 0}], region]] &, virtualHard];
loopNormalization


(* ::Section:: *)
(*Counterterms*)


(* ::Text:: *)
(*All counterterms multiply D-dimensional Born terms (through order eps) and carry the MSbar factor S_eps = (4 Pi)^eps e^(-eps EulerGamma) and the scale ratio (muR^2/mu^2)^eps of the subtraction: each gives a pole tuple K b0 and a finite tuple K (MSbarConstant b0 + b1 + log(muR^2/mu^2) b0), in the basis {regular, D0, delta}.  Splitting functions in units of alphaS/(2 Pi) are the kernels {regular(z), c0, cdelta} = regular(z) + c0 [1/(1-z)]_+ + cdelta delta(1-z) of FeynFacet/Factorization.wl.  The convolution formulas are those derived in Processes/Hadronic/4-Unpolarized.wl:*)
(**)
(*Coupling renormalization: -(alphaS/(2 Pi)) (beta0/eps) S_eps x Born, beta0 = (11 CA - 2 Nf)/3.*)
(**)
(*Incoming leg a (polarized): + (alphaS/(2 Pi)) (1/eps) S_eps (muR^2/muF^2)^eps Delta_T P_qq convoluted with the transfer Born of a b -> c, integral dz/z P(z) sigma_Born(z s, z t, u); the Born delta function fixes z = w, and P(w)/w = [regular(w) + c0]/w + c0 [1/(1-w)]_+ + cdelta delta(1-w).  Only the species of a itself contributes (a quark stays a quark of the same flavour under the transversity kernel; an antiquark an antiquark).*)
(**)
(*Incoming leg b (unpolarized): the same with the unpolarized kernels P_{c' <- b} of every species c' and the transfer Born of a c' -> c, integral dz/z P(z) sigma_Born(z s, t, z u); the Born delta function fixes z = z0 = (1-v)/(1-vw), and as a distribution in w the convolution is regular(z0) h(w) + c0 [g(w) - g(1)]/(1-w) + c0 g(1) D0(w) + [cdelta + c0 log(v/(1-v))] b0 delta(1-w), h(w) = |M|^2(vw)/(16 Pi^2 s (1-v)), g(w) = h(w) (1-vw)/v, g(1) = b0, with |M|^2(vw) the transfer Born square at v -> vw.  A quark of flavour different from i (qp, qpb) stands for the Nf - 1 other flavours when its flavour is summed over.*)
(**)
(*Fragmentation (polarized): + (alphaS/(2 Pi)) (1/eps) S_eps (muR^2/muFF^2)^eps Delta_T P_qq convoluted with the transfer Born of a b -> c at the rescaled final momentum, z = 1 - v + v w with the Jacobian v w/(1-v+vw)^2; for j(w) = the Born at the shifted point times the Jacobian, (1/v) { c0 [j(w) - j(1)]/(1-w) + v regular(z) j(w),  c0 j(1),  j(1) [cdelta + c0 log v] }, and the order-eps part of the Jacobian, z^(2 eps), adds 2 log(z) j0 to the finite part.  Only the species of c itself contributes (a gluon cannot fragment into a transversely polarized quark).*)


(* ::Input:: *)
quarkQ[x_] := x =!= "g";
(* P_{c <- a}: the unpolarized kernels for leg b *)
kernel[c_, a_] := Which[
  !quarkQ[c] && !quarkQ[a], SplittingKernel["gg", z] /. {TR -> 1/2, nf -> Nf},
  quarkQ[c] && quarkQ[a], If[c === a, SplittingKernel["qq", z], {0, 0, 0}],
  !quarkQ[c] && quarkQ[a], SplittingKernel["gq", z],
  True, SplittingKernel["qg", z] /. TR -> 1/2];
transversityKernel = SplittingKernel["Transverse qq", z];               (* CF {-2, 2, 3/2}: 2z/(1-z) = 2/(1-z) - 2 *)
species = {"q", "qb", "qp", "qpb", "g"};
(* a quark of another flavour (qp, qpb) in a counterterm stands for Nf - 1 flavours when it is summed over, i.e. when the flavour k is not one of the
   external partons of the class (then it is the fixed flavour of that external parton); this matches the weights of the real-emission processes *)
multiplicity[x_, leg_] := If[MemberQ[{"qp", "qpb"}, x] && FreeQ[{a, b, c}, "qp" | "qpb"], Nf - 1, 1];
{a, b, c} = If[spec["Born"] =!= None, spec["Born"][[1 ;; 3]], spec["Real"][[1, 1]][[1 ;; 3]]];
(* {b0, b1} of a Born square |M|^2/gs^4 (given in s, t) at the point v = vValue: |M|^2 (4 Pi alphaS)^2/(16 Pi^2 s vValue) with t = -s (1 - vValue), u = -s vValue *)
bornPieces[square_, vValue_] := Module[{ser = Series[square (4 Pi alphaS)^2/(16 Pi^2 s vValue) /. D -> 4 - 2 eps /. {t -> -s (1 - vValue), u -> -s vValue}, {eps, 0, 1}]},
  {SeriesCoefficient[ser, {eps, 0, 0}], SeriesCoefficient[ser, {eps, 0, 1}]}];
bornAt[square_] := bornPieces[square, v];
beta0 = (11 CA - 2 Nf)/3;
finiteOf[K_List, {b0_, b1_}, logScale_] := K (MSbarConstant b0 + b1 + logScale b0);
zero = <|"Pole" -> {0, 0, 0}, "Finite" -> {0, 0, 0}|>;

(* renormalization *)
renormalization[spin_] := With[{bb = If[spec["Born"] === None, {0, 0}, bornAt[algebra["Born"][spin]/gs^4]]},
  <|"Pole" -> {0, 0, -beta0 bb[[1]]}, "Finite" -> {0, 0, -beta0 (MSbarConstant bb[[1]] + bb[[2]])}|>];

(* leg a: the transversity kernel, the species of a only *)
legA[spin_] := With[{k = transversityKernel, square = bornOf[spin][a, b, c]},
  If[square === 0, zero,
    With[{bb = bornAt[square], K = {((k[[1]] /. z -> w) + k[[2]])/w, k[[2]], k[[3]]}},
      <|"Pole" -> K bb[[1]], "Finite" -> finiteOf[K, bb, Log[muR^2/muF^2]]|>]]];

(* leg b: the unpolarized kernels, the transfer Born of a c' -> c at v -> v w, with the weights h and g *)
z0 = (1 - v)/(1 - v w);
legB[spin_] := Total[Table[With[{k = kernel[cp, b], square = bornOf[spin][a, cp, c]},
  If[k === {0, 0, 0} || square === 0, zero,
    Module[{b0, b1, h0, h1, g0, g1, m = multiplicity[cp, b]},
      {h0, h1} = bornPieces[square, v w] (v w)/(1 - v) // Simplify;            (* h(w) = |M|^2(v w)/(16 Pi^2 s (1-v)): the Born hard part at v w, with s v w replaced by s (1-v) *)
      {g0, g1} = {h0, h1} (1 - v w)/v // Simplify;                             (* g(w); g(1) = b0 *)
      {b0, b1} = bornAt[square];
      Require[Simplify[(g0 /. w -> 1) - b0] === 0, "g(1) differs from the Born", {g0, b0}];
      <|"Pole" -> m {(k[[1]] /. z -> z0) h0 + k[[2]] (g0 - b0)/(1 - w), k[[2]] b0, (k[[3]] + k[[2]] Log[v/(1 - v)]) b0},
        "Finite" -> m {(k[[1]] /. z -> z0) (MSbarConstant h0 + h1 + Log[muR^2/muF^2] h0) + k[[2]] ((MSbarConstant g0 + g1 + Log[muR^2/muF^2] g0) - (MSbarConstant b0 + b1 + Log[muR^2/muF^2] b0))/(1 - w),
          k[[2]] (MSbarConstant b0 + b1 + Log[muR^2/muF^2] b0), (k[[3]] + k[[2]] Log[v/(1 - v)]) (MSbarConstant b0 + b1 + Log[muR^2/muF^2] b0)}|>]]], {cp, species}]];

(* fragmentation: the transversity kernel, the species of c only *)
dw = 1 - v + v w;
jFunction[bornTerm_] := v w (bornTerm /. v -> v w/dw)/dw^2;
fragmentationTuple[j_, k_] := Module[{jEnd = Simplify[j /. w -> 1, region]},
  (1/v) {k[[2]] (j - jEnd)/(1 - w) + v (k[[1]] /. z -> dw) j, k[[2]] jEnd, jEnd (k[[3]] + k[[2]] Log[v])}];
fragmentation[spin_] := With[{k = transversityKernel, square = bornOf[spin][a, b, c]},
  If[square === 0, zero,
    With[{bb = bornAt[square]},
      <|"Pole" -> fragmentationTuple[jFunction[bb[[1]]], k],
        "Finite" -> fragmentationTuple[MSbarConstant jFunction[bb[[1]]] + jFunction[bb[[2]]] + 2 jFunction[bb[[1]]] Log[dw] + Log[muR^2/muFF^2] jFunction[bb[[1]]], k]|>]]];
counterterms = Association @ Table[spin -> <|"renormalization" -> renormalization[spin], "incoming a" -> legA[spin], "incoming b" -> legB[spin], "fragmentation" -> fragmentation[spin]|>, {spin, spins}];
Map[Simplify[#, region] &, Map[#["Pole"] &, counterterms, {2}], {3}]


(* ::Section:: *)
(*Pole cancellation*)


(* ::Input:: *)
colourRelation = CF -> (CA^2 - 1)/(2 CA);                      (* SU(N): the cancellations hold only with it *)
Do[Require[Simplify[(realExpanded[spin]["Pole2"] + virtualPole2[spin]) /. colourRelation, region] === 0, "Double poles do not cancel", spin], {spin, spins}];
poleTable = Association @ Table[spin -> Join[{realExpanded[spin]["Pole1"], {0, 0, virtualPole1[spin]}}, alphaS/(2 Pi) Map[#["Pole"] &, Values[counterterms[spin]]]], {spin, spins}];
singlePoles = Map[Total[#] /. reduced["MasterFunctions"] &, poleTable];        (* the pole coefficient of the mixed masters, F3[-1][z] = -2/(1-z), is a function of w *)
Do[Require[And @@ (ExactDifference[# /. colourRelation, 0, region] === 0 & /@ singlePoles[spin]), "Single poles do not cancel", {spin, Simplify[singlePoles[spin] /. colourRelation, region]}], {spin, spins}];
Column[Table[Grid[Prepend[MapThread[Prepend[Map[Simplify[#, region] &, #2], #1] &,
  {Join[{"real emission", "virtual"}, Keys[counterterms[spin]], {"sum"}], Append[poleTable[spin], singlePoles[spin]]}],
  {spin <> ": 1/eps coefficient of", "regular", "D0", "delta"}], Alignment -> Left, Dividers -> {{2 -> True}, {2 -> True, -2 -> True}}], {spin, spins}]]


(* ::Section:: *)
(*The finite hard coefficients*)


(* ::Input:: *)
finiteCounterterms = Map[alphaS/(2 Pi) {#["Finite"][[1]], #["Finite"][[2]], 0, #["Finite"][[3]]} &, counterterms, {2}];
finite = Association @ Table[spin -> ((realExpanded[spin]["Finite"] + {0, 0, 0, virtualFinite[spin]} + Total[Values[finiteCounterterms[spin]]]) /. reduced["MasterFunctions"]), {spin, spins}];
prefactor = alphaS^3/(Pi s v w);
prefactorEndpoint = prefactor /. w -> 1;
coefficients = Map[Map[Collect[Expand[#], {_Log, _PolyLog}, Factor] &, #/{prefactor, prefactorEndpoint, prefactorEndpoint, prefactorEndpoint}] &, finite];
Column[Table[TableForm[Transpose[{ThresholdBasis, Map[Short[#, 6] &, coefficients[spin]]}], TableHeadings -> {None, {spin, ""}}], {spin, spins}]]


(* ::Text:: *)
(*Scale dependence.  The derivatives of the finite coefficients with respect to the three scales must reproduce the structure of the factorized cross section (with the SU(N) colour relation, as the pole cancellation): muF^2 d/dmuF^2 is minus the single-pole tuple of the two incoming-leg counterterms (the transversity convolution on a plus the unpolarized ones on b), muFF^2 d/dmuFF^2 minus that of the fragmentation counterterm (the transversity convolution on c), and muR^2 d/dmuR^2 is + (alphaS/(2 Pi)) beta0 times the Born term (the running of the two powers of alphaS in the Born term), all in the basis {regular, D0, D1, delta}.*)


(* ::Input:: *)
asTuple[{r_, d0_, de_}] := {r, d0, 0, de};
scaleDerivative[expr_, mu_] := mu/2 D[expr, mu];                                   (* mu^2 d/dmu^2 *)
(* logarithms of integers are written through their prime factors (log 4 = 2 log 2) before simplifying, so that such combinations cancel *)
primeLogs = Log[n_Integer?Positive] :> Total[#[[2]] Log[#[[1]]] & /@ FactorInteger[n]];
(* In the physical region 0 < v < 1, 1-v is positive, so
   Log[(v-1)^2] = Log[(1-v)^2] = 2 Log[1-v].  Normalize this
   identity after Simplify; do not expand logarithms without their domain. *)
physicalScaleDifference[expr_] := Module[{difference},
  difference = Simplify[expr /. colourRelation /. primeLogs, region];
  Simplify[difference /. Log[(v - 1)^2] -> 2 Log[1 - v], region]];
scaleChecks = Association @ Table[spin -> Map[physicalScaleDifference, {
    scaleDerivative[finite[spin], muF] + alphaS/(2 Pi) asTuple[counterterms[spin]["incoming a"]["Pole"] + counterterms[spin]["incoming b"]["Pole"]],
    scaleDerivative[finite[spin], muFF] + alphaS/(2 Pi) asTuple[counterterms[spin]["fragmentation"]["Pole"]],
    scaleDerivative[finite[spin], muR] - alphaS/(2 Pi) beta0 {0, 0, 0, born0[spin]}}], {spin, spins}];
Do[Require[scaleChecks[spin] === {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, "A scale derivative is not the corresponding splitting convolution", {spin, scaleChecks[spin]}], {spin, spins}];
"scale dependence verified for " <> StringRiffle[spins, " and "]


(* ::Text:: *)
(*The finite XX and NN coefficients agree at NLO after renormalization and factorization, with CF = (CA^2 - 1)/(2 CA). They must remain separate before taking the finite limit because their evanescent terms differ. The displayed difference retains CF and CA as separate symbols; impose the SU(N) relation to obtain zero.*)


(* ::Input:: *)
componentDifference = Map[Simplify[#, region] &, coefficients["XX"] - coefficients["NN"]];
TableForm[Transpose[{ThresholdBasis, Map[Short[#, 6] &, componentDifference]}]]


(* ::Section:: *)
(*Checks across classes*)


(* ::Text:: *)
(*Classes 1 (q qp -> q) and 3 (q qpb -> q) proceed through the t-channel gluon exchange only, so their Born terms must be identical: the tensor of the unpolarized line is the same for a quark and an antiquark.  At NLO they differ, as the unpolarized classes 1 and 3 do: the interference of the gluon emissions from the two lines and the box graphs distinguish a quark from an antiquark on the unpolarized line.  When the other class of the pair has already been run, its Born term is compared with this one (exactly), and the differences of the four NLO coefficients are recorded with the result.*)


(* ::Input:: *)
partnerClass = <|"1" -> "3", "3" -> "1"|>;
partnerFile = If[KeyExistsQ[partnerClass, class], FileNameJoin[{ResultsDirectory[route <> "/class-" <> partnerClass[class]], "Transfer.wl"}], None];
crossClassCheck = If[partnerFile =!= None && FileExistsQ[partnerFile],
  Module[{other = Get[partnerFile]["Finite"], otherBorn = Get[partnerFile]["Born"], points = Map[Select[#, FreeQ[#[[1]], muR | muF | muFF] &] &, PhysicalPoints], differences, bornDifference},   (* s, v, w, Nf rational; the scales symbolic *)
    bornDifference = Table[Simplify[born0[spin] - otherBorn[spin]], {spin, spins}];
    Require[bornDifference === {0, 0}, "The Born terms of classes 1 and 3 differ", bornDifference];
    differences = Table[Together[(finite[spin] - other[spin]) /. point], {spin, spins}, {point, points}];
    <|"BornAgreement" -> True, "NLODifferences" -> differences,
      "Summary" -> "Born terms of classes 1 and 3 identical; NLO coefficients differ (recorded)"|>],
  "the partner class has not been run yet: nothing to compare"]


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
result = <|"Basis" -> ThresholdBasis, "Prefactor" -> prefactor, "Coefficients" -> coefficients, "Finite" -> finite, "Born" -> born0, "Virtual" -> virtualFinite,
  "Real" -> Map[#["Finite"] &, realExpanded], "Counterterms" -> finiteCounterterms, "ComponentDifference" -> componentDifference,
  "CrossClassCheck" -> crossClassCheck, "Route" -> route, "Class" -> class, "Spins" -> spins|>;
SaveResult["Transfer", result, route <> "/class-" <> class];
Keys[result]

