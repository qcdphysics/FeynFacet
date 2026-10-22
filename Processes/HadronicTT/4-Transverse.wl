(* ::Title:: *)
(*Transversely polarized proton-proton classes (A_TT), step 4: the NLO hard coefficient of one class*)


(* ::Text:: *)
(*Input: Results/<route>/class-<n>/Algebra.wl and Reduced.wl (steps 2 and 3), whose Born table of transversely polarized 2 -> 2 squares replaces the table of Processes/PartonicBorn.  Output: Results/<route>/class-<n>/Transverse.wl with the finite hard coefficient of A_TT in the basis {regular, D0, D1, delta} and its comparison with the Fortran code of Mukherjee, Stratmann and Vogelsang (Reference/pp-transverse-fortran.wl).  The structure is that of Processes/Hadronic/4-Unpolarized.wl (read that first): the same threshold expansion, virtual normalization, coupling renormalization, convolution formulas for the two incoming legs and the fragmenting parton, pole cancellation and finite coefficient.*)
(**)
(*The hard coefficient is sigma(s,v,w) = s E dsigma/d^3p of the tagged parton (docs/CONVENTIONS.md), t = -s(1-v), u = -s v w, s4 = s v (1-w), here for the transverse double-spin difference: delta sigma_TT = (1/2)[sigma(up,up) - sigma(up,down)] projected on cos(2 Phi) (step 2: the projection TT = (XX - NN)/2),*)
(**)
(*   delta sigma_TT = alphaS^2 b(v) delta(1-w) + (alphaS^3/Pi) [ R(v,w)/w + B(v) D0(w) + C(v) D1(w) + A(v) delta(1-w) ] + O(alphaS^4),*)
(**)
(*with the Born b = |M|^2_TT/(gs^4 s v) [at D = 4], the distributions D0 = [1/(1-w)]_+ and D1 = [log(1-w)/(1-w)]_+ defined at w = 1.  Scales: muR (renormalization), muF (transversity distributions of both incoming hadrons), muFF (fragmentation).*)
(**)
(*What differs from the spin average (docs/TRANSVERSE-PP.md): (1) the HVBM scheme, whose hatted invariants step 3 has turned into integrals in D + 2 dimensions, and the four-dimensional spin vectors; (2) the collinear counterterms of the two incoming legs use the transversity kernel Delta_T P_qq(z) = CF [2 z/(1-z)_+ + 3/2 delta(1-z)] (SplittingKernel["Transverse qq"]) with no mixing (a gluon has no transversity, and the kernel connects a quark only to the same quark), so each incoming counterterm is the kernel convoluted with the class's own TT Born; (3) the fragmenting parton is unpolarized: the unpolarized kernels with the TT Born table (class 12: P_gg x Born(q qb -> g)_TT and P_gq x [Born(q qb -> q)_TT + Born(q qb -> qb)_TT + 2 (Nf-1) Born(q qb -> q')_TT]; class 7, without Born term: P_gq x Born(q q -> q)_TT only); (4) there is no finite helicity-restoring term (it belongs to the longitudinal case); (5) the evanescent Phi-independent projection "Symmetric" of step 2 is carried through the same steps as a check: its poles must cancel too, and its finite part is displayed.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
class = "5";
spin = "TT";
algebra = LoadResult["Algebra", route <> "/class-" <> class];
reduced = LoadResult["Reduced", route <> "/class-" <> class];
spec = algebra["Spec"];
Get[FileNameJoin[{NotebookDirectory[], "classes.wl"}]];
region = s > 0 && muR > 0 && muF > 0 && muFF > 0 && 0 < v < 1 && 0 < w < 1 && Nf > 0;
toVW = {t -> -s (1 - v), s4 -> s v (1 - w), u -> -s v w};
couplings = {gs -> Sqrt[4 Pi alphaS]};
(* the transversely polarized Born table of step 2 (|M|^2/gs^4 in s, t, D for the projection being assembled); a channel absent from the table is looked up as
   its charge conjugate (every quark <-> antiquark, the same transverse square), as in Processes/Hadronic; a channel absent in both forms is zero *)
conjugateName[name_] := StringReplace[name, {"qpb" -> "qp", "qp" -> "qpb", "qb" -> "q", "q" -> "qb"}];
bornOf[a_, b_, c_, projection_: spin] := Module[{name = a <> " " <> b <> " -> " <> c, table = algebra["BornTable"][projection], candidates},
  candidates = Select[{name, conjugateName[name]}, KeyExistsQ[table, #] &];
  If[candidates === {}, 0, table[First[candidates]]]];
hard[square_] := square (4 Pi alphaS)^2/(16 Pi^2 s v);        (* the hard cross section of a Born square |M|^2/gs^4 *)
spec["Name"]


(* ::Section:: *)
(*Born*)


(* ::Text:: *)
(*The TT Born hard part through order eps (the order-eps term multiplies the poles of the counterterms).  For the symmetric projection the order-eps^0 part vanishes (step 2) and only the order-eps part remains.*)


(* ::Input:: *)
bornHard = If[spec["Born"] === None, 0, Series[hard[algebra["Born"][spin]/gs^4] /. D -> 4 - 2 eps /. toVW, {eps, 0, 1}]];
born0 = If[bornHard === 0, 0, Simplify[SeriesCoefficient[bornHard, {eps, 0, 0}]]];
born1 = If[bornHard === 0, 0, Simplify[SeriesCoefficient[bornHard, {eps, 0, 1}]]];
symmetricBorn = If[spec["Born"] === None, {0, 0}, With[{ser = Series[hard[algebra["Born"]["Symmetric"]/gs^4] /. D -> 4 - 2 eps /. toVW, {eps, 0, 1}]},
  Simplify /@ {SeriesCoefficient[ser, {eps, 0, 0}], SeriesCoefficient[ser, {eps, 0, 1}]}]];
Require[symmetricBorn[[1]] === 0, "The symmetric projection of the Born term does not vanish at eps = 0", symmetricBorn];
{born0, born1, "symmetric Born: eps x " -> symmetricBorn[[2]]}


(* ::Section:: *)
(*Real emission and its threshold expansion*)


(* ::Text:: *)
(*As in Processes/Hadronic: the integrated real contribution contains s4^(-eps) and, from the unmixed families, z^(-1-eps) and z^(-eps) with z proportional to s4 = s v (1-w); ThresholdSplit makes every power of (1-w) explicit and EndpointExpand expands the branches (1-w)^(-1-lambda eps) into the distributions.  The measure of the two unobserved partons gives 1/(32 Pi^3); gs^6 -> (4 Pi alphaS)^3 with one factor muR^(2 eps) kept for the extra coupling.  The same expansion is done for the symmetric projection.*)


(* ::Input:: *)
realHardOf[projection_] := muR^(2 eps) (4 Pi alphaS)^3 reduced["Real"][projection]/gs^6/(32 Pi^3) /. toVW;
expandReal[projection_] := Module[{split = ThresholdSplit[realHardOf[projection], w, eps]},
  EndpointExpand[split["Regular"], split["Branches"], eps, w, region] /. reduced["ThresholdValues"]];   (* the inert master functions at w = 1 are numbers *)
realExpanded = expandReal[spin];
symmetricReal = expandReal["Symmetric"];
Simplify[realExpanded["Pole2"], region]


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
(*The same normalization multiplies the symmetric projection.*)


(* ::Input:: *)
virtualRawOf[projection_] := If[spec["Virtual"], (4 Pi alphaS)^3 muR^(2 eps) reduced["Virtual"][projection]/gs^6/(16 Pi^2 s v) /. toVW, 0];
virtualRaw = virtualRawOf[spin];
virtualSeries = If[virtualRaw === 0, 0, Series[virtualRaw, {eps, 0, 0}]];
virtualPole2Raw = If[virtualRaw === 0, 0, Simplify[SeriesCoefficient[virtualSeries, {eps, 0, -2}], region]];
loopNormalization = If[virtualRaw === 0, 0, 2/(4 Pi)^2 (4 Pi)^eps];          (* 2/(4 Pi)^(D/2) with D/2 = 2 - eps, phase +1 (derivation above) *)
(* the ratio of the two double poles ("Pole2" is the coefficient of delta(1-w)/eps^2) must be the constant 2/(4 Pi)^2 of the loop measure, up to a phase;
   the SU(N) relation CF = (CA^2-1)/(2 CA) is needed to see it.  The phase itself is tested by the cancellation of the double poles below *)
poleRatio = If[virtualRaw === 0, 0, Simplify[-realExpanded["Pole2"]/virtualPole2Raw /. CF -> (CA^2 - 1)/(2 CA), region]];
Require[virtualRaw === 0 || FreeQ[poleRatio, v | w | s] && MemberQ[{1, -1, I, -I}, Simplify[poleRatio/(2 (4 Pi)^(-2))]],
  "The virtual normalization is not 2/(4 Pi)^2 up to a phase", poleRatio];
virtualHard = loopNormalization virtualSeries;
virtualPole2 = If[virtualRaw === 0, 0, SeriesCoefficient[virtualHard, {eps, 0, -2}]];
(* the universal double pole of a one-loop amplitude: -(alphaS/(2 Pi)) (sum of the Casimirs of the four partons) x Born *)
casimirSum = If[spec["Born"] === None, 0, Total[If[# === "g", CA, CF] & /@ spec["Born"]]];
Require[virtualRaw === 0 || Simplify[(virtualPole2 + alphaS/(2 Pi) casimirSum born0) /. CF -> (CA^2 - 1)/(2 CA), region] === 0,
  "The virtual double pole is not -(alphaS/(2 Pi)) (C_a + C_b + C_c + C_d) times the Born term", {virtualPole2, casimirSum born0}];
virtualPole1 = If[virtualRaw === 0, 0, Simplify[SeriesCoefficient[virtualHard, {eps, 0, -1}], region]];
virtualFinite = If[virtualRaw === 0, 0, Simplify[SeriesCoefficient[virtualHard, {eps, 0, 0}], region]];
symmetricVirtual = If[virtualRaw === 0, {0, 0, 0}, With[{ser = loopNormalization Series[virtualRawOf["Symmetric"], {eps, 0, 0}]},
  Table[Simplify[SeriesCoefficient[ser, {eps, 0, n}], region], {n, -2, 0}]]];
loopNormalization


(* ::Section:: *)
(*Counterterms*)


(* ::Text:: *)
(*All counterterms multiply D-dimensional Born terms (through order eps) and carry the MSbar factor S_eps = (4 Pi)^eps e^(-eps EulerGamma) and the scale ratio (muR^2/mu^2)^eps of the subtraction: each gives a pole tuple K b0 and a finite tuple K (MSbarConstant b0 + b1 + log(muR^2/mu^2) b0), in the basis {regular, D0, delta}.  Splitting functions in units of alphaS/(2 Pi) are the kernels {regular(z), c0, cdelta} of FeynFacet/Factorization.wl.*)
(**)
(*Coupling renormalization: -(alphaS/(2 Pi)) (beta0/eps) S_eps x Born with beta0 = (11 CA - 2 Nf)/3.*)
(**)
(*Incoming leg a: + (alphaS/(2 Pi)) (1/eps) S_eps (muR^2/muF^2)^eps Delta_T P_qq convoluted with the TT Born of a b -> c (the transversity kernel connects the parton a only to itself: no sum over species), integral dz/z P(z) sigma_Born(z s, z t, u); the Born delta function fixes z = w and the convolution is P(w)/w times the Born hard part at (s, v), as a distribution [regular(w) + c0]/w + c0 [1/(1-w)]_+ + cdelta delta(1-w).*)
(**)
(*Incoming leg b: the same with the Born of a b -> c at the rescaled parton b, integral dz/z P(z) sigma_Born(z s, t, z u); the Born delta function fixes z = z0 = (1-v)/(1-vw) and, as a distribution in w,*)
(**)
(*   regular(z0) h(w) + c0 [g(w) - g(1)]/(1-w) + c0 g(1) D0(w) + [cdelta + c0 log(v/(1-v))] b0 delta(1-w),   h(w) = |M|^2(vw)/(16 Pi^2 s (1-v)),  g(w) = h(w) (1-vw)/v,  g(1) = b0,*)
(**)
(*where |M|^2(vw) is the Born square at v -> vw (Processes/Hadronic derives and checks this pullback numerically).  Transversity enters the leg-b formula only through the kernel and the TT Born: the leg-exchanged (a <-> b) Born entries of the table are the t <-> u images of the originals (checked in step 2) and the pullback uses the class's own entry at the rescaled kinematics.*)
(**)
(*Fragmentation: + (alphaS/(2 Pi)) (1/eps) S_eps (muR^2/muFF^2)^eps Sum over d of the unpolarized P_{c <- d} convoluted with the TT Born of a b -> d at the rescaled final momentum, z = 1 - v + v w with the Jacobian v w/(1-v+vw)^2; for a kernel {regular, c0, cdelta} and j(w) the Born at the shifted point times the Jacobian,*)
(**)
(*   (1/v) { c0 [j(w) - j(1)]/(1-w) + v regular(z) j(w),  c0 j(1),  j(1) [cdelta + c0 log v] },*)
(**)
(*the order-eps part of the Jacobian, z^(-2 eps), adding 2 log(z) j0 to the finite part.  A quark of flavour different from i (qp, qpb) stands for the Nf - 1 other flavours when it is summed over (class 12), exactly as in Processes/Hadronic.*)


(* ::Input:: *)
quarkQ[x_] := x =!= "g";
(* P_{c <- a} for the incoming legs: the transversity kernel, quark -> same quark only *)
incomingKernel[cp_, a_] := If[quarkQ[a] && cp === a, SplittingKernel["Transverse qq", z], {0, 0, 0}];
(* P_{c <- d} for the fragmenting parton: the unpolarized kernels *)
kernel[c_, d_] := Which[
  !quarkQ[c] && !quarkQ[d], SplittingKernel["gg", z] /. {TR -> 1/2, nf -> Nf},
  quarkQ[c] && quarkQ[d], If[c === d, SplittingKernel["qq", z], {0, 0, 0}],
  !quarkQ[c] && quarkQ[d], SplittingKernel["gq", z],
  True, SplittingKernel["qg", z] /. TR -> 1/2];
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
z0 = (1 - v)/(1 - v w);
dw = 1 - v + v w;
jFunction[bornTerm_] := v w (bornTerm /. v -> v w/dw)/dw^2;
fragmentationTuple[j_, k_] := Module[{jEnd = Simplify[j /. w -> 1, region]},
  (1/v) {k[[2]] (j - jEnd)/(1 - w) + v (k[[1]] /. z -> dw) j, k[[2]] jEnd, jEnd (k[[3]] + k[[2]] Log[v])}];
zeroTuple = <|"Pole" -> {0, 0, 0}, "Finite" -> {0, 0, 0}|>;

(* the counterterms of one projection (TT, or Symmetric for the check) *)
countertermsOf[projection_] := Module[{renormalization, legA, legB, fragmentation},
  renormalization = With[{bb = If[spec["Born"] === None, {0, 0}, bornAt[algebra["Born"][projection]/gs^4]]},
    <|"Pole" -> {0, 0, -beta0 bb[[1]]}, "Finite" -> {0, 0, -beta0 (MSbarConstant bb[[1]] + bb[[2]])}|>];
  (* leg a *)
  legA = Total[Table[With[{k = incomingKernel[cp, a], square = bornOf[cp, b, c, projection]},
    If[k === {0, 0, 0} || square === 0, zeroTuple,
      With[{bb = bornAt[square], K = multiplicity[cp, a] {((k[[1]] /. z -> w) + k[[2]])/w, k[[2]], k[[3]]}},
        <|"Pole" -> K bb[[1]], "Finite" -> finiteOf[K, bb, Log[muR^2/muF^2]]|>]]], {cp, species}]];
  (* leg b: the Born of a c' -> c at v -> v w, with the weights h and g *)
  legB = Total[Table[With[{k = incomingKernel[cp, b], square = bornOf[a, cp, c, projection]},
    If[k === {0, 0, 0} || square === 0, zeroTuple,
      Module[{b0, b1, h0, h1, g0, g1, m = multiplicity[cp, b]},
        {h0, h1} = bornPieces[square, v w] (v w)/(1 - v) // Simplify;            (* h(w) = |M|^2(v w)/(16 Pi^2 s (1-v)): the Born hard part at v w, with s v w replaced by s (1-v) *)
        {g0, g1} = {h0, h1} (1 - v w)/v // Simplify;                             (* g(w); g(1) = b0 *)
        {b0, b1} = bornAt[square];
        Require[Simplify[(g0 /. w -> 1) - b0] === 0, "g(1) differs from the Born", {g0, b0}];
        <|"Pole" -> m {(k[[1]] /. z -> z0) h0 + k[[2]] (g0 - b0)/(1 - w), k[[2]] b0, (k[[3]] + k[[2]] Log[v/(1 - v)]) b0},
          "Finite" -> m {(k[[1]] /. z -> z0) (MSbarConstant h0 + h1 + Log[muR^2/muF^2] h0) + k[[2]] ((MSbarConstant g0 + g1 + Log[muR^2/muF^2] g0) - (MSbarConstant b0 + b1 + Log[muR^2/muF^2] b0))/(1 - w),
            k[[2]] (MSbarConstant b0 + b1 + Log[muR^2/muF^2] b0), (k[[3]] + k[[2]] Log[v/(1 - v)]) (MSbarConstant b0 + b1 + Log[muR^2/muF^2] b0)}|>]]], {cp, species}]];
  (* fragmentation *)
  fragmentation = Total[Table[With[{k = kernel[c, d], square = bornOf[a, b, d, projection]},
    If[k === {0, 0, 0} || square === 0, zeroTuple,
      With[{bb = bornAt[square]},
        <|"Pole" -> multiplicity[d, c] fragmentationTuple[jFunction[bb[[1]]], k],
          "Finite" -> multiplicity[d, c] fragmentationTuple[MSbarConstant jFunction[bb[[1]]] + jFunction[bb[[2]]] + 2 jFunction[bb[[1]]] Log[dw] + Log[muR^2/muFF^2] jFunction[bb[[1]]], k]|>]]], {d, species}]];
  <|"renormalization" -> renormalization, "incoming a" -> legA, "incoming b" -> legB, "fragmentation" -> fragmentation|>];
counterterms = countertermsOf[spin];
symmetricCounterterms = countertermsOf["Symmetric"];
Map[Simplify[#, region] &, Map[#["Pole"] &, counterterms]]


(* ::Section:: *)
(*Pole cancellation*)


(* ::Text:: *)
(*The double poles of the real emission and the virtual correction cancel; the single poles of real, virtual and the four counterterms cancel in each of the three slots (regular, D0, delta).  SU(N) is needed: CF = (CA^2-1)/(2 CA).*)


(* ::Input:: *)
colourRelation = CF -> (CA^2 - 1)/(2 CA);                      (* SU(N): the cancellations hold only with it *)
Require[Simplify[(realExpanded["Pole2"] + virtualPole2) /. colourRelation, region] === 0, "Double poles do not cancel"];
poleTable = Join[{realExpanded["Pole1"], {0, 0, virtualPole1}}, alphaS/(2 Pi) Map[#["Pole"] &, Values[counterterms]]];
singlePoles = Total[poleTable] /. reduced["MasterFunctions"];        (* the pole coefficient of the mixed masters, F3[-1][z] = -2/(1-z), is a function of w *)
Require[And @@ (ExactDifference[# /. colourRelation, 0, region] === 0 & /@ singlePoles), "Single poles do not cancel", Simplify[singlePoles /. colourRelation, region]];
Grid[Prepend[MapThread[Prepend[Map[Simplify[#, region] &, #2], #1] &,
  {Join[{"real emission", "virtual"}, Keys[counterterms], {"sum"}], Append[poleTable, singlePoles]}],
  {"1/eps coefficient of", "regular", "D0", "delta"}], Alignment -> Left, Dividers -> {{2 -> True}, {2 -> True, -2 -> True}}]


(* ::Text:: *)
(*The same for the symmetric projection.  Its Born term is of order eps, so the counterterms contribute no poles (their pole tuples are K b0 with b0 = 0), and the real and virtual parts must be finite by themselves after their double poles cancel: the whole structure of the evanescent part is tested here.  Its finite part (displayed in the next section) is a Phi-independent remainder of the HVBM scheme and does not belong to the physical cross section; it is not part of the result.*)


(* ::Input:: *)
Require[Simplify[(symmetricReal["Pole2"] + symmetricVirtual[[1]]) /. colourRelation, region] === 0, "Double poles of the symmetric projection do not cancel"];
symmetricPoleTable = Join[{symmetricReal["Pole1"], {0, 0, symmetricVirtual[[2]]}}, alphaS/(2 Pi) Map[#["Pole"] &, Values[symmetricCounterterms]]];
symmetricSinglePoles = Total[symmetricPoleTable] /. reduced["MasterFunctions"];
Require[And @@ (ExactDifference[# /. colourRelation, 0, region] === 0 & /@ symmetricSinglePoles), "Single poles of the symmetric projection do not cancel", Simplify[symmetricSinglePoles /. colourRelation, region]];
"poles of the symmetric projection cancel"


(* ::Section:: *)
(*The finite hard coefficient*)


(* ::Input:: *)
finiteCounterterms = alphaS/(2 Pi) Map[{#["Finite"][[1]], #["Finite"][[2]], 0, #["Finite"][[3]]} &, counterterms];
finite = realExpanded["Finite"] + {0, 0, 0, virtualFinite} + Total[Values[finiteCounterterms]];
finite = finite /. reduced["MasterFunctions"];
prefactor = alphaS^3/(Pi s v w);
prefactorEndpoint = prefactor /. w -> 1;
coefficients = Map[Collect[Expand[#], {_Log, _PolyLog}, Factor] &, finite/{prefactor, prefactorEndpoint, prefactorEndpoint, prefactorEndpoint}];
TableForm[Transpose[{ThresholdBasis, Map[Short[#, 6] &, coefficients]}]]


(* ::Input:: *)
symmetricFinite = (symmetricReal["Finite"] + {0, 0, 0, symmetricVirtual[[3]]} + Total[Values[alphaS/(2 Pi) Map[{#["Finite"][[1]], #["Finite"][[2]], 0, #["Finite"][[3]]} &, symmetricCounterterms]]]) /. reduced["MasterFunctions"];
symmetricCoefficients = Map[Simplify[# /. colourRelation, region] &, symmetricFinite/{prefactor, prefactorEndpoint, prefactorEndpoint, prefactorEndpoint}];
"finite Phi-independent remainder of the HVBM scheme (not part of the result), coefficients {regular, D0, D1, delta}" -> symmetricCoefficients


(* ::Section:: *)
(*Scale dependence*)


(* ::Text:: *)
(*As in Processes/Collins/5-Transverse.wl: the derivatives of the finite coefficient with respect to the three scales are fixed by the renormalization group, mu^2 d/dmu^2 of the finite coefficient equals minus the corresponding pole tuple (in units of alphaS/(2 Pi)): for muF the transversity kernel convoluted with the Born on both incoming legs, for muFF the unpolarized kernels convoluted with the Born table for the fragmenting parton, for muR the beta0 term, (alphaS/(2 Pi)) beta0 x Born in the delta(1-w) slot.  The muR check uses the muR^(2 eps) of the real and virtual parts together with the counterterm logarithms and is a consequence of the pole cancellation; all three are required.*)


(* ::Input:: *)
scaleDerivative[mu_] := mu/2 D[finite, mu];                                   (* mu^2 d/dmu^2 *)
tuple4[t3_List] := {t3[[1]], t3[[2]], 0, t3[[3]]};
scaleChecks = <|
  "muF" -> Simplify[scaleDerivative[muF] + alphaS/(2 Pi) tuple4[counterterms["incoming a"]["Pole"] + counterterms["incoming b"]["Pole"]], region],
  "muFF" -> Simplify[scaleDerivative[muFF] + alphaS/(2 Pi) tuple4[counterterms["fragmentation"]["Pole"]], region],
  "muR" -> Simplify[(scaleDerivative[muR] + alphaS/(2 Pi) tuple4[counterterms["renormalization"]["Pole"]]) /. colourRelation, region]|>;
Require[And @@ (ExactDifference[#, 0, region] === 0 & /@ Flatten[Values[scaleChecks]]), "A scale derivative is not the corresponding splitting convolution (or the beta0 term)", scaleChecks];
"scale dependence verified: muF (transversity kernel on both legs), muFF (unpolarized kernels), muR (beta0)"


(* ::Section:: *)
(*Comparison with the reference code*)


(* ::Text:: *)
(*Reference/pp-transverse-fortran.wl transcribes the pieces of the Fortran code of Mukherjee, Stratmann and Vogelsang (me-transpol-ms.f) for this class as exact expressions in v, w, s, the three scales, CA, CF and Nf, and PPTransverseFortranCoefficients[J0] maps them onto this basis (the map, fixed on the Born term of class 5, is explained in the header of that file and in Reference/fortran/README.md: the code's slots hold one half of the cos(2 Phi) coefficient, so {R, B, C, A} = 2 alphaS^3/(Pi s CC) {STRUV, AVWPL(1), AVLO(1), AVDEL} and Born = 2 (1-v) F0 alphaS^2/Pi).  Both sides being exact, they are compared at exact rational points with 30 digits: several (v, w), distinct values of the three scales, two values of Nf (the dependence is linear), at SU(3) and at other values of CA with CF = (CA^2-1)/(2 CA).  The relative difference of each coefficient is taken with respect to its reference value plus 1/100 of the largest of the four at that point (a coefficient that nearly vanishes at a point is compared at absolute precision, as in Processes/Hadronic).  The SU(3) comparison is required; the general-colour one is reported (the code was written for SU(3) and the transcription cannot decide whether its CA, CF dependence is the general one).*)


(* ::Input:: *)
fortranFile = FileNameJoin[{NotebookDirectory[], "..", "..", "Reference", "pp-transverse-fortran.wl"}];
comparisonPoints = Flatten[Table[Join[{v -> vw[[1]], w -> vw[[2]], s -> 7/5, muR -> Sqrt[9/10], muF -> Sqrt[7/10], muFF -> Sqrt[13/10], Nf -> nf, alphaS -> 1}, colour],
  {vw, {{3/10, 1/5}, {11/20, 3/5}, {4/5, 9/10}, {1/2, 1/2}, {2/5, 7/10}, {19/20, 1/20}}}, {nf, {4, 5}},
  {colour, {{CA -> 3, CF -> 4/3}, {CA -> 2, CF -> 3/4}, {CA -> 5, CF -> 12/5}}}], 2];
relative[o_, t_] := With[{scale = Max[Abs[t]]/100}, MapThread[If[Abs[#1] < 10^-25 && Abs[#2] < 10^-25, 0, (#1 - #2)/(Abs[#2] + scale)] &, {o, t}]];
If[!FileExistsQ[fortranFile],
  Print["Reference/pp-transverse-fortran.wl is not available yet: the comparison with the reference code is skipped"]; comparison = None,
  Get[fortranFile];
  theirs = PPTransverseFortranCoefficients[ToExpression[class]];
  theirCoefficients = {theirs["Regular"], theirs["D0"], theirs["D1"], theirs["Delta"]};
  comparison = Table[With[{o = N[coefficients /. point, 30], th = N[theirCoefficients /. point, 30]},
    <|"v" -> (v /. point), "w" -> (w /. point), "Nf" -> (Nf /. point), "CA" -> (CA /. point), "ours" -> o, "theirs" -> th, "relative" -> relative[o, th],
      "Born" -> N[(born0 - theirs["Born"]) /. point, 30]|>], {point, comparisonPoints}];
  su3 = Select[comparison, #["CA"] == 3 &]; general = Select[comparison, #["CA"] != 3 &];
  Require[Max[Abs[Lookup[su3, "Born"]]] < 10^-20, "The Born term differs from the reference code", Lookup[su3, "Born"]];
  Require[Max[Abs[Flatten[Lookup[su3, "relative"]]]] < 10^-20, "The hard coefficient differs from the reference code at SU(3)", Grid[Map[{#["v"], #["w"], #["Nf"], NumberForm[#["relative"], 3]} &, su3]]];
  Print["largest relative difference at SU(3): ", Max[Abs[Flatten[Lookup[su3, "relative"]]]], "; at CA = 2 and 5: ", Max[Abs[Flatten[Lookup[general, "relative"]]]]]];
If[comparison =!= None, Grid[Prepend[Map[{#["v"], #["w"], #["Nf"], #["CA"], NumberForm[#["relative"], 3]} &, comparison], {"v", "w", "Nf", "CA", "differences of {R, B, C, A} (relative, see above)"}], Frame -> All]]


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
result = <|"Basis" -> ThresholdBasis, "Prefactor" -> prefactor, "Coefficients" -> coefficients, "Finite" -> finite, "Born" -> born0, "Virtual" -> virtualFinite,
  "Real" -> realExpanded["Finite"], "Counterterms" -> finiteCounterterms, "SymmetricRemainder" -> symmetricCoefficients, "ScaleChecks" -> scaleChecks,
  "Comparison" -> comparison, "Route" -> route, "Class" -> class, "Spin" -> spin,
  "Observable" -> "delta sigma_TT = (1/2)[sigma(up,up) - sigma(up,down)] projected on cos(2 Phi); sigma = s E dsigma/d^3p; scales muR, muF (both transversity densities), muFF"|>;
SaveResult["Transverse", result, route <> "/class-" <> class];
Keys[result]
