(* ::Title:: *)
(*Lepton-quark and lepton-gluon scattering with a different tagged parton, step 2 (FeynCalc route): spin and colour algebra*)


(* ::Text:: *)
(*Input: Results/Amplitudes.wl from step 1.  Output: Results/FeynCalc/Algebra.wl with the squared matrix elements of*)
(*  Born:      l q -> l q (needed by the collinear counterterms of both channels),*)
(*  RealQG:    l q -> l g q, the gluon p observed, the quark r unobserved,*)
(*  RealGQ:    l g -> l q qbar, the quark p observed, the antiquark r unobserved,*)
(*  PhotonQG:  gamma q -> g q and  PhotonGQ:  gamma g -> q qbar, the photon Born processes of the photon-in-lepton subtraction,*)
(*each for the spin average ("Unpolarized") and the double-longitudinal helicity difference ("Polarized") of the incoming lepton (or photon) and parton.  Step 3 reads either route's Algebra.wl; the three routes must agree here.*)
(**)
(*Conventions, as in Processes/LeptonQuark: s = 2 l.k, t = -2 k.p, u = -2 l.p; real emission: a = (k - r)^2 = -2 k.r, b = (p + r)^2 = 2 p.r, tau = (k - p - r)^2 = a + b + t (the photon virtuality) and s4 = (lp + r)^2 = s + t + u; the photon processes use sh = 2 q.k, th = -2 k.p, uh = -2 q.p.  D = 4 - 2 eps.  Averages: 1/2 per incoming fermion, 1/(D-2) per incoming photon or gluon, 1/Nc for a quark and 1/(Nc^2 - 1) for a gluon in colour.  The helicity difference Delta sigma = (sigma(++) - sigma(+-))/2 uses the density p-slash gamma5/2 for an incoming fermion and i eps^{mu nu l r}/(2 l.r) for an incoming vector boson of momentum l (r the other incoming momentum), with gamma5 in the HVBM scheme ($LeviCivitaSign = -1): l, k, p and the recoil pRec = lp + r are four-dimensional, lp and r keep their (D-4)-dimensional components with spe3 = rhat.rhat = lphat.lphat = -lphat.rhat.*)
(**)
(*Which gluon has a polarization sum.  In l q -> l g q the gluon is the observed parton and in the photon processes it is an external parton as well; each process has a single external gluon, so its polarization sum may be taken as -g^{mu nu} (D-dimensional): the amplitude is transverse to the gluon momentum by itself, there is no ghost contribution.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[];
amplitudes = LoadAmplitudes[]["FeynCalc"];
Keys[amplitudes]


(* ::Section:: *)
(*Born: l q -> l q*)


(* ::Text:: *)
(*Exactly as in Processes/LeptonQuark (the same checks): 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - eps) for the spin average and 2 ee^4 eq^2 (s^2 - u^2)/t^2 for the helicity difference, with u = -s - t.*)


(* ::Input:: *)
FCSetDiracGammaScheme["BMHV"];
$LeviCivitaSign = -1;
couplingNames = {SMP["e"] -> ee, SMP["g_s"] -> gs};          (* FeynCalc's names of the couplings -> ours *)
FCClearScalarProducts[];
SetMandelstam[s, t, u, l, k, -lp, -p, 0, 0, 0, 0];
Scan[(Momentum[#, D - 4] = 0) &, {l, k, lp, p}];
bornAmplitude = Total[amplitudes["Born"]];
bornUnpolarized = bornAmplitude ComplexConjugate[bornAmplitude]/SUNN //
  FeynAmpDenominatorExplicit // SUNSimplify[#, Explicit -> True] & //
  FermionSpinSum[#, ExtraFactor -> 1/4] & // DiracSimplify[#, DiracTraceEvaluate -> True] & //
  Contract // TrickMandelstam[#, {s, t, u, 0}] & // Factor;
bornUnpolarized = bornUnpolarized /. couplingNames;
helicityTag = {Spinor[Momentum[l, D], m_, o_] :> Spinor[Momentum[lHel, D], m, o],
  Spinor[Momentum[k, D], m_, o_] :> Spinor[Momentum[kHel, D], m, o]};
helicityDensity = {DiracGamma[Momentum[lHel, D], D] :> GSD[l].GA[5], DiracGamma[Momentum[kHel, D], D] :> GSD[k].GA[5]};
bornPolarized = (bornAmplitude ComplexConjugate[bornAmplitude]/SUNN /. helicityTag) //
  FeynAmpDenominatorExplicit // SUNSimplify[#, Explicit -> True] & //
  FermionSpinSum[#, ExtraFactor -> 1/4] & // ReplaceAll[#, helicityDensity] & //
  DiracSimplify[#, DiracTraceEvaluate -> True] & // Contract // ExpandScalarProduct //
  TrickMandelstam[Cancel[Together[#]], {s, t, u, 0}] & // Factor;
bornPolarized = bornPolarized /. couplingNames;
Require[Simplify[(bornUnpolarized - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];
Require[Simplify[(bornPolarized - 2 ee^4 eq^2 (s^2 - u^2)/t^2) /. u -> -s - t] === 0, "Polarized Born square differs from the expected result"];
born = <|"Unpolarized" -> bornUnpolarized, "Polarized" -> bornPolarized|>


(* ::Section:: *)
(*Real emission, spin average*)


(* ::Text:: *)
(*The five-particle kinematics is set through all scalar products.  Both channels are squared with FeynCalc's generic tools: FermionSpinSum closes the two fermion lines (lepton and quark; for l g -> l q qbar the quark line runs from the antiquark to the quark), the gluon is summed with -g^{mu nu}, the colour sums are explicit.  Averages: 1/4 and 1/Nc for l q, 1/2 x 1/(D-2) and 1/(Nc^2-1) for l g.  The colour factors are Tr[T^a T^a]/Nc = CF and Tr[T^a T^a]/(Nc^2-1) = 1/2; the latter is written as the symbol TR (FeynCalc would evaluate it to 1/2).*)


(* ::Input:: *)
FCSetDiracGammaScheme["NDR"];
realKinematics := (FCClearScalarProducts[];
  SPD[l] = 0; SPD[k] = 0; SPD[lp] = 0; SPD[p] = 0; SPD[r] = 0;
  SPD[l, k] = s/2; SPD[k, p] = -t/2; SPD[l, p] = -u/2; SPD[k, r] = -a/2; SPD[p, r] = b/2;
  SPD[l, lp] = -(a + b + t)/2; SPD[l, r] = (s + t + u + a + b)/2; SPD[k, lp] = (s + t + a)/2;
  SPD[lp, p] = -(t + u + b)/2; SPD[lp, r] = (s + t + u)/2);
realKinematics;
squareUnpolarized[amps_, gluon_, fermionAverage_, colourAverage_] := Module[{amp},
  amp = Total[SUNSimplify[#, SUNNToCACF -> True] & /@ amps] // DiracSimplify // Contract;
  amp ComplexConjugate[amp] colourAverage // FeynAmpDenominatorExplicit //
    SUNSimplify[#, Explicit -> True, SUNNToCACF -> True] & // FermionSpinSum[#, ExtraFactor -> fermionAverage] & //
    DiracSimplify // Contract // DoPolarizationSums[#, gluon, 0, VirtualBoson -> True] & // Contract // Simplify //
    SUNSimplify[#, SUNNToCACF -> True] & // Factor // ReplaceAll[#, couplingNames] &];
realUnpolarized = <|
  "QG" -> squareUnpolarized[amplitudes["RealQG"], p, 1/4, 1/SUNN],
  "GQ" -> squareUnpolarized[amplitudes["RealGQ"], k, 1/2, 2 TR/(D - 2)/(SUNN^2 - 1)]|>;     (* colour: Tr[T^a T^a]/(Nc^2-1) = 1/2, written as TR *)
Require[FreeQ[realUnpolarized, SUNN | SUNIndex | Polarization | Spinor], "The unpolarized real squares are not scalar"];
realUnpolarized["QG"]


(* ::Section:: *)
(*Real emission, helicity difference*)


(* ::Text:: *)
(*As in the lepton-quark process the traces are written out explicitly (FeynCalc's spin sum would insert p-slash, not the helicity densities).  All gamma matrices and Lorentz indices are D-dimensional; the four-dimensional objects are the physical momenta and the Levi-Civita tensors of the densities (in a trace with gamma5 the (D-4)-dimensional components of the photon index drop out anyway, so this agrees with the four-dimensional photon index used in Processes/LeptonQuark).  The lepton trace is Tr[lp-slash gamma^mu l-slash gamma5 gamma^nu]/2.  The quark line of each diagram is a string W_i between the external spinors, with the photon index mu and the gluon index rho (W-bar_i is the reversed string with nu, sigma):*)
(**)
(*  l q -> l g q:    W_1 = gamma^mu (k-slash - p-slash) gamma^rho / t,     W_2 = gamma^rho (r-slash + p-slash) gamma^mu / b,     trace Tr[r-slash W_i k-slash gamma5 W-bar_j] / 2;*)
(*  l g -> l q qbar: W_1 = gamma^mu (k-slash - r-slash) gamma^rho / a,     W_2 = gamma^rho (p-slash - k-slash) gamma^mu / t,     trace Tr[p-slash W_i r-slash W-bar_j],*)
(**)
(*the gluon is summed with -g_{rho sigma} in D dimensions or, when it is the polarized incoming gluon, replaced by the circular density i eps^{rho sigma k l}/(2 k.l).  The common factors are the photon propagators 1/tau^2, the couplings ee^4 eq^2 gs^2 and the colour factors CF (Tr[T^a T^a]/Nc) and TR (Tr[T^a T^a]/(Nc^2-1)).  The same function with gamma5 removed and the spin averages inserted (1/2 for the quark, -g/(D-2) for the gluon) must reproduce the squares of the previous cell: this checks the strings, propagators and signs.  The unobserved momenta lp and r are D-dimensional (GSD), the observed ones four-dimensional (GS); the (D-4)-dimensional products are spe3.*)


(* ::Input:: *)
FCSetDiracGammaScheme["BMHV"];
hvbmKinematics := (realKinematics;
  Scan[(Momentum[#, D - 4] = 0) &, {l, k, p, pRec}];
  SPE[lp] = spe3; SPE[r] = spe3; SPE[lp, r] = -spe3;
  SP[l] = 0; SP[k] = 0; SP[p] = 0; SP[lp] = -spe3; SP[r] = -spe3;
  SP[l, k] = s/2; SP[k, p] = -t/2; SP[l, p] = -u/2; SP[k, r] = -a/2; SP[p, r] = b/2;
  SP[l, lp] = -(a + b + t)/2; SP[l, r] = (s + t + u + a + b)/2; SP[k, lp] = (s + t + a)/2;
  SP[lp, p] = -(t + u + b)/2; SP[lp, r] = (s + t + u)/2 + spe3);
hvbmKinematics;
reverse[string_] := (string /. {mu -> nu, rho -> sigma}) /. Dot -> (Dot @@ Reverse[{##}] &);
simplifyTrace[expr_] := Factor[Together[Simplify[Contract[ExpandScalarProduct[Contract[DiracSimplify[expr, DiracTraceEvaluate -> True]]]]]]];
(* density matrices of the incoming particles: polarized or averaged *)
leptonDensity[spin_] := If[spin === "Polarized", GSD[l].GA[5]/2, GSD[l]/2];
quarkDensity[spin_] := If[spin === "Polarized", GS[k].GA[5]/2, GS[k]/2];
gluonDensity[spin_] := If[spin === "Polarized", Contract[I LC[rho, sigma, al, be] FV[k, al] FV[l, be]/(2 SP[k, l])], -MTD[rho, sigma]/(D - 2)];
leptonTrace[spin_] := DiracTrace[GSD[lp].GAD[mu].leptonDensity[spin].GAD[nu]];
realByTraces["QG", spin_] := Module[{strings, traces},
  strings = {GAD[mu].(GSD[k] - GSD[p]).GAD[rho]/t, GAD[rho].(GSD[r] + GSD[p]).GAD[mu]/b};
  traces = Table[simplifyTrace[leptonTrace[spin] DiracTrace[GSD[r].strings[[i]].quarkDensity[spin].reverse[strings[[j]]]] (-MTD[rho, sigma])], {i, 2}, {j, 2}];
  Factor[Together[CF ee^4 eq^2 gs^2/(a + b + t)^2 Total[traces, 2]]]];
realByTraces["GQ", spin_] := Module[{strings, traces},
  strings = {GAD[mu].(GSD[k] - GSD[r]).GAD[rho]/a, GAD[rho].(GSD[p] - GSD[k]).GAD[mu]/t};
  traces = Table[simplifyTrace[leptonTrace[spin] DiracTrace[GS[p].strings[[i]].GSD[r].reverse[strings[[j]]]] gluonDensity[spin]], {i, 2}, {j, 2}];
  Factor[Together[TR ee^4 eq^2 gs^2/(a + b + t)^2 Total[traces, 2]]]];
Do[Require[Together[(realByTraces[channel, "Unpolarized"] /. spe3 -> 0) - realUnpolarized[channel]] === 0,
  "The explicit traces do not reproduce the unpolarized square", channel], {channel, {"QG", "GQ"}}];
realPolarized = <|"QG" -> realByTraces["QG", "Polarized"], "GQ" -> realByTraces["GQ", "Polarized"]|>;
real = <|"QG" -> <|"Unpolarized" -> realUnpolarized["QG"], "Polarized" -> realPolarized["QG"]|>,
  "GQ" -> <|"Unpolarized" -> realUnpolarized["GQ"], "Polarized" -> realPolarized["GQ"]|>|>;
Factor[realPolarized["GQ"] /. spe3 -> 0]


(* ::Section:: *)
(*Photon Born processes: gamma q -> g q and gamma g -> q qbar*)


(* ::Text:: *)
(*Needed through order eps for the photon-in-lepton subtraction.  Invariants sh = 2 q.k, th = -2 k.p, uh = -2 q.p with sh + th + uh = 0; all momenta are four-dimensional.  Spin average: FeynCalc's generic tools with the averages 1/2 (quark), 1/(D-2) (photon, and the gluon of gamma g -> q qbar), 1/Nc or 1/(Nc^2-1).  Helicity difference: explicit traces with the circular photon density i eps^{mu nu q k}/(2 q.k), the quark density k-slash gamma5/2 or the circular gluon density i eps^{rho sigma k q}/(2 k.q), the strings*)
(**)
(*  gamma q -> g q:    W_1 = gamma^rho (k-slash + q-slash) gamma^mu / sh,   W_2 = gamma^mu (k-slash - p-slash) gamma^rho / th,   trace Tr[r-slash W_i k-slash gamma5 W-bar_j]/2, gluon sum -g;*)
(*  gamma g -> q qbar: W_1 = gamma^mu (p-slash - q-slash) gamma^rho / uh,   W_2 = gamma^rho (p-slash - k-slash) gamma^mu / th,   trace Tr[p-slash W_i r-slash W-bar_j],*)
(**)
(*and the same consistency check against the generic spin average.*)


(* ::Input:: *)
photonKinematics := (FCClearScalarProducts[];
  Scan[(Momentum[#, D - 4] = 0) &, {q, k, p, r}];
  SPD[q] = 0; SPD[k] = 0; SPD[p] = 0; SPD[r] = 0;
  SPD[q, k] = sh/2; SPD[k, p] = -th/2; SPD[q, p] = -uh/2; SPD[p, r] = sh/2; SPD[k, r] = -uh/2; SPD[q, r] = -th/2;
  SP[q] = 0; SP[k] = 0; SP[p] = 0; SP[r] = 0;
  SP[q, k] = sh/2; SP[k, p] = -th/2; SP[q, p] = -uh/2; SP[p, r] = sh/2; SP[k, r] = -uh/2; SP[q, r] = -th/2);
FCSetDiracGammaScheme["NDR"];
photonKinematics;
squarePhoton[amps_, bosons_, fermionAverage_, colourAverage_] := Module[{amp, expr},
  amp = Total[SUNSimplify[#, SUNNToCACF -> True] & /@ amps] // DiracSimplify // Contract;
  expr = amp ComplexConjugate[amp] colourAverage // FeynAmpDenominatorExplicit //
    SUNSimplify[#, Explicit -> True, SUNNToCACF -> True] & // FermionSpinSum[#, ExtraFactor -> fermionAverage] & // DiracSimplify // Contract;
  Do[expr = DoPolarizationSums[expr, boson, 0, VirtualBoson -> True] // Contract, {boson, bosons}];
  TrickMandelstam[Factor[SUNSimplify[Simplify[expr], SUNNToCACF -> True]], {sh, th, uh, 0}] /. couplingNames];
photonUnpolarized = <|
  "QG" -> squarePhoton[amplitudes["PhotonQG"], {q, p}, 1/2, 1/SUNN]/(D - 2),
  "GQ" -> squarePhoton[amplitudes["PhotonGQ"], {q, k}, 1, 2 TR/(SUNN^2 - 1)]/(D - 2)^2|>;
Require[FreeQ[photonUnpolarized, SUNN | SUNIndex | Polarization | Spinor], "The unpolarized photon squares are not scalar"];
photonUnpolarized


(* ::Input:: *)
FCSetDiracGammaScheme["BMHV"];
photonKinematics;
simplifyPhoton[expr_] := Factor[TrickMandelstam[Cancel[Together[Simplify[Contract[ExpandScalarProduct[Contract[DiracSimplify[expr, DiracTraceEvaluate -> True]]]]]]], {sh, th, uh, 0}]];
photonDensity[spin_] := If[spin === "Polarized", Contract[I LC[mu, nu, al, be] FV[q, al] FV[k, be]/(2 SP[q, k])], -MTD[mu, nu]/(D - 2)];
gluonDensityPhoton[spin_] := If[spin === "Polarized", Contract[I LC[rho, sigma, al, be] FV[k, al] FV[q, be]/(2 SP[k, q])], -MTD[rho, sigma]/(D - 2)];
photonByTraces["QG", spin_] := Module[{strings, traces},
  strings = {GAD[rho].(GSD[k] + GSD[q]).GAD[mu]/sh, GAD[mu].(GSD[k] - GSD[p]).GAD[rho]/th};
  traces = Table[simplifyPhoton[photonDensity[spin] DiracTrace[GS[r].strings[[i]].quarkDensity[spin].reverse[strings[[j]]]] (-MTD[rho, sigma])], {i, 2}, {j, 2}];
  simplifyPhoton[CF ee^2 eq^2 gs^2 Total[traces, 2]]];
photonByTraces["GQ", spin_] := Module[{strings, traces},
  strings = {GAD[mu].(GSD[p] - GSD[q]).GAD[rho]/uh, GAD[rho].(GSD[p] - GSD[k]).GAD[mu]/th};
  traces = Table[simplifyPhoton[photonDensity[spin] DiracTrace[GS[p].strings[[i]].GS[r].reverse[strings[[j]]]] gluonDensityPhoton[spin]], {i, 2}, {j, 2}];
  simplifyPhoton[TR ee^2 eq^2 gs^2 Total[traces, 2]]];
Do[Require[Simplify[(photonByTraces[channel, "Unpolarized"] - photonUnpolarized[channel]) /. uh -> -sh - th] === 0,
  "The explicit traces do not reproduce the unpolarized photon square", channel], {channel, {"QG", "GQ"}}];
photonPolarized = <|"QG" -> photonByTraces["QG", "Polarized"], "GQ" -> photonByTraces["GQ", "Polarized"]|>;
photon = <|"QG" -> <|"Unpolarized" -> photonUnpolarized["QG"], "Polarized" -> photonPolarized["QG"]|>,
  "GQ" -> <|"Unpolarized" -> photonUnpolarized["GQ"], "Polarized" -> photonPolarized["GQ"]|>|>;
photonPolarized


(* ::Section:: *)
(*Save*)


(* ::Text:: *)
(*Algebra.wl holds everything step 3 needs: the Born square (function of D), the real squares of both channels (functions of D, a, b, s, t, u and spe3) and the photon squares (functions of D, sh, th, uh).*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Photon" -> photon|>;
SaveResult["Algebra", algebra, "FeynCalc"];
Keys[algebra]
