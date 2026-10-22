(* ::Title:: *)
(*Lepton-quark scattering, step 2 (FeynCalc route): spin and colour algebra*)


(* ::Text:: *)
(*Input: Results/Amplitudes.wl from step 1.  Output: Results/FeynCalc/Algebra.wl, the squared matrix elements as scalar functions of the invariants, for the spin average ("Unpolarized") and the incoming double-helicity difference ("Polarized").  Step 3 reads either route's Algebra.wl; the three routes must agree here.*)
(**)
(*Conventions.  Massless kinematics with s = 2 l.k, t = -2 k.p, u = -2 l.p.  For real emission the extra invariants are a = (k - r)^2 = -2 k.r, b = (p + r)^2 = 2 p.r, tau = (k - p - r)^2 = a + b + t (photon virtuality) and s4 = (lp + r)^2 = s + t + u.  Dimensional regularization with D = 4 - 2 eps.*)
(**)
(*Spin.  The spin average is 1/2 per incoming fermion, the colour average 1/Nc for the quark; final spins and colours are summed.  The double-longitudinal difference Delta sigma = (sigma(++) - sigma(+-))/2 is obtained by replacing the density matrix p-slash/2 of each incoming fermion by p-slash gamma5/2; gamma5 is treated in the HVBM (BMHV) scheme, in which the external momenta l, k, p are four-dimensional while the unobserved momenta lp, r and the loop momentum keep their (D-4)-dimensional ("hatted") components.  spe3 denotes the hatted invariant rhat.rhat = lphat.lphat = -lphat.rhat of the unobserved pair.*)


(* ::Section:: *)
(*Load*)


(* ::Text:: *)
(*The package must be loaded in a cell of its own, before any cell that uses its functions is parsed.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[];
amplitudes = LoadAmplitudes[]["FeynCalc"];
Keys[amplitudes]


(* ::Section:: *)
(*Born: l q -> l q*)


(* ::Text:: *)
(*Unpolarized square: sum over the two fermion lines with FermionSpinSum, average 1/4 over the incoming spins and 1/Nc over the incoming quark colour.  The D-dimensional result keeps its (D-4) term, which multiplies collinear poles later.*)


(* ::Input:: *)
FCClearScalarProducts[];
SetMandelstam[s, t, u, l, k, -lp, -p, 0, 0, 0, 0];
bornAmplitude = Total[amplitudes["Born"]];
bornUnpolarized = bornAmplitude ComplexConjugate[bornAmplitude]/SUNN //
  FeynAmpDenominatorExplicit // SUNSimplify[#, Explicit -> True] & //
  FermionSpinSum[#, ExtraFactor -> 1/4] & // DiracSimplify[#, DiracTraceEvaluate -> True] & //
  Contract // TrickMandelstam[#, {s, t, u, 0}] & // Factor;
bornUnpolarized = bornUnpolarized /. {SMP["e"] -> ee, SMP["g_s"] -> gs}


(* ::Text:: *)
(*Expected: 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - eps) with D = 4 - 2 eps (TrickMandelstam has used s + t + u = 0, which holds for the 2 -> 2 kinematics).*)


(* ::Input:: *)
Require[Simplify[(bornUnpolarized - 2 ee^4 eq^2 ((s^2 + u^2)/t^2 - (4 - D)/2)) /. u -> -s - t] === 0, "Born square differs from the textbook result"];


(* ::Text:: *)
(*Helicity difference.  FeynCalc's FermionSpinSum inserts p-slash for every spinor; to get p-slash gamma5 for the two incoming fermions we temporarily rename their momenta (lHel, kHel) and replace the slashes it produces by GSD[l].GA[5] and GSD[k].GA[5].  FCSetDiracGammaScheme["BMHV"] selects the HVBM scheme; $LeviCivitaSign = -1 fixes Tr[gamma^a gamma^b gamma^c gamma^d gamma5] = -4 i eps^{abcd}, the convention in which Tracer and FORM below are written.  Momentum[x, D-4] = 0 declares x four-dimensional.*)


(* ::Input:: *)
FCSetDiracGammaScheme["BMHV"];
$LeviCivitaSign = -1;
FCClearScalarProducts[];
SetMandelstam[s, t, u, l, k, -lp, -p, 0, 0, 0, 0];
Scan[(Momentum[#, D - 4] = 0) &, {l, k, lp, p}];
helicityTag = {Spinor[Momentum[l, D], m_, o_] :> Spinor[Momentum[lHel, D], m, o],
  Spinor[Momentum[k, D], m_, o_] :> Spinor[Momentum[kHel, D], m, o]};
helicityDensity = {DiracGamma[Momentum[lHel, D], D] :> GSD[l].GA[5],
  DiracGamma[Momentum[kHel, D], D] :> GSD[k].GA[5]};
bornPolarized = (bornAmplitude ComplexConjugate[bornAmplitude]/SUNN /. helicityTag) //
  FeynAmpDenominatorExplicit // SUNSimplify[#, Explicit -> True] & //
  FermionSpinSum[#, ExtraFactor -> 1/4] & // ReplaceAll[#, helicityDensity] & //
  DiracSimplify[#, DiracTraceEvaluate -> True] & // Contract // ExpandScalarProduct //
  TrickMandelstam[Cancel[Together[#]], {s, t, u, 0}] & // Factor;
bornPolarized = bornPolarized /. {SMP["e"] -> ee, SMP["g_s"] -> gs}


(* ::Text:: *)
(*Expected: 2 ee^4 eq^2 (s^2 - u^2)/t^2, with no (D-4) term.*)


(* ::Input:: *)
Require[Simplify[(bornPolarized - 2 ee^4 eq^2 (s^2 - u^2)/t^2) /. u -> -s - t] === 0, "Polarized Born square differs from the expected result"];
born = <|"Unpolarized" -> bornUnpolarized, "Polarized" -> bornPolarized|>;


(* ::Section:: *)
(*Virtual: the one-loop vertex correction*)


(* ::Text:: *)
(*The interference 2 Re(M_vertex M_Born^*) of the vertex graph with the Born amplitude is a vector-current form factor times the Born square, so it is the same for both spin cases.  TID expresses the tensor loop integral through the scalar integrals B0 and C0 (Passarino-Veltman functions); dividing by the Born square leaves a scalar expression.  We write it in terms of the scalar integral family*)
(**)
(*   VirtualIntegral[n0, n1, n2] = integral over the loop momentum of 1/((loop^2)^n0 ((loop-k)^2)^n1 ((loop-p)^2)^n2),*)
(**)
(*which Kira reduces to a single master integral in step 3.  The on-shell bubble B0[0,0,0] is scaleless and vanishes; keeping it visible ($KeepLogDivergentScalelessIntegrals) documents where it would enter.*)


(* ::Input:: *)
FCSetDiracGammaScheme["NDR"];
FCClearScalarProducts[];
SetMandelstam[s, t, u, l, k, -lp, -p, 0, 0, 0, 0];
$KeepLogDivergentScalelessIntegrals = True;
vertexAmplitude = Total[amplitudes["VirtualVertex"]] // DiracSimplify // Contract // SUNSimplify[#, SUNNToCACF -> True] &;
vertexScalar = vertexAmplitude // TID[#, loop, ToPaVe -> True] & // DiracSimplify // Contract //
  SUNSimplify[#, SUNNToCACF -> True] & // FeynAmpDenominatorExplicit;
bornReduced = bornAmplitude // DiracSimplify // Contract // FeynAmpDenominatorExplicit;
square[x_, y_] := x ComplexConjugate[y] // FermionSpinSum[#, ExtraFactor -> 1/4] & // DiracSimplify //
  SUNSimplify[#, SUNNToCACF -> True] & // TrickMandelstam[#, {s, t, u, 0}] & // Simplify;
virtualRatio = Factor[square[2 vertexScalar, bornReduced]/square[bornReduced, bornReduced]]


(* ::Text:: *)
(*The loop measure of TID is d^D loop/(2 pi)^D; the family is defined with measure d^D loop/(i pi^(D/2)), hence the factor 1/(2 CF Pi^2 gs^2) that also removes the colour factor CF and the coupling (the factor 2 is the 2 Re).  The master integral normalization is fixed in step 3.*)


(* ::Input:: *)
virtualFamily = Factor[virtualRatio/(2 CF Pi^2 SMP["g_s"]^2)] /. {
  B0[0, 0, 0] -> VirtualIntegral[1, 1, 0], B0[t, 0, 0] -> VirtualIntegral[0, 1, 1],
  C0[0, 0, t, 0, 0, 0] -> VirtualIntegral[1, 1, 1]}


(* ::Section:: *)
(*Real emission: l q -> l q g*)


(* ::Text:: *)
(*Unpolarized: square the sum of the two diagrams, sum the gluon polarizations with -g^{mu nu} (DoPolarizationSums with VirtualBoson -> True), average over spins and colour.  The five-particle kinematics is specified through all scalar products.*)


(* ::Input:: *)
FCSetDiracGammaScheme["NDR"];
FCClearScalarProducts[];
SPD[l] = 0; SPD[k] = 0; SPD[lp] = 0; SPD[p] = 0; SPD[r] = 0;
SPD[l, k] = s/2; SPD[k, p] = -t/2; SPD[l, p] = -u/2; SPD[k, r] = -a/2; SPD[p, r] = b/2;
SPD[l, lp] = -(a + b + t)/2; SPD[l, r] = (s + t + u + a + b)/2; SPD[k, lp] = (s + t + a)/2;
SPD[lp, p] = -(t + u + b)/2; SPD[lp, r] = (s + t + u)/2;
realAmplitude = Total[SUNSimplify[#, SUNNToCACF -> True] & /@ amplitudes["Real"]] // DiracSimplify // Contract;
realUnpolarized = realAmplitude ComplexConjugate[realAmplitude]/SUNN // FeynAmpDenominatorExplicit //
  SUNSimplify[#, Explicit -> True, SUNNToCACF -> True] & // FermionSpinSum[#, ExtraFactor -> 1/4] & //
  DiracSimplify // Contract // DoPolarizationSums[#, r, 0, VirtualBoson -> True] & // Contract // Simplify //
  SUNSimplify[#, SUNNToCACF -> True] & // Factor;
realUnpolarized = realUnpolarized /. {SMP["e"] -> ee, SMP["g_s"] -> gs}


(* ::Text:: *)
(*Helicity difference in the HVBM scheme.  Here l, k, p and the recoil momentum pRec = lp + r are four-dimensional, while lp and r individually have (D-4)-dimensional components with rhat.rhat = lphat.lphat = -lphat.rhat = spe3.  FeynCalc distinguishes D-dimensional products SPD, four-dimensional ones SP and (D-4)-dimensional ones SPE; all three are set.*)
(**)
(*Rather than the generic spin sum we write the two traces explicitly: the lepton trace Tr[lp-slash gamma^mu l-slash gamma5 gamma^nu]/2 and, for each pair (i, j) of the two quark currents, Tr[p-slash W_i k-slash gamma5 Wbar_j]/2, where W_1 = gamma^mu (k-slash - r-slash) gamma^rho and W_2 = gamma^rho (p-slash + r-slash) gamma^mu are the quark strings of the two diagrams (photon index mu, gluon index rho) and Wbar their reversed strings.  The factors 1/2 are the two helicity density matrices.  The photon indices mu, nu are four-dimensional (GA), the gluon indices D-dimensional (GAD); the gluon sum is -g_{rho sigma}.*)


(* ::Input:: *)
FCSetDiracGammaScheme["BMHV"];
FCClearScalarProducts[];
Scan[(Momentum[#, D - 4] = 0) &, {l, k, p, pRec}];
SPD[l] = 0; SPD[k] = 0; SPD[lp] = 0; SPD[p] = 0; SPD[r] = 0;
SPD[l, k] = s/2; SPD[k, p] = -t/2; SPD[l, p] = -u/2; SPD[k, r] = -a/2; SPD[p, r] = b/2;
SPD[l, lp] = -(a + b + t)/2; SPD[l, r] = (s + t + u + a + b)/2; SPD[k, lp] = (s + t + a)/2;
SPD[lp, p] = -(t + u + b)/2; SPD[lp, r] = (s + t + u)/2;
SPE[lp] = spe3; SPE[r] = spe3; SPE[lp, r] = -spe3;
SP[l] = 0; SP[k] = 0; SP[p] = 0; SP[lp] = -spe3; SP[r] = -spe3;
SP[l, k] = s/2; SP[k, p] = -t/2; SP[l, p] = -u/2; SP[k, r] = -a/2; SP[p, r] = b/2;
SP[l, lp] = -(a + b + t)/2; SP[l, r] = (s + t + u + a + b)/2; SP[k, lp] = (s + t + a)/2;
SP[lp, p] = -(t + u + b)/2; SP[lp, r] = (s + t + u)/2 + spe3;
quarkStrings = {GA[mu].(GSD[k] - GSD[r]).GAD[rho], GAD[rho].(GSD[p] + GSD[r]).GA[mu]};
quarkStringsBar = {GAD[sigma].(GSD[k] - GSD[r]).GA[nu], GA[nu].(GSD[p] + GSD[r]).GAD[sigma]};
leptonTrace = DiracTrace[GSD[lp].GA[mu].GSD[l].GA[5].GA[nu]]/2;
quarkTraces = Table[DiracTrace[GS[p].quarkStrings[[i]].GS[k].GA[5].quarkStringsBar[[j]]]/2, {i, 2}, {j, 2}];
interference = Table[Factor[Together[Simplify[Contract[ExpandScalarProduct[Contract[
  DiracSimplify[leptonTrace quarkTraces[[i, j]], DiracTraceEvaluate -> True] (-MTD[rho, sigma])]]]]]], {i, 2}, {j, 2}];
MatrixForm[interference]


(* ::Text:: *)
(*Multiply by the propagators of the two diagrams, 1/(a tau) and 1/(b tau), the colour factor CF (from Tr[T^a T^a]/Nc), and the couplings ee^4 eq^2 gs^2 (the generated prefactors, 3 eq/2 times the charge 2/3 times i's, square to one).*)


(* ::Input:: *)
realPolarized = Factor[Together[CF ee^4 eq^2 gs^2/(a + b + t)^2 (
  interference[[1, 1]]/a^2 + interference[[2, 2]]/b^2 + (interference[[1, 2]] + interference[[2, 1]])/(a b))]];
real = <|"Unpolarized" -> realUnpolarized, "Polarized" -> realPolarized|>;
Factor[realPolarized /. spe3 -> 0]


(* ::Section:: *)
(*Photon-quark Born: gamma q -> q g*)


(* ::Text:: *)
(*Needed through order eps for the photon-in-lepton subtraction.  Invariants sh = 2 q.k, th = -2 k.p, uh = -2 q.p.  Unpolarized: average 1/2 over the quark spin and 1/(D-2) over the photon polarizations, 1/Nc over colour.*)


(* ::Input:: *)
FCSetDiracGammaScheme["NDR"];
FCClearScalarProducts[];
SPD[q] = 0; SPD[k] = 0; SPD[p] = 0; SPD[r] = 0;
SPD[q, k] = sh/2; SPD[k, p] = -th/2; SPD[q, p] = -uh/2; SPD[p, r] = sh/2; SPD[k, r] = -uh/2; SPD[q, r] = -th/2;
photonAmplitude = Total[SUNSimplify[#, SUNNToCACF -> True] & /@ amplitudes["Photon"]] // DiracSimplify // Contract;
photonUnpolarized = photonAmplitude ComplexConjugate[photonAmplitude]/SUNN // FeynAmpDenominatorExplicit //
  SUNSimplify[#, Explicit -> True, SUNNToCACF -> True] & // FermionSpinSum[#, ExtraFactor -> 1/2] & //
  DiracSimplify // Contract // DoPolarizationSums[#, r, 0, VirtualBoson -> True] & //
  DoPolarizationSums[#, q, 0, VirtualBoson -> True] & // Contract // Simplify // SUNSimplify[#, SUNNToCACF -> True] & // Factor;
photonUnpolarized = photonUnpolarized/(D - 2) /. {SMP["e"] -> ee, SMP["g_s"] -> gs}


(* ::Text:: *)
(*Polarized: the circularly polarized photon has the density difference i eps^{mu nu q k}/(2 q.k) (with $LeviCivitaSign = -1) and the quark the density k-slash gamma5/2; all momenta of this 2 -> 2 process are four-dimensional, the final gluon is summed with -g in D dimensions.*)


(* ::Input:: *)
FCSetDiracGammaScheme["BMHV"];
FCClearScalarProducts[];
Scan[(Momentum[#, D - 4] = 0) &, {q, k, p, r}];
SPD[q] = 0; SPD[k] = 0; SPD[p] = 0; SPD[r] = 0;
SPD[q, k] = sh/2; SPD[k, p] = -th/2; SPD[q, p] = -uh/2; SPD[p, r] = sh/2; SPD[k, r] = -uh/2; SPD[q, r] = -th/2;
SP[q] = 0; SP[k] = 0; SP[p] = 0; SP[r] = 0;
SP[q, k] = sh/2; SP[k, p] = -th/2; SP[q, p] = -uh/2; SP[p, r] = sh/2; SP[k, r] = -uh/2; SP[q, r] = -th/2;
circularPhoton = Contract[I LC[mu, nu, alpha, beta] FV[q, alpha] FV[k, beta]/(2 SP[q, k])];
photonStrings = {GAD[rho].(GSD[k] + GSD[q]).GA[mu], GA[mu].(GSD[k] - GSD[r]).GAD[rho]};
photonStringsBar = {GA[nu].(GSD[k] + GSD[q]).GAD[sigma], GAD[sigma].(GSD[k] - GSD[r]).GA[nu]};
photonTraces = Table[DiracTrace[GS[p].photonStrings[[i]].GS[k].GA[5].photonStringsBar[[j]]]/2, {i, 2}, {j, 2}];
photonInterference = Table[Factor[TrickMandelstam[Cancel[Together[Simplify[Contract[ExpandScalarProduct[Contract[
  DiracSimplify[circularPhoton photonTraces[[i, j]], DiracTraceEvaluate -> True] (-MTD[rho, sigma])]]]]]], {sh, th, uh, 0}]], {i, 2}, {j, 2}];
photonPolarized = CF ee^2 eq^2 gs^2 Factor[TrickMandelstam[Cancel[Together[
  photonInterference[[1, 1]]/sh^2 + photonInterference[[2, 2]]/uh^2 + (photonInterference[[1, 2]] + photonInterference[[2, 1]])/(sh uh)]], {sh, th, uh, 0}]]


(* ::Input:: *)
photon = <|"Unpolarized" -> photonUnpolarized, "Polarized" -> photonPolarized|>;


(* ::Section:: *)
(*Save*)


(* ::Text:: *)
(*Algebra.wl holds everything step 3 needs: the Born and photon squares (functions of D), the real squares (functions of D, a, b, s, t, u and spe3) and the virtual scalar family.*)


(* ::Input:: *)
algebra = <|"Born" -> born, "Real" -> real, "Photon" -> photon, "VirtualFamily" -> virtualFamily|>;
SaveResult["Algebra", algebra, "FeynCalc"];
Keys[algebra]
