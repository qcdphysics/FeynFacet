(* ::Title:: *)
(*TMD beam and soft functions, step 3: phase space, rapidity regulator and master integrals*)


(* ::Text:: *)
(*This step does not depend on the algebra route.  It takes the structure of the squared matrix elements established in step 2 (real emission: gs^2 K(x,eps)/kT2; soft: 4 gs^2 C_i/((n.k)(nb.k))), performs the phase-space integration with the measurement of the TMD functions, applies the rapidity regulator, and evaluates the resulting transverse integrals.  Output: Results/Integrals.wl, used by steps 4 and 5.*)
(**)
(*Measurement.  The beam function at fixed momentum fraction x and transverse separation b is*)
(**)
(*   B^(1)(x, b) = Integral d^D k/(2 Pi)^(D-1) delta+(k^2) delta(1 - x - n.k/n.p) e^{i b.kT} R(k) M2(x, kT2),*)
(**)
(*with M2 from step 2 (it already contains the normalization 1/(2 n.p) or 1/(x (n.p)^2) of the operator) and the exponential rapidity regulator R(k) = exp[-tau e^{-EulerGamma} (n.k + nb.k)] of Li, Neill and Zhu [arXiv:1604.00392], nu = 1/tau, with the limit tau -> 0 taken before eps -> 0.  The soft function has the same form without the delta function in x.  In light-cone coordinates d^D k = (1/2) d(n.k) d(nb.k) d^{D-2}kT and delta+(k^2) fixes nb.k = kT2/(n.k), so that*)
(**)
(*   beam:  Integral d^D k/(2 Pi)^(D-1) delta+(k^2) delta(1 - x - n.k/n.p) ... = Integral d^{D-2}kT/((2 Pi)^(D-1) 2 (1-x)) ... at n.k = (1-x) n.p,*)
(*   soft:  Integral d^D k/(2 Pi)^(D-1) delta+(k^2) ...                      = Integral d^{D-2}kT/((2 Pi)^(D-1) 2) Integral_0^Infinity dy/y ...  with y = n.k.*)
(**)
(*The couplings are renormalized in the MSbar scheme: gs^2 -> 4 Pi alphaS muR^(2 eps) (e^EulerGamma/(4 Pi))^eps, and as = alphaS/(4 Pi).  Collecting the factors, every contribution is a_s times a transverse Fourier integral of the form*)
(**)
(*   T[f] = 4 Pi muR^(2 eps) (e^EulerGamma/(4 Pi))^eps Integral d^{D-2}kT/(2 Pi)^(D-2) e^{i b.kT} f(kT)/kT2,*)
(**)
(*namely (step 5)  B^(1) = as K(x,eps)/(1-x) T[R]  and (step 4)  S^(1) = 4 C_i as T[ Integral dy/y R ].  The results are expressed through Lp = Log[b^2 muR^2/b0^2] with b0 = 2 e^{-EulerGamma}, Lv = Log[nu^2/muR^2] and, for the beam function, LQ = 2 Log[x n.p/nu].*)
(**)
(*Two transverse integrals occur: the one with f = 1 (the regulator drops out away from x = 1), whose value is the Gamma function F_perp(eps, Lp) = e^{-EulerGamma eps} Gamma[-eps] e^{eps Lp}, and the one with the regulated longitudinal integral, f = Integral dy/y exp[-tau e^{-EulerGamma}(y + kT2/y)] = 2 K0(2 tau e^{-EulerGamma} kT), a Bessel function.  Both are written as one-dimensional Euler integrals, reduced with Kira to integrals that are finite at eps = 0, and those are evaluated with SubTropica; the closed forms serve as checks.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadSubTropica[]


(* ::Section:: *)
(*The transverse master F_perp*)


(* ::Text:: *)
(*The Fourier transform of 1/kT2 in D - 2 = 2 - 2 eps dimensions is a standard integral,*)
(**)
(*   Integral d^{2-2eps}kT/(2 Pi)^{2-2eps} e^{i b.kT} (kT2)^(-1) = Gamma[-eps] (Pi b^2)^eps/(4 Pi),*)
(**)
(*so that T[1] = e^{eps EulerGamma} Gamma[-eps] (muR^2 b^2/4)^eps = e^{-EulerGamma eps} Gamma[-eps] e^{eps Lp} = F_perp.  To stay within the same tool chain as the other processes we nevertheless write Gamma[-eps] Gamma[1+eps] as the Euler integral T0 = Integral_0^Infinity du u^(-1-eps)/(1+u) (a Beta function), derive the integration-by-parts identity by differentiating u^(-eps)/(1+u),*)
(**)
(*   0 = Integral du d/du[u^(-eps)/(1+u)] = -eps T0 - T1,   T1 = Integral du u^(-eps)/(1+u)^2,*)
(**)
(*solve it with Kira, and evaluate T1, which is finite at eps = 0, with SubTropica.  The surface terms vanish for -1 < eps < 0.  Kira's "user-defined system" takes the identities as a text file; the integrals are labelled by integers, large labels for the integrals to be reduced (100: T0, 200: E0 below) and small ones for the masters (1: T1, 2: E1, 3: E2).*)


(* ::Input:: *)
ibpCoefficients[weight_, basis_, u_] := Module[{unknowns, derivative, polynomial},
  unknowns = Array[cc, Length[basis]];
  derivative = D[weight/(1 + u), u]/weight;
  polynomial = Numerator[Together[derivative - unknowns.basis]];
  unknowns /. First[Solve[Thread[CoefficientList[polynomial, u] == 0], unknowns]]];
transverseIBP = ibpCoefficients[u^(-eps), {1/(u (1 + u)), 1/(1 + u)^2}, u]     (* coefficients of T0 and T1 in the identity *)


(* ::Text:: *)
(*The soft master below gives a second identity.  Its integrand carries the parameter c = 1 + b^2/(b0^2 tau^2) > 1 (derived in the soft section):*)
(**)
(*   E0 = Integral du u^(-1-eps) (1+c u)^eps/(1+u),   E1 = Integral du u^(-eps) (1+c u)^(eps-1)/(1+u),   E2 = Integral du u^(-eps) (1+c u)^eps/(1+u)^2,*)
(**)
(*and differentiating u^(-eps)(1+c u)^eps/(1+u) gives -eps E0 + eps c E1 - E2 = 0.*)


(* ::Input:: *)
softIBP = ibpCoefficients[u^(-eps) (1 + c u)^eps, {1/(u (1 + u)), 1/((1 + u) (1 + c u)), 1/(1 + u)^2}, u]


(* ::Input:: *)
kiraDir = PrepareKira[FileNameJoin[{NotebookDirectory[], "kira", "measured"}], FileNameJoin[{ResultsDirectory[], "kira", "measured"}]];
equation[labels_, coefficients_] := StringJoin[MapThread[ToString[#1] <> "*(" <> StringReplace[ToString[#2 /. eps -> ep, InputForm], " " -> ""] <> ")\n" &, {labels, coefficients}]];
equations = equation[{100, 1}, transverseIBP] <> "\n" <> equation[{200, 2, 3}, softIBP];
Export[FileNameJoin[{kiraDir, "equations.kira"}], equations, "Text"];
WriteKiraTargets[kiraDir, {100, 200}, {}];
equations


(* ::Input:: *)
RunKira[kiraDir];
kiraRules = KiraRules[kiraDir, "Tuserweight", {}] /. ep -> eps;
CheckReductionCoverage[Tuserweight /@ {100, 200}, kiraRules, {Tuserweight}, Tuserweight /@ {1, 2, 3}];
kiraRules


(* ::Text:: *)
(*As expected: T0 = -T1/eps and E0 = c E1 - E2/eps.  Now the finite integrals with SubTropica (EulerIntegral[prefactor, integrand, variables, parameters, order, eps] expands to the given order in eps).  T1 is needed to order eps^2: the soft master multiplies F_perp by a 1/eps pole, so F_perp is needed to order eps.*)


(* ::Input:: *)
T1 = EulerIntegral[1, u^(-eps)/(1 + u)^2, {u}, {}, 2, eps]


(* ::Input:: *)
T0 = Tuserweight[100] /. kiraRules /. Tuserweight[1] -> T1;
Fperp = Collect[Normal[Series[Exp[-EulerGamma eps + eps Lp] T0/Gamma[1 + eps], {eps, 0, 1}]], eps, Simplify]


(* ::Text:: *)
(*Check against the closed form e^{-EulerGamma eps} Gamma[-eps] e^{eps Lp}: the difference must vanish through the computed order.*)


(* ::Input:: *)
Require[Simplify[Normal[Series[Fperp - Exp[-EulerGamma eps + eps Lp] Gamma[-eps], {eps, 0, 1}]]] === 0, "F_perp differs from the Gamma function"];


(* ::Section:: *)
(*The rapidity regulator at the endpoint x -> 1*)


(* ::Text:: *)
(*The diagonal kernels behave like Gamma_i/(1-x) for x -> 1 (Gamma_i = 4 C_i, the cusp anomalous dimension).  There the regulator matters: with y = 1 - x, P = n.p and a = tau e^{-EulerGamma},*)
(**)
(*   R/(1-x) = e^{-a(y P + kT2/(y P))}/y  ->  [1/y]_+  +  delta(y) Integral_0^1 dy/y e^{-a(y P + kT2/(y P))}      (as a -> 0),*)
(**)
(*because for a test function f(y), Integral dy [f(y) - f(0)] e^{...}/y -> Integral dy [f(y) - f(0)]/y.  The coefficient of delta(y) is, up to terms that vanish with a,*)
(**)
(*   Integral_0^1 dy/y e^{-a(y P + kT2/(y P))} = Integral_0^Infinity - Integral_1^Infinity = 2 K0(2 a kT) + EulerGamma + Log[a P]   (Integral_1^Infinity dy/y e^{-a P y} = -EulerGamma - Log[a P] + O(a)),*)
(**)
(*in which the Bessel function is exactly the longitudinal integral of the soft function, and EulerGamma + Log[a P] = Log[tau P] = Log[P/nu].  So the endpoint of the bare beam function contains the bare soft function (this is the overlap removed by the zero-bin subtraction) plus Log[P/nu] = LQ/2 times the transverse master.  We verify the identity numerically.*)


(* ::Input:: *)
endpointCheck = Table[With[{a = 10^-m, P = 3, kT = 1},
    NIntegrate[Exp[-a (y P + kT^2/(y P))]/y, {y, 0, 1}, WorkingPrecision -> 20, PrecisionGoal -> 12] - N[2 BesselK[0, 2 a kT] + EulerGamma + Log[a P], 20]],
  {m, {3, 4, 5}}];
Require[Abs[endpointCheck[[3]]] < 10^-4 && Abs[endpointCheck[[3]]/endpointCheck[[1]]] < 10^-1, "The endpoint identity fails", endpointCheck];
endpointCheck      (* the differences are O(a): they shrink by a factor 10 with a *)


(* ::Section:: *)
(*The soft master: Bessel function in transverse momentum, Fourier transformed*)


(* ::Text:: *)
(*The soft function requires T[2 K0(2 a kT)].  With the integral representation K0(z) = (1/2) Integral_0^Infinity dt/t exp[-t - z^2/(4t)] and 1/kT2 = Integral_0^Infinity ds e^{-s kT2}, the kT integral is Gaussian, Integral d^{2-2eps}kT e^{i b.kT} e^{-sigma kT2} = (Pi/sigma)^{1-eps} e^{-b^2/(4 sigma)} with sigma = s + a^2/t; the s integral gives an incomplete Gamma function and the t integral turns it into*)
(**)
(*   T[2 K0(2 a kT)] = F_perp z^(-eps) E(z, eps),   E(z, eps) = Integral_0^1 du u^(-1-eps) (1 + z u)^eps,   z = b^2/(4 a^2) = b^2/(b0^2 tau^2),*)
(**)
(*and the substitution u -> u/(1+u) maps E onto the Euler integral E0 above with c = 1 + z.  We first check the representation numerically against the direct Bessel-Fourier integral (angular average of e^{i b.kT} in 2 - 2 eps dimensions: Gamma[1-eps] (2/(b kT))^(-eps) J_{-eps}(b kT)) at eps = -1/4.*)


(* ::Input:: *)
directSoft[epsValue_, bValue_, aValue_] := Module[{eps = epsValue, b = bValue, a = aValue},
  (* T[2 K0] without the MSbar factors: 4 Pi Integral d^{2-2eps}kT/(2 Pi)^{2-2eps} e^{i b.kT} 2 K0(2 a kT)/kT2 *)
  4 Pi/(2 Pi)^(2 - 2 eps) 2 Pi^(1 - eps)/Gamma[1 - eps] NIntegrate[
    kT^(1 - 2 eps) Gamma[1 - eps] (b kT/2)^eps BesselJ[-eps, b kT] 2 BesselK[0, 2 a kT]/kT^2, {kT, 0, Infinity},
    WorkingPrecision -> 20, PrecisionGoal -> 10, MaxRecursion -> 20]];
eulerSoft[epsValue_, bValue_, aValue_] := Module[{eps = epsValue, b = bValue, a = aValue, z},
  z = b^2/(4 a^2);
  Gamma[-eps] (Pi b^2)^eps z^(-eps) NIntegrate[u^(-1 - eps) (1 + z u)^eps, {u, 0, 1}, WorkingPrecision -> 20, PrecisionGoal -> 10]];
softRepresentationCheck = {directSoft[-1/4, 2, 1/3], eulerSoft[-1/4, 2, 1/3]};
Require[Abs[softRepresentationCheck[[1]]/softRepresentationCheck[[2]] - 1] < 10^-8, "The Euler representation of the soft master is wrong", softRepresentationCheck];
softRepresentationCheck


(* ::Text:: *)
(*Now the finite masters E1 and E2 with SubTropica, as functions of c (the result contains logarithms and dilogarithms of c), and E0 from the Kira identity.  E0 has a 1/eps pole (the u -> 0 endpoint) and F_perp another, so for the finite part of the soft function E0 is needed to order eps, i.e. E1 to order eps and E2 to order eps^2.*)


(* ::Input:: *)
E1 = EulerIntegral[1, u^(-eps) (1 + c u)^(eps - 1)/(1 + u), {u}, {c}, 1, eps];
E2 = EulerIntegral[1, u^(-eps) (1 + c u)^eps/(1 + u)^2, {u}, {c}, 2, eps];
E0 = Tuserweight[200] /. kiraRules /. {Tuserweight[2] -> E1, Tuserweight[3] -> E2};
E0 = Collect[Normal[Series[E0, {eps, 0, 1}]], eps, Simplify]


(* ::Input:: *)
softMasterExact = Collect[Normal[Series[Fperp (c - 1)^(-eps) E0, {eps, 0, 0}]], eps, Simplify]


(* ::Text:: *)
(*Check of the reconstructed E0 against the direct Euler integral at c = 7: the truncation error of the series is of order eps^2, so halving eps must reduce the difference by about a factor of four.*)


(* ::Input:: *)
e0Direct[epsValue_, cValue_] := NIntegrate[u^(-1 - epsValue) (1 + cValue u)^epsValue/(1 + u), {u, 0, Infinity}, WorkingPrecision -> 25, PrecisionGoal -> 15];
e0Check = Table[(E0 /. {eps -> epsValue, c -> 7}) - e0Direct[epsValue, 7], {epsValue, {-1/5, -1/10}}];
Require[Abs[e0Check[[2]]/e0Check[[1]]] < 1/3, "E0 disagrees with its direct evaluation", e0Check];
e0Check


(* ::Text:: *)
(*Finally the limit tau -> 0, i.e. c = 1 + b^2/(b0^2 tau^2) -> Infinity at fixed b, keeping the logarithms of c and dropping the power corrections O(1/c): Log[c - 1] -> Log[c] with Log[c] = Log[b^2/(b0^2 tau^2)] = Lp + Lv.  The expansion of the dilogarithms at infinity is done by Mathematica.  The result is the soft master*)
(**)
(*   T[2 K0] -> F_perp (Lp + Lv - EulerGamma - PolyGamma[0, -eps]) + O(tau),*)
(**)
(*which is also what one gets by expanding 2 K0(2 a kT) = -Log[a^2 kT2] - 2 EulerGamma + O(a^2) under the integral and using T[Log kT2] = -d/d alpha T[(kT2)^(-alpha)] at alpha = 0; this closed form is the check.*)


(* ::Input:: *)
asymptotic[expr_] := Normal[Series[expr /. Log[c - 1] -> Log[c] /. Log[-1 + c] -> Log[c], {c, Infinity, 0}, Assumptions -> c > 1]] /. Log[c] -> Lp + Lv;
softMaster = Collect[asymptotic[softMasterExact], eps, Simplify]


(* ::Input:: *)
softMasterClosedForm = Normal[Series[Exp[-EulerGamma eps + eps Lp] Gamma[-eps] (Lp + Lv - EulerGamma - PolyGamma[0, -eps]), {eps, 0, 0}]];
Require[Simplify[softMaster - softMasterClosedForm] === 0, "The soft master differs from its closed form", {softMaster, softMasterClosedForm}];


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
integrals = <|"Fperp" -> Fperp, "SoftMaster" -> softMaster, "SoftMasterExact" -> softMasterExact,
  "Endpoint" -> "R/(1-x) -> [1/(1-x)]_+ + delta(1-x) (T[2K0] + F_perp LQ/2), see step 3",
  "KiraRules" -> kiraRules, "T1" -> T1, "E1" -> E1, "E2" -> E2|>;
SaveResult["Integrals", integrals];
Keys[integrals]
