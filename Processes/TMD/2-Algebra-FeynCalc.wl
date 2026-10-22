(* ::Title:: *)
(*TMD beam and soft functions, step 2 (FeynCalc route): squared matrix elements*)


(* ::Text:: *)
(*Input: the amplitudes of step 1 (Results/Amplitudes.wl).  Output: Results/FeynCalc/Algebra.wl with the squared operator matrix elements, summed over spins and colours, as functions of the momentum fraction x and the transverse momentum kT2 of the emitted parton.  The FORM and Tracer routes (2-Algebra-FORM.wl, 2-Algebra-Tracer.wl) produce the same file from the same amplitudes.*)
(**)
(*What is squared.  An operator matrix element <X| O |p> is the open amplitude of step 1 with the operator leg "O" carrying the operator's open index.  Its square, summed over the spins and colours of X and averaged over those of p, is closed with the operator's own structure:*)
(**)
(*   quark operator:  M2_q = 1/(2 n.p) * (1/2)(1/Nc) Sum Tr[ n-slash . Gamma . rho(p) . Gammabar ]       rho(p) = p-slash (incoming quark) or the gluon polarization sum,*)
(*   gluon operator:  M2_g = 1/(x (n.p)^2) * (1/(D-2))(1/(Nc^2-1)) Sum (-g_{mu nu}) M^mu M^nu*,*)
(**)
(*where n-slash is the gamma^+ of the operator qbar gamma^+ q, -g_{mu nu} contracts the two F^{+ mu} of the gluon operator, and the prefactors 1/(2 n.p) and 1/(x (n.p)^2) come from the Fourier transform in b^- that defines the beam function (docs/CONVENTIONS.md).  With these normalizations the Born matrix elements are exactly 1, which is the statement B = delta(1-x) at lowest order; the real-emission matrix elements, integrated over the one-parton phase space with the longitudinal constraint n.k = (1-x) n.p (step 5), give the NLO beam functions.*)
(**)
(*Kinematics.  n^2 = nb^2 = 0, n.nb = 2; p^2 = k^2 = 0; n.p = pn; n.k = (1-x) pn; and k_T, the transverse momentum of k with respect to p and n, satisfies kT2 = (n.k)(nb.k) and p.k = kT2/(2(1-x)).  External gluons are summed with the physical polarization sum d^{mu nu}(k) = -g^{mu nu} + (k^mu n^nu + n^mu k^nu)/(n.k) (light-cone gauge, reference vector n); the Ward identities of step 1 guarantee that the result does not depend on this choice, which is checked below for one channel.  The soft matrix element is a function of n.k and nb.k only.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
LoadFeynCalc[];
amplitudes = LoadAmplitudes[]["FeynCalc"];
Keys[amplitudes]


(* ::Section:: *)
(*Kinematics and polarization sums*)


(* ::Input:: *)
FCClearScalarProducts[];
SPD[p] = 0; SPD[k] = 0; SPD[n] = 0; SPD[nb] = 0; SPD[n, nb] = 2;
SPD[p, n] = pn; SPD[k, n] = (1 - x) pn; SPD[p, k] = kT2/(2 (1 - x)); SPD[k, nb] = kT2/((1 - x) pn); SPD[p, nb] = 0;
eikonal = FAEikonal[mom_, v_] :> 1/ExpandScalarProduct[SPD[v, mom]];
(* physical polarization sum of a gluon with momentum l, reference vector n *)
polarizationSum[l_, mu_, nu_] := -MTD[mu, nu] + (FVD[l, mu] FVD[n, nu] + FVD[n, mu] FVD[l, nu])/SPD[l, n];
polarizationSum[k, mu, nu]


(* ::Text:: *)
(*square[amplitude, conjugate, ...] multiplies an amplitude by the conjugate of another, replaces every polarization vector by an index and inserts the polarization sums, closes the Dirac chains with the spin sums, contracts the specified polarization tensors before expanding the traces, and performs the colour sum. The polarization tensors and their dimensional prescription are unchanged; only the order of the contractions is different.  The operator leg is treated as follows: for the quark operator the placeholder spinor ubar(q) ... u(q) of the two factors is replaced by n-slash (operatorSpinSum); for the gluon operator the placeholder polarization vectors epsilon(q) epsilon*(q) are replaced by -g_{mu nu}.*)


(* ::Input:: *)
(* the operator placeholder spinor ubar(q) is renamed, so that the slash that
   FermionSpinSum produces for it can be recognized and replaced by n-slash;
   the external massless fermions get their momentum-slash as usual *)
operatorQuark = Spinor[Momentum[q, D], 0, 1] -> Spinor[Momentum[qOperator, D], 0, 1];
operatorSpinSum = DiracGamma[Momentum[qOperator, D], D] :> GSD[n];
square[amp_, conj_, polarizationRules_, sums_:1] := Module[{expr},
  expr = (Total[amp] /. operatorQuark) ComplexConjugate[Total[conj] /. operatorQuark] /. eikonal;
  expr = expr /. polarizationRules // FeynAmpDenominatorExplicit;
  expr = FermionSpinSum[expr] /. operatorSpinSum;
  expr = Contract[expr sums] // ExpandScalarProduct;
  expr = DiracSimplify[expr, DiracTraceEvaluate -> True] // Contract // ExpandScalarProduct;
  expr = SUNSimplify[expr, Explicit -> True, SUNNToCACF -> False] // Contract;
  Factor[Together[expr /. SUNN -> Nc]]];
(* the colour sums are written through CF = (Nc^2-1)/(2 Nc) and CA = Nc at the
   end; the gluon normalization Tr(T^a T^b) = TF delta^{ab} with TF = 1/2 is FeynCalc's *)
colourFactors[expr_] := Factor[expr] /. (-1 + Nc) (1 + Nc) -> 2 Nc CF /. Nc -> CA;


(* ::Section:: *)
(*Born matrix elements: normalization*)


(* ::Text:: *)
(*Quark: Tr[n-slash p-slash]/(2 n.p) with the spin average 1/2 and the colour average 1/Nc of delta_ii = Nc.  Gluon: the amplitude is gOp (n.p eps^mu - p^mu n.eps); contracted with -g_{mu nu}, summed over the physical polarizations (which gives (D-2) (n.p)^2) and averaged with 1/((D-2)(Nc^2-1)), divided by (n.p)^2 (x = 1).*)


(* ::Input:: *)
bornQ = square[amplitudes["BornQ"], amplitudes["BornQ"], {}] /. x -> 1;
bornQ = Factor[bornQ (1/2) (1/Nc)/(2 pn)/gOp^2];
bornG = square[amplitudes["BornG"], amplitudes["BornG"],
  {Momentum[Polarization[q, -I], D] -> LorentzIndex[mu, D], Momentum[Polarization[q, I], D] -> LorentzIndex[nu, D],
   Momentum[Polarization[p, I], D] -> LorentzIndex[al, D], Momentum[Polarization[p, -I], D] -> LorentzIndex[be, D]}];
bornG = Factor[Contract[bornG (-MTD[mu, nu]) polarizationSum[p, al, be]] /(D - 2)/(Nc^2 - 1)/pn^2/gOp^2];
Require[bornQ === 1 && bornG === 1, "The Born matrix elements are not normalized to 1", {bornQ, bornG}];
{bornQ, bornG}


(* ::Section:: *)
(*Real emission: the four beam channels*)


(* ::Text:: *)
(*q -> q (gluon emitted): incoming quark averaged with 1/2 and 1/Nc, the emitted gluon summed with the physical polarization sum.  The chains are closed with n-slash (operator) and p-slash (incoming quark).*)


(* ::Input:: *)
gluonK = {Momentum[Polarization[k, -I], D] -> LorentzIndex[mu, D], Momentum[Polarization[k, I], D] -> LorentzIndex[nu, D]};
gluonP = {Momentum[Polarization[p, I], D] -> LorentzIndex[al, D], Momentum[Polarization[p, -I], D] -> LorentzIndex[be, D]};
operatorG = {Momentum[Polarization[q, -I], D] -> LorentzIndex[rho, D], Momentum[Polarization[q, I], D] -> LorentzIndex[sigma, D]};

realQQ = square[amplitudes["BeamQQ"], amplitudes["BeamQQ"], gluonK, polarizationSum[k, mu, nu]] (1/2) (1/Nc)/(2 pn)/gOp^2 // Factor


(* ::Text:: *)
(*g -> q (the antiquark is emitted): incoming gluon averaged with 1/(D-2) and 1/(Nc^2-1), summed with the physical polarization sum; the outgoing antiquark's spin sum is k-slash.*)


(* ::Input:: *)
realQG = square[amplitudes["BeamQG"], amplitudes["BeamQG"], gluonP, polarizationSum[p, al, be]]/(D - 2)/(Nc^2 - 1)/(2 pn)/gOp^2 // Factor


(* ::Text:: *)
(*q -> g (the quark is emitted; gluon operator): incoming quark averaged with 1/2 and 1/Nc; the operator indices are contracted with -g_{rho sigma}; the prefactor is 1/(x (n.p)^2).*)


(* ::Input:: *)
realGQ = square[amplitudes["BeamGQ"], amplitudes["BeamGQ"], operatorG, -MTD[rho, sigma]] (1/2) (1/Nc)/(x pn^2)/gOp^2 // Factor


(* ::Text:: *)
(*g -> g: both gluons with physical polarization sums, the incoming one averaged with 1/((D-2)(Nc^2-1)).*)


(* ::Input:: *)
realGG = square[amplitudes["BeamGG"], amplitudes["BeamGG"], Join[gluonK, gluonP, operatorG],
  polarizationSum[k, mu, nu] polarizationSum[p, al, be] (-MTD[rho, sigma])]/(D - 2)/(Nc^2 - 1)/(x pn^2)/gOp^2 // Factor


(* ::Text:: *)
(*All four are proportional to gs^2/kT2 times a function of x and D only: the dependence on n.p has cancelled, as it must by boost invariance along the beam.  Written as gs^2 K_ij(x, D)/kT2, the K_ij are the D-dimensional real-emission splitting functions times (1-x):*)
(**)
(*   K_qq/(1-x) = 2 CF [(1+x^2)/(1-x) - eps (1-x)],   K_qg/(1-x) = 2 TF [1 - 2x(1-x)/(1-eps)],   K_gq/(1-x) = 2 CF [(1+(1-x)^2)/x - eps x],   K_gg/(1-x) = 4 CA [x/(1-x) + (1-x)/x + x(1-x)],*)
(**)
(*with D = 4 - 2 eps and TF = 1/2.  At eps = 0 they are the LO splitting functions in the normalization alpha_s/(4 Pi) (twice the usual ones).*)


(* ::Input:: *)
real = colourFactors /@ <|"qq" -> realQQ, "qg" -> realQG, "gq" -> realGQ, "gg" -> realGG|>;
Require[FreeQ[real, pn], "The real matrix elements depend on n.p", real];
kernels = Map[Factor[# kT2/gs^2] &, real];
Require[FreeQ[kernels, kT2], "The real matrix elements are not proportional to 1/kT2"];
expected = <|"qq" -> 2 CF ((1 + x^2)/(1 - x) - eps (1 - x)), "qg" -> 2 (1/2) (1 - 2 x (1 - x)/(1 - eps)),
  "gq" -> 2 CF ((1 + (1 - x)^2)/x - eps x), "gg" -> 4 CA (x/(1 - x) + (1 - x)/x + x (1 - x))|>;
Require[And @@ Table[Simplify[kernels[c]/(1 - x) /. D -> 4 - 2 eps] === Simplify[expected[c]], {c, Keys[kernels]}],
  "The real-emission kernels differ from the D-dimensional splitting functions", kernels];
TableForm[Normal[kernels], TableHeadings -> {None, {"channel", "kT2 M2/gs^2"}}]


(* ::Text:: *)
(*Check of gauge invariance: the qq matrix element computed with -g^{mu nu} as polarization sum of the emitted gluon (allowed for a single external gluon when the Ward identity holds) is the same.*)


(* ::Input:: *)
realQQFeynman = square[amplitudes["BeamQQ"], amplitudes["BeamQQ"], gluonK, -MTD[mu, nu]] (1/2) (1/Nc)/(2 pn)/gOp^2 // Factor;
Require[Simplify[realQQFeynman - realQQ] === 0, "The qq matrix element depends on the gluon polarization sum", {realQQ, realQQFeynman}];


(* ::Section:: *)
(*The soft matrix element*)


(* ::Text:: *)
(*The soft current J^mu = gs T^a [n^mu/(n.k) - nb^mu/(nb.k)] (fundamental) or the same with the adjoint generator.  Its square with -g_{mu nu} (one external gluon), the colour trace Tr(T^a T^a) = CF Nc divided by Nc for the normalization Tr(1)/Nc = 1 of the soft function, gives 4 gs^2 CF/((n.k)(nb.k)); the adjoint line gives CA instead of CF.  Here the invariants n.k = kn and nb.k = kb are kept separate (the soft function has no longitudinal constraint).*)


(* ::Input:: *)
FCClearScalarProducts[];
SPD[k] = 0; SPD[n] = 0; SPD[nb] = 0; SPD[n, nb] = 2; SPD[k, n] = kn; SPD[k, nb] = kb;
softQ = square[amplitudes["SoftQ"], amplitudes["SoftQ"], gluonK, -MTD[mu, nu]]/Nc/gSoft^2 // Factor;
softG = square[amplitudes["SoftG"], amplitudes["SoftG"], gluonK, -MTD[mu, nu]]/(Nc^2 - 1)/gSoft^2 // Factor;
soft = colourFactors /@ <|"q" -> softQ, "g" -> softG|>


(* ::Input:: *)
Require[Simplify[soft["q"] - 4 gs^2 CF/(kn kb)] === 0 && Simplify[soft["g"] - 4 gs^2 CA/(kn kb)] === 0, "Soft matrix elements differ from 4 gs^2 C_i/((n.k)(nb.k))", soft];


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
algebra = <|"Born" -> <|"q" -> bornQ, "g" -> bornG|>, "Real" -> real, "Soft" -> soft|>;
SaveResult["Algebra", algebra, "FeynCalc"];
Keys[algebra]
