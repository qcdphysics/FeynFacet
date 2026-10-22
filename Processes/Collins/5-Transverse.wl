(* ::Title:: *)
(*Lepton-quark scattering with transverse spin transfer, step 5: the NLO transverse-spin-transfer hard coefficient*)


(* ::Text:: *)
(*Same structure as 4-Unpolarized.wl (and Processes/LeptonQuark/4-Unpolarized.wl, read that first), for the transverse spin transfer T of Section "NLO transverse-spin-transfer hard part" of the paper: the full difference sigma(k up -> p up) - sigma(k up -> p down), projected onto the Collins combination (T_XX + T_NN)/2 ("Transverse" entries of Algebra.wl and Reduced.wl).  Output: Results/<route>/Transverse.wl and the comparison with Eqs. (transverse-finite-hard), (transverse-endpoints) and (transverse-regular) of the paper (reference/Draft-transverse.wl).*)
(**)
(*What differs from the unpolarized case.  (1) The HVBM scheme: the hatted invariant spe3 of step 2 became integrals in D + 2 dimensions in step 3; the normal component of the unobserved gluon momentum was reduced in step 2.  (2) The collinear counterterms of the incoming and of the fragmenting quark both use the transversity kernel Delta_T P_qq(z) = CF [2 z/(1-z)_+ + 3/2 delta(1-z)] (the transverse spin is transferred along the quark line; there is no gluon transversity).  (3) Three factorization scales are kept apart, muI for the incoming quark, muO for the fragmenting quark and muL for the photon in the lepton, as in the paper; the photon in the lepton is unpolarized, P_photon/lepton(y) = (1 + (1-y)^2)/y.  (4) There is no finite helicity-restoring term (that term belongs to the longitudinal case).  (5) The virtual correction is the vector form factor times the transverse Born term.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
spin = "Transverse";
algebra = LoadResult["Algebra", route];
reduced = LoadResult["Reduced", route];
region = s > 0 && muR > 0 && muI > 0 && muO > 0 && muL > 0 && 0 < v < 1 && 0 < w < 1;
toVW = {t -> -s (1 - v), s4 -> s v (1 - w), u -> -s v w};
couplings = {ee -> Sqrt[4 Pi alphaEM], gs -> Sqrt[4 Pi alphaS]};
Keys[reduced]


(* ::Section:: *)
(*Born terms*)


(* ::Text:: *)
(*The transverse Born hard part, born0 + eps born1, multiplies the collinear counterterms; the Collins projection has no order-eps term (step 2).  Expected: born0 = 4 alphaEM^2 eq^2/(s (1-v)^2), the Born term of Eq. (transverse-finite-hard).*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "Draft-transverse.wl"}]];
paper = DraftTransverse[s, v, w, muI, muO, muL, alphaEM, alphaS, eq, CF];
bornHard = Series[algebra["Born"][spin]/(16 Pi^2 s v) /. D -> 4 - 2 eps /. couplings /. {t -> -s (1 - v), u -> -s v}, {eps, 0, 1}];
born0 = SeriesCoefficient[bornHard, {eps, 0, 0}] // Simplify;
born1 = SeriesCoefficient[bornHard, {eps, 0, 1}] // Simplify;
Require[Simplify[born0 - paper["Born"]] === 0 && born1 === 0, "Transverse Born hard part differs from the paper", {born0, born1}];
{born0, born1}


(* ::Text:: *)
(*Photon-quark Born with the transverse projection, normalized as in the unpolarized case: the photon carries the fraction y = (1-v)/(1-vw) of the lepton momentum, sh = y s, th = -s(1-v), uh = -y s v w, and photonHard = |M|^2/(16 Pi^2 y s (1 - v w)) in units of alphaS/Pi, through order eps.*)


(* ::Input:: *)
yPhoton = (1 - v)/(1 - v w);
photonHard = Series[(Pi/alphaS) algebra["Photon"][spin]/(16 Pi^2 yPhoton s (1 - v w)) /. D -> 4 - 2 eps /. couplings /.
  {sh -> yPhoton s, th -> -s (1 - v), uh -> -yPhoton s v w}, {eps, 0, 1}];
photon0 = Simplify[SeriesCoefficient[photonHard, {eps, 0, 0}], region];
photon1 = Simplify[SeriesCoefficient[photonHard, {eps, 0, 1}], region];
{photon0, photon1}


(* ::Section:: *)
(*Real emission and its threshold expansion*)


(* ::Text:: *)
(*As in the unpolarized case: the exact threshold factors s4^(-eps), zAB^(-1-eps) and zAB^(-eps) give the branches (1-w)^(-1-2eps) and (1-w)^(-1-eps); ThresholdSplit separates them and EndpointExpand turns each into plus distributions.  The scale muR^(2 eps) of the strong coupling is attached to the real emission; the three factorization scales enter through the counterterms.*)


(* ::Input:: *)
realHard = muR^(2 eps) reduced["Real"][spin]/(32 Pi^3) /. couplings /. toVW;
split = ThresholdSplit[realHard, w, eps];
realExpanded = EndpointExpand[split["Regular"], split["Branches"], eps, w, region];
Simplify[realExpanded["Pole2"], region]


(* ::Section:: *)
(*Virtual correction*)


(* ::Input:: *)
virtualHard = Series[alphaS CF/(2 Pi) (4 Pi muR^2)^eps reduced["Virtual"] (born0 + eps born1) /. t -> -s (1 - v), {eps, 0, 0}];
virtualHard = Map[Simplify[#, region] &, virtualHard];
virtualPole2 = SeriesCoefficient[virtualHard, {eps, 0, -2}];
virtualPole1 = SeriesCoefficient[virtualHard, {eps, 0, -1}];
virtualFinite = SeriesCoefficient[virtualHard, {eps, 0, 0}];
Require[Simplify[realExpanded["Pole2"] + virtualPole2, region] === 0, "Double poles do not cancel"];
"double poles cancel"


(* ::Section:: *)
(*Collinear counterterms and the single poles*)


(* ::Text:: *)
(*With the convolution formulas of the lepton-quark and lepton-off-diagonal processes, for a kernel {R(z), c0, cdelta} = R(z) + c0 [1/(1-z)]_+ + cdelta delta(1-z) in units of alphaS/(2 Pi):*)
(**)
(*  incoming quark (scale muI):   + (alphaS/(2 Pi)) (1/eps) b0 {(R(w) + c0)/w, c0, cdelta}   in the basis {regular, D0, delta},*)
(*  fragmenting quark (scale muO): + (alphaS/(2 Pi)) (1/eps) (1/v) {c0 (j - j(1))/(1-w) + v R(dw) j, c0 j(1), j(1) (cdelta + c0 log v)},  dw = 1 - v + v w,  j = v w b0(v -> v w/dw)/dw^2,*)
(*  photon in the lepton (scale muL): + (alphaEM/(2 Pi)) (1/eps) P_photon/lepton(y) (alphaS/Pi) photon0.*)
(**)
(*Here R(z) = -2 CF, c0 = 2 CF, cdelta = 3 CF/2 (the transversity kernel) on both quark legs.  The finite parts carry the MSbar constant, the logarithm of muR^2 over the respective factorization scale, the order-eps Born terms and, for the fragmentation, the order-eps Jacobian 2 j0 log(dw).  The table lists every contribution to the single pole; the sum must vanish in each column.*)


(* ::Input:: *)
incomingTuple[kernel_, bornTerm_] := bornTerm {(kernel[[1]] + kernel[[2]])/w, kernel[[2]], kernel[[3]]};
dw = 1 - v + v w;
jFunction[bornTerm_] := v w (bornTerm /. v -> v w/dw)/dw^2;
fragmentationTuple[kernel_, j_] := Module[{jEnd = Simplify[j /. w -> 1, region]},
  (1/v) {kernel[[2]] (j - jEnd)/(1 - w) + v kernel[[1]] j, kernel[[2]] jEnd, jEnd (kernel[[3]] + kernel[[2]] Log[v])}];
kernelIn = SplittingKernel["Transverse qq", w];
kernelOut = SplittingKernel["Transverse qq", dw];
photonKernel = SplittingKernel["photon/lepton", yPhoton][[1]];
pdfPole = alphaS/(2 Pi) incomingTuple[kernelIn, born0];
ffPole = alphaS/(2 Pi) fragmentationTuple[kernelOut, jFunction[born0]];
photonPole = {alphaEM/(2 Pi) photonKernel alphaS/Pi photon0, 0, 0};
poleTable = {realExpanded["Pole1"], {0, 0, virtualPole1}, pdfPole, ffPole, photonPole};
singlePoles = Total[poleTable];
Require[And @@ (ExactDifference[#, 0, region] === 0 & /@ singlePoles), "Single poles do not cancel", singlePoles];
Grid[Prepend[MapThread[Prepend[Map[Simplify[#, region] &, #2], #1] &,
  {{"real emission", "virtual", "incoming quark", "outgoing quark", "photon in lepton", "sum"}, Append[poleTable, singlePoles]}],
  {"1/eps coefficient of", "regular", "D0", "delta"}], Alignment -> Left, Dividers -> {{2 -> True}, {2 -> True, 7 -> True}}]


(* ::Section:: *)
(*The finite hard coefficient*)


(* ::Input:: *)
pdfFinite = alphaS/(2 Pi) ((MSbarConstant + Log[muR^2/muI^2]) incomingTuple[kernelIn, born0] + incomingTuple[kernelIn, born1]);
ffFinite = alphaS/(2 Pi) ((MSbarConstant + Log[muR^2/muO^2]) fragmentationTuple[kernelOut, jFunction[born0]] +
  fragmentationTuple[kernelOut, jFunction[born1] + 2 jFunction[born0] Log[dw]]);
photonFinite = {alphaEM/(2 Pi) photonKernel alphaS/Pi ((MSbarConstant + Log[muR^2/muL^2]) photon0 + photon1), 0, 0};
counterterms = {#[[1]], #[[2]], 0, #[[3]]} &[pdfFinite + ffFinite + photonFinite];
finite = realExpanded["Finite"] + {0, 0, 0, virtualFinite} + counterterms;
finiteTable = {realExpanded["Finite"], {0, 0, 0, virtualFinite}, {#[[1]], #[[2]], 0, #[[3]]} &[pdfFinite],
  {#[[1]], #[[2]], 0, #[[3]]} &[ffFinite], {#[[1]], #[[2]], 0, #[[3]]} &[photonFinite]};
finiteNames = {"real emission", "virtual", "incoming quark", "outgoing quark", "photon in lepton"};
finite = finite /. reduced["MasterFunctions"];
prefactor = CF alphaEM^2 eq^2 alphaS/(Pi s v w);
prefactorEndpoint = prefactor /. w -> 1;
coefficients = Map[Collect[Expand[#], {_Log, _PolyLog}, Factor] &, finite/{prefactor, prefactorEndpoint, prefactorEndpoint, prefactorEndpoint}];
Require[FreeQ[Simplify[D[finite, muR], region], muR] && Simplify[D[finite, muR], region] === {0, 0, 0, 0}, "The finite coefficient depends on the renormalization scale (no strong coupling in the Born term: it must not)"];
TableForm[Transpose[{ThresholdBasis, coefficients}]]


(* ::Text:: *)
(*The scale muR cancels (checked above): the Born process has no strong coupling, so no renormalization logarithm appears at this order.  The derivatives with respect to the three factorization scales are the convolutions of the splitting kernels with the Born terms: muI^2 d/dmuI^2 of the finite coefficient equals minus the incoming-quark pole tuple, and likewise for muO and muL.*)


(* ::Input:: *)
scaleDerivative[mu_] := mu/2 D[finite, mu];                                   (* mu^2 d/dmu^2 *)
scaleChecks = {Simplify[scaleDerivative[muI] + {#[[1]], #[[2]], 0, #[[3]]} &[pdfPole], region],
  Simplify[scaleDerivative[muO] + {#[[1]], #[[2]], 0, #[[3]]} &[ffPole], region],
  Simplify[scaleDerivative[muL] + {#[[1]], #[[2]], 0, #[[3]]} &[photonPole], region]};
Require[scaleChecks === {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, "A factorization-scale derivative is not the corresponding splitting convolution", scaleChecks];
"scale dependence verified"


(* ::Section:: *)
(*Comparison with the paper*)


(* ::Text:: *)
(*reference/Draft-transverse.wl transcribes the paper's result as written: its prefactor and its coefficients R_T, B_T, C_T, A_T of 1, D0, D1 and delta(1-w).  The paper's prefactor contains 1/w and multiplies the distributions as well, so CanonicalDistribution brings it to the basis used here.  The difference must vanish identically in the physical region; it is also evaluated numerically at three points.*)


(* ::Input:: *)
literature = CanonicalDistribution[paper["Prefactor"], paper["Coefficients"], w];
difference = Table[ExactDifference[finite[[i]], literature[[i]], region], {i, 4}];
numeric = NumericDifference[finite, literature, PhysicalPoints /. {muR -> 1, muI -> 1, muO -> 1, muL -> 1} /. {muI -> 7/10, muO -> 13/10, muL -> 9/10}, 20];
{difference, numeric}


(* ::Section:: *)
(*Save*)


(* ::Input:: *)
result = <|"Basis" -> ThresholdBasis, "Prefactor" -> prefactor, "Coefficients" -> coefficients, "Finite" -> finite,
  "Born" -> born0, "Virtual" -> virtualFinite, "Real" -> realExpanded["Finite"], "Counterterms" -> counterterms,
  "LiteratureDifference" -> difference, "Route" -> route, "Spin" -> spin|>;
SaveResult["Transverse", result, route];
Require[difference === {0, 0, 0, 0}, "Finite coefficient differs from the paper", difference];
coefficients[[4]]
