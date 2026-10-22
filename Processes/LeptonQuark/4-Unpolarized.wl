(* ::Title:: *)
(*Lepton-quark scattering, step 4: the unpolarized NLO hard coefficient*)


(* ::Text:: *)
(*Input: Results/<route>/Algebra.wl (step 2) and Results/<route>/Reduced.wl (step 3).  Output: Results/<route>/Unpolarized.wl with the finite coefficient in the basis {regular, D0, D1, delta}, and its comparison with Hinderer, Schlegel and Vogelsang, arXiv:1505.06415.*)
(**)
(*Hard coefficient.  Following that paper we use sigma(s,v,w) = s E dsigma/d^3p for the partonic cross section, with t = -s(1-v), u = -s v w, s4 = s v (1-w).  Hadronic cross sections follow by convoluting sigma with the parton distribution and fragmentation function (see docs/CONVENTIONS.md).  At Born level the recoiling lepton fixes w = 1:*)
(**)
(*   sigma_Born = |M|^2 / (16 Pi^2 s v) delta(1-w),*)
(**)
(*where |M|^2 is the spin- and colour-averaged squared matrix element.  Real emission with the two-body phase space Phi2 of the unobserved pair gives sigma_real = (integral of |M|^2 over Phi2)/(32 Pi^3).  The couplings carry the scale muR^(2 eps) in D dimensions (MSbar), and the four-dimensional measure of the observed parton is common to every term.*)
(**)
(*Steps: (1) Born and photon-quark Born through order eps; (2) real emission: separate the threshold branch (1-w)^(-1-2eps) and expand into plus distributions; (3) virtual correction; (4) check that the double and single poles cancel against the MSbar counterterms of the incoming quark, the outgoing quark and the photon in the lepton; (5) collect the finite coefficient and compare with the literature.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
algebra = LoadResult["Algebra", route];
reduced = LoadResult["Reduced", route];
region = s > 0 && muR > 0 && 0 < v < 1 && 0 < w < 1;          (* physical region; fixes every logarithm's branch *)
toVW = {t -> -s (1 - v), s4 -> s v (1 - w), u -> -s v w};
couplings = {ee -> Sqrt[4 Pi alphaEM], gs -> Sqrt[4 Pi alphaS]};
Keys[reduced]


(* ::Section:: *)
(*Born terms*)


(* ::Text:: *)
(*The D-dimensional Born hard part, bornHard = born0 + eps born1 + ..., multiplies the collinear counterterms, so its order-eps term contributes to the finite result.  At Born level u = -s v (w = 1).*)


(* ::Input:: *)
bornHard = Series[algebra["Born"]["Unpolarized"]/(16 Pi^2 s v) /. D -> 4 - 2 eps /. couplings /. {t -> -s (1 - v), u -> -s v}, {eps, 0, 1}];
born0 = SeriesCoefficient[bornHard, {eps, 0, 0}] // Simplify;
born1 = SeriesCoefficient[bornHard, {eps, 0, 1}] // Simplify;
Require[Simplify[born0 - 2 alphaEM^2 eq^2 (1 + v^2)/(s v (1 - v)^2)] === 0, "Born hard part differs from Eq. (22) of HSV"];
{born0, born1}


(* ::Text:: *)
(*Photon-quark Born gamma(q) + q(k) -> q(p) + g(r).  The photon carries the fraction y of the lepton momentum, so that sh = y s, th = -s(1-v), uh = -y s v w, and its on-shell condition sh + th + uh = 0 fixes y = (1-v)/(1-vw).  Its hard part photonHard = photon0 + eps photon1 is normalized like sigma_Born, with the Jacobian for the photon fraction: photonHard = |M|^2/(16 Pi^2 y s (1 - v w)) in units of alphaS/Pi.*)


(* ::Input:: *)
yPhoton = (1 - v)/(1 - v w);
photonHard = Series[(Pi/alphaS) algebra["Photon"]["Unpolarized"]/(16 Pi^2 yPhoton s (1 - v w)) /. D -> 4 - 2 eps /. couplings /.
  {sh -> yPhoton s, th -> -s (1 - v), uh -> -yPhoton s v w}, {eps, 0, 1}];
photon0 = Simplify[SeriesCoefficient[photonHard, {eps, 0, 0}], region];
photon1 = Simplify[SeriesCoefficient[photonHard, {eps, 0, 1}], region];
Require[Simplify[photon0 - 2 Pi CF alphaEM eq^2 (1 + v^2 w^2)/(s (1 - v) v w)] === 0, "Photon-quark Born differs from the expected form"];
{photon0, photon1}


(* ::Section:: *)
(*Real emission and its threshold expansion*)


(* ::Text:: *)
(*The integrated real contribution from step 3, in the variables v, w.  It contains the exact factors s4^(-eps) = (s v (1-w))^(-eps) and, from the AB family, zAB^(-1-eps) and zAB^(-eps) with zAB = (1-v)(1-w)/(1-v+vw).  Together they give the threshold branches (1-w)^(-1-2eps) and (1-w)^(-1-eps) that produce the plus distributions.  ThresholdSplit makes every power of (1-w) explicit (all other factors are positive in the physical region) and writes the result as*)
(**)
(*   sigma_real = A(w, eps) + Sum over lambda of (1-w)^(-1-lambda eps) B_lambda(w, eps),*)
(**)
(*with A and B_lambda regular at w = 1.*)


(* ::Input:: *)
realHard = muR^(2 eps) reduced["Real"]["Unpolarized"]/(32 Pi^3) /. couplings /. toVW;
split = ThresholdSplit[realHard, w, eps];
Keys[split["Branches"]]


(* ::Text:: *)
(*EndpointExpand applies (1-w)^(-1-lambda eps) = -delta(1-w)/(lambda eps) + D0 - lambda eps D1 + ... to every branch, with D0 = [1/(1-w)]_+ and D1 = [log(1-w)/(1-w)]_+, and splits w-dependent coefficients into their value at w = 1 times the distribution plus an ordinary remainder.  The result has the double pole, the single poles in the basis {regular, D0, delta} and the finite part in {regular, D0, D1, delta}.*)


(* ::Input:: *)
realExpanded = EndpointExpand[split["Regular"], split["Branches"], eps, w, region];
Simplify[realExpanded["Pole2"], region]


(* ::Section:: *)
(*Virtual correction*)


(* ::Text:: *)
(*The one-loop vertex correction is the quark form factor times the Born term: 2 Re(M1 M0*) / |M0|^2 = (alphaS CF/(2 Pi)) (4 Pi muR^2)^eps x virtual, where virtual is the reduced expression of step 3 with the massless bubble inserted, (-t)^(-eps) Gamma(eps) Gamma(1-eps)^2/Gamma(2-2eps).  Its Laurent expansion is -2/eps^2 - 3/eps - 8 + O(eps) times (4 Pi muR^2/(-t))^eps Gamma(1+eps) Gamma(1-eps)^2/Gamma(1-2eps), the standard result.  The external-leg graphs are scaleless and vanish; in the MSbar scheme their UV pole is cancelled by the wave-function counterterm and their IR pole is part of the collinear subtraction below.*)


(* ::Input:: *)
virtualHard = Series[alphaS CF/(2 Pi) (4 Pi muR^2)^eps reduced["Virtual"] (born0 + eps born1) /. t -> -s (1 - v), {eps, 0, 0}];
virtualHard = Map[Simplify[#, region] &, virtualHard];
virtualPole2 = SeriesCoefficient[virtualHard, {eps, 0, -2}];
virtualPole1 = SeriesCoefficient[virtualHard, {eps, 0, -1}];
virtualFinite = SeriesCoefficient[virtualHard, {eps, 0, 0}];
Require[Simplify[virtualPole2 + 2 alphaS CF/(2 Pi) born0] === 0, "Virtual double pole is not -2 times the Born term"];
virtualPole2


(* ::Section:: *)
(*Pole cancellation*)


(* ::Text:: *)
(*The double poles of real and virtual terms (both proportional to delta(1-w)) must cancel exactly.*)


(* ::Input:: *)
Require[Simplify[realExpanded["Pole2"] + virtualPole2, region] === 0, "Double poles do not cancel"];
"double poles cancel"


(* ::Text:: *)
(*The remaining single poles are removed by the MSbar collinear counterterms; the table at the end of this cell lists every contribution separately, and their sum must vanish in each column.  With the hadronic convolution measure dx/x (docs/CONVENTIONS.md) the counterterm for the incoming quark is*)
(**)
(*   + (alphaS/(2 Pi)) (1/eps) [(4 Pi)^eps / Gamma(1-eps)] integral dz/z Pqq(z) sigma_Born(z s, z t, u),*)
(**)
(*and the Born delta function fixes z = w, so the convolution becomes Pqq(w)/w times the Born hard part at (s, v).  Pqq(w)/w = CF [(1-w)/w + 2/(1-w)]_+ with the delta coefficient 3/2.  The outgoing quark fragments with momentum fraction z: sigma_Born(s, t/z, u/z) has v -> vF = v w/(1 - v + v w) and its delta function fixes z = 1 - v + v w, with the Jacobian v w/(1 - v + v w)^2 for the convolution variable; the plus prescription at z = 1 maps to w = 1 but the lower limit of z produces the log v in the delta coefficient.  The photon in the lepton gives + (alphaEM/(2 Pi)) (1/eps) P_photon/lepton(y) times the photon-quark hard part at y = (1-v)/(1-vw).  Each counterterm multiplies the D-dimensional Born term, whose order-eps part, together with the MSbar constant log(4 Pi) - EulerGamma, feeds the finite part.*)


(* ::Input:: *)
incomingTuple[bornTerm_] := CF bornTerm {(1 - w)/w, 2, 3/2};                  (* {regular, D0, delta} *)
dw = 1 - v + v w;
jFunction[bornTerm_] := v w (bornTerm /. v -> v w/dw)/dw^2;              (* Born at the shifted point, with the Jacobian *)
fragmentation[j_] := Module[{jEnd = Simplify[j /. w -> 1, region]},
  CF/v {2 (j - jEnd)/(1 - w) - v j (1 + dw), 2 jEnd, jEnd (2 Log[v] + 3/2)}];
photonKernel = (1 + (1 - yPhoton)^2)/yPhoton;                               (* P_photon/lepton *)
pdfPole = alphaS/(2 Pi) incomingTuple[born0];
ffPole = alphaS/(2 Pi) fragmentation[jFunction[born0]];
photonPole = {alphaEM/(2 Pi) photonKernel alphaS/Pi photon0, 0, 0};
poleTable = {realExpanded["Pole1"], {0, 0, virtualPole1}, pdfPole, ffPole, photonPole};
singlePoles = Total[poleTable];
Require[And @@ (ExactDifference[#, 0, region] === 0 & /@ singlePoles), "Single poles do not cancel", singlePoles];
Grid[Prepend[MapThread[Prepend[Map[Simplify[#, region] &, #2], #1] &,
  {{"real emission", "virtual", "incoming quark", "outgoing quark", "photon in lepton", "sum"}, Append[poleTable, singlePoles]}],
  {"1/eps coefficient of", "regular", "D0", "delta"}], Alignment -> Left, Dividers -> {{2 -> True}, {2 -> True, 7 -> True}}]


(* ::Section:: *)
(*The finite hard coefficient*)


(* ::Text:: *)
(*The regular parts still contain the inert symbols F1[n][zAB] etc. of the hypergeometric functions; they are replaced by their explicit expressions (logarithms and dilogarithms) here.  Finite parts of the counterterms: the MSbar constant times the pole coefficient, plus the order-eps Born terms (for the fragmentation also the order-eps Jacobian 2 j0 log(1-v+vw) from z^(-2eps)).  The finite coefficient is the sum of real, virtual and counterterm pieces in the basis {regular, D0, D1, delta}.*)


(* ::Input:: *)
pdfFinite = alphaS/(2 Pi) (MSbarConstant incomingTuple[born0] + incomingTuple[born1]);
ffFinite = alphaS/(2 Pi) (MSbarConstant fragmentation[jFunction[born0]] + fragmentation[jFunction[born1] + 2 jFunction[born0] Log[dw]]);
photonFinite = {alphaEM/(2 Pi) photonKernel alphaS/Pi (MSbarConstant photon0 + photon1), 0, 0};
counterterms = {#[[1]], #[[2]], 0, #[[3]]} &[pdfFinite + ffFinite + photonFinite];
finite = realExpanded["Finite"] + {0, 0, 0, virtualFinite} + counterterms;
finiteTable = {realExpanded["Finite"], {0, 0, 0, virtualFinite}, {#[[1]], #[[2]], 0, #[[3]]} &[pdfFinite],
  {#[[1]], #[[2]], 0, #[[3]]} &[ffFinite], {#[[1]], #[[2]], 0, #[[3]]} &[photonFinite]};
finiteNames = {"real emission", "virtual", "incoming quark", "outgoing quark", "photon in lepton"};
finite = finite /. reduced["MasterFunctions"];                 (* insert the explicit F1, F2 functions *)
prefactor = CF alphaEM^2 eq^2 alphaS/(Pi s v w);
prefactorEndpoint = prefactor /. w -> 1;          (* the D0, D1 and delta coefficients are defined at w = 1 *)
coefficients = Map[Collect[Expand[#], {_Log, _PolyLog}, Factor] &,
  finite/{prefactor, prefactorEndpoint, prefactorEndpoint, prefactorEndpoint}];
TableForm[Transpose[{ThresholdBasis, coefficients}]]


(* ::Text:: *)
(*The same coefficients split into their contributions (in units of the prefactor; regular parts with the explicit F1, F2 functions inserted):*)


(* ::Input:: *)
Grid[Prepend[MapThread[Prepend[Map[Simplify[#, region] &, #2 /. reduced["MasterFunctions"]], #1] &,
  {finiteNames, Map[#/{prefactor, prefactorEndpoint, prefactorEndpoint, prefactorEndpoint} &, finiteTable]}],
  {"finite part of", "regular", "D0", "D1", "delta"}], Alignment -> Left, Dividers -> {{2 -> True}, {2 -> True}}]


(* ::Text:: *)
(*So sigma_NLO = born0 delta(1-w) + CF alphaEM^2 eq^2 alphaS/(Pi s v) [R/w + B D0 + C D1 + A delta(1-w)] with the four coefficients above: the regular coefficient R is quoted with the factor 1/w of the papers, the distributions are defined at w = 1.*)


(* ::Section:: *)
(*Comparison with the literature*)


(* ::Text:: *)
(*reference/HSV2015-unpolarized.wl transcribes the paper's result exactly as printed: its prefactor and its coefficients of 1, D0, D1 and delta(1-w).  The paper's prefactor contains 1/w and multiplies the distributions as well (and its coefficient of D0 depends on w), so the two forms must be compared as distributions: CanonicalDistribution turns the paper's form into the basis used here with the identities f(w) D0 = f(1) D0 + (f(w) - f(1))/(1-w) and f(w) D1 = f(1) D1 + (f(w) - f(1)) log(1-w)/(1-w).  The difference must vanish identically in the physical region; it is also evaluated numerically at three points.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "HSV2015-unpolarized.wl"}]];
paper = HSV2015Unpolarized[s, v, w, muR, alphaEM, alphaS, eq, CF];
literature = CanonicalDistribution[paper["Prefactor"], paper["Coefficients"], w];
difference = Table[ExactDifference[finite[[i]], literature[[i]], region], {i, 4}];
numeric = NumericDifference[finite, literature, PhysicalPoints, 20];
{difference, numeric}


(* ::Section:: *)
(*Save*)


(* ::Text:: *)
(*The result is saved before the final check, so that a failed comparison still leaves every ingredient on disk for inspection.*)


(* ::Input:: *)
result = <|"Basis" -> ThresholdBasis, "Prefactor" -> prefactor, "Coefficients" -> coefficients, "Finite" -> finite,
  "Born" -> born0, "Virtual" -> virtualFinite, "Real" -> realExpanded["Finite"], "Counterterms" -> counterterms,
  "LiteratureDifference" -> difference, "Route" -> route|>;
SaveResult["Unpolarized", result, route];
Require[difference === {0, 0, 0, 0}, "Finite coefficient differs from HSV 2015", difference];
coefficients[[4]]
