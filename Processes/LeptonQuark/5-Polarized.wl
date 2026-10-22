(* ::Title:: *)
(*Lepton-quark scattering, step 5: the double-longitudinal NLO hard coefficient*)


(* ::Text:: *)
(*Same structure as 4-Unpolarized.wl (read that first), for the helicity difference Delta sigma = (sigma(++) - sigma(+-))/2 of the incoming lepton and quark.  Input: Results/<route>/Algebra.wl and Reduced.wl.  Output: Results/<route>/Polarized.wl and the comparison with Hinderer, Schlegel and Vogelsang, arXiv:1703.10872.*)
(**)
(*Three things differ from the unpolarized case.  (1) The HVBM scheme: the polarized matrix elements of step 2 contain the hatted invariant spe3, which step 3 turned into integrals in D + 2 dimensions; they are already included in reduced["Real"]["Polarized"].  (2) The splitting functions: the incoming quark and the photon in the lepton are polarized, Delta Pqq = Pqq and Delta P_photon/lepton(y) = 2 - y; the outgoing quark fragments unpolarized.  (3) The HVBM scheme breaks helicity conservation at the quark-gluon vertex by a finite amount; the standard remedy (Ref. Vogelsang 1996, used in arXiv:1703.10872) is a finite renormalization of the polarized quark density with the kernel Delta f_qq = -4 CF (1-z), i.e. the hard part receives the additional term -(alphaS/(2 Pi)) Delta f_qq (x) Delta sigma_Born.*)
(**)
(*All factorization scales (incoming quark, outgoing quark, photon in lepton) are set equal to muR, as in that paper.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
algebra = LoadResult["Algebra", route];
reduced = LoadResult["Reduced", route];
region = s > 0 && muR > 0 && 0 < v < 1 && 0 < w < 1;
toVW = {t -> -s (1 - v), s4 -> s v (1 - w), u -> -s v w};
couplings = {ee -> Sqrt[4 Pi alphaEM], gs -> Sqrt[4 Pi alphaS]};
Keys[reduced]


(* ::Section:: *)
(*Born terms*)


(* ::Text:: *)
(*The polarized Born square has no order-eps term; born0 = 2 alphaEM^2 eq^2 (1+v)/(s v (1-v)).*)


(* ::Input:: *)
bornHard = Series[algebra["Born"]["Polarized"]/(16 Pi^2 s v) /. D -> 4 - 2 eps /. couplings /. {t -> -s (1 - v), u -> -s v}, {eps, 0, 1}];
born0 = SeriesCoefficient[bornHard, {eps, 0, 0}] // Simplify;
born1 = SeriesCoefficient[bornHard, {eps, 0, 1}] // Simplify;
Require[Simplify[born0 - 2 alphaEM^2 eq^2 (1 + v)/(s v (1 - v))] === 0, "Polarized Born hard part differs from the expected form"];
{born0, born1}


(* ::Text:: *)
(*Polarized photon-quark Born, normalized as in the unpolarized case; its order-eps part (the evanescent term of the HVBM scheme) enters the finite photon-in-lepton subtraction.*)


(* ::Input:: *)
yPhoton = (1 - v)/(1 - v w);
photonHard = Series[(Pi/alphaS) algebra["Photon"]["Polarized"]/(16 Pi^2 yPhoton s (1 - v w)) /. D -> 4 - 2 eps /. couplings /.
  {sh -> yPhoton s, th -> -s (1 - v), uh -> -yPhoton s v w}, {eps, 0, 1}];
photon0 = Simplify[SeriesCoefficient[photonHard, {eps, 0, 0}], region];
photon1 = Simplify[SeriesCoefficient[photonHard, {eps, 0, 1}], region];
{photon0, photon1}


(* ::Section:: *)
(*Real emission and its threshold expansion*)


(* ::Text:: *)
(*As in the unpolarized case: the exact threshold factors s4^(-eps), zAB^(-1-eps) and zAB^(-eps) give the branches (1-w)^(-1-2eps) and (1-w)^(-1-eps); ThresholdSplit separates them and EndpointExpand turns each into plus distributions.  The polarized result includes the dimension-shifted (D+2) masters of the hatted terms.*)


(* ::Input:: *)
realHard = muR^(2 eps) reduced["Real"]["Polarized"]/(32 Pi^3) /. couplings /. toVW;
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
(*The vector form factor multiplies the polarized Born term; helicity is conserved at the photon vertex.*)


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
(*As in the unpolarized case, with Delta Pqq(w)/w = CF [(1-w)/w + 2/(1-w)]_+ + 3/2 delta(1-w) for the incoming polarized quark, the unpolarized Pqq for the fragmenting quark, and Delta P_photon/lepton(y) = 2 - y.*)


(* ::Input:: *)
incomingTuple[bornTerm_] := CF bornTerm {(1 - w)/w, 2, 3/2};
dw = 1 - v + v w;
jFunction[bornTerm_] := v w (bornTerm /. v -> v w/dw)/dw^2;              (* Born at the shifted point, with the Jacobian *)
fragmentation[j_] := Module[{jEnd = Simplify[j /. w -> 1, region]},
  CF/v {2 (j - jEnd)/(1 - w) - v j (1 + dw), 2 jEnd, jEnd (2 Log[v] + 3/2)}];
photonKernel = 2 - yPhoton;                                                  (* Delta P_photon/lepton *)
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
(*The regular parts still contain the inert symbols F1[n][zAB] etc. of the hypergeometric functions; they are replaced by their explicit expressions (logarithms and dilogarithms) here.  Finite parts of the counterterms as before, plus the helicity-restoring term -(alphaS/(2 Pi)) Delta f_qq (x) Delta sigma_Born with Delta f_qq(z) = -4 CF (1-z) on the incoming quark leg (HelicityRestorationKernel), which contributes +(alphaS/(2 Pi)) 4 CF (1-w)/w born0 to the regular part.*)


(* ::Input:: *)
pdfFinite = alphaS/(2 Pi) (MSbarConstant incomingTuple[born0] + incomingTuple[born1]);
ffFinite = alphaS/(2 Pi) (MSbarConstant fragmentation[jFunction[born0]] + fragmentation[jFunction[born1] + 2 jFunction[born0] Log[dw]]);
photonFinite = {alphaEM/(2 Pi) photonKernel alphaS/Pi (MSbarConstant photon0 + photon1), 0, 0};
helicityRestoration = {-alphaS/(2 Pi) HelicityRestorationKernel[w]/w born0, 0, 0};
counterterms = {#[[1]], #[[2]], 0, #[[3]]} &[pdfFinite + ffFinite + photonFinite + helicityRestoration];
finite = realExpanded["Finite"] + {0, 0, 0, virtualFinite} + counterterms;
finiteTable = {realExpanded["Finite"], {0, 0, 0, virtualFinite}, {#[[1]], #[[2]], 0, #[[3]]} &[pdfFinite],
  {#[[1]], #[[2]], 0, #[[3]]} &[ffFinite], {#[[1]], #[[2]], 0, #[[3]]} &[photonFinite], {#[[1]], #[[2]], 0, #[[3]]} &[helicityRestoration]};
finiteNames = {"real emission", "virtual", "incoming quark", "outgoing quark", "photon in lepton", "helicity restoration"};
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


(* ::Section:: *)
(*Comparison with the literature*)


(* ::Text:: *)
(*reference/HSV2017-polarized.wl transcribes the paper's result exactly as printed: its prefactor and its coefficients of 1, D0, D1 and delta(1-w).  The paper's prefactor contains 1/w and multiplies the distributions as well, so the two forms must be compared as distributions: CanonicalDistribution turns the paper's form into the basis used here with the identities f(w) D0 = f(1) D0 + (f(w) - f(1))/(1-w) and f(w) D1 = f(1) D1 + (f(w) - f(1)) log(1-w)/(1-w).  The difference must vanish identically in the physical region; it is also evaluated numerically at three points.*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "HSV2017-polarized.wl"}]];
paper = HSV2017Polarized[s, v, w, muR, alphaEM, alphaS, eq, CF];
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
SaveResult["Polarized", result, route];
Require[difference === {0, 0, 0, 0}, "Finite coefficient differs from HSV 2017", difference];
coefficients[[4]]
