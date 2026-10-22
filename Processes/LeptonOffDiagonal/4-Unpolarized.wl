(* ::Title:: *)
(*Lepton-quark and lepton-gluon scattering with a different tagged parton, step 4: the unpolarized NLO coefficients*)


(* ::Text:: *)
(*Input: Results/<route>/Algebra.wl (step 2) and Results/<route>/Reduced.wl (step 3).  Output: Results/<route>/Unpolarized.wl with the finite coefficients of the two channels l q -> g + X and l g -> q + X and their comparison with Hinderer, Schlegel and Vogelsang, arXiv:1505.06415, Eqs. (27) and (28).*)
(**)
(*Hard coefficient, as in Processes/LeptonQuark: sigma(s,v,w) = s E dsigma/d^3p with t = -s(1-v), u = -s v w, s4 = s v (1-w).  Neither channel has a Born term: the real emission (integrated phase space of the unobserved pair, divided by 32 Pi^3) is the whole order alphaEM^2 alphaS result, its collinear pole is removed by one MSbar counterterm built from the Born process l q -> l q (the fragmentation q -> g of the Born quark for the first channel, the parton distribution g -> q of the incoming quark for the second), and the photon-in-lepton subtraction with the photon Born processes gamma q -> g q and gamma g -> q qbar.  There are no virtual corrections and, since no soft gluon can be emitted (the gluon is observed or absent), no plus distributions: the coefficients are ordinary functions of v and w.*)
(**)
(*Steps: (1) Born l q -> l q and photon Born terms through order eps; (2) real emission of both channels, with the threshold expansion machinery used only to verify that no distribution arises; (3) the counterterms and the pole cancellation; (4) the finite coefficients and the comparison.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
spin = "Unpolarized";
algebra = LoadResult["Algebra", route];
reduced = LoadResult["Reduced", route];
region = s > 0 && muR > 0 && 0 < v < 1 && 0 < w < 1;          (* physical region; fixes every logarithm's branch *)
toVW = {t -> -s (1 - v), s4 -> s v (1 - w), u -> -s v w};
couplings = {ee -> Sqrt[4 Pi alphaEM], gs -> Sqrt[4 Pi alphaS]};
channels = {"QG", "GQ"};
Keys[reduced]


(* ::Section:: *)
(*Born terms*)


(* ::Text:: *)
(*The Born hard part of l q -> l q, bornHard = born0 + eps born1 + ..., multiplies the collinear counterterms; at Born level u = -s v (w = 1).  Expected: born0 = 2 alphaEM^2 eq^2 (1 + v^2)/(s v (1-v)^2).*)


(* ::Input:: *)
bornHard = Series[algebra["Born"][spin]/(16 Pi^2 s v) /. D -> 4 - 2 eps /. couplings /. {t -> -s (1 - v), u -> -s v}, {eps, 0, 1}];
born0 = SeriesCoefficient[bornHard, {eps, 0, 0}] // Simplify;
born1 = SeriesCoefficient[bornHard, {eps, 0, 1}] // Simplify;
Require[Simplify[born0 - 2 alphaEM^2 eq^2 (1 + v^2)/(s v (1 - v)^2)] === 0, "Born hard part differs from Eq. (22) of HSV"];
{born0, born1}


(* ::Text:: *)
(*Photon Born terms gamma(q) + q(k) -> g(p) + q(r) and gamma(q) + g(k) -> q(p) + qbar(r).  The photon carries the fraction y of the lepton momentum, so that sh = y s, th = -s(1-v), uh = -y s v w, and sh + th + uh = 0 fixes y = (1-v)/(1-vw).  Their hard parts photonHard = photon0 + eps photon1 are normalized like a Born term with the Jacobian of the photon fraction, |M|^2/(16 Pi^2 y s (1 - v w)), in units of alphaS/Pi.  Expected (Eq. 29 of arXiv:1505.06415): 2 Pi CF alphaEM eq^2 (1 + (1-vw)^2)/(s (1-v) (1-vw)) and 2 Pi TR alphaEM eq^2 (v^2 w^2 + (1-vw)^2)/(s (1-v) v w (1-vw)).*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "HSV2015-offdiagonal.wl"}]];
Get[FileNameJoin[{NotebookDirectory[], "reference", "HSV2017-offdiagonal.wl"}]];
yPhoton = (1 - v)/(1 - v w);
photonHard = Association @ Table[channel -> Series[(Pi/alphaS) algebra["Photon"][channel][spin]/(16 Pi^2 yPhoton s (1 - v w)) /. D -> 4 - 2 eps /. couplings /.
  {sh -> yPhoton s, th -> -s (1 - v), uh -> -yPhoton s v w}, {eps, 0, 1}], {channel, channels}];
photon0 = Map[Simplify[SeriesCoefficient[#, {eps, 0, 0}], region] &, photonHard];
photon1 = Map[Simplify[SeriesCoefficient[#, {eps, 0, 1}], region] &, photonHard];
expectedPhoton = <|"QG" -> HSV2015PhotonBorn["gamma q -> g", s, v, w, alphaEM, eq, CF], "GQ" -> HSV2015PhotonBorn["gamma g -> q", s, v, w, alphaEM, eq, TR]|>;
Do[Require[Simplify[photon0[channel] - expectedPhoton[channel], region] === 0, "A photon Born term differs from Eq. (29) of HSV 2015", channel], {channel, channels}];
{photon0, photon1}


(* ::Section:: *)
(*Real emission*)


(* ::Text:: *)
(*The integrated real contributions of step 3 in the variables v, w, with the couplings and the scale muR^(2 eps) of the MSbar coupling.  ThresholdSplit and EndpointExpand are the tools of the lepton-quark process: they would turn any power (1-w)^(-1-lambda eps) into plus distributions.  Here the only non-integer powers of 1-w come from s4^(-eps) and (1-z)^eps of the angular masters with no 1/(1-w) in front, so the branch coefficients vanish at w = 1 and the result has no double pole and no D0, D1 or delta term: this is verified.  The single pole is an ordinary function of w.*)


(* ::Input:: *)
realHard = Association @ Table[channel -> (muR^(2 eps) reduced["Real"][channel][spin]/(32 Pi^3) /. couplings /. toVW), {channel, channels}];
split = Map[ThresholdSplit[#, w, eps] &, realHard];
realExpanded = Map[EndpointExpand[#["Regular"], #["Branches"], eps, w, region] &, split];
Do[Require[Simplify[realExpanded[channel]["Pole2"], region] === 0 && Simplify[realExpanded[channel]["Pole1"][[2 ;; 3]], region] === {0, 0} &&
    Simplify[realExpanded[channel]["Finite"][[2 ;; 4]], region] === {0, 0, 0},
  "A distribution appeared in an off-diagonal channel", channel], {channel, channels}];
realPole = Map[Simplify[#["Pole1"][[1]] /. reduced["MasterFunctions"], region] &, realExpanded]


(* ::Section:: *)
(*Collinear counterterms and the pole cancellation*)


(* ::Text:: *)
(*With the hadronic convolutions of docs/CONVENTIONS.md (the same as in the lepton-quark process), a kernel P(z) = R(z) + c0 [1/(1-z)]_+ + cdelta delta(1-z) in units of alphaS/(2 Pi) gives*)
(**)
(*  incoming parton:   + (alphaS/(2 Pi)) (1/eps) integral dz/z P(z) sigma_Born(z s, z t, u)  =  (alphaS/(2 Pi)) (1/eps) b0 {(R(w) + c0)/w, c0, cdelta}   in the basis {regular, D0, delta},*)
(*  fragmentation:     + (alphaS/(2 Pi)) (1/eps) integral dz/z^2 P(z) sigma_Born(s, t/z, u/z)  =  (alphaS/(2 Pi)) (1/eps) (1/v) {c0 (j - j(1))/(1-w) + v R(dw) j, c0 j(1), j(1) (cdelta + c0 log v)},  dw = 1 - v + v w,  j = v w b0(v -> v w/dw)/dw^2,*)
(**)
(*where the Born delta function fixes z = w for the incoming leg and z = dw for the fragmenting one (the kernels here have c0 = cdelta = 0).  The channel l q -> g + X needs the fragmentation q -> g of the Born quark with P_gq(z) = CF (1 + (1-z)^2)/z; the channel l g -> q + X needs the parton distribution g -> q with P_qg(z) = TR (z^2 + (1-z)^2).  The photon in the lepton gives + (alphaEM/(2 Pi)) (1/eps) P_photon/lepton(y) times the photon hard part at y = (1-v)/(1-vw).  Every counterterm multiplies a D-dimensional Born term, whose order-eps part, with the MSbar constant log(4 Pi) - EulerGamma (and the order-eps Jacobian 2 j0 log(dw) of the fragmentation), feeds the finite part.  The table at the end lists each contribution to the single pole, whose sum must vanish.*)


(* ::Input:: *)
kernelPrefix = If[spin === "Polarized", "Delta ", ""];                    (* the incoming partons are polarized in the LL case, the fragmenting one is not *)
incomingTuple[kernel_, bornTerm_] := bornTerm {(kernel[[1]] + kernel[[2]])/w, kernel[[2]], kernel[[3]]};
dw = 1 - v + v w;
jFunction[bornTerm_] := v w (bornTerm /. v -> v w/dw)/dw^2;
fragmentationTuple[kernel_, j_] := Module[{jEnd = Simplify[j /. w -> 1, region]},
  (1/v) {kernel[[2]] (j - jEnd)/(1 - w) + v kernel[[1]] j, kernel[[2]] jEnd, jEnd (kernel[[3]] + kernel[[2]] Log[v])}];
photonKernel = SplittingKernel[kernelPrefix <> "photon/lepton", yPhoton][[1]];
partonPole = <|
  "QG" -> alphaS/(2 Pi) fragmentationTuple[SplittingKernel["gq", dw], jFunction[born0]],
  "GQ" -> alphaS/(2 Pi) incomingTuple[SplittingKernel[kernelPrefix <> "qg", w], born0]|>;
partonFinite = <|
  "QG" -> alphaS/(2 Pi) (MSbarConstant fragmentationTuple[SplittingKernel["gq", dw], jFunction[born0]] +
    fragmentationTuple[SplittingKernel["gq", dw], jFunction[born1] + 2 jFunction[born0] Log[dw]]),
  "GQ" -> alphaS/(2 Pi) (MSbarConstant incomingTuple[SplittingKernel[kernelPrefix <> "qg", w], born0] + incomingTuple[SplittingKernel[kernelPrefix <> "qg", w], born1])|>;
photonPole = Map[{alphaEM/(2 Pi) photonKernel alphaS/Pi #, 0, 0} &, photon0];
photonFinite = Association @ Table[channel -> {alphaEM/(2 Pi) photonKernel alphaS/Pi (MSbarConstant photon0[channel] + photon1[channel]), 0, 0}, {channel, channels}];
counterName = <|"QG" -> "fragmentation q -> g", "GQ" -> "parton distribution g -> q"|>;
poleTable = Association @ Table[channel -> {realExpanded[channel]["Pole1"], partonPole[channel], photonPole[channel]}, {channel, channels}];
singlePoles = Map[Total, poleTable];
Do[Require[And @@ (ExactDifference[#, 0, region] === 0 & /@ singlePoles[channel]), "Single poles do not cancel", {channel, singlePoles[channel]}], {channel, channels}];
Column[Table[Grid[Prepend[MapThread[Prepend[Map[Simplify[#, region] &, #2], #1] &,
  {{"real emission", counterName[channel], "photon in lepton", "sum"}, Append[poleTable[channel], singlePoles[channel]]}],
  {channel <> ": 1/eps coefficient of", "regular", "D0", "delta"}], Alignment -> Left, Dividers -> {{2 -> True}, {2 -> True, 5 -> True}}], {channel, channels}]]


(* ::Section:: *)
(*The finite coefficients*)


(* ::Text:: *)
(*The regular parts still contain the inert symbols F1[n][z], F2[n][z] of the hypergeometric functions, replaced here by their explicit forms.  The finite coefficient is the sum of the real emission and the two counterterms; it is an ordinary function (its D0, D1 and delta components vanish), quoted in units of the prefactor alphaEM^2 eq^2 CF alphaS/(Pi s v w) (TR for the second channel) of the paper.*)


(* ::Input:: *)
colourOf = <|"QG" -> CF, "GQ" -> TR|>;
finite = Association @ Table[channel -> ((realExpanded[channel]["Finite"] + {#[[1]], #[[2]], 0, #[[3]]} &[partonFinite[channel] + photonFinite[channel]]) /. reduced["MasterFunctions"]), {channel, channels}];
Do[Require[Simplify[finite[channel][[2 ;; 4]], region] === {0, 0, 0}, "The finite coefficient has a distribution part", channel], {channel, channels}];
prefactor = Association @ Table[channel -> colourOf[channel] alphaEM^2 eq^2 alphaS/(Pi s v w), {channel, channels}];
coefficients = Association @ Table[channel -> Collect[Expand[finite[channel][[1]]/prefactor[channel]], {_Log, _PolyLog}, Factor], {channel, channels}];
coefficients["GQ"]


(* ::Section:: *)
(*Comparison with the literature*)


(* ::Text:: *)
(*reference/HSV2015-offdiagonal.wl transcribes Eqs. (27), (28) and Appendix A of arXiv:1505.06415 exactly as printed (the coefficient of alphaS/Pi, with a single scale mu for the renormalization and both factorizations, as here).  The difference must vanish identically in the physical region; it is also evaluated numerically at three points.*)


(* ::Input:: *)
literature = <|"QG" -> alphaS/Pi HSV2015OffDiagonal["q -> g", s, v, w, muR, alphaEM, eq, CF], "GQ" -> alphaS/Pi HSV2015OffDiagonal["g -> q", s, v, w, muR, alphaEM, eq, TR]|>;
difference = Association @ Table[channel -> ExactDifference[finite[channel][[1]], literature[channel], region], {channel, channels}];
numeric = Association @ Table[channel -> NumericDifference[{finite[channel][[1]]}, {literature[channel]}, PhysicalPoints, 20], {channel, channels}];
{difference, numeric}


(* ::Section:: *)
(*Save*)


(* ::Text:: *)
(*The result is saved before the final check, so that a failed comparison still leaves every ingredient on disk for inspection.*)


(* ::Input:: *)
result = <|"Basis" -> ThresholdBasis, "Prefactor" -> prefactor, "Coefficients" -> coefficients, "Finite" -> finite, "Born" -> born0,
  "Real" -> Map[#["Finite"] &, realExpanded], "Counterterms" -> <|"parton" -> partonFinite, "photon" -> photonFinite|>,
  "LiteratureDifference" -> difference, "Route" -> route, "Spin" -> spin|>;
SaveResult[spin, result, route];
Require[Values[difference] === {0, 0}, "Finite coefficients differ from HSV 2015", difference];
coefficients
