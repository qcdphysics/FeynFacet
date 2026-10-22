(* ::Title:: *)
(*TMD beam and soft functions, step 5: the NLO beam-function matching coefficients and the TMD PDF*)


(* ::Text:: *)
(*Input: Results/<route>/Algebra.wl (step 2), Results/Integrals.wl (step 3) and Results/<route>/Soft.wl (step 4).  Output: Results/<route>/Beam.wl with the matching coefficients I_ij (i, j = q, g) of the beam functions and the TMD PDF coefficients f_ij = I_ij Sqrt[S_i], both compared with Luo, Yang, Zhu and Zhu, arXiv:2012.03256.*)
(**)
(*Assembly, channel by channel, all in units of as = alpha_s/(4 Pi):*)
(*  1. bare real emission: as K(x,eps)/(1-x) T[R] (steps 2 and 3).  Away from x = 1 the regulator drops out and T[R] = F_perp.  For the diagonal channels K/(1-x) = Reg(x,eps) + Gamma_i/(1-x) with Gamma_i = K(1,eps); at the endpoint the regulator gives (step 3)*)
(*       Gamma_i/(1-x) T[R] -> Gamma_i { F_perp [1/(1-x)]_+ + delta(1-x) [ T[2K0] + F_perp LQ/2 ] };*)
(*  2. zero-bin subtraction: the overlap with the soft region is removed by dividing by the bare soft function, at this order by subtracting delta(1-x) S_i^(1),bare = delta(1-x) 4 C_i as T[2K0]: since Gamma_i = 4 C_i, this removes the T[2K0] term of the endpoint exactly;*)
(*  3. collinear (PDF) counterterm: + (as/eps) P^(0)_ij(x), the LO splitting functions in the alpha_s/(4 Pi) normalization, i.e. twice SplittingKernel of FeynFacet/Factorization.wl;*)
(*  4. the remaining poles must be proportional to delta(1-x) (local); their negative is the MSbar renormalization of the beam-function operator.  What is left is I_ij^(1).*)
(*Finally f_ij = I_ij + (1/2) delta_ij delta(1-x) S_i^(1); the dependence on nu cancels when LQ is written as Lh - Lv with Lh = Log[(x n.p)^2/muR^2].*)
(**)
(*Distributions are written with the inert symbols plusD[0, 1-x] = [1/(1-x)]_+ and delta[1-x], as in the reference.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
algebra = LoadResult["Algebra", route];
integrals = LoadResult["Integrals"];
softResult = LoadResult["Soft", route];
channels = {"qq", "qg", "gq", "gg"};
colour = <|"qq" -> CF, "gg" -> CA|>;        (* diagonal channels and their Casimir *)
Fperp = integrals["Fperp"];
softMaster = integrals["SoftMaster"];
Keys[algebra]


(* ::Section:: *)
(*1. Bare real emission*)


(* ::Input:: *)
kernel = Map[Factor[# kT2/gs^2/(1 - x) /. D -> 4 - 2 eps] &, algebra["Real"]];          (* K(x,eps)/(1-x) *)
cusp = AssociationMap[Factor[Limit[(1 - x) kernel[#], x -> 1]] &, Keys[colour]];
regular = Association @ Table[channel -> Factor[kernel[channel] - If[KeyExistsQ[cusp, channel], cusp[channel]/(1 - x), 0]], {channel, channels}];
Require[cusp === <|"qq" -> 4 CF, "gg" -> 4 CA|>, "The endpoint coefficients are not 4 C_i", cusp];
Require[And @@ Table[FreeQ[Limit[regular[channel], x -> 1], ComplexInfinity | Indeterminate | DirectedInfinity], {channel, channels}], "A regular part is singular at x = 1"];
TableForm[Normal[regular]]


(* ::Input:: *)
bare = Association @ Table[channel -> as (Fperp regular[channel] +
    If[KeyExistsQ[cusp, channel], cusp[channel] (Fperp plusD[0, 1 - x] + delta[1 - x] (softMaster + Fperp LQ/2)), 0]), {channel, channels}];
Collect[bare["qq"], {plusD[0, 1 - x], delta[1 - x]}, Simplify]


(* ::Section:: *)
(*2. Zero-bin subtraction*)


(* ::Input:: *)
zeroBin = Association @ Table[channel -> If[KeyExistsQ[colour, channel], -delta[1 - x] softResult["SoftBare"][StringTake[channel, 1]], 0], {channel, channels}];
subtracted = Association @ Table[channel -> Collect[bare[channel] + zeroBin[channel], {plusD[0, 1 - x], delta[1 - x]}, Simplify], {channel, channels}];
Require[FreeQ[subtracted, Lv], "The rapidity scale nu survives the zero-bin subtraction", subtracted];
subtracted["qq"]


(* ::Section:: *)
(*3. Collinear counterterm*)


(* ::Text:: *)
(*P^(0)_ij in the alpha_s/(4 Pi) normalization: 2 x SplittingKernel (which is in alpha_s/(2 Pi) units and returns {regular, D0, delta}), with TF = 1/2 as in step 2.  The endpoint terms, 3 CF and (11 CA - 2 Nf)/3, follow from quark-number and momentum conservation.*)


(* ::Input:: *)
splitting[channel_] := With[{k = 2 SplittingKernel[channel, x] /. TR -> 1/2}, k[[1]] + k[[2]] plusD[0, 1 - x] + k[[3]] delta[1 - x]];
pdfCounterterm = Association @ Table[channel -> as splitting[channel]/eps, {channel, channels}];
TableForm[Normal[Map[Collect[# eps/as, {plusD[0, 1 - x], delta[1 - x]}] &, pdfCounterterm]]]


(* ::Section:: *)
(*4. Operator renormalization and the finite coefficients*)


(* ::Input:: *)
withCounterterm = Association @ Table[channel -> Collect[Normal[Series[subtracted[channel] + pdfCounterterm[channel], {eps, 0, 0}]], {eps, plusD[0, 1 - x], delta[1 - x]}, Simplify], {channel, channels}];
poles = Map[{Simplify[Coefficient[#, eps, -2]], Simplify[Coefficient[#, eps, -1]]} &, withCounterterm];
Require[And @@ Table[poles[channel][[1]] === 0 && FreeQ[Simplify[poles[channel][[2]] /. delta[1 - x] -> 0], x | plusD], {channel, channels}],
  "After the PDF counterterm the poles are not local", poles];
Grid[Prepend[Table[{channel, poles[channel][[1]], poles[channel][[2]]}, {channel, channels}], {"channel", "1/eps^2", "1/eps"}], Frame -> All]


(* ::Text:: *)
(*The remaining single pole is proportional to delta(1-x): it is removed by the renormalization constant of the beam-function operator (its coefficient, -2 C_i LQ + gamma_i with gamma_q = 3 CF and gamma_g = (11 CA - 2 Nf)/3, is the one-loop anomalous dimension of the beam function).  What is left at order eps^0 is the matching coefficient.*)


(* ::Input:: *)
uvCounterterm = Map[-(#[[1]]/eps^2 + #[[2]]/eps) &, poles];
beam = Association @ Table[channel -> Collect[Coefficient[withCounterterm[channel], eps, 0], {plusD[0, 1 - x], delta[1 - x]}, Simplify], {channel, channels}];
TableForm[Normal[beam]]


(* ::Section:: *)
(*Comparison with the literature, and the TMD PDF*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "Luo2020.wl"}]];
reference = Luo2020NLO[x, Lp, LQ, Lv, Lh, CF, CA, 1/2] /. zeta[2] -> Zeta[2];       (* TF = 1/2 *)
compare[ours_, theirs_] := Simplify[Collect[ours - theirs, {plusD[0, 1 - x], delta[1 - x]}, Simplify]];
beamDifference = Association @ Table[channel -> compare[beam[channel]/as, reference["Beam"][channel]], {channel, channels}]


(* ::Input:: *)
tmdpdf = Association @ Table[channel -> Collect[(beam[channel] + If[KeyExistsQ[colour, channel], delta[1 - x] softResult["Soft"][StringTake[channel, 1]]/2, 0]) /. LQ -> Lh - Lv,
  {plusD[0, 1 - x], delta[1 - x]}, Simplify], {channel, channels}];
Require[FreeQ[tmdpdf, Lv], "The TMD PDF coefficients depend on the rapidity scale nu", tmdpdf];
tmdpdfDifference = Association @ Table[channel -> compare[tmdpdf[channel]/as, reference["TMDPDF"][channel]], {channel, channels}]


(* ::Input:: *)
result = <|"Beam" -> beam, "TMDPDF" -> tmdpdf, "Bare" -> bare, "ZeroBin" -> zeroBin, "PDFCounterterm" -> pdfCounterterm,
  "UVCounterterm" -> uvCounterterm, "Kernels" -> kernel, "BeamDifference" -> beamDifference, "TMDPDFDifference" -> tmdpdfDifference, "Route" -> route|>;
SaveResult["Beam", result, route];
Require[And @@ (# === 0 & /@ Values[beamDifference]) && And @@ (# === 0 & /@ Values[tmdpdfDifference]),
  "The beam or TMD PDF coefficients differ from arXiv:2012.03256", {beamDifference, tmdpdfDifference}];
TableForm[Normal[tmdpdf]]
