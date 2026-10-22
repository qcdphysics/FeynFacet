(* ::Title:: *)
(*TMD beam and soft functions, step 4: the NLO soft function*)


(* ::Text:: *)
(*Input: Results/<route>/Algebra.wl (step 2) and Results/Integrals.wl (step 3).  Output: Results/<route>/Soft.wl with the bare and the renormalized NLO soft function for quarks and gluons, compared with Luo, Yang, Zhu and Zhu, arXiv:2012.03256.*)
(**)
(*The soft function S_i(b, mu, nu) = (1/N_i) <0| Tr[...] |0> is 1 at lowest order.  At NLO the real emission of one gluon from the two Wilson lines gives, with the squared current of step 2, M2_s = 4 gs^2 C_i/((n.k)(nb.k)), and the phase space and regulator of step 3,*)
(**)
(*   S_i^(1),bare = 4 C_i as T[2 K0(2 a kT)],*)
(**)
(*where T[2 K0] is the soft master integral of step 3.  Virtual corrections are scaleless and vanish.  The renormalized soft function is obtained by minimal subtraction of the poles: S_i^(1) = finite part of S_i^(1),bare.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Input:: *)
route = "FeynCalc";
algebra = LoadResult["Algebra", route];
integrals = LoadResult["Integrals"];
algebra["Soft"]


(* ::Section:: *)
(*From the squared current to the bare soft function*)


(* ::Text:: *)
(*Phase space (step 3): d^D k/(2 Pi)^(D-1) delta+(k^2) = d^{D-2}kT d(n.k)/(2 (n.k) (2 Pi)^(D-1)) with (n.k)(nb.k) = kT2.  Hence*)
(**)
(*   S^(1) = Integral d^{D-2}kT/(2 Pi)^(D-2) e^{i b.kT} [M2_s kT2/gs^2] gs^2/(4 Pi kT2) Integral dy/y R(y, kT)   (y = n.k),*)
(**)
(*and with gs^2 = 16 Pi^2 as muR^(2 eps) (e^EulerGamma/(4 Pi))^eps and the definition of T[f] (which contains 4 Pi muR^(2 eps) (e^EulerGamma/(4 Pi))^eps): S^(1) = as [M2_s kT2/gs^2] T[Integral dy/y R].  The bracket is the colour factor 4 C_i, and the regulated longitudinal integral is 2 K0(2 a kT).*)


(* ::Input:: *)
colourFactor = Map[Simplify[# kn kb/gs^2] &, algebra["Soft"]];
Require[colourFactor === <|"q" -> 4 CF, "g" -> 4 CA|>, "Unexpected soft colour factors", colourFactor];
softBare = Map[Collect[# as integrals["SoftMaster"], {as, eps}, Simplify] &, colourFactor]


(* ::Text:: *)
(*Poles: 4 C_i/eps^2 - 4 C_i Lv/eps (the double pole is the cusp, the single pole carries the rapidity logarithm).*)


(* ::Input:: *)
poles = Map[{Coefficient[#, eps, -2], Coefficient[#, eps, -1]} &, softBare];
Require[Simplify[poles["q"] - {4 CF as, -4 CF as Lv}] === {0, 0} && Simplify[poles["g"] - {4 CA as, -4 CA as Lv}] === {0, 0},
  "Unexpected poles of the bare soft function", poles];
poles


(* ::Section:: *)
(*Renormalized soft function and comparison*)


(* ::Input:: *)
soft = Map[Simplify[Coefficient[#, eps, 0]] &, softBare]


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "reference", "Luo2020.wl"}]];
reference = Luo2020NLO[x, Lp, LQ, Lv, Lh, CF, CA, TF] /. zeta[2] -> Zeta[2];
difference = Map[Simplify[soft[#]/as - reference["Soft"][#]] &, <|"q" -> "q", "g" -> "g"|>]


(* ::Input:: *)
result = <|"SoftBare" -> softBare, "Soft" -> soft, "LiteratureDifference" -> difference, "Route" -> route|>;
SaveResult["Soft", result, route];
Require[difference === <|"q" -> 0, "g" -> 0|>, "The soft function differs from arXiv:2012.03256", difference];
soft
