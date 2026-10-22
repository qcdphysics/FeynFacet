(* ::Package:: *)

(* Factorization.wl: the universal ingredients of MSbar collinear
   subtraction.  Conventions: alpha_s/(2 pi) P(z) is the one-loop splitting
   function; a plus-distribution kernel is stored as {regular, D0, delta}
   meaning regular(z) + D0 [1/(1-z)]_+ + delta delta(1-z). *)

BeginPackage["FeynFacet`"];

MSbarConstant::usage = "MSbarConstant = Log[4 Pi] - EulerGamma, the finite part \
of the MSbar pole 1/eps - EulerGamma + Log[4 Pi].";

SplittingKernel::usage = "SplittingKernel[\"qq\" | \"gq\" | \"qg\" | \"gg\", z] gives the \
unpolarized LO splitting function as a {regular, D0, delta} kernel, in units of \
alpha_s/(2 pi).  SplittingKernel[\"Delta qq\" ...] gives the longitudinally \
polarized kernels Delta P, SplittingKernel[\"Transverse qq\", z] the transversity kernel.  SplittingKernel[\"photon/lepton\", y] and \
SplittingKernel[\"Delta photon/lepton\", y] are the photon-in-lepton kernels in \
units of alpha_em/(2 pi).";

HelicityRestorationKernel::usage = "HelicityRestorationKernel[z] = -4 CF (1-z): the \
finite kernel that restores helicity conservation for an incoming longitudinally \
polarized quark in the HVBM scheme (subtracted once per polarized quark leg).";

Begin["`Private`"];

(* colour factors and Nf are whatever the notebook (or FeynCalc) calls CF, CA, TR, Nf *)
cf := Symbol["CF"]; ca := Symbol["CA"]; tr := Symbol["TR"]; nf := Symbol["Nf"];

MSbarConstant = Log[4 Pi] - EulerGamma;

SplittingKernel["qq", z_] := {-cf (1 + z), 2 cf, 3 cf/2};
SplittingKernel["gq", z_] := {cf (1 + (1 - z)^2)/z, 0, 0};
SplittingKernel["qg", z_] := {tr (z^2 + (1 - z)^2), 0, 0};
(* z/(1-z) = 1/(1-z) - 1: the regular part of 2 CA [z/(1-z)]_+ is -2 CA *)
SplittingKernel["gg", z_] := {2 ca ((1 - z)/z + z (1 - z) - 1), 2 ca,
  (11 ca - 4 tr nf)/6};
SplittingKernel["Delta qq", z_] := SplittingKernel["qq", z];
SplittingKernel["Delta gq", z_] := {cf (2 - z), 0, 0};
SplittingKernel["Delta qg", z_] := {tr (2 z - 1), 0, 0};
SplittingKernel["Delta gg", z_] := {2 ca (1 - 2 z), 2 ca,
  (11 ca - 4 tr nf)/6};
(* transversity: Delta_T P_qq = CF [2 z/(1-z)_+ + 3/2 delta(1-z)], 2z/(1-z) = 2/(1-z) - 2 *)
SplittingKernel["Transverse qq", z_] := {-2 cf, 2 cf, 3 cf/2};
SplittingKernel["photon/lepton", y_] := {(1 + (1 - y)^2)/y, 0, 0};
SplittingKernel["Delta photon/lepton", y_] := {2 - y, 0, 0};

HelicityRestorationKernel[z_] := -4 cf (1 - z);

End[];
EndPackage[];
