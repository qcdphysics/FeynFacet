(* Separate collinear factorization scales for single-inclusive lepton
   scattering. All hard parts use s E_c d sigma/d^3 p_c. Vectors have basis
   {regular, [1/(1-w)]_+, [log(1-w)/(1-w)]_+, delta(1-w)}.
   The residues include their couplings; they are not dimensionless tables. *)
BeginPackage["FeynFacet`LeptonScales`"];

LeptonScaleResidues::usage =
 "LeptonScaleResidues[bornD,photonD,spin,channel,{v,w,eps},{alphaS,alphaEM,CF,TR}] " <>
 "returns the PDF, FF and photon-in-lepton collinear residues. bornD is the " <>
 "dimensionally continued lepton-quark Born hard part without its delta " <>
 "function. photonD is the physical photon Born hard part including alphaS " <>
 "and its photon-fraction Jacobian. spin is Unpolarized or Polarized; " <>
 "channel is QQ, QG or GQ (incoming and observed partons).";

LeptonHardAtScales::usage =
 "LeptonHardAtScales[common,residues,{muR,muPDF,muFF,muPhoton}] restores " <>
 "independent collinear scales to the common-scale finite hard-part vector. " <>
 "common must include all overall prefactors, in the same basis as residues.";

Begin["`Private`"];

threshold[{regular_, plus_, endpoint_}] := {regular, plus, 0, endpoint};

LeptonScaleResidues[bornD_, photonD_, spin_String, channel_String,
    {vv_Symbol, ww_Symbol, epsilon_Symbol}, {as_, ae_, cf_, tr_}] := Module[
  {born, photon, z, y, j, endpoint, pdf = {0,0,0,0}, ff = {0,0,0,0},
   kernel, singular, delta, photonKernel},

  If[!MemberQ[{"Unpolarized","Polarized"}, spin] ||
     !MemberQ[{"QQ","QG","GQ"}, channel],
    Return[Failure["LeptonScaleChannel", <|"Spin"->spin,"Channel"->channel|>]]];
  born = SeriesCoefficient[bornD, {epsilon,0,0}];
  photon = SeriesCoefficient[photonD, {epsilon,0,0}];
  z = 1-vv+vv ww;
  y = (1-vv)/(1-vv ww);

  (* Incoming convolution: the Born delta fixes the splitting fraction to w.
     Multiplication of a plus distribution by 1/w adds its regular 1/w term. *)
  If[channel === "QQ",
    pdf = as/(2 Pi) cf born {(1-ww)/ww, 2, 0, 3/2}];
  If[channel === "GQ",
    kernel = If[spin === "Polarized", tr (2 ww-1), tr (ww^2+(1-ww)^2)];
    pdf = {as/(2 Pi) born kernel/ww,0,0,0}];

  (* Fragmentation fixes z=1-v+v w and vBorn=v w/z. Pulling back the original
     z-plus distribution supplies delta(1-w) log(v)/v as well as D0(w)/v. *)
  If[MemberQ[{"QQ","QG"},channel],
    j = vv ww (born /. vv -> vv ww/z)/z^2;
    endpoint = FullSimplify[j /. ww -> 1, 0<vv<1];
    {kernel,singular,delta} = If[channel === "QQ",
      {-cf (1+z),2 cf,3 cf/2}, {cf (1+(1-z)^2)/z,0,0}];
    ff = as/(2 Pi) threshold[{
      kernel j + singular (j-endpoint)/(vv (1-ww)),
      singular endpoint/vv,
      endpoint (delta+singular Log[vv])/vv}]];

  photonKernel = If[spin === "Polarized",2-y,(1+(1-y)^2)/y];
  <|"PDF"->pdf,"FF"->ff,
    "Photon"->{ae/(2 Pi) photonKernel photon,0,0,0}|>
];

LeptonHardAtScales[common_, residues_Association,
    {muR_, muPDF_, muFF_, muPhoton_}] := (
  common - residues["PDF"] Log[muPDF^2/muR^2]
         - residues["FF"] Log[muFF^2/muR^2]
         - residues["Photon"] Log[muPhoton^2/muR^2]);

End[];
EndPackage[];
