(* ::Package:: *)

(* Amplitudes.wl: load FeynCalc with FeynArts, and turn an FCFAConvert
   amplitude into a plain record that any algebra program can read.

   A record is one term of an amplitude written the way one writes it by hand:

     <| "Coefficient" -> I ee^2 eq gs ColorT[Glu5, Col4, Col2] Den[k - r] Den[k - p - r]
                         Metric[Lor2, Lor3] PolVector[r, Lor1, -1],
        "Chains" -> { Chain[Ubar[lp], GammaMatrix[Lor2], Uspinor[l]],
                      Chain[Ubar[p], GammaMatrix[Lor3], Slash[k - r], GammaMatrix[Lor1], Uspinor[k]] } |>

   Vocabulary (all symbols live in Global` so that they print without context):
     Chain[Ubar[p], ..., Uspinor[k]]   open fermion line, read left to
                                 right as ubar(p) ... u(k); Vbar/Vspinor
                                 for antifermions
     FermionLoop[GammaMatrix[mu], Slash[loop], ...]   closed fermion loop (the Dirac
                                 trace of one whole gamma word), listed under
                                 "Traces"; its commuting factors (couplings,
                                 generators, signs) are in "Coefficient".
                                 Legacy records list those
                                 factors as extra one-element entries
                                 FermionLoop[gs^2], ... under "Traces"; the
                                 algebra routes multiply such entries into the
                                 coefficient, so both forms give the same result
     GammaMatrix[mu]             gamma^mu            Slash[p]   p-slash
     GammaFive                   gamma_5
     Den[q]                      1/q^2 for a massless propagator with momentum q
     Metric[mu, nu]              g^{mu nu}           MomentumComponent[p, mu]  p^mu
     MomentumProduct[p, q]       p.q (appears only after contractions)
     PolVector[p, mu, sign]  polarization vector of the external boson with
                                 momentum p (sign -1: conjugated, outgoing)
     ColorT[{a, b, ..}, i, j]    (T^a T^b ...)_{ij}  ColorF[a, b, c]  f^{abc}
     ColorDelta[i, j]            delta_{ij} (fundamental)
     GluonDelta[a, b]            delta^{ab} (adjoint)
     Epsilon4[a, b, c, d]      epsilon tensor (four dimensional)
     Eikonal[q, n]               1/(n.q), the denominator of a Wilson-line
                                 (eikonal) propagator with momentum q along
                                 the line in the direction n (TMD process)
     ee, gs, eq                  electromagnetic and strong couplings, quark
                                 charge in units of ee
   All Dirac matrices and momenta are D dimensional unless the notebook says
   otherwise. *)

BeginPackage["FeynFacet`"];

LoadFeynCalc::usage = "LoadFeynCalc[] loads FeynCalc together with FeynArts and \
silences their start-up output.  Evaluate it in a cell of its own.";

AmplitudeRecords::usage = "AmplitudeRecords[amplitude] converts one FCFAConvert \
amplitude (obtained with Contract -> False) into a list of records \
<|\"Coefficient\" -> ..., \"Chains\" -> {Chain[...], ...}|>, one per term, after \
setting all gauge parameters to one (Feynman gauge).  See Amplitudes.wl for the \
vocabulary.";

ShowRecords::usage = "ShowRecords[records] displays a list of amplitude records as a table.";

SaveAmplitudes::usage = "SaveAmplitudes[association] writes Results/Amplitudes.wl with \
full symbol contexts, so that a kernel without FeynCalc can read it.";

LoadAmplitudes::usage = "LoadAmplitudes[] reads Results/Amplitudes.wl; LoadAmplitudes[\"class-1\"] reads \
Results/class-1/Amplitudes.wl (written by SaveAmplitudes[data, \"class-1\"]).";

Begin["`Private`"];

LoadFeynCalc[] := (
  Global`$FeynCalcStartupMessages = False;
  Global`$LoadAddOns = {"FeynArts"};
  Needs["FeynCalc`"];
  FeynArts`$FAVerbose = 0;
  "FeynCalc " <> FeynCalc`$FeynCalcVersion <> ", FeynArts " <> FeynArts`$FeynArtsVersion);

(* ---- scalars ------------------------------------------------------------ *)
scalar[x_?NumericQ] := x;
scalar[x_Symbol] := x;
scalar[x_Plus] := scalar /@ x;
scalar[x_Times] := scalar /@ x;
scalar[Power[x_, n_]] := Power[scalar[x], n];
scalar[FeynCalc`SMP["e"]] := Global`ee;
scalar[FeynCalc`SMP["g_s"]] := Global`gs;
scalar[FeynCalc`SMP[name_]] := Symbol["Global`" <> StringReplace[name, "_" -> ""]];
scalar[FeynCalc`Pair[FeynCalc`LorentzIndex[a_, ___], FeynCalc`LorentzIndex[b_, ___]]] := Global`Metric[a, b];
scalar[FeynCalc`Pair[FeynCalc`LorentzIndex[a_, ___], FeynCalc`Momentum[FeynCalc`Polarization[p_, sign_, ___], ___]]] :=
  Global`PolVector[p, a, sign/I];
scalar[FeynCalc`Pair[FeynCalc`Momentum[FeynCalc`Polarization[p_, sign_, ___], ___], FeynCalc`LorentzIndex[a_, ___]]] :=
  Global`PolVector[p, a, sign/I];
scalar[FeynCalc`Pair[FeynCalc`LorentzIndex[a_, ___], FeynCalc`Momentum[p_, ___]]] := Global`MomentumComponent[p, a];
scalar[FeynCalc`Pair[FeynCalc`Momentum[p_, ___], FeynCalc`LorentzIndex[a_, ___]]] := Global`MomentumComponent[p, a];
scalar[FeynCalc`Pair[FeynCalc`Momentum[p_, ___], FeynCalc`Momentum[q_, ___]]] := Global`MomentumProduct[p, q];
scalar[FeynCalc`FeynAmpDenominator[ds___]] := Times @@ (denominator /@ {ds});
scalar[FeynCalc`SUNTF[a_List, FeynCalc`SUNFIndex[i_], FeynCalc`SUNFIndex[j_]]] :=
  Global`ColorT[a /. FeynCalc`SUNIndex[x_] :> x, i, j];
scalar[FeynCalc`SUNFDelta[FeynCalc`SUNFIndex[i_], FeynCalc`SUNFIndex[j_]]] := Global`ColorDelta[i, j];
scalar[FeynCalc`SUNDelta[FeynCalc`SUNIndex[a_], FeynCalc`SUNIndex[b_]]] := Global`GluonDelta[a, b];
scalar[FeynCalc`SUNF[FeynCalc`SUNIndex[a_], FeynCalc`SUNIndex[b_], FeynCalc`SUNIndex[c_]]] := Global`ColorF[a, b, c];
scalar[FeynCalc`Eps[args__]] := Global`Epsilon4 @@ (Replace[{args},
  {FeynCalc`Momentum[p_, ___] :> p, FeynCalc`LorentzIndex[a_, ___] :> a}, {1}]);
scalar[Global`FAEikonal[q_, v_]] := Global`Eikonal[Expand[q], v];
scalar[Global`Eikonal[q_, v_]] := Global`Eikonal[Expand[q], v];
scalar[x_] := (Message[AmplitudeRecords::unknown, x]; Throw[$Failed, "records"]);
AmplitudeRecords::unknown = "Unrecognised object in the amplitude: `1`";

denominator[FeynCalc`PropagatorDenominator[c_. FeynCalc`Momentum[q_, ___], 0]] := Global`Den[canonical[c q]];
denominator[x_] := scalar[x];
(* Den[q] = Den[-q]: write the momentum with a positive leading coefficient *)
canonical[q_] := Module[{e = Expand[q], lead},
  lead = If[Head[e] === Plus, First[e], e];
  If[MatchQ[lead, Times[_?Negative, ___] | _?Negative], -e, e]];

(* ---- fermion chains ----------------------------------------------------- *)
gamma[FeynCalc`DiracGamma[FeynCalc`LorentzIndex[a_, ___], ___]] := Global`GammaMatrix[a];
gamma[FeynCalc`DiracGamma[FeynCalc`Momentum[p_, ___], ___]] := Global`Slash[p];
gamma[FeynCalc`DiracGamma[5]] := Global`GammaFive;
gamma[x_] := scalar[x];
spinor[FeynCalc`Spinor[c_. FeynCalc`Momentum[p_, ___], 0, _], side_] :=
  Switch[{side, c}, {"left", 1}, Global`Ubar[p], {"right", 1}, Global`Uspinor[p],
    {"left", -1}, Global`Vbar[p], {"right", -1}, Global`Vspinor[p]];
chain[d_Dot] := With[{f = List @@ d},
  Global`Chain @@ Join[{spinor[First[f], "left"]}, gamma /@ f[[2 ;; -2]], {spinor[Last[f], "right"]}]];

(* ---- one amplitude -> records --------------------------------------------- *)
termList[e_] := If[Head[e] === Plus, List @@ e, {e}];
loopTerm[t_] := Module[{factors = If[Head[t] === Times, List @@ t, {t}], word},
  word = Select[factors, !FreeQ[#, FeynCalc`DiracGamma] &];
  (Times @@ Select[factors, FreeQ[#, FeynCalc`DiracGamma] &]) FeynCalc`DiracTrace[Times @@ word]];
trace[FeynCalc`DiracTrace[d_Dot]] := Global`FermionLoop @@ (gamma /@ (List @@ d));
trace[FeynCalc`DiracTrace[d_]] := Global`FermionLoop[gamma[d]];
record[term_] := Module[{factors, chains, traces, rest},
  factors = If[Head[term] === Times, List @@ term, {term}];
  chains = Select[factors, !FreeQ[#, FeynCalc`Spinor] &];
  traces = Select[factors, MatchQ[#, _FeynCalc`DiracTrace] &];
  rest = Times @@ Select[factors, FreeQ[#, FeynCalc`Spinor] && !MatchQ[#, _FeynCalc`DiracTrace] &];
  If[traces === {}, <|"Coefficient" -> scalar[rest], "Chains" -> chain /@ chains|>,
    <|"Coefficient" -> scalar[rest], "Chains" -> chain /@ chains, "Traces" -> trace /@ traces|>]];

AmplitudeRecords[amplitude_] := Catch[Module[{expr, terms},
  expr = amplitude /. FeynCalc`GaugeXi[_] -> 1;
  (* a closed fermion loop comes as DiracTrace[sum of chains]: expand it into one trace per term.  The trace is
     linear over the terms of a sum only; a single term c (g1.g2...gn) keeps its whole gamma word g1.g2...gn
     inside one trace and its commuting coefficient c outside (mapping over the Dot itself would trace each
     gamma matrix separately) *)
  expr = expr /. FeynCalc`DiracTrace[x_] :> Total[loopTerm /@
      termList[Expand[FeynCalc`DotSimplify[x, FeynCalc`Expanding -> True]]]];
  expr = Expand[FeynCalc`DotSimplify[expr, FeynCalc`Expanding -> True]];
  terms = If[Head[expr] === Plus, List @@ expr, {expr}];
  record /@ terms], "records"];

ShowRecords[records_List] := Grid[Prepend[
  Table[{i, records[[i, "Coefficient"]], Column[records[[i, "Chains"]]]}, {i, Length[records]}],
  {"term", "coefficient", "fermion chains"}], Alignment -> Left, Frame -> All];

SaveAmplitudes[data_Association] := (
  Block[{$ContextPath = {"System`"}},
    Put[data, FileNameJoin[{ResultsDirectory[], "Amplitudes.wl"}]]];
  FileNameJoin[{ResultsDirectory[], "Amplitudes.wl"}]);
(* processes with several classes keep one file per class under Results/<subdirectory>/ *)
SaveAmplitudes[data_Association, subdirectory_String] := (
  Block[{$ContextPath = {"System`"}},
    Put[data, FileNameJoin[{ResultsDirectory[subdirectory], "Amplitudes.wl"}]]];
  FileNameJoin[{ResultsDirectory[subdirectory], "Amplitudes.wl"}]);

LoadAmplitudes[] := LoadResult["Amplitudes"];
LoadAmplitudes[subdirectory_String] := LoadResult["Amplitudes", subdirectory];

End[];
EndPackage[];
