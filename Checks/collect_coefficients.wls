(* collect_coefficients.wls: write the final coefficients of every process to
   Coefficients/, in a form that can be used without any of the algebra tools.

     wolframscript -file Checks/collect_coefficients.wls

   To refresh only the public file/input manifest of the existing tables:

     wolframscript -file Checks/collect_coefficients.wls manifest-only

   It reads the results that the process notebooks saved (Processes/<P>/Results/
   <route>/...), and for every coefficient
     1. checks that the three routes (FeynCalc, FORM, Tracer) give the same value
        at the six test points below;
     2. rewrites it in a canonical form:
          - every logarithm is split into the logarithms of the positive factors of
            its argument (in the physical region 0 < v, w < 1, s > 0);
          - the constants EulerGamma, log 2, log 3 and log pi, which must cancel in
            MSbar coefficients, are checked to cancel (at the six test points) and
            removed;
          - the scales appear only through LR = log(muR^2/s), LF = log(muF^2/s),
            LFF = log(muFF^2/s) (and LI, LO, LL for muI, muO, muL in Collins);
            the leftover log s is checked to cancel;
          - colour factors are written in Nc (CA = Nc, CF = (Nc^2 - 1)/(2 Nc));
          - the result is collected in the logarithms, with factored coefficients;
     3. checks that the canonical form equals the saved coefficient at the six
        test points;
   and writes one file per process.  Nothing in Processes/ is modified.
   The test points are random rational points of the physical region with
   random scales and N_f, three with N_c = 3 and three with N_c = 5, evaluated
   with 40-digit numerics: checks 1-3 are sampled numerical comparisons (the
   difference must be below 10^-30 at every point), not symbolic identities.
   The colour dependence enters through N_c (C_A = N_c, C_F = (N_c^2-1)/(2 N_c)),
   so two values of N_c test it beyond the physical SU(3).

   The transversely polarized proton-proton processes (Processes/HadronicTT and
   Processes/HadronicTransfer) are collected class by class: a class is written
   only when its final file (Transverse.wl, Transfer.wl) exists for the
   FeynCalc route and every other route is either complete as well or
   certified at step 2.  Steps 3 and 4 are route-independent (they act on the
   step-2 output Algebra.wl, which the FORM and Tracer files of step 2 compare
   at exact-rational sample points with the FeynCalc one). When a route has
   its own final file, that finite result is compared directly.  A route whose final file is missing is therefore
   accepted when its step-2 Algebra.wl exists: the comparison of the end of
   step 2 is repeated here (Born, Virtual, Real for every key and every entry
   of the Born table, exactly at two random rational points of s, t, u, A, B,
   D, CA, spe3, Nf, gs, with every scalar integral replaced by its own random
   rational number and HatPower[n] by hp^n), and the class records
   "StepsThreeFour" -> "run once (FeynCalc)" and, per route, what was used
   ("Routes").  The classes still missing are printed and listed in the file
   under "NotYetCollected", so the script can simply be rerun when they are
   done.
   To rewrite only some of the files, name them in the environment variable
   COLLECT, e.g.

     COLLECT=ProtonProtonTT,ProtonProtonTransfer wolframscript -file Checks/collect_coefficients.wls *)

$HistoryLength = 0;
root = DirectoryName[ExpandFileName[$InputFileName], 2];
Get[FileNameJoin[{root, "Checks", "ScriptArguments.wl"}]];
processes = FileNameJoin[{root, "Processes"}];
outDir = FileNameJoin[{root, "Coefficients"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir]];
routes = {"FeynCalc", "FORM", "Tracer"};
currentName = "";

(* The public manifest records table hashes, repository-relative input hashes
   and the distributed source hashes. Execution receipts remain with the local
   process results; they are not embedded in the distributed manifest. *)
sha[f_] := IntegerString[FileHash[f, "SHA256"], 16, 64];
rootPrefix = StringReplace[ExpandFileName[root], "/" ~~ EndOfString -> ""] <> "/";
relative[f_] := StringReplace[ExpandFileName[f], StartOfString ~~ rootPrefix -> ""];
readLog = <||>; writtenTables = {};
Unprotect[Get, Put];
Get[f_String, rest___] /; !TrueQ[$inLog] := Block[{$inLog = True},
  If[StringStartsQ[ExpandFileName[f], processes] && StringContainsQ[f, "/Results/"], readLog[relative[f]] = sha[f]]; Get[f, rest]];
Put[args__, f_String] /; !TrueQ[$inLog] && StringStartsQ[ExpandFileName[f], outDir] := Block[{$inLog = True},
  With[{r = Put[args, f]}, AppendTo[writtenTables, f];
    manifest[FileNameTake[f]] = <|"sha256" -> sha[f], "inputs" -> readLog|>; readLog = <||>; r]];
Protect[Get, Put];
manifestFile = FileNameJoin[{outDir, "MANIFEST.json"}];
manifest = If[FileExistsQ[manifestFile], Lookup[Import[manifestFile, "RawJSON"], "files", <||>], <||>];
(* COLLECT=<file>,<file>,... restricts the run to these output files (names without .wl); default: all *)
manifestOnly = MemberQ[ScriptArguments[], "manifest-only"];
selected = If[manifestOnly, {}, With[{e = Environment["COLLECT"]},
  If[StringQ[e] && StringTrim[e] =!= "", StringTrim /@ StringSplit[e, ","], All]]];
want[name_] := selected === All || MemberQ[selected, name];
fail[msg_, detail_: Null] := (Print["FAILED (", currentName, "): ", msg]; If[detail =!= Null, Print["  ", Short[detail, 10]]]; Quit[1]);

(* ---- physical test points (0 < v, w < 1, positive scales): three with Nc = 3, three with Nc = 5 ---- *)
SeedRandom[20261011];
testPoints = Join @@ Table[Table[<|v -> RandomInteger[{5, 95}]/100, w -> RandomInteger[{5, 95}]/100, s -> RandomInteger[{3, 40}],
   muR -> RandomInteger[{1, 9}]/3, muF -> RandomInteger[{1, 9}]/4, muFF -> RandomInteger[{1, 9}]/5,
   muI -> RandomInteger[{1, 9}]/6, muO -> RandomInteger[{1, 9}]/7, muL -> RandomInteger[{1, 9}]/2,
   Nf -> RandomInteger[{1, 6}], CA -> nc, CF -> (nc^2 - 1)/(2 nc), Nc -> nc, TR -> 1/2|>, {3}], {nc, {3, 5}}];
reference = <|v -> 3/10, w -> 1/5, s -> 7, muR -> 2/3, muF -> 3/4, muFF -> 4/5, muI -> 1/2, muO -> 5/7, muL -> 3/2|>;
value[expr_, point_] := Quiet[N[expr /. Normal[point], 40], N::meprec];   (* an exact cancellation cannot reach 40 relative digits *)

(* ---- 1. logarithms of positive factors ---- *)
scaleSymbols = {muR, muF, muFF, muI, muO, muL};
logOf[arg_] := Module[{t = Together[arg], parts, terms, sign},
  parts = Join[FactorList[Numerator[t]], {#[[1]], -#[[2]]} & /@ FactorList[Denominator[t]]];
  sign = 1;
  terms = Table[With[{f = part[[1]], n = part[[2]]},
      Which[
        IntegerQ[f] || Head[f] === Rational, (sign *= Sign[f]^n; n Total[#[[2]] Log[#[[1]]] & /@ FactorInteger[Abs[f]]]),
        f === Pi, n Log[Pi],
        True, With[{g = If[(f /. reference) > 0, Expand[f], sign *= (-1)^n; Expand[-f]]}, n Log[g]]]],
    {part, parts}];
  If[sign =!= 1, fail["a logarithm has a negative argument in the physical region", arg]];
  Total[terms]];
splitLogs[expr_] := expr /. Log[a_] :> logOf[a];

(* ---- 2. constants and scales ---- *)
scaleLog = <|muR -> LR, muF -> LF, muFF -> LFF, muI -> LI, muO -> LO, muL -> LL|>;
logValues[pt_] := Normal[Association[KeyValueMap[#2 -> Log[(#1 /. Normal[pt])^2/(s /. Normal[pt])] &, scaleLog]]];
canonical[expr0_] := Module[{expr = splitLogs[expr0], probe},
  expr = expr /. {EulerGamma -> cEuler, Log[2] -> cLog2, Log[3] -> cLog3, Log[Pi] -> cLogPi, Log[s] -> cLogS} /.
    Log[m_Symbol] /; MemberQ[scaleSymbols, m] :> (scaleLog[m] + cLogS)/2;
  Do[If[!FreeQ[expr, c],
      probe = Table[value[D[expr, c] /. {cEuler -> 0, cLog2 -> 0, cLog3 -> 0, cLogPi -> 0, cLogS -> 0} /. logValues[pt], pt], {pt, testPoints}];
      If[Max[Abs[probe]] > 10^-30, fail["the constant " <> ToString[c] <> " does not cancel", probe]]],
    {c, {cEuler, cLog2, cLog3, cLogPi, cLogS}}];
  expr = expr /. {cEuler -> 0, cLog2 -> 0, cLog3 -> 0, cLogPi -> 0, cLogS -> 0};
  expr = expr /. {CF -> (Nc^2 - 1)/(2 Nc), CA -> Nc};
  Collect[Expand[expr], Union[Cases[expr, (Log | PolyLog)[__], Infinity], Values[scaleLog]], Factor]];

(* ---- 3. check against the saved form ---- *)
checkCanonical[original_, new_] := Module[{d},
  d = Table[With[{a = value[original /. {CF -> (Nc^2 - 1)/(2 Nc), CA -> Nc}, pt], b = value[new /. logValues[pt], pt]},
      If[Abs[a] < 10^-30, Abs[b], Abs[a - b]/Abs[a]]], {pt, testPoints}];
  If[Max[d] > 10^-30, fail["the canonical form differs from the saved coefficient", d]];
  new];
sameInRoutes[name_, values_List] := Module[{d},
  (* every route minus the first, coefficient by coefficient (values is a list of routes, each a list of four) *)
  d = Table[Max[Abs[Flatten[value[(# - values[[1]]) /. {CF -> (Nc^2 - 1)/(2 Nc), CA -> Nc}, pt] & /@ values]]], {pt, testPoints}];
  If[!VectorQ[d, NumericQ], fail["the route comparison did not evaluate to numbers for " <> name, d]];
  If[Max[d] > 10^-30, fail["the routes differ for " <> name, d]]];
(* a result without plus distributions (the off-diagonal lepton channels) is stored as its regular term alone: pad it to {regular, 0, 0, 0} *)
asFour[c_] := If[ListQ[c], c, {c, 0, 0, 0}];
canonicalList[name_, perRoute0_List] := With[{perRoute = asFour /@ perRoute0}, currentName = name; Print["  ", name]; sameInRoutes[name, perRoute];
  Table[checkCanonical[perRoute[[1, k]], canonical[perRoute[[1, k]]]], {k, 4}]];
load[path_] := If[FileExistsQ[path], Get[path], fail["missing result file " <> path]];
basisText = "{regular, D0, D1, delta} with D0 = [1/(1-w)]_+, D1 = [log(1-w)/(1-w)]_+, delta = delta(1-w)";
scaleText = "LR = log(muR^2/s), LF = log(muF^2/s), LFF = log(muFF^2/s)";

(* ---- proton-proton classes ---- *)
If[want["ProtonProton"],
Get[FileNameJoin[{processes, "Hadronic", "classes.wl"}]];
pp = Association @ Table[Print["pp class ", n]; Module[{files, born},
    files = Table[load[FileNameJoin[{processes, "Hadronic", "Results", r, "class-" <> ToString[n], spinFile <> ".wl"}]], {spinFile, {"Unpolarized", "Polarized"}}, {r, routes}];
    n -> <|"Process" -> HadronicClasses[n]["Name"],
      "Born" -> <|"UU" -> Factor[files[[1, 1]]["Born"] /. {CF -> (Nc^2 - 1)/(2 Nc), CA -> Nc}], "LL" -> Factor[files[[2, 1]]["Born"] /. {CF -> (Nc^2 - 1)/(2 Nc), CA -> Nc}]|>,
      "UU" -> canonicalList["pp class " <> ToString[n] <> " UU", #["Coefficients"] & /@ files[[1]]],
      "LL" -> canonicalList["pp class " <> ToString[n] <> " LL", #["Coefficients"] & /@ files[[2]]]|>], {n, 1, 16}];
Put[<|"Description" -> "NLO hard coefficients of single-inclusive parton production a b -> c + X in proton-proton collisions, the 16 channel classes of Processes/Hadronic (see Coefficients/README.md)",
  "Cross section" -> "sigma = s E dsigma/d^3p = Born delta(1-w) + Prefactor (c_regular + c_D0 D0 + c_D1 D1 + c_delta delta), with the distributions' coefficients taken at w = 1 in the prefactor",
  "Prefactor" -> alphaS^3/(Pi s v w), "Basis" -> basisText, "Scales" -> scaleText,
  "Variables" -> "s = (pa+pb)^2, v = 1 + t/s, w = -u/(s+t); UU = spin average, LL = helicity difference (sigma(++) - sigma(+-))/2; colour SU(Nc)",
  "Classes" -> pp|>, FileNameJoin[{outDir, "ProtonProton.wl"}]]];

(* ---- lepton processes ---- *)
leptonFile[process_, file_] := load /@ Table[FileNameJoin[{processes, process, "Results", r, file}], {r, routes}];
If[want["LeptonQuark"],
lq = Association @ Table[spin -> With[{f = leptonFile["LeptonQuark", spin <> ".wl"]},
    <|"Born" -> Factor[f[[1]]["Born"]], "Prefactor" -> f[[1]]["Prefactor"], "Coefficients" -> canonicalList["lepton-quark " <> spin, #["Coefficients"] & /@ f]|>], {spin, {"Unpolarized", "Polarized"}}];
Put[<|"Description" -> "NLO hard coefficients of l q -> q + X (tagged quark), Processes/LeptonQuark", "Basis" -> basisText,
  "Scales" -> "LR = log(muR^2/s) (one common scale for renormalization, the incoming PDF and the fragmentation)",
  "Cross section" -> "sigma = s E dsigma/d^3p = Born delta(1-w) + Prefactor (c_regular + c_D0 D0 + c_D1 D1 + c_delta delta), distributions' coefficients with the prefactor at w = 1",
  "Variables" -> "s = (l+k)^2, v = 1 + t/s, w = -u/(s+t); Unpolarized = spin average, Polarized = helicity difference",
  "Results" -> lq|>, FileNameJoin[{outDir, "LeptonQuark.wl"}]]];
If[want["LeptonOffDiagonal"],
lod = Association @ Table[spin -> With[{f = leptonFile["LeptonOffDiagonal", spin <> ".wl"]},
    Association @ Table[ch -> <|"Prefactor" -> f[[1]]["Prefactor"][ch], "Coefficients" -> canonicalList["off-diagonal " <> ch <> " " <> spin, (#["Coefficients"][ch] & /@ f)]|>, {ch, {"QG", "GQ"}}]],
  {spin, {"Unpolarized", "Polarized"}}];
Put[<|"Description" -> "NLO hard coefficients of the off-diagonal lepton channels l q -> g + X (QG) and l g -> q + X (GQ), Processes/LeptonOffDiagonal; they start at order alphaS, without Born term",
  "Basis" -> basisText, "Scales" -> "LR = log(muR^2/s) (common scale)", "Results" -> lod|>, FileNameJoin[{outDir, "LeptonOffDiagonal.wl"}]]];
If[want["Collins"],
collins = Association @ Table[part -> With[{f = leptonFile["Collins", part <> ".wl"]},
    <|"Born" -> Factor[f[[1]]["Born"]], "Prefactor" -> f[[1]]["Prefactor"], "Coefficients" -> canonicalList["Collins " <> part, #["Coefficients"] & /@ f]|>], {part, {"Unpolarized", "Transverse"}}];
Put[<|"Description" -> "NLO hard coefficients of l q -> q + X with transverse spin transfer from the incoming to the tagged quark (Transverse), and the spin average (Unpolarized), Processes/Collins",
  "Basis" -> basisText, "Scales" -> "LR = log(muR^2/s); LI = log(muI^2/s), LO = log(muO^2/s), LL = log(muL^2/s) for the incoming quark, the outgoing quark and the photon-in-lepton subtraction",
  "Results" -> collins, "ResolvedBorn" -> load[FileNameJoin[{processes, "Collins", "Results", "ResolvedBorn.wl"}]]|>, FileNameJoin[{outDir, "Collins.wl"}]]];

(* ---- TMD and Born tables: already in final variables; the routes must agree ---- *)
If[want["TMD"],
tmd = Table[Join[load[FileNameJoin[{processes, "TMD", "Results", r, "Beam.wl"}]], load[FileNameJoin[{processes, "TMD", "Results", r, "Soft.wl"}]]], {r, routes}];
If[!And @@ Flatten[Table[
    Keys[tmd[[1]][k]] === Keys[tmd[[j]][k]] &&
      Union[Simplify[Values[tmd[[1]][k]] - Values[tmd[[j]][k]]]] === {0},
    {j, 2, Length[routes]}, {k, {"Beam", "TMDPDF", "Soft"}}]], fail["the TMD routes differ"]];
Put[<|"Description" -> "NLO TMD beam functions, TMD PDFs and soft functions, Processes/TMD (notation in its README: as = alphaS/(4 Pi), Lp, Lv, LQ, Lh)",
  "Beam" -> tmd[[1]]["Beam"], "TMDPDF" -> tmd[[1]]["TMDPDF"], "Soft" -> tmd[[1]]["Soft"]|>, FileNameJoin[{outDir, "TMD.wl"}]]];
If[want["PartonicBorn"],
born = Get[FileNameJoin[{processes, "PartonicBorn", "Results", "Born.wl"}]];
Put[<|"Description" -> "D-dimensional Born squares |M|^2/gs^4 of all 2 -> 2 partonic channels, spin average (UU) and helicity difference (LL), Processes/PartonicBorn", "Born" -> born["Born"]|>,
  FileNameJoin[{outDir, "PartonicBorn.wl"}]]];

(* ---- transversely polarized proton-proton collisions ---- *)
(* a class enters when its final file exists in the FeynCalc route and every other route has either its
   final file or an identical step-2 result (see the header); the routes with a final file must agree on
   the Born term (exactly) and on every coefficient (numerically, as above) *)
colourToNc = {CF -> (Nc^2 - 1)/(2 Nc), CA -> Nc};
classFile[process_, route_, n_, file_] := FileNameJoin[{processes, process, "Results", route, "class-" <> ToString[n], file}];
(* the comparison of the end of the FORM and Tracer files of step 2: two step-2 results (Algebra.wl) must be
   identical; exact, at two random rational points (s > 0; t, u, A, B < 0; C = -(s+t+u) - A - B;
   CF = (CA^2 - 1)/(2 CA)), every scalar integral K...[n1, ..., n4] replaced by its own random rational
   number and HatPower[n] by hp^n, so that the coefficient of every integral is compared *)
sameAlgebra[name_, fileA_, fileB_] := Module[{a = Get[fileA], b = Get[fileB], integrals, values, points, evaluate, parts, bad},
  currentName = name;
  If[Keys[a["BornTable"]] =!= Keys[b["BornTable"]] || Or @@ Table[Keys[a["BornTable"][k]] =!= Keys[b["BornTable"][k]], {k, Keys[a["BornTable"]]}],
    fail["the two step-2 results have different Born tables"]];
  Do[If[Keys[a[part]] =!= Keys[b[part]], fail["the two step-2 results have different keys in " <> part]], {part, {"Born", "Virtual", "Real"}}];
  integrals = Union[Cases[{a["Virtual"], b["Virtual"]}, (h_Symbol)[__Integer] /; StringStartsQ[SymbolName[h], "KVirtual"], Infinity]];
  BlockRandom[SeedRandom[20261013];
    values = Table[RandomInteger[{1, 10^9}, Length[integrals]]/RandomInteger[{1, 10^9}, Length[integrals]], {2}];
    If[!And @@ (DuplicateFreeQ /@ values), fail["the random values of the integrals are not distinct"]];
    points = Table[Join[Thread[{s, t, u, A, B, D, CA, spe3, Nf, gs, hp} -> RandomInteger[{1, 10^6}, 11]/RandomInteger[{1, 10^6}, 11] {1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1}],
      Thread[integrals -> values[[k]]]], {k, 2}]];
  evaluate[x_, point_] := Together[x /. HatPower[n_] :> hp^n /. point /. C -> -(s + t + u) - A - B /. CF -> (CA^2 - 1)/(2 CA) /. point];
  parts = Join[Flatten[Table[{part, key} -> {a[part][key], b[part][key]}, {part, {"Born", "Virtual", "Real"}}, {key, Keys[a[part]]}]],
    Flatten[Table[{"BornTable", key, ch} -> {a["BornTable"][key][ch], b["BornTable"][key][ch]}, {key, Keys[a["BornTable"]]}, {ch, Keys[a["BornTable"][key]]}]]];
  bad = Select[parts, !(And @@ Table[evaluate[#[[2, 1]] - #[[2, 2]], point] === 0, {point, points}]) &];
  If[bad === {}, Print["    ", name, ": step-2 results identical (", Length[parts], " expressions, 2 exact points)"]; Return[True, Module]];
  (* the one accepted difference: a virtual interference may be written differently in the two routes, by a multiple
     of the four-dimensional Gram identity of the loop momentum, which integrates to zero (see the comparison cell of
     2-Algebra-FORM.wl of HadronicTransfer).  The step-2 file compares the two routes after the loop-angle projection
     and stops when they differ; it is accepted here only when that file completed (its receipt) after this
     Algebra.wl was written *)
  If[!MatchQ[bad[[All, 1]], {{"Virtual", _} ..}] || !receiptCertified[fileA],
    fail["the step-2 results differ", bad[[All, 1]]]];
  Print["    ", name, ": step-2 results identical (", Length[parts] - Length[bad], " expressions, 2 exact points); ",
    "the virtual of ", StringRiffle[bad[[All, 1, 2]], ", "], " identical after the loop-angle projection (required by the completed step-2 file)"];
  True];
(* a step-2 result whose own step file completed after it was written: the receipt (Checks/run_process.py) lists the
   class among the overrides, says "completed": true, and is newer than the Algebra.wl file *)
receiptCertified[algebraFile_] := Module[{route = FileNameTake[algebraFile, {-3}], n = StringDrop[FileNameTake[algebraFile, {-2}], 6],
    results = FileNameDrop[algebraFile, -3], receipts},
  receipts = Select[FileNames["2-Algebra-" <> route <> "-*.receipt.json", results], Module[{r = Import[#, "RawJSON"]},
      TrueQ[r["completed"]] && MemberQ[r["overrides"], {"class", n}] && MemberQ[r["overrides"], {"route", route}] &&
        AbsoluteTime[FileDate[#]] > AbsoluteTime[FileDate[algebraFile]]] &];
  receipts =!= {}];
(* per class and route: "steps 2-4" (final file there), "step 2" (final file missing, step-2 result there and
   identical to the FeynCalc one: steps 3-4 run once), or missing *)
routeStatus[process_, file_, n_] := AssociationMap[Which[
    FileExistsQ[classFile[process, #, n, file]], "steps 2-4",
    # =!= "FeynCalc" && FileExistsQ[classFile[process, "FeynCalc", n, file]] && FileExistsQ[classFile[process, #, n, "Algebra.wl"]] && FileExistsQ[classFile[process, "FeynCalc", n, "Algebra.wl"]], "step 2",
    True, "missing"] &, routes];
completeClasses[process_, file_, all_List] := Module[{status, done, missing},
  status = AssociationMap[routeStatus[process, file, #] &, all];
  done = Select[all, FreeQ[Values[status[#]], "missing"] &];
  missing = Complement[all, done];
  Do[Print["  ", process, " class ", n, ": not collected, ", file, " missing for ",
     StringRiffle[Keys[Select[status[n], # === "missing" &]], ", "], " (rerun this script when it is there)"], {n, missing}];
  {done, missing, status}];
(* the routes whose final file is used, after the step-2 comparison of the others; the record of the class *)
usedRoutes[process_, n_, status_] := (
  Do[If[status[r] === "step 2", sameAlgebra[process <> " class " <> ToString[n] <> " " <> r <> " vs FeynCalc (step 2)",
      classFile[process, r, n, "Algebra.wl"], classFile[process, "FeynCalc", n, "Algebra.wl"]]], {r, routes}];
  Keys[Select[status, # === "steps 2-4" &]]);
routeRecord[status_] := <|"Routes" -> Map[Replace[{"steps 2-4" -> "steps 2-4", "step 2" -> "step 2 only: Algebra.wl identical to the FeynCalc route (compared in exact rational arithmetic at two random points; a virtual may differ in representation and is then compared after the loop-angle projection, in the step-2 file)"}], status],
   "StepsThreeFour" -> If[MemberQ[Values[status], "step 2"], "run once (FeynCalc): steps 3-4 act on the step-2 result, which is identical in the routes listed as step 2 only", "run in every route"]|>;
sameBorn[name_, borns_List] := (currentName = name;
  If[!(And @@ (Simplify[(# - borns[[1]]) /. colourToNc] === 0 & /@ borns)), fail["the routes give different Born terms", borns]];
  Factor[borns[[1]] /. colourToNc]);
ttScaleText = "LR = log(muR^2/s), LF = log(muF^2/s) (both transversity densities), LFF = log(muFF^2/s) (fragmentation)";

If[want["ProtonProtonTT"],
Get[FileNameJoin[{processes, "HadronicTT", "classes.wl"}]];
{ttDone, ttMissing, ttStatus} = completeClasses["HadronicTT", "Transverse.wl", Keys[HadronicTTClasses]];
tt = Association @ Table[Print["TT class ", n]; Module[{files = load /@ Table[classFile["HadronicTT", r, n, "Transverse.wl"], {r, usedRoutes["HadronicTT", n, ttStatus[n]]}]},
    If[Union[#["Spin"] & /@ files] =!= {"TT"}, fail["unexpected spin entry in Transverse.wl of class " <> ToString[n]]];
    n -> Join[<|"Process" -> HadronicTTClasses[n]["Name"],
      "Born" -> sameBorn["TT class " <> ToString[n] <> " Born", #["Born"] & /@ files],
      "TT" -> canonicalList["TT class " <> ToString[n], #["Coefficients"] & /@ files]|>, routeRecord[ttStatus[n]]]], {n, ttDone}];
Put[<|"Description" -> "NLO hard coefficients of single-inclusive parton production a b -> c + X with both incoming partons transversely polarized (the double-spin asymmetry A_TT), Processes/HadronicTT: the coefficient of cos(2 Phi) (see Coefficients/README.md and docs/TRANSVERSE-PP.md)",
  "Cross section" -> "delta sigma_TT = s E d delta sigma/d^3p projected on cos(2 Phi) = Born delta(1-w) + Prefactor (c_regular + c_D0 D0 + c_D1 D1 + c_delta delta), with the distributions' coefficients taken at w = 1 in the prefactor",
  "Prefactor" -> alphaS^3/(Pi s v w), "Basis" -> basisText, "Scales" -> ttScaleText,
  "Variables" -> "s = (pa+pb)^2, v = 1 + t/s, w = -u/(s+t); TT = (1/2)[sigma(up,up) - sigma(up,down)], the coefficient of cos(2 Phi), Phi the azimuth of p about the common spin axis; HVBM gamma5, MSbar; colour SU(Nc)",
  "Classes" -> tt, "NotYetCollected" -> ttMissing|>, FileNameJoin[{outDir, "ProtonProtonTT.wl"}]]];

If[want["ProtonProtonTransfer"],
Get[FileNameJoin[{processes, "HadronicTransfer", "classes.wl"}]];
{trDone, trMissing, trStatus} = completeClasses["HadronicTransfer", "Transfer.wl", Keys[HadronicTransferClasses]];
tr = Association @ Table[Print["transfer class ", n]; Module[{files = load /@ Table[classFile["HadronicTransfer", r, n, "Transfer.wl"], {r, usedRoutes["HadronicTransfer", n, trStatus[n]]}]},
    n -> Join[<|"Process" -> HadronicTransferClasses[n]["Name"],
      "Born" -> Association @ Table[comp -> sameBorn["transfer class " <> ToString[n] <> " Born " <> comp, #["Born"][comp] & /@ files], {comp, {"XX", "NN"}}]|>,
      Association @ Table[comp -> canonicalList["transfer class " <> ToString[n] <> " " <> comp, #["Coefficients"][comp] & /@ files], {comp, {"XX", "NN"}}],
      routeRecord[trStatus[n]]]], {n, trDone}];
Put[<|"Description" -> "NLO hard coefficients of the transverse spin transfer a^ b -> c^ + X in proton-proton collisions (incoming parton a and observed parton c transversely polarized), Processes/HadronicTransfer: the components XX (both spins in the scattering plane) and NN (both along the normal)",
  "Cross section" -> "sigma_comp = s E d sigma/d^3p, component comp = XX or NN, = Born delta(1-w) + Prefactor (c_regular + c_D0 D0 + c_D1 D1 + c_delta delta), with the distributions' coefficients taken at w = 1 in the prefactor; sigma(s_a, s_c) = XX cos(phi_a) cos(phi_c) + NN sin(phi_a) sin(phi_c)",
  "Prefactor" -> alphaS^3/(Pi s v w), "Basis" -> basisText,
  "Scales" -> "LR = log(muR^2/s), LF = log(muF^2/s) (both incoming partons), LFF = log(muFF^2/s) (fragmentation of the polarized parton c)",
  "Variables" -> "s = (pa+pb)^2, v = 1 + t/s, w = -u/(s+t); density pa-slash gamma5 sa-slash/2 for a, projection p-slash gamma5 sc-slash for c (sigma(s_c) - sigma(-s_c)); HVBM gamma5, MSbar; colour SU(Nc).  The finite XX and NN coefficients are equal in every class (in this form, after CF -> (Nc^2-1)/(2 Nc)); both are kept",
  "Classes" -> tr, "NotYetCollected" -> trMissing|>, FileNameJoin[{outDir, "ProtonProtonTransfer.wl"}]]];

If[writtenTables =!= {}, Print["written: ", writtenTables]];

(* Only source files enter this snapshot, never generated Results/ trees.
   The input hashes describe the calculation products used for each table;
   the sources block describes the files distributed with those tables. *)
sourceExtensions = {"wl", "wls", "m", "mod", "gen", "frm", "h", "yaml", "py"};
sourceFiles[dir_] := Module[{entries = FileNames["*", dir], files, folders},
  files = Select[entries, !DirectoryQ[#] && MemberQ[sourceExtensions, FileExtension[#]] &];
  folders = Select[entries, DirectoryQ[#] && !StringStartsQ[FileNameTake[#], "."] &&
    !MemberQ[{"Results", "__pycache__"}, FileNameTake[#]] &];
  Join[files, Flatten[sourceFiles /@ folders]]];
distributedSources = Sort @ Flatten[sourceFiles /@ (FileNameJoin[{root, #}] & /@ {"FeynFacet", "Processes", "Checks"})];
manifest = Association @ KeyValueMap[Function[{name, record},
  If[!FileExistsQ[FileNameJoin[{outDir, name}]] || sha[FileNameJoin[{outDir, name}]] =!= record["sha256"],
    fail["Coefficient file differs from its calculation-input binding", name]];
  name -> <|"sha256" -> record["sha256"], "inputs" -> KeySort @ Association @
    KeyValueMap[Function[{path, value}, path -> <|"sha256" -> If[AssociationQ[value], value["sha256"], value]|>],
      record["inputs"]]|>], manifest];
Block[{$inLog = True}, Export[manifestFile,
  <|"description" -> "SHA-256 bindings of the distributed coefficient tables, scientific source files, and relative calculation inputs. Intermediate calculation products are generated by running the numbered workflows and are not distributed.",
    "generator" -> <|"file" -> "Checks/collect_coefficients.wls", "sha256" -> sha[FileNameJoin[{root, "Checks", "collect_coefficients.wls"}]]|>,
    "sources" -> Association @ Table[relative[file] -> sha[file], {file, distributedSources}],
    "files" -> KeySort[manifest]|>, "RawJSON"]];
Print["manifest: ", Length[manifest], " coefficient tables, ",
  Total[Length /@ Values[manifest[[All, "inputs"]]]], " calculation inputs, ",
  Length[distributedSources], " distributed source files"];
