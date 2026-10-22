(* Compare the independently assembled direct and crossed transfer results.
   Both inputs are NEW independent.wl outputs, never production Transfer.wl.
   No production crossing, distribution, or factorization helper is loaded. *)
$HistoryLength = 0;
directDirectory = Environment["INDEPENDENT_DIRECT_RESULTS"];
crossedDirectory = Environment["INDEPENDENT_CROSSED_RESULTS"];
outputDirectory = Environment["INDEPENDENT_CROSSING_OUTPUT"];
If[!And @@ (StringQ /@ {directDirectory, crossedDirectory, outputDirectory}),
  Print["Set INDEPENDENT_DIRECT_RESULTS, INDEPENDENT_CROSSED_RESULTS and INDEPENDENT_CROSSING_OUTPUT."]; Quit[2]];
If[!DirectoryQ[outputDirectory], CreateDirectory[outputDirectory, CreateIntermediateDirectories -> True]];
assumptions = s > 0 && 0 < v < 1 && 0 < w < 1 && CA > 1 && Nf > 0 &&
  muR > 0 && muF > 0 && muFF > 0;
(* Expand logarithms only when every real rational factor has a certified
   nonzero sign and the complete argument is positive. Unknown and complex
   factors preserve the original logarithm, including its phase. *)
Clear[certifiedFactorSign, positiveRationalLog, certifiedLogOf, certifiedLogs];
certifiedFactorSign[factor_, assumptions_] := certifiedFactorSign[factor, assumptions] = Which[
  !TrueQ[FullSimplify[Element[factor, Reals], assumptions]], 0,
  TrueQ[FullSimplify[factor > 0, assumptions]], 1,
  TrueQ[FullSimplify[factor < 0, assumptions]], -1,
  True, 0];

positiveRationalLog[value_] := If[(IntegerQ[value] || Head[value] === Rational) && value > 0,
  Total[#[[2]] Log[#[[1]]] & /@ FactorInteger[Numerator[value]]] -
    Total[#[[2]] Log[#[[1]]] & /@ FactorInteger[Denominator[value]]], Log[value]];

certifiedLogOf[argument_, assumptions_] := certifiedLogOf[argument, assumptions] = Module[
  {piPlaceholder, factored, factor, exponent, sign, productSign = 1, terms = {}, valid = True},
  (* Treat Pi as a formal positive factor, so a numeric coefficient such as
     16 Pi^2/9 is decomposed rather than left as one indivisible number. *)
  factored = Quiet[Check[FactorList[Together[argument /. Pi -> piPlaceholder]], $Failed]];
  If[factored === $Failed || !ListQ[factored] ||
    !AllTrue[factored, MatchQ[#, {_, _Integer}] &], Return[Log[argument]]];
  Do[
    factor = pair[[1]] /. piPlaceholder -> Pi;
    exponent = pair[[2]];
    sign = certifiedFactorSign[factor, assumptions];
    If[sign === 0, valid = False; Break[]];
    If[sign === -1, factor = Expand[-factor]; productSign *= (-1)^exponent];
    AppendTo[terms, exponent positiveRationalLog[factor]], {pair, factored}];
  If[!valid, Log[argument], If[productSign === 1, Total[terms], Log[argument]]]];

certifiedLogs[expression_, assumptions_] :=
  expression /. Log[argument_] :> certifiedLogOf[argument, assumptions];

certifiedZeroQ[expression_, assumptions_] := Module[{rational},
  rational = Together[Expand[certifiedLogs[expression, assumptions]]];
  TrueQ[rational === 0 || FullSimplify[rational, assumptions] === 0]];

exactZero[expression_] := certifiedZeroQ[expression /. CF -> (CA^2 - 1)/(2 CA), assumptions];
rows = {}; inputs = <||>; differences = <||>;
started = AbsoluteTime[];
Do[
  directFile = FileNameJoin[{directDirectory, "class-" <> number, "independent.wl"}];
  crossedFile = FileNameJoin[{crossedDirectory, "class-" <> number <> "x", "independent.wl"}];
  Do[inputs[file] = IntegerString[FileHash[file, "SHA256"], 16, 64], {file, {directFile, crossedFile}}];
  direct = Get[directFile]; crossed = Get[crossedFile];
  If[direct["Failures"] =!= {} || crossed["Failures"] =!= {} ||
    Sort[Keys[direct["Assembly"]]] =!= {"NN", "XX"} ||
    Sort[Keys[crossed["Assembly"]]] =!= {"NN", "XX"},
    Print["An input independent assembly did not pass in full: class ", number]; Quit[1]];
  componentDifferences = <||>;
  Do[
    difference = direct["Assembly"][spin]["Finite"] - crossed["Assembly"][spin]["Finite"];
    componentDifferences[spin] = difference;
    Do[
      phaseStarted = AbsoluteTime[];
      passed = exactZero[difference[[slot]]];
      (* These finite changes leave poles and explicit scale derivatives alone. *)
      perturbation = alphaS^3/(Pi s v) v (1 - v) If[slot === 1, w (1 - w), 1];
      controlPassed = !exactZero[difference[[slot]] + perturbation];
      AppendTo[rows, <|"class" -> number, "spin" -> spin, "slot" -> slot,
        "exactIdentity" -> passed, "nonzeroControlDetected" -> controlPassed,
        "seconds" -> N[AbsoluteTime[] - phaseStarted]|>];
      Print["CROSSING ", number, " ", spin, " ", slot, " exact=", passed,
        " control=", controlPassed], {slot, 4}], {spin, {"XX", "NN"}}];
  differences[number] = componentDifferences,
  {number, {"1", "3", "6", "10", "11", "13", "17"}}];
passed = And @@ (TrueQ[#["exactIdentity"]] && TrueQ[#["nonzeroControlDetected"]] & /@ rows);
Put[<|"Differences" -> differences, "Checks" -> rows, "Inputs" -> inputs|>,
  FileNameJoin[{outputDirectory, "crossing.wl"}]];
Export[FileNameJoin[{outputDirectory, "summary.json"}], <|"passed" -> passed,
  "exactIdentities" -> Count[Lookup[rows, "exactIdentity"], True],
  "detectedControls" -> Count[Lookup[rows, "nonzeroControlDetected"], True],
  "checks" -> rows, "inputs" -> inputs, "seconds" -> N[AbsoluteTime[] - started],
  "checkerSHA256" -> IntegerString[FileHash[$InputFileName, "SHA256"], 16, 64]|>, "RawJSON"];
Quit[If[passed, 0, 1]];
