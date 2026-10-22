(* Analytic q g -> q endpoint benchmark.
   F. Aversa et al., Nucl. Phys. B327 (1989) 105-143,
   Eqs. (22), (23), (41), (43), printed pp. 113-114 and 127-128.
   Eq. (43) is used at CQ=0; its quark finite PDF/FF terms are retained
   by the paper and must be converted to the MSbar convention of the tables.
   Only Coefficients/ProtonProton.wl is loaded; no process code is called. *)
$HistoryLength = 0;
sourceRoot = Environment["FEYNFACET_SOURCE"];
outputRoot = Environment["FEYNFACET_BENCHMARK_OUTPUT"];
If[!And @@ (StringQ /@ {sourceRoot, outputRoot}),
  Print["Set FEYNFACET_SOURCE and FEYNFACET_BENCHMARK_OUTPUT."]; Quit[2]];
coefficientFile = FileNameJoin[{sourceRoot, "Coefficients", "ProtonProton.wl"}];
inputHash = IntegerString[FileHash[coefficientFile, "SHA256"], 16, 64];
coefficients = Get[coefficientFile];
class = coefficients["Classes"][13];
checks = {}; failures = {};
check[name_, condition_] := (AppendTo[checks, <|"name" -> name, "passed" -> TrueQ[condition]|>];
  If[!TrueQ[condition], AppendTo[failures, name]; Print["FAIL ", name]]; TrueQ[condition]);
(* The endpoint expressions contain only Log[v], Log[1-v] and the three
   independent scale logs. Treating these as algebraically independent
   indeterminates establishes an identity without numerical sampling or
   logarithmic continuation. The physical region is 0<v<1. *)
canonical[expression_] := Factor[Together[expression /.
  {Log[v] -> logV, Log[1 - v] -> logOneMinusV}]];
zeroQ[expression_] := TrueQ[canonical[expression] === 0];
check["tagged-channel", class["Process"] === "q g -> q"];
check["exact-input", FreeQ[{class["Born"]["UU"], class["UU"][[2 ;; 4]]}, _Real]];

(* Born normalization: sigma=s E dsigma/d^3p=(1/(pi v)) dsigma/dv.
   The coefficient convention is sigma_NLO=alphaS^3/(pi s v) times
   {cD0,cD1,cDelta}; Eq. (22) factors Born*alphaS/(2pi). *)
p = Nc^2 (1 + v^2) - (1 - v)^2;
publishedBorn = Pi alphaS^2 (1 + v^2) p/(2 Nc^2 s v (1 - v)^2);
check["Born/Eq41", zeroQ[Pi v class["Born"]["UU"] - publishedBorn]];
bornUnit = canonical[s class["Born"]["UU"]/alphaS^2];
normalizedMSbar = canonical /@ (2 class["UU"][[2 ;; 4]]/(v bornUnit));

(* Eq. (23): fqq=CF[(1+z^2)D1-3/2 D0+r_f+d_f delta],
   dqq=CF[(1+z^2)D1-3/2 D0+r_d+d_d delta], with
   d_f=-9/2-pi^2/3, d_d=2pi^2/3-9/2. Regular r_f,r_d contribute no
   endpoint distributions. Their kernel maps have slopes 1 and v.
   Pullback of j D1 gives endpoint rows {j1/k log k, j1/k,
   j1/(2k) log^2 k}; that of j D0 gives {j1/k,0,j1/k log k}.
   j1/k equals Born on each leg; the (1+z^2) factor supplies two.
   Consequently the summed finite quark-scheme shift, in the normalized
   {D0,D1,delta} basis, is the expression below. Its subtraction converts
   the released MSbar coefficient to Eq. (43)'s CQ=0 convention. *)
quarkCF = (Nc^2 - 1)/(2 Nc);
quarkSchemeShift = quarkCF {2 Log[v] - 3, 4,
  Log[v]^2 - 3 Log[v]/2 + Pi^2/3 - 9};
normalizedAversa = canonical /@ (normalizedMSbar - quarkSchemeShift);

(* Literal Eq. (43) transcription at CQ=0, grouped by powers of the two
   kinematic logarithms. Tr=Nf/2. The printed logs ln(s/M^2),
   ln(s/M_f^2), ln(s/mu^2) equal -LF, -LFF, -LR, respectively.
   The paper prints 2 C_i; the expressions below are C_i. *)
tr = Nf/2;
a = 2 Nc^2 (4 v^2 - 7 v + 4) - Nc^4 (3 v^2 + 7) - (v - 1)^2;
b = 2 Nc^2 (v^2 - 3 v + 1) - Nc^4 (v^2 + 1) + (v - 1)^2;
c = Nc^4 (-2 v^4 - 2 v^3 - 11 v^2 - 9) +
  Nc^2 (v - 1) (8 v^3 - 8 v^2 + 9 v - 10) - (v - 1)^2 (2 v^2 - 2 v + 1);
d = Nc^4 (13 v^4 + 38 v^2 + 6 v - 5) -
  2 Nc^2 (2 v^4 - 19 v^3 + v^2 + 14 v - 16) -
  3 (v - 1)^2 (3 v^2 - 4 v + 9) - 8 tr Nc^3 (v^2 + 1)^2 +
  8 tr Nc (v - 1)^2 (v^2 + 1);
e = Nc^2 (2 v^4 - 6 v^3 + 5 v^2 - 9 v + 4) -
  Nc^4 (v^4 + 2 v^3 + 5 v^2 + 4) - (v - 2) (v - 1)^2 v;
f = Nc^2 (12 Pi^2 v^4 - 103 v^4 - 18 Pi^2 v^3 + 179 v^3 +
    15 Pi^2 v^2 - 188 v^2 + 9 Pi^2 v + 179 v - 6 Pi^2 - 103) +
  Nc^4 (9 Pi^2 v^4 + 85 v^4 + 18 Pi^2 v^3 + 9 v^3 +
    9 Pi^2 v^2 + 188 v^2 + 9 v + 85) +
  3 (v - 1)^2 (5 Pi^2 v^2 + 6 v^2 - 6 Pi^2 v + 2 Pi^2 + 6) -
  20 tr Nc^3 (v^2 + 1)^2 + 20 tr Nc (v - 1)^2 (v^2 + 1);
g = -(v - 1)^2 (3 v^2 - 4 v + 3) - Nc^4 (v + 1)^2 (v^2 + 1) + Nc^2 (v - 1)^2 v;
h = Nc^4 (v^2 + 10 v + 1) - Nc^2 (v^2 + 10 v + 1) + 4 (v - 1)^2;
publishedD1 = 2 (2 Nc^2 - 1)/Nc;
publishedD0 = -(6 a Log[v] - 12 b Log[1 - v] +
  6 (Nc^2 - 1) p LFF + 6 (3 Nc^2 - 1) p LF +
  (2 Nc^2 - 4 tr Nc + 9) p)/(6 Nc p);
publishedDelta = -(18 c Log[v]^2 +
  3 LF (v^2 + 1) p (24 Nc^2 Log[v] - 24 Nc^2 Log[1 - v] + 31 Nc^2 - 8 tr Nc - 9) +
  9 (Nc^2 - 1) LFF (v^2 + 1) p (4 Log[v] + 3) + 3 d Log[v] -
  36 Log[1 - v] e Log[v] - 2 f + 18 Log[1 - v]^2 g -
  18 Log[1 - v] v h - 12 Nc (11 Nc - 4 tr) LR (v^2 + 1) p)/
  (36 Nc (v^2 + 1) p);
published = {publishedD0, publishedD1, publishedDelta};
residuals = canonical /@ (normalizedAversa - published);
labels = {"D0", "D1", "delta"};
Do[check["endpoint/Eq43/" <> labels[[slot]], zeroQ[residuals[[slot]]]], {slot, 3}];

(* Finite-only perturbations leave pole cancellation and scale derivatives
   unchanged. Each is tested against the actual reconstructed endpoint. *)
Do[perturbed = ReplacePart[normalizedAversa,
    slot -> normalizedAversa[[slot]] + v (1 - v)];
  check["control/finite-only/" <> labels[[slot]],
    !zeroQ[perturbed[[slot]] - published[[slot]]]], {slot, 3}];
check["control/omit-scheme-conversion", !zeroQ[normalizedMSbar[[2]] - publishedD1]];
check["input-unchanged", inputHash === IntegerString[FileHash[coefficientFile, "SHA256"], 16, 64]];
If[!DirectoryQ[outputRoot], CreateDirectory[outputRoot, CreateIntermediateDirectories -> True]];
Put[<|"Source" -> "Aversa et al., NPB327 (1989)105, Eqs22,23,41,43; CQ=0", "Basis" -> labels,
  "MSbar" -> normalizedMSbar, "QuarkSchemeShift" -> quarkSchemeShift,
  "Aversa" -> normalizedAversa, "Published" -> published, "Residuals" -> residuals,
  "Checks" -> checks|>, FileNameJoin[{outputRoot, "endpoints.wl"}]];
Export[FileNameJoin[{outputRoot, "summary.json"}], <|"passed" -> (failures === {}),
  "checkCount" -> Length[checks], "failures" -> failures, "exactEndpointIdentities" -> Count[residuals, 0],
  "genericVariables" -> {"v", "Nc", "Nf", "LR", "LF", "LFF"},
  "coefficientFileSHA256" -> inputHash,
  "checkerSHA256" -> IntegerString[FileHash[$InputFileName, "SHA256"], 16, 64]|>, "RawJSON"];
Print["AVERSA ENDPOINTS: ", Length[checks], " checks, ", Length[failures], " failures; residuals=", residuals];
Quit[If[failures === {}, 0, 1]];
