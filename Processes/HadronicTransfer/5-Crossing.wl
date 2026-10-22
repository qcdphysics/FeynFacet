(* ::Title:: *)
(*Transverse spin transfer in proton-proton collisions, step 5: the generalized crossing from p^ p^ -> c + X*)


(* ::Text:: *)
(*Generalized crossing symmetry relates the transverse spin transfer cross sections to the double-transverse ones of Processes/HadronicTT: knowing a^ b^ -> c + X, the process a^ cbar -> bbar^ + X follows by s <-> t, v -> -v/(1-v), w unchanged, including the phase-space integration.  The implementation must account for log(v) -> log(-v/(1-v)), the pi^2 terms, the spin averaging factors and the interchange of the incoming-parton and fragmentation counterterms.  This step checks the crossing against the results of Processes/HadronicTT at LO (exactly) and at NLO (the complete finite coefficients), for every class of this process.  Input: Results/<route>/class-<n>/Algebra.wl of both processes (and Transfer.wl of this process for the NLO comparison); output: Results/<route>/crossing/class-<n>.wl per pair, and the pseudo-classes Results/<route>/class-<n>x/ described below.*)
(**)
(*The crossing.  The TT process a(pa) b(pb) -> c(p) + d(k1) [+ e(k2)] is crossed by b <-> cbar: pb -> -p and p -> -pb (the momenta of the transfer process), with the unobserved momenta k1, k2 unchanged.  The invariants: s <-> t, u unchanged, A = (pa-k1)^2 unchanged, B = (pb-k1)^2 <-> C = (p+k1)^2, spe3 (the hatted k1) unchanged; v = 1 + t/s -> -v/(1-v), w = -u/(s+t) unchanged.  The pairs (transfer class <- TT class; the unobserved parton k1 of the two class entries must coincide, which is checked on the species lists):*)
(**)
(*   3 <- 5:    q qb -> qp (k1 = qpb)  -->  q qpb -> q + X;    11 <- 11:  q qb -> q (k1 = qb)  -->  q qb -> q + X;    17 <- 6:  q q -> q  -->  q qb -> qb (= qb q -> q by charge conjugation, class 17);*)
(*   13 <- 12:  q qb -> g  -->  q g -> q + X;    10 <- 7:  q q -> g (NLO only)  -->  q g -> qb + X;*)
(*   1 <- 5 and 6 <- 11 observe the antiquark of the TT pair: the TT process with the observed antiquark, q qb -> qpb (q qb -> qb), equals by charge conjugation the class's own square with a and b exchanged, pa <-> pb (t <-> u, A <-> B; checked in step 2 of HadronicTT on the Born table), so these two pairs use the composed map pa -> -p, pb -> pa, p -> -pb: s -> t, t -> u, u -> s, A -> C, B -> A, C -> B.*)
(**)
(*The constant factor.  Per pair the transfer square is a constant times the crossed TT square, component by component:*)
(**)
(*   XX_transfer = - sigma g x XX_TT(crossed),   NN_transfer = + sigma g x NN_TT(crossed),   sigma = +1, g = 1 (b a quark or antiquark),   sigma = -1, g = 2 CA/((D-2)(CA^2-1)) = 1/((D-2) CF) (b a gluon).*)
(**)
(*Derivation.  (i) Spinor structures: the TT density of b, pb-slash gamma5 sb-slash/2, becomes -p-slash gamma5 sb-slash/2 at pb -> -p, while the transfer projection of the observed quark is p-slash gamma5 sc-slash (no 1/2); the TT spin sum of an observed quark, p-slash, becomes -pb-slash at p -> -pb, while the transfer's unpolarized quark b carries pb-slash/2: the product (-1/2)(-2) = 1 contains the two fermion crossing signs.  When the TT observed parton is a gluon (the pairs 13 <- 12 and 10 <- 7, whose transfer parton b is that gluon), its polarization sum is even in its momentum and the second sign is absent: sigma = -1 for these pairs, +1 otherwise.  (ii) Colour: the TT average 1/Nc^2 against the transfer's 1/Nc x 1/N_b with N_b = Nc (quark) or Nc^2 - 1 (gluon).  (iii) A TT gluon c is summed over its physical polarizations with the reference vector pa, and the transfer's incoming gluon b is averaged with the same sum, reference pa, and 1/(D-2); the sums coincide (the reference momentum enters as a ratio).  Together: sigma = 1, g = 1 for a quark b and sigma = -1, g = 2 Nc/((D-2)(Nc^2-1)) for a gluon b.  (iv) The axes.  The TT in-plane axis aX = p + (u/s) pa + (t/s) pb, crossed, is -pb + (u/t) pa - (s/t) p (transfer momenta), orthogonal to pa and p.  As the spin vector of a it equals (-s/t) eX + (2u/t) pa, where eX = (u/s) pa + (t/s) pb + p is the transfer's in-plane axis of a and the multiple of pa is annihilated by pa-slash gamma5; as the spin vector of c it equals -((t+u)/t) eC - 2 s u/(t (t+u)) p with the transfer's eC = -u/(t+u) pa + t/(t+u) pb - s (u-t)/(t+u)^2 p, the multiple of p annihilated by p-slash gamma5.  With the normalizations (TT divides by -aX^2 -> s u/t; the transfer by |eX| |eC| = -t u/(t+u)) the product of the two coefficients gives exactly -1: the crossed in-plane axis of b points along -eC (the azimuth map phi_b -> pi - phi_c).  The normal n ~ eps^{mu pa pb p} is reflected under the crossing, but both TT densities carry it: NN keeps its sign.  (v) Projections: TT's cos(2 Phi) coefficient TT = (XX-NN)/2 therefore crosses into -(XX+NN)/(2 sigma g) of the transfer (the coefficient of cos(phi_a - phi_c)), and TT's evanescent Symmetric = (XX+NN)/2 into -(XX-NN)/(2 sigma g), the difference of the two transfer components (not into their sum).  All signs and factors are fixed by this derivation and verified exactly on the D-dimensional Born terms of the six pairs with a Born term; for the NLO-only pair 10 <- 7 (a gluon b, like 13 <- 12 whose Born term verifies sigma and g including the D dependence) the same formula is used and is tested by the pole cancellation of Part B: its real-emission poles must cancel against the transfer counterterm P_{qb <- g} x Born(q qb -> qb)_transfer, which fixes the normalization including the gluon average.*)
(**)
(*The one-loop families.  The scalar integrals of the TT virtual interference are defined by their propagators; under the crossing VirtualST (s,t box) maps onto itself with the indices permuted (s <-> t), and VirtualSU <-> VirtualTU; for the composed map all three permute cyclically.  The map is found on the propagator sets by a shift and a reflection of the loop momentum (the integrals are invariant under both, the hatted loop momentum too) and verified to be a permutation of the family's propagators.  The integrals are then those of the transfer process with the transfer's Feynman prescription: the crossing is done before any integral is evaluated, and the real parts are taken in the transfer's physical region by step 3.*)
(**)
(*Part A (exact, algebra level) compares, for every pair, the crossed TT objects of step 2 with those of this process: the Born squares and the complete Born table (rational-function identities in s, t, D, CA), the virtual XX and NN after the reduction to master integrals (term by term the scalar-integral coefficients agree only when both files route every graph onto the same family members, which is recorded; in general the two sides differ by integral identities, such as the same bubble written in two families), XX and NN modulo the integration identities that relate the two treatments of the normal component of the loop momentum (HadronicTT uses the exact identity (n.loop)^2 = loop_par^2 - loop^2 + loopHat2, this process the angular average -loop_perp^2/(D-3); both sides are brought to a common form by replacing every hatted power by its angular average, loopHat2^h -> ((D-4)/2)_h/((D-3)/2)_h (loop_perp^2)^h, and reduced with Kira in D dimensions: the coefficients of the master integrals must agree exactly), and the real-emission squares (rational functions with tens of thousands of terms, compared exactly at random rational points of all invariants).*)
(**)
(*Part B (the complete NLO check) writes the crossed TT algebra as a pseudo-class "<n>x" of this process (the transfer class's entry of classes.wl as its specification, the crossed TT Born table as its counterterm table, which Part A has shown to equal the transfer's) and runs steps 3 and 4 of this process on it, exactly as on a class of its own.  This implements generalized crossing symmetry at the level of phase-space integrals: the crossed integrands are integrated by step 3 in the transfer region (the phase-space masters at s4 = s v (1-w) of the transfer process, the loop masters continued in the transfer's s channel, with the real part taken there), so no finite TT expression is continued and no log(v) -> log(-v/(1-v)) or pi^2 term has to be tracked by hand; the counterterms are those of step 4 of this process (transversity kernels on leg a and on the fragmenting quark, unpolarized kernels on leg b with the lower transfer Born terms): this is the interchange of the incoming-parton and fragmentation counterterms, and the cancellation of the poles of the crossed integrals against them tests the crossing factors once more, including the pair without a Born term.  The finite coefficients {regular, D0, D1, delta} of both components of the pseudo-class must then equal those of the directly computed class, exactly at rational points of s, v, w, Nf with symbolic scales; the real-emission and virtual finite parts are compared separately as well.*)


(* ::Section:: *)
(*Load*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "..", "..", "FeynFacet", "FeynFacet.wl"}]]


(* ::Text:: *)
(*Choices: the route of the step-2 files; the transfer classes to check ("all", or a comma-separated list such as "3,17"; a class set through class = "<n>", as RunSteps[..., "class=<n>"] does, selects that one pair); the stage ("A": Part A and the pseudo-class algebra only; "all": also steps 3 and 4 of the pseudo-class, run here in fresh kernels when their results are missing, and the NLO comparison); and whether pseudo-class results of an earlier run of this file may be reused ("no" reruns steps 3 and 4).  Every pair reads the TT file Processes/HadronicTT/Results/<route>/class-<m>/Algebra.wl and this process's Algebra.wl (and, for the NLO comparison, Transfer.wl) of class <n>.*)


(* ::Input:: *)
route = "FeynCalc";
class = "";                                  (* set by RunSteps[..., "class=<n>"]: the pair of that transfer class only *)
pairs = "all";
stage = "all";
reuse = "yes";
Get[FileNameJoin[{NotebookDirectory[], "classes.wl"}]];
ttDirectory = FileNameJoin[{NotebookDirectory[], "..", "HadronicTT"}];
crossingTable = <|
  "1" -> <|"TT" -> "5", "Relabel" -> True|>, "3" -> <|"TT" -> "5", "Relabel" -> False|>, "6" -> <|"TT" -> "11", "Relabel" -> True|>,
  "10" -> <|"TT" -> "7", "Relabel" -> False|>, "11" -> <|"TT" -> "11", "Relabel" -> False|>, "13" -> <|"TT" -> "12", "Relabel" -> False|>,
  "17" -> <|"TT" -> "6", "Relabel" -> False|>|>;
pairList = Which[class =!= "", {class}, pairs === "all", Keys[crossingTable], True, StringTrim /@ StringSplit[pairs, ","]];
Require[SubsetQ[Keys[crossingTable], pairList], "Unknown transfer class in the pairs", pairList];
crossingDirectory = ResultsDirectory[route <> "/crossing"];
ttAlgebraFile[m_] := FileNameJoin[{ttDirectory, "Results", route, "class-" <> m, "Algebra.wl"}];
transferFile[n_, name_] := FileNameJoin[{ResultsDirectory[route <> "/class-" <> n], name <> ".wl"}];
pairList


(* ::Section:: *)
(*The crossing map*)


(* ::Text:: *)
(*The momentum map of a pair (direct: pb -> -p, p -> -pb; composed with a <-> b for the pairs that observe the antiquark: pa -> -p, pb -> pa, p -> -pb) is turned into rules for the invariants by evaluating the squares (pa+pb)^2, (pa-p)^2, (pb-p)^2, (pa-k1)^2, (pb-k1)^2, (p+k1)^2 with the transfer's scalar products (the D-dimensional k1^2 = 0; spe3 is k1's own and unchanged).  The constant factor g and the sign sigma of the pair are read from the species of the transfer's parton b, and the signs of the components are those of the introduction.  The TT files hold the projections TT and Symmetric of the virtual and real parts; XX = TT + Symmetric and NN = Symmetric - TT are rebuilt (the Born holds all four).  The SU(N) relation CF = (CA^2-1)/(2 CA) is used in every exact comparison, since the colour factors of the two processes are written through CA and CF differently.*)


(* ::Input:: *)
momentumMap[relabel_] := If[relabel, {pa -> -p, pb -> pa, p -> -pb}, {pa -> pa, pb -> -p, p -> -pb}];
products = <|{pa, pb} -> s/2, {p, pa} -> -t/2, {p, pb} -> -u/2, {k1, pa} -> -A/2, {k1, pb} -> -B/2, {k1, p} -> C/2,
  {loop, loop} -> L2, {loop, pa} -> Lpa, {loop, pb} -> Lpb, {loop, p} -> Lp|>;
scalarProduct[x_, y_] := Lookup[products, Key[Sort[{x, y}]], 0];              (* pa^2 = pb^2 = p^2 = k1^2 = 0 *)
squareOf[q_] := Module[{e = Expand[q], vectors = {pa, pb, p, k1, loop}},
  Expand[Sum[Coefficient[e, x] Coefficient[e, y] scalarProduct[x, y], {x, vectors}, {y, vectors}]]];
crossingRules[map_] := {s -> squareOf[(pa + pb) /. map], t -> squareOf[(pa - p) /. map], u -> squareOf[(pb - p) /. map],
  A -> squareOf[(pa - k1) /. map], B -> squareOf[(pb - k1) /. map], C -> squareOf[(p + k1) /. map]};
colourRelation = CF -> (CA^2 - 1)/(2 CA);
signs = <|"XX" -> -1, "NN" -> 1|>;
gluonFactor = 2 CA/((D - 2) (CA^2 - 1));
crossingSign[b_] := If[b === "g", -1, 1];                     (* sigma: one fermion crossing sign is absent when the TT observed parton (the transfer's b) is a gluon *)
crossingFactor[b_] := If[b === "g", gluonFactor, 1];          (* g *)
conjugate = <|"q" -> "qb", "qb" -> "q", "qp" -> "qpb", "qpb" -> "qp", "g" -> "g"|>;
partonsOf[spec_] := (spec["Born"] /. None -> spec["Real"][[1, 1]])[[1 ;; 4]];        (* {a, b, c, d}: d the unobserved parton k1 *)
ttComponents[object_] := <|"XX" -> object["TT"] + object["Symmetric"], "NN" -> object["Symmetric"] - object["TT"]|>;
Column[Table[n <> " <- " <> crossingTable[n]["TT"] <> ":  " <> ToString[crossingRules[momentumMap[crossingTable[n]["Relabel"]]], InputForm], {n, Keys[crossingTable]}]]


(* ::Text:: *)
(*One pair: the two algebra files, the species check (the crossed TT species a, cbar, bbar, d must be the transfer's a, b, c, d; for the composed map the TT species are first exchanged a <-> b and conjugated, which is the charge-conjugation argument above), the invariant rules and the factor.*)


(* ::Input:: *)
loadPair[n_] := Module[{entry = crossingTable[n], tt, transfer, spec, map, ttPartons, crossedPartons, b, g},
  tt = Get[ttAlgebraFile[entry["TT"]]];
  transfer = Get[transferFile[n, "Algebra"]];
  spec = HadronicTransferClasses[ToExpression[n]];
  Require[spec[{"Name", "Born", "Virtual", "Real"}] === transfer["Spec"][{"Name", "Born", "Virtual", "Real"}], "The class entry differs from the one used by step 2", n];
  map = momentumMap[entry["Relabel"]];
  ttPartons = partonsOf[tt["Spec"]];
  If[entry["Relabel"], ttPartons = conjugate /@ ttPartons[[{2, 1, 3, 4}]]];                    (* a <-> b and charge conjugation: the same TT square *)
  crossedPartons = {ttPartons[[1]], conjugate[ttPartons[[3]]], conjugate[ttPartons[[2]]], ttPartons[[4]]};
  (* class 17 is written as the charge conjugate (qb q -> q) of the crossed process (q qb -> qb): the same square, checked in step 2 on the Born table *)
  Require[crossedPartons === partonsOf[spec] || (conjugate /@ crossedPartons) === partonsOf[spec], "The crossed TT species are not those of the transfer class (up to charge conjugation)", {n, crossedPartons, partonsOf[spec]}];
  b = partonsOf[spec][[2]];
  <|"Transfer" -> n, "TT" -> entry["TT"], "Relabel" -> entry["Relabel"], "Map" -> map, "Rules" -> crossingRules[map], "Factor" -> crossingFactor[b], "Sign" -> crossingSign[b], "Spec" -> spec,
    "TTAlgebra" -> tt, "TransferAlgebra" -> transfer, "Name" -> spec["Name"] <> "  <-  " <> tt["Spec"]["Name"] <> If[entry["Relabel"], " (observed antiquark, a <-> b)", ""]|>];
crossObject[pair_, expr_, component_] := signs[component] pair["Sign"] pair["Factor"] (expr /. pair["Rules"]);
loaded = Association @ Table[n -> loadPair[n], {n, pairList}];
Column[Table[loaded[n]["Name"] <> ":  sigma = " <> ToString[loaded[n]["Sign"]] <> ", g = " <> ToString[loaded[n]["Factor"], InputForm], {n, pairList}]]


(* ::Section:: *)
(*Part A: the Born terms and the Born table*)


(* ::Text:: *)
(*The crossed TT Born square of the pair must equal the transfer's Born square in D dimensions, for XX and NN (classes without a Born term: both vanish), and the projections cross as stated: -sigma g x TT -> (XX+NN)/2, -sigma g x Symmetric -> (XX-NN)/2.  The complete transfer Born table follows from the TT table by the direct crossing of every channel, a b -> c (k1 = d) into a cbar -> bbar (k1 = d), with sigma and g of the species of cbar: the seven transfer channels of classes.wl are all obtained (qb q -> q from the TT channel qb qb -> qb).  At LO the relation can be written in terms of the hard cross sections b(v) = |M|^2 (4 Pi alphaS)^2/(16 Pi^2 s v) of step 4: since s' v' = s v under the crossing (s' = t = -s(1-v), v' = -v/(1-v)) and b is proportional to 1/s,*)
(**)
(*   b_transfer,XX(v) = + sigma g b_TT,XX(-v/(1-v))/(1-v),   b_transfer,NN(v) = - sigma g b_TT,NN(-v/(1-v))/(1-v),   and at D = 4 both equal sigma g b_TT(-v/(1-v))/(1-v) with b_TT the cos(2 Phi) coefficient of HadronicTT (b_TT,XX = -b_TT,NN = b_TT there),*)
(**)
(*an identity of rational functions of v, checked here with the D-dimensional Born squares of the two files (for the pairs 1 <- 5 and 6 <- 11, b_TT is that of the TT process with the observed antiquark, the class's own with t <-> u).*)


(* ::Input:: *)
exactZero[expr_] := Factor[expr /. u -> -s - t /. colourRelation] === 0;
hard[square_] := square (4 Pi alphaS)^2/(16 Pi^2 s v);
bornCheck[pair_] := Module[{tt = pair["TTAlgebra"], transfer = pair["TransferAlgebra"], differences, projections, slide, table, tableDifferences},
  differences = Association @ Table[comp -> exactZero[crossObject[pair, tt["Born"][comp], comp] - transfer["Born"][comp]], {comp, {"XX", "NN"}}];
  Require[And @@ Values[differences], "The crossed TT Born square differs from the transfer's", {pair["Name"], differences}];
  projections = {exactZero[-pair["Sign"] pair["Factor"] (tt["Born"]["TT"] /. pair["Rules"]) - (transfer["Born"]["XX"] + transfer["Born"]["NN"])/2],
    exactZero[-pair["Sign"] pair["Factor"] (tt["Born"]["Symmetric"] /. pair["Rules"]) - (transfer["Born"]["XX"] - transfer["Born"]["NN"])/2]};
  Require[projections === {True, True}, "The TT projections do not cross into -(XX+NN)/2 and -(XX-NN)/2", pair["Name"]];
  (* the slide at LO, with the hard cross sections as functions of v at fixed s *)
  (* for the pairs that observe the antiquark the TT square is the class's own with a <-> b, i.e. t <-> u (t -> -s - t with u eliminated) *)
  slide = Table[With[{ttBorn = If[pair["Relabel"], tt["Born"][comp] /. u -> -s - t /. t -> -s - t, tt["Born"][comp] /. u -> -s - t]},
    Factor[(hard[transfer["Born"][comp]] /. t -> -s (1 - v) /. colourRelation) -
      (-signs[comp]) pair["Sign"] pair["Factor"] ((hard[ttBorn] /. t -> -s (1 - v)) /. v -> -v/(1 - v))/(1 - v) /. colourRelation] === 0], {comp, {"XX", "NN"}}];
  Require[slide === {True, True}, "The LO crossing relation of the hard cross sections fails", pair["Name"]];
  (* the Born table: every transfer channel from a TT channel by the direct crossing *)
  table = Association @ Table[name -> Module[{ch = StringSplit[StringReplace[name, "->" -> ""]], ttName, ttChannel, g},
      ttChannel = {ch[[1]], conjugate[ch[[3]]], conjugate[ch[[2]]]};
      ttName = ttChannel[[1]] <> " " <> ttChannel[[2]] <> " -> " <> ttChannel[[3]];
      Require[KeyExistsQ[tt["BornTable"]["XX"], ttName], "The TT Born table lacks the channel needed for " <> name, ttName];
      g = crossingSign[ch[[2]]] crossingFactor[ch[[2]]];
      <|"TT channel" -> ttName, "Factor" -> g, "Crossed" -> Association @ Table[comp -> signs[comp] g (tt["BornTable"][comp][ttName] /. crossingRules[momentumMap[False]]), {comp, {"XX", "NN"}}]|>],
    {name, Keys[transfer["BornTable"]["XX"]]}];
  tableDifferences = Association @ Table[{name, comp} -> exactZero[table[name]["Crossed"][comp] - transfer["BornTable"][comp][name]], {name, Keys[table]}, {comp, {"XX", "NN"}}];
  Require[And @@ Values[tableDifferences], "A crossed TT Born table entry differs from the transfer's", Select[tableDifferences, # === False &]];
  <|"Born" -> differences, "Projections" -> projections, "SlideLO" -> slide, "BornTable" -> Map[#["TT channel"] &, table], "CrossedTable" -> Map[#["Crossed"] &, table]|>];
bornResults = Association @ Table[n -> bornCheck[loaded[n]], {n, pairList}];
Column[Table[loaded[n]["Name"] <> ": Born XX, NN and the table of " <> ToString[Length[bornResults[n]["BornTable"]]] <> " channels cross exactly", {n, pairList}]]


(* ::Section:: *)
(*Part A: the virtual interference*)


(* ::Text:: *)
(*The family map: for each TT family, the crossed propagator momenta are matched, up to a shift and a reflection of the loop momentum, with a permutation of the propagators of one family of this process (the families are the same in both files, which is checked); the indices of KVirtual..[n1,n2,n3,n4] are permuted accordingly, the hat powers are unchanged.  The XX component is then compared coefficient by coefficient (the coefficient of every distinct integral x hat power), after dropping the scaleless integrals (those whose propagators with positive index have only vanishing differences of momenta); this term-by-term identity holds when both files route the graphs alike and is recorded, not required: the required comparison is the one after the reduction (next cell).*)


(* ::Input:: *)
familyHeads = <|"VirtualST" -> KVirtualST, "VirtualSU" -> KVirtualSU, "VirtualTU" -> KVirtualTU|>;
headNames = Association[Reverse /@ Normal[familyHeads]];
integralHeads = Alternatives @@ Values[familyHeads];
shifts = Flatten[Table[c1 pa + c2 pb + c3 p, {c1, -1, 1}, {c2, -1, 1}, {c3, -1, 1}]];
sameUpToSign[q1_, q2_] := Expand[q1 - q2] === 0 || Expand[q1 + q2] === 0;
(* {target family, sigma, delta, permutation}: propagator i of the crossed family is propagator permutation[[i]] of the target *)
familyMap[families_, map_] := Association @ Table[name -> Module[{props = families[name] /. map, found},
    found = Catch[Do[Module[{shifted = props /. loop -> sigma loop + delta, perm},
        perm = Table[FirstPosition[families[target], f_ /; sameUpToSign[f, q], {0}][[1]], {q, shifted}];
        If[FreeQ[perm, 0] && Sort[perm] === Range[4], Throw[{target, sigma, delta, perm}]]],
      {target, Keys[families]}, {sigma, {1, -1}}, {delta, shifts}]; None];
    Require[found =!= None, "No family of this process matches the crossed propagators of " <> name, props]; found],
  {name, Keys[families]}];
integralRules[fm_] := Table[With[{h = familyHeads[name], g = familyHeads[fm[name][[1]]], perm = fm[name][[4]]},
  h[n__] :> g @@ Permute[{n}, InversePermutation[perm]]], {name, Keys[fm]}];
(* the terms of a virtual expression grouped by integral and hat power: {integral, HatPower[h]} -> coefficient *)
monomials[expr_] := If[expr === 0, <||>, Module[{terms = With[{e = Expand[expr]}, Which[e === 0, {}, Head[e] === Plus, List @@ e, True, {e}]]},   (* a single product is one term *)
  Association @ KeyValueMap[#1 -> Total[#2] &, GroupBy[terms,
    Function[term, {FirstCase[term, integralHeads[__], None, {0, Infinity}], FirstCase[term, _HatPower, None, {0, Infinity}]}],
    Function[term, term/(FirstCase[term, integralHeads[__], 1, {0, Infinity}] FirstCase[term, _HatPower, 1, {0, Infinity}])]]]]];
scalelessQ[families_, (head_)[n__]] := Module[{q = families[headNames[head]], present = Flatten[Position[{n}, _?Positive]]},
  And @@ Flatten[Table[squareOf[q[[i]] - q[[j]]] === 0, {i, present}, {j, present}]]];
compareMonomials[families_, a_, b_] := Module[{ma = monomials[a], mb = monomials[b], keys},
  keys = Select[Union[Keys[ma], Keys[mb]], !scalelessQ[families, #[[1]]] &];
  Select[Association @ Table[k -> Together[(Lookup[ma, Key[k], 0] - Lookup[mb, Key[k], 0]) /. colourRelation], {k, keys}], # =!= 0 &]];
crossVirtual[pair_, component_] := With[{tt = pair["TTAlgebra"]},
  crossObject[pair, ttComponents[tt["Virtual"]][component] /. integralRules[familyMap[tt["Families"], pair["Map"]]], component] /. u -> -s - t];   (* 2 -> 2: u eliminated again (the composed map sends t to u) *)
virtualXXCheck[pair_] := Module[{tt = pair["TTAlgebra"], transfer = pair["TransferAlgebra"], fm, bad},
  Require[tt["Families"] === transfer["Families"], "The one-loop families of the two processes differ"];
  fm = familyMap[tt["Families"], pair["Map"]];
  bad = compareMonomials[tt["Families"], crossVirtual[pair, "XX"], transfer["Virtual"]["XX"]];
  (* identical term by term when both files route every graph onto the same family members; otherwise equal only up to integral identities (next cell) *)
  <|"FamilyMap" -> fm, "Integrals" -> Length[monomials[transfer["Virtual"]["XX"]]], "TermByTerm" -> (bad === <||>), "DifferingTerms" -> Length[bad]|>];
virtualXXResults = Association @ Table[n -> If[loaded[n]["Spec"]["Virtual"], virtualXXCheck[loaded[n]], <|"FamilyMap" -> None, "TermByTerm" -> "no virtual correction"|>], {n, pairList}];
Column[Table[loaded[n]["Name"] <> ": virtual XX identical term by term: " <> ToString[virtualXXResults[n]["TermByTerm"]] <> "; family map " <>
  ToString[If[virtualXXResults[n]["FamilyMap"] === None, None, Map[#[[1]] &, virtualXXResults[n]["FamilyMap"]]], InputForm], {n, pairList}]]


(* ::Text:: *)
(*The NN component contains the normal component of the loop momentum, which the two step-2 files reduce differently, and cannot be compared term by term.  Both sides are brought to a common form: every hatted power is replaced by its angular average over the (D-3)-dimensional space orthogonal to the scattering plane, loopHat2^h -> ((D-4)/2)_h/((D-3)/2)_h (loop_perp^2)^h with loop_perp^2 = loop^2 - d.G^-1.d written through the inverse propagators of the family (so that a term becomes a sum of integrals with numerators); the union of the integrals of both sides is reduced by Kira (the families of kira/virtual, copied to Results/<route>/crossing/kira/class-<n>/virtual, with the seed ranges enlarged if the numerators need it) in D dimensions, the masters are identified by their values (bubble of one channel, box of a pair of channels, as in step 3, since the same bubble can be a master of two families), and their coefficients must agree exactly as rational functions of s, t, D, CA, Nf.  The XX component is passed through the same reduction: this is the required check for XX as well.*)


(* ::Input:: *)
loopSymbols = {L2, Lpa, Lpb, Lp};
gram = {{0, s/2, -t/2}, {s/2, 0, -u/2}, {-t/2, -u/2, 0}} /. u -> -s - t;
inverseProps[families_, name_] := inverseProps[families, name] = Module[{sol},
  sol = Solve[Table[Dj[j] == (squareOf[families[name][[j]]] /. u -> -s - t), {j, 4}], loopSymbols];
  Require[Length[sol] == 1, "Cannot express the loop scalar products through the propagators of " <> name]; First[sol]];
loopPerp2[families_, name_] := loopPerp2[families, name] = Factor[(L2 - {Lpa, Lpb, Lp} . Inverse[gram] . {Lpa, Lpb, Lp}) /. inverseProps[families, name]];
perpMoment[h_] := Pochhammer[(D - 4)/2, h]/Pochhammer[(D - 3)/2, h];
(* every integral x HatPower[h] -> moment x sum of integrals with the numerator (loop_perp^2)^h *)
dehat[families_, expr_] := Total[KeyValueMap[Function[{key, coefficient}, With[{k = key[[1]], h = key[[2]][[1]]},
  If[h == 0, coefficient k,
    coefficient perpMoment[h] Total[(#[[2]] (Head[k] @@ (List @@ k - #[[1]])) &) /@ CoefficientRules[Expand[loopPerp2[families, headNames[Head[k]]]^h], Array[Dj, 4]]]]]],
  monomials[expr]]];
virtualHeads = {KVirtualST -> "virtual_ST", KVirtualSU -> "virtual_SU", KVirtualTU -> "virtual_TU"};
reduceBoth[pair_, components_] := Module[{families = pair["TTAlgebra"]["Families"], sides, targets, rMax, sMax, kiraDir, jobs, rules, masters, keys},
  sides = Association @ Table[comp -> {dehat[families, crossVirtual[pair, comp]], dehat[families, pair["TransferAlgebra"]["Virtual"][comp]]}, {comp, components}];
  targets = Union[Cases[sides, integralHeads[__], Infinity]];
  rMax = Max[Total[Select[List @@ #, Positive]] & /@ targets]; sMax = Max[-Total[Select[List @@ #, Negative]] & /@ targets];
  kiraDir = PrepareKira[FileNameJoin[{NotebookDirectory[], "kira", "virtual"}], FileNameJoin[{crossingDirectory, "kira", "class-" <> pair["Transfer"], "virtual"}]];
  jobs = FileNameJoin[{kiraDir, "jobs.yaml"}];
  Export[jobs, StringReplace[Import[jobs, "Text"], "r: 7, s: 4" -> "r: " <> ToString[Max[7, rMax]] <> ", s: " <> ToString[Max[4, sMax]]], "Text"];
  WriteKiraTargets[kiraDir, targets, virtualHeads];
  RunKira[kiraDir];
  rules = RestoreScale[Join @@ Table[KiraRules[kiraDir, name, virtualHeads], {name, Select[Values[virtualHeads], FileExistsQ[FileNameJoin[{kiraDir, "results", #, "kira_targets.m"}]] &]}], s, {t -> t/s}, Keys[virtualHeads]];
  sides = Map[Collect[Expand[# /. Dispatch[rules] /. d -> D], integralHeads[__], Together] &, sides, {2}];
  masters = Union[Cases[sides, integralHeads[__], Infinity]];
  (* a target that is itself a master has no rule: what remains must be bubbles and boxes (indices 0 and 1 only), as in step 3 *)
  Require[And @@ (SubsetQ[{0, 1}, List @@ #] & /@ masters), "Some integrals were not reduced by Kira", Select[masters, !SubsetQ[{0, 1}, List @@ #] &]];
  (* masters of different families with the same value (a bubble of one channel, the box of one pair of channels) are identified by their value, as in step 3 *)
  sides = Map[Collect[# /. Dispatch[Table[m -> masterValue[families, m], {m, masters}]], _bubble | _box | _triangle, Together] &, sides, {2}];
  keys = Union[Cases[sides, _bubble | _box | _triangle, Infinity]];
  <|"Targets" -> Length[targets], "Ranges" -> {rMax, sMax}, "Masters" -> masters, "MasterValues" -> keys,
    "Differences" -> Association @ Table[comp -> Select[Association @ Table[m -> Together[(Coefficient[sides[comp][[1]], m] - Coefficient[sides[comp][[2]], m]) /. colourRelation], {m, keys}], # =!= 0 &], {comp, components}]|>];
(* the value of a master (indices 0, 1): bubble[P^2] of the momentum through its two propagators, box[{S, T}] of its two channels (symmetric), with u = -s - t *)
masterValue[families_, (head_)[n__]] := Module[{q = families[headNames[head]], present = Flatten[Position[{n}, 1]], inv},
  inv[x_] := Factor[squareOf[x] /. u -> -s - t];
  Switch[Length[present],
    2, bubble[inv[q[[present[[1]]]] - q[[present[[2]]]]]],
    4, box[Sort[{inv[q[[3]] - q[[1]]], inv[q[[4]] - q[[2]]]}]],
    _, triangle[Sort[Table[inv[q[[present[[i]]]] - q[[present[[Mod[i, Length[present]] + 1]]]]], {i, Length[present]}]]]]];
virtualNNResults = Association @ Table[n -> If[loaded[n]["Spec"]["Virtual"], reduceBoth[loaded[n], {"NN", "XX"}], "no virtual correction in this class"], {n, pairList}];
Do[Require[!AssociationQ[virtualNNResults[n]] || And @@ (# === <||> & /@ Values[virtualNNResults[n]["Differences"]]),
  "The crossed TT virtual NN (or XX) differs from the transfer's after the reduction to master integrals", {loaded[n]["Name"], virtualNNResults[n]["Differences"]}], {n, pairList}];
Column[Table[loaded[n]["Name"] <> ": " <> If[AssociationQ[virtualNNResults[n]],
  "virtual NN and XX agree exactly on " <> ToString[Length[virtualNNResults[n]["MasterValues"]]] <> " distinct master integrals " <> ToString[virtualNNResults[n]["MasterValues"], InputForm] <> " (" <> ToString[virtualNNResults[n]["Targets"]] <> " targets)", virtualNNResults[n]], {n, pairList}]]


(* ::Section:: *)
(*Part A: the real emission*)


(* ::Text:: *)
(*The crossed TT real-emission square (XX and NN, and the two projections) against the transfer's, exactly at random rational points of s, t, u, A, B (C = -s4 - A - B), D, CA, spe3, Nf: the expressions have tens of thousands of terms and are not brought to a common denominator.  The reference vectors of the gluon polarization sums differ between the two files (the TT's k2 uses pb, i.e. -p after the crossing; the transfer's uses pb), which is immaterial for the complete squares (gauge invariance, checked in both step-2 files).*)


(* ::Input:: *)
realPoints = Table[With[{r = RandomInteger[{1, 97}, 8]/RandomInteger[{98, 197}, 8]},
  {s -> 1, t -> -r[[1]], u -> -r[[2]] (1 - r[[1]]), A -> -r[[3]], B -> -r[[4]], D -> 4 - r[[5]], CA -> 2 + 10 r[[6]], spe3 -> r[[7]], Nf -> 5}], {3}];
atPoint[expr_, point_] := Together[expr /. C -> -(s + t + u) - A - B /. colourRelation /. point];
realCheck[pair_] := Module[{tt = ttComponents[pair["TTAlgebra"]["Real"]], transfer = pair["TransferAlgebra"]["Real"], differences, projections},
  differences = Association @ Table[comp -> Table[atPoint[crossObject[pair, tt[comp], comp] - transfer[comp], point], {point, realPoints}], {comp, {"XX", "NN"}}];
  Require[And @@ (# === {0, 0, 0} & /@ Values[differences]), "The crossed TT real-emission square differs from the transfer's", {pair["Name"], differences}];
  projections = {Table[atPoint[-pair["Sign"] pair["Factor"] (pair["TTAlgebra"]["Real"]["TT"] /. pair["Rules"]) - (transfer["XX"] + transfer["NN"])/2, point], {point, realPoints}],
    Table[atPoint[-pair["Sign"] pair["Factor"] (pair["TTAlgebra"]["Real"]["Symmetric"] /. pair["Rules"]) - (transfer["XX"] - transfer["NN"])/2, point], {point, realPoints}]};
  Require[projections === {{0, 0, 0}, {0, 0, 0}}, "The TT projections of the real emission do not cross into -(XX+NN)/2 and -(XX-NN)/2", pair["Name"]];
  <|"Real" -> differences, "Projections" -> projections, "LeafCounts" -> LeafCount /@ transfer|>];
realResults = Association @ Table[n -> realCheck[loaded[n]], {n, pairList}];
Column[Table[loaded[n]["Name"] <> ": real emission XX and NN cross exactly at " <> ToString[Length[realPoints]] <> " rational points", {n, pairList}]]


(* ::Text:: *)
(*Part A is saved per pair (the map, the factor, the family map, the checks), before Part B.*)


(* ::Input:: *)
partA = Association @ Table[n -> <|"Pair" -> loaded[n]["Name"], "Transfer" -> n, "TT" -> loaded[n]["TT"], "Map" -> loaded[n]["Map"], "Rules" -> loaded[n]["Rules"], "Factor" -> loaded[n]["Factor"], "Sign" -> loaded[n]["Sign"],
  "Signs" -> signs, "Born" -> bornResults[n], "VirtualXX" -> virtualXXResults[n], "VirtualReduced" -> virtualNNResults[n], "Real" -> realResults[n], "Route" -> route,
  "Projections" -> "-sigma g x TT -> (XX+NN)/2, -sigma g x Symmetric -> (XX-NN)/2 of the transfer process"|>, {n, pairList}];
Do[SaveResult["class-" <> n, <|"PartA" -> partA[n], "PartB" -> "not run"|>, route <> "/crossing"], {n, pairList}];
"Part A saved for " <> StringRiffle[pairList, ", "]


(* ::Section:: *)
(*Part B: the pseudo-classes*)


(* ::Text:: *)
(*The crossed TT algebra is written as Results/<route>/class-<n>x/Algebra.wl in the format of step 2 of this process: Born, Virtual and Real with the components XX and NN (the crossed TT objects with their factors; the virtual with the families relabelled), the Born table crossed from the TT table (equal to the transfer's by Part A), the class entry of classes.wl as the specification.  Steps 3 and 4 of this process are then run on the pseudo-class in fresh kernels (RunSteps, with the override class=<n>x; skipped when their results exist and reuse = "yes"): step 3 integrates the crossed integrands in the transfer region, step 4 assembles them with the transfer counterterms and requires the cancellation of the poles.*)


(* ::Input:: *)
pseudoName[n_] := n <> "x";
writePseudoClass[pair_] := Module[{tt = pair["TTAlgebra"], n = pair["Transfer"], algebra},
  algebra = <|"Born" -> Association @ Table[comp -> Factor[crossObject[pair, tt["Born"][comp], comp] /. u -> -s - t], {comp, {"XX", "NN"}}],
    "Virtual" -> Association @ Table[comp -> Collect[crossVirtual[pair, comp], {integralHeads[__], _HatPower}, Factor], {comp, {"XX", "NN"}}],
    "Real" -> Association @ Table[comp -> crossObject[pair, ttComponents[tt["Real"]][comp], comp], {comp, {"XX", "NN"}}],
    "BornTable" -> Association @ Table[comp -> Map[#[comp] &, bornResults[n]["CrossedTable"]], {comp, {"XX", "NN"}}],
    "Class" -> pseudoName[n], "Spec" -> pair["Spec"], "Spins" -> {"XX", "NN"}, "Families" -> tt["Families"],
    "Variables" -> tt["Variables"] <> "; crossed from HadronicTT class " <> pair["TT"] <> " with " <> ToString[pair["Rules"], InputForm] <> ", factor " <> ToString[pair["Sign"] pair["Factor"], InputForm] <> " and the signs XX: -1, NN: +1",
    "Origin" -> "Processes/HadronicTransfer/5-Crossing.wl"|>;
  SaveResult["Algebra", algebra, route <> "/class-" <> pseudoName[n]]];
Do[writePseudoClass[loaded[n]], {n, pairList}];
"pseudo-classes written: " <> StringRiffle[pseudoName /@ pairList, ", "]


(* ::Input:: *)
pseudoResultsValidQ[n_, started_: None] := Module[{choices, files, receipts, algebraFile, reducedFile, resultFile,
    algebraHash, reducedHash, resultHash, root, relative, bound},
  choices = {{"route", route}, {"class", pseudoName[n]}};
  files = FileNameJoin[{NotebookDirectory[], #}] & /@ {"3-Reduce.wl", "4-Transfer.wl"};
  If[!And @@ (StepReceiptValidQ[#, choices, started] & /@ files), Return[False]];
  receipts = Import[FileNameJoin[{ResultsDirectory[], FileBaseName[#] <> "-" <> route <> "-" <> pseudoName[n] <> ".receipt.json"}], "RawJSON"] & /@ files;
  {algebraFile, reducedFile, resultFile} = transferFile[pseudoName[n], #] & /@ {"Algebra", "Reduced", "Transfer"};
  If[!And @@ (FileExistsQ /@ {algebraFile, reducedFile, resultFile}), Return[False]];
  If[NumberQ[started] && !And @@ (AbsoluteTime[FileDate[#]] >= Floor[started] & /@ {reducedFile, resultFile}), Return[False]];
  {algebraHash, reducedHash, resultHash} = IntegerString[FileHash[#, "SHA256"], 16, 64] & /@ {algebraFile, reducedFile, resultFile};
  root = StringReplace[ExpandFileName[$FeynFacetDirectory], $PathnameSeparator ~~ EndOfString -> ""];
  relative[file_] := StringDrop[ExpandFileName[file], StringLength[root] + 1];
  bound[receipt_, kind_, file_, hash_] := Lookup[receipt["provenance"][kind], relative[file], None] === hash;
  bound[receipts[[1]], "inputs", algebraFile, algebraHash] &&
    bound[receipts[[1]], "outputs", reducedFile, reducedHash] &&
    bound[receipts[[2]], "inputs", algebraFile, algebraHash] &&
    bound[receipts[[2]], "inputs", reducedFile, reducedHash] &&
    bound[receipts[[2]], "outputs", resultFile, resultHash]];
runPseudoClass[n_] := Module[{summary, rows, started},
  If[reuse === "yes" && pseudoResultsValidQ[n], Return["reused"]];
  started = AbsoluteTime[];
  summary = RunSteps["HadronicTransfer", route, {3, 4}, "class=" <> pseudoName[n]];
  rows = If[Head[summary] === Grid, Rest[First[summary]], {}];
  Require[MatchQ[rows, {{"3-Reduce", "passed", _}, {"4-Transfer", "passed", _}}],
    "Steps 3 and 4 of the pseudo-class must both finish successfully", {pseudoName[n], summary}];
  Require[pseudoResultsValidQ[n, started],
    "The pseudo-class needs fresh receipts binding Algebra, Reduced and Transfer", pseudoName[n]];
  summary];
pseudoRuns = If[stage === "all", Association @ Table[n -> runPseudoClass[n], {n, pairList}], "stage A: steps 3 and 4 of the pseudo-classes not run"];
pseudoRuns


(* ::Section:: *)
(*Part B: the NLO comparison*)


(* ::Text:: *)
(*The hard coefficients of the pseudo-class (Transfer.wl of class <n>x) against those of the directly computed class <n>, both components: the Born term (exactly), the finite real-emission tuple and the finite virtual term separately, the counterterm tuples, and the complete finite coefficients {regular, D0, D1, delta} (with the SU(N) colour relation; the comparison is exact at rational points of s, v, w, Nf with the three scales symbolic, ExactDifference; a numerical difference to 30 digits is recorded as well).  The difference of the two components of the pseudo-class equals, by construction, -2 sigma g times the crossed evanescent Symmetric projection of HadronicTT assembled with the transfer counterterms; whether it vanishes (as the TT Symmetric remainder does) is reported with the result, evaluated to 30 digits at the physical points (the expression holds logarithms of rational numbers that Simplify does not combine).  The pole cancellation of the pseudo-class was required by step 4 itself; its receipt is recorded.*)


(* ::Input:: *)
scaleAssumptions = muR > 0 && muF > 0 && muFF > 0 && s > 0;
exactPoints = Map[Select[#, FreeQ[#[[1]], muR | muF | muFF] &] &, PhysicalPoints];            (* s, v, w, Nf rational; the scales symbolic *)
(* the difference of two scalars or tuples at one point: Together, PowerExpand of the scale logarithms, Simplify, and FullSimplify only if needed *)
differenceAt[a_, b_, point_] := Module[{d = Together[(a - b) /. colourRelation /. point]},
  d = Simplify[PowerExpand[d, Assumptions -> scaleAssumptions], scaleAssumptions];
  If[d =!= 0, d = FullSimplify[FunctionExpand[d], scaleAssumptions]]; d];
exactDifferences[a_, b_] := Table[differenceAt[x[[1]], x[[2]], point], {point, exactPoints}, {x, Transpose[{Flatten[{a}], Flatten[{b}]}]}];
nloCheck[n_] := Module[{ours, theirs, receipt, born, real, virtual, counterterms, finite, numeric, componentDifference},
  Require[FileExistsQ[transferFile[n, "Transfer"]], "Transfer.wl of class " <> n <> " (step 4 of this process) is needed for the NLO comparison"];
  ours = Get[transferFile[pseudoName[n], "Transfer"]]; theirs = Get[transferFile[n, "Transfer"]];
  receipt = With[{file = FileNameJoin[{ResultsDirectory[], "4-Transfer-" <> route <> "-" <> pseudoName[n] <> ".receipt.json"}]}, If[FileExistsQ[file], Import[file, "RawJSON"]["completed"], "no receipt"]];
  born = Association @ Table[spin -> Simplify[(ours["Born"][spin] - theirs["Born"][spin]) /. colourRelation], {spin, {"XX", "NN"}}];
  Require[And @@ (# === 0 & /@ Values[born]), "The Born terms differ", {n, born}];
  real = Association @ Table[spin -> exactDifferences[ours["Real"][spin], theirs["Real"][spin]], {spin, {"XX", "NN"}}];
  virtual = Association @ Table[spin -> exactDifferences[ours["Virtual"][spin], theirs["Virtual"][spin]], {spin, {"XX", "NN"}}];
  counterterms = Association @ Table[spin -> Association @ Table[key -> exactDifferences[ours["Counterterms"][spin][key], theirs["Counterterms"][spin][key]], {key, Keys[theirs["Counterterms"][spin]]}], {spin, {"XX", "NN"}}];
  finite = Association @ Table[spin -> exactDifferences[ours["Finite"][spin], theirs["Finite"][spin]], {spin, {"XX", "NN"}}];
  numeric = Association @ Table[spin -> NumericDifference[ours["Finite"][spin] /. colourRelation /. CA -> 3, theirs["Finite"][spin] /. colourRelation /. CA -> 3, PhysicalPoints, 30], {spin, {"XX", "NN"}}];
  (* XX - NN of the pseudo-class (= that of the class, by the comparison above), at 30 digits: its logarithms of numbers are not combined by Simplify *)
  componentDifference = NumericDifference[ours["ComponentDifference"] /. colourRelation /. CA -> 3, 0 ours["ComponentDifference"], PhysicalPoints, 30];
  Require[And @@ (# === 0 & /@ Flatten[Values[real]]), "The finite real-emission parts differ", {n, real}];
  Require[And @@ (# === 0 & /@ Flatten[Values[virtual]]), "The finite virtual parts differ", {n, virtual}];
  Require[And @@ (# === 0 & /@ Flatten[Values[Map[Values, counterterms]]]), "The counterterms differ", {n, counterterms}];
  Require[And @@ (# === 0 & /@ Flatten[Values[finite]]), "The finite hard coefficients of the pseudo-class differ from the directly computed ones", {n, finite}];
  <|"PseudoClass" -> pseudoName[n], "Step4Completed" -> receipt, "Born" -> born, "Real" -> real, "Virtual" -> virtual, "Counterterms" -> counterterms, "Finite" -> finite,
    "NumericFinite" -> numeric, "ComponentDifference" -> componentDifference, "Points" -> exactPoints|>];
partB = If[stage === "all", Association @ Table[n -> nloCheck[n], {n, pairList}], "stage A: no NLO comparison"];
If[AssociationQ[partB], Column[Table[loaded[n]["Name"] <> ": NLO finite coefficients of XX and NN agree exactly at " <> ToString[Length[exactPoints]] <> " points; poles of the pseudo-class cancelled (step 4 completed: " <>
  ToString[partB[n]["Step4Completed"]] <> "); XX - NN of the transfer class at three points (30 digits): " <> ToString[partB[n]["ComponentDifference"], InputForm], {n, pairList}]], partB]


(* ::Section:: *)
(*Save and summary*)


(* ::Input:: *)
Do[SaveResult["class-" <> n, <|"PartA" -> partA[n], "PartB" -> If[AssociationQ[partB], partB[n], partB], "Stage" -> stage|>, route <> "/crossing"], {n, pairList}];
summary = Association @ Table[With[{file = FileNameJoin[{crossingDirectory, "class-" <> n <> ".wl"}]}, n -> If[FileExistsQ[file], With[{r = Get[file]},
  <|"Pair" -> r["PartA"]["Pair"], "Factor" -> r["PartA"]["Sign"] r["PartA"]["Factor"], "Part A" -> "passed", "Part B" -> If[AssociationQ[r["PartB"]], "passed (finite XX, NN equal)", r["PartB"]]|>], "not run"]], {n, Keys[crossingTable]}];
Dataset[summary]
