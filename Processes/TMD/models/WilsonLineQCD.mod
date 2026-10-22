(* WilsonLineQCD.mod: classes model "massless QCD with light-like Wilson lines
   and TMD operators", to be used with the generic model WilsonLine.gen.

   Fields
     F[1]  q      massless quark (one flavour; Nf enters only through the
                  gluon PDF counterterm, which the notebooks write by hand)
     V[1]  g      gluon
     U[1]  u_g    gluon ghost (not needed at the order computed here)
     S[1]  W_n    Wilson line in the direction n,  fundamental representation
     S[2]  W_nb   Wilson line in the direction nb, fundamental representation
     S[3]  W_n^A  Wilson line in the direction n,  adjoint representation
     S[4]  W_nb^A Wilson line in the direction nb, adjoint representation
     F[2]  O_q    quark TMD operator  (external sink, Dirac index open)
     V[2]  O_g    gluon TMD operator  (external sink, Lorentz index open)
     S[5]  O_S    soft operator       (external sink)
   The light-cone vectors n and nb (n^2 = nb^2 = 0, n.nb = 2) are the Global`
   symbols n and nb of the notebooks.

   Couplings (gs = strong coupling; gOp and gSoft mark the operator
   insertions and are divided out at the end; they are not interactions):
     quark-gluon           -i gs T^a gamma^mu                  (standard QCD)
     3 gluons, 4 gluons, ghosts                               (standard QCD)
     Wilson line-gluon     -i gs T^a n^mu  (fundamental),  -gs f^{abc} n^mu (adjoint)
     O_q q W_n             i gOp delta_{ij}:  the operator O_q-bar(x) W_n-bar(x) q(x)
     O_g g W_n^A           gOp (n.q g^{mu rho} - q^mu n^rho) delta^{ab}  and
     O_g g g W_n^A         -i gOp gs f^{abc} (n^rho g^{mu sigma} - n^sigma g^{mu rho}):
                           the operator O_g^mu(x) W_n^A(x)^{ab} n_nu F^{nu mu, b}(x)
                           with F^{nu mu} = d^nu A^mu - d^mu A^nu + gs f A^nu A^mu
     O_S W_n W_nb          i gSoft delta_{ij}:  the operator W_n-bar(x) W_nb(x)
   The sign of the eikonal vertex relative to the quark-gluon vertex, and the
   signs of the gluon-operator couplings, are fixed by gauge invariance: step 1
   of the TMD process checks the Ward identities of every generated current. *)

IndexRange[Index[Colour]] = NoUnfold[Range[3]];
IndexRange[Index[Gluon]] = NoUnfold[Range[8]];

M$ClassesDescription = {
  F[1] == {SelfConjugate -> False, Indices -> {Index[Colour]}, Mass -> 0,
    PropagatorLabel -> "q", PropagatorType -> Straight, PropagatorArrow -> Forward},
  V[1] == {SelfConjugate -> True, Indices -> {Index[Gluon]}, Mass -> 0,
    PropagatorLabel -> "g", PropagatorType -> Cycles, PropagatorArrow -> None},
  U[1] == {SelfConjugate -> False, Indices -> {Index[Gluon]}, Mass -> 0,
    PropagatorLabel -> ComposedChar["u", "g"], PropagatorType -> GhostDash, PropagatorArrow -> Forward},

  S[1] == {SelfConjugate -> False, Indices -> {Index[Colour]}, Mass -> 0,
    PropagatorLabel -> ComposedChar["W", "n"], PropagatorType -> ScalarDash, PropagatorArrow -> Forward},
  S[2] == {SelfConjugate -> False, Indices -> {Index[Colour]}, Mass -> 0,
    PropagatorLabel -> ComposedChar["W", "n\:0304"], PropagatorType -> ScalarDash, PropagatorArrow -> Forward},
  S[3] == {SelfConjugate -> False, Indices -> {Index[Gluon]}, Mass -> 0,
    PropagatorLabel -> ComposedChar["W", "n", "A"], PropagatorType -> ScalarDash, PropagatorArrow -> Forward},
  S[4] == {SelfConjugate -> False, Indices -> {Index[Gluon]}, Mass -> 0,
    PropagatorLabel -> ComposedChar["W", "n\:0304", "A"], PropagatorType -> ScalarDash, PropagatorArrow -> Forward},

  F[2] == {SelfConjugate -> False, Indices -> {}, Mass -> 0,
    PropagatorLabel -> ComposedChar["O", "q"], PropagatorType -> Straight, PropagatorArrow -> Forward},
  V[2] == {SelfConjugate -> True, Indices -> {}, Mass -> 0,
    PropagatorLabel -> ComposedChar["O", "g"], PropagatorType -> Sine, PropagatorArrow -> None},
  S[5] == {SelfConjugate -> True, Indices -> {}, Mass -> 0,
    PropagatorLabel -> ComposedChar["O", "S"], PropagatorType -> Straight, PropagatorArrow -> None}
};

FAGaugeXi[V[1, ___]] = FAGaugeXi[g];
FAGaugeXi[U[1, ___]] = FAGaugeXi[g];

(* direction of each Wilson line, and the orientation of a line on an
   internal propagator: +1 if the arrow runs along the propagator momentum,
   -1 (field -S) if against; see WilsonLine.gen *)
Global`WilsonDirection[S[1 | 3, ___]] = Global`n;
Global`WilsonDirection[S[2 | 4, ___]] = Global`nb;
Global`WilsonDirection[-field_] := Global`WilsonDirection[field];
Global`WilsonLineSign[S[_Integer, ___]] = 1;
Global`WilsonLineSign[-field_] := -Global`WilsonLineSign[field];

(* The index names col1, col2, glu1 ... glu4 in the couplings are pattern
   symbols of FeynArts; do not assign values to them. *)
M$CouplingMatrices = {

  (* ---- QCD (as in FeynArts' SMQCD.mod, massless, no counterterms) ---- *)
  C[-F[1, {col1}], F[1, {col2}], V[1, {glu1}]] == -I gs FASUNT[glu1, col1, col2] {{1}, {1}},

  C[V[1, {glu1}], V[1, {glu2}], V[1, {glu3}]] == gs FASUNF[glu1, glu2, glu3] {{1}},

  C[V[1, {glu1}], V[1, {glu2}], V[1, {glu3}], V[1, {glu4}]] == -I gs^2 *
    {(FASUNF[glu1, glu3, glu2, glu4] - FASUNF[glu1, glu4, glu3, glu2]) {1},
     (FASUNF[glu1, glu2, glu3, glu4] + FASUNF[glu1, glu4, glu3, glu2]) {1},
     (-FASUNF[glu1, glu2, glu3, glu4] - FASUNF[glu1, glu3, glu2, glu4]) {1}},

  C[-U[1, {glu1}], U[1, {glu2}], V[1, {glu3}]] == gs FASUNF[glu1, glu2, glu3] {{1}, {0}},

  (* ---- Wilson lines ---- *)
  C[-S[1, {col1}], S[1, {col2}], V[1, {glu1}]] == -I gs FASUNT[glu1, col1, col2] {{1}},
  C[-S[2, {col1}], S[2, {col2}], V[1, {glu1}]] == -I gs FASUNT[glu1, col1, col2] {{1}},
  C[-S[3, {glu1}], S[3, {glu2}], V[1, {glu3}]] == -gs FASUNF[glu3, glu1, glu2] {{1}},
  C[-S[4, {glu1}], S[4, {glu2}], V[1, {glu3}]] == -gs FASUNF[glu3, glu1, glu2] {{1}},

  (* ---- operators ---- *)
  (* quark operator: the quark enters, the Wilson line W_n leaves *)
  C[-F[2], F[1, {col1}], -S[1, {col2}]] == I gOp IndexDelta[col1, col2] {{1}, {1}},
  C[-F[1, {col1}], F[2], S[1, {col2}]] == I gOp IndexDelta[col1, col2] {{1}, {1}},

  (* gluon operator: n_nu F^{nu mu}, abelian and non-abelian parts *)
  C[V[1, {glu1}], V[2], -S[3, {glu2}]] == gOp IndexDelta[glu1, glu2] {{1}, {0}},
  C[V[1, {glu1}], V[2], S[3, {glu2}]] == gOp IndexDelta[glu1, glu2] {{1}, {0}},
  C[V[1, {glu1}], V[1, {glu2}], V[2], -S[3, {glu3}]] == -I gOp gs FASUNF[glu1, glu2, glu3] {{0}, {-1}, {1}},
  C[V[1, {glu1}], V[1, {glu2}], V[2], S[3, {glu3}]] == -I gOp gs FASUNF[glu1, glu2, glu3] {{0}, {-1}, {1}},

  (* soft operator: the colour charge enters along nb and leaves along n *)
  C[S[5], -S[1, {col1}], S[2, {col2}]] == I gSoft IndexDelta[col1, col2] {{1}},
  C[S[5], -S[3, {glu1}], S[4, {glu2}]] == I gSoft IndexDelta[glu1, glu2] {{1}}
};

M$LastModelRules = {};
