(* classes.wl: the 16 tagged channel classes of single-inclusive parton
   production in proton-proton collisions, a(pa) b(pb) -> c(p) + X, with the
   numbering of the reference codes (J0) and of the paper.

   Particle names: q = quark of flavour i, qp = quark of a different flavour
   k, qb and qpb their antiquarks, g = gluon.  Each class lists
     "Born"     the 2 -> 2 process a b -> c d (None for the classes that start
                at order alpha_s^3), with the unobserved parton d,
     "Real"     the 2 -> 3 processes a b -> c d e with their weights:
                1/2 for two identical unobserved partons (gluons); 1 when the
                tagged parton has an identical unobserved partner (the
                identical-particle factor 1/2 is cancelled by the two ways
                of tagging); Nf - 1 for an unobserved pair of a quark flavour
                different from i (one representative generated); Nf for an
                unobserved quark pair when no flavour is singled out,
     "Virtual"  whether the one-loop correction to the Born process is needed,
     "Average"  the colour normalization CC of the reference codes (for the
                comparison only).
   The flavour sums of the hadronic cross section (over i, over k != i,
   quark <-> antiquark, incoming-leg exchange) are outside these classes. *)

HadronicClasses = <|
  1 -> <|"Name" -> "q qp -> q", "Born" -> {"q", "qp", "q", "qp"}, "Virtual" -> True,
    "Real" -> {{"q", "qp", "q", "qp", "g"} -> 1}|>,
  2 -> <|"Name" -> "q qp -> g", "Born" -> None, "Virtual" -> False,
    "Real" -> {{"q", "qp", "g", "q", "qp"} -> 1}|>,
  3 -> <|"Name" -> "q qpb -> q", "Born" -> {"q", "qpb", "q", "qpb"}, "Virtual" -> True,
    "Real" -> {{"q", "qpb", "q", "qpb", "g"} -> 1}|>,
  4 -> <|"Name" -> "q qpb -> g", "Born" -> None, "Virtual" -> False,
    "Real" -> {{"q", "qpb", "g", "q", "qpb"} -> 1}|>,
  5 -> <|"Name" -> "q qb -> qp", "Born" -> {"q", "qb", "qp", "qpb"}, "Virtual" -> True,
    "Real" -> {{"q", "qb", "qp", "qpb", "g"} -> 1}|>,
  6 -> <|"Name" -> "q q -> q", "Born" -> {"q", "q", "q", "q"}, "Virtual" -> True,
    "Real" -> {{"q", "q", "q", "q", "g"} -> 1}|>,
  7 -> <|"Name" -> "q q -> g", "Born" -> None, "Virtual" -> False,
    "Real" -> {{"q", "q", "g", "q", "q"} -> 1/2}|>,
  8 -> <|"Name" -> "q g -> qp", "Born" -> None, "Virtual" -> False,
    "Real" -> {{"q", "g", "qp", "qpb", "q"} -> 1}|>,
  9 -> <|"Name" -> "q g -> qpb", "Born" -> None, "Virtual" -> False,
    "Real" -> {{"q", "g", "qpb", "qp", "q"} -> 1}|>,
  10 -> <|"Name" -> "q g -> qb", "Born" -> None, "Virtual" -> False,
    "Real" -> {{"q", "g", "qb", "q", "q"} -> 1/2}|>,
  11 -> <|"Name" -> "q qb -> q", "Born" -> {"q", "qb", "q", "qb"}, "Virtual" -> True,
    "Real" -> {{"q", "qb", "q", "qb", "g"} -> 1}|>,
  12 -> <|"Name" -> "q qb -> g", "Born" -> {"q", "qb", "g", "g"}, "Virtual" -> True,
    "Real" -> {{"q", "qb", "g", "g", "g"} -> 1/2, {"q", "qb", "g", "q", "qb"} -> 1, {"q", "qb", "g", "qp", "qpb"} -> Nf - 1}|>,
  13 -> <|"Name" -> "q g -> q", "Born" -> {"q", "g", "q", "g"}, "Virtual" -> True,
    "Real" -> {{"q", "g", "q", "g", "g"} -> 1/2, {"q", "g", "q", "q", "qb"} -> 1, {"q", "g", "q", "qp", "qpb"} -> Nf - 1}|>,
  14 -> <|"Name" -> "q g -> g", "Born" -> {"q", "g", "g", "q"}, "Virtual" -> True,
    "Real" -> {{"q", "g", "g", "q", "g"} -> 1}|>,
  15 -> <|"Name" -> "g g -> g", "Born" -> {"g", "g", "g", "g"}, "Virtual" -> True,
    "Real" -> {{"g", "g", "g", "g", "g"} -> 1/2, {"g", "g", "g", "q", "qb"} -> Nf}|>,
  16 -> <|"Name" -> "g g -> q", "Born" -> {"g", "g", "q", "qb"}, "Virtual" -> True,
    "Real" -> {{"g", "g", "q", "qb", "g"} -> 1}|>|>;
