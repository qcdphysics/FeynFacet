(* classes.wl: the channel classes of single-inclusive parton production in
   collisions of two transversely polarized protons, p^ p^ -> h + X (the
   asymmetry A_TT), a(pa) b(pb) -> c(p) + X, with the numbering of
   Processes/Hadronic/classes.wl and of the reference codes (J0).

   The entries have the same meaning as in Processes/Hadronic/classes.wl:
     "Born"     the 2 -> 2 process a b -> c d (None for class 7, which starts
                at order alpha_s^3), with the unobserved parton d,
     "Real"     the 2 -> 3 processes a b -> c d e with their weights (1/2 for
                two identical unobserved partons, 1 when the tagged parton has
                an identical unobserved partner, Nf - 1 for an unobserved pair
                of another flavour),
     "Virtual"  whether the one-loop correction to the Born process is needed,
     "Average"  the colour normalization CC of the reference codes.

   Why only five of the sixteen classes contribute.  Both incoming partons
   carry a transverse spin, which enters through the chiral-odd densities
   pa-slash gamma5 sa-slash/2 and pb-slash gamma5 sb-slash/2 (docs/TRANSVERSE-PP.md).
   A gluon has no transversity (helicity flip of two units), so every
   gluon-initiated class (8, 9, 10, 13, 14, 15, 16) vanishes.  For two incoming
   quarks or antiquarks, each density flips the helicity of its line; along a
   massless fermion line with vector couplings the helicity is conserved, so
   the squared amplitude is non-zero only if the two polarized partons sit on
   one closed fermion line of the cut diagram: the s-channel annihilation
   line (q qb -> q' qb', q qb -> q qb, q qb -> g g, with the two gluons
   fragmenting, class 12, or the quark, classes 5 and 11) or the exchange
   interference of two identical quarks (q q -> q q, class 6, and at NLO
   q q -> g q q, class 7).  Classes 1-4 (q q' and q qb' scattering, different
   flavours) have the two polarized partons on two distinct lines connected
   only by t-channel gluons, and vanish identically at every order; this is
   verified on the Born table of step 2 (control channels q q' -> q and
   q qb' -> q).  Class 2 of Processes/Hadronic (q q' -> g) and class 4 vanish
   for the same reason. *)

HadronicTTClasses = <|
  5 -> <|"Name" -> "q qb -> qp", "Born" -> {"q", "qb", "qp", "qpb"}, "Virtual" -> True,
    "Real" -> {{"q", "qb", "qp", "qpb", "g"} -> 1}|>,
  6 -> <|"Name" -> "q q -> q", "Born" -> {"q", "q", "q", "q"}, "Virtual" -> True,
    "Real" -> {{"q", "q", "q", "q", "g"} -> 1}|>,
  7 -> <|"Name" -> "q q -> g", "Born" -> None, "Virtual" -> False,
    "Real" -> {{"q", "q", "g", "q", "q"} -> 1/2}|>,
  11 -> <|"Name" -> "q qb -> q", "Born" -> {"q", "qb", "q", "qb"}, "Virtual" -> True,
    "Real" -> {{"q", "qb", "q", "qb", "g"} -> 1}|>,
  12 -> <|"Name" -> "q qb -> g", "Born" -> {"q", "qb", "g", "g"}, "Virtual" -> True,
    "Real" -> {{"q", "qb", "g", "g", "g"} -> 1/2, {"q", "qb", "g", "q", "qb"} -> 1, {"q", "qb", "g", "qp", "qpb"} -> Nf - 1}|>|>;
