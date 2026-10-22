(* classes.wl: the channel classes of the transverse spin transfer in
   proton-proton collisions, a(pa)^ b(pb) -> c(p)^ + X: the incoming parton a
   transversely polarized, the transverse spin of the observed parton c
   measured (docs/TRANSVERSE-PP.md).  The numbering is that of
   Processes/Hadronic/classes.wl (and of the reference codes), plus class 17.

   Particle names and the meaning of "Born", "Real", "Virtual" and the
   real-emission weights are those of Processes/Hadronic/classes.wl: q = quark
   of flavour i, qp = quark of a different flavour k, qb and qpb their
   antiquarks, g = gluon; weight 1/2 for two identical unobserved partons, 1
   when the tagged parton has an identical unobserved partner, Nf - 1 for an
   unobserved pair of another flavour.

   Why only these seven of the 16 + 1 classes contribute.  Transversity is a
   chirality-odd quantity: a massless quark line conserves helicity, and the
   transverse spin of a is transferred to c only if a and c sit on one fermion
   line (the t-channel exchange classes 1, 3, 6, 11, 13) or if they sit on two
   lines that are joined by the interference of an annihilation graph with an
   exchange graph (classes 6, 11 and 17, where a and c are on the same line in
   one graph and on different lines in the other; the NLO-only class 10, where
   the observed antiquark's line is linked to the polarized quark's line by the
   interference of graphs in which the two identical unobserved quarks are
   exchanged).  A gluon carries no
   transversity (helicity flip 2 cannot feed a quark's helicity flip 1), so
   every class with an incoming gluon on leg a (14, 15, 16 and the g-initiated
   halves of the others) or with an observed gluon (2, 4, 7, 12, 14, 15)
   vanishes identically; a quark of another flavour observed (5, 8, 9) vanishes
   because a and c are then on different lines without an annihilation graph
   between them (q qb -> qp qpb has no exchange graph to interfere with).
   The same reasoning gives the classes of the double-transverse asymmetry
   (Processes/HadronicTT); the two sets are related by the generalized
   crossing of docs/TRANSVERSE-PP.md:

     transfer class   process               crossing partner in HadronicTT
     1                q qp  -> q + X        class 5  (q qb -> qp, observing qpb)
     3                q qpb -> q + X        class 5  (observing qp); Born term equal to class 1
     6                q q   -> q + X        class 11 with the observed antiquark
     10               q g   -> qb + X       class 7  (NLO only)
     11               q qb  -> q + X        class 11
     13               q g   -> q + X        class 12
     17               qb q  -> q + X        class 6

   Class 17 is q qb -> qb with the quark polarized and the antiquark observed,
   written by charge conjugation with the polarized antiquark on leg a and the
   observed quark, so that the observed parton is a quark as in the other
   classes (the counterterm bookkeeping of step 4 reads the species of a and c).

   The flavour sums of the hadronic cross section (over i, over k != i,
   quark <-> antiquark, incoming-leg exchange) are outside these classes. *)

HadronicTransferClasses = <|
  1 -> <|"Name" -> "q qp -> q", "Born" -> {"q", "qp", "q", "qp"}, "Virtual" -> True,
    "Real" -> {{"q", "qp", "q", "qp", "g"} -> 1}|>,
  3 -> <|"Name" -> "q qpb -> q", "Born" -> {"q", "qpb", "q", "qpb"}, "Virtual" -> True,
    "Real" -> {{"q", "qpb", "q", "qpb", "g"} -> 1}|>,
  6 -> <|"Name" -> "q q -> q", "Born" -> {"q", "q", "q", "q"}, "Virtual" -> True,
    "Real" -> {{"q", "q", "q", "q", "g"} -> 1}|>,
  10 -> <|"Name" -> "q g -> qb", "Born" -> None, "Virtual" -> False,
    "Real" -> {{"q", "g", "qb", "q", "q"} -> 1/2}|>,
  11 -> <|"Name" -> "q qb -> q", "Born" -> {"q", "qb", "q", "qb"}, "Virtual" -> True,
    "Real" -> {{"q", "qb", "q", "qb", "g"} -> 1}|>,
  13 -> <|"Name" -> "q g -> q", "Born" -> {"q", "g", "q", "g"}, "Virtual" -> True,
    "Real" -> {{"q", "g", "q", "g", "g"} -> 1/2, {"q", "g", "q", "q", "qb"} -> 1, {"q", "g", "q", "qp", "qpb"} -> Nf - 1}|>,
  17 -> <|"Name" -> "qb q -> q", "Born" -> {"qb", "q", "q", "qb"}, "Virtual" -> True,
    "Real" -> {{"qb", "q", "q", "qb", "g"} -> 1}|>|>;

(* The 2 -> 2 transfer channels a^ b -> c^ whose D-dimensional squares the
   collinear counterterms of step 4 need (the lower Born of every class after
   a splitting on leg a, on leg b or of the fragmenting parton), computed in
   step 2 with the same machinery as the class squares.  The polarized parton
   is the first one, the observed parton the one after the arrow. *)
TransferBornChannels = {{"q", "qp", "q"}, {"q", "qpb", "q"}, {"q", "q", "q"}, {"q", "qb", "q"}, {"q", "qb", "qb"}, {"q", "g", "q"}, {"qb", "q", "q"}};
