(* SSV-born.wl: the lowest-order transversely polarized 2 -> 2 cross sections,
   transcribed verbatim from J. Soffer, M. Stratmann and W. Vogelsang,
   "Accessing transversity in double-spin asymmetries at the BNL-RHIC",
   Phys. Rev. D 65 (2002) 114024 [hep-ph/0204058], Table I and Eq. (6).

   Eq. (6) of the paper:

       d^2 delta sigma-hat/(dt dPhi) = (2 alpha_s^2(mu_R)/(9 s^2)) (t u/s^2) cos(2 Phi) delta|M|^2,

   with Phi the azimuthal angle of the observed parton with respect to the
   (common) initial spin axis, and t = (p_a - p_c)^2 for the first-listed
   final-state parton c.  Table I (jets and inclusive hadrons):

       q q   -> q q      :  2 s^2/(3 t u)
       q q'  -> q q'     :  --  (vanishes)
       q qb  -> q qb     :  2 - 2 s/(3 t)
       q qb  -> q' qb'   :  2
       q qb' -> q qb'    :  --  (vanishes)
       q qb  -> g g      :  16 s^2/(3 t u) - 12

   The entries are for N_c = 3 (the 2/9 of the prefactor contains the colour average).

   The association SSVMatrixElement holds delta|M|^2 keyed by our channel names
   (a b -> c, c the observed parton); the two entries of q qb -> q qb are
   related by t <-> u (observing the antiquark instead of the quark), and the
   a <-> b relabelled channels likewise.  SSVTransverseBorn[name, s, t, u] is
   the cos(2 Phi) coefficient of d^2 delta sigma-hat/(dt dPhi) divided by
   alpha_s^2, Eq. (6) without the coupling. *)

SSVMatrixElement = <|
  "q q -> q" -> 2 s^2/(3 t u),
  "q qp -> q" -> 0, "q qp -> qp" -> 0,
  "q qb -> q" -> 2 - 2 s/(3 t), "q qb -> qb" -> 2 - 2 s/(3 u),
  "q qb -> qp" -> 2, "q qb -> qpb" -> 2,
  "q qpb -> q" -> 0, "q qpb -> qpb" -> 0,
  "q qb -> g" -> 16 s^2/(3 t u) - 12,
  (* a <-> b relabelled (t <-> u) *)
  "qb q -> q" -> 2 - 2 s/(3 u), "qb q -> qb" -> 2 - 2 s/(3 t),
  "qb q -> qp" -> 2, "qb q -> qpb" -> 2,
  "qb q -> g" -> 16 s^2/(3 t u) - 12,
  "qb qb -> qb" -> 2 s^2/(3 t u)|>;

SSVTransverseBorn[name_String, sIn_, tIn_, uIn_] := (2/(9 sIn^2)) (tIn uIn/sIn^2) (SSVMatrixElement[name] /. {s -> sIn, t -> tIn, u -> uIn});
