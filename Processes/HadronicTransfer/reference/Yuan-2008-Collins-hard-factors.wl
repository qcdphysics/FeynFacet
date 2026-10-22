(* Yuan-2008-Collins-hard-factors.wl: the leading-order hard factors of the
   Collins asymmetry in p^ p -> jet + X, which are the leading-order transverse
   spin transfer cross sections a^ b -> c^ of the partonic channels, transcribed
   verbatim from

     Feng Yuan, "Azimuthal asymmetric distribution of hadrons inside a jet at
     hadron collider", Phys. Rev. Lett. 100 (2008) 032003 [arXiv:0709.3272],
     Eq. (6) (arXiv source collins1.tex, lines 190-200):

   \begin{eqnarray}
   H_{qq'\to qq'}^{\rm Collins}&=&H_{q\bar q'\to q\bar q'}^{\rm
   Collins}=\frac{\alpha_s^2\pi}{\hat
   s^2}\frac{N_c^2-1}{4N_c^2}\frac{4\hat s\hat u}{-\hat t^2} \ ,
   \nonumber\\
   H_{qq\to qq}^{\rm Collins}&=&\frac{\alpha_s^2\pi}{\hat
   s^2}\frac{N_c^2-1}{4N_c^2}\left[\frac{4\hat s\hat u}{-\hat
   t^2}-\frac{1}{N_c}\frac{4\hat s}{-\hat t}\right] \ ,
   \nonumber\\
   H_{qg\to qg}^{\rm Collins}&=&\frac{\alpha_s^2\pi}{\hat
   s^2}\left[\frac{N_c^2-1}{N_c^2}+\frac{1}{2}\frac{4\hat s\hat
   u}{-\hat t^2}\right] \ ,
   \end{eqnarray}

   "where $\hat s$, $\hat t$, and $\hat u$ are the usual partonic Mandelstam
   variables. The hard factors for the partonic channels associated with the
   antiquark transversity are the same as the above. Not surprisingly, these
   hard factors are exactly the same as those calculated for the transverse
   spin transfer in the hard partonic processes [Stratmann and Vogelsang,
   Phys. Lett. B295 (1992) 277; Collins, Heppelmann and Ladinsky, Nucl. Phys.
   B420 (1994) 565]."  The unpolarized hard factors of the same paper are the
   differential cross sections, H^{uu} = d sigma/d t = alphaS^2 Pi/s^2 x
   |M|^2/gs^4 (averaged over the initial spins and colours).

   YuanCollinsHardFactors[s, t, u, Nc] returns h(s, t, u) with
   H = alphaS^2 Pi/s^2 x h for the channels, named as in classes.wl (first
   parton polarized, observed parton after the arrow). *)

YuanCollinsHardFactors[s_, t_, u_, Nc_] := <|
  "q qp -> q" -> (Nc^2 - 1)/(4 Nc^2) (4 s u/(-t^2)),
  "q qpb -> q" -> (Nc^2 - 1)/(4 Nc^2) (4 s u/(-t^2)),
  "q q -> q" -> (Nc^2 - 1)/(4 Nc^2) (4 s u/(-t^2) - (1/Nc) (4 s/(-t))),
  "q g -> q" -> ((Nc^2 - 1)/Nc^2 + (1/2) (4 s u/(-t^2)))|>;
