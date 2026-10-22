# Reference values from the Fortran codes

`pp-fortran-values.wl` (one level up) holds the NLO coefficient pieces of the
two reference codes for single-inclusive hadron production in pp collisions,
evaluated at a few rational points:

- `kang2.f` — unpolarized, F. Aversa, P. Chiappetta, M. Greco and
  J.-P. Guillet, Nucl. Phys. B327 (1989) 105 (the version used in
  arXiv:1505.06415);
- `ourh_p.f` — longitudinally polarized, B. Jäger, A. Schäfer, M. Stratmann
  and W. Vogelsang, Phys. Rev. D67 (2003) 054005.

The codes themselves are not part of this repository.  `driver.f` is the
program that produced the table; with gfortran it is built and run as

```sh
gfortran -std=legacy -o kang2 driver.f kang2.f -Wl,--unresolved-symbols=ignore-all
echo "1 0.3 0.6 1 0.7 0.9 1.3 4" | ./kang2      # J0 v w s muF^2 muR^2 muFF^2 Nf
```

(the unresolved symbols are the parton-distribution routines of the full
program, which the coefficient functions do not call).  Each output line is
`J0 F0 AVDELMU AVDEL AVWPL(w) AVWPL(1) AVLO(w) AVLO(1) Regular CC`, in the
notation of the codes: F0 the Born term, AVDEL the δ(1−w) coefficient, AVDELMU
the coefficient of ln(s/μR²) δ(1−w) (printed by the code already divided by
CC, unlike the other pieces), AVWPL and AVLO the coefficients of 1/(1−w) and
ln(1−w)/(1−w) (their values at w = 1 are the subtraction terms of the plus
distributions), Regular the rest, CC the colour normalization (8 Nc² for two
incoming quarks, 8 Nc (Nc²−1) for quark–gluon, 8 (Nc²−1)² for gluon–gluon).
The code's hadronic integrand carries the factor 1/w of the partonic
normalization σ = s E dσ/d³p in the Jacobian of its variables, so the pieces
compare with our coefficients without it.  `Processes/Hadronic/4-Unpolarized.wl` explains how the
generated coefficients are mapped onto these pieces.

## The transversely polarized code and the normalization map

`pp-transverse-fortran.wl` (one level up) holds, as Mathematica expressions,
the pieces of the third reference code,

- `me-transpol-ms.f` with `drive-transpol.f` — transversely polarized
  (A_TT), A. Mukherjee, M. Stratmann and W. Vogelsang, Phys. Rev. D72 (2005)
  034011 [hep-ph/0506315] (M. Stratmann, status 03/11/04),

transcribed verbatim by `transcribe_transpol.py` (this directory): the script
parses the fixed-form Fortran, resolves the precomputed quantities of the
common blocks to their definitions, writes the file and checks it by
evaluating the Fortran expressions in Python/mpmath and the Mathematica file
with wolframscript at five exact rational points (agreement to 1e-37), and by
parsing the raw FortranForm files `proc5/6/7/11/12-transpol.f` shipped with
the code (the same expressions before they were pasted into the routines;
not compiled by the makefile) as an independent second source.  The code
itself is not part of the repository; the script takes its directory as
argument.  The header of the .wl file documents the code's routines, line
ranges, the hadronic assembly `DPLUS` and the check.  The code's formulas are
written in CA and CF (the driver sets the SU(3) values CF = 4/3, CA = 3,
NF = 4); whether they hold for general CA is not decidable from the code
(the comparison below is done at general CA with CF = (CA²−1)/(2 CA)).

For the five classes with transverse polarization (J0 = 5, 6, 7, 11, 12 of
the code's numbering; the same as `Processes/Hadronic/classes.wl`) the pieces
are `F0` (Born, routine FBOR, including the phase-space factor 1/(v(1−v))),
`AVDEL` (δ(1−w), with ln(μR²/s) inside, no separate AVDELMU), `AVWPL`
(1/(1−w)₊), `AVLO` ([ln(1−w)/(1−w)]₊), `STRUV` (regular), and `CC` = 8 Nc².
AVWPL and AVLO depend on v only, so the plus-distribution subtraction terms of
the canonical form vanish and the regular part is STRUV alone.  The structure
is that of the paper's Eq. (10), s dδσ̂⁽¹⁾/(dv dw dΦ) = cos(2Φ)(αs/π)²
[(A0 δ(1−w) + B0/(1−w)₊ + C0) ln(μ²/s) + A δ(1−w) + B/(1−w)₊ + C
+ D [ln(1−w)/(1−w)]₊ + E ln w + F ln v + G ln(1−v) + H ln(1−w) + I ln(1−vw)
+ J ln(1−v+vw) + K ln w/(1−w) + L ln((1−v)/(1−vw))/(1−w) + M ln(1−v+vw)/(1−w)]:
AVDEL ↔ A, A0; AVWPL ↔ B, B0; AVLO ↔ D; STRUV ↔ C ... M (the code's LVW and
L1VW are the L and M structures).  The paper's overall (αs/π)² is not the
code's normalization; the latter is fixed by DPLUS as follows.

### Derivation of the map

1. Kinematics.  With t = −s(1−v), u = −s v w one has dt du = s² v dv dw and,
   for a massless parton, d³p/E = p_T dp_T dy dΦ = dt du dΦ/(2s) (so that
   dp_T² dy = s v dv dw).  Hence s E dσ/d³p = (2/v) dσ/(dv dw dΦ); for an
   azimuth-independent cross section s E dσ/d³p = (1/(π v)) dσ/(dv dw), and
   at Born level dσ̂/dv = s dσ̂/dt = |M|²/(16π s) gives the Born coefficient
   b = αs² |M|²/(g⁴ s v) of our σ = s E dσ̂/d³p.

2. The hadronic formula of the unpolarized and longitudinal codes.  From
   E dσ/d³p = Σ ∫dz/z² ∫dx_a dx_b f_a f_b D_c (E dσ̂/d³p)(ŝ, v, w) with
   x_a = VW/(v w z), x_b = (1−V)/((1−v) z) at fixed z, dx_a dx_b
   = x_a x_b dv dw/(w (1−v)), ŝ = x_a x_b S, one obtains
   E dσ/d³p = (1/(π S)) Σ ∫dz/z² ∫dv/(v(1−v)) ∫dw/w f_a f_b D_c
   [dσ̂⁽⁰⁾/dv δ(1−w) + (αs/π) dσ̂⁽¹⁾/(dv dw)] — the formula coded in DPLUS
   (`DPLUS = HC2/(π S) × Σ ...`), whose Born slot is therefore
   F0 = (1/(v(1−v))) dσ̂⁽⁰⁾/dv /αs² = PRELO × |M|²/g⁴ = π b/(1−v),
   PRELO = π/(ŝ v (1−v)) (FBOR: "incl. phase space factor 1/v/(1−v)").  This is
   the established map of `Processes/Hadronic/4-Unpolarized.wl`, together with
   {R, B, C, A} = αs³/(π s CC) × {Regular, AVWPL1, AVLO1, AVDEL(+CC AVDELMU
   ln(s/μR²))} for the NLO slots, whose hadronic weight equals the Born slot's
   (the explicit 1/W of FHOREST cancels against XJAC = W BXJAC).

3. The transverse code.  Its DPLUS is identical (same slots, same Jacobians,
   same 1/(πS), same plus prescription on [0, 1] in w subtracted at w = 1 at
   fixed hadronic variables), with two changes: the final factor 4/π ("include
   cos(2φ) dependence, four quadrants") and ALPASHO = αs × IORD × 0 (the
   shipped driver switches the NLO weight off).  For a Φ-dependent cross
   section the same change of variables gives
   d³δσ/(dp_T dη dΦ) = p_T E dδσ/d³p = (2 p_T/S) Σ ∫dz/z² ∫dv/(v(1−v)) ∫dw/w
   δf_a δf_b D_c [dδσ̂⁽⁰⁾/(dv dΦ) δ(1−w) + ...] with the true Φ-differential
   partonic cross sections, whereas the paper's Eq. (3) keeps the prefactor
   p_T/(π S) of the azimuth-integrated formula.  The slots of Eq. (3), (10)
   and of the code therefore hold some fixed multiple of the Φ-differential
   cross section, and the multiple is read off the code's Born:
   F0(5)/PRELO = (CF/Nc) v (1−v).

4. The literature LO cross section.  Soffer, Stratmann and Vogelsang, Phys.
   Rev. D65 (2002) 114024 [hep-ph/0204058], Eq. (prefactor) and Table I
   (which agree with Jaffe and Saito, Phys. Lett. B382 (1996) 165, up to an
   overall 2π in that paper's normalization): d²δσ̂/(dt dΦ) = (2αs²/(9 s²))
   (t u/s²) cos(2Φ) δ|M|², δ|M|²(q q̄ → q′ q̄′) = 2, with dδσ = ½[dσ(↑↑) −
   dσ(↑↓)] and δf = f↑ − f↓ — the same definitions as ours (the partonic
   densities ½ p̸ γ5 s̸ give ½[σ̂(s_a) − σ̂(−s_a)] per leg).  With t u/s² =
   v(1−v) and 2/9 = (CF/Nc)/2 at Nc = 3 the coefficient of cos(2Φ) in
   s E dδσ̂/d³p = 2s² d²δσ̂/(dt dΦ) δ(s+t+u) is 2αs² (CF/Nc)(1−v)/s δ(1−w):
   our b_TT(5) = 2αs² (CF/Nc)(1−v)/s, i.e. δ|M|²_TT/g⁴ = 2 (CF/Nc) v(1−v) in
   the normalization b = αs² |M|²/(g⁴ s v) of step 1.  The code's slot holds
   (CF/Nc) v(1−v): one half.  (Classes 6, 11, 12 agree with the table entries
   2s²/(3tu), 2 − 2s/(3t), 16s²/(3tu) − 12 with the same half.)  So the code's
   pieces are the cos(2Φ) moments ∫dΦ cos(2Φ) dδσ̂/dΦ = π × coefficient in
   the slots where the unpolarized code holds ∫dΦ dσ̂/dΦ = 2π × dσ̂/dΦ.  This
   is consistent with the 4/π: the code's convolution is then one half of the
   coefficient C of cos(2Φ) in E dδσ/d³p, and DPLUS (ISIGM = 3) = 2π p_T ×
   (C/2) × 4/π = 4 p_T C = the alternating-sign integral over the four
   quadrants of d³δσ/(dp_T dy dΦ) (∫ ± cos(2Φ) dΦ = 4, SSV's prescription),
   which divided by the full-azimuth unpolarized dσ/(dp_T dy) = 2π p_T
   E dσ/d³p is the asymmetry A_TT.

5. The map.  Born_TT = 2 (1−v) F0/π × αs², and, since the NLO slots carry
   the same weight as the Born slot (step 2, unchanged in this code), the
   same factor 2 for the NLO pieces:

       {R, B, C, A}_TT = 2 αs³/(π s CC) × {STRUV (+ vanishing subtraction
                         terms), AVWPL(1), AVLO(1), AVDEL},

   or, in the coefficient basis σ = Born δ(1−w) + P c_regular + P(v,1)
   [c_D0 D0 + c_D1 D1 + c_δ δ(1−w)], P = αs³/(π s v w):
   c_regular = 2 v w STRUV/CC, c_D0 = 2 v AVWPL(1)/CC, c_D1 = 2 v AVLO(1)/CC,
   c_δ = 2 v AVDEL/CC, with CC = 8 Nc².  `PPTransverseFortranCoefficients[J0]`
   in `pp-transverse-fortran.wl` returns <|"Born", "Regular", "D0", "D1",
   "Delta", "Prefactor", "Pieces" (the raw F0, AVDEL, AVWPL, AVLO, STRUV, CC)|>
   in the variables v, w, s, muR, muF, muFF (scale logs Log[mu²/s]), CA, CF,
   Nf, alphaS — the API for `Processes/HadronicTT/4-Transverse.wl`.

6. Checks of the factor 2 that do not use the Born of the code:
   (a) renormalization scale: c_δ ∋ (β0 v/2) × (s b/αs²) LR, β0 = (11CA − 2Nf)/3;
   the code's 2 v × (4/3) CA CF (1−v)(11CA − 2Nf)/CC = CF v(1−v) β0/CA
   equals it with b_TT(5) of step 4;
   (b) factorization scale: the transversity kernel δP_qq = CF[2z/(1−z)₊
   + (3/2) δ(1−z)] on both legs gives c_D0 ∋ −2 CF v (s b/αs²) LF; the
   code's 2 v (−16 CA CF²(1−v))/CC = −4 CF² v(1−v)/CA equals it;
   (c) the ratio D0 P(v,1)/Born (the soft-gluon structure, spin independent
   for the single colour structure of class 5) of the mapped Fortran pieces
   equals exactly the unpolarized and longitudinal ratios of
   `Coefficients/ProtonProton.wl`, which were verified against the other two
   codes.

### The leg-exchanged ordering

For its hadronic convolution the code also builds each class with the two
incoming partons exchanged, a ↔ b, i.e. t ↔ u at fixed s and observed
momentum (`DPLUS`, slots FHODEL2/FHOREST2, for the classes of its lists
IA2/IA3, among them J0 = 5 and 11), from the same pieces:

    F02    = FBOR(1−v)
    FDEL2  = [AVDEL(1−v) + ln(v/(1−v)) AVWPL(1, 1−v) + ½ AVLO(1, 1−v) ln²((1−v)/v)] (1−v)/v
    FVWPL2 = [AVWPL(wx, vx) + AVLO(wx, vx) ln(v/vx)] vx/v
    FVLO2  = AVLO(wx, vx) vx/v
    FRESC2 = STRUV(wx, vx),       vx = 1 − v w,  wx = (1 − v)/(1 − v w),

assembled exactly as the first ordering (plus prescription on [0, 1] in w,
subtracted at w = 1, where vx = 1 − v and wx = 1; the same weights BXJAC and
XJAC).  Since FVWPL2 and FVLO2 depend on w, the regular part also holds their
remainders (FVWPL2(w) − FVWPL2(1))/(1−w) and (FVLO2(w) − FVLO2(1)) ln(1−w)/(1−w),
as `DPLUS` assembles them.  `PPTransverseFortranExchangedCoefficients[J0]` in
`pp-transverse-fortran.wl` (and in `transcribe_transpol.py`, which writes it)
returns this ordering in our basis, with the line numbers of the code's
functions.  The two orderings are different functions of (v, w): D1 is the same
for class 5, regular, D0 and δ differ.  Our classes observe the parton of the
code's direct ordering (`Processes/HadronicTT/4-Transverse.wl`).
