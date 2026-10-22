# Transversely polarized proton–proton collisions: conventions

Two observables with transversity, both built from the pieces of
`Processes/Hadronic` (the 16 channel classes, the squaring, the reduction, the
counterterms, the comparison with the reference codes) and of
`Processes/Collins` (transverse spin densities, in-plane and normal axes, the
reduction of the normal component of an unobserved momentum):

| Process folder | Observable | Reference |
|---|---|---|
| `Processes/HadronicTT` | p↑ p↑ → h + X: both incoming partons transversely polarized, the hadron unpolarized (the asymmetry A_TT) | Mukherjee, Stratmann, Vogelsang, Phys. Rev. D 72 (2005) 034011 [hep-ph/0506315]; the Fortran code `me-transpol-ms.f` of the authors, transcribed in `Reference/pp-transverse-fortran.wl` |
| `Processes/HadronicTransfer` | p↑ p → h↑ + X: one incoming parton transversely polarized, the hadron transversely polarized (transverse spin transfer D_TT) | new NLO results; checked against the double-transverse process by generalized crossing symmetry |

## Kinematics (as in Processes/Hadronic)

a(p_a) b(p_b) → c(p) + X, with s = (p_a+p_b)², t = (p_a−p)², u = (p_b−p)²,
v = 1 + t/s, w = −u/(s+t), s4 = s + t + u = s v (1−w).  The partonic cross
section is σ = s E dσ/d³p, written as

    σ = Born δ(1−w) + P(v,w) c_regular + P(v,1) [ c_D0 D0 + c_D1 D1 + c_δ δ(1−w) ],
    P = α_s³/(π s v w),  D0 = [1/(1−w)]_+,  D1 = [ln(1−w)/(1−w)]_+ ,

the same basis {regular, D0, D1, delta} and the same scale logarithms
LR = ln(μ_R²/s), LF = ln(μ_F²/s) (incoming partons), LFF = ln(μ_FF²/s)
(fragmentation) as `Coefficients/ProtonProton.wl`.

## Spin densities (HVBM γ5, four-dimensional spin vectors)

- Incoming transversely polarized quark a: u ū → p̸_a γ5 s̸_a / 2, with
  s_a·p_a = 0, s_a² = −1.  This is the density whose partonic cross section
  multiplies the transversity distribution δf_a = f_{↑↑} − f_{↓↑} in the
  hadronic cross section, i.e. the partonic quantity is
  ½[σ̂(s_a) − σ̂(−s_a)].  An incoming antiquark: v v̄ → p̸_a γ5 s̸_a / 2 as well
  (the sign is the same for transversity; HVBM: γ5 anticommutes with the
  four-dimensional γ's only).
- Observed transversely polarized quark c (HadronicTransfer): the projection
  u ū → p̸ γ5 s̸_c (no ½): the partonic quantity is σ̂(c: s_c) − σ̂(c: −s_c) at
  fixed s_a, as in Processes/Collins.
- Unpolarized legs as in Processes/Hadronic (p̸/2 incoming quark, physical
  polarization sum with 1/(D−2) for an incoming gluon, p̸ for an outgoing quark).

## Axes

All spin vectors are written in the basis of the scattering plane of p_a, p_b, p:

- in-plane axis for the incoming partons (transverse to both beams, along the
  transverse momentum of p):  e_X = sqrt(s/(t u)) [ p + (u/s) p_a + (t/s) p_b ],
  e_X·p_a = e_X·p_b = 0, e_X² = −1, e_X·p = −sqrt(tu/s) = −p_T
  (in Processes/Collins this is `aX` with l → p_b, k → p_a: aX = (t/s) l + (u/s) k + p
  up to the normalization);
- in-plane axis for the observed parton (HadronicTransfer), transverse to p
  and to p_a + p_b: Processes/Collins' `bX` with the relabeling l → p_b, k → p_a;
- the normal n ∝ ε^{μ p_a p_b p}, n² = −1, orthogonal to p_a, p_b, p; its
  products with an unobserved momentum k1 (and k2 = p_a+p_b−p−k1 = −k1 on n)
  occur only in even powers; (n·k1)² = spe3 + k1∥², with k1∥² the square of the
  projection of k1 onto the plane (Gram matrix of p_a, p_b, p; Processes/Collins,
  `normalReduction`).

The spin vectors are s_a = cos φ_a e_X + sin φ_a n, s_b = cos φ_b e_X + sin φ_b n
(HadronicTT) and s_c = cos φ_c e_C + sin φ_c n (HadronicTransfer).  The
squared matrix elements are computed as the four components XX, NN, XN, NX
(the coefficients of cos cos, sin sin, cos sin, sin cos); the mixed
components must vanish (parity).

## HadronicTT

The transversely polarized cross section has the azimuthal dependence
cos(2Φ), Φ the azimuth of p relative to the common spin axis of the two
protons (both along x̂).  With both spins along e_X one has Φ = −φ_a = −φ_b,
so the cross section is T_XX cos φ_a cos φ_b + T_NN sin φ_a sin φ_b
= ½(T_XX − T_NN) cos(φ_a + φ_b) + ½(T_XX + T_NN) cos(φ_a − φ_b), and the
absence of a Φ-independent term (the paper's statement that the Φ-integral
vanishes) means T_NN = −T_XX in four dimensions.  The hard part is the
coefficient of cos(2Φ):

    δσ̂_TT = (T_XX − T_NN)/2   (the paper's projector
                    F(p, s_a, s_b) = (s/(tu)) [2 (p·s_a)(p·s_b) + (tu/s) s_a·s_b],
                    which equals +1 for s_a = s_b = e_X and −1 for s_a = s_b = n).

In the HVBM scheme the symmetric part (T_XX + T_NN)/2 does not vanish in D
dimensions: already the Born squares carry an evanescent term proportional
to D − 4 (from the hatted γ matrices between the two densities).  It is
carried through the reduction and the assembly as the projection
`Symmetric`: its poles cancel by themselves and its finite remainder after
renormalization and factorization vanishes in every class (checked in step 4
of `Processes/HadronicTT`), so the scheme leaves no Φ-independent artifact
and needs no finite renormalization.

Channel classes (numbering of Processes/Hadronic/classes.wl and of the
reference codes).  Both polarized partons must be quarks or antiquarks on
different fermion lines connected by annihilation or by an exchange graph; a
t-channel-only process (q q′ → q X, q q̄′ → q X) and every gluon-initiated
process vanish identically (helicity conservation along a massless line):

| Class | Process | Born |
|---|---|---|
| 5 | q q̄ → q′ + X | yes |
| 6 | q q → q + X | yes |
| 7 | q q → g + X | no (NLO only) |
| 11 | q q̄ → q + X | yes |
| 12 | q q̄ → g + X | yes |

Counterterms: transversity kernel δP_qq(z) = C_F [2z/(1−z)_+ + (3/2) δ(1−z)]
(`SplittingKernel["Transverse qq", z]`) on both incoming legs (no gluon
mixing), the unpolarized kernels for the fragmenting parton, the Born table of
the TT cross sections in D dimensions (its O(ε) terms kept), MSbar coupling
renormalization.  No helicity-restoring finite term (that term belongs to the
longitudinal case only).

The Born terms and four finite coefficients of all five classes are
checked against the reference code at sampled SU(3) points. Comparisons
for C_A = 2, 5 are reported separately; the reference code's general-colour
validity is not assumed. All five classes are reduced and assembled separately in each of the three
algebra routes. The checks are described in
`Processes/HadronicTT/README.md`, and the coefficients are in
`Coefficients/ProtonProtonTT.wl`.

Reduction design (both processes). Steps 3 and 4 share their mathematical
implementation and act on each route's own `Algebra.wl`. The step-2 comparison
uses exact rational arithmetic at two random points, assigning every scalar
integral its own value: this is a sampled comparison, not a symbolic identity.
Each route supplies its own finite coefficients, which the collector compares
at six sampled kinematic points.

## HadronicTransfer

δσ̂_T(a↑ b → c↑) is the cross section with the densities above; its
components T_XX (cos φ_a cos φ_c) and T_NN (sin φ_a sin φ_c) are both kept.
In D dimensions they differ by evanescent terms, but after renormalization
and factorization the finite coefficients of XX and NN are equal in every
class (checked in step 4 and by the crossing of step 5): the spin transfer
depends only on cos(φ_a − φ_c), as for the Born term at D = 4.

Channel classes (a is the polarized incoming parton, c the observed polarized
parton; a and c must be on one fermion line, or connected through an
annihilation–exchange interference):

| Class | Process | Born | Obtained from HadronicTT by crossing |
|---|---|---|---|
| 1 | q q′ → q + X | yes | class 5 (observing q̄′) |
| 3 | q q̄′ → q + X | yes | class 5 (observing q′); its Born term equals that of class 1, the NLO terms differ (as for the unpolarized classes 1 and 3) |
| 6 | q q → q + X | yes | class 11 with the observed antiquark |
| 10 | q g → q̄ + X | no (NLO only) | class 7 |
| 11 | q q̄ → q + X | yes | class 11 |
| 13 | q g → q + X | yes | class 12 |
| 17 | q̄ q → q + X (the antiquark polarized, the quark observed; = q q̄ → q̄ by charge conjugation) | yes | class 6 |

Counterterms: transversity kernel on leg a and on the fragmenting quark c
(δP_qq, no mixing), unpolarized kernels on leg b (with the lower Born of the
same transfer type), the D-dimensional Born table of the transfer cross
sections.

All seven classes are reduced and assembled separately in the three algebra
routes, with finite coefficients compared at six sampled points.
The calculation is described in `Processes/HadronicTransfer/README.md`;
the coefficient tables are in `Coefficients/ProtonProtonTransfer.wl`.

## Generalized crossing symmetry

Knowing a⃗ b⃗ → c + X, the process a⃗ c̄ → b̄⃗ + X follows by s ↔ t:
v → −v/(1−v), w unchanged, with attention to log(v) → log(−v/(1−v)), the π²
terms (the virtual part must be crossed before its real part is taken), the
spin averaging factors, and the interchange of the incoming-parton and
fragmentation counterterms.  `Processes/HadronicTransfer/5-Crossing.wl` checks
it at LO and NLO against the results of `Processes/HadronicTT`.

The map is b ↔ c̄: p_b → −p, p → −p_b, the unobserved momenta unchanged
(s ↔ t, u unchanged; for the pairs that observe the antiquark of the TT
process, 1 ← 5 and 6 ← 11, the composed map s → t, t → u, u → s).  With the
conventions above the components cross as

    XX_transfer = −σ g XX_TT(crossed),   NN_transfer = +σ g NN_TT(crossed),

σ = +1, g = 1 when the parton b of the transfer process is a quark or
antiquark, σ = −1, g = 2C_A/((D−2)(C_A²−1)) when it is a gluon (g from the
colour average 1/(C_A²−1) and the polarization average 1/(D−2) of the
incoming gluon against the quark averages of the TT process; σ because only
one fermion line is crossed; `Processes/HadronicTransfer/README.md` derives
both).  Hence the TT projection (XX−NN)/2 crosses
into the transfer's (XX+NN)/2 and the evanescent Symmetric projection into
(XX−NN)/2; the vanishing finite Symmetric remainder of the TT classes is the
statement XX = NN for the finite transfer coefficients.  At LO and D = 4,
b_transfer(v) = σ g b_TT(−v/(1−v))/(1−v).

Outcome, for all seven pairs: Part A holds (Born and Born table
identically, the virtual identically after reduction, the real-emission
squares at three random exact rational points), and so does Part B (the
crossed TT algebra run through steps 3 and 4 of the transfer process: poles
cancel, and the finite coefficients equal the direct ones in exact arithmetic
at three rational points of s, v, w, N_f with symbolic scales).

## The reference code

`me-transpol-ms.f` (Stratmann, 2004) has the structure of the unpolarized and
longitudinal codes (`Reference/fortran/README.md`): per class J0 the pieces
FBOR (Born), AVDEL (δ(1−w), with ln(μ_R²/s) inside), AVWPL (1/(1−w)),
AVLO (ln(1−w)/(1−w)), STRUV (regular), the colour normalization CC(J0), and
the hadronic assembly of `DPLUS`, which multiplies by 4/π for the cos(2Φ)
projection over four quadrants.  Its coefficients are those of cos(2Φ) in
s dδσ̂/(dv dw dΦ) (the paper's Eq. (10)); the normalization map to our
σ = s E dσ/d³p is fixed on the Born term of class 5 and must then hold for
every class and coefficient.  It is (`Reference/fortran/README.md` derives
it): Born = 2(1−v) F0 α_s²/π and {regular, D0, D1, δ} = 2 α_s³/(π s CC)
{STRUV, AVWPL(1), AVLO(1), AVDEL}, CC = 8 N_c², i.e. the code's slots hold
one half of the cos(2Φ) coefficient.  For the classes of its lists IA2 and IA3
(among the transverse classes, J0 = 5 and 11) `DPLUS` also builds the ordering with the two
incoming partons exchanged, (v, w) → (1 − vw, (1 − v)/(1 − vw)) in the
regular part with the distributions re-expanded (FDEL2, FVWPL2, FVLO2, FRESC2,
FBOR(1−v)); `PPTransverseFortranExchangedCoefficients[J0]` of
`Reference/pp-transverse-fortran.wl` maps it onto our basis.  Reference comparisons evaluate the analytic coefficient routines directly.
