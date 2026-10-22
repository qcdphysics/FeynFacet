# Conventions

These conventions are common to all processes; each process README
restates the ones that matter for it.

## Kinematics and hard coefficients

Massless partons, D = 4 − 2ε dimensions.  For lepton scattering with the
incoming lepton l, the incoming parton k and the observed parton p,

    s = 2 l·k,   t = −2 k·p,   u = −2 l·p,   v = 1 + t/s,   w = −u/(s+t),   s4 = s + t + u = s v (1−w).

For proton–proton scattering the incoming partons are pa, pb and the
observed parton pc, with s = 2 pa·pb, t = −2 pa·pc, u = −2 pb·pc and the
same v, w.  The physical region is s > 0, 0 < v < 1, 0 < w < 1; Born
kinematics is w = 1.

The partonic hard coefficient is σ(s, v, w) = s E dσ/d³p.  The hadronic
cross section is

    E_h dσ/d³P_h = (1/S) Σ ∫ (dx/x) (dz/z²) f(x) D(z) σ(s, t, u)       (lepton–nucleon)

with one more ∫ dx/x for proton–proton collisions.  The measure dx/x
fixes the incoming-parton convolution.  With D = 4 − 2ε, the collinear
counterterm added to the bare hard coefficient in the notebooks is

    +(α_s/2π)(S_ε/ε)(μ_R²/μ_F²)^ε ∫ dz/z P(z) σ_Born(z s, z t, u),

where S_ε = (4π)^ε exp(−εγ_E).  This specifies the hard-coefficient
subtraction convention used in the finite assembly.

At NLO every result is written as

    σ = σ_Born δ(1−w) + (prefactor) [ R(v,w) + B(v) D0 + C(v) D1 + A(v) δ(1−w) ],
    D0 = [1/(1−w)]₊,   D1 = [ln(1−w)/(1−w)]₊,

the basis {regular, D0, D1, delta} of the package.  The plus prescription
is ∫₀¹ dw Dₙ(w) φ(w) = ∫₀¹ dw lnⁿ(1−w)/(1−w) [φ(w) − φ(1)].

## Spin

"Unpolarized" (U): average over the spins of the incoming particles, sum
over the final ones; density matrix p̸/2 for an incoming massless fermion.
"Polarized" (LL): the double-longitudinal helicity difference
Δσ = (σ₊₊ − σ₊₋)/2, obtained with the density matrix p̸γ₅/2 for each
incoming fermion; incoming gluons use the helicity density
i ε^{μνρσ} p_ρ n_σ/(2 p·n) with a light-like reference vector n.  Final
fragmentation is always unpolarized.

γ₅ is treated in the HVBM ('t Hooft–Veltman–Breitenlohner–Maison) scheme:
the external momenta of observed and incoming particles are four-dimensional,
the unobserved and loop momenta are D-dimensional and may have
(D−4)-dimensional ("hatted") components, which we denote by spe3 = r̂·r̂ and
similar symbols.  The trace convention is
Tr[γ^a γ^b γ^c γ^d γ₅] = −4i ε^{abcd} (FeynCalc with $LeviCivitaSign = −1);
Tracer's ε tensor has the opposite sign and the Tracer notebooks account for
it.  The HVBM scheme violates helicity conservation at the quark–gluon
vertex by a finite amount; the polarized quark density is therefore
renormalized with the additional finite kernel Δf_qq(z) = −4 CF (1−z)
(one subtraction per polarized incoming quark leg).

"Transverse" (Processes/Collins): the transverse spin transfer
T = σ(k↑ → p↑) − σ(k↑ → p↓) at fixed incoming transverse spin, from the
density k̸γ₅s̸ₖ/2 for the incoming quark and the projection p̸γ₅s̸ₚ for the
observed one, with the normal axis (T_NN) and the in-plane transverse axes of
the lepton–quark centre-of-mass frame (T_XX); the Collins asymmetry uses
T_C = (T_XX + T_NN)/2.  The spin vectors are four-dimensional; the component
of an unobserved momentum along the normal is reduced to invariants with the
Gram determinant of the physical momenta.  Both quark legs are subtracted
with the transversity kernel Δ_T P_qq(z) = CF [2z/(1−z)₊ + 3/2 δ(1−z)], the
photon in the lepton with the unpolarized kernel, and no helicity-restoring
term enters.

Proton–proton classes (Processes/Hadronic): a(pa) b(pb) → c(p) + X with
s = (pa+pb)², t = (pa−p)², u = (pb−p)², the 16 tagged channel classes and
real-emission weights of `Processes/Hadronic/classes.wl`; polarization sums
of external gluons with a light-like reference vector; closed quark loops
counted with Nf; a quark of another flavour (qp) stands for Nf − 1 flavours
when it is summed over.  The collinear counterterms of the two incoming legs
and of the fragmenting parton are built from the Born table of
Processes/PartonicBorn; their convolution formulas (legs a and b,
fragmentation) are derived in `4-Unpolarized.wl` and were checked by direct
numerical integration.

## Couplings, colour, scales

e² = 4π α, g_s² = 4π α_s; e_q is the quark charge in units of e.
Tr(T^a T^b) = δ^{ab}/2, CF = (Nc² − 1)/(2 Nc), CA = Nc, TR = 1/2; colour
averages 1/Nc per incoming quark and 1/(Nc² − 1) per incoming gluon.  The
coupling in D dimensions carries μ^{2ε}; the MSbar subtraction removes
poles together with ln 4π − γ_E.  muR is the renormalization scale, muF and
muFF the initial- and final-state factorization scales (set equal to muR
in the lepton-quark notebooks, as in the papers compared with; kept apart
in the Collins and proton–proton notebooks, with muL for the photon in the
lepton).  The coupling renormalization counterterm of an order-α_s² Born
term is −(α_s/2π)(β₀/ε) S_ε × Born with β₀ = (11 CA − 2 Nf)/3.  Splitting
functions in units of α_s/(2π) are listed in `FeynFacet/Factorization.wl`.

## Integrals

Real emission: the unobserved two-particle phase space is written with the
cut propagators of reverse unitarity, so that Kira can reduce the integrals;
each process notebook defines its families.  The two-body phase-space
measure is d^{D}r/(2π)^{D−1} δ₊(r²) 2π δ₊((pRec − r)²), giving the volume
(4π)^ε s4^{−ε}/(8π) Γ(1−ε)/Γ(2−2ε).  A hatted invariant spe3ⁿ under the
integral is replaced by (4π)ⁿ (−1)ⁿ (−ε)ₙ times the integral in D + 2n
dimensions.  Loop integrals use the measure d^D loop/(i π^{D/2}); with it a
power (l̂·l̂)ⁿ of the hatted loop momentum is replaced by (−1)ⁿ (−ε)ₙ times the
integral in D + 2n dimensions (for n = 1: ε I_{D+2}).  In both cases only the
dimension of the integral changes: the coefficients produced by the trace keep
D = 4 − 2ε.  Shifting their D as well is a tempting mistake (the reduction
rules carry a dimension symbol of their own) and changes the δ(1−w)
coefficient of the polarized pp classes by a Born-proportional constant.
Master integrals are one-dimensional Euler integrals evaluated with
SubTropica (HyperIntica backend) and checked against closed forms.

## TMD operators (Processes/TMD)

Light-cone vectors n, nb with n² = nb² = 0 and n.nb = 2; the incoming parton
p is along nb, so n.p is its large light-cone component and γ⁺ = n-slash.
Beam functions: B_q(x, b) = ∫ db⁻/(4π) e^{−i x n.p b⁻} ⟨p| q̄(b) W_n†(b) n-slash W_n(0) q(0) |p⟩
and B_g = ∫ db⁻/(2π x n.p) e^{…} ⟨p| F^{+μ a}(b) W^A… F^{+}_μ{}^{a}(0) |p⟩ (sum over the
transverse components, i.e. contraction with −g_{μν}), both equal to
δ(1−x) at lowest order; the Wilson lines run from the operator point to +∞
along n.  The soft function is the vacuum expectation value of the trace of
the Wilson lines along n and nb, divided by the dimension of the
representation, equal to 1 at lowest order.  Squared operator matrix
elements are normalized by 1/(2 n.p) (quark) and 1/(x (n.p)²) (gluon), which
are the factors of the b⁻ Fourier transforms.

Real emission: the emitted parton k has n.k = (1−x) n.p and transverse
momentum kT2 = (n.k)(nb.k), and p.k = kT2/(2(1−x)).  External gluons are
summed with the physical polarization sum −g^{μν} + (k^μ n^ν + n^μ k^ν)/(n.k).
Rapidity regulator: the real phase space is weighted by
exp[−τ e^{−γ_E}(n.k + nb.k)] with ν = 1/τ (Li, Neill, Zhu), the limit τ → 0
being taken before ε → 0; the zero-bin subtraction divides by the bare soft
function.  Fourier transform to b space: B(x, b) = ∫ d^{D−2}kT e^{i b.kT} B(x, kT).
Results in as = αs/(4π) with Lp = ln(b² μ²/b0²), b0 = 2 e^{−γ_E},
LQ = 2 ln(x n.p/ν), Lν = ln(ν²/μ²), Lh = LQ + Lν; the TMD PDF coefficient is
f_ij = I_ij √S_i.  These are the conventions of Luo, Yang, Zhu and Zhu,
arXiv:2012.03256, whose ancillary tables are the reference.
