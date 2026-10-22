# Lepton–quark scattering with transverse spin transfer: ℓ q↑ → q↑ + X (Collins)

The NLO hard part for single-inclusive quark production in lepton–quark
scattering with a transversely polarized incoming quark and the transverse
spin of the observed quark measured, the quantity that feeds the Collins
asymmetry of single-inclusive jets at the EIC (Section "Application to EIC
phenomenology: the Collins effect" of the paper).  The diagrams are those of
`Processes/LeptonQuark`; the new ingredients are the transverse spin
projections, the transversity splitting kernel on both quark legs, and the
three separate factorization scales.  The notebooks end by comparing with the
result written in the paper (Eqs. transverse-finite-hard, transverse-endpoints
and the appendix coefficient), and the resolved-photon Born table with the
paper's Eq. (resolved-transverse-born).

## The observable

T = σ(k↑ → p↑) − σ(k↑ → p↓) at fixed incoming transverse spin, obtained from
the spin average by inserting the density k̸ γ₅ s̸ₖ/2 for the incoming quark and
the projection p̸ γ₅ s̸ₚ for the observed one.  Two axes are kept: the normal n
to the scattering plane (T_NN) and the in-plane transverse directions of the
lepton–quark centre-of-mass frame (T_XX); the Collins asymmetry needs the
projection T_C = (T_XX + T_NN)/2, whose Born value is −4su/t² (times e⁴e_q²)
with no evanescent term.  The HVBM γ₅ is used; the component of the
unobserved gluon momentum along the normal is reduced to the invariants with
the Gram determinant of l, k, p.

## Files

| File | What it does | Kernel time |
|---|---|---|
| `1-Amplitudes.wl` | FeynArts diagrams and amplitudes (as in the lepton–quark process) | 10 s |
| `2-Algebra-FeynCalc.wl` | Squares with FeynCalc: spin average, T_NN, T_XX and their Collins projection (Born, real, photon Born), virtual form factor | 1 min |
| `2-Algebra-FORM.wl` | The same with FORM; the programs are in `form/` | 13 min (the XX real-emission trace) |
| `2-Algebra-Tracer.wl` | The same with Tracer | 1 min |
| `3-Reduce.wl` | Reverse unitarity, Kira reduction (`kira/`, numerator indices allowed), master integrals | 1 min |
| `4-Unpolarized.wl` | The unpolarized coefficient re-derived from this process's files and compared with arXiv:1505.06415 (consistency check of the pipeline) | 2 min |
| `5-Transverse.wl` | Threshold expansion, virtual term, transversity counterterms with scales μ_I, μ_O, μ_L, scale-derivative checks, finite coefficient, comparison with the paper | 3 min |
| `6-ResolvedBorn.wl` | Transverse spin transfer at Born level in parton–parton scattering (the resolved-photon channels), compared with the paper's table | 2 min |
| `form/`, `kira/`, `reference/` | FORM programs, Kira families (those of the lepton–quark process), the paper's result and arXiv:1505.06415 transcribed for the comparisons |

## How to run

As for the lepton–quark process: steps 1, 2 (one route), 3, 4, 5 in fresh
kernels; step 6 is independent.  Headless:

```sh
python Checks/run_process.py Processes/Collins --route FeynCalc
python Checks/run_process.py Processes/Collins --route FORM
python Checks/run_process.py Processes/Collins --route Tracer
```

## What is checked, and where

| Step | Check |
|---|---|
| 2 (all routes) | unpolarized Born square equals the textbook form; T_NN and T_XX at Born level equal −4su/t² ∓ (D−4); the in-plane axes have the stated orthogonality and lengths; the two transverse components of the real emission coincide in four dimensions; the photon Born transverse tensor equals 4 − 2ε²t²/((1−ε)su); the explicit-trace construction reproduces FeynCalc's generic spin average (FeynCalc route); the three routes agree |
| 3 | every term rebuilt from its family rows; the family terms add up to the matrix element; integrals within the Kira ranges; master integrals as in the lepton–quark process |
| 4 | the unpolarized coefficient equals arXiv:1505.06415 |
| 5 | Born hard part equals the paper's; double poles cancel; single poles cancel between real, virtual, the two transversity counterterms and the photon-in-lepton subtraction (table displayed); no μ_R dependence; the derivatives with respect to μ_I, μ_O and μ_L equal the corresponding splitting convolutions; the finite coefficient equals the paper's, symbolically and at three numerical points |
| 6 | the two transverse components coincide in four dimensions; the Collins projections equal the paper's table; control channels (transfer to a different flavour) vanish; charge-conjugate channels equal their representatives |

## Conventions in brief

σ(s, v, w) = s E dσ/d³p, t = −s(1−v), u = −s v w, s4 = s v (1−w).  The result
is written as 4α²e_q²/(s(1−v)²) δ(1−w) + α²e_q²C_F α_s/(π s v w)[R_T + B_T D0 +
C_T D1 + A_T δ(1−w)] with L_I = ln(μ_I²/s), L_O = ln(μ_O²/s), L_L = ln(μ_L²/s)
for the incoming-quark, fragmentation and lepton–photon scales; the
transversity kernel is C_F[2z/(1−z)₊ + 3/2 δ(1−z)] on both legs, the photon
in the lepton is unpolarized, and no helicity-restoring term enters.
