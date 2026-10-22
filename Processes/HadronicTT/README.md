# Transversely polarized proton–proton collisions: p↑ p↑ → h + X (A_TT) at NLO

The NLO QCD corrections to single-inclusive parton production in collisions
of two transversely polarized protons, a(pa) b(pb) → c(p) + X with both
incoming partons transversely polarized and the observed parton unpolarized:
the hard part of the double-spin asymmetry A_TT = dδσ/dσ,
dδσ = ½[dσ(↑↑) − dσ(↑↓)].  Five of the sixteen channel classes of
`Processes/Hadronic` contribute (`classes.wl` explains why the others vanish):
5 (q q̄ → q′), 6 (q q → q), 7 (q q → g, NLO only), 11 (q q̄ → q) and
12 (q q̄ → g).  The diagrams, the squaring machinery, the reduction and the
counterterm formulas are those of `Processes/Hadronic`; the transverse spin
projections and the normal-component reduction are those of
`Processes/Collins` (docs/TRANSVERSE-PP.md states the shared conventions).
The results are compared with the Fortran code of Mukherjee, Stratmann and
Vogelsang, Phys. Rev. D 72 (2005) 034011 [hep-ph/0506315], transcribed
verbatim in `Reference/pp-transverse-fortran.wl` (map onto our normalization
in `Reference/fortran/README.md`), and the Born table with the lowest-order
cross sections of Soffer, Stratmann and Vogelsang, Phys. Rev. D 65 (2002)
114024 [hep-ph/0204058], transcribed in `reference/SSV-born.wl`.

## The observable

Each incoming quark or antiquark carries the chiral-odd density
p̸ γ5 s̸/2 (HVBM γ5, four-dimensional spin vectors; the same expression for an
antiquark), whose partonic cross section multiplies the transversity
δf = f↑↑ − f↓↑ of its proton.  The spin vectors are written in the basis of
the scattering plane, s = cos φ e_X + sin φ n, with e_X the in-plane axis
along the transverse momentum of the observed parton and n the normal, and
the squared matrix element is computed as the four components XX, NN, XN, NX
(coefficients of cos φ_a cos φ_b, sin sin, cos sin, sin cos).  Parity demands
XN = NX = 0, and the absence of a Φ-independent term (two helicity flips
carry two units of angular momentum along the beam) demands XX + NN = 0 in
four dimensions; both are checked for the Born, the virtual and the real
emission of every class.  In the HVBM scheme XX + NN does not vanish in D
dimensions (an evanescent term ∝ (D−4)(s_a·s_b) from the hatted γ matrices
between the two densities), so the hard part of A_TT is defined as the
coefficient of cos(2Φ),

    δσ̂_TT = (XX − NN)/2,

which is the covariant projector F = (s/(tu))[2(p·s_a)(p·s_b) + (tu/s) s_a·s_b]
of the paper and equals XX in four dimensions.  The evanescent symmetric
part (XX + NN)/2 is carried through the reduction and the assembly as a
check: its poles cancel by themselves (its Born term is of order ε) and its
finite remainder vanishes for every class, so the HVBM scheme leaves no
Φ-independent artifact and no finite renormalization is needed (unlike the
helicity-restoring term of the longitudinal case).

## Files

| File | What it does |
|---|---|
| `classes.wl` | the five classes (entries as in `Processes/Hadronic/classes.wl`) and why the other eleven vanish |
| `1-Amplitudes.wl` | FeynArts diagrams and amplitudes of one class (as in `Processes/Hadronic`) |
| `2-Algebra-FeynCalc.wl` | squares with FeynCalc (Hadronic's parallel, colour-separated, i ≤ j, grouped-projection machinery) with the two transverse densities, for the four components XX, NN, XN, NX: Born, one-loop interference (scalar-integral families with the normal component of the loop momentum reduced through the Gram matrix, (n·l)² = l∥² − l² + l̂²), real emission (the normal component of k1 reduced as in `Processes/Collins`: odd powers dropped, (n·k1)² = k1∥² + spe3); the checks XN = NX = 0 and XX + NN = 0 at D = 4; gauge invariance of the real emission; the D-dimensional TT Born table of all channels needed by the counterterms, checked against the class's own Born, under a ↔ b relabelling, charge conjugation, vanishing control channels (q q′ → q, q q̄′ → q) and the literature at D = 4, N_c = 3 |
| `2-Algebra-FORM.wl` | the same step with FORM doing the spin and colour algebra (`form/common.h`, `form/square.frm`, `form/checks.frm`: the two transverse densities with the HVBM γ5 written as epsilon tensors, the axes as dot-product substitutions, the colour library of FORM; one FORM job per diagram pair and component); the normal components are reduced in Mathematica as in the FeynCalc route; the same Born table and checks; ends with the exact comparison with the file of the FeynCalc route |
| `2-Algebra-Tracer.wl` | the same step with the classic Tracer package for the Dirac algebra and the explicit SU(N) colour reduction of `Processes/Hadronic/2-Algebra-Tracer.wl`; the transverse densities enter as the strings `q, G5, s` when the fermion lines are closed, with the axes aX and nn as Tracer vectors; the Born table is squared from hand-written tree records; the same checks, and at the end the exact comparison of Born, Virtual, Real (TT and Symmetric) and the Born table with the FeynCalc file |
| `3-Reduce.wl` | `Processes/Hadronic/3-Reduce.wl` for the projections TT and Symmetric: partial fractions, Kira (`kira/`, Hadronic's families), master integrals, continuation |
| `4-Transverse.wl` | threshold expansion, virtual normalization 2/(4π)^{D/2} with phase +1, derived from the FeynArts loop measure and the master-integral measure (checked by the universal double pole −(α_s/2π) Σ C_i × Born and, separately, by the cancellation against the real emission), MSbar coupling renormalization, collinear counterterms (transversity kernel on both incoming legs with the TT Born; unpolarized kernels for the fragmenting parton with the TT Born table; no helicity-restoring term), pole cancellation for TT and for the symmetric projection, scale-derivative checks (μ_F, μ_FF, μ_R), finite coefficient, comparison with the Fortran transcription |
| `reference/SSV-born.wl` | Table I and Eq. (6) of Soffer, Stratmann and Vogelsang, transcribed |
| `kira/` | the Kira families of `Processes/Hadronic` (copied, unchanged) |

## How to run

Each step file begins with `class = "5"`; set the class (5, 6, 7, 11, 12) and
run steps 1, 2, 3, 4 in fresh kernels.  Headless, one class at a time:

```sh
python Checks/run_process.py Processes/HadronicTT --set class=5
FEYNFACET_SET="route=FeynCalc;class=5;workers=8" wolframscript -file Checks/run_step.wls Processes/HadronicTT/2-Algebra-FeynCalc.wl
```

Step 2 runs the products in parallel (`workers = "auto"`: one per core, at
most 8).  Three routes are implemented for step 2: FeynCalc
(`2-Algebra-FeynCalc.wl`), FORM (`2-Algebra-FORM.wl`) and Tracer
(`2-Algebra-Tracer.wl`).  The override `route=FORM` or `route=Tracer` selects
the file of step 2, and steps 3 and 4 run with the same override read its
output; the FORM and Tracer routes end with a comparison with the FeynCalc result
(exact rational arithmetic at two random rational points, see "Run-time
design"), which is skipped, with a message, when the FeynCalc result of the
class is missing: run the FeynCalc route of a class first.  One class
through all steps of one route:

```sh
for step in 1-Amplitudes 2-Algebra-FORM 3-Reduce 4-Transverse; do
  FEYNFACET_SET="route=FORM;class=5" wolframscript -file Checks/run_step.wls Processes/HadronicTT/$step.wl
done
for step in 1-Amplitudes 2-Algebra-Tracer 3-Reduce 4-Transverse; do
  FEYNFACET_SET="route=Tracer;class=5" wolframscript -file Checks/run_step.wls Processes/HadronicTT/$step.wl
done
```

(`1-Amplitudes.wl` is the same for all routes and needs to run only once per
class.)  Results go to `Results/<route>/class-<n>/` (`Algebra.wl`,
`Reduced.wl`, `Transverse.wl`, the product files in `pairs/`).  Class 12
(q q̄ → g) is the largest.  Step 2 uses parallel algebra kernels; reduce
`workers` when memory is limited.  Steps 3 and 4 use one Wolfram kernel, with
external Kira processes during reduction.  For
class 12, sampled process-tree resident memory reaches approximately
10.9 GiB during finite assembly (step 4). This periodic measurement is a
lower bound on the peak, not a sufficient-memory requirement. Reducing the
algebra worker count does not limit the memory of scalar reduction or finite
assembly. Measure these stages separately; see
[resource guidance](../../docs/RESOURCES.md).

Reduction design. Steps 3 and 4 use the same mathematical implementation
for each route, acting on that route's own `Algebra.wl`. FORM and Tracer
compare their algebra with FeynCalc at two exact-rational points, assigning
every scalar integral its own value. This is a sampled comparison.
All five classes are reduced and assembled separately in all three routes.
The collector compares their finite coefficients at six sampled points.

## Results

The rows summarize the comparisons performed by steps 2–4 in the three routes.  "Squares" is
the comparison at the end of step 2 (exact rational arithmetic at two random
rational points) (Born, virtual and real-emission squares for
the projections TT and Symmetric, and the Born table; FORM and Tracer against
FeynCalc).  "Fortran" is the comparison of step 4 of every route with the
transcribed reference code (`Reference/pp-transverse-fortran.wl`, mapped as
in `Reference/fortran/README.md`): the four finite coefficients
{regular, D0, D1, δ} and the Born term, at six rational (v, w) points,
N_f = 4, 5 and distinct fixed scales, evaluated with 30-digit arithmetic.
The required SU(3) residual is below 10⁻²⁰. Additional C_A = 2, 5 comparisons,
with C_F = (C_A²−1)/(2C_A), are reported separately; the reference code's
general-colour validity is not assumed. "Symmetric" is the finite
remainder of the evanescent projection (XX + NN)/2 after renormalization and
factorization.  Poles cancel and the three scale-derivative checks pass in
every row.

| Class | Process | Born | Squares FeynCalc = FORM = Tracer | Finite coefficients vs Fortran (SU(3)) | Symmetric remainder | Status |
|---|---|---|---|---|---|---|
| 5 | q q̄ → q′ | yes | agree (2 points) | agree (all four, each route) | 0 | complete |
| 6 | q q → q | yes | agree (2 points) | agree (all four, each route) | 0 | complete |
| 7 | q q → g | no | agree (2 points) | agree (all four, each route) | 0 | complete |
| 11 | q q̄ → q | yes | agree (2 points) | agree (all four, each route) | 0 | complete |
| 12 | q q̄ → g | yes | agree (2 points) | agree (all four, each route) | 0 | complete |

The coefficients of the three routes agree for every class, at
the six random points (N_c = 3 and 5) of the numerical comparison of
`Checks/collect_coefficients.wls` (they are collected, with this check, in
`Coefficients/ProtonProtonTT.wl`).  The Born terms equal the lowest-order
cross sections of Soffer, Stratmann and Vogelsang (`reference/SSV-born.wl`).

HVBM note.  The Born squares of step 2 have an evanescent symmetric part:
(XX + NN)/2 is proportional to D − 4 (class 5: −(D−4) C_F g⁴/(2C_A); classes
6 and 11 likewise; class 12: zero); the
hard part is TT = (XX − NN)/2, the coefficient of cos(2Φ), and the finite
remainder of the symmetric projection (XX + NN)/2 vanishes in every completed
class (column "Symmetric").

The crossing of these results into the transverse spin transfer
p↑ p → h↑ + X is `Processes/HadronicTransfer/5-Crossing.wl`.

## Conventions in brief

a(pa) b(pb) → c(p) + X; s = (pa+pb)², t = (pa−p)², u = (pb−p)², v = 1 + t/s,
w = −u/(s+t), s4 = s v (1−w).  δσ̂_TT = s E dδσ/d³p projected on cos(2Φ), Φ
the azimuth of p relative to the common spin axis of the two protons;
δσ̂_TT = α_s² b δ(1−w) + α_s³/(π s v w) [R + B D0 + C D1 + A δ(1−w)] in the
basis {regular, D0, D1, delta} of the package, with the scales μ_R
(coupling), μ_F (both transversity densities), μ_FF (fragmentation).  The
transversity kernel is C_F[2z/(1−z)₊ + (3/2) δ(1−z)] on both incoming legs
(no gluon mixing); the fragmenting parton is subtracted with the unpolarized
kernels; MSbar throughout.  The Fortran code's slots hold one half of the
cos(2Φ) coefficient: Born = 2(1−v) F0 α_s²/π and
{R, B, C, A} = 2 α_s³/(π s CC) {STRUV, AVWPL(1), AVLO(1), AVDEL}
(`Reference/fortran/README.md`).
