# Transverse spin transfer in proton–proton collisions: p↑ p → h↑ + X at NLO

The NLO QCD corrections to single-inclusive parton production in proton–proton
collisions with the incoming parton a transversely polarized and the transverse
spin of the observed parton c measured, a↑ b → c↑ + X, the partonic hard parts
of the transverse spin transfer D_TT (docs/TRANSVERSE-PP.md).  The diagrams,
the squaring machinery, the integrals and the counterterm formulas are those of
`Processes/Hadronic`; the new ingredients are the transverse spin projections
of `Processes/Collins` (the density p̸_a γ5 s̸_a/2 for a, the projection
p̸ γ5 s̸_c for c, the in-plane axes e_X, e_C and the normal n), the reduction of
the normal components of the unobserved and loop momenta, the transversity
kernel on both polarized legs with a D-dimensional transfer Born table, and the
seven channel classes in which transversity can be transferred (`classes.wl`
explains why the other classes vanish and lists the crossing partner of each
class in `Processes/HadronicTT`).  There is no published NLO result; the
checks are the internal ones (pole cancellation, scale derivatives, gauge
invariance, the vanishing of the mixed spin components, agreement
of the three algebra routes at sampled exact-rational points), the generalized crossing from the
double-transverse results of `Processes/HadronicTT`, which agree with the
reference code of Mukherjee, Stratmann and Vogelsang (step 5).

## The observable

σ(s_a, s_c) = s E dσ/d³p with the spin vectors s_a = cos φ_a e_X + sin φ_a n
and s_c = cos φ_c e_C + sin φ_c n, where e_X = (u/s) pa + (t/s) pb + p
(normalized; transverse to both beams, along the transverse momentum of p),
e_C = −u/(t+u) pa + t/(t+u) pb − s(u−t)/(t+u)² p (normalized; transverse to p
and to pa + pb) and n ∝ ε^{μ pa pb p}.  The squared matrix elements are the four
components XX, NN, XN, NX (the coefficients of cos φ_a cos φ_c, sin sin, cos sin,
sin cos); XN and NX vanish by parity (checked at every stage), XX and NN are
both carried through to the finite coefficients and saved.  In D dimensions
they differ by evanescent terms; the finite NLO coefficients of XX and NN are
equal in every class (to 30 digits), so the spin transfer depends only on
cos(φ_a − φ_c).  The transverse spin transfer D_TT of a hadron
polarized along the in-plane direction of the parent proton uses XX.

## Files

| File | What it does |
|---|---|
| `classes.wl` | the seven classes 1, 3, 6, 10, 11, 13, 17 with their Born, virtual and real-emission processes and weights (the numbering of `Processes/Hadronic`; 17 = q̄ q → q with the antiquark polarized), why the others vanish, the crossing partners; the list of 2 → 2 transfer channels of the Born table |
| `1-Amplitudes.wl` | FeynArts diagrams and amplitudes of one class (as `Processes/Hadronic`), plus the Born amplitudes of the transfer channels of the table |
| `2-Algebra-FeynCalc.wl` | the squaring machinery of `Processes/Hadronic` (parallel kernels, colour separated and cached, products i ≤ j, grouped polarization sums) with the transverse projections of both polarized legs; the D-dimensional transfer Born table (checked against the literature, against charge conjugation and against `Processes/Collins/6-ResolvedBorn.wl`); Born, virtual (families with HatPower, the normal component of the loop momentum reduced by symmetric integration over the (D−3)-dimensional orthogonal space) and real emission ((n·k1)² = k1∥² + spe3 through the Gram matrix) in the four components; XN = NX = 0 required; gauge invariance of the real emission checked |
| `2-Algebra-FORM.wl` | the same step with FORM doing the spin and colour algebra (`form/common.h`, `form/square.frm`, `form/checks.frm`: the density of a and the projection of c with the HVBM gamma5 as two epsilon tensors, the axes as dot-product substitutions, one FORM job per diagram pair and component); the same Born table and checks; the normal components reduced in Mathematica; ends with the exact comparison with the file of the FeynCalc route |
| `2-Algebra-Tracer.wl` | the same step with the classic Tracer package for the Dirac algebra and the explicit SU(N) colour reduction of `Processes/Hadronic/2-Algebra-Tracer.wl`; the density of a and the projection of c enter as the string `q, G5, s` when a fermion line is closed, the axes eX, eC and the normal as Tracer vectors; the same Born table and checks; ends with the exact comparison with the file of the FeynCalc route |
| `3-Reduce.wl` | `Processes/Hadronic`'s step 3 for the two components XX and NN: partial fractions, Kira (`kira/`), master integrals |
| `4-Transfer.wl` | threshold expansion, virtual normalization 2/(4π)^{D/2} with phase +1, derived from the FeynArts loop measure and the master-integral measure (checked by the universal double pole −(α_s/2π) Σ C_i × Born and, separately, by the cancellation against the real emission), MSbar renormalization, counterterms (transversity kernel on a and on the fragmenting c with the transfer Born table, unpolarized kernels on b with the lower transfer Borns), pole cancellation, scale-derivative checks (μ_F, μ_FF, μ_R), finite coefficients of both components, the equality of the class-1 and class-3 Born terms |
| `5-Crossing.wl` | the generalized crossing from `Processes/HadronicTT` (section below): Part A compares the crossed TT algebra of step 2 with this process's, exactly; Part B runs steps 3 and 4 on the crossed algebra as a pseudo-class and compares the NLO coefficients |
| `reference/Yuan-2008-Collins-hard-factors.wl` | the leading-order hard factors of arXiv:0709.3272, Eq. (6), transcribed verbatim (= the LO transverse spin transfer of q q′ → q, q q̄′ → q, q q → q, q g → q) |
| `kira/` | the Kira families of `Processes/Hadronic` (real and virtual) |

## How to run

Each notebook begins with `class = "1"`; set the class number (1, 3, 6, 10, 11,
13, 17), then run steps 1, 2, 3, 4 in fresh kernels.  Headless, one class at a
time, in one of the three routes FeynCalc, FORM and Tracer:

```sh
for step in 1-Amplitudes 2-Algebra-FeynCalc 3-Reduce 4-Transfer; do
  FEYNFACET_SET="route=FeynCalc;class=1" wolframscript -file Checks/run_step.wls Processes/HadronicTransfer/$step.wl
done
for step in 2-Algebra-FORM 3-Reduce 4-Transfer; do
  FEYNFACET_SET="route=FORM;class=1" wolframscript -file Checks/run_step.wls Processes/HadronicTransfer/$step.wl
done
for step in 2-Algebra-Tracer 3-Reduce 4-Transfer; do
  FEYNFACET_SET="route=Tracer;class=1" wolframscript -file Checks/run_step.wls Processes/HadronicTransfer/$step.wl
done
```

(`1-Amplitudes.wl` is common to the routes; the FORM and Tracer files of
step 2 end with a comparison with the FeynCalc file, in exact rational
arithmetic at two random rational points, which is skipped, with a message,
when that file is missing: run the FeynCalc route of a class first.)  Step 2 runs the products on parallel
kernels (`workers = "auto"`: one per processor core, at most 8); with
`route=FORM` the same products run as FORM jobs, with `route=Tracer` on
Tracer kernels, and steps 3 and 4 with the same override read the files of
that route.  Results in
`Results/<route>/class-<n>/`:
`Algebra.wl` (squares and the transfer Born table), `Reduced.wl`, `Transfer.wl`
(the finite coefficients with the keys "XX" and "NN", in the format of
`Processes/Hadronic`'s `Unpolarized.wl`).

Step 2 uses parallel algebra kernels; reduce `workers` when memory is
limited.  Step 3, including its use by the crossing calculation in step 5,
uses one Wolfram kernel with external Kira processes. For class 13,
sampled process-tree resident memory reaches approximately 13.4–14.4 GiB
in the crossing calculation, including its finite assembly. These periodic
measurements are lower bounds on the true peaks, not sufficient-memory
requirements. Reducing the algebra worker count does not limit the memory
of these stages. Measure them separately; see
[resource guidance](../../docs/RESOURCES.md).

Reduction design. Steps 3 and 4 use the same mathematical implementation
for each route, acting on that route's own `Algebra.wl`. FORM and Tracer
compare their algebra with FeynCalc at two exact-rational points, assigning
every scalar integral its own value. This is a sampled comparison.
All seven classes are reduced and assembled separately in all three routes.
The collector compares their finite coefficients at six sampled points.

## The generalized crossing from p↑ p↑ → c + X (step 5)

The transverse spin transfer cross sections and the double-transverse ones of
`Processes/HadronicTT` are related by generalized crossing symmetry:
knowing a↑ b↑ → c + X, the process a↑ c̄ → b̄↑ + X follows
by s ↔ t, v → −v/(1−v), w unchanged.  `5-Crossing.wl` checks it for every
class, at LO exactly and at NLO on the complete finite coefficients, against
the results of HadronicTT.

The map.  b ↔ c̄ is pb → −p, p → −pb with the unobserved momenta unchanged:
s ↔ t, u, A = (pa−k1)², spe3 unchanged, B = (pb−k1)² ↔ C = (p+k1)²; the
one-loop family VirtualST maps onto itself with s ↔ t, VirtualSU ↔ VirtualTU
(found on the propagator sets).  The pairs (transfer ← TT): 3 ← 5 (q q̄ → q′,
k1 = q̄′), 11 ← 11, 17 ← 6 (q q → q crosses into q q̄ → q̄, class 17 by charge
conjugation), 13 ← 12 (q q̄ → g), 10 ← 7 (q q → g, NLO only); 1 ← 5 and
6 ← 11 observe the antiquark of the TT pair, which by charge conjugation is
the TT square with a ↔ b, so they use the composed map pa → −p, pb → pa,
p → −pb (s → t, t → u, u → s, A → C, B → A, C → B).

The factor.  XX_transfer = −σ g XX_TT(crossed), NN_transfer = +σ g NN_TT(crossed)
with σ = +1, g = 1 when the transfer's parton b is a quark or antiquark and
σ = −1, g = 2C_A/((D−2)(C_A²−1)) = 1/((D−2)C_F) when it is a gluon.  The pieces:
the TT density p̸_b γ5 s̸_b/2 → −p̸ γ5 s̸/2 against the transfer projection
p̸ γ5 s̸_c, and the TT spin sum p̸ of the observed quark → −p̸_b against the
transfer's p̸_b/2 (two fermion crossing signs; only one when the TT observed
parton is a gluon, hence σ); the colour averages 1/N_c² against 1/N_c × 1/N_b;
the gluon polarization sum (reference vector pa on both sides) with the
transfer's 1/(D−2); and the axes: the crossed TT in-plane axis equals
(−s/t) e_X + (2u/t) pa as the spin vector of a and −((t+u)/t) e_C + (multiple of p)
as the spin vector of c, which with the normalizations gives −1 for XX (the
azimuth map φ_b → π − φ_c) and +1 for NN (the normal is reflected, but both
TT densities carry it).  Consequently the TT projections cross as
−σg × TT → (XX+NN)/2 (the coefficient of cos(φ_a − φ_c)) and
−σg × Symmetric → (XX−NN)/2, the difference of the two transfer components,
not their sum.  At LO, with b(v) = |M|²(4πα_s)²/(16π² s v):
b_transfer,XX(v) = σg b_TT,XX(−v/(1−v))/(1−v), b_transfer,NN(v) = −σg b_TT,NN(−v/(1−v))/(1−v),
and at D = 4 both equal σg b_TT(−v/(1−v))/(1−v) with b_TT the cos(2Φ)
coefficient of HadronicTT.  All of this is verified exactly on the
D-dimensional Born terms and on the complete Born table (all seven transfer
channels follow from the TT table).

Part A (exact, algebra level, per pair): Born and Born table (rational-function
identities); virtual XX and NN after both sides are reduced to master
integrals (term by term the scalar-integral coefficients agree only when the
two files route the graphs alike, which is recorded) and brought to a common form (the two
step-2 files reduce the normal component of the loop momentum differently:
every hatted power is replaced by its angular average
loopHat2^h → ((D−4)/2)_h/((D−3)/2)_h (loop⊥²)^h and both sides are reduced
with Kira in D dimensions; the masters are identified by their values, bubble
of a channel or box of two channels, and their coefficients must agree exactly); the
real-emission squares exactly at random rational points of all invariants.

Part B (the complete NLO check): the crossed TT algebra is written as a
pseudo-class `Results/<route>/class-<n>x/Algebra.wl` (the transfer class entry
of `classes.wl` as its specification, the crossed TT Born table as its
counterterm table) and steps 3 and 4 of this process are run on it in fresh
kernels (`RunSteps[..., "class=<n>x"]`, skipped when the results exist).  Step 3
integrates the crossed integrands in the transfer region, so no finite TT
expression is continued (no log(v) → log(−v/(1−v)) or π² bookkeeping); step 4
assembles them with the transfer counterterms (the interchange of the
incoming-parton and fragmentation counterterms) and requires the
poles to cancel, which tests the factor once more, including the pair 10 ← 7
without a Born term.  The finite coefficients {regular, D0, D1, delta} of both
components, and separately the finite real, virtual and counterterm pieces,
must equal those of the directly computed class, exactly at rational points of
s, v, w, N_f with symbolic scales.

For all seven pairs, the crossed and direct calculations use the same
`3-Reduce.wl` and `4-Transfer.wl` implementations, including their transfer
counterterm construction.  Part B checks their consistency for independently
constructed direct and crossed algebra; it does not provide an independent
implementation of the reduction or finite subtraction.

A separate finite-assembly implementation is provided in
[`Checks/transfer/`](../../Checks/transfer/README.txt). It independently
constructs the endpoint distributions, renormalization and collinear
counterterms from the existing scalar integrals and D-dimensional Born
tensors. The checks retain generic kinematics, color, flavor and separate
scales for all seven direct and crossed classes, including class 13.
They compare every regular, D0, D1 and delta coefficient for both XX and NN,
and check the four XX=NN identities per class. The 56 direct-crossed
coefficient comparisons check assembly consistency: the reduced real
integrals, Born tensors and master values already agree between the inputs;
the virtual expressions can differ in form. The nontrivial crossing evidence
comes from step 5 above. Polynomial perturbations test the comparison routine.
Four class-1 tests remove actual Born O(epsilon) or fragmentation-measure
terms, separately for XX and NN, and verify the resulting finite mismatch.
This supplies an independent subtraction and finite-assembly check using
shared amplitude algebra, reduction and master integrals. The
[derivation](../../Checks/transfer/DERIVATION.txt) gives the Jacobian,
evanescent Born and fragmentation-dimensional factors.

```sh
FEYNFACET_SET="route=FeynCalc;pairs=3" wolframscript -file Checks/run_step.wls Processes/HadronicTransfer/5-Crossing.wl   # one pair
FEYNFACET_SET="route=FeynCalc;pairs=all" wolframscript -file Checks/run_step.wls Processes/HadronicTransfer/5-Crossing.wl # all seven, sequentially
```

The evanescent part.  Since TT's Symmetric projection crosses into (XX−NN)/2,
the vanishing finite Symmetric remainder of HadronicTT (classes 5, 6, 7, 11)
predicts XX = NN for the finite transfer coefficients; indeed XX − NN of
`Transfer.wl` vanishes (with C_F = (C_A²−1)/(2C_A), checked to 30 digits at
the physical points by step 5) for every class checked, although the two
components differ in D dimensions.

`pairs` takes a transfer class, a comma-separated list or `all`; `stage=A`
stops after Part A and the pseudo-class algebra (steps 3 and 4 can then be
run by hand with `class=<n>x`, and a later run with `reuse=yes` reads them).
Inputs: `Processes/HadronicTT/Results/<route>/class-<m>/Algebra.wl` and this
process's `Algebra.wl` and `Transfer.wl` of the class.  Results in
`Results/<route>/crossing/class-<n>.wl` (the map, the factor, the family map,
every check of Parts A and B, the numerical differences of the finite
coefficients).

The class-13 native comparison cells can also be run separately using
[`Checks/native/`](../../Checks/native/README.txt). They retain the original
FORM angular projection and Tracer comparison, save the exact sample points,
and check every required tensor and Born-table entry without rerunning algebra.

## Results

The rows summarize the comparisons performed by steps 2–4 in the three routes
and by step 5.
"Squares" is the comparison at the end of step 2 (Born, virtual, real
emission and the transfer Born table, FORM and Tracer against FeynCalc, in
exact rational arithmetic at two random rational points); "Routes" the
comparison of the finite coefficients {regular, D0, D1, δ} of both components
from `Transfer.wl` (numerical, at the six random points of
`Checks/collect_coefficients.wls`); "XX = NN" the equality of the finite coefficients of the
two components; "Crossing" the check of step 5 (Part A: Born, Born table,
virtual after reduction and real emission exactly; Part B: the crossed TT
algebra through steps 3 and 4 of this process, poles cancel and the finite
coefficients equal those of the direct calculation exactly).  Poles cancel
and the scale-derivative checks pass in every completed row.

| Class | Process | Born | TT partner | Squares FeynCalc = FORM = Tracer | Routes (finite coefficients) | XX = NN | Crossing | Status |
|---|---|---|---|---|---|---|---|---|
| 1 | q q′ → q | yes | 5 | agree (2 points) | agree | yes | Part A exact; Part B equal | complete |
| 3 | q q̄′ → q | yes (= class 1) | 5 | agree (2 points) | agree | yes | Part A exact; Part B equal | complete; NLO terms differ from class 1 |
| 6 | q q → q | yes | 11 | agree (2 points) | agree | yes | Part A exact; Part B equal | complete |
| 10 | q g → q̄ | no | 7 | agree (2 points) | agree | yes | Part A exact; Part B equal (no Born term: the poles cancel between the real emission and the counterterms alone) | complete |
| 11 | q q̄ → q | yes | 11 | agree (2 points) | agree | yes | Part A exact; Part B equal | complete |
| 13 | q g → q | yes | 12 | agree (2 points; FORM's NN virtual after the loop-angle projection, see the comparison cell of `2-Algebra-FORM.wl`) | agree | yes | Part A exact; Part B equal | complete |
| 17 | q̄ q → q | yes | 6 | agree (2 points) | agree | yes | Part A exact; Part B equal | complete |

The finite coefficients of the completed classes are collected in
`Coefficients/ProtonProtonTransfer.wl`.

The crossing factors found in step 5, for every pair:
XX_transfer = −σ g XX_TT(crossed) and NN_transfer = +σ g NN_TT(crossed),
with σ = +1, g = 1 when the transfer's parton b is a quark or antiquark
(pairs 1 ← 5, 3 ← 5, 6 ← 11, 11 ← 11, 17 ← 6) and σ = −1,
g = 2C_A/((D−2)(C_A²−1)) when it is a gluon (10 ← 7, 13 ← 12).  At LO and
D = 4: b_transfer(v) = σ g b_TT(−v/(1−v))/(1−v), with b_TT the cos(2Φ)
Born coefficient of HadronicTT.  Since the TT coefficients of the partner
classes agree with the Fortran code of Mukherjee, Stratmann and Vogelsang,
the crossing ties this process to a published calculation.  Part A (Born and
Born table identically, the virtual identically after reduction, the
real-emission squares at three random exact rational points) and the
end-to-end Part B (the crossed TT algebra run through steps 3 and 4 of this
process: poles cancel, finite coefficients equal to the direct ones in exact
arithmetic at three rational points with symbolic scales) hold for all seven
pairs.  Step 5 needs steps 1–4 of the transfer class and step 2 of its TT
partner; one pair is run with
`FEYNFACET_SET="route=FeynCalc;pairs=13" wolframscript -file Checks/run_step.wls Processes/HadronicTransfer/5-Crossing.wl`
(Part B runs steps 3 and 4 of the pseudo-class `13x` itself, through `RunSteps`).

Classes 1 and 3 have equal Born terms (the t-channel graph alone) but
different NLO coefficients, as for the unpolarized classes 1 and 3 of
`Processes/Hadronic`: at one loop the box graphs and the real emission
distinguish q′ from q̄′.  (`classes.wl` describes class 3 as "equal to
class 1"; that holds for the Born term only.)

## Conventions in brief

a(pa) b(pb) → c(p) + X; s = (pa+pb)², t = (pa−p)², u = (pb−p)², v = 1 + t/s,
w = −u/(s+t), s4 = s v (1−w).  σ = s E dσ/d³p; NLO coefficient in the basis
{regular, D0, D1, delta} with the prefactor α_s³/(π s v w) and the Born
α_s² b δ(1−w), for each of the components XX and NN; scales μ_R (coupling),
μ_F (both incoming legs), μ_FF (fragmentation).  HVBM γ5 with four-dimensional
spin vectors; transversity kernel δP_qq(z) = C_F [2z/(1−z)_+ + 3/2 δ(1−z)] on
leg a and on the fragmenting quark, unpolarized kernels on leg b, MSbar coupling
renormalization, no helicity-restoring term.
