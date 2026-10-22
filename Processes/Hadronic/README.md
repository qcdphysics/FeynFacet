# Proton–proton channel classes: a b → c + X at NLO

The NLO QCD corrections to single-inclusive parton production in
proton–proton collisions, for the spin average (UU) and the double-longitudinal
helicity difference of the incoming partons (LL), organized in the 16 tagged
channel classes of the reference codes (`classes.wl` lists them: the Born
process, whether a virtual correction exists, and the real-emission processes
with their weights).  The results are compared point by point with the
reference codes of Aversa, Chiappetta, Greco and Guillet (unpolarized) and
Jäger, Schäfer, Stratmann and Vogelsang (polarized), whose values at nine kinematic
points per class are stored in `Reference/pp-fortran-values.wl` (see
`Reference/fortran/README.md` for how they were produced).

## Files

| File | What it does |
|---|---|
| `classes.wl` | the 16 classes: Born and real-emission processes with their weights |
| `1-Amplitudes.wl` | FeynArts diagrams and amplitudes of one class (Born, one-loop, real emission with closed quark loops counted by N_f) |
| `2-Algebra-FeynCalc.wl` | squares with FeynCalc: Born (checked against `Processes/PartonicBorn`), one-loop interference written in three scalar-integral families, real-emission squares (gauge invariance checked by changing the gluon reference vectors).  Each amplitude is simplified once, colour is separated and evaluated once per colour network, only the products i ≤ j are formed, gluon polarization sums are applied one gluon at a time, and the products run on parallel kernels with one file per product (`Results/FeynCalc/class-<n>/pairs/`) |
| `2-Algebra-FORM.wl`, `form/` | the same squares with FORM (route `FORM`, results in `Results/FORM/class-<n>/`): Mathematica writes every diagram as colour structures times Lorentz parts (`work/<process>/amp<i>.inc`, `bar<j>.inc`, `process.h`); `form/square.frm` forms one product M_i M_j^* per job (transversality of the gluon polarizations, colour with the package `color.h` before the Dirac algebra, traces with the HVBM γ5 of `form/common.h`, polarization sums and circular densities, kinematics), at most `workers` jobs at a time. The last cell compares the file with the FeynCalc route's in exact rational arithmetic at two random rational points of all variables, every scalar integral given its own random value (a sampled test, not a symbolic identity); it passes for all 16 classes.  If the FeynCalc file of the class is missing, the comparison is skipped and the cell says so ("the FeynCalc route has not been run for this class"). |
| `2-Algebra-Tracer.wl` | the same squares with the Tracer package (route `Tracer`, results in `Results/Tracer/class-<n>/`, run in a kernel without FeynCalc): the amplitude records become Tracer expressions (fermion chains as inert gamma words), colour is reduced by an explicit SU(N) algorithm written in the file (f through traces, Fierz identity, closed loops; checked on known colour sums and on the Born table), the fermion lines are closed into traces, the gluon polarization sums are contracted into the words before the traces are taken, and the products run on parallel kernels with one file per product. The last cell compares the file with the FeynCalc route's in exact rational arithmetic at two random rational points of all variables, every scalar integral given its own random value (a sampled test, not a symbolic identity); it passes for all 16 classes.  If the FeynCalc file of the class is missing, the comparison is skipped and the cell says so ("the FeynCalc route has not been run for this class"). |
| `3-Reduce.wl` | partial fractions into six cut families and the loop families, Kira reduction (`kira/`), master integrals: phase-space volumes, angular integrals (including the mixed families with their regular hypergeometric branch), bubbles and boxes from SubTropica, analytic continuation to the physical region |
| `4-Unpolarized.wl` | threshold expansion, virtual normalization 2/(4π)^{D/2} with phase +1, derived from the FeynArts loop measure and the master-integral measure (checked by the universal double pole −(α_s/2π) Σ C_i × Born and, separately, by the cancellation against the real emission), coupling renormalization, collinear counterterms of both incoming legs and of the fragmenting parton from the Born table, pole cancellation, finite coefficient, comparison with the reference code |
| `5-Polarized.wl` | the same for the helicity difference (HVBM scheme, polarized kernels, helicity-restoring term) |
| `kira/` | Kira families for the real (`real_AB` … `real_BCp`) and virtual (`VirtualST`, `VirtualSU`, `VirtualTU`) integrals |

## How to run

Each notebook begins with `class = "13"`; set the class number (1–16), then
run steps 1, 2, 3, 4 and 5 in fresh kernels (step 4 is the spin average,
step 5 the helicity difference).  Headless, one class at a time:

```sh
python Checks/run_process.py Processes/Hadronic --set class=1                  # FeynCalc route
python Checks/run_process.py Processes/Hadronic --set class=1 --route FORM     # FORM route
python Checks/run_process.py Processes/Hadronic --set class=1 --route Tracer   # Tracer route
```

Each route keeps its own results in `Results/<route>/class-<n>/` (step 2
writes one file per diagram product in `pairs/`).  Run `Processes/PartonicBorn`
first: its Born table is read by steps 2 and 4.  Step 2 runs the products in
parallel (`workers = "auto"`: one per processor core, at most 8; set "1" for
none).  Runtime and memory use depend on the channel, algebra route, and
worker count; the larger classes can take hours.  With limited memory, start
with two workers and increase the count only when memory permits.  The
command-line runner records each cell's time in its execution receipt.

## Results

The three routes give the same squares for every class: each route's
step 2 ends with a comparison with the FeynCalc file, in exact rational arithmetic at two random rational points of all variables, every scalar integral given its own random value (a sampled test, not a symbolic identity); it is
skipped, with a message, when the FeynCalc file of the class is missing, so
run the FeynCalc route first.  The finite coefficients of the three routes
agree at the six random points (N_c = 3 and 5) at which
`Checks/collect_coefficients.wls` compares them numerically.  Every class agrees with the reference codes at the nine kinematic points, in
all four coefficients (regular, D0, D1, δ(1−w)), for the spin average and the
helicity difference.  The relative differences are at the level of the
double-precision reference values (10⁻¹⁵ to 2·10⁻⁷).

| Class | Spin average (Aversa et al.) | Helicity difference (Jäger, Schäfer, Stratmann, Vogelsang) |
|---|---|---|
| 1  q q′ → q, 3  q q̄′ → q, 5  q q̄ → q′, 6  q q → q, 11  q q̄ → q | agree (≤ 2·10⁻⁷, class 1; ≤ 3·10⁻⁸ the others) | agree (≤ 4·10⁻¹⁵) |
| 2, 4, 7, 8, 9, 10 (no Born term) | agree (10⁻¹⁶) | agree (10⁻¹⁶) |
| 13  q g → q | agrees (≤ 10⁻⁸) | agrees (10⁻¹⁵) |
| 14  q g → g | agrees (10⁻¹⁵) | agrees (10⁻¹⁵) |
| 15  g g → g | agrees (10⁻¹⁵) | agrees (≤ 5·10⁻¹⁵) |
| 16  g g → q | agrees (≤ 10⁻⁸); at v = 4/5 the D0 coefficient nearly vanishes (3·10⁻⁵ against 160 for the regular term) and agrees within the absolute precision of the double-precision reference (see the comparison cell of step 4) | agrees (≤ 10⁻¹¹) |
| 12  q q̄ → g | agrees (10⁻¹⁵) | agrees (≤ 3·10⁻¹⁵) |

Two lessons of this comparison are built into the step files:

- In the HVBM scheme a hatted loop momentum shifts only the dimension of the
  integral (I_D[l̂² F] = ε I_{D+2}[F] for the normalized loop measure), never
  the dimension D of the coefficients that the trace produced.  Shifting both
  changes the δ(1−w) coefficient of the helicity difference by a rational,
  Born-proportional amount (7/3 (α_s/2π) × Born for the t-channel classes);
  `3-Reduce.wl` therefore sets D = 4 − 2ε in the coefficients before the
  reduction rules, which carry their own dimension symbol d.
- A gluon whose momentum is eliminated by momentum conservation needs its
  transversality condition in the eliminated form, and a cached contraction
  must not be reused between the Born and the real kinematics (the scalar
  products of the same symbols differ).

## Conventions in brief

a(pa) b(pb) → c(p) + X; s = (pa+pb)², t = (pa−p)², u = (pb−p)², v = 1 + t/s,
w = −u/(s+t), s4 = s v (1−w).  σ = s E dσ/d³p; NLO coefficient in the basis
{regular, D0, D1, delta} with the prefactor α_s³/(π s v w) and the Born
α_s² b δ(1−w); scales μ_R (coupling), μ_F (both incoming legs), μ_FF
(fragmentation).  The comparison with the reference codes is done in the
reference codes' own normalization (`4-Unpolarized.wl`, section "Comparison"):
{R, B, C, A} = α_s³/(π s CC) {Regular + (AVWPL − AVWPL1)/(1−w) + (AVLO − AVLO1) ln(1−w)/(1−w),
AVWPL1, AVLO1, AVDEL + CC·AVDELMU ln(s/μ_R²)}, with the Born F0 = π b/(1−v).
