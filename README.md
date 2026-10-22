# FeynFacet: analytic NLO QCD calculations

*One Feynman-diagram calculation, many facets.*

Executable Mathematica notebooks that derive next-to-leading-order QCD
coefficient functions from Feynman diagrams, with three independent
implementations of the spin and colour algebra (FeynCalc, FORM, Tracer),
integration-by-parts reduction with Kira, and master integrals evaluated with
SubTropica. Known results are checked against the literature; the new
transverse spin-transfer results are checked across the three algebra routes
and, for proton-proton scattering, by generalized crossing. This is the
source code accompanying the paper
[*FeynFacet: AI agents for analytic NLO QCD*](https://arxiv.org/abs/2610.31745)
(Kang, Vogelsang, Zhang), arXiv:2610.31745.

## Quick start

An AI coding agent installs and checks the tools for you; you do not need to
set them up one by one.

1. Install Mathematica and clone this repository.
2. Start Codex (`codex --approve-for-me`) or Claude Code
   (`claude --permission-mode auto`) in the repository folder, and paste:

   ```text
   Set up FeynFacet on this machine. Install the tools it needs, following
   docs/INSTALL.md (docs/INSTALL-macOS.md on a Mac; on Apple silicon, build
   Kira with tools/install_kira.py as described there), and run
   python3 feynfacet.py check --all-routes until it reports READY.
   Ask me before anything that needs an administrator password.
   ```

   The agent asks before anything that needs your administrator password
   (on a Mac, for example, Homebrew and Apple's Command Line Tools).
   To install the tools by hand instead, follow `docs/INSTALL.md`
   (`docs/INSTALL-macOS.md` for a Mac) and run `python3 feynfacet.py check`.
3. Open `Start.nb` in Mathematica. Evaluate the loading and installation
   cells, then the lepton–quark example, which runs the process end to end
   and compares with the published coefficients. Select subsequent examples
   individually; the larger hadronic classes can require hours.
4. Open the notebooks in `Processes/LeptonQuark/` to read and evaluate the
   calculation step by step.

To compute a new process with an agent, see the Example on the project website.

To use the results without running anything, load the files in
`Coefficients/` (its README gives the formulas and an example).

The paper's comparison settings are supplied in
[Checks/validation-settings.wl](Checks/validation-settings.wl); see the
[checks documentation](Checks/README.md#comparison-settings).

See [memory and parallel runs](docs/RESOURCES.md) before running the larger
hadronic classes. Resource requirements depend strongly on the process and stage.

Project website: https://qcdphysics.github.io/FeynFacet/  
Source: https://github.com/qcdphysics/FeynFacet

No configuration file is needed when the tools are in their usual places;
`Paths.example.wl` explains the one-line override if a program is not found.

## Layout

```
FeynFacet/        shared package: tool lookup, amplitude records, Kira and FORM
                  interfaces, master integrals, threshold expansion, comparisons
Processes/        one folder per process, with numbered notebook files
  LeptonQuark/    l q -> q + X, unpolarized and double-longitudinal (NLO)
  LeptonOffDiagonal/  l q -> g + X and l g -> q + X (NLO, no Born term)
  Collins/        l q -> q + X with transverse spin transfer (NLO), resolved Born table
  PartonicBorn/   every 2 -> 2 Born channel in D dimensions (UU, LL)
  Hadronic/       the 16 proton-proton channel classes (NLO, UU and LL)
  HadronicTT/     p^ p^ -> h + X, both protons transversely polarized (A_TT, NLO)
  HadronicTransfer/  p^ p -> h^ + X, transverse spin transfer (NLO), and the
                  generalized crossing from HadronicTT
  TMD/            NLO TMD beam-function matching coefficients and soft function,
                  from a FeynArts model with Wilson-line and operator fields
Coefficients/     the final coefficients of every process as plain Mathematica
                  expressions, usable without any of the algebra tools
Checks/           headless runner used to test every notebook, and the script
                  that writes Coefficients/ from the notebooks' results
Reference/        transcribed reference codes used by the comparisons
docs/             conventions (CONVENTIONS.md; TRANSVERSE-PP.md for the
                  transversely polarized proton-proton processes) and
                  installation notes
Paths.example.wl  local tool configuration (copy to Paths.wl)
```

Each process folder contains numbered steps:

1. `1-Amplitudes.wl` – FeynArts diagrams and amplitudes, written as plain
   records (one coefficient and the open fermion chains per term);
2. `2-Algebra-<route>.wl` – squared matrix elements with FeynCalc, FORM or
   Tracer; all three write a common output format;
3. `3-Reduce.wl` – reverse unitarity, sparse scalar-integral construction,
   Kira reduction, and master integrals;
4. `4-Unpolarized.wl`, `5-Polarized.wl` (or `4-Transverse.wl`,
   `4-Transfer.wl`) – threshold expansion, counterterms, finite coefficient,
   comparison with the literature.

The numbered workflows read preceding results from `Results/` and write
their own. The runners use fresh kernels for successive calculations; a
workflow such as the crossing check also invokes further calculation steps.
The files are Wolfram Language package files with cell markers: open them in
Mathematica and they appear as notebooks with text and input cells; evaluate
the cells from top to bottom. They are also readable as plain text on GitHub;
`.nb` copies generated from them are provided next to each file.

## Requirements

Mathematica (tested with 15.0), FeynCalc 10.2 with FeynArts 3.12, Tracer 1.1,
FORM 5.0 with the colour library, Kira 3.1 with Fermat, SubTropica 1.2.10
(bundles HyperIntica) with polymake.  See `docs/INSTALL.md` for the
installation of each tool and for `Paths.wl`, and `docs/INSTALL-macOS.md`
for an Apple-silicon Mac. **Apple silicon requires the documented Kira
arithmetic patch**; unpatched Kira 3.1 release builds can silently give wrong
reductions. See the [correctness notice](docs/KIRA-CORRECTNESS.md).
After installing, open a notebook and evaluate

```wl
Get["/path/to/FeynFacet/FeynFacet/FeynFacet.wl"];   (* the clone of the repository *)
CheckInstallation[]
```

which lists every tool with its location and version.  To test that the tools
actually work together, run from a terminal

```sh
python3 feynfacet.py check               # the FeynCalc route
python3 feynfacet.py check --all-routes  # also FORM and Tracer
```

To use the kernel directly, add `--wolfram /path/to/WolframKernel`.

Each test performs a small real calculation in a fresh Wolfram kernel (a
FeynCalc trace, eight exact Kira identities through Fermat, a native polymake
calculation, SubTropica integrals including a one-loop box that needs polymake,
and FORM and Tracer traces), prints what failed with a suggested fix, and writes
its report to its own folder `Results/toolchain-check/run-<time>-<process>-<id>/`,
with an atomic copy of the latest completed report in `Results/toolchain-check/report.json`.
The ten most recently completed runs are retained; unfinished runs are preserved.
It installs nothing and runs no process
workflow. Each test has a 120-second timeout; `--timeout SECONDS` overrides it.

## Running

Interactively: open `Processes/LeptonQuark/1-Amplitudes.wl` in Mathematica
and evaluate it; continue with `2-Algebra-FeynCalc.wl`, `3-Reduce.wl` and
`4-Unpolarized.wl`, each in a fresh kernel.  The process README explains the
physics of each step.

All steps of one route from inside Mathematica, each in a fresh kernel:

```wl
RunSteps["LeptonQuark", "FeynCalc"]
```

or, from a terminal, `python Checks/run_process.py Processes/LeptonQuark --route FeynCalc`.
The TMD process runs the same way with `"TMD"` in place of `"LeptonQuark"`
and writes regenerated figures to `Processes/TMD/Results/figures/`.

A run stops at the first failed check and prints the cell that failed;
`Results/<step>-<route>[-<choice values>].log` contains the complete kernel output
for the selected route and settings.

## Calculations and checks

| Process | Routes | Spins | Literature comparison |
|---|---|---|---|
| LeptonQuark | FeynCalc, FORM, Tracer | U, LL | exact agreement with arXiv:1505.06415 and arXiv:1703.10872 |
| TMD | FeynCalc, FORM, Tracer | U | beam matching coefficients I_qq, I_qg, I_gq, I_gg, soft functions S_q, S_g and soft-subtracted TMDPDF matching coefficients: exact agreement with the ancillary files of arXiv:2012.03256 |
| LeptonOffDiagonal | FeynCalc, FORM, Tracer | U, LL | ℓq → g + X and ℓg → q + X: exact agreement with arXiv:1505.06415 Eqs. (27), (28) and arXiv:1703.10872 Eqs. (13), (14) |
| Collins | FeynCalc, FORM, Tracer | U, transverse spin transfer | ℓq↑ → q↑ + X (the Collins hard part) and the resolved-photon Born table: comparison with the paper; the unpolarized case re-derived as a check |
| PartonicBorn | FeynCalc, FORM, Tracer | UU, LL | all 2 → 2 Born channels in D dimensions, checked against the textbook table and used by Hadronic and Collins; the FORM and Tracer tables equal the FeynCalc one exactly (any D and N_c) |
| Hadronic | FeynCalc, FORM, Tracer | UU, LL | the 16 proton–proton channel classes: all four coefficients (regular, D0, D1, δ(1−w)) of every class, for both spin cases, agree with the reference codes of Aversa et al. (unpolarized) and Jäger, Schäfer, Stratmann, Vogelsang (polarized) at nine kinematic points each (`Reference/pp-fortran-values.wl`); the FORM and Tracer squares agree with the FeynCalc ones for every class (exact rational arithmetic at two random rational points), and so do their coefficients (numerically, at six random points) |
| HadronicTT | FeynCalc, FORM, Tracer | TT (both incoming partons transversely polarized) | the five contributing classes of p↑p↑ → h + X: the Born terms and all four finite coefficients are checked against Mukherjee, Stratmann, Vogelsang (hep-ph/0506315; `Reference/pp-transverse-fortran.wl`) at sampled SU(3) kinematic points; additional C_A = 2, 5 comparisons are reported separately. All five classes are reduced and assembled separately in every route |
| HadronicTransfer | FeynCalc, FORM, Tracer | transverse spin transfer a↑ b → c↑ | new NLO results for all seven classes (1, 3, 6, 10, 11, 13, 17). All seven classes are reduced and assembled separately in every route; the collector compares their coefficients at six sampled points. Step 5 checks generalized crossing from HadronicTT. An additional [independent finite assembly](Checks/transfer/README.txt) checks the subtraction and finite coefficients for all seven classes, including class 13, and the finite XX = NN identities, using the same reduced scalar integrals. |

## Algebra and reduction

The reductions keep additive fractions separate and combine scalar
coefficients multiplying identical integration monomials before constructing
IBP targets. Partial-fraction and family maps are reused for identical monomials.
Every contribution retains its reconstruction check, integral-family indices
and dimension-shift label. Kira targets are generated from the current rows;
coverage is checked before the reduced result is expressed in master integrals.
Equivalent partial-fraction representations can require different target
integrals. Their reduced coefficients are compared in a common master basis.

The FeynCalc implementation groups identical algebra structures and uses bounded
caches for repeated contractions. The four transverse-spin projections remain
separate. The trace dimension, HVBM prescription, polarization projectors and
master-integral formulas are unchanged by these representation choices. Runtime
improvements depend on the channel; amplitude algebra, reduction and final
assembly should be timed separately.

## Threshold preparation and endpoint expansion

This procedure applies to NLO single-inclusive partonic cross sections with
recoil invariant `s4 = s v (1-w)` and threshold `w = 1`. Reuse it when
developing new processes with the same kinematics and threshold structure.
After evaluating the real-emission integrals, prepare and expand the
threshold factors with the shared routines in
`FeynFacet/Distributions.wl`, used by all three algebra routes:

1. **Prepare the threshold branches.** Write each master integral as a sum
   of terms with a single exact threshold power. Keep factors such as
   `(1-w)^(-1-lambda*eps)` unexpanded, and use `ThresholdSplit` to separate
   the regular part from these branches. Retain inert master functions
   where the process template does so. Determine the required epsilon
   orders from the full pole structure: higher-order epsilon terms can
   contribute to the finite answer when multiplied by poles.
2. **Expand into distributions.** Use `EndpointExpand` to obtain the
   epsilon poles and finite `{regular, D0, D1, delta}` coefficients, with
   `D0 = [1/(1-w)]_+`, `D1 = [log(1-w)/(1-w)]_+` and
   `delta = delta(1-w)`. It groups identical epsilon-dependent structures,
   expands each once per required order, and reuses endpoint series of
   identical w-dependent structures, including logarithms and master
   functions. Reuse this implementation to avoid expanding a complete
   bulky branch or forming its global common denominator.
3. **Validate the result.** Check exact pole cancellation after combining
   real, virtual and counterterm contributions. When changing the expansion
   implementation, run `Checks/test_distributions.wls` and compare affected
   coefficients exactly where tractable. These helper tests complement the
   process-specific physics checks.

In these single-inclusive workflows, `3-Reduce.wl` prepares the threshold
branches of the real-emission master integrals,
and the assembly steps (`4-...`, `5-...`) apply `ThresholdSplit` and
`EndpointExpand`. The virtual correction is evaluated separately at Born
kinematics and contributes to the `delta(1-w)` coefficient. The assembly then
combines the virtual correction, the expanded real contribution and the
renormalization and factorization counterterms and checks cancellation of all
poles. The smallest worked example is LeptonQuark: see the master-integral
section of
[3-Reduce.wl](Processes/LeptonQuark/3-Reduce.wl), followed by “Real emission
and its threshold expansion” in
[4-Unpolarized.wl](Processes/LeptonQuark/4-Unpolarized.wl). Function signatures
and coefficient ordering are documented in
[Distributions.wl](FeynFacet/Distributions.wl), with further details in
[Checks/README.md](Checks/README.md#distribution-expansion).

Small exact helper checks can be run without regenerating amplitudes:

```sh
wolframscript -file Checks/test_preprocessing.wls
wolframscript -file Checks/test_sparse.wls
wolframscript -file Checks/test_algebra_reuse.wls
wolframscript -file Checks/test_distributions.wls
wolframscript -file Checks/test_script_arguments.wls
```

## License

The original FeynFacet code is released under the [MIT license](LICENSE).
Third-party software and reference material retain their respective terms
and attribution. Please cite the accompanying paper,
[arXiv:2610.31745](https://arxiv.org/abs/2610.31745), when using FeynFacet in your research.
