# Checks: running the notebooks without the front end

## Toolchain correctness

`python3 feynfacet.py check --all-routes` performs small functional checks,
including eight exact Kira identities. It does not run a physics workflow or
rebuild tools. See [the Kira fixtures and derivations](kira/README.md) and
[the Apple Silicon correctness notice](../docs/KIRA-CORRECTNESS.md).

The native polymake test exercises both `polytope` and `fan`; the SubTropica
test checks a one-loop box against its exact pole coefficients and requires
a recorded polymake call. Every calculation also checks external launches
against complete per-call records. Nonzero exits, reported errors, missing
records and unfinished calls abort and retain the integral's workspace.
SubTropica's tropical-data cache is isolated for each integral, so a failed
call cannot poison a later retry. A legitimate zero or an integral needing
no polymake call is accepted.

Each check gets a unique directory. Reports are published atomically; cleanup
retains the ten most recently completed checks and leaves unfinished runs
alone. Python regressions exercise overlapping checker processes and the
installer's dependency selection:

```sh
python3 Checks/test_toolchain_runner.py
python3 Checks/test_kira_installer.py
```

After a successful toolchain check, run the backend failure regressions with
a new scratch directory. They deliberately break polymake after a healthy
box, test missing records and same-kernel retries, and retain the expected
failure diagnostics:

```sh
wolframscript -file Checks/test_masters.wls /tmp/feynfacet-master-regression-1
```

The `Kira arithmetic` GitHub workflow tests the actual upstream arithmetic
source before and after the patch on Linux x86-64 and macOS ARM64. This fast
source-level test needs no Wolfram or Fermat installation; it does not test
an installed Kira, so the exact reduction checks of `feynfacet.py check`
remain necessary on every machine.

## Process runner

For the `.wls` commands below, `WolframKernel -script` can replace
`wolframscript -file`. The scripts use `ScriptArguments.wl` to read user
arguments with either launcher and `$InputFileName` to locate their sources.
`test_script_arguments.wls` checks both command-line formats, including
optional arguments and paths containing spaces.

`run_process.py` evaluates the numbered step files of a process one after
another, each in a fresh Wolfram kernel, exactly as a user would in the
notebook interface:

```sh
python Checks/run_process.py Processes/LeptonQuark                 # FeynCalc route
python Checks/run_process.py Processes/LeptonQuark --route FORM
python Checks/run_process.py Processes/LeptonQuark --route Tracer
python Checks/run_process.py Processes/LeptonQuark --steps 3,4      # only some steps
python Checks/run_process.py Processes/LeptonQuark --copy /tmp/fresh --paths ~/Paths.wl
python Checks/run_process.py Processes/LeptonQuark --steps 4,5 --plan   # list what would run
```

`--route` selects which `2-Algebra-*.wl` file runs and sets the variable
`route` in the later steps.  `--set name=value` overrides any top-level
string choice in the files (for instance `--set spin=Polarized` where a
notebook has `spin = "Unpolarized"`).  `--copy DIR` first makes a copy of
the repository without any `Results/` directory and runs there: the
cleanest test that the shipped sources are complete (the process folder is
given relative to the repository root, and with `--copy` the folder of the
copy is run).  `--plan` lists the steps that would run, with the receipt
each must write, and runs nothing.

The runner fails (non-zero exit code) when the process folder does not
exist, when `--steps` selects no step file, or when a step does not end with
a completed receipt written by this run: every step must exit with code 0
and leave `Results/<step>-<override values>.receipt.json` with
`"completed": true`, written after the run started and naming that step file.

`run_step.wls` is the per-file driver: it splits the file into its input
cells, evaluates them in order, and stops at the first cell that produces
a message, a `Failure`, `$Failed`, `False` or an abort (the `Require`
checks in the notebooks abort with an explanation).  The kernel output of
each step goes to `Results/<step>.log`, and a JSON receipt with the time
and result of every cell to `Results/<step>-<override values>.receipt.json`
(`"completed": true` once the last cell has passed).  The receipt's
`"provenance"` records the Git commit (and which sources differ from it) and
the SHA-256 of the sources (the step file, the package `FeynFacet/`, the other
files of the process folder except notebooks and Markdown, `run_step.wls`
and its `ScriptArguments.wl` helper),
of every repository file the cells read (`Get`, `Import`; hashed when read)
and of writes observed by its `Put`, `Export` and `SaveResult` hooks;
the working files of Kira (`kira/` directories) are left out.  These hooks
do not capture every write by external programs or parallel kernels.
Use `run_recorded_step.py` below when a complete scoped output inventory
and an observed process exit are needed.

## Distribution expansion

`EndpointExpand` groups identical epsilon-dependent factors and expands each
distinct factor once per required order. Before taking endpoint limits it
also collects identical w-dependent structures, including logarithms and
master functions, and reuses their series. The caches belong to one call
and are shared by all its threshold
branches; the same routine serves FeynCalc, FORM and Tracer outputs. Scalar
coefficients and dimensional powers are kept separate until the expansion,
avoiding a large common denominator for the complete branch.

The distribution identities and required epsilon orders are unchanged.
Intermediate endpoint and epsilon poles are retained until their exact sum
cancels. An uncancelled pole outside the documented NLO form stops the
calculation. Run the exact regression tests, which compare against direct
whole-branch expansions of small expressions, with Mathematica only:

```sh
wolframscript -file Checks/test_distributions.wls
```

The tests include multiple branches, dimensional powers, Gamma functions,
mixed denominators, removable endpoint singularities, finite delta terms
from O(epsilon) coefficients, and rejection of unsupported inputs.

## Comparison settings

[validation-settings.wl](validation-settings.wl) collects the fixed comparison
points, scales, colour and flavour settings, precision and acceptance criteria
used in the paper, together with the scripts defining the algebra probes.
From the repository root, load the association with
`Get["Checks/validation-settings.wl"]`. It documents the existing checks;
the comparison scripts remain executable independently.

## Recording a fresh algebra step

`run_recorded_step.py` preserves the Wolfram cell receipt and separately
records the observed process exit and SHA-256 inventory of every new or
changed file in one output scope, including FORM and parallel-kernel
outputs.  The scoped record complements the Wolfram receipt; it does not
replace the cell checks or turn reused files into newly computed outputs.

After generating the class-13 amplitudes in step 1, run its FeynCalc algebra
from the repository root with an empty `Results/FeynCalc/class-13/` scope:

```sh
python3 Checks/run_recorded_step.py Processes/HadronicTransfer/2-Algebra-FeynCalc.wl \
  --scope Processes/HadronicTransfer/Results/FeynCalc/class-13 \
  --set route=FeynCalc --set class=13 --set workers=2 --set resume=no \
  --require-empty --record ../execution-records/transfer13-FeynCalc
```

The evidence directory must be new and outside the source repository.
`--require-empty` rejects pre-existing products and requires `resume=no`;
use a clean calculation copy when preserving an earlier result.  Concurrent
recorded runs must use disjoint scopes.  The wrapper locks overlapping
scopes and shared receipt names; its file locking is supported on Linux
and macOS.  Records include local paths and execution output: keep them
outside a source release.

When a process needs the results of another one, `run_process.py` runs that
one first if its results are missing (`Processes/Hadronic` reads the Born
table of `Processes/PartonicBorn`).

## Native transverse-spin algebra comparisons

[native/README.txt](native/README.txt) gives standalone FORM/FeynCalc and
Tracer/FeynCalc comparisons for hadronic transverse-transfer class 13.
They replay the production comparison expressions and test 20 scalar quantities
at each of two saved exact-rational points: 40 evaluations per route,
including all seven Born-table channels. The held-source guard covers the final
comparison block; FORM's auxiliary hat-projection definitions are copied
separately, as described in [native/README.txt](native/README.txt).
FORM records both the literal and
the required angularly projected virtual comparisons; its NN literal result
can differ while the projected result agrees. These algebra checks complement
the separate [finite-assembly checks](transfer/README.txt).

## Independent transverse-spin-transfer assembly

[transfer/README.txt](transfer/README.txt) explains a separate implementation
of the endpoint expansion, renormalization and collinear factorization for
the seven hadronic transverse-spin-transfer classes. The two standalone
Wolfram programs use the existing Born, algebra and reduced scalar-integral
inputs; they load no production finite-assembly or distribution helper.

The full check reconstructs regular, D0, D1 and delta coefficients for XX
and NN, including class 13, at generic kinematics and separate scales.
It checks the four XX=NN finite-coefficient identities for each class and also
compares the direct and crossed assemblies in 56 coefficient slots. The
reduced real integrals, Born tensors and master values already agree between
those inputs, so this last comparison checks assembly consistency. The
nontrivial crossing tests are in production step 5. Polynomial perturbations
test the comparison routine; four class-1 omission tests instead remove actual
Born O(epsilon) or fragmentation-measure terms and check that the finite result
changes. See [transfer/DERIVATION.txt](transfer/DERIVATION.txt) for the
distribution and dimensional factors. Amplitude algebra, scalar reduction and
master integrals remain shared mathematical inputs.

## Independent lepton factorization scales

After running the lepton process steps, `check_lepton_scales.wls` derives
PDF, fragmentation and photon-in-lepton scale residues for all six UU/LL
QQ, QG and GQ channels.  It checks the saved collinear poles and derivatives
of the computed results and original coefficient tables.  To check all
three algebra routes, with a new output directory outside the source tree:

```sh
wolframscript -file Checks/check_lepton_scales.wls . ../lepton-scale-check All
```

The separate-scale module is loaded explicitly with
`Get["FeynFacet/LeptonScaleDependence.wl"]`.  Its coupling-inclusive residues
are supplied in `Coefficients/LeptonScaleResidues.wl`, an auxiliary dataset
with its own schema and `AUXILIARY-MANIFEST.json`; the original eight
coefficient tables are unchanged.  See
[LEPTON-SCALES.txt](../docs/LEPTON-SCALES.txt) for the derivation and usage.

## Analytic qg endpoint benchmark

`AversaEndpoints.wls` checks the unpolarized `q g -> q` Born term and all
three endpoint coefficients against Aversa et al., Nucl. Phys. B327 (1989)
105, Eqs. (41) and (43), with the finite quark PDF/FF scheme conversion of
Eq. (23).  The identities retain generic kinematics, colour, flavour count
and separate scale logarithms.  It reads `Coefficients/ProtonProton.wl`
and needs only Mathematica; no process calculation or Fortran execution
is required.  From the repository root:

```sh
FEYNFACET_SOURCE="$PWD" FEYNFACET_BENCHMARK_OUTPUT="../aversa-check" \
  math -noprompt -script Checks/AversaEndpoints.wls
```

See [AVERSA-ENDPOINTS.txt](AVERSA-ENDPOINTS.txt) for normalization, the
published scheme convention and finite-only sensitivity controls.
Keep the generated `endpoints.wl` and `summary.json` outside the source tree.

## Soft, collinear and identical-gluon checks

`UnresolvedLimits.wls` reads the FeynCalc algebra outputs for ordinary
hadronic classes 1 and 15 and transverse-transfer class 1. It proves the
four-dimensional soft and incoming-quark collinear limits at generic
kinematic ratios for UU, LL and transfer XX/NN. It also compares exchange
of the two unobserved gluons in `g g -> g g g` at three saved exact-rational
points for UU/LL; that exchange comparison is sampled. Run after the
corresponding algebra steps, keeping outputs outside the source tree:

```sh
FEYNFACET_SOURCE="$PWD" FEYNFACET_BENCHMARK_OUTPUT="../unresolved-check" \
  math -noprompt -script Checks/UnresolvedLimits.wls
```

See [UNRESOLVED-LIMITS.txt](UNRESOLVED-LIMITS.txt) for the colour derivation,
normalizations, points and controls. This check does not require reduction
or master-integral evaluation.

`QGExchange.wls` separately reconstructs `q g -> q g g` from the class-13
FeynCalc interference caches. It checks all 136 upper-triangle products per
UU/LL spin projection (12 exact zeros and 124 proved nonzero expressions),
and tests identical-gluon exchange at the same three exact-rational points.
It requires those pair caches rather than the combined real-emission sum:

```sh
FEYNFACET_SOURCE="$PWD" FEYNFACET_BENCHMARK_OUTPUT="../qg-exchange-check" \
  math -noprompt -script Checks/QGExchange.wls
```

See [QG-EXCHANGE.txt](QG-EXCHANGE.txt) for the interference weights, map,
points and finite negative control.

`make_notebooks.wls` converts the `.wl` step files into `.nb` notebooks:

```sh
wolframscript -file Checks/make_notebooks.wls                 # all processes
wolframscript -file Checks/make_notebooks.wls Processes/LeptonQuark
```

`collect_coefficients.wls` writes the final coefficients of all processes to
`Coefficients/` after the notebooks have run in all three routes: it checks
that the routes agree, rewrites every coefficient in a canonical form
(logarithms of positive factors, scales as log(μ²/s), colour in N_c) and
checks the canonical form against the computed one.  These checks are
numerical, with 40-digit arithmetic at six random physical points (three with
N_c = 3, three with N_c = 5): sampled comparisons, not symbolic identities.
It also writes `Coefficients/MANIFEST.json`, which binds every coefficient
file to the relative paths and hashes of its calculation inputs and records
the hashes of the distributed scientific source files. Execution receipts
stay with the local `Results/` products.

Native hadronic FORM and Tracer comparisons use unseeded random exact-rational
probes whose coordinates and integral assignments are not saved; a rerun uses
new probes. The collector uses fixed seeds for its six finite-coefficient points
and its algebra probes for routes with shared scalar stages. These are
reproducible additional comparisons; they do not replay the original native
probes.

```sh
wolframscript -file Checks/collect_coefficients.wls
```

To refresh the manifest of the existing tables without collecting new
calculation products:

```sh
wolframscript -file Checks/collect_coefficients.wls manifest-only
```

This retains each table's calculation-input hashes and verifies its file hash;
it refreshes the distributed-source hashes and the collector's own hash.

## Comparing coefficient tables and file integrity

`compare_coefficient_tables.wls` compares every mathematical entry of the
eight process tables in two directories using exact symbolic arithmetic.
The `Routes` and `StepsThreeFour` execution metadata are excluded. For example:

```sh
wolframscript -file Checks/compare_coefficient_tables.wls tables-a tables-b ../comparison.json
```

The JSON report gives the individual comparisons; a nonzero exit denotes a
failed comparison. This compares the supplied tables; it does not regenerate
their amplitudes or master integrals.

`SHA256SUMS` covers the distributed source and data files. From the repository
root, check it with `python3 Checks/file_checksums.py --check`. After
intentionally changing those files, regenerate it with:

```sh
python3 Checks/file_checksums.py --write
```

Locally generated `Results/`, website output, logs and `Paths.wl` are excluded.

## Representation helpers

`test_preprocessing.wls` checks termwise denominator mapping, Laurent grouping
and reduction coverage. `test_sparse.wls` checks reuse of partial-fraction and
integral-family maps, exact reconstruction and failure propagation.
`test_algebra_reuse.wls` checks ordered-structure grouping, bounded-cache reuse
and failure handling. These small exact tests complement the amplitude and
coefficient comparisons; they do not replace a process calculation.
