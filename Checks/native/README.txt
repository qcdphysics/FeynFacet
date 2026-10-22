Transverse spin-transfer algebra comparisons

After computing HadronicTransfer class 13 through FeynCalc and the route to
be compared, run from the repository root:

 export FEYNFACET_SOURCE="$PWD"
 export FEYNFACET_BENCHMARK_OUTPUT="/path/to/new/form-comparison"
 wolframscript -file Checks/native/TransferRouteFORM.wls

For TRACER, choose another new output directory and run
Checks/native/TransferRouteTracer.wls. The output directory must not already exist.
The scripts use built-in Wolfram algebra; they do not load FeynCalc or TRACER.

Each benchmark replays the final comparison cell of its corresponding
2-Algebra-<route>.wl file. It checks both XX and NN, the Born, virtual and real
parts, and all seven Born-table channels: 20 scalar quantities, each evaluated
at two exact-rational points (40 evaluations). Missing channels, missing inputs and skipped comparisons
are errors. A held-expression check verifies that the copied comparison
block matches the source after removing the one point-saving statement.
For FORM, this guard covers the final comparison block, not the auxiliary
hat-projection definitions preceding it. Those definitions are copied from
the notebook and match this source version; the guard does not check them.

The fixed seeds are 2026101701 (FORM) and 2026101702 (TRACER). The exact points
are written to exact-points.wl. comparison.json contains the verdicts, file
hashes, sample-point hash and Wolfram version. Thus the samples can be checked
and regenerated without relying on a summary of a previous run.

FORM reports both literal and angularly projected virtual comparisons. The
NN literal comparison can be false because the two integrands differ by a
term whose angular projection vanishes. The projected comparison is required
to pass; both verdicts remain visible. These are sampled algebra comparisons,
not a symbolic equality claim or separate scalar integrations.
