# Working with FeynFacet

Before any calculation, run `python3 feynfacet.py check` (add `--all-routes`
when FORM or Tracer is used). Resolve failures from its report, rerun it, and
consider the toolchain ready only when it reports READY. The Wolfram kernel,
Kira and Fermat write to their user and temporary directories; if a sandbox
blocks them, ask the user for permission rather than changing the source.

Start with [README.md](README.md), [the conventions](docs/CONVENTIONS.md),
and the README and numbered `.wl` files of the closest existing process.
Use its shared routines, runner and applicable physics checks when developing
a new calculation. The `.wl` workflows are authoritative; regenerate their
`.nb` copies with `Checks/make_notebooks.wls` after editing them.

For NLO single-inclusive partonic cross sections with recoil invariant
`s4 = s v (1-w)`, follow
[Threshold preparation and endpoint expansion](README.md#threshold-preparation-and-endpoint-expansion)
in the main README. It explains the branch preparation, shared expansion
routines, reuse of repeated structures and required validation. Follow both
the real contribution's threshold preparation in `3-Reduce.wl` and its use of
`ThresholdSplit` and `EndpointExpand` in the assembly steps. Keep the separate
virtual contribution and combine it with the real term and counterterms before
checking pole cancellation. Use the linked LeptonQuark files as the smallest
complete example. Preserve the distinction between the trace dimension and the
shifted integration dimension.

In new step files, load the FeynFacet package in a cell of its own, before
any cell that calls its functions such as `LoadFeynCalc[]`. Mathematica reads
a whole cell before evaluating it, so an unqualified call in the same cell
can create a `Global` symbol in a fresh kernel that shadows the package symbol.

Before a check that searches an expression for forbidden terms, rewrite
logarithms under the physical-region assumptions. For example, with
`s > 0` and `0 < v < 1`, rewrite `Log[s - s v]` as
`Log[s] + Log[1 - v]`. Use `ExactDifference` to compare expressions; it does
not rewrite the expression being exported.

After changing distribution handling, run
`wolframscript -file Checks/test_distributions.wls` and compare affected
coefficients exactly where tractable. The helper tests complement the
process checks; they do not replace them. Use fresh result locations when
validating changes, and identify any reused calculation inputs explicitly.

Keep local settings, execution logs, progress reports and internal reviews
outside the distributed source. After intentional release-file changes,
refresh and verify `SHA256SUMS` using `Checks/file_checksums.py`.
