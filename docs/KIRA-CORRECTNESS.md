# Kira correctness notice: Apple Silicon

Unpatched Kira 3.1 release builds using the floating-point `mod_mul` path
can silently produce wrong IBP reductions on Apple Silicon. The algorithm
requires at least 64 significant binary bits in `long double`; this macOS
platform provides 53. An executable can exit successfully while returning
wrong coefficients or incorrectly declaring integrals zero.

FeynFacet supplies the reviewed [upstream MR !61 fix](https://gitlab.com/kira-pyred/kira/-/merge_requests/61)
with pinned source revisions, a patch checksum and an installation procedure
in [INSTALL-macOS.md](INSTALL-macOS.md#kira-built-from-source).
This is a locally patched Kira release, not an assertion that an official
release already contains the fix. The precision guard also leaves platforms
with sufficient `long double` precision on their original arithmetic path.

## Earlier calculations

If a calculation used an affected binary, install and validate the corrected
Kira, regenerate its reductions in fresh result directories, and rerun the
dependent integration/assembly steps and their physics checks. Do not reuse
cached reductions from that binary. Preserve old results separately for
comparison and record which binary generated each set.

Successful pole cancellation alone does not certify every reduction.
Conversely, this notice does not imply that results generated with a known
patched binary or a separately validated unaffected build are wrong.

## What FeynFacet checks

- The macOS recipe applies the reviewed patch before building and runs the
  modular-arithmetic test and Kira's upstream regression tests before
  installing.
- The installer forces the bundled dependencies, including FLINT 2.8.4
  and yaml-cpp. Kira's checked FLINT source archive and Meson patch hashes
  are recorded in `installation.json`. `--force-fallback-for=flint` alone
  does not control FireFly's direct library search; the installer uses
  `--wrap-mode=forcefallback`.
- The arithmetic test is a target in Kira's own Meson build, using its
  configuration, compiler, `libpyred` and dependency objects. It does not
  reconstruct linker flags from a possibly different system yaml-cpp.
  After installation, `ldd` (Linux) or `otool -L` (macOS) checks Kira and
  FireFly against the installed libraries recorded by Meson. On macOS,
  relative install names are resolved using the loader location and
  `otool -l` run paths, including their search order. Static
  libraries are recorded explicitly instead of being mistaken for missing
  shared libraries. The loader output is retained in `installation.json`;
  an unexpected or missing shared library prevents a successful report.
- `python3 feynfacet.py check` compares the full rational coefficients of
  eight reduction identities, including explicit scaleless zeros, in a one-
  and a two-loop family, and reports the Kira executable it tested.
- `Checks/check_kira.py` runs the same exact identities without Mathematica
  (it needs SymPy); see [Checks/kira/README.md](../Checks/kira/README.md).

These checks detect the reproduced failure and several kinds of incomplete or
corrupted output, but no finite set of identities proves every reduction
correct: the physics checks of each process remain necessary.
