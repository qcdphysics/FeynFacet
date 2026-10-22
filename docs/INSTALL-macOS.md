# Installing the toolchain on macOS (Apple silicon)

Most users do not need to follow these steps by hand: an AI coding agent
(Codex or Claude Code) can install everything from them; see the Quick start
in [README.md](../README.md#quick-start).  They document each installation
step for the agent and for anyone installing manually.

Companion to `INSTALL.md`.  The calculations were developed on Linux; these
are the corresponding steps for an Apple-silicon Mac with Homebrew.
The Kira recipe below includes the required modular-arithmetic fix and tests
the build before installing it. Run the functional toolchain check before
starting a calculation; the processes' own physics checks remain required.

## Homebrew packages

```sh
brew install form                      # FORM and TFORM
brew tap apaffenholz/polymake && brew install polymake
brew install meson ninja pkg-config ginac cln yaml-cpp sqlite zlib gmp mpfr   # for building Kira
```

### polymake

SubTropica needs a working polymake for multi-variable integrals such as the
one-loop boxes. A polymake that prints its version may still be unable to
compute, so `python3 feynfacet.py check` runs a native polymake calculation
and a box integral, and every SubTropica integral in a calculation stops if a
polymake call fails. **If polymake passes the check, keep the installation
as it is**, whatever the macOS or Perl version.

Two problems with the Homebrew tap `apaffenholz/polymake` (polymake
4.15_2, October 2026) are known:

- **macOS 26 (Tahoe) and newer:** the formula is disabled ("does not
  build"), so no package can be installed there.
- **A newer Homebrew Perl:** the package's Perl modules were built for Perl
  5.42. After Homebrew upgrades Perl to 5.44, polymake fails with
  `Perl API version v5.42.0 of IO.c does not match v5.44.0`, and a rebuild
  from source refuses Perl 5.44 as well.

For an existing installation with the Perl mismatch, a repair using private
Perl 5.42.2 plus `SVG` and `XML::Writer` has passed native polymake and box
integral checks on one machine. It leaves Homebrew's Perl unchanged. The
available transcript does not contain the complete configure command or
launcher script, so this is **not yet a reproducible fresh-install recipe**.
Keep a working private runtime and launcher until their configuration has
been recorded and a replacement passes the same checks.

To recover the exact settings from that known working installation, read
its launcher and Perl's recorded configure arguments (adjust the prefix to
the existing installation):

```sh
cat "$HOME/.local/bin/polymake"
"$HOME/.local/opt/feynfacet-perl-5.42.2/bin/perl" -V:config_args -V:archname
```

The visible configure flags were `-Duseshrplib -Dusethreads`; they are not
claimed to be the complete command. A replacement recipe needs the full
record, a verified source archive, the compatible Perl modules, and a
launcher that selects that runtime. Select an existing working launcher
with `$ToolPaths["Polymake"]` in `Paths.wl`, then run
`python3 feynfacet.py check --all-routes`. Do not proceed with calculations
when the polymake or SubTropica test fails. On a fresh affected Mac, resolve
the tap's installation problem or use a separately validated compatible
installation first; a version banner alone is insufficient.

The platform restriction can be checked in the
[polymake Homebrew formula](https://github.com/apaffenholz/homebrew-polymake/blob/master/Formula/polymake.rb).

FORM's colour library:

```sh
mkdir -p ~/.local/share/form/color
curl -fL https://www.nikhef.nl/~form/maindir/packages/color/color.h -o ~/.local/share/form/color/color.h
```

## Mathematica packages

Install the tested FeynCalc revision
`881b6cec662665ca19c606be5b12fcda89501cdd`, whose package version is
10.2.1.  The commit pin below specifies the tested source precisely; it
is not a reference to the 10.2.1 release tag.  Evaluate the following
two inputs in separate cells of the same fresh kernel:

```wl
Import["https://raw.githubusercontent.com/FeynCalc/feyncalc/881b6cec662665ca19c606be5b12fcda89501cdd/install.m"]
```

```wl
InstallFeynCalc[AutoInstallFeynArts -> True,
  FeynCalcStableVersionLink -> "https://codeload.github.com/FeynCalc/feyncalc/zip/881b6cec662665ca19c606be5b12fcda89501cdd"]
```

Install SubTropica in a fresh kernel:

```wl
PacletInstall["https://github.com/SubTropica/SubTropica/releases/download/v1.2.10/SubTropica-1.2.10.paclet"]
```

Tracer: put `Tracer.m` into the directory printed by
`FileNameJoin[{$UserBaseDirectory, "Applications"}]`
(normally `~/Library/Wolfram/Applications`).

## Fermat

Download the Apple-silicon archive from https://home.bway.net/lewis/ and
extract it, for instance to `~/Applications/Fermat/`.  The executable is
`Ferm7a/fer64`.  macOS may refuse to run a downloaded binary the first time;
allow it in System Settings → Privacy & Security, or run
`xattr -d com.apple.quarantine ~/Applications/Fermat/Ferm7a/fer64`.

## Kira (built from source)

**Unpatched Kira 3.1 release builds on Apple silicon can silently produce
incorrect IBP reductions.** The floating-point modular multiplication assumes
more precision than this platform's `long double` provides. Do not wait for
a failed pole-cancellation check: use the patched build below. See the
[correctness notice](KIRA-CORRECTNESS.md), including guidance for earlier runs.

If Kira is already installed, run the strengthened `python3 feynfacet.py check`
first. Keep an existing validated patched build; installation is only needed
when the tool is missing or its build needs correction.

From the FeynFacet repository, with the packages above and a native Fermat
executable already installed:

```sh
python3 tools/install_kira.py \
  --workdir "$HOME/Library/Caches/feynfacet-kira-build-1" \
  --fermat "$HOME/Applications/Fermat/Ferm7a/fer64"
```

Adjust `--fermat` to its actual location (for example `~/.local/bin/fer64`).
The work directory must be new. The installer verifies the pinned Kira source
and bundled patch, pins FireFly and the FLINT source/build archives, and
forces bundled dependencies, including yaml-cpp. Its modular-arithmetic
probe is a Meson target using Kira's own compiler settings and dependencies.
After that probe and the two full upstream regressions pass, it installs under
`~/.local`. It does not use `sudo` or change shell/Wolfram configuration.
An existing `~/.local/bin/kira` is protected: explicitly pass `--replace` to
replace it, or choose a separate `--prefix /absolute/path` for evaluation.
The prefix also receives the bundled libraries and headers. The installer
checks the installed Kira/FireFly library references with `otool -L`
(`ldd` on Linux) and records their paths in `installation.json`.

This is a one-time build, not part of routine toolchain checking. It uses at
most six build workers by default (`--jobs N` overrides that) and a total
20-minute deadline (`--timeout SECONDS`). On failure it stops, retains
`install.log` and `installation.json` in the work directory, and does not
silently retry. Downloads, compilation and tests count toward the deadline.
If a dependency or network failure occurs, resolve the reported cause before
retrying in a new work directory. Do not rebuild an already validated Kira
just to run another calculation.

The build is **Kira 3.1 plus the reviewed MR !61 patch**, not an official
fixed Kira release. The exact revisions and patch checksum are recorded in
[the patch provenance](../tools/patches/README.md). A version banner alone
does not prove that a binary contains the fix. Rosetta has not been validated
by this recipe and is not a recommended workaround.

## Paths.wl

```wl
$ToolPaths["Fermat"] = FileNameJoin[{$HomeDirectory, "Applications", "Fermat", "Ferm7a", "fer64"}];
```

Homebrew's programs live in `/opt/homebrew/bin`, which the package searches;
Kira built as above is found in `~/.local/bin`.  Add `$ToolPaths["Kira"]`
etc. only if `CheckInstallation[]` reports something missing.

## Check

```wl
Get["/path/to/FeynFacet/FeynFacet/FeynFacet.wl"];   (* the clone of the repository *)
CheckInstallation[]
```

then test that the tools work together, with small real calculations
(a Kira reduction through Fermat, a native polymake calculation, SubTropica
integrals including a one-loop box, FeynCalc, FORM and Tracer traces):

```sh
python3 feynfacet.py check --all-routes
```

It reports each failure with a suggested fix. The Kira test checks eight
exact identities in two small standalone integral families and reports the
selected executable. It does not rebuild tools or run a physics process.
The per-test timeout is 120 seconds; increase `--timeout` only when needed.
For a custom installation prefix, set `$ToolPaths["Kira"]` in `Paths.wl`
to that prefix's `bin/kira` before checking.
Then run `Processes/LeptonQuark` as described in its
README.  Runtime depends on the process, toolchain and machine.
