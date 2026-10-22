# Installing the toolchain

Most users do not need to follow these steps by hand: an AI coding agent
(Codex or Claude Code) can install everything from them; see the Quick start
in [README.md](../README.md#quick-start).  They document each installation
step for the agent and for anyone installing manually.

Everything below is installed once, outside this repository.  The notebooks
find the programs automatically on the usual search path (`PATH`,
`~/.local/bin`, `/usr/local/bin`, `/opt/homebrew/bin`); anything else is
declared in `Paths.wl`.  The calculations were run on Linux x86_64 with the
versions in brackets.

For large hadronic calculations, see [memory and parallel runs](RESOURCES.md)
before choosing the number of simultaneous jobs.

## Mathematica (15.0.1)

A normal desktop installation.  In a fresh notebook, evaluate
`FileNameJoin[{$UserBaseDirectory, "Applications"}]`: this is the directory
where FeynCalc and Tracer go below.

## FeynCalc (10.2.1) with FeynArts (3.12)

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

Accept the FeynArts patch prompt.  Then quit the kernel and check

```wl
$LoadAddOns = {"FeynArts"}; Needs["FeynCalc`"];
{$FeynCalcVersion, FeynArts`$FeynArtsVersion}
```

## Tracer (1.1.1)

Download `Tracer.m` from the Wolfram Library Archive
(https://library.wolfram.com/infocenter/MathSource/2987/) and put it into
the `Applications` directory found above.  `FindFile["Tracer.m"]` must then
return its path.  The package manual, `Tracer.ps`, is distributed with it
and explains the notation used in the Tracer notebooks.

## FORM (5.0.1) and the colour library

Download the binary release from https://github.com/form-dev/form/releases
(or build from source), and put `form` and `tform` into `~/.local/bin`.
Download the colour package `color.h` from
https://www.nikhef.nl/~form/maindir/packages/color/color.html into
`~/.local/share/form/color/color.h`.  `form -v` must print the version.

## Kira (3.1) and Fermat (7.9b)

Kira: https://kira.hepforge.org/ provides a static Linux executable;
install it as `~/.local/bin/kira`. On Apple silicon use the
[pinned patched build](INSTALL-macOS.md#kira-built-from-source): unpatched
Kira 3.1 release builds can silently produce incorrect reductions; see the
[correctness notice](KIRA-CORRECTNESS.md). Fermat: download from https://home.bway.net/lewis/ and
extract; the executable is `Ferl7/fer64` on Linux, `Ferm7i/fer64` on Intel
Mac, `Ferm7a/fer64` on Apple silicon.  Kira is told where Fermat is through
`Paths.wl`, see below.

## SubTropica (1.2.10) and polymake (4.4)

In Mathematica:

```wl
PacletInstall["https://github.com/SubTropica/SubTropica/releases/download/v1.2.10/SubTropica-1.2.10.paclet"]
```

HyperIntica, the analytic integrator used here, is bundled with it.
polymake: `brew install polymake` (macOS, after `brew tap apaffenholz/polymake`)
or `sudo apt-get install polymake` (Debian/Ubuntu).

## Paths.wl

Copy `Paths.example.wl` to `Paths.wl` in the repository root and set the
Fermat executable:

```wl
$ToolPaths["Fermat"] = "/path/to/Fermat/Ferl7/fer64";
```

Add `$ToolPaths["Kira"] = ...`, `$ToolPaths["FORM"] = ...` and so on only if
a program is not found automatically.  `Paths.wl` is private to your
computer and is not part of the repository.

## Check

```wl
Get["/path/to/FeynFacet/FeynFacet/FeynFacet.wl"];   (* the clone of the repository *)
CheckInstallation[]
```

Every row must show a location (and a version where applicable).  Then test
that the tools work together, with small real calculations:

```sh
python3 feynfacet.py check --all-routes
```

and run the first process headlessly:

```sh
python Checks/run_process.py Processes/LeptonQuark --route FeynCalc
```
