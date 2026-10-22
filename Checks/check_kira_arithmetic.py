#!/usr/bin/env python3
"""Fast source-level diagnostic; complements tests of the installed executable.

Compiles actual coeff_helper.cpp with a minimal release configuration.
Requires only Python and a C++ compiler, not Wolfram, Fermat or SymPy.
"""
import argparse
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="Kira source checkout")
    parser.add_argument("--workdir", required=True, type=Path, help="new scratch directory")
    args = parser.parse_args()
    source, work = args.source.resolve(), args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    (work / "pyred").mkdir()
    # No VALGRIND: test the release arithmetic path implicated in MR !61.
    (work / "pyred/config.h").write_text("#pragma once\n")
    # Other methods refer to Config::finished(); this field is irrelevant
    # to mod_mul. No arithmetic is mocked or copied into the test.
    state = work / "state.cpp"
    state.write_text('#include "pyred/defs.h"\nbool pyred::Config::s_finished = false;\n')
    executable = work / "mod-mul-check"
    subprocess.run(["c++", "-std=c++17", "-O2", "-I" + str(work),
                    "-I" + str(source / "src"),
                    str(source / "src/pyred/coeff_helper.cpp"), str(state),
                    str(Path(__file__).with_name("kira_mod_mul.cpp")), "-o", str(executable)],
                   check=True, timeout=60)
    return subprocess.run([str(executable)], timeout=15).returncode


if __name__ == "__main__":
    raise SystemExit(main())
