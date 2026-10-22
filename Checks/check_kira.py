#!/usr/bin/env python3
"""Standalone exact Kira regression for release/CI (requires SymPy).

The normal feynfacet.py check uses Wolfram and does not require SymPy.
This script shares its fixtures and expected coefficients with that check.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from sympy import Tuple, cancel
from sympy.parsing.mathematica import parse_mathematica

FIXTURES = Path(__file__).resolve().parent / "kira"


def verify(text, spec):
    # Kira emits a newline before each unary '+'. Normalize those line breaks
    # before the Mathematica parser treats them as expression separators.
    try:
        parsed = parse_mathematica(" ".join(text.split()))
    except (IndexError, SyntaxError, TypeError) as error:
        raise ValueError("Malformed Kira output") from error
    if not isinstance(parsed, Tuple):
        raise ValueError("Kira result is not a list of rules")
    rules = {}
    for rule in parsed:
        if str(rule.func) != "Rule" or len(rule.args) != 2:
            raise ValueError("Malformed Kira rule")
        lhs, rhs = rule.args
        if lhs in rules:
            raise ValueError(f"Duplicate rule: {lhs}")
        rules[lhs] = rhs
    master = parse_mathematica(spec["master"])
    for case in spec["cases"]:
        target = parse_mathematica(case["target"])
        expected = parse_mathematica(case["coefficient"]) * master
        if target not in rules:
            raise ValueError(f"Missing rule: {target}")
        difference = cancel(rules[target] - expected)
        if difference != 0:
            raise ValueError(f"Incorrect reduction for {target}: difference {difference}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kira", required=True, type=Path)
    parser.add_argument("--fermat", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path,
                        help="new directory; existing directories are never overwritten")
    parser.add_argument("--timeout", default=60, type=int, help="seconds per reduction")
    args = parser.parse_args()
    kira, fermat = args.kira.resolve(), args.fermat.resolve()
    for path in (kira, fermat):
        if not path.is_file() or not os.access(path, os.X_OK):
            parser.error(f"Not an executable: {path}")
    args.workdir.mkdir(parents=True, exist_ok=False)
    specs = json.loads((FIXTURES / "expectations.json").read_text())
    report = {"executable": str(kira), "fermat": str(fermat), "tests": [], "passed": False}
    try:
        for spec in specs:
            work = args.workdir / spec["directory"]
            shutil.copytree(FIXTURES / spec["directory"], work)
            (work / "targets").write_text("\n".join(c["target"] for c in spec["cases"]) + "\n")
            start = time.monotonic()
            with (work / "kira.log").open("w") as log:
                process = subprocess.Popen([str(kira), "jobs.yaml"], cwd=work,
                                           env={**os.environ, "FERMATPATH": str(fermat)},
                                           stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    process.wait(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
                    raise
                if process.returncode:
                    raise subprocess.CalledProcessError(process.returncode, str(kira))
            result = work / "results" / spec["family"] / "kira_targets.m"
            verify(result.read_text(), spec)
            row = {"family": spec["family"], "exactIdentities": len(spec["cases"]),
                   "seconds": round(time.monotonic() - start, 2)}
            report["tests"].append(row)
            print("PASS", row, flush=True)
        report["passed"] = True
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        report["error"] = str(error)
        print("FAIL:", error, flush=True)
    finally:
        (args.workdir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
