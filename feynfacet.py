#!/usr/bin/env python3
"""Check that the FeynFacet toolchain works on this machine.

    python3 feynfacet.py check               # tools needed by the FeynCalc route
    python3 feynfacet.py check --all-routes  # also FORM and Tracer

Each test runs Checks/check_toolchain.wls in a fresh Wolfram kernel and performs
a small real calculation: a FeynCalc trace and FeynArts topologies, a Kira
reduction through Fermat, a native polymake calculation, SubTropica integrals
including a one-loop box that needs polymake, and with --all-routes a FORM
trace and a Tracer trace. Nothing is installed or changed outside Results/.
Each run works in its own folder Results/toolchain-check/run-<time>-<process>-<id>/,
so that checks started at the same time do not interfere; the ten most recently
completed runs are kept. Unfinished runs are preserved for diagnosis. The report is
printed, written to that folder, and copied to
Results/toolchain-check/report.json for the most recently finished run.
The exit status is 0 only when every required test passes.
"""

import argparse
import fcntl
import glob
import json
import os
import platform
import signal
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "Checks" / "check_toolchain.wls"
CHECKS = ROOT / "Results" / "toolchain-check"
KEEP_RUNS = 10

TESTS = [  # name, required for the FeynCalc route, required with --all-routes
    ("kernel", True, True),
    ("tools", True, True),
    ("feyncalc", True, True),
    ("kira", True, True),
    ("polymake", True, True),
    ("subtropica", True, True),
    ("form", False, True),
    ("tracer", False, True),
]
TOOLS_FEYNCALC = ["FeynCalc", "SubTropica", "Kira", "Fermat", "Polymake"]
TOOLS_ALL = TOOLS_FEYNCALC + ["FORM", "FORMColor", "Tracer"]

HINTS = {
    "kernel": "The Wolfram kernel did not start or could not write its files. Check the Mathematica "
              "licence, give the path with --wolfram, and when an agent runs this check, allow it to run "
              "Mathematica outside a restrictive sandbox: the kernel writes to its user and temporary directories.",
    "tools": "Install the missing programs (docs/INSTALL.md, docs/INSTALL-macOS.md) or give their locations in "
             "Paths.wl (see Paths.example.wl).",
    "feyncalc": "Install FeynCalc with FeynArts (docs/INSTALL.md). If it is installed outside the standard "
                "locations, load it from Paths.wl with PacletDirectoryLoad or AppendTo[$Path, ...].",
    "kira": "Read the kira.log named above. Common causes: a Fermat binary for another platform (Ferl7: Linux, "
            "Ferm7i: Intel Mac, Ferm7a: Apple silicon), a quarantined binary on macOS "
            "(xattr -d com.apple.quarantine .../fer64), or a sandbox that blocks Fermat's temporary files. "
            "On Apple silicon, unpatched release builds of Kira 3.1 can silently return wrong coefficients "
            "because of insufficient long-double precision (upstream MR !61). Use the pinned patched build "
            "in docs/INSTALL-macOS.md; see docs/KIRA-CORRECTNESS.md. A version banner alone is insufficient.",
    "polymake": "polymake starts but cannot compute; SubTropica needs it for multi-variable integrals. Read the "
                "errors above. On macOS with Homebrew a common cause is a Perl that no longer matches polymake's "
                "modules; see the polymake section of docs/INSTALL-macOS.md.",
    "subtropica": "Install SubTropica with PacletInstall (docs/INSTALL.md) and polymake. If SubTropica is "
                  "installed elsewhere, add PacletDirectoryLoad[\"/path/to/SubTropica\"] to Paths.wl. If a "
                  "polymake call failed, fix the polymake test first.",
    "form": "Install FORM and its color.h package (docs/INSTALL.md); the FORM log is in the folder named above.",
    "tracer": "Put Tracer.m into the Applications folder of $UserBaseDirectory, or give its path as "
              "$ToolPaths[\"Tracer\"] in Paths.wl.",
}


def find_wolfram(explicit):
    """wolframscript, or a WolframKernel that accepts -script."""
    if explicit:
        return [explicit] if explicit.endswith("wolframscript") else [explicit, "-script"]
    found = shutil.which("wolframscript")
    if found:
        return [found]
    patterns = ["/Applications/*.app/Contents/MacOS/wolframscript",
                "/usr/local/Wolfram/*/*/Executables/wolframscript",
                "/opt/Wolfram/*/*/Executables/wolframscript"]
    for pattern in patterns:
        for path in sorted(glob.glob(pattern), reverse=True):
            if os.access(path, os.X_OK):
                return [path]
    for pattern in ["/Applications/*.app/Contents/MacOS/WolframKernel",
                    "/usr/local/Wolfram/*/*/Executables/WolframKernel"]:
        for path in sorted(glob.glob(pattern), reverse=True):
            if os.access(path, os.X_OK):
                return [path, "-script"]
    return None


def run_test(wolfram, name, timeout, scratch):
    command = wolfram + (["-file"] if wolfram[-1].endswith("wolframscript") else []) + \
        [str(SCRIPT), name, str(scratch / name)]
    start = time.time()
    try:
        process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, start_new_session=True)
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            # Stop the kernel and its Kira/Fermat children, not just the launcher.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = process.communicate()
            return {"test": name, "status": "fail", "seconds": round(time.time() - start),
                    "detail": f"no result after {timeout} s (process group stopped)",
                    "output": (stdout + stderr)[-2000:]}
        output = stdout + stderr
    except OSError as error:
        return {"test": name, "status": "fail", "seconds": 0, "detail": str(error)}
    result = None
    for line in output.splitlines():
        if line.startswith("FEYNFACET_CHECK "):
            result = json.loads(line[len("FEYNFACET_CHECK "):])
    if result is None:
        result = {"test": name, "status": "fail", "detail": "the kernel stopped without a result",
                  "output": output[-2000:]}
    if process.returncode != 0 and result["status"] == "pass":
        result["status"] = "fail"
        result["detail"] = f"test printed pass but exited with status {process.returncode}"
    if name == "kira":
        for line in output.splitlines():
            if line.startswith("Checking Kira executable: "):
                result["executable"] = line.removeprefix("Checking Kira executable: ")
    result["seconds"] = round(time.time() - start)
    return result


def check(args):
    required = {name: (all_routes if args.all_routes else feyncalc) for name, feyncalc, all_routes in TESTS}
    CHECKS.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix=f"run-{time.strftime('%Y%m%d-%H%M%S')}-{os.getpid()}-", dir=CHECKS))
    report = {"system": {"os": platform.system(), "release": platform.release(),
                         "machine": platform.machine(), "python": platform.python_version()},
              "routes": "all" if args.all_routes else "FeynCalc", "tests": [],
              "reportPath": str(scratch / "report.json")}
    wolfram = find_wolfram(args.wolfram)
    report["wolfram"] = wolfram[0] if wolfram else None
    print(f"FeynFacet toolchain check ({report['routes']} route{'s' if args.all_routes else ''})")
    print(f"  system: {report['system']['os']} {report['system']['machine']}, Python {report['system']['python']}")
    print(f"  Wolfram: {report['wolfram'] or 'NOT FOUND'}")
    print(f"  Report: {scratch / 'report.json'}\n", flush=True)
    if not wolfram:
        report["tests"].append({"test": "kernel", "status": "fail", "detail": "no wolframscript or WolframKernel found"})
    else:
        for name, _, _ in TESTS:
            if name in ("form", "tracer") and not args.all_routes and not args.optional:
                continue
            result = run_test(wolfram, name, args.timeout, scratch)
            if name == "tools" and result["status"] == "pass":
                needed = TOOLS_ALL if args.all_routes else TOOLS_FEYNCALC
                missing = [tool for tool in needed if not result["detail"].get(tool)]
                if missing:
                    result["status"] = "fail"
                    result["missing"] = missing
            report["tests"].append(result)
            print(f"  {result['status'].upper():4}  {name:11} {result['seconds']:>4} s  "
                  f"{summary(result)}{'' if required[name] else '  (optional)'}")
            if name == "kernel" and result["status"] != "pass":
                break
    failures = [r for r in report["tests"] if r["status"] != "pass" and required.get(r["test"], True)]
    report["ready"] = not failures and bool(wolfram)
    for result in report["tests"]:
        if result["status"] != "pass":
            result["hint"] = HINTS.get(result["test"], "")
    print()
    for result in failures:
        print(f"{result['test']}: {result.get('hint', '')}")
        detail = result.get("detail")
        print("  detail: " + (json.dumps(detail) if isinstance(detail, dict) else str(detail)))
        if result.get("executable"):
            print("  Kira executable: " + result["executable"])
        log = Path(detail["log"]) if isinstance(detail, dict) and detail.get("log") else scratch / "kira" / "triangle" / "kira.log"
        if result["test"] == "kira" and log.exists():
            result["log"] = log.read_text(errors="replace")[-4000:]
            print("  end of kira.log:\n    " + "\n    ".join(result["log"].strip().splitlines()[-12:]))
        if result.get("output"):
            print("  last output:\n    " + "\n    ".join(result["output"].strip().splitlines()[-8:]))
    finish_run(report, scratch)
    print(("READY: every required test passed." if report["ready"] else
           "NOT READY: fix the failures above and run the check again.") +
          f"  Report: {scratch / 'report.json'}")
    return 0 if report["ready"] else 1


def atomic_report(path, text):
    """Publish complete JSON in one rename, including when a writer is interrupted."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", prefix=".report-",
                                         suffix=".json", dir=CHECKS, delete=False) as output:
            temporary = Path(output.name)
            output.write(text)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def finish_run(report, scratch):
    """Serialize publication and cleanup only; calculations never hold this lock."""
    with (CHECKS / ".report.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        text = json.dumps(report, indent=1)
        atomic_report(scratch / "report.json", text)
        atomic_report(CHECKS / "report.json", text)
        prune_runs(scratch)


def prune_runs(current):
    """Keep the newest completed runs. Called with the report lock held.

    A directory without a report may belong to a running or interrupted check;
    never infer that it is safe to remove from its age or process ID.
    """
    completed = []
    for path in CHECKS.glob("run-*"):
        if not path.is_dir() or path.is_symlink():
            continue
        try:
            completed.append(((path / "report.json").stat().st_mtime_ns, path))
        except FileNotFoundError:
            continue
    for _, path in sorted(completed)[:-KEEP_RUNS]:
        if path != current:
            shutil.rmtree(path, ignore_errors=True)


def summary(result):
    detail = result.get("detail")
    if result.get("missing"):
        return "missing: " + ", ".join(result["missing"])
    if isinstance(detail, dict):
        if result["test"] == "kernel":
            return detail.get("version", "")
        if result["test"] == "tools":
            return "all located"
        if result["test"] == "kira" and result["status"] == "pass":
            count = sum(row["exactIdentities"] for row in detail["reductions"])
            return f"{count} exact identities; {detail['executable']}"
        return "; ".join(f"{key} {value}" for key, value in detail.items())[:90]
    return str(detail)[:90]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    checker = commands.add_parser("check", help="run the functional toolchain checks")
    checker.add_argument("--all-routes", action="store_true", help="also require FORM and Tracer")
    checker.add_argument("--optional", action="store_true", help="also run the FORM and Tracer tests without requiring them")
    checker.add_argument("--wolfram", help="path to wolframscript or WolframKernel")
    checker.add_argument("--timeout", type=int, default=120, help="seconds allowed for each test (default 120; increase on slow machines)")
    args = parser.parse_args()
    return check(args)


if __name__ == "__main__":
    sys.exit(main())
