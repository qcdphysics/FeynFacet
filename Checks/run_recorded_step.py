#!/usr/bin/env python3
"""Run a notebook step and record every new or changed file in its output scope.

The Wolfram runner keeps its cell-level receipt. This wrapper adds the actual
process exit and an inventory that includes external programs and subkernels.
For a clean algebra run, use --require-empty and --set resume=no. Keep records
outside the source tree; they contain local paths and are not release files.
"""
import argparse
from contextlib import ExitStack
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def inventory(directory):
    """Hash contents, rather than using timestamps to infer changes."""
    result = {}
    for p in sorted(directory.rglob("*")):
        if ".recording-locks" in p.relative_to(directory).parts:
            continue
        if p.is_symlink():
            raise RuntimeError("output scope contains a symbolic link: " + str(p))
        if p.is_file():
            result[str(p.relative_to(directory))] = digest(p)
    return result


def sources(root):
    result = subprocess.run(["git", "-C", str(root), "ls-files", "-z"],
                            capture_output=True, check=False)
    if result.returncode == 0:
        paths = [root / os.fsdecode(s) for s in result.stdout.split(b"\0") if s]
    else:
        paths = []
    # Include local calculation dependencies, even when not tracked by Git.
    # Prune generated directories before descending: Results can contain many
    # thousands of pair products and Kira files, none of which are source files.
    # Tracked files above remain included, including tracked hidden files.
    # Ignore untracked hidden working directories during local discovery.
    excluded = {"Results", "__pycache__", "public"}
    extensions = {".wl", ".wls", ".m", ".frm", ".h", ".inc",
                  ".yaml", ".yml", ".py", ".json", ".mod", ".gen"}
    for directory, folders, filenames in os.walk(root):
        folders[:] = [name for name in folders
                      if name not in excluded and not name.startswith(".")]
        paths += [Path(directory) / name for name in filenames
                  if Path(name).suffix in extensions]
    if (root / "Paths.wl").is_file():
        paths.append(root / "Paths.wl")
    return {str(p.relative_to(root)): digest(p) for p in sorted(set(paths))}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("step", help="step file relative to the repository root")
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--scope", type=Path, required=True,
                    help="exclusive output directory relative to the root; repeat runs must not overlap")
    ap.add_argument("--record", type=Path, required=True,
                    help="new private directory for this execution's evidence")
    ap.add_argument("--bound-files", type=Path,
                    help="JSON map of immutable absolute file paths to SHA256")
    ap.add_argument("--wolfram", default="wolframscript")
    ap.add_argument("--require-empty", action="store_true")
    ap.add_argument("--timeout", type=float, help="optional child wall time in seconds")
    args = ap.parse_args()
    root, record = args.root.resolve(), args.record.resolve()
    step = (root / args.step).resolve()
    scope = (root / args.scope).resolve()
    step_relative = step.relative_to(root)
    scope_relative = scope.relative_to(root)
    if step_relative.parts[:1] != ("Processes",) or "Results" not in scope_relative.parts:
        ap.error("step must be in Processes/ and scope must be in its Results/")
    if scope_relative.parts[:2] != step_relative.parts[:2]:
        ap.error("step and output scope must belong to the same process")
    if root == record or root in record.parents:
        ap.error("--record must be outside the source repository")
    if any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*=[A-Za-z0-9_,.+-]+", item)
           for item in args.set):
        ap.error("each --set must be a distinct NAME=VALUE")
    choices = [item.split("=", 1) for item in args.set]
    if len({key for key, value in choices}) != len(choices):
        ap.error("duplicate override names")
    suffix = "-" + "-".join(value for key, value in choices) if choices else ""
    receipt = step.parent / "Results" / (step.stem + suffix + ".receipt.json")
    runner = root / "Checks/run_step.wls"
    if not step.is_file() or not runner.is_file():
        ap.error("step or Wolfram runner is missing")
    record.mkdir(parents=True, exist_ok=False)
    scope.mkdir(parents=True, exist_ok=True)
    lock_dir = step.parent / "Results/.recording-locks"
    lock_dir.mkdir(exist_ok=True)
    with ExitStack() as stack:
        receipt_lock_name = "receipt-" + hashlib.sha256(str(receipt).encode()).hexdigest() + ".lock"
        receipt_lock = stack.enter_context((lock_dir / receipt_lock_name).open("a+"))
        fcntl.flock(receipt_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        results_root = step.parent / "Results"
        hierarchy = [results_root]
        for part in scope.relative_to(results_root).parts:
            hierarchy.append(hierarchy[-1] / part)
        for directory in hierarchy:
            lock_name = hashlib.sha256(str(directory).encode()).hexdigest() + ".lock"
            lock = stack.enter_context((lock_dir / lock_name).open("a+"))
            mode = fcntl.LOCK_EX if directory == scope else fcntl.LOCK_SH
            fcntl.flock(lock, mode | fcntl.LOCK_NB)
        before = inventory(scope)
        if args.require_empty and before:
            raise RuntimeError("clean execution requires an empty output scope")
        if args.require_empty and dict(choices).get("resume") != "no":
            raise RuntimeError("clean algebra execution must explicitly set resume=no")
        previous_receipt = None
        if receipt.exists():
            previous_receipt = digest(receipt)
            shutil.copy2(receipt, record / "previous-receipt.json")
            receipt.unlink()
        source_before = sources(root)
        bound = json.loads(args.bound_files.read_text()) if args.bound_files else {}
        for name, expected in bound.items():
            if digest(name) != expected:
                raise RuntimeError("immutable input changed: " + name)
        command = [args.wolfram, "-file", str(runner), str(step)]
        env = dict(os.environ, FEYNFACET_SET=";".join(args.set))
        started = time.time()
        evidence = {"format": 1, "step": str(step), "root": str(root),
                    "scope": str(scope), "command": command, "overrides": choices,
                    "started_unix": started, "wrapper_sha256": digest(__file__),
                    "sources": source_before, "immutable_bound_files": bound,
                    "scope_before": before, "require_empty": args.require_empty}
        evidence["previous_receipt_sha256"] = previous_receipt
        write_json(record / "launch.json", evidence)
        log = stack.enter_context((record / "stdout.log").open("xb"))
        proc = subprocess.Popen(command, cwd=step.parent, env=env, stdout=log,
                                stderr=subprocess.STDOUT, start_new_session=True,
                                close_fds=True)
        write_json(record / "running.json", {"pid": proc.pid,
                   "started_unix": started, "owner_pid": os.getpid()})
        timed_out = False
        try:
            code = proc.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                code = proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                code = proc.wait()
        # Write observed exit before subsequent validation can fail.
        write_json(record / "actual-exit.json", {"pid": proc.pid,
                   "exit_code": code, "exit_status_observed": True,
                   "timed_out": timed_out, "ended_unix": time.time()})
        after = inventory(scope)
        new_or_changed = {name: value for name, value in after.items()
                          if before.get(name) != value}
        removed = sorted(set(before) - set(after))
        errors = []
        if args.require_empty and ("Algebra.wl" not in new_or_changed or not any(
                name.startswith("pairs/") and name.endswith((".m", ".wl"))
                for name in new_or_changed)):
            errors.append("clean algebra run did not produce Algebra.wl and pair products")
        receipt_data = None
        if receipt.is_file() and receipt.stat().st_mtime >= started - 1:
            try:
                receipt_data = json.loads(receipt.read_text())
            except ValueError:
                errors.append("invalid Wolfram receipt JSON")
        else:
            errors.append("missing fresh Wolfram receipt")
        if receipt_data is not None:
            if receipt_data.get("completed") is not True:
                errors.append("Wolfram receipt incomplete")
            if Path(receipt_data.get("file", "")).resolve() != step:
                errors.append("Wolfram receipt names another step")
            if receipt_data.get("overrides") != choices:
                errors.append("Wolfram receipt overrides differ")
            cells = receipt_data.get("cells", [])
            if not cells or any(cell.get("failed") is not False or
                                cell.get("cell") != index
                                for index, cell in enumerate(cells, 1)):
                errors.append("Wolfram receipt needs nonempty sequential successful cells")
            provenance = receipt_data.get("provenance", {})
            if not provenance.get("sources"):
                errors.append("Wolfram receipt has no source bindings")
            for kind in ("sources", "inputs", "outputs"):
                for name, expected in provenance.get(kind, {}).items():
                    path = root / name
                    if not path.is_file() or digest(path) != expected:
                        errors.append("Wolfram " + kind + " hash differs: " + name)
        if sources(root) != source_before:
            errors.append("source files changed during execution")
        for name, expected in bound.items():
            if not Path(name).is_file() or digest(name) != expected:
                errors.append("immutable dependency changed during execution: " + name)
        if code != 0 or timed_out:
            errors.append("child did not exit successfully")
        evidence.update({"actual_exit_code": code, "exit_status_observed": True,
                         "ended_unix": time.time(), "scope_after": after,
                         "outputs": new_or_changed, "removed": removed,
                         "receipt": str(receipt),
                         "receipt_sha256": digest(receipt) if receipt.is_file() else None,
                         "source_and_bound_files_unchanged": not any(
                             "changed during execution" in e for e in errors),
                         "errors": errors, "passed": not errors})
        write_json(record / "execution.json", evidence)
        print(json.dumps({"record": str(record), "passed": not errors,
                          "exit_code": code, "outputs": len(new_or_changed),
                          "errors": errors}), flush=True)
        return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
