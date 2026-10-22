#!/usr/bin/env python3
"""Exercise report publication and overlapping checks without external physics tools.

Run with python3 Checks/test_toolchain_runner.py. The stand-in kernel tests the
Python runner only; real toolchain calculations are checked by feynfacet.py check.
"""

import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]

KERNEL = r'''
import json, os, pathlib, sys, time
name, folder = sys.argv[-2:]
folder = pathlib.Path(folder)
folder.mkdir(parents=True, exist_ok=True)
if name == "kernel" and os.environ.get("FF_TEST_GATE"):
    (folder / "started").write_text("running")
    deadline = time.monotonic() + 30
    parent = os.getppid()
    while not pathlib.Path(os.environ["FF_TEST_GATE"]).exists():
        if os.getppid() != parent:
            raise SystemExit("test runner stopped")
        if time.monotonic() > deadline:
            raise SystemExit("test gate timed out")
        time.sleep(.01)
    if not (folder / "started").exists():
        raise SystemExit("another check deleted this run")
details = "stand-in calculation"
if name == "kernel":
    details = {"version": "13.0" if os.environ.get("FF_TEST_OLD") else "15.0.1"}
elif name == "tools":
    details = {tool: "/stand-in/" + tool for tool in
               ("FeynCalc", "SubTropica", "Kira", "Fermat", "Polymake", "FORM", "FORMColor", "Tracer")}
elif name == "kira":
    details = {"executable": "/stand-in/kira", "reductions": [{"exactIdentities": 8}]}
print("FEYNFACET_CHECK " + json.dumps({"test": name, "status": "pass", "detail": details}))
'''


class ToolchainRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="feynfacet-runner-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        shutil.copy2(ROOT / "feynfacet.py", self.root / "feynfacet.py")
        self.kernel = self.root / "wolframscript"
        self.kernel.write_text("#!" + sys.executable + "\n" + KERNEL)
        self.kernel.chmod(0o755)
        self.checks = self.root / "Results" / "toolchain-check"

    def start(self, **environment):
        process = subprocess.Popen(
            [sys.executable, str(self.root / "feynfacet.py"), "check", "--all-routes",
             "--wolfram", str(self.kernel)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            env={**os.environ, **environment})
        self.addCleanup(self.stop, process)
        return process

    @staticmethod
    def stop(process):
        if process.poll() is None:
            process.kill()
        process.communicate()

    def finish(self, process):
        output, _ = process.communicate(timeout=15)
        self.assertEqual(process.returncode, 0, output)
        self.assertIn("READY: every required test passed.", output)
        report_line = next(line for line in output.splitlines() if line.startswith("  Report: "))
        report_path = Path(report_line.removeprefix("  Report: "))
        report = json.loads(report_path.read_text())
        self.assertTrue(report["ready"])
        self.assertEqual(report["reportPath"], str(report_path))
        return report_path, output

    def load_module(self):
        spec = importlib.util.spec_from_file_location("feynfacet_test", self.root / "feynfacet.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_overlapping_checks_preserve_active_run_and_complete_reports(self):
        gate = self.root / "release-slow-check"
        slow = self.start(FF_TEST_GATE=str(gate))
        deadline = time.monotonic() + 10
        while not list(self.checks.glob("run-*/kernel/started")):
            self.assertIsNone(slow.poll(), "slow check exited before reaching its gate")
            self.assertLess(time.monotonic(), deadline, "slow check did not start")
            time.sleep(.01)
        started = next(self.checks.glob("run-*/kernel/started"))
        active = started.parents[1]
        # These names sort after the active run. The former age-based pruning
        # deleted the active directory before either fast check could complete.
        for i in range(10):
            historical = self.checks / f"run-99999999-history-{i:02}"
            historical.mkdir()
            report = historical / "report.json"
            report.write_text(json.dumps({"ready": True, "history": i}))
            os.utime(report, (1000 + i, 1000 + i))
        latest = self.checks / "report.json"
        latest.write_text(json.dumps({"ready": True, "history": "latest"}))
        fast = [self.start(), self.start()]
        while any(process.poll() is None for process in fast):
            # A concurrent reader must always see complete JSON.
            self.assertTrue(json.loads(latest.read_text())["ready"])
            self.assertTrue(started.exists(), "cleanup removed an active run")
            self.assertLess(time.monotonic(), deadline + 10, "fast checks did not finish")
            time.sleep(.005)
        paths = [self.finish(process)[0] for process in fast]
        self.assertNotEqual(paths[0], paths[1])
        self.assertFalse((active / "report.json").exists())
        self.assertEqual(len(list(self.checks.glob("run-*/report.json"))), 10)
        gate.touch()
        slow_path, _ = self.finish(slow)
        self.assertEqual(slow_path.parent, active)
        self.assertEqual(json.loads(latest.read_text())["reportPath"], str(slow_path))
        self.assertEqual(len(list(self.checks.glob("run-*/report.json"))), 10)
        self.assertTrue(all(path.exists() for path in paths))

    def test_passing_older_kernel_is_not_rejected_by_version(self):
        report_path, output = self.finish(self.start(FF_TEST_OLD="1"))
        self.assertNotIn("WARNING:", output)
        report = json.loads(report_path.read_text())
        kernel = next(item for item in report["tests"] if item["test"] == "kernel")
        self.assertEqual(kernel["status"], "pass")
        self.assertEqual(kernel["detail"]["version"], "13.0")
        self.assertTrue(report["ready"])

    def test_repeated_runs_same_process_and_second_have_distinct_folders(self):
        module = self.load_module()
        args = type("Args", (), {"all_routes": True, "optional": False,
                                  "wolfram": str(self.kernel), "timeout": 10})()
        with patch.object(module.time, "strftime", return_value="20000101-000000"), \
                contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(module.check(args), 0)
            self.assertEqual(module.check(args), 0)
        reports = list(self.checks.glob("run-*/report.json"))
        self.assertEqual(len(reports), 2)
        self.assertNotEqual(reports[0].parent, reports[1].parent)

    def test_failed_publication_leaves_previous_report_complete(self):
        module = self.load_module()
        self.checks.mkdir(parents=True)
        latest = self.checks / "report.json"
        latest.write_text('{"previous": true}')
        with patch.object(module.os, "replace", side_effect=OSError("injected publication failure")):
            with self.assertRaises(OSError):
                module.atomic_report(latest, '{"new": true}')
        self.assertEqual(json.loads(latest.read_text()), {"previous": True})
        self.assertEqual(list(self.checks.glob(".report-*.json")), [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
