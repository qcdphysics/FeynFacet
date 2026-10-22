#!/usr/bin/env python3
"""Dependency-selection regressions; no installed tools or libraries are changed."""
import importlib.util
import contextlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("install_kira", ROOT / "tools/install_kira.py")
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


class InstallerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="kira-installer-test-")
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.libraries = {name: [{"type": "shared library", "installed": str(self.work / f"lib{name}.so")}]
                          for name in ("flint", "yaml-cpp", "firefly")}

    def test_pin_requires_expected_revision(self):
        wrap = self.work / "dependency.wrap"
        wrap.write_text("[wrap-git]\nrevision = head")
        installer.pin_wrap(wrap, "head", "1234")
        self.assertEqual(wrap.read_text(), "[wrap-git]\nrevision = 1234")
        with self.assertRaises(ValueError):
            installer.pin_wrap(wrap, "head", "5678")

    def test_flint_archive_pin_rejects_changed_source(self):
        (self.work / "subprojects").mkdir()
        wrap = self.work / "subprojects/flint2.wrap"
        wrap.write_text("[wrap-file]\ndirectory = flint-2.8.4\nsource_hash = "
                        + installer.FLINT_SOURCE_SHA256 + "\npatch_hash = " + installer.FLINT_PATCH_SHA256)
        self.assertEqual(installer.check_flint_wrap(self.work)["directory"], "flint-2.8.4")
        wrap.write_text(wrap.read_text().replace(installer.FLINT_SOURCE_SHA256, "0" * 64))
        with self.assertRaises(ValueError):
            installer.check_flint_wrap(self.work)

    def test_linux_records_exact_libraries(self):
        output = "\n".join(f"lib{name}.so => {items[0]['installed']} (0x1234)"
                           for name, items in self.libraries.items())
        self.assertEqual(set(installer.check_library_references(output, self.libraries, "Linux")),
                         set(self.libraries))

    def test_prefix_with_spaces(self):
        path = str(self.work / "prefix with spaces/libflint.so")
        self.libraries["flint"][0]["installed"] = path
        for system, line in (("Linux", f"libflint.so => {path} (0x1234)"),
                             ("Darwin", f"\t{path} (compatibility version 0.0.0, current version 0.0.0)")):
            self.assertEqual(installer.check_library_references(line, self.libraries, system), {"flint": path})

    def test_linux_rejects_system_flint_and_yaml(self):
        for name in ("flint", "yaml-cpp"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                installer.check_library_references(f"lib{name}.so => /usr/lib/lib{name}.so (0x123)",
                                                   self.libraries, "Linux")

    def test_missing_library_is_a_failure(self):
        with self.assertRaises(ValueError):
            installer.check_library_references("libgmp.so => not found", self.libraries, "Linux")

    def test_macos_references(self):
        path = self.libraries["flint"][0]["installed"]
        output = f"kira:\n\t{path} (compatibility version 0.0.0, current version 0.0.0)"
        self.assertEqual(installer.check_library_references(output, self.libraries, "Darwin"),
                         {"flint": path})
        with self.assertRaises(ValueError):
            installer.check_library_references(output.replace(path, "/opt/homebrew/lib/libflint.dylib"),
                                               self.libraries, "Darwin")

    def mac_fixture(self):
        prefix = self.work / "Mac prefix"
        (prefix / "lib").mkdir(parents=True)
        for name in self.libraries:
            path = prefix / "lib" / f"lib{name}.dylib"
            path.touch()
            self.libraries[name][0]["installed"] = str(path)
        return prefix

    def test_macos_rpath_and_loader_origins(self):
        prefix = self.mac_fixture()
        executable = prefix / "bin/kira"
        firefly = prefix / "lib/libfirefly.dylib"
        version = " (compatibility version 0.0.0, current version 0.0.0)"
        def run(command, capture=False):
            if command[1] == "-l":
                if Path(command[2]) == executable:
                    return "Load command 12\n          cmd LC_RPATH\n      cmdsize 48\n         path @executable_path/../lib (offset 12)\n"
                return "Load command 8\n          cmd LC_RPATH\n      cmdsize 32\n         path @loader_path (offset 12)\n"
            if Path(command[2]) == executable:
                return "kira:\n\t@rpath/libfirefly.dylib" + version + "\n\t@executable_path/../lib/libyaml-cpp.dylib" + version
            self.assertEqual(Path(command[2]), firefly)
            return "firefly:\n\t@loader_path/libflint.dylib" + version
        report = installer.inspect_linkage(executable, self.libraries, "Darwin", run)
        self.assertEqual(report["references"], {name: products[0]["installed"]
                                               for name, products in self.libraries.items()})
        self.assertEqual(len(report["loadCommands"]), 2)

    def test_macos_build_linkage_uses_build_library_paths(self):
        prefix = self.mac_fixture()
        build = self.work / "build with spaces"
        executable = build / "src/kira"
        for name, products in self.libraries.items():
            path = build / "subprojects" / name / f"lib{name}.dylib"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
            products[0]["build"] = str(path)
        firefly = Path(self.libraries["firefly"][0]["build"])
        version = " (compatibility version 0.0.0, current version 0.0.0)"
        calls = []

        def run(command, capture=False):
            calls.append(command)
            binary = Path(command[2])
            if command[1] == "-l":
                paths = (["@executable_path/../subprojects/firefly",
                          "@executable_path/../subprojects/yaml-cpp"]
                         if binary == executable else ["@loader_path/../flint"])
                return "".join("Load command 12\n          cmd LC_RPATH\n      cmdsize 80\n"
                               f"         path {path} (offset 12)\n" for path in paths)
            if binary == executable:
                return ("kira:\n\t@rpath/libfirefly.dylib" + version
                        + "\n\t@rpath/libyaml-cpp.dylib" + version)
            self.assertEqual(binary, firefly)
            return "firefly:\n\t@rpath/libflint.dylib" + version

        report = installer.inspect_linkage(executable, self.libraries, "Darwin", run,
                                           location="build")
        self.assertEqual(report["references"], {name: products[0]["build"]
                                               for name, products in self.libraries.items()})
        self.assertTrue(all(str(prefix) not in str(command[2]) for command in calls))

    def test_macos_rpath_search_order_rejects_external_library(self):
        prefix = self.mac_fixture()
        external = self.work / "external"
        external.mkdir()
        (external / "libflint.dylib").touch()
        output = "\t@rpath/libflint.dylib (compatibility version 0.0.0, current version 0.0.0)"
        with self.assertRaises(ValueError):
            installer.check_library_references(output, self.libraries, "Darwin",
                prefix / "bin/kira", prefix / "bin/kira", [external, prefix / "lib"])
        found = installer.check_library_references(output, self.libraries, "Darwin",
                prefix / "bin/kira", prefix / "bin/kira", [prefix / "lib", external])
        self.assertEqual(found["flint"], str(prefix / "lib/libflint.dylib"))

    def test_macos_dylib_identity_does_not_prove_kira_loads_it(self):
        prefix = self.mac_fixture()
        executable = prefix / "bin/kira"
        version = " (compatibility version 0.0.0, current version 0.0.0)"
        def run(command, capture=False):
            if command[1] == "-l":
                return ""
            if Path(command[2]) == executable:
                return "kira:\n\t" + str(prefix / "lib/libyaml-cpp.dylib") + version
            # otool -L includes FireFly's LC_ID_DYLIB as well as its loads.
            return ("firefly:\n\t" + str(prefix / "lib/libfirefly.dylib") + version
                    + "\n\t" + str(prefix / "lib/libflint.dylib") + version)
        with self.assertRaisesRegex(ValueError, "Kira does not reference.*FireFly"):
            installer.inspect_linkage(executable, self.libraries, "Darwin", run)

    def test_macos_rpath_needs_an_existing_file(self):
        prefix = self.mac_fixture()
        with self.assertRaises(ValueError):
            installer.check_library_references(
                "\t@rpath/libflint.dylib (compatibility version 0.0.0, current version 0.0.0)",
                self.libraries, "Darwin", prefix / "bin/kira", prefix / "bin/kira", [self.work / "missing"])

    def test_static_libraries_are_explicit(self):
        targets, installed = [], {}
        for name in self.libraries:
            filename = str(self.work / f"build/lib{name}.a")
            targets.append({"name": name, "type": "static library", "subproject": name,
                            "filename": [filename]})
            installed[filename] = str(self.work / f"lib{name}.a")
        libraries = installer.bundled_libraries(targets, installed)
        calls = []
        def run(command, capture=False):
            calls.append(command)
            return "libc.so => /usr/lib/libc.so (0x123)"
        report = installer.inspect_linkage(self.work / "kira", libraries, "Linux", run)
        self.assertEqual(report["references"], {})
        self.assertEqual(len(calls), 1)

    def test_external_library_cannot_count_as_bundled(self):
        with self.assertRaises(ValueError):
            installer.bundled_libraries([], {})

    def test_missing_shared_library_is_not_assumed_static(self):
        with self.assertRaises(ValueError):
            installer.inspect_linkage(self.work / "kira", self.libraries, "Linux",
                                      lambda *args, **kwargs: "libc.so => /usr/lib/libc.so (0x123)")

    def run_install_fixture(self, name, wrong_library=None, old_prefix=False):
        """Run main and the real linkage checks; replace only external commands.

        The fake install deliberately overwrites sentinel files. A check moved
        after installation will therefore fail the preservation assertion.
        """
        fixture = self.work / name
        prefix = fixture / "validated installation"
        source = fixture / "work/kira"
        build = source / "builddir"
        executable = build / "src/kira"
        fermat = fixture / "fermat"
        (prefix / "bin").mkdir(parents=True)
        (prefix / "lib").mkdir()
        existing = [prefix / "bin/kira"] + [prefix / "lib" / f"lib{lib}.so"
                                              for lib in ("flint", "yaml-cpp", "firefly")]
        for path in existing:
            path.write_bytes(("previously validated " + path.name).encode())
        before = {path: path.read_bytes() for path in existing}
        fermat.write_text("#!/bin/sh\nexit 0\n")
        fermat.chmod(0o755)
        targets = [{"name": "kira", "type": "executable", "filename": [str(executable)]}]
        installed = {str(executable): str(prefix / "bin/kira")}
        for lib in ("flint", "yaml-cpp", "firefly"):
            path = build / "subprojects" / lib / f"lib{lib}.so"
            targets.append({"name": lib, "type": "shared library", "subproject": lib,
                            "filename": [str(path)]})
            installed[str(path)] = str(prefix / "lib" / path.name)
        calls = []

        def external(command, **kwargs):
            command = tuple(map(str, command))
            calls.append(command)
            output = ""
            if command[:2] == ("git", "clone"):
                (source / "subprojects/firefly").mkdir(parents=True)
                (source / "tests").mkdir()
                (source / "tests/meson.build").write_text("")
                (source / "subprojects/firefly.wrap").write_text("[wrap-git]\nrevision = kira-2\n")
                (source / "subprojects/flint2.wrap").write_text(
                    "[wrap-file]\ndirectory = flint-2.8.4\nsource_hash = "
                    + installer.FLINT_SOURCE_SHA256 + "\npatch_hash = " + installer.FLINT_PATCH_SHA256)
            elif command[:2] == ("git", "rev-parse"):
                output = (installer.FIREFLY_COMMIT if Path(kwargs["cwd"]).name == "firefly"
                          else installer.KIRA_COMMIT)
            elif command[:2] == ("git", "apply"):
                pass
            elif command[:2] == ("meson", "setup"):
                build.mkdir()
            elif command[:3] == ("meson", "introspect", "--targets"):
                output = json.dumps(targets)
            elif command[:3] == ("meson", "introspect", "--installed"):
                output = json.dumps(installed)
            elif command[0] == "ninja":
                for target in targets:
                    path = Path(target["filename"][0])
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(("new build " + path.name).encode())
            elif command[:2] == ("meson", "test"):
                pass
            elif command[:2] == ("meson", "install"):
                for built, destination in installed.items():
                    Path(destination).write_bytes(Path(built).read_bytes())
            elif command[0] == "ldd":
                binary = Path(command[1])
                in_build = binary == executable or binary.is_relative_to(build)
                libraries = ("flint",) if binary.name == "libfirefly.so" else ("flint", "yaml-cpp", "firefly")
                lines = []
                for lib in libraries:
                    expected = (build / "subprojects" / lib / f"lib{lib}.so" if in_build
                                else prefix / "lib" / f"lib{lib}.so")
                    if in_build and lib == wrong_library:
                        expected = ((prefix / "lib") if old_prefix else Path("/usr/lib")) / f"lib{lib}.so"
                    lines.append(f"lib{lib}.so => {expected} (0x1234)")
                output = "\n".join(lines)
            elif command == (str(prefix / "bin/kira"), "--version"):
                output = "Kira 3.1 fixture"
            else:
                raise AssertionError(f"Unexpected external command: {command}")

            class Completed:
                returncode = 0

                def communicate(self, timeout=None):
                    return output, None

            return Completed()

        argv = ["install_kira.py", "--workdir", str(fixture / "work"), "--prefix", str(prefix),
                "--fermat", str(fermat), "--replace"]
        with patch.object(sys, "argv", argv), patch.object(installer.platform, "system", return_value="Linux"), \
                patch.object(installer.platform, "platform", return_value="Linux fixture"), \
                patch.object(installer.shutil, "which", side_effect=lambda tool: "/fixture/bin/" + tool), \
                patch.object(installer.subprocess, "Popen", side_effect=external), \
                contextlib.redirect_stdout(io.StringIO()):
            status = installer.main()
        report = json.loads((fixture / "work/installation.json").read_text())
        return status, report, calls, before, executable, prefix

    def test_wrong_build_linkage_preserves_existing_installation(self):
        for library in ("flint", "yaml-cpp"):
            for old_prefix in (False, True):
                with self.subTest(library=library, old_prefix=old_prefix):
                    status, report, calls, before, executable, _ = self.run_install_fixture(
                        f"bad-{library}-{old_prefix}", wrong_library=library, old_prefix=old_prefix)
                    self.assertEqual(status, 1)
                    self.assertFalse(report["completed"])
                    self.assertIn(f"unexpected {library}", report["error"])
                    self.assertIn(("ldd", str(executable)), calls)
                    self.assertFalse(any(command[:2] == ("meson", "install") for command in calls))
                    self.assertEqual({path: path.read_bytes() for path in before}, before)

    def test_success_checks_build_before_install_and_installed_paths_after(self):
        status, report, calls, before, executable, prefix = self.run_install_fixture("success")
        self.assertEqual(status, 0)
        self.assertTrue(report["completed"])
        install_at = next(i for i, command in enumerate(calls) if command[:2] == ("meson", "install"))
        build_firefly = executable.parent.parent / "subprojects/firefly/libfirefly.so"
        for path in (executable, build_firefly):
            self.assertLess(calls.index(("ldd", str(path))), install_at)
        for path in (prefix / "bin/kira", prefix / "lib/libfirefly.so"):
            self.assertGreater(calls.index(("ldd", str(path))), install_at)
        self.assertEqual(set(report["buildLinkage"]["references"]), {"flint", "yaml-cpp", "firefly"})
        self.assertTrue(all(str(executable.parent.parent) in path
                            for path in report["buildLinkage"]["references"].values()))
        self.assertTrue(all(str(prefix / "lib") in path for path in report["linkage"]["references"].values()))
        self.assertTrue(all(path.read_bytes() != content for path, content in before.items()))

    @unittest.skipUnless(shutil.which("meson") and shutil.which("ninja") and shutil.which("c++"),
                         "Meson, Ninja and C++ are needed for the compiled dependency test")
    def test_probe_uses_meson_configuration_and_bundled_dependency(self):
        # The probe calls pyRed from a static archive that needs the selected YAML
        # library. No pkg-config yaml-cpp is supplied: rebuilding flags independently
        # would fail. The reference arithmetic remains Checks/kira_mod_mul.cpp.
        source = self.work / "source with spaces"
        (source / "tests").mkdir(parents=True)
        (source / "src/pyred").mkdir(parents=True)
        (source / "src/pyred/coeff_helper.h").write_text(
            "#include <cstdint>\nnamespace pyred { struct CoeffHelper {\n"
            "static uint64_t mod_mul(uint64_t,uint64_t,uint64_t); }; }\n")
        (source / "src/pyred.cpp").write_text(
            '#include "pyred/coeff_helper.h"\n#include "config.h"\n'
            '#ifndef MESON_SELECTED_CONFIG\n#error missing shared compiler configuration\n#endif\n'
            'extern int yaml_marker();\nuint64_t pyred::CoeffHelper::mod_mul(uint64_t a,uint64_t b,uint64_t m) {'
            'return yaml_marker() == 73 && CONFIG_VALUE == 19 ? (__uint128_t(a)*b)%m : 0; }\n')
        (source / "yaml.cpp").write_text("int yaml_marker() { return 73; }\n")
        (source / "tests/meson.build").write_text("")
        (source / "meson.build").write_text("""
project('probe-graph-check', 'cpp', default_options : ['cpp_std=c++17'])
add_project_arguments('-DMESON_SELECTED_CONFIG', language : 'cpp')
conf = configuration_data()
conf.set('CONFIG_VALUE', 19)
config_h = configure_file(output : 'config.h', configuration : conf)
yaml = static_library('yaml-selected-by-meson', 'yaml.cpp')
deps = [declare_dependency(link_with : yaml)]
deps_kira = []
include_dirs = include_directories('src', '.')
pyred = static_library('pyred', ['src/pyred.cpp', config_h], include_directories : include_dirs)
rpath = ''
subdir('tests')
""")
        installer.add_arithmetic_test(source)
        build = self.work / "build with spaces"
        for command in (["meson", "setup", str(build), str(source)],
                        ["meson", "test", "-C", str(build), "--suite", "feynfacet", "--print-errorlogs"]):
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    timeout=60)
            self.assertEqual(result.returncode, 0, result.stdout)
        log = (build / "meson-logs/testlog.txt").read_text()
        self.assertIn("mismatches=0/100072", log)


if __name__ == "__main__":
    unittest.main()
