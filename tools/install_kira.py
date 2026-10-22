#!/usr/bin/env python3
"""Build the reviewed Kira 3.1 + modular arithmetic fix in a fresh directory.

No sudo, Homebrew installation, shell-profile edits, or physics calculations.
Requires git, meson, ninja, pkg-config, a C/C++ compiler, Kira's libraries,
and a platform-appropriate Fermat executable. See docs/INSTALL-macOS.md.
"""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parent.parent
KIRA_COMMIT = "5760f42f8decd7eb5474c780db75b5e5e7b450a0"
FIREFLY_COMMIT = "0e155a9d879362521e9157e8fdac10b9553fcacb"
FLINT_VERSION = "2.8.4"
FLINT_SOURCE_SHA256 = "61df92ea8c8e9dc692d46c71d7f50aaa09a33d4ba08d02a1784730a445e5e4be"
FLINT_PATCH_SHA256 = "7fdfead7d94f60bf4093ee4e16fa499178be03054895ccf1b1430adc472f714f"
PATCH = ROOT / "tools/patches/kira-3.1-mod-mul.patch"
PATCH_SHA256 = "ac74e82e9f1e5ebf3945d6738159d2dee251789f252900fafae3de2ec8ad03ec"


def pin_wrap(path, old_revision, revision):
    """Reject unexpected dependency metadata instead of editing an unreviewed wrap."""
    text = path.read_text()
    old = "revision = " + old_revision
    if text.splitlines().count(old) != 1:
        raise ValueError(f"Unexpected dependency revision in {path}")
    path.write_text(re.sub(r"^" + re.escape(old) + r"$", "revision = " + revision,
                           text, flags=re.MULTILINE))


def add_arithmetic_test(source):
    """Let Meson use the exact compiler, configuration and dependencies of Kira."""
    shutil.copy2(ROOT / "Checks/kira_mod_mul.cpp", source / "tests/feynfacet_mod_mul.cpp")
    with (source / "tests/meson.build").open("a") as stream:
        stream.write("""
# FeynFacet's independent modular arithmetic reference, built with Kira's graph.
feynfacet_mod_mul = executable(
  'feynfacet-mod-mul', ['feynfacet_mod_mul.cpp'] + [config_h],
  include_directories : include_dirs,
  link_with : pyred,
  dependencies : deps + deps_kira,
  install_rpath : rpath,
  install : false)
test('feynfacet-mod-mul', feynfacet_mod_mul, suite : ['feynfacet'])
""")


def check_flint_wrap(source):
    # Kira's root wrap takes precedence over FireFly's nested flint2.wrap.
    # Meson verifies these archive checksums when it extracts the subproject.
    config = configparser.ConfigParser()
    config.read(source / "subprojects/flint2.wrap")
    if not config.has_section("wrap-file"):
        raise ValueError("Expected the reviewed archive-based FLINT wrap")
    wrap = config["wrap-file"]
    expected = {"directory": "flint-" + FLINT_VERSION,
                "source_hash": FLINT_SOURCE_SHA256, "patch_hash": FLINT_PATCH_SHA256}
    if any(wrap.get(key) != value for key, value in expected.items()):
        raise ValueError("FLINT source or build-definition pin differs from the reviewed version")
    return expected


def bundled_libraries(targets, installed):
    """Record the actual build products, including explicitly static libraries."""
    result = {}
    for name in ("flint", "yaml-cpp", "firefly"):
        rows = [target for target in targets if target["name"] == name
                and target["type"] in ("static library", "shared library")]
        if not rows or any(not target.get("subproject") for target in rows):
            raise ValueError(f"Meson did not build the bundled {name} library")
        result[name] = []
        for target in rows:
            for filename in target["filename"]:
                destination = installed.get(filename)
                if not destination:
                    raise ValueError(f"No installation destination recorded for {filename}")
                result[name].append({"type": target["type"], "build": filename,
                                     "installed": destination})
    return result


def darwin_path(path, binary, executable):
    for marker, origin in (("@loader_path", binary), ("@executable_path", executable)):
        if path == marker or path.startswith(marker + "/"):
            if origin is None:
                raise ValueError(f"Cannot resolve {path} without its loader location")
            return (Path(origin).parent / path[len(marker):].lstrip("/")).resolve()
    if Path(path).is_absolute():
        return Path(path).resolve()
    raise ValueError(f"Unsupported relative Mach-O library path: {path}")


def darwin_rpaths(output, binary, executable):
    paths = re.findall(r"^\s*cmd LC_RPATH\s*\n\s*cmdsize \d+\s*\n\s*path (.*?) \(offset \d+\)",
                       output, re.MULTILINE)
    return [darwin_path(path, binary, executable) for path in paths]


def resolve_darwin_reference(reference, binary, executable, rpaths):
    if reference.startswith("@rpath/"):
        for directory in rpaths:
            candidate = directory / reference[len("@rpath/"):]
            if candidate.is_file():
                return candidate.resolve()
        raise ValueError(f"Library cannot be found through LC_RPATH: {reference}")
    return darwin_path(reference, binary, executable)


def check_library_references(output, libraries, system, binary=None, executable=None, rpaths=(),
                             location="installed"):
    """Reject missing or external FLINT/yaml-cpp/FireFly in loader-tool output."""
    found = {}
    for line in output.splitlines():
        if system == "Darwin":
            match = re.match(r"\s+(.+?) \(compatibility version", line)
            if not match:
                continue
            reference = match.group(1)
        else:
            match = re.match(r"\s*(\S+)\s+=>\s+(.+?)(?:\s+\(0x[0-9a-fA-F]+\))?$", line)
            if not match:
                continue
            reference = match.group(2)
            if reference == "not found":
                raise ValueError(f"Kira ({location}) has a missing library: {line.strip()}")
        for name, products in libraries.items():
            if Path(match.group(1)).name.startswith("lib" + name + "."):
                allowed = {Path(item[location]).resolve() for item in products
                           if item["type"] == "shared library"}
                resolved = (resolve_darwin_reference(reference, binary, executable, rpaths)
                            if system == "Darwin" else Path(reference))
                if not resolved.is_absolute() or resolved.resolve() not in allowed:
                    raise ValueError(f"Kira ({location}) uses an unexpected {name}: {reference}")
                found[name] = str(resolved.resolve())
    return found


def inspect_linkage(executable, libraries, system, run, location="installed"):
    """Inspect Kira and FireFly using Meson's build or installed library paths."""
    if location not in ("build", "installed"):
        raise ValueError(f"Unknown library location: {location}")
    if system not in ("Darwin", "Linux"):
        raise ValueError("Library verification is supported on macOS and Linux")
    binaries = [executable]
    binaries += [Path(item[location]) for item in libraries["firefly"]
                 if item["type"] == "shared library"]
    found, output, load_commands = {}, {}, {}
    executable_rpaths = []
    if system == "Darwin":
        load_commands[str(executable)] = run(["otool", "-l", executable], capture=True)
        executable_rpaths = darwin_rpaths(load_commands[str(executable)], executable, executable)
    for binary in binaries:
        rpaths = executable_rpaths
        if system == "Darwin" and binary != executable:
            load_commands[str(binary)] = run(["otool", "-l", binary], capture=True)
            # FireFly is loaded by Kira: search its loader paths, then the
            # inherited executable paths. Never accept only a matching basename.
            rpaths = darwin_rpaths(load_commands[str(binary)], binary, executable) + executable_rpaths
        command = ["otool", "-L", binary] if system == "Darwin" else ["ldd", binary]
        listing = run(command, capture=True)
        output[str(binary)] = listing
        references = check_library_references(listing, libraries, system, binary, executable, rpaths,
                                              location=location)
        # On macOS a dylib's -L listing includes its own LC_ID_DYLIB. Its
        # existence cannot establish the executable -> FireFly dependency.
        if binary == executable and any(item["type"] == "shared library"
                                        for item in libraries["firefly"]) and "firefly" not in references:
            raise ValueError("Kira does not reference the expected FireFly shared library")
        found.update(references)
    for name, products in libraries.items():
        if any(item["type"] == "shared library" for item in products) and name not in found:
            raise ValueError(f"Cannot confirm the {location} {name} shared library")
    return {"references": found, "loaderOutput": output, "loadCommands": load_commands}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True, type=Path, help="new build directory; retained for diagnosis")
    parser.add_argument("--prefix", type=Path, default=Path.home() / ".local")
    parser.add_argument("--fermat", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=min(6, os.cpu_count() or 2))
    parser.add_argument("--timeout", type=int, default=1200, help="total time budget in seconds (default 1200)")
    parser.add_argument("--replace", action="store_true", help="explicitly permit replacing an existing prefix/bin/kira")
    args = parser.parse_args()
    prefix, work, fermat = args.prefix.expanduser().resolve(), args.workdir.expanduser().resolve(), args.fermat.expanduser().resolve()
    if args.jobs < 1 or args.timeout < 1:
        parser.error("--jobs and --timeout must be positive")
    if work.exists():
        parser.error("--workdir must not exist; use a fresh location")
    if not fermat.is_file() or not os.access(fermat, os.X_OK):
        parser.error("--fermat must name an executable for this platform")
    if (prefix / "bin/kira").exists() and not args.replace:
        parser.error("Kira already exists at the prefix; choose another --prefix or explicitly use --replace")
    if hashlib.sha256(PATCH.read_bytes()).hexdigest() != PATCH_SHA256:
        parser.error("Bundled patch checksum differs; stop and restore the reviewed patch")
    loader_tool = "otool" if platform.system() == "Darwin" else "ldd"
    for tool in ("git", "meson", "ninja", "pkg-config", "c++", loader_tool):
        if not shutil.which(tool):
            parser.error(f"Missing prerequisite: {tool}; see docs/INSTALL-macOS.md")
    work.mkdir(parents=True)
    started = time.monotonic()
    deadline = started + args.timeout
    env = {**os.environ, "FERMATPATH": str(fermat)}
    source, build = work / "kira", work / "kira/builddir"
    report = {"kiraCommit": KIRA_COMMIT, "fireflyCommit": FIREFLY_COMMIT,
              "patchSHA256": PATCH_SHA256, "prefix": str(prefix),
              "system": platform.platform(), "completed": False}

    with (work / "install.log").open("w") as log:
        def run(command, cwd=work, capture=False):
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError("Installation time budget exhausted")
            message = "+ " + shlex.join(map(str, command))
            print(message, flush=True)
            log.write(message + "\n"); log.flush()
            process = subprocess.Popen(list(map(str, command)), cwd=cwd, env=env,
                                       stdout=subprocess.PIPE if capture else log,
                                       stderr=log, text=True, start_new_session=True)
            try:
                output, _ = process.communicate(timeout=remaining)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.communicate()
                raise TimeoutError("Installation time budget exhausted; build process group stopped")
            if process.returncode:
                raise subprocess.CalledProcessError(process.returncode, command)
            if capture:
                log.write(output); log.flush()
                return output.strip()

        try:
            run(["git", "clone", "--branch", "3.1", "--depth", "1",
                 "https://gitlab.com/kira-pyred/kira.git", source])
            if run(["git", "rev-parse", "HEAD"], source, True) != KIRA_COMMIT:
                raise ValueError("Kira tag does not match reviewed source commit")
            run(["git", "apply", "--check", PATCH], source)
            run(["git", "apply", PATCH], source)
            # Upstream's kira-2 branch is mutable. Pin the reviewed FireFly source.
            pin_wrap(source / "subprojects/firefly.wrap", "kira-2", FIREFLY_COMMIT)
            report["flintSource"] = check_flint_wrap(source)
            add_arithmetic_test(source)
            options = ["--prefix=" + str(prefix), "--libdir=lib", "--buildtype=release",
                       "--wrap-mode=forcefallback", "-Dtest-level=full"]
            if platform.system() == "Darwin":
                if not shutil.which("brew"):
                    raise ValueError("Homebrew is required for this macOS build recipe")
                libs = [run(["brew", "--prefix", name], capture=True) + "/lib" for name in ("gmp", "mpfr")]
                flags = " ".join("-L" + path for path in libs)
                options += ["-Dc_link_args=" + flags, "-Dcpp_link_args=" + flags]
            run(["meson", "setup", *options, build], source)
            if run(["git", "rev-parse", "HEAD"], source / "subprojects/firefly", True) != FIREFLY_COMMIT:
                raise ValueError("FireFly checkout does not match the reviewed commit")
            targets = json.loads(run(["meson", "introspect", "--targets", build], capture=True))
            installed = json.loads(run(["meson", "introspect", "--installed", build], capture=True))
            report["libraries"] = bundled_libraries(targets, installed)
            run(["ninja", "-C", build, "-j", str(args.jobs)])
            report["testLogs"] = {"arithmetic": str(build / "meson-logs/feynfacet-arithmetic.txt"),
                                  "upstream": str(build / "meson-logs/kira-regressions.txt")}
            run(["meson", "test", "-C", build, "--suite", "feynfacet", "--logbase",
                 "feynfacet-arithmetic", "--print-errorlogs"])
            run(["meson", "test", "-C", build, "--suite", "kira", "--logbase", "kira-regressions",
                 "--timeout-multiplier", "10", "--print-errorlogs"])
            executables = [target for target in targets
                           if target["name"] == "kira" and target["type"] == "executable"]
            if len(executables) != 1 or len(executables[0]["filename"]) != 1:
                raise ValueError("Meson did not identify a unique Kira build executable")
            report["buildLinkage"] = inspect_linkage(
                Path(executables[0]["filename"][0]), report["libraries"], platform.system(), run,
                location="build")
            # A failed build check stops before any installation writes, including
            # with --replace. Installation itself is not a rollback mechanism.
            run(["meson", "install", "-C", build, "--no-rebuild"])
            report["linkage"] = inspect_linkage(prefix / "bin/kira", report["libraries"], platform.system(), run)
            report["version"] = run([prefix / "bin/kira", "--version"], capture=True)
            report["binarySHA256"] = hashlib.sha256((prefix / "bin/kira").read_bytes()).hexdigest()
            report["completed"] = True
        except (OSError, ValueError, configparser.Error, subprocess.SubprocessError, TimeoutError) as error:
            report["error"] = str(error)
            print(f"Installation failed: {error}\nDetails: {work / 'install.log'}", flush=True)
        finally:
            report["seconds"] = round(time.monotonic() - started, 2)
            (work / "installation.json").write_text(json.dumps(report, indent=2) + "\n")
    if report["completed"]:
        print(f"Installed and tested: {prefix / 'bin/kira'}\nRecord: {work / 'installation.json'}")
        print("Next: select this executable in Paths.wl if necessary, then run python3 feynfacet.py check --all-routes")
    return 0 if report["completed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
