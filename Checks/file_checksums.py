#!/usr/bin/env python3
"""Check SHA256SUMS, or refresh it after an intentional source/data change."""

import argparse
import hashlib
from pathlib import Path, PurePosixPath
import re

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "SHA256SUMS"


def source_files(directory):
    """Leave locally generated calculation and website products out."""
    for path in sorted(directory.iterdir()):
        if path.is_symlink():
            raise ValueError(f"Symbolic link in source tree: {path.relative_to(ROOT)}")
        if path.is_dir():
            if (path.name.startswith(".") and path.relative_to(ROOT).as_posix() != ".github") or path.name in {"Results", "__pycache__"}:
                continue
            if path.relative_to(ROOT).as_posix() in {"site/public", "site/artifact"}:
                continue
            yield from source_files(path)
        elif path.is_file():
            if path == MANIFEST or path.name in {"Paths.wl", ".DS_Store"}:
                continue
            if path.suffix in {".log", ".pyc"} or path.match("xform*.str"):
                continue
            if path.name.startswith(".") and path.name != ".gitignore":
                continue
            yield path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--write", action="store_true", help="regenerate SHA256SUMS")
    action.add_argument("--check", action="store_true", help="verify SHA256SUMS (default)")
    args = parser.parse_args()
    files = {p.relative_to(ROOT).as_posix(): p for p in source_files(ROOT)}
    if args.write:
        MANIFEST.write_text("".join(f"{sha(path)}  {name}\n" for name, path in sorted(files.items())))
        print(f"Wrote checksums for {len(files)} source and data files.")
        return 0
    expected = {}
    for line in MANIFEST.read_text().splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match:
            raise ValueError("Invalid checksum record")
        value, name = match.groups()
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or name in expected:
            raise ValueError(f"Invalid or repeated checksum path: {name}")
        expected[name] = value
    errors = [f"Missing: {name}" for name in sorted(set(expected) - set(files))]
    errors += [f"Not checksummed: {name}" for name in sorted(set(files) - set(expected))]
    errors += [f"Hash differs: {name}" for name in sorted(set(files) & set(expected))
               if sha(files[name]) != expected[name]]
    if errors:
        print("\n".join(errors))
        return 1
    print(f"Verified {len(expected)} source and data files.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        raise SystemExit(str(error))
