#!/usr/bin/env python3
"""Redraw the two website TMD illustrations with TikZ-Feynman.

Run `python3 site/make_tmd_figures.py` from the repository root.
Requires pdflatex, the standalone and tikz-feynman packages, and Poppler
(pdftocairo and pdftoppm). Building the website uses the checked-in SVGs
and does not require TeX.
"""
from pathlib import Path
import shutil
import subprocess
import tempfile

SITE = Path(__file__).resolve().parent


def main():
    for command in ("pdflatex", "pdftocairo", "pdftoppm"):
        if shutil.which(command) is None:
            raise SystemExit(f"Missing {command}; install TeX Live and Poppler.")

    output = SITE / "figures" / "diagrams"
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="feynfacet-tmd-") as temporary:
        work = Path(temporary)
        for name in ("tmd-igg", "tmd-sq"):
            source = SITE / "figures" / "tikz" / f"{name}.tex"
            result = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                 f"-output-directory={work}", str(source)],
                cwd=work, capture_output=True, text=True,
            )
            if result.returncode:
                raise SystemExit(result.stdout + result.stderr)
            pdf = work / f"{name}.pdf"
            subprocess.run(
                ["pdftocairo", "-svg", str(pdf), str(output / f"{name}.svg")],
                check=True,
            )
            subprocess.run(
                ["pdftoppm", "-png", "-singlefile", "-scale-to", "900",
                 str(pdf), str(output / name)],
                check=True,
            )
            print(f"Written {name}.svg and {name}.png")


if __name__ == "__main__":
    main()
