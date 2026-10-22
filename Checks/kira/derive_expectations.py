#!/usr/bin/env python3
"""Verify the nonzero fixture coefficients from Feynman-parameter Gamma factors.

Optional maintainer check: Python + SymPy, no Kira or physics workflow.
The derivation and the reasons for the scaleless zeros are in README.md.
"""
import json
from pathlib import Path
from sympy import Symbol, gamma, expand_func, simplify, cancel
from sympy.parsing.mathematica import parse_mathematica

d = Symbol("d")


def G(a, b):
    return gamma(a+b-d/2) * gamma(d/2-a) * gamma(d/2-b) / (
        gamma(a) * gamma(b) * gamma(d-a-b))


def J(a, b, c):
    return G(b, c) * G(a, b+c-d/2)


specs = json.loads(Path(__file__).with_name("expectations.json").read_text())
count = 0
for spec in specs:
    for case in spec["cases"]:
        expected = parse_mathematica(case["coefficient"])
        if expected == 0:
            continue  # No scale; these vanish by dimensional regularization.
        indices = parse_mathematica(case["target"]).args
        if spec["family"] == "virtual":
            derived = gamma(3-d/2) * gamma(d/2-2)**2 / gamma(d-3) / G(1, 1)
        else:
            derived = J(*indices[:3]) / J(1, 1, 1)
        derived = cancel(simplify(expand_func(derived)))
        if cancel(derived - expected) != 0:
            raise AssertionError((case["target"], derived, expected))
        print("PASS", case["target"], "coefficient =", derived)
        count += 1
print(f"Derived {count} nonzero rational coefficients independently of Kira.")
