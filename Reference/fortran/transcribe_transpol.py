#!/usr/bin/env python3
"""Transcribe the transversely polarized reference code into Mathematica.

The Fortran code of A. Mukherjee, M. Stratmann and W. Vogelsang, Phys. Rev.
D 72 (2005) 034011 [hep-ph/0506315] for single-inclusive hadron production in
collisions of transversely polarized protons (file me-transpol-ms.f, with the
driver drive-transpol.f; the code itself is NOT part of this repository) holds,
for every channel class J0, the pieces of the NLO hard cross section as
machine-generated Fortran expressions (Mathematica FortranForm output pasted
into the code):

    FBOR      F0(J0)          Born term, with the phase-space factor 1/(v(1-v))
    AVDEL     AVDEL(V,S)      coefficient of delta(1-w)   (ln(muR^2/s) inside)
    AVWPL     AVWPL(W,V,S)    coefficient of 1/(1-w)_+
    AVLO      AVLO(W,V,S)     coefficient of [ln(1-w)/(1-w)]_+
    STRUVn    STRUV(W,V,X3,S) the regular remainder
    CC(J0)    colour normalization (driver)

This script parses the fixed-form Fortran (continuation lines, D exponents,
** powers, the precomputed quantities of the common blocks PREV, PREW, PRELOG,
PRECOLOR and the local definitions of each routine) and writes the pieces for
J0 in {5, 6, 7, 11, 12} verbatim -- no simplification -- as Mathematica
expressions in v, w, s, muR, muF, muFF, CA, CF, Nf into

    Reference/pp-transverse-fortran.wl

together with the normalization map onto the basis of this package
(see Reference/fortran/README.md, section "The transversely polarized code").

Checks performed (and recorded in the header of the .wl file):

 1. every precomputed quantity (L1V, LVW, V2, CACF2, PRELO, ...) is resolved
    from a table that is verified numerically against the definitions parsed
    from the code itself (PRECALC, FBOR, AVWPL, AVDEL, AVLO, STRUVn);
 2. the Fortran expressions are evaluated in Python (mpmath, 40 digits) at
    several exact rational points, and the generated Mathematica file is
    evaluated at the same points with wolframscript (30 digits); the two must
    agree to 1e-20 relative;
 3. the raw FortranForm files proc5/6/7/11/12-transpol.f shipped with the
    code (the same expressions before they were pasted into the routines:
    blocks tagged dw*(..) = delta term, (..)*wp = 1/(1-w)_+ term,
    lwp*(..) = [ln(1-w)/(1-w)]_+ term, the remainder = regular term; they are
    not referenced by the makefile and are not compiled) are parsed
    independently and compared numerically with the routines' expressions.

Usage:
    python3 Reference/fortran/transcribe_transpol.py /path/to/ATT [--out FILE]
                                                     [--no-wolfram]
    The directory must hold me-transpol-ms.f, drive-transpol.f and the
    proc*-transpol.f files.  Temporary files go to --tmp (default: a new
    temporary directory).
"""

import argparse
import os
import tempfile
import re
import subprocess
import sys
from fractions import Fraction

try:
    import mpmath
    from mpmath import mp, mpf
except ImportError:  # pragma: no cover
    sys.exit("mpmath is required (pip install mpmath)")

mp.dps = 40

CLASSES = [5, 6, 7, 11, 12]
CLASS_NAMES = {5: "q qb -> qp", 6: "q q -> q", 7: "q q -> g", 11: "q qb -> q", 12: "q qb -> g"}
PIECES = ["F0", "AVDEL", "AVWPL", "AVLO", "STRUV"]

# ---------------------------------------------------------------------------
# Fixed-form Fortran reading
# ---------------------------------------------------------------------------

def logical_lines(path):
    """Join continuation lines; drop comments.  Returns [(first, last, text)].

    Fixed form: a comment has C, c, * or ! in column 1; columns 1-5 hold the
    label, column 6 a continuation mark (anything but blank or 0), the
    statement starts in column 7.  Text after an exclamation mark is a comment.
    Characters beyond column 72 are ignored by a fixed-form compiler; the
    script reports any non-blank character there (there are none in these
    files, only trailing blanks).
    """
    out = []
    cur = None
    long_lines = []
    with open(path, encoding="latin-1") as fh:
        for no, raw in enumerate(fh, 1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line[0] in "Cc*!":
                continue
            if "!" in line:
                line = line[: line.index("!")]
            if len(line) > 72 and line[72:].strip():
                long_lines.append((no, line[72:]))
            body = line[:72]
            if len(body) > 5 and body[5] not in " 0" and not body[:5].strip():
                if cur is None:
                    raise ValueError("%s:%d: continuation without a statement" % (path, no))
                cur[2] += " " + body[6:]
                cur[1] = no
            else:
                if cur is not None:
                    out.append(tuple(cur))
                cur = [no, no, body[6:] if len(body) > 6 else body.strip()]
    if cur is not None:
        out.append(tuple(cur))
    if long_lines:
        raise ValueError("%s: non-blank text beyond column 72: %r" % (path, long_lines[:3]))
    return out


def split_units(lines):
    """Split logical lines into program units {name: [(first,last,text),...]}."""
    units = {}
    name = None
    start = re.compile(r"^\s*(?:DOUBLE\s+PRECISION\s+|REAL\*8\s+)?(FUNCTION|SUBROUTINE|PROGRAM)\s+([A-Za-z0-9_]+)", re.I)
    for first, last, text in lines:
        m = start.match(text)
        if m:
            name = m.group(2).upper()
            units[name] = []
            continue
        if name is None:
            continue
        if re.match(r"^\s*END\s*$", text, re.I):
            units[name].append((first, last, text))
            name = None
            continue
        units[name].append((first, last, text))
    return units


# ---------------------------------------------------------------------------
# Tokenizer and converters
# ---------------------------------------------------------------------------

NUM = re.compile(r"(\d+\.\d*|\.\d+|\d+)([DdEe][+-]?\d+)?")
IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def tokenize(expr):
    tokens = []
    i = 0
    n = len(expr)
    while i < n:
        c = expr[i]
        if c.isspace():
            i += 1
            continue
        if c.isdigit() or (c == "." and i + 1 < n and expr[i + 1].isdigit()):
            m = NUM.match(expr, i)
            tokens.append(("num", m.group(0)))
            i = m.end()
            continue
        m = IDENT.match(expr, i)
        if m:
            tokens.append(("id", m.group(0).lower()))
            i = m.end()
            continue
        if expr.startswith("**", i):
            tokens.append(("op", "**"))
            i += 2
            continue
        if c in "+-*/(),":
            tokens.append(("op", c))
            i += 1
            continue
        raise ValueError("cannot tokenize %r at %d: %r" % (expr, i, expr[i:i + 20]))
    return tokens


def number(text):
    """Exact rational value of a Fortran constant ('3D0', '1.', '.5', '2')."""
    return Fraction(re.sub(r"[DdEe]", "e", text))


# The precomputed quantities and local variables of the code, resolved to
# Mathematica (strings) and to Python values (functions of a point).
# Every entry is verified numerically against the definition parsed from the
# code (verify_definitions below).
POWERS = {("v%d" % k): ("(v^%d)" % k) for k in range(2, 13)}
POWERS.update({("w%d" % k): ("(w^%d)" % k) for k in range(2, 13)})

MMA = {
    "v": "v", "w": "w", "s": "s", "shd": "s", "ca": "CA", "nc": "CA", "cf": "CF", "nf": "Nf", "pi": "Pi",
    "l1v": "(Log[1 - v])", "lv": "(Log[v])", "l1w": "(Log[1 - w])", "lw": "(Log[w])",
    "lvw": "(Log[(1 - v)/(1 - v*w)]/(1 - w))", "l1vw": "(Log[1 - v + v*w]/(1 - w))",
    "lmu": "(Log[muR^2/s])", "lms": "(Log[muF^2/s])", "lmss": "(Log[muFF^2/s])",
    "ca2": "(CA^2)", "ca4": "(CA^4)", "cacf": "(CA*CF)", "cacf2": "(CA*CF^2)",
    # FBOR locals
    "vc": "(CA^2 - 1)", "nc2": "(CA^2)", "vm": "(1 - v)", "vm2": "((1 - v)^2)", "prelo": "(Pi/(s*v*(1 - v)))",
    # tags of the raw FortranForm files (proc*-transpol.f)
    "dw": "1", "wp": "1", "lwp": "1",
}
MMA.update(POWERS)

# what the code defines these as (lower case, blanks removed); checked
EXPECTED_DEFINITIONS = {
    "PRECALC": {
        "l1v": "dlog(1.d0-v)", "lv": "dlog(v)", "l1w": "dlog(1.d0-w)", "lw": "dlog(w)",
        "lvw": "dlog((1.d0-v)/(1.d0-v*w))/(1.d0-w)", "l1vw": "dlog(1.d0-v+v*w)/(1.d0-w)",
        "lmu": "dlog(q2mu/s)", "lms": "dlog(q2fac/s)", "lmss": "dlog(q2frag/s)",
        "cacf": "ca*cf", "ca2": "ca**2", "ca4": "ca**4",
        **{"v%d" % k: "v**%d" % k for k in range(2, 13)}, **{"w%d" % k: "w**%d" % k for k in range(2, 13)},
    },
    "FBOR": {"v2": "v**2", "v3": "v**3", "v4": "v**4", "vc": "(nc**2-1.d0)", "nc2": "nc**2",
             "vm": "(1.d0-v)", "vm2": "vm**2", "prelo": "pi/shd/v/vm"},
    "AVWPL": {"m": "dsqrt(q2fac)", "mp": "dsqrt(q2frag)", "l1v": "dlog(1.-v)", "lv": "dlog(v)",
              "lms": "dlog(m**2/s)", "lmss": "dlog(mp**2/s)", "nf": "2.*gtr"},
    "AVDEL": {"m": "dsqrt(q2fac)", "mp": "dsqrt(q2frag)", "delw": "1.", "l1v": "dlog(1.-v)", "lv": "dlog(v)",
              "lmu": "dlog(q2mu/s)", "lms": "dlog(m**2/s)", "lmss": "dlog(mp**2/s)", "nf": "2.*gtr"},
    "AVLO": {"m": "dsqrt(q2fac)", "mp": "dsqrt(q2frag)", "loplus": "1.", "l1v": "dlog(1.-v)", "lv": "dlog(v)",
             "lmu": "dlog(q2mu/s)", "lms": "dlog(m**2/s)", "lmss": "dlog(mp**2/s)", "nf": "2.*gtr"},
    "STRUV": {"cacf2": "ca*cf**2", "nf": "2.d0*gtr"},
}


def point_values(pt):
    """Python (mpmath) values of all symbols at a point given as Fractions."""
    F = lambda x: mpf(x.numerator) / x.denominator
    v, w, s = F(pt["v"]), F(pt["w"]), F(pt["s"])
    CA, CF, Nf = F(pt["CA"]), F(pt["CF"]), F(pt["Nf"])
    muR2, muF2, muFF2 = F(pt["muR2"]), F(pt["muF2"]), F(pt["muFF2"])
    log = mpmath.log
    S = {
        "v": v, "w": w, "s": s, "shd": s, "ca": CA, "nc": CA, "cf": CF, "nf": Nf, "pi": mp.pi,
        "l1v": log(1 - v), "lv": log(v), "l1w": log(1 - w), "lw": log(w),
        "lvw": log((1 - v) / (1 - v * w)) / (1 - w), "l1vw": log(1 - v + v * w) / (1 - w),
        "lmu": log(muR2 / s), "lms": log(muF2 / s), "lmss": log(muFF2 / s),
        "ca2": CA ** 2, "ca4": CA ** 4, "cacf": CA * CF, "cacf2": CA * CF ** 2,
        "vc": CA ** 2 - 1, "nc2": CA ** 2, "vm": 1 - v, "vm2": (1 - v) ** 2, "prelo": mp.pi / (s * v * (1 - v)),
        "dw": mpf(1), "wp": mpf(1), "lwp": mpf(1), "delw": mpf(1), "loplus": mpf(1),
        # quantities of the code's common blocks used only in the local definitions
        "q2mu": muR2, "q2fac": muF2, "q2frag": muFF2, "gtr": Nf / 2,
        "m": mpmath.sqrt(muF2), "mp": mpmath.sqrt(muFF2),
    }
    for k in range(2, 13):
        S["v%d" % k] = v ** k
        S["w%d" % k] = w ** k
    return S


FUNCS_PY = {"dlog": "mpmath.log", "dsqrt": "mpmath.sqrt", "log": "mpmath.log", "sqrt": "mpmath.sqrt"}
FUNCS_MMA = {"dlog": "Log", "dsqrt": "Sqrt", "log": "Log", "sqrt": "Sqrt"}


def to_python(tokens):
    parts = []
    for kind, val in tokens:
        if kind == "num":
            fr = number(val)
            parts.append("(mpf(%d)/%d)" % (fr.numerator, fr.denominator))
        elif kind == "id":
            if val in FUNCS_PY:
                parts.append(FUNCS_PY[val])
            else:
                parts.append("S[%r]" % val)
        else:
            parts.append(val)
    return "".join(parts)


def evaluate(tokens, S):
    return eval(to_python(tokens), {"S": S, "mpf": mpf, "mpmath": mpmath})


def to_mathematica(tokens, table=MMA):
    """Verbatim conversion: ** -> ^, constants exact, symbols resolved."""
    out = []
    stack = []
    prev = None
    for kind, val in tokens:
        if kind == "num":
            fr = number(val)
            out.append(str(fr.numerator) if fr.denominator == 1 else "(%d/%d)" % (fr.numerator, fr.denominator))
        elif kind == "id":
            if val in FUNCS_MMA:
                out.append(FUNCS_MMA[val])
                stack.append("[")
                continue
            if val not in table:
                raise KeyError("unknown Fortran symbol %r" % val)
            out.append(table[val])
        else:
            if val == "(":
                if stack and stack[-1] == "[" and prev == ("fn",):
                    out.append("[")
                    stack[-1] = "]"
                else:
                    out.append("(")
                    stack.append(")")
            elif val == ")":
                out.append(stack.pop())
            elif val == "**":
                out.append("^")
            elif val in "+-":
                unary = prev is None or prev == ("op", "(") or (prev[0] == "op" and prev[1] not in ")")
                out.append(val if unary else " %s " % val)
            else:
                out.append(val)
        prev = ("fn",) if (kind == "id" and val in FUNCS_MMA) else (kind, val)
    if stack:
        raise ValueError("unbalanced parentheses")
    return "".join(out)


def wrap(text, width=96, indent="      "):
    """Break a Mathematica expression after binary +/- operators."""
    lines = []
    cur = ""
    i = 0
    while i < len(text):
        j = len(text)
        for op in (" + ", " - "):
            k = text.find(op, i)
            if k != -1:
                j = min(j, k + len(op))
        piece = text[i:j]
        if cur and len(cur) + len(piece) > width:
            lines.append(cur.rstrip())
            cur = indent + piece
        else:
            cur += piece
        i = j
    if cur:
        lines.append(cur.rstrip())
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

ASSIGN = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$")
F0ASSIGN = re.compile(r"^\s*F0\((\d+)\)\s*=\s*(.*?)\s*$", re.I)
IFJ0 = re.compile(r"^\s*(?:ELSE\s*)?IF\s*\(\s*J0\s*\.EQ\.\s*(\d+)\s*\)\s*THEN\s*$", re.I)


def extract_pieces(units):
    """Return pieces[J0][name] = (tokens, (first,last)) and definitions."""
    pieces = {j: {} for j in CLASSES}
    definitions = {}
    # PRECALC
    definitions["PRECALC"] = {}
    for first, last, text in units["PRECALC"]:
        m = ASSIGN.match(text)
        if m and "(" not in m.group(1):
            definitions["PRECALC"][m.group(1).lower()] = m.group(2)
    # FBOR
    definitions["FBOR"] = {}
    for first, last, text in units["FBOR"]:
        m = F0ASSIGN.match(text)
        if m:
            j = int(m.group(1))
            if j in pieces:
                pieces[j]["F0"] = (tokenize(m.group(2)), (first, last))
            continue
        m = ASSIGN.match(text)
        if m:
            definitions["FBOR"][m.group(1).lower()] = m.group(2)
    # AVWPL, AVDEL, AVLO
    for name in ("AVWPL", "AVDEL", "AVLO"):
        definitions[name] = {}
        current = None
        for first, last, text in units[name]:
            m = IFJ0.match(text)
            if m:
                current = int(m.group(1))
                continue
            m = ASSIGN.match(text)
            if not m:
                continue
            lhs, rhs = m.group(1), m.group(2)
            if lhs.upper() == name:
                if current in pieces:
                    pieces[current][name] = (tokenize(rhs), (first, last))
            elif current is None:
                definitions[name][lhs.lower()] = rhs
    # STRUVn
    definitions["STRUV"] = {}
    for j in CLASSES:
        unit = units["STRUV%d" % j]
        for first, last, text in unit:
            m = ASSIGN.match(text)
            if not m:
                continue
            if m.group(1).upper() == "STRUV%d" % j:
                pieces[j]["STRUV"] = (tokenize(m.group(2)), (first, last))
            else:
                definitions["STRUV"].setdefault(m.group(1).lower(), m.group(2))
    for j in CLASSES:
        missing = [p for p in PIECES if p not in pieces[j]]
        if missing:
            raise ValueError("class %d: pieces not found: %s" % (j, missing))
    return pieces, definitions


def extract_cc(driver_lines):
    """CC(J0) of the driver: three branches selected by J0."""
    branches = []
    cond = None
    for first, last, text in driver_lines:
        m = re.match(r"^\s*(?:ELSE\s*)?IF\s*\((.*)\)\s*THEN\s*$", text, re.I)
        if m and "J0" in m.group(1).upper():
            cond = [int(x) for x in re.findall(r"J0\.EQ\.(\d+)", m.group(1), re.I)]
            continue
        m = re.match(r"^\s*ELSE\s*$", text, re.I)
        if m and cond is not None:
            cond = "else"
            continue
        m = re.match(r"^\s*CC\(J0\)\s*=\s*(.*?)\s*$", text, re.I)
        if m and cond is not None:
            branches.append((cond, m.group(1)))
            if cond == "else":
                cond = None
    table = {"vc": "(CA^2 - 1)", "nc": "CA"}
    cc = {}
    for j in CLASSES:
        for cond, expr in branches:
            if cond == "else" or j in cond:
                cc[j] = (tokenize(expr), to_mathematica(tokenize(expr), table))
                break
    return cc, branches


def verify_definitions(definitions, points):
    """Every local definition of the code agrees with the resolution table."""
    report = []
    for unit, defs in definitions.items():
        expected = EXPECTED_DEFINITIONS[unit]
        got = {k: re.sub(r"\s+", "", v).lower() for k, v in defs.items()}
        if got != expected:
            raise ValueError("definitions of %s changed: %r vs expected %r" % (unit, got, expected))
        for name, rhs in defs.items():
            toks = tokenize(rhs)
            for pt in points:
                S = point_values(pt)
                lhs = S[name]
                val = evaluate(toks, S)
                if abs(val - lhs) > mpf(10) ** -35 * (1 + abs(lhs)):
                    raise ValueError("%s: %s = %s does not match the table" % (unit, name, rhs))
        report.append("%s (%d definitions)" % (unit, len(defs)))
    return report


# ---------------------------------------------------------------------------
# The raw FortranForm files proc*-transpol.f
# ---------------------------------------------------------------------------

def proc_blocks(path):
    """Blocks of a proc*-transpol.f file: expressions separated by blank lines.

    These files are FortranForm output (continuation marks in column 6, a
    stray backslash at the end of one line in proc7); the tags dw, wp, lwp
    multiply the delta, 1/(1-w)_+ and [ln(1-w)/(1-w)]_+ terms.
    """
    blocks = []
    cur = []
    with open(path, encoding="latin-1") as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r").rstrip("\\")
            if not line.strip():
                if cur:
                    blocks.append(" ".join(cur))
                    cur = []
                continue
            body = line[6:] if len(line) > 6 else ""
            cur.append(body)
    if cur:
        blocks.append(" ".join(cur))
    out = {}
    for b in blocks:
        toks = tokenize(b)
        ids = {t[1] for t in toks if t[0] == "id"}
        tag = "AVDEL" if "dw" in ids else "AVWPL" if "wp" in ids else "AVLO" if "lwp" in ids else "STRUV"
        if tag in out:
            raise ValueError("%s: two blocks tagged %s" % (path, tag))
        out[tag] = toks
    return out


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

POINTS = [
    {"v": Fraction(3, 10), "w": Fraction(1, 5), "s": Fraction(1), "muR2": Fraction(9, 10), "muF2": Fraction(7, 10),
     "muFF2": Fraction(13, 10), "Nf": Fraction(4), "CA": Fraction(3), "CF": Fraction(4, 3)},
    {"v": Fraction(11, 20), "w": Fraction(3, 5), "s": Fraction(7, 5), "muR2": Fraction(1, 2), "muF2": Fraction(2),
     "muFF2": Fraction(3, 2), "Nf": Fraction(5), "CA": Fraction(3), "CF": Fraction(4, 3)},
    {"v": Fraction(4, 5), "w": Fraction(9, 10), "s": Fraction(1), "muR2": Fraction(9, 10), "muF2": Fraction(7, 10),
     "muFF2": Fraction(13, 10), "Nf": Fraction(3), "CA": Fraction(2), "CF": Fraction(3, 4)},
    {"v": Fraction(1, 2), "w": Fraction(1, 2), "s": Fraction(1), "muR2": Fraction(1), "muF2": Fraction(1),
     "muFF2": Fraction(1), "Nf": Fraction(4), "CA": Fraction(3), "CF": Fraction(4, 3)},
    {"v": Fraction(2, 5), "w": Fraction(7, 10), "s": Fraction(2), "muR2": Fraction(1), "muF2": Fraction(3),
     "muFF2": Fraction(1, 3), "Nf": Fraction(6), "CA": Fraction(5, 2), "CF": Fraction(7, 3)},
]


def mma_point(pt):
    def fr(x):
        return "%d" % x.numerator if x.denominator == 1 else "%d/%d" % (x.numerator, x.denominator)
    return ("<|v -> %s, w -> %s, s -> %s, muR -> Sqrt[%s], muF -> Sqrt[%s], muFF -> Sqrt[%s], Nf -> %s, CA -> %s, CF -> %s|>"
            % (fr(pt["v"]), fr(pt["w"]), fr(pt["s"]), fr(pt["muR2"]), fr(pt["muF2"]), fr(pt["muFF2"]), fr(pt["Nf"]), fr(pt["CA"]), fr(pt["CF"])))


HEADER = r'''(* pp-transverse-fortran.wl: the NLO pieces of the Fortran code of
   A. Mukherjee, M. Stratmann and W. Vogelsang, Phys. Rev. D 72 (2005) 034011
   [hep-ph/0506315], for single-inclusive hadron production in collisions of
   two transversely polarized protons (the double-spin asymmetry A_TT), as
   Mathematica expressions.  Generated by Reference/fortran/transcribe_transpol.py
   from the files me-transpol-ms.f (routines FBOR, AVDEL, AVWPL, AVLO, STRUV5,
   STRUV6, STRUV7, STRUV11, STRUV12, PRECALC, DPLUS) and drive-transpol.f
   (CC(J0), the colour constants); the Fortran code is not part of this
   repository.  Nothing was simplified: every expression is the code's
   expression, token by token, with ** written as ^, the D exponents of the
   constants removed (3D0 = 3, 1. = 1, .5 = 1/2) and the precomputed quantities
   of the code resolved to their definitions:

     PRECALC (common blocks PREV, PREW, PRELOG, PRECOLOR):
       V2 ... V12 = v^2 ... v^12,  W2 ... W12 = w^2 ... w^12,
       L1V = Log[1 - v],  LV = Log[v],  L1W = Log[1 - w],  LW = Log[w],
       LVW = Log[(1 - v)/(1 - v w)]/(1 - w),  L1VW = Log[1 - v + v w]/(1 - w),
       LMU = Log[muR^2/s],  LMS = Log[muF^2/s],  LMSS = Log[muFF^2/s]
         (Q2MU, Q2FAC, Q2FRAG are the squared renormalization, initial-state
          and final-state factorization scales; S the partonic s),
       CACF = CA CF,  CA2 = CA^2,  CA4 = CA^4;
     local to AVWPL, AVDEL, AVLO: l1v, lv, lmu, lms, lmss as above,
       Nf = 2 GTR (GTR = NF/2 of the driver);
     local to STRUVn: CACF2 = CA CF^2, NF = 2 GTR;
     local to FBOR: PRELO = Pi/(s v (1 - v)), VM = 1 - v, VM2 = (1 - v)^2,
       NC2 = CA^2, VC = CA^2 - 1, V2, V3, V4 (NC of the common block PARAMS is
       the CA of the driver).

   The resolution table was verified numerically against the definitions
   parsed from the code (routines %(verified)s).

   Colour: the formulas are written in terms of CA and CF (common block CONS:
   CA = NC, CF = CF of PARAMS); the driver sets CF = 4/3, CA = 3 (SU(3)
   numbers) and NF = 4.  Whether the expressions hold for general CA with
   CF = (CA^2 - 1)/(2 CA), or only at SU(3), cannot be decided from the code;
   step 4 of Processes/HadronicTT compares them with our classes at CA = 2, 3, 5.
   CC(J0) = 8 NC^2 for all five classes here (the driver's assignment for
   classes other than 8, 9, 10, 13, 14 [8 VC NC] and 15, 16 [8 VC^2]).

   Classes (J0 of the code; the five that do not vanish for transverse
   polarization, driver's PROCI flags): 5 = QI QBI -> QK, 6 = QI QI -> QI,
   7 = QI QI -> G (no Born term), 11 = QI QBI -> QI, 12 = QI QBI -> G.
   The code's AVWPL and AVLO depend only on v (no w dependence), so the
   regular part is STRUV alone and the plus-distribution subtraction terms
   (AVWPL(w) - AVWPL(1))/(1-w) of the earlier codes vanish identically.
   AVDEL holds the ln(muR^2/s) term (no separate AVDELMU).

   The hadronic assembly (function DPLUS of me-transpol-ms.f, lines %(dplus)s),
   identical to that of the longitudinally polarized code ourh_p.f up to the
   two points marked NEW:
     - Born: BORN(J0) = (F01(J0) GRRT(J0) + F02(J0) GRRC(J0)) BXJAC with
       F01 = FBOR(V), F02 = FBOR(1-V) (the leg-exchanged class), GRRT/GRRC the
       products of transversity densities and fragmentation functions of STRU,
       BXJAC = (X3MAX-X3MIN)(VMAX-VMIN)/(BX1 BX2 X3^2) the Monte-Carlo
       Jacobian; FBOR includes the phase-space factor 1/(v(1-v)) of the
       hadronic formula (PRELO = Pi/(SHD V (1-V))).
     - delta(1-w) and plus-distribution terms (FHODEL1, classes IA1):
       [FDEL1 + FVWPL1(1) ln(1-WMIN) + FVLO1(1) ln(1-WMIN)^2/2]/(1-V)
        - (WMAX-WMIN)/(1-V) [FVWPL1(1) + FVLO1(1) ln(1-W)]/(1-W),
       divided by CC(J0), times BXJAC; FDEL1 = AVDEL(V,SHD)/SHD,
       FVWPL1(W) = AVWPL(W,V,SH)/SH, FVLO1(W) = AVLO(W,V,SH)/SH (SH the
       partonic s at the NLO momentum fractions): the plus prescription on
       [0, 1] in w, the subtraction at w = 1 (FVWPL1(UN, ...), UN = 1).
     - regular and w-dependent terms (FHOREST1):
       (WMAX-WMIN)/(1-V) [FVWPL1(W)/(W(1-W)) + FVLO1(W) ln(1-W)/(W(1-W))
        + FRESC1(W)/W] / CC(J0) times XJAC, FRESC1 = STRUV(W,V,X3,SH)/SH;
       XJAC = W BXJAC, so the explicit 1/W cancels against the Jacobian: the
       NLO slot carries the same weight as the Born slot (the factor 1/w of the
       partonic normalization sits in the Jacobian, as in the earlier codes).
     - the leg-exchanged pieces (FHODEL2, FHOREST2, classes IA2/IA3) use
       VX = 1 - V W, WX = (1-V)/(1-V W) with the corresponding logarithms.
     - sum: GHD = sum_J0 [(FHODEL1 GRRT + FHODEL2 GRRC) ALPASHO + BORN] PROCI,
       GHE = sum_J0 (FHOREST1 GPPT + FHOREST2 GPPC) ALPASHO PROCI,
       DPLUS = (GHD + GHE) ALPORD^2, ALPORD = 4 Pi ALPHAS(Q2MU) = alpha_s
       (ALPHAS returns alpha_s/(4 Pi)); NEW: ALPASHO = ALPORD IORD * 0., i.e.
       the shipped driver multiplies the whole NLO weight by zero (a debugging
       leftover; the analytic pieces are unaffected).
     - DPLUS = FACIN DPLUS HC2/(Pi S) XJETA XJPT with FACIN = Pi
       (dsigma/dy/dpT^2), 1 (E dsigma/d^3p), 2 Pi PT (dsigma/dpT/dy, the
       driver's ISIGM = 3) or 2 Pi PT^4; HC2 = (hbar c)^2 in pb GeV^2;
       NEW: finally DPLUS = DPLUS 4/Pi, "include cos(2 phi) dependence (four
       quadrants)".
   Hence the TT code computes, in the slots where the longitudinal code holds
   E dsigma/d^3p, (1/(Pi S)) x the convolution of the pieces, times 4/Pi.

   Paper: Eq. (3) writes d^3 dsigma/(dpT deta dPhi) = pT/(Pi S) sum int
   dz/z^2 dv/(v(1-v)) dw/w deltaf deltaf D [d dsigma^(0)/(dv dPhi) delta(1-w)
   + (alpha_s/Pi) d dsigma^(1)/(dv dw dPhi)], and Eq. (10)
     s d dsigma^(1)/(dv dw dPhi) = cos(2 Phi) (alpha_s/Pi)^2 [
       (A0 delta(1-w) + B0/(1-w)_+ + C0) ln(mu^2/s) + A delta(1-w)
       + B/(1-w)_+ + C + D [ln(1-w)/(1-w)]_+ + E ln w + F ln v + G ln(1-v)
       + H ln(1-w) + I ln(1-vw) + J ln(1-v+vw) + K ln w/(1-w)
       + L ln((1-v)/(1-vw))/(1-w) + M ln(1-v+vw)/(1-w) ].
   The code's pieces are these functions in the same slots: AVDEL <-> A (and
   A0 through lmu, with lms, lmss the separate scale logarithms), AVWPL <-> B
   (B0), AVLO <-> D, STRUV <-> C, E ... M (the LVW and L1VW of PRECALC are the
   L and M structures); the overall normalization of the code's pieces is not
   the paper's (alpha_s/Pi)^2 but the one fixed by DPLUS, see "Map" below.

   Normalization map onto the basis of this package (derivation in
   Reference/fortran/README.md): with sigma = s E dsigma/d^3p =
   Born delta(1-w) + P c_regular + P(v,1) [c_D0 D0 + c_D1 D1 + c_delta delta(1-w)],
   P = alphaS^3/(Pi s v w), our deltasigma_TT = T_XX the coefficient of
   cos(2 Phi) in sigma, the code's slots hold ONE HALF of the cos(2 Phi)
   coefficient in the normalization in which the unpolarized code's slots hold
   the azimuth-integrated cross section (equivalently: the cos(2 Phi)-weighted
   azimuthal integral, int dPhi cos(2 Phi) ..., of the Phi-differential cross
   section where the unpolarized code has int dPhi ...).  Therefore
     Born_TT      = alphaS^2 2 (1 - v) F0(J0)/Pi      [unpolarized: Pi b/(1-v) = F0]
     c_regular    = 2 v w STRUV/CC   (+ the vanishing subtraction terms)
     c_D0         = 2 v AVWPL(w=1)/CC,   c_D1 = 2 v AVLO(w=1)/CC,
     c_delta      = 2 v AVDEL/CC,
   i.e. {R, B, C, A}_TT = 2 alphaS^3/(Pi s CC) {STRUV + ..., AVWPL1, AVLO1, AVDEL}
   against {R, B, C, A} = alphaS^3/(Pi s CC) {...} of Processes/Hadronic.
   PPTransverseFortranCoefficients[J0] returns these.
   PPTransverseFortranExchangedCoefficients[J0] returns the same map for the
   ordering with the incoming partons exchanged (t <-> u), built with the
   code's FDEL2, FVWPL2, FVLO2, FRESC2 and FBOR(1-v) (see the end of the file).

   Transcription check (transcribe_transpol.py): %(check)s
*)
'''


# The leg-exchanged ordering (DPLUS slots FHODEL2/FHOREST2), written after the
# pieces; the line numbers refer to me-transpol-ms.f (status 03/11/04).
EXCHANGED_WL = r'''(* The leg-exchanged ordering of class J0 (DPLUS, slots FHODEL2/FHOREST2, used
   for the classes IA2 = {1, 3, 5, 11, 13, 14} and IA3 = {1, 3, 5, 8, 9, 10, 11,
   13, 14}): the same class with the two incoming partons interchanged, t <-> u.
   For J0 = 5 this is qbar q -> q' with the antiquark from hadron A, i.e. (by
   charge conjugation) q qbar -> qbar' with the antiquark qbar' observed.  The
   code builds it from the same pieces with its functions (me-transpol-ms.f):
     FDEL2   lines 417-429:  (AVDEL(1-v) + Log[v/(1-v)] AVWPL(1, 1-v)
                              + 1/2 AVLO(1, 1-v) Log[(1-v)/v]^2) (1-v)/v
     FVWPL2  lines 442-453:  (AVWPL(wx, vx) + AVLO(wx, vx) Log[v/vx]) vx/v
     FVLO2   lines 466-477:  AVLO(wx, vx) vx/v
     FRESC2  lines 490-501:  STRUV(wx, vx)   (after CALL PRECALC(VX, WX), line 247)
     F02     line 147: CALL FBOR(1-V, SHD, F02)
   with vx = 1 - v w, wx = (1-v)/(1-v w) (lines 245-246), each divided by
   the same partonic s and assembled by DPLUS exactly as the first ordering
   (FHODEL2, lines 219-223: plus prescription on [0, 1] in w subtracted at w = 1, i.e. at
   vx = 1-v, wx = 1; FHOREST2, lines 251-253: FVWPL2(w)/(w(1-w)) + FVLO2(w) Log[1-w]/(w(1-w))
   + FRESC2(w)/w, same weights BXJAC and XJAC as the first ordering).  The map
   onto our basis is therefore the one of PPTransverseFortranCoefficients with
   these combinations in place of AVDEL, AVWPL, AVLO, STRUV, F0; FVWPL2 and
   FVLO2 depend on w, so the regular part contains their remainders
   (FVWPL2(w) - FVWPL2(1))/(1-w) and (FVLO2(w) - FVLO2(1)) Log[1-w]/(1-w). *)
PPTransverseFortranExchangedCoefficients[J0_Integer] := Module[{p = PPTransverseFortran["Classes"][J0], vx = 1 - v w, wx = (1 - v)/(1 - v w), at, del2, wpl2, lo2, resc2, f02},
  at[piece_, vv_, ww_] := p[piece] /. {v -> vv, w -> ww};
  del2 = (at["AVDEL", 1 - v, 1] + Log[v/(1 - v)] at["AVWPL", 1 - v, 1] + 1/2 at["AVLO", 1 - v, 1] Log[(1 - v)/v]^2) (1 - v)/v;
  wpl2 = (at["AVWPL", vx, wx] + at["AVLO", vx, wx] Log[v/vx]) vx/v;
  lo2 = at["AVLO", vx, wx] vx/v;
  resc2 = at["STRUV", vx, wx];
  f02 = at["F0", 1 - v, w];
  <|"Class" -> J0, "Name" -> p["Name"] <> " (legs exchanged)",
    "Born" -> alphaS^2 2 (1 - v) f02/Pi,
    "Regular" -> 2 v w (resc2 + (wpl2 - (wpl2 /. w -> 1))/(1 - w) + (lo2 - (lo2 /. w -> 1)) Log[1 - w]/(1 - w))/p["CC"],
    "D0" -> 2 v (wpl2 /. w -> 1)/p["CC"],
    "D1" -> 2 v (lo2 /. w -> 1)/p["CC"],
    "Delta" -> 2 v del2/p["CC"],
    "Prefactor" -> alphaS^3/(Pi s v w), "Basis" -> {"regular", "D0", "D1", "delta"},
    "Scales" -> {Log[muR^2/s], Log[muF^2/s], Log[muFF^2/s]}, "Pieces" -> <|"F02" -> f02, "FDEL2" -> del2, "FVWPL2" -> wpl2, "FVLO2" -> lo2, "FRESC2" -> resc2, "CC" -> p["CC"]|>|>];
'''


def write_wl(path, pieces, cc, lines_info, check_text, verified):
    out = []
    out.append(HEADER % {"verified": verified, "dplus": lines_info["DPLUS"], "check": check_text})
    out.append("PPTransverseFortran = <|")
    out.append('  "Provenance" -> <|"Paper" -> "A. Mukherjee, M. Stratmann, W. Vogelsang, Phys. Rev. D 72 (2005) 034011 [hep-ph/0506315]",')
    out.append('    "Code" -> "me-transpol-ms.f (status 03/11/04, M. Stratmann) with drive-transpol.f; not shipped",')
    out.append('    "Lines" -> <|' + ", ".join('"%s" -> {%d, %d}' % (k, a, b) for k, (a, b) in lines_info.items()) + "|>|>,")
    out.append('  "Variables" -> "v = 1 + t/s, w = -u/(s + t), s the partonic s; muR, muF, muFF the renormalization, initial- and final-state factorization scales; CA, CF, Nf; Born in units of alphaS^2, NLO pieces in units of alphaS^3 (the code multiplies by ALPORD^2 = alphaS^2 and ALPASHO = alphaS)",')
    out.append('  "Pieces" -> "F0 (FBOR, includes Pi/(s v (1-v))), AVDEL (delta(1-w), with Log[muR^2/s]), AVWPL (1/(1-w)_+), AVLO ([Log[1-w]/(1-w)]_+), STRUV (regular), CC (colour normalization)",')
    out.append('  "Map" -> "Born = alphaS^2 2 (1 - v) F0/Pi; {regular, D0, D1, delta} = 2 v w {STRUV, AVWPL(1)/w, AVLO(1)/w, AVDEL/w}/CC in the basis sigma = Born delta(1-w) + alphaS^3/(Pi s v w) regular + alphaS^3/(Pi s v) (D0 [1/(1-w)]_+ + D1 [Log[1-w]/(1-w)]_+ + delta delta(1-w))",')
    out.append('  "Classes" -> <|')
    for n, j in enumerate(CLASSES):
        out.append('    %d -> <|"Name" -> "%s",' % (j, CLASS_NAMES[j]))
        for p in PIECES:
            toks, (a, b) = pieces[j][p]
            expr = to_mathematica(toks)
            out.append('      "%s" -> (* lines %d-%d *)' % (p, a, b))
            out.append(wrap("        " + expr + ","))
        out.append('      "CC" -> %s|>%s' % (cc[j][1], "," if n < len(CLASSES) - 1 else ""))
    out.append("  |>")
    out.append("|>;")
    out.append("")
    out.append("(* The code's pieces of class J0 in the basis of this package (see the header). *)")
    out.append("PPTransverseFortranCoefficients[J0_Integer] := Module[{p = PPTransverseFortran[\"Classes\"][J0]},")
    out.append("  <|\"Class\" -> J0, \"Name\" -> p[\"Name\"],")
    out.append("    \"Born\" -> alphaS^2 2 (1 - v) p[\"F0\"]/Pi,")
    out.append("    \"Regular\" -> 2 v w (p[\"STRUV\"] + (p[\"AVWPL\"] - (p[\"AVWPL\"] /. w -> 1))/(1 - w) + (p[\"AVLO\"] - (p[\"AVLO\"] /. w -> 1)) Log[1 - w]/(1 - w))/p[\"CC\"],")
    out.append("    \"D0\" -> 2 v (p[\"AVWPL\"] /. w -> 1)/p[\"CC\"],")
    out.append("    \"D1\" -> 2 v (p[\"AVLO\"] /. w -> 1)/p[\"CC\"],")
    out.append("    \"Delta\" -> 2 v p[\"AVDEL\"]/p[\"CC\"],")
    out.append("    \"Prefactor\" -> alphaS^3/(Pi s v w), \"Basis\" -> {\"regular\", \"D0\", \"D1\", \"delta\"},")
    out.append("    \"Scales\" -> {Log[muR^2/s], Log[muF^2/s], Log[muFF^2/s]}, \"Pieces\" -> p|>];")
    out.append("")
    out.append(EXCHANGED_WL)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))


def run_wolfram_check(wl_path, pieces, cc, points, tmpdir):
    """Evaluate the .wl at the points with wolframscript; return the values."""
    script = os.path.join(tmpdir, "check_transpol.wls")
    keys = PIECES + ["CC"]
    with open(script, "w", encoding="utf-8") as fh:
        fh.write('SetOptions[$Output, PageWidth -> Infinity];\n')
        fh.write('Get["%s"];\n' % wl_path)
        fh.write("pts = {" + ", ".join(mma_point(pt) for pt in points) + "};\n")
        fh.write('Do[Print[j0, " ", key, " ", i, " ", ToString[N[PPTransverseFortran["Classes"][j0][key] /. Normal[pts[[i]]], 30], InputForm]],\n')
        fh.write("  {j0, %s}, {key, %s}, {i, Length[pts]}];\n" % ("{" + ", ".join(map(str, CLASSES)) + "}", "{" + ", ".join('"%s"' % k for k in keys) + "}"))
    res = subprocess.run(["wolframscript", "-file", script], capture_output=True, text=True, timeout=600)
    if res.returncode != 0:
        raise RuntimeError("wolframscript failed:\n%s\n%s" % (res.stdout, res.stderr))
    values = {}
    for line in res.stdout.splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue
        j0, key, i, val = int(parts[0]), parts[1], int(parts[2]), parts[3]
        val = re.sub(r"`[0-9.]*", "", val).replace("*^", "e")
        values[(j0, key, i)] = mpf(val)
    return values


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("srcdir", help="directory with me-transpol-ms.f, drive-transpol.f, proc*-transpol.f")
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--out", default=os.path.join(here, "..", "pp-transverse-fortran.wl"))
    ap.add_argument("--tmp", default=None, help="directory for temporary files (default: a new temporary directory)")
    ap.add_argument("--no-wolfram", action="store_true", help="skip the wolframscript evaluation")
    args = ap.parse_args()
    if args.tmp is None:
        args.tmp = tempfile.mkdtemp(prefix="transpol-")
    os.makedirs(args.tmp, exist_ok=True)
    out_path = os.path.abspath(args.out)

    main_f = os.path.join(args.srcdir, "me-transpol-ms.f")
    lines = logical_lines(main_f)
    units = split_units(lines)
    pieces, definitions = extract_pieces(units)
    driver = split_units(logical_lines(os.path.join(args.srcdir, "drive-transpol.f")))
    cc, branches = extract_cc(driver["HADRON"])
    print("CC branches of the driver:", branches)
    lines_info = {}
    for name in ["DPLUS", "PRECALC", "FBOR", "AVWPL", "AVDEL", "AVLO", "STRUV"] + ["STRUV%d" % j for j in CLASSES]:
        lines_info[name] = (units[name][0][0], units[name][-1][1])

    ids = sorted({t[1] for j in CLASSES for p in PIECES for t in pieces[j][p][0] if t[0] == "id"})
    print("Fortran symbols used by the pieces:", " ".join(ids))

    # 1. the resolution table against the code's definitions
    verified = ", ".join(verify_definitions(definitions, POINTS))
    print("definitions verified:", verified)

    # 2a. Python values of the Fortran expressions
    py = {}
    for j in CLASSES:
        for p in PIECES:
            for i, pt in enumerate(POINTS, 1):
                py[(j, p, i)] = evaluate(pieces[j][p][0], point_values(pt))
        for i, pt in enumerate(POINTS, 1):
            py[(j, "CC", i)] = evaluate(cc[j][0], point_values(pt))

    # 3. the raw FortranForm files
    proc_report = []
    for j in CLASSES:
        blocks = proc_blocks(os.path.join(args.srcdir, "proc%d-transpol.f" % j))
        expected = {"STRUV"} if j == 7 else {"AVDEL", "AVWPL", "AVLO", "STRUV"}
        if set(blocks) != expected:
            raise ValueError("proc%d: blocks %s, expected %s" % (j, sorted(blocks), sorted(expected)))
        worst = mpf(0)
        for tag, toks in blocks.items():
            for i, pt in enumerate(POINTS, 1):
                a = evaluate(toks, point_values(pt))
                b = py[(j, tag, i)]
                worst = max(worst, abs(a - b) / (1 + abs(b)))
        proc_report.append("proc%d (%s): %s" % (j, ", ".join(sorted(blocks)), mpmath.nstr(worst, 3)))
    print("raw FortranForm files against the routines (max |difference|/(1+|value|)):", "; ".join(proc_report))
    worst_proc = max(mpf(r.split(": ")[-1]) for r in proc_report)
    if worst_proc > mpf(10) ** -30:
        raise ValueError("the proc*.f files differ from the routines")

    # 2b. write, evaluate with Mathematica, compare
    write_wl(out_path, pieces, cc, lines_info, "pending", verified)
    if args.no_wolfram:
        check_text = "Mathematica evaluation skipped (--no-wolfram); Python evaluation and proc*.f comparison done: " + "; ".join(proc_report)
    else:
        mma = run_wolfram_check(out_path, pieces, cc, POINTS, args.tmp)
        worst = mpf(0)
        worst_key = None
        for key, a in py.items():
            if key not in mma:
                raise RuntimeError("missing Mathematica value for %s" % (key,))
            d = abs(mma[key] - a) / (1 + abs(a))
            if d > worst:
                worst, worst_key = d, key
        print("Mathematica file against the Fortran expressions: max |difference|/(1+|value|) = %s at %s" % (mpmath.nstr(worst, 3), worst_key))
        if worst > mpf(10) ** -20:
            raise ValueError("the Mathematica file does not reproduce the Fortran expressions")
        pts_text = "; ".join("(v, w, s, muR^2, muF^2, muFF^2, Nf, CA, CF) = (%s)" % ", ".join(str(pt[k]) for k in ["v", "w", "s", "muR2", "muF2", "muFF2", "Nf", "CA", "CF"]) for pt in POINTS)
        check_text = ("the 26 pieces x %d points, evaluated from the Fortran text in Python/mpmath (40 digits) and from this file in Mathematica (30 digits), "
                      "agree to %s relative (largest |difference|/(1+|value|)); the raw FortranForm files proc*-transpol.f, parsed separately, agree with the routines to %s. "
                      "Points: %s." % (len(POINTS), mpmath.nstr(worst, 2), mpmath.nstr(worst_proc, 2), pts_text))
    write_wl(out_path, pieces, cc, lines_info, check_text, verified)
    print("written", out_path)


if __name__ == "__main__":
    main()
