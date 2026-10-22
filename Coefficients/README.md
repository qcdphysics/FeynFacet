# Coefficients: the final results, ready to use

The eight process tables below hold the final NLO coefficients of every process in this
repository, as plain Mathematica expressions.  They can be used without
FeynCalc, FORM, Tracer, Kira or SubTropica; reading them needs only
`Get`.  They are produced by `Checks/collect_coefficients.wls` from the
results of the process notebooks.  For every coefficient that script checks,
numerically with 40-digit arithmetic at six random physical points (three
with N_c = 3, three with N_c = 5), that the three algebra routes (FeynCalc,
FORM, Tracer) agree and that the form written here equals the computed one:
sampled checks, not symbolic identities.  The physics checks (pole
cancellation, scale dependence, comparison with the literature and the
reference codes) are in the process notebooks.

| File | Content |
|---|---|
| `ProtonProton.wl` | a b → c + X in proton–proton collisions: the 16 channel classes of `Processes/Hadronic`, spin average (UU) and helicity difference (LL), with their Born terms |
| `ProtonProtonTT.wl` | p↑ p↑ → c + X, both incoming partons transversely polarized (A_TT): the coefficient of cos(2Φ) (`TT`) of the classes 5, 6, 7, 11, 12 of `Processes/HadronicTT`, with their Born terms |
| `ProtonProtonTransfer.wl` | a↑ b → c↑ + X, transverse spin transfer from the incoming parton a to the observed parton c: the components `XX` and `NN` of the classes 1, 3, 6, 10, 11, 13, 17 of `Processes/HadronicTransfer`, with their Born terms |
| `LeptonQuark.wl` | ℓ q → q + X, unpolarized and longitudinally polarized |
| `LeptonOffDiagonal.wl` | ℓ q → g + X (`QG`) and ℓ g → q + X (`GQ`), unpolarized and polarized |
| `Collins.wl` | ℓ q → q + X with transverse spin transfer (`Transverse`), the spin average (`Unpolarized`), and the resolved-photon Born table (`ResolvedBorn`) |
| `TMD.wl` | NLO beam and soft-subtracted TMDPDF matching coefficients, and quark/gluon soft functions |
| `PartonicBorn.wl` | the D-dimensional Born squares of all 2 → 2 partonic channels |

Each of these eight files is an association; its `"Description"`, `"Variables"`,
`"Scales"` and `"Basis"` entries say what it contains.

## Auxiliary lepton scale residues

`LeptonScaleResidues.wl` is an additional dataset for restoring separate
PDF, fragmentation and photon-in-lepton scales in the UU/LL lepton
coefficients.  It is an association with six channel keys such as
`"Unpolarized QQ"`; each contains `"PDF"`, `"FF"` and `"Photon"` vectors
in the basis `{regular,D0,D1,delta}`.  These vectors include couplings and
all hard-part prefactors.  They do not have the dimensionless coefficient
schema of the eight process tables above.

`Checks/check_lepton_scales.wls` derives this auxiliary dataset from the
process results and checks it against their pole residues and logarithmic
scale derivatives.  Its SHA-256 and schema are documented separately in
`AUXILIARY-MANIFEST.json`; `MANIFEST.json` and
`Checks/collect_coefficients.wls` cover the original eight process tables.
See [LEPTON-SCALES.txt](../docs/LEPTON-SCALES.txt) for the normalization,
independent-scale module and an executable example.

## The transversely polarized proton–proton files

`ProtonProtonTT.wl` and `ProtonProtonTransfer.wl` have the form of
`ProtonProton.wl`: `pp["Classes"][n]` holds `"Process"`, `"Born"` and the
coefficient lists `{c_regular, c_D0, c_D1, c_δ}` with the same prefactor
α_s³/(π s v w), basis and scale logarithms (LF for both incoming partons, LFF
for the fragmenting parton).  The spin conventions are those of
`docs/TRANSVERSE-PP.md` (HVBM γ5 with four-dimensional spin vectors, MSbar,
transversity kernel δP_qq on every transversely polarized leg, no finite
renormalization).

- `ProtonProtonTT.wl`, key `"TT"`: the coefficient of cos(2Φ) of
  ½[σ(↑↑) − σ(↑↓)], Φ the azimuth of the observed parton about the common
  spin axis; the partonic densities are those whose cross sections multiply
  the transversity distributions δf_a δf_b.
- `ProtonProtonTransfer.wl`, keys `"XX"` and `"NN"` (and `"Born"` →
  `<|"XX", "NN"|>`): σ(s_a, s_c) = XX cos φ_a cos φ_c + NN sin φ_a sin φ_c,
  with the spin of a and of c in the scattering plane (XX) or along its normal
  (NN), for σ(s_c) − σ(−s_c) at the observed parton.  The two components are
  equal for every class in these files (they differ only by evanescent terms
  before renormalization); both are kept.

A class is collected when the required calculation products exist and the
route comparisons pass. The `"Routes"` and `"StepsThreeFour"` metadata identify
which products were used. The collector also supports an algebra-comparison
fallback for a route lacking a finite result; this fallback is distinguished
explicitly from a route with its own reduction and finite assembly. A missing
class appears under `"NotYetCollected"`.

## File checksums and calculation inputs

`MANIFEST.json` binds each coefficient file to its SHA-256 and the relative
calculation-input paths and hashes used by `Checks/collect_coefficients.wls`.
Its `sources` block hashes the distributed scientific source files, and its
`generator` entry identifies the collector. The collector writes this public
format directly.
Intermediate `Results/` trees and machine-specific execution records are not
part of the source release. Running the numbered process steps produces new
local calculation products and receipts; running the collector then binds the
new tables to those products.

To refresh the file/source manifest without recollecting the tables, run
`wolframscript -file Checks/collect_coefficients.wls manifest-only` from the
repository root. Existing calculation-input bindings are retained; a table
whose file hash no longer matches is rejected.

`SHA256SUMS` in the repository root covers the distributed source and data files.
It can be checked from the repository root with `sha256sum -c SHA256SUMS`
(on macOS, `shasum -a 256 -c SHA256SUMS`). Hash agreement checks file integrity;
the process notebooks and comparison scripts provide the physics checks.
The portable alternative is `python3 Checks/file_checksums.py --check`;
use `--write` to regenerate it after intentional source or data changes.

## The form of the hadronic and leptonic coefficients

For a measured parton of momentum p, the partonic cross section
σ = s E dσ/d³p is written in the variables

    v = 1 + t/s,   w = −u/(s + t)        (0 < v, w < 1)

with s, t, u the partonic invariants (`docs/CONVENTIONS.md` gives the
momenta of each process), as

    σ = Born · δ(1−w)
      + P(v, w) c_regular(v, w)
      + P(v, 1) [ c_D0(v) D0 + c_D1(v) D1 + c_δ(v) δ(1−w) ],

    D0 = [1/(1−w)]_+ ,   D1 = [log(1−w)/(1−w)]_+ ,

where P is the `"Prefactor"` (α_s³/(π s v w) for proton–proton) and the
coefficients are the four entries of each list `{c_regular, c_D0, c_D1, c_δ}`
in that order.  The scales enter through

    LR = log(μ_R²/s),   LF = log(μ_F²/s),   LFF = log(μ_FF²/s)

(μ_R renormalization, μ_F incoming parton densities, μ_FF fragmentation;
the lepton files use one common scale, `LR`, and `Collins.wl` the separate
scales described in its `"Scales"` entry).  All coefficients are in the
MSbar scheme; the helicity differences use the HVBM scheme for γ5 with the
standard finite renormalization that restores helicity conservation.  Colour
is SU(N_c): C_A = N_c, C_F = (N_c² − 1)/(2N_c); N_f is the number of light
flavours.

### Translating the manuscript notation

| Code symbol | Manuscript notation | Meaning |
| --- | --- | --- |
| `muFF`, `LFF` in the hadronic calculation | μ_F′, L_F′ | Outgoing fragmentation scale and log(μ_F′²/s) |
| `muI`, `LI` in `Collins.wl` and its calculation | μ_F, L_F | Incoming-quark factorization scale and logarithm |
| `muO`, `LO` in `Collins.wl` and its calculation | μ_F′, L_F′ | Outgoing-quark factorization scale and logarithm |
| `muL`, `LL` in `Collins.wl` and its calculation | μ_L, L_L | Photon-in-lepton factorization scale and logarithm |
| `x` in `TMD.wl` | ξ | Beam-matching convolution variable |

In the TMD tables, `LQ` is a separate external-momentum logarithm,
2 log(xP⁺/ν). Hold it fixed when convolving over the kernel variable `x`;
that variable represents ξ in the paper, rather than the external fraction.

## Example

```wl
pp = Get["Coefficients/ProtonProton.wl"];
c  = pp["Classes"][13]["LL"];                (* q g -> q, helicity difference *)
born = pp["Classes"][13]["Born"]["LL"];
(* the four coefficients at a point, for SU(3), Nf = 5 and all scales equal to the partonic s *)
c /. {Nc -> 3, Nf -> 5, LR -> 0, LF -> 0, LFF -> 0, v -> 3/10, w -> 1/5} // N
```

The plus distributions act on a test function f(w) as
∫₀¹ dw D_n(w) f(w) = ∫₀¹ dw logⁿ(1−w)/(1−w) [f(w) − f(1)];
with a lower limit w₀ > 0 add f(1) logⁿ⁺¹(1−w₀)/(n+1).

## Numerical evaluation without Mathematica: `evaluate.py`

`evaluate.py` is a stand-alone Python module (Python ≥ 3.8, needs only
`mpmath`) that evaluates every coefficient of `ProtonProton.wl`,
`ProtonProtonTT.wl`, `ProtonProtonTransfer.wl`, `LeptonQuark.wl`,
`LeptonOffDiagonal.wl` and `Collins.wl` at arbitrary precision.  It is generated by `Checks/make_evaluator.wls`, which translates
the expressions of the `.wl` files into Python (one function per result,
each coefficient a separate expression); do not edit it by hand, rerun

```sh
wolframscript -file Checks/make_evaluator.wls      # PYTHON=... selects the python used for the final check
```

The coefficient functions return the bare entries
`{c_regular, c_D0, c_D1, c_δ}` as a dict with keys `regular`, `D0`, `D1`,
`delta`, **without** the `"Prefactor"`.  Born terms and prefactors are separate
functions, evaluated at s = 1 with all couplings (α_s, α_em, e_q) set to 1
(in the lepton prefactors C_F = (N_c²−1)/(2N_c), T_R = 1/2):

```python
import sys; sys.path.insert(0, 'Coefficients')
from mpmath import mp
import evaluate as ev
mp.dps = 30
ev.pp(13, 'LL', '3/10', '1/5', Nc=3, Nf=5, LR=0, LF=0, LFF=0)   # {'regular': ..., 'D0': ..., 'D1': ..., 'delta': ...}
ev.born(13, 'LL', '3/10')                     # Born * s / alphaS^2
ev.pp_prefactor('3/10', '1/5')                # 1/(pi v w)
ev.tt(5, '3/10', '1/5', Nc=3, Nf=5)           # q qb -> q', cos(2 Phi) coefficient
ev.tt_born(5, '3/10')
ev.transfer(1, 'XX', '3/10', '1/5')           # q^ q' -> q^, in-plane component
ev.transfer_born(1, 'XX', '3/10')
ev.tt_classes(), ev.transfer_classes()        # the classes present
ev.lepton_quark('Polarized', 0.3, 0.2, LR=0)
ev.lepton_offdiagonal('QG', 'Unpolarized', 0.3, 0.2)
ev.collins('Transverse', 0.3, 0.2, LR=0, LI=0, LO=0, LL=0)
```

Arguments may be numbers, `mpf`, `Fraction` or strings (`'0.3'`, `'3/10'`).
c_D0, c_D1 and c_δ do not depend on w; c_regular needs 0 < w < 1.
From the command line:

```sh
python3 Coefficients/evaluate.py pp 13 LL 0.3 0.2 --digits 30          # options --Nc --Nf --LR --LF --LFF
python3 Coefficients/evaluate.py tt 5 0.3 0.2                          # same options
python3 Coefficients/evaluate.py transfer 1 XX 0.3 0.2                 # same options
python3 Coefficients/evaluate.py lepton_quark Polarized 0.3 0.2 --LR 0.5
python3 Coefficients/evaluate.py lepton_offdiagonal GQ Unpolarized 0.3 0.2
python3 Coefficients/evaluate.py collins Transverse 0.3 0.2 --LI -0.3   # options --LR --LI --LO --LL
python3 Coefficients/evaluate.py --check                                # built-in checks
```

`--check` evaluates every function (32 proton–proton coefficient sets with
their 32 Born terms and the prefactor; 5 TT and 14 transfer coefficient sets
with their Born terms; all lepton coefficients, Born terms and prefactors:
123 functions) at three exact rational points, with nonzero scale logarithms
and N_c = 3 and 4, at `mp.dps = 60`, and compares with 50-digit reference
values that `make_evaluator.wls` computed directly from the `.wl` expressions
and embedded in the file.  Agreement to 40 significant digits is required; the
result is `900 checks passed`, with the largest relative deviation about
1e-49.
The reference values are computed from the same saved `.wl` expressions, so
`--check` tests the translation of those expressions into Python, not the
physics of the coefficients: that is tested in the process notebooks.
