# Lepton–quark and lepton–gluon scattering with a different tagged parton: ℓ q → g + X and ℓ g → q + X

The two off-diagonal channels of single-inclusive parton production in
lepton–proton collisions, which start at order α² α_s: the gluon is observed
in ℓ q → g + X and the quark in ℓ g → q + X.  Both the spin average and the
double-longitudinal helicity difference of the incoming lepton and parton are
computed.  The results are the coefficients of Hinderer, Schlegel and
Vogelsang, Phys. Rev. D 92 (2015) 014001 [arXiv:1505.06415], Eqs. (27), (28)
and Appendix A, and Phys. Rev. D 96 (2017) 014002 [arXiv:1703.10872],
Eqs. (13), (14) and Appendix A, and the notebooks end by comparing with them.

The structure is that of `Processes/LeptonQuark` (read that README first);
the same momenta, invariants, phase-space families and master integrals are
used.  What differs: there is no Born term and no virtual correction, the real
emission is ℓ q → ℓ g q or ℓ g → ℓ q q̄, the collinear counterterm is the
fragmentation q → g of the Born quark (first channel) or the parton
distribution g → q of the incoming quark (second channel), both built from
the Born process ℓ q → ℓ q, and the photon-in-lepton subtraction uses the
photon Born processes γ q → g q and γ g → q q̄.  No soft gluon can be emitted,
so there are no plus distributions: the coefficients are ordinary functions of
v and w.

## Files

| File | What it does | Kernel time |
|---|---|---|
| `1-Amplitudes.wl` | FeynArts diagrams and amplitudes (Born ℓq → ℓq, the two real-emission processes, the two photon Born processes), converted to plain records | 5 s |
| `2-Algebra-FeynCalc.wl` | Squared matrix elements with FeynCalc; the helicity differences are built from explicit traces, checked against the generic spin average | 10 s |
| `2-Algebra-FORM.wl` | The same with FORM; the programs are in `form/` | 5 s |
| `2-Algebra-Tracer.wl` | The same with Tracer | 10 s |
| `3-Reduce.wl` | Reverse unitarity, Kira reduction (`kira/`), master integrals with SubTropica | 20 s |
| `4-Unpolarized.wl` | Real emission, collinear and photon counterterms, pole cancellation, finite coefficients, comparison with arXiv:1505.06415 | 1.5 min |
| `5-Polarized.wl` | The same for the helicity difference, comparison with arXiv:1703.10872 | 1.5 min |
| `form/` | FORM programs: `common.h`, `Born.frm`, `RealQG.frm`, `RealGQ.frm`, `PhotonQG.frm`, `PhotonGQ.frm` |
| `kira/` | Kira families `real_a_tau` and `real_b_tau` (one per channel) |
| `reference/` | The published coefficients and photon Born terms, transcribed for the final comparison only |
| `Results/` | Created when the notebooks run (not part of the repository) |

## How to run

Open the files in Mathematica in the order 1, 2 (any one of the three
routes), 3, 4, 5, each in a fresh kernel, and evaluate the cells from top to
bottom.  Steps 3, 4 and 5 begin with `route = "FeynCalc"`; change it to
`"FORM"` or `"Tracer"` to continue from another route's step 2.

Headless, from the repository root:

```sh
python Checks/run_process.py Processes/LeptonOffDiagonal --route FeynCalc
python Checks/run_process.py Processes/LeptonOffDiagonal --route FORM
python Checks/run_process.py Processes/LeptonOffDiagonal --route Tracer
```

A completed run means that every check below passed, in particular
`LiteratureDifference -> <|"QG" -> 0, "GQ" -> 0|>` in
`Results/<route>/Unpolarized.wl` and `Polarized.wl`.

## What is checked, and where

| Step | Check |
|---|---|
| 1 | number of diagrams per process (1, 2, 2, 2, 2) |
| 2 (all routes) | Born squares equal 2 e⁴ e_q² [(s²+u²)/t² − ε] and 2 e⁴ e_q² (s²−u²)/t²; the three routes write the same file and agree exactly (Born, real squares, photon squares) |
| 2 (FeynCalc) | the explicit-trace construction with the helicity densities reproduces the generic spin average when the densities are replaced by the averages (checks strings, propagators and signs) |
| 2 (FORM) | the explicit HVBM γ₅ reproduces Tr[γᵃγᵇγᶜγᵈγ₅] = −4i ε^{abcd}; every propagator recognized |
| 3 | the family terms reproduce the matrix element; every integral reduced; two-body volumes equal their Γ-function forms; the angular masters agree with their hypergeometric closed forms and with direct numerical integration |
| 4, 5 | Born hard part equals the textbook form; the photon Born terms equal Eq. (29) of the 2015 paper and Eq. (15) of the 2017 paper; no double pole and no distribution arises; the single pole cancels between real emission, the collinear counterterm and the photon-in-lepton subtraction (table displayed); the finite coefficients equal the published ones, symbolically and at three numerical points |

## Conventions in brief

As in the lepton–quark process: σ(s, v, w) = s E dσ/d³p, t = −s(1−v),
u = −s v w, s4 = s v (1−w); the NLO coefficient is written as
α² e_q² C_F α_s/(π s v w) × R(v, w) (T_R instead of C_F for ℓ g → q + X),
with one scale μ_R for the renormalization and all factorizations, as in the
papers.  The colour factor of the gluon-initiated channel,
Tr[TᵃTᵃ]/(N_c²−1) = 1/2, is kept as the symbol TR.  HVBM γ₅ for the
helicity difference; the polarized kernels ΔP_qg(z) = T_R(2z−1) and
ΔP_{γ/ℓ}(y) = 2 − y enter the incoming legs, the fragmenting quark is
unpolarized, and no finite helicity-restoring term appears (it belongs to the
ℓ q → q + X channel).  See `docs/CONVENTIONS.md` for the complete list.
