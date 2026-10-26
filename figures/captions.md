# Figure captions

The coefficient plots and FeynArts diagrams are produced by `site/make_figures.wls` (run
`wolframscript -file site/make_figures.wls` from the repository root; the
arguments `diagrams` or `plots` regenerate one part only).  Each figure exists
as SVG and as PNG at twice the nominal size.

The website's two TMD illustrations, `diagrams/tmd-igg` and `diagrams/tmd-sq`,
use TikZ-Feynman instead. Regenerate their SVG and 900-pixel PNG files with
`python3 site/make_tmd_figures.py` (TeX Live and Poppler required).
Editable sources are in `tikz/`. They show one representative real cut each,
with labels $I_{gg}$ and $S_q$; they are not complete NLO diagram sets.
The soft diagram uses double Wilson lines, solid attachment dots and no arrows.

## Feynman-diagram tiles (`diagrams/`)

One diagram per file, drawn with FeynArts `Paint` from the same
`CreateTopologies`/`InsertFields` calls and diagram selections as
`Processes/*/1-Amplitudes.wl`.  The labels of FeynArts' SMQCD model are
renamed to the generic ones: e → ℓ (lepton), u → q, c → q′ (a second quark
flavour), with γ the photon and g the gluon.  In the Collins tile the quark is
written q↑ (transversely polarized).  The caption of every tile is in
`diagrams/index.json` (and the same list in `diagrams/index.wl`), with the fields
`id`, `process`, `group`, `svg`, `png`, `caption`.

| Files | Content |
|---|---|
| `leptonquark-born` | ℓ q → ℓ q, one-photon exchange |
| `leptonquark-vertex` | the one-loop QCD vertex correction (the self-energy graphs on the external legs vanish) |
| `leptonquark-real-1`, `-2` | ℓ q → ℓ q g, gluon from the incoming / outgoing quark |
| `collins-born` | ℓ q↑ → ℓ q↑; the NLO diagrams are those of lepton–quark scattering |
| `leptonoffdiagonal-qg-1`, `-2` | ℓ q → ℓ q g with the gluon observed (ℓ q → g + X) |
| `leptonoffdiagonal-gq-1`, `-2` | ℓ g → ℓ q q̄ with the quark observed (ℓ g → q + X), photon–gluon fusion |
| `pp-class1-born-1` | q q′ → q q′ (class 1), t-channel gluon exchange |
| `pp-class13-born-1` … `-3` | q g → q g (class 13), the three Born diagrams |
| `pp-class15-born-1` … `-4` | g g → g g (class 15), the four Born diagrams (including the four-gluon vertex) |
| `pp-ggggg-real-1`, `-13`, `-15`, `-16` | four of the 25 real-emission diagrams of g g → g g g (class 15); the number is FeynArts' diagram number |
| `tmd-*` | TMD beam and soft functions from the model `Processes/TMD/models/WilsonLineQCD.mod` (Wilson lines W as dashed lines, operator legs O_q, O_g, O_S), with the diagram layouts of `Processes/TMD/models/ShapeData`: lowest order (`bornq`, `borng`), real emission for I_qq, I_qg, I_gq, I_gg (`beamqq`, `beamqg`, `beamgq`, `beamgg`), soft function (`softq`, `softg`), and the scaleless one-loop virtual graph (`virtualq`) |

## Coefficient plots (`plots/`)

Notation as in `Coefficients/README.md`: v = 1 + t/s, w = −u/(s + t), and

σ = Born · δ(1−w) + P(v, w) c_regular(v, w) + P(v, 1) [c_D0(v) D0 + c_D1(v) D1 + c_δ(v) δ(1−w)],

with D0 = [1/(1−w)]₊, D1 = [log(1−w)/(1−w)]₊.  The plotted quantities are the
coefficients c themselves, without the prefactor P (α_s³/(π s v w) for
proton–proton, α_em² α_s C_F e_q²/(π s v w) for lepton–quark).

### `pp-regular-v05` — proton–proton, regular coefficient

c_regular(v = 0.5, w) as a function of w for 0.1 ≤ w ≤ 0.99, for the channel
classes 1 (q q′ → q), 13 (q g → q), 15 (g g → g) and 16 (g g → q) (rows), spin
average UU (left, solid) and helicity difference LL = (σ(++) − σ(+−))/2 (right,
dashed).  The two panels of a row share the vertical scale.  SU(3) color
(N_c = 3), N_f = 5, all scale logarithms zero (μ_R = μ_F = μ_FF = √s, i.e.
LR = LF = LFF = 0).  Source: `Coefficients/ProtonProton.wl`, entry
`"Classes" → n → "UU"/"LL"`, first element.

### `pp-threshold` — proton–proton, threshold coefficients

The coefficients of D0 (left column), D1 (middle) and δ(1−w) (right) as functions
of v for 0.1 ≤ v ≤ 0.9, for the classes 1, 13, 15 and 16 (rows); UU solid, LL
dashed.  N_c = 3, N_f = 5, LR = LF = LFF = 0.  Source: elements 2–4 of the same
lists in `Coefficients/ProtonProton.wl`.  For class 16 the three threshold
coefficients are exactly opposite for UU and LL (c_LL = −c_UU identically, for
any scales), as are the Born terms of g g → q q̄.

### `leptonquark-regular-v05` — lepton–quark, regular coefficient

c_regular(v = 0.5, w) for ℓ q → q + X as a function of w for 0.01 ≤ w ≤ 0.99:
unpolarized (solid), helicity difference LL (dashed; `Coefficients/LeptonQuark.wl`,
`"Polarized"`) and transverse spin transfer, the T_C projection (dotted;
`Coefficients/Collins.wl`, `"Transverse"`).  All scale logarithms zero
(LR = 0 for the unpolarized and helicity cases; LI = LO = LL = 0 for the
transverse case).  The unpolarized coefficients of `Collins.wl` and
`LeptonQuark.wl` are identical (the script checks this).

### `leptonquark-spin-transfer` — lepton–quark, spin-transfer ratios at threshold

Ratio polarized / unpolarized as a function of v of the Born cross sections
(lines) and of the coefficients of D0, D1 and δ(1−w) (open markers at
v = 0.05, 0.10, …, 0.95), for the helicity difference (dashed; ratio of `"Polarized"`
to `"Unpolarized"` in `LeptonQuark.wl`) and transverse spin transfer (dotted;
`"Transverse"` to `"Unpolarized"` in `Collins.wl`).  Scale logarithms zero.
At threshold every coefficient has exactly the Born ratio, (1 − v²)/(1 + v²)
for helicity and 2v/(1 + v²) for transverse spin transfer (the script verifies
this symbolically); the spin dependence of the NLO correction sits entirely in
the regular part.  The same proportionality holds for any common value of the
scale logarithms.
