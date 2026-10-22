#!/usr/bin/env python3
"""Build the static results website of FeynFacet into site/public/.

    python3 site/build.py                      # writes site/public/
    REPO_URL=https://github.com/<owner>/<repo> python3 site/build.py   # another repository

Plain Python 3 (no packages).  The pages are written from the content below;
figures are taken from site/figures/ (made by site/make_figures.wls), the
coefficient files and the evaluator from Coefficients/.  Formulas are typeset
in the browser with MathJax.  site/publish.sh publishes site/public/ to SITE_BRANCH
(default gh-pages) in SITE_REPO, preserving the source-code branch.
"""
import html
import json
import os
import pathlib
import shutil

ROOT = pathlib.Path(__file__).resolve().parent.parent
SITE = ROOT / "site"
ARTIFACT = os.environ.get("ARTIFACT") == "1"     # embedded entry page without its own <html> skeleton
OUT = SITE / ("artifact" if ARTIFACT else "public")
REPO = os.environ.get("REPO_URL", "https://github.com/qcdphysics/FeynFacet")
BRANCH = os.environ.get("REPO_BRANCH", "main")
PAPER = os.environ.get("PAPER_URL", "https://arxiv.org/abs/2610.31745")
# NOINDEX=1 asks search engines not to index the pages (for test deployments)
NOINDEX = os.environ.get("NOINDEX", "0") == "1"


def src(path):
    """Link to a file (blob) or folder (tree) of the repository on GitHub."""
    kind = "tree" if (ROOT / path).is_dir() else "blob"
    return f"{REPO}/{kind}/{BRANCH}/{path}"


# ---------------------------------------------------------------------------
# page frame
# ---------------------------------------------------------------------------
NAV = [("index.html", "Overview"), ("results.html", "Results"), ("method.html", "Method"),
       ("reproduce.html", "Reproduce"), ("example.html", "Example"), ("references.html", "References"),
       ("paper.html", "Paper")]

STYLE = """
/* a reading column for physicists: serif text, slate-blue accent, light and dark themes */
:root { --ink:#1d2433; --muted:#5a6275; --accent:#3d4f8a; --accent2:#6b4fa0; --rule:#d9dce4;
        --panel:#f5f6f9; --ok:#21734a; --new:#8a3d6b; --bg:#ffffff; --tile:#ffffff; --tileink:#555b66; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) { --ink:#e3e6ee; --muted:#a3a9b8; --accent:#9fb0ea; --accent2:#c3a8ef; --rule:#343a48;
          --panel:#1b1f29; --ok:#6fcf9b; --new:#e39ac5; --bg:#12151c; --tile:#f4f5f8; --tileink:#3d4250; color-scheme:dark; } }
:root[data-theme="dark"] { --ink:#e3e6ee; --muted:#a3a9b8; --accent:#9fb0ea; --accent2:#c3a8ef; --rule:#343a48;
          --panel:#1b1f29; --ok:#6fcf9b; --new:#e39ac5; --bg:#12151c; --tile:#f4f5f8; --tileink:#3d4250; color-scheme:dark; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink); font:17px/1.6 Georgia,"DejaVu Serif",serif; }
header { border-bottom:1px solid var(--rule); }
.bar { max-width:980px; margin:0 auto; padding:14px 20px; display:flex; flex-wrap:wrap; align-items:baseline; gap:6px 22px; }
.brand { font-size:22px; font-weight:bold; color:var(--accent); text-decoration:none; letter-spacing:.3px; }
.brand span { color:var(--accent2); }
.brand img { height:16px; width:auto; display:block; }
.bar { align-items:center; }
h1.logo { margin:34px 0 20px; }
h1.logo + p.lead { font-size:17.5px; }
h1.logo img { display:block; width:min(100%, 360px); height:auto; margin-bottom:12px; }
h1.logo span { display:block; font-size:26px; font-weight:normal; color:var(--muted); text-wrap:balance; }
nav a { color:var(--muted); text-decoration:none; font-size:15px; margin-right:16px; }
nav a.on, nav a:hover { color:var(--accent); }
main { max-width:980px; margin:0 auto; padding:10px 20px 60px; min-width:0; }
.wide { overflow-x:auto; }
h1 { font-size:32px; line-height:1.25; margin:28px 0 6px; }
h2 { font-size:23px; margin:38px 0 10px; padding-bottom:4px; border-bottom:1px solid var(--rule); }
h3 { font-size:19px; margin:24px 0 6px; }
p.lead { font-size:19px; color:var(--muted); margin-top:4px; }
a { color:var(--accent); }
table { border-collapse:collapse; width:100%; margin:12px 0 18px; font-size:15px; }
th, td { border-bottom:1px solid var(--rule); padding:7px 8px; text-align:left; vertical-align:top; }
th { font-weight:normal; color:var(--muted); font-size:14px; }
tr.group td { font-weight:bold; background:var(--panel); }
.tag { display:inline-block; font-size:12.5px; padding:1px 8px; border-radius:10px; border:1px solid currentColor; white-space:nowrap; }
.tag.new { color:var(--new); } .tag.known { color:var(--ok); }
.panel { background:var(--panel); border-radius:6px; padding:14px 18px; margin:16px 0; }
.standard { counter-reset:n; list-style:none; padding:0; }
.standard li { counter-increment:n; margin:10px 0 10px 38px; position:relative; }
.standard li::before { content:counter(n); position:absolute; left:-38px; top:0; width:26px; height:26px; border-radius:50%;
  background:var(--accent); color:var(--bg); text-align:center; font:bold 14px/26px sans-serif; }
.tiles { display:grid; grid-template-columns:repeat(auto-fill,minmax(150px,1fr)); gap:14px; margin:14px 0; }
.tiles figure { margin:0; background:var(--tile); border:1px solid var(--rule); border-radius:6px; padding:6px; }
.tiles .diagram { aspect-ratio:6 / 5; }
.tiles img { width:100%; height:100%; object-fit:contain; object-position:center; display:block; }
.tiles.tmd { grid-template-columns:repeat(2,minmax(0,240px)); }
.tiles figcaption { font-size:12.5px; color:var(--tileink); text-align:center; margin-top:4px; }
figure.plot { margin:22px 0; } figure.plot img { width:100%; background:var(--tile); border-radius:4px; }
figure.plot figcaption { font-size:14.5px; color:var(--muted); margin-top:6px; }
pre, code { font-family:"DejaVu Sans Mono",Menlo,monospace; font-size:14px; }
pre { background:var(--panel); padding:12px 14px; border-radius:6px; overflow-x:auto; }
.refs li { margin:8px 0; } .refs .note { color:var(--muted); font-size:14.5px; }
.prompt { position:relative; } .prompt pre { white-space:pre-wrap; padding-top:40px; }
.prompt button { position:absolute; top:8px; right:8px; font:13px sans-serif; padding:4px 12px; border-radius:4px;
  border:1px solid var(--rule); background:var(--bg); color:var(--accent); cursor:pointer; }
.steps { padding-left:22px; } .steps li { margin:6px 0; }
footer { border-top:1px solid var(--rule); color:var(--muted); font-size:14px; }
footer .bar { display:block; }
"""

KATEX = """
<script>window.MathJax = {tex: {inlineMath: [['$', '$']], displayMath: [['$$', '$$']]}, svg: {fontCache: 'global'}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>
"""


def page(filename, title, body, depth=0):
    up = "../" * depth
    links = " ".join(f'<a href="{up}{href}"{" class=on" if href == filename else ""}>{label}</a>'
                     for href, label in NAV)
    links += f' <a href="{REPO}">Code</a>'
    robots = '<meta name="robots" content="noindex, nofollow">\n' if NOINDEX else ""
    head = robots + f"""<title>{"FeynFacet" if filename == "index.html" else html.escape(title) + " — FeynFacet"}</title>
<style>{STYLE}</style>{KATEX}"""
    bare = ARTIFACT and filename == "index.html"     # the embedding publisher supplies the document skeleton
    doc = (head + "\n" if bare else f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
{head}</head>
<body>
""") + f"""<header><div class="bar"><a class="brand" href="{up}index.html"><img src="{up}figures/logo/feynfacet-logo.svg" alt="FeynFacet"></a><nav>{links}</nav></div></header>
<main>
{body}
</main>
<footer><div class="bar">FeynFacet — Z.-B. Kang, W. Vogelsang, C. Zhang.
Source code and coefficients: <a href="{REPO}">{REPO}</a>.
<a href="{up}LICENSE.txt">MIT license</a>.</div></footer>
""" + ("" if bare else "</body></html>\n")
    doc = doc.replace('<a href="https://', '<a target="_blank" rel="noopener" href="https://')
    target = OUT / filename
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(doc, encoding="utf-8")


# ---------------------------------------------------------------------------
# figures (produced by site/make_figures.wls)
# ---------------------------------------------------------------------------
def figure_index():
    for name in ("figures/index.json", "figures/diagrams/index.json"):
        path = SITE / name
        if path.exists():
            try:
                return json.loads(path.read_text())
            except json.JSONDecodeError:
                pass
    return None


def tiles(prefix, depth=0, names=None, css_class=""):
    """Diagram tiles whose file name starts with prefix (site/figures/diagrams/<prefix>-*.svg), or the given names."""
    folder = SITE / "figures" / "diagrams"
    if not folder.exists():
        return ""
    captions = {}
    index = figure_index()
    if isinstance(index, list):
        for item in index:
            if isinstance(item, dict) and ("svg" in item or "file" in item):
                captions[pathlib.Path(item.get("svg", item.get("file"))).name] = item.get("caption", "")
    elif isinstance(index, dict):
        for key, value in index.items():
            captions[pathlib.Path(key).name] = value if isinstance(value, str) else value.get("caption", "")
    files = [folder / (n + ".svg") for n in names] if names else sorted(folder.glob(prefix + "*.svg"))
    files = [f for f in files if f.exists()]
    if not files:
        return ""
    up = "../" * depth
    # the two TikZ cut diagrams of the TMD page: short captions of their own
    tmd_labels = {
        "tmd-igg.svg": "Gluon TMD PDF: the beam function I_gg at NLO (real emission).",
        "tmd-sq.svg": "Soft function S_q at NLO: a real gluon between the Wilson lines along n and n̄.",
    }
    cells = []
    for f in files:
        alt = tmd_labels.get(f.name, captions.get(f.name, f.stem))
        caption = tmd_labels.get(f.name, captions.get(f.name, ""))
        label = f'<figcaption>{html.escape(caption)}</figcaption>' if caption else ""
        cells.append(f'<figure><div class="diagram"><img src="{up}figures/diagrams/{f.name}" '
                     f'alt="{html.escape(alt)}" loading="lazy"></div>{label}</figure>')
    return f'<div class="tiles {html.escape(css_class)}">{"".join(cells)}</div>'


def plot(name, caption, depth=0):
    folder = SITE / "figures"
    for ext in ("svg", "png"):
        matches = sorted(folder.rglob(f"{name}.{ext}"))
        if matches:
            rel = matches[0].relative_to(SITE).as_posix()
            return (f'<figure class="plot"><img src="{"../" * depth}{rel}" alt="{html.escape(caption)}" loading="lazy">'
                    f'<figcaption>{caption}</figcaption></figure>')
    return ""


# ---------------------------------------------------------------------------
# content
# ---------------------------------------------------------------------------
NEW = '<span class="tag new">new</span>'
KNOWN = '<span class="tag known">known, rederived</span>'

RESULTS = [
    # (page, title, process, status, compared with)
    ("collins", "Transverse spin transfer in lepton–quark scattering",
     r"$\ell\, q^{\uparrow} \to q^{\uparrow} + X$ at NLO; resolved-photon Born channels", NEW,
     "new"),
    ("transfer", "Transverse spin transfer in proton–proton collisions",
     r"$a^{\uparrow} b \to c^{\uparrow} + X$ at NLO, seven channel classes", NEW,
     r"related to the double-transverse result by a generalized crossing symmetry found in this calculation"),
    ("pp", "Proton–proton collisions: the 16 channel classes",
     r"$a\,b \to c + X$, spin average and double-longitudinal helicity difference", KNOWN,
     "Aversa, Chiappetta, Greco, Guillet (1989); Jäger, Schäfer, Stratmann, Vogelsang (2003)"),
    ("tt", "Double-transverse proton–proton collisions",
     r"$a^{\uparrow} b^{\uparrow} \to c + X$ ($A_{TT}$), five channel classes", KNOWN,
     "Mukherjee, Stratmann, Vogelsang (2005), their Fortran code"),
    ("leptonquark", "Lepton–quark scattering",
     r"$\ell\, q \to q + X$, unpolarized and double-longitudinal", KNOWN,
     "Hinderer, Schlegel, Vogelsang (2015, 2017)"),
    ("offdiagonal", "Lepton–quark and lepton–gluon, other tagged parton",
     r"$\ell\, q \to g + X$ and $\ell\, g \to q + X$", KNOWN,
     "Hinderer, Schlegel, Vogelsang (2015, 2017)"),
    ("tmd", "TMD beam functions and soft function",
     r"quark and gluon TMD PDFs at NLO, exponential rapidity regulator", KNOWN,
     "Luo, Yang, Zhu, Zhu (2021)"),
    ("born", "Born table",
     r"all $2\to2$ partonic channels in $D$ dimensions, unpolarized and helicity difference", KNOWN,
     "Ellis, Stirling, Webber (Table 7.1); helicity and crossing relations"),
]


def results_table(depth=0):
    up = "../" * depth
    rows = "".join(f'<tr><td><a href="{up}results/{p}.html">{t}</a></td><td>{proc}</td><td>{st}</td><td>{cmp}</td></tr>'
                   for p, t, proc, st, cmp in RESULTS)
    return ('<div class="wide"><table><tr><th>Result</th><th>Process</th><th>Status</th><th>Compared with</th></tr>'
            + rows + '</table></div>')


def overview():
    body = f"""
<h1 class="logo"><img src="figures/logo/feynfacet-logo.svg" alt="FeynFacet"><span>One Feynman-diagram calculation, many facets</span></h1>
<p class="lead">FeynFacet derives perturbative results in quantum field theory from Feynman diagrams.  Every
calculation follows one chain of steps: diagrams, spin and color algebra, integration-by-parts reduction, master
integrals, renormalization, and the factorization or matching that the observable requires.  A new observable
reuses the chain and changes only what it needs.  The chain is written as executable Mathematica notebooks that a
physicist can read cell by cell, with three independent implementations of the spin and color algebra.</p>

<p>The calculations released here give the NLO partonic hard functions for single-inclusive hadron and jet production in
proton–proton and lepton–proton collisions, for unpolarized and longitudinally polarized beams, the double-transverse spin asymmetry
$A_{{TT}}$ in proton–proton collisions, NLO TMDPDF matching coefficients and soft functions, and two new results: the NLO
transverse spin transfer in lepton–quark scattering, which enters transverse-spin measurements at the
Electron-Ion Collider, and the NLO transverse spin transfer in proton–proton collisions.  Every published result is rederived, not
copied, and compared with the literature: identically where the published coefficients are analytic expressions,
otherwise with the reference codes at many points to their numerical precision.</p>
{tiles("", 0, ["leptonquark-born", "leptonquark-vertex", "leptonquark-real-1", "leptonoffdiagonal-gq-1", "pp-class15-born-4", "pp-ggggg-real-13", "tmd-igg", "tmd-sq"])}
{getting_started()}
<p>Then reproduce a calculation (<a href="reproduce.html">Reproduce</a>) or let the agent compute a new process
(<a href="example.html">Example</a>).</p>

<h2>The standard</h2>
<p>A result is shown here only when it meets all four requirements:</p>
<ol class="standard">
<li><b>Three algebra routes.</b> The squared matrix elements are computed independently with
<a href="https://feyncalc.github.io">FeynCalc</a>, <a href="https://form-dev.github.io/form/">FORM</a> and
<a href="https://library.wolfram.com/infocenter/MathSource/2987/">Tracer</a>.  For proton–proton collisions the three
squares must agree before any integration (in exact rational arithmetic at two random rational points of all variables,
every scalar integral given its own value); the integration and assembly act on each route's squares, and the finite
coefficients are compared across routes.  For the lepton–proton and TMD calculations the finite coefficients of the
three routes agree symbolically.</li>
<li><b>Published benchmarks and crossing.</b> For known results, the finite coefficients agree with the published
coefficients or the reference codes: identically where analytic expressions are available, otherwise to the
precision of the reference at many points.  The new transverse spin transfer in proton–proton collisions is tied to
the published double-transverse result by a generalized crossing symmetry found in this calculation (verified
symbolically at Born level and at exact test points at NLO).</li>
<li><b>Every step visible.</b> Diagram generation (<a href="https://www.feynarts.de/">FeynArts</a>), traces, reverse unitarity, <a href="https://kira.hepforge.org/">Kira</a> reduction, master
integrals (<a href="https://subtropi.ca/">SubTropica</a>), threshold expansion, renormalization and factorization are cells of the notebooks; every
check is a visible <code>Require</code> that stops the calculation if it fails (pole cancellation, Born tables,
scale dependence, and gauge invariance where it is tested).</li>
<li><b>Usable without the tools.</b> The final coefficients are plain expressions in
<a href="{src('Coefficients/README.md')}"><code>Coefficients/</code></a>, with a stand-alone evaluator that
checks itself against 50-digit reference values.</li>
</ol>

<h2>The results</h2>
{results_table(0)}

<h2>How it was made</h2>
<p>The calculations were carried out by AI agents working with physicists: the agents wrote and ran the
process programs and used physics checks to locate and correct inconsistencies; the physicists specified
the problems and reviewed the results.  The paper describes the workflow{f' (<a href="{PAPER}">arXiv</a>)' if PAPER else ''};
these pages present the calculation itself.</p>
{COPY_SCRIPT}
"""
    page("index.html", "Overview", body)


def results_overview():
    body = f"""<h1>Results</h1>
<p class="lead">Each page states what is computed, how, which checks it passes, and where the coefficients are.</p>
{results_table(0)}
<p>All coefficients are in the MSbar scheme; helicity differences use the HVBM scheme for $\\gamma_5$ with
the standard finite renormalization that restores helicity conservation for incoming quarks; the
transverse results use the HVBM scheme with four-dimensional spin vectors and need no finite
renormalization.</p>"""
    page("results.html", "Results", body)


PP_INVARIANTS = r"$s=(k_a+k_b)^2$, $t=(k_a-p)^2$, $u=(k_b-p)^2$ for incoming partons of momenta $k_a$, $k_b$"
LEPTON_INVARIANTS = r"$s=(l+k)^2$, $t=(k-p)^2$, $u=(l-p)^2$ for the lepton momentum $l$ and the incoming parton momentum $k$"
PP_SCALES = (r"$L_R=\ln(\mu_R^2/s)$, $L_F=\ln(\mu_F^2/s)$ for the incoming partons and $L_{FF}=\ln(\mu_{FF}^2/s)$ "
             r"for the fragmenting one ($\mu_{FF}$ is the $\mu_F'$ of the paper)")
LEPTON_SCALES = r"$L_R=\ln(\mu^2/s)$, with one common scale $\mu_R=\mu_F=\mu_F'=\mu_L=\mu$"
COLLINS_SCALES = r"the separate scales $\mu_I$, $\mu_O$, $\mu_L$ of the incoming quark, the outgoing quark and the photon-in-lepton subtraction ($\mu_F$, $\mu_F'$, $\mu_L$ in the paper; logarithms $L_I$, $L_O$, $L_L$ in <code>Collins.wl</code>); the unpolarized entry of the same file, a check of the pipeline, uses one common scale as for lepton–quark scattering"


def distribution_form(invariants, scales):
    return r"""
<p>For a measured parton of momentum $p$ the partonic cross section $\sigma = s\,E\,d\sigma/d^3p$ is written in
$v = 1 + t/s$ and $w = -u/(s+t)$, with """ + invariants + r""" and $0 &lt; v,w &lt; 1$:</p>
$$\sigma = \sigma_{\rm Born}\,\delta(1-w) + P(v,w)\,c_{\rm reg}(v,w) + P(v,1)\big[c_{0}(v)\,\big[\tfrac{1}{1-w}\big]_+ + c_{1}(v)\,\big[\tfrac{\ln(1-w)}{1-w}\big]_+ + c_{\delta}(v)\,\delta(1-w)\big],$$
<p>Here $\sigma_{\rm Born}$ is the lowest-order term, and $c_{\rm reg}$, $c_0$, $c_1$, $c_\delta$ are the coefficients
$R$, $B$, $C$, $A$ of the paper, with $P$ its normalization $\mathcal N$.  The scales enter through """ + scales + ".</p>"


def process_page(key, title, body):
    page(f"results/{key}.html", title, f"<h1>{title}</h1>\n{body}", depth=1)


PP_CLASSES = [
    (1, r"q\,q'\to q"), (2, r"q\,q'\to g"), (3, r"q\,\bar q'\to q"), (4, r"q\,\bar q'\to g"),
    (5, r"q\,\bar q\to q'"), (6, r"q\,q\to q"), (7, r"q\,q\to g"), (8, r"q\,g\to q'"),
    (9, r"q\,g\to \bar q'"), (10, r"q\,g\to \bar q"), (11, r"q\,\bar q\to q"), (12, r"q\,\bar q\to g"),
    (13, r"q\,g\to q"), (14, r"q\,g\to g"), (15, r"g\,g\to g"), (16, r"g\,g\to q"),
]


def pp_page():
    rows = "".join(f"<tr><td>{n}</td><td>${p}$</td><td>agrees</td><td>agrees</td><td>agree (2 random exact points)</td></tr>"
                   for n, p in PP_CLASSES)
    body = f"""
<p class="lead">Single-inclusive parton production $a\\,b\\to c+X$ at NLO in the 16 tagged channel classes,
for the spin average (UU) and the helicity difference of the incoming partons (LL).</p>
{tiles("pp-", 1)}
<h2>What is computed</h2>
{distribution_form(PP_INVARIANTS, PP_SCALES)}
<p>Here $P = \\alpha_s^3/(\\pi s v w)$.  The fragmenting parton $c$ is unpolarized; the helicity difference is
$\\tfrac12[\\sigma(++)-\\sigma(+-)]$, in the HVBM scheme with the finite term
$-(\\alpha_s/2\\pi)\\Delta f_{{qq}}\\otimes\\sigma_{{LL}}^{{(0)}}$ on each polarized incoming quark or antiquark
(see <a href="../method.html">Method</a>).  Classes 2, 4, 7, 8, 9, 10 start at order $\\alpha_s^3$ (no Born term).</p>

<h2>The calculation</h2>
<p>Step 1 generates the Born, virtual (one-loop) and real-emission diagrams with FeynArts.  Step 2 squares them in
the three routes; for $g\\,g\\to g\\,g\\,g$ (25 diagrams, 325 products), every amplitude is simplified once with the transversality conditions of the gluon
polarizations, color is separated and each color network evaluated once, only the products $i\\le j$ are
formed, gluon polarization sums are applied one gluon at a time, and the products run in parallel.
Step 3 maps the real emission onto six reverse-unitarity families and the virtual graphs onto three box
families, reduces them with Kira, and evaluates the masters with SubTropica.  Steps 4 and 5 expand at
threshold, renormalize, subtract the collinear singularities of both incoming partons and the
fragmenting one, and compare with the reference codes.</p>

<h2>Checks</h2>
<div class="wide"><table><tr><th>Class</th><th>Process</th><th>UU vs Aversa et al.</th><th>LL vs Jäger et al.</th><th>FeynCalc = FORM = Tracer</th></tr>{rows}</table></div>
<p>The reference codes were evaluated at nine points per class; all four coefficients agree to the precision of
the double-precision reference values ($10^{{-15}}$ to $2\\cdot10^{{-7}}$).  Inside the notebooks:
double and single poles cancel, the Born squares equal the Born table, for selected classes the real-emission
squares do not depend on the gluon reference vectors, and the squares of the three routes agree (exact arithmetic at two random rational points).</p>
{plot("pp-regular-v05", "The regular coefficient $c_{\\rm reg}(v=0.5,w)$ for the classes 1 ($q\\,q'\\to q$), 13 ($q\\,g\\to q$), 15 ($g\\,g\\to g$) and 16 ($g\\,g\\to q$); spin average (left) and helicity difference (right).  $N_c=3$, $N_f=5$, all scales $\\sqrt s$.", 1)}
{plot("pp-threshold", "The coefficients of $[1/(1-w)]_+$, $[\\ln(1-w)/(1-w)]_+$ and $\\delta(1-w)$ as functions of $v$, same classes; spin average solid, helicity difference dashed.  For $g\\,g\\to q$ the three are exactly opposite for the two spin cases, as are the Born terms.", 1)}

<h2>Files</h2>
<p>Notebooks: <a href="{src('Processes/Hadronic')}">Processes/Hadronic</a>.
Coefficients: <a href="{src('Coefficients/ProtonProton.wl')}">Coefficients/ProtonProton.wl</a>.</p>
<pre>python3 Coefficients/evaluate.py pp 13 LL 0.3 0.2 --digits 30</pre>
"""
    process_page("pp", "Proton–proton collisions", body)


TT_CLASSES = [(5, r"q\,\bar q\to q'", "yes"), (6, r"q\,q\to q", "yes"), (7, r"q\,q\to g", "no"),
              (11, r"q\,\bar q\to q", "yes"), (12, r"q\,\bar q\to g", "yes")]


def tt_page():
    rows = "".join(
        f"<tr><td>{n}</td><td>${p}$</td><td>{b}</td><td>agree (2 random exact points)</td>"
        + "<td>agrees (all four coefficients, each route"
        + ")</td><td>0</td></tr>"
        for n, p, b in TT_CLASSES)
    body = f"""
<p class="lead">The NLO hard part of the double-spin asymmetry $A_{{TT}}$ of single-inclusive hadron production in
collisions of two transversely polarized protons, $a^\\uparrow b^\\uparrow\\to c+X$.  {KNOWN}</p>
<h2>What is computed</h2>
<p>Each incoming quark or antiquark carries the chiral-odd density $\\tfrac12\\not k_a\\gamma_5\\not s_a$ (and likewise
for $b$), whose partonic
cross section multiplies its transversity distribution $\\delta f$; the observed parton is unpolarized.  The cross
section depends on the azimuth $\\Phi$ of the observed parton about the common spin axis as $\\cos 2\\Phi$, and the hard
part is the coefficient of $\\cos 2\\Phi$ in $\\tfrac12[\\sigma(\\uparrow\\uparrow)-\\sigma(\\uparrow\\downarrow)]$, written in the distribution form below.</p>
{distribution_form(PP_INVARIANTS, PP_SCALES)}
<p>Here $P = \\alpha_s^3/(\\pi s v w)$; $L_F$ refers to both transversity densities.
Only five of the sixteen channel classes contribute: a gluon has no transversity, and the two helicity flips
must sit on one fermion line of the cut diagram (annihilation, or the exchange interference of identical
quarks).  Class 7 starts at order $\\alpha_s^3$.</p>

<h2>The calculation</h2>
<p>The diagrams, the squaring, the reduction and the counterterm formulas are those of the unpolarized and
longitudinal classes.  The spin vectors are written in the basis of the scattering plane,
$s=\\cos\\varphi\\,e_X+\\sin\\varphi\\,e_N$, with $e_X$ in the plane along the transverse momentum of the observed
parton and $e_N$ the normal, and the squares are computed as the components $XX$ and $NN$ (the mixed ones vanish
by parity, which is checked).  The normal component of the unobserved momentum is reduced through the Gram
matrix of the three external momenta, the normal component of the loop momentum by symmetric integration.
In the HVBM scheme the symmetric part $(XX+NN)/2$ is evanescent but not zero in $D$ dimensions; the hard part is
$(XX-NN)/2$, and the symmetric part is carried along as a check: its finite remainder vanishes in every
class, so no finite renormalization is needed.  Counterterms: the transversity kernel
$P^T_{{qq}}(z)=C_F[2z/(1-z)_+ + \\tfrac32\\delta(1-z)]$ ($\\delta P_{{qq}}$ in the code) on both incoming legs, the unpolarized kernels for the
fragmenting parton, with the $D$-dimensional TT Born table.</p>

<h2>Checks</h2>
<div class="wide"><table><tr><th>Class</th><th>Process</th><th>Born</th><th>Squares: FeynCalc = FORM = Tracer</th>
<th>vs Mukherjee–Stratmann–Vogelsang code</th><th>Symmetric remainder</th></tr>{rows}</table></div>
<p>The reference is the Fortran code of Mukherjee, Stratmann and Vogelsang (2005), transcribed verbatim into
Mathematica; its slots hold one half of the $\\cos 2\\Phi$ coefficient, a factor fixed on the Born term and
confirmed by the renormalization- and factorization-scale terms. The four coefficients are evaluated with
30-digit arithmetic at six rational kinematic points and $N_f=4,5$ (12 SU(3) points), using distinct fixed
renormalization, initial-state and fragmentation scales. The required SU(3) scaled residual is below
$10^{{-20}}$. Additional $C_A=2,5$ diagnostics, with $C_F=(C_A^2-1)/(2C_A)$, are reported separately; the
reference code's general-color validity is not assumed. The Born terms equal the lowest-order cross sections of Soffer, Stratmann and Vogelsang (2002).  Inside the
notebooks: poles cancel, the scale derivatives are those of the evolution equations, the real-emission squares
are gauge invariant, and the squares of the three routes agree (exact arithmetic at two random rational points).</p>
<p>All five classes are reduced and assembled separately in each of the three algebra routes.
Their squared matrix elements agree at two exact-rational points, and their finite coefficients are compared at six sampled points.</p>

<h2>Files</h2>
<p>Notebooks: <a href="{src('Processes/HadronicTT')}">Processes/HadronicTT</a>; conventions:
<a href="{src('docs/TRANSVERSE-PP.md')}">docs/TRANSVERSE-PP.md</a>; reference code:
<a href="{src('Reference/pp-transverse-fortran.wl')}">Reference/pp-transverse-fortran.wl</a>.
Coefficients: <a href="{src('Coefficients/ProtonProtonTT.wl')}">Coefficients/ProtonProtonTT.wl</a>
(classes 5, 6, 7, 11, 12).</p>
<pre>python3 Coefficients/evaluate.py tt 5 0.3 0.2 --digits 30</pre>
"""
    process_page("tt", "Double-transverse proton–proton collisions", body)


TRANSFER_CLASSES = [(1, r"q\,q'\to q", "yes", 5), (3, r"q\,\bar q'\to q", "yes", 5), (6, r"q\,q\to q", "yes", 11),
                    (10, r"q\,g\to \bar q", "no", 7), (11, r"q\,\bar q\to q", "yes", 11), (13, r"q\,g\to q", "yes", 12),
                    (17, r"\bar q\,q\to q", "yes", 6)]

def transfer_page():
    def row(n, p, b, partner):
        crossing = "Born and virtual identical, real squares at 3 exact points; NLO coefficients at 3 exact points"
        return (f"<tr><td>{n}</td><td>${p}$</td><td>{b}</td><td>{partner}</td><td>agree (2 random exact points)</td>"
                f"<td>{crossing}</td></tr>")
    rows = "".join(row(*c) for c in TRANSFER_CLASSES)
    body = f"""
<p class="lead">The NLO hard parts of the transverse spin transfer in proton–proton collisions,
$a^\\uparrow b\\to c^\\uparrow+X$: a transversely polarized parton of the polarized proton scatters and the transverse
spin of the observed parton is measured, as for a transversely polarized hyperon, the spin asymmetry of a
single hadron or a dihadron inside a jet, or one- and two-point energy correlators in the jet.  The results are related
to the double-transverse cross sections by a <a href="#crossing">generalized crossing symmetry</a> found in this
calculation.  {NEW}</p>
<h2>What is computed</h2>
<p>The incoming parton $a$ carries the density $\\tfrac12\\not k_a\\gamma_5\\not s_a$, the observed parton $c$ is
projected with $\\not p\\,\\gamma_5\\not s_c$ ($\\sigma(s_c)-\\sigma(-s_c)$), and the cross section is
$\\sigma = XX\\,\\cos\\varphi_a\\cos\\varphi_c + NN\\,\\sin\\varphi_a\\sin\\varphi_c$, with both spins in the scattering
plane ($XX$) or along its normal ($NN$).  Both components are written in the distribution form below.</p>
{distribution_form(PP_INVARIANTS, PP_SCALES)}
<p>Here $P = \\alpha_s^3/(\\pi s v w)$; their finite coefficients turn out equal in every class (they differ in $D$
dimensions by evanescent terms), so the transfer depends on $\\cos(\\varphi_a-\\varphi_c)$ only.  Seven classes
contribute; class 10 starts at order $\\alpha_s^3$.  Classes 1 ($q\\,q'\\to q$) and 3 ($q\\,\\bar q'\\to q$) have the
same Born term but different NLO coefficients.</p>

<h2>The calculation</h2>
<p>As for the double-transverse process, with the transversity kernel $P^T_{{qq}}$ on the polarized incoming parton
and on the fragmenting quark or antiquark, the unpolarized kernels on the other incoming parton, and a $D$-dimensional Born table
of all seven transfer channels.  The squares are computed with FeynCalc, FORM and Tracer.</p>

<h2 id="crossing">A generalized crossing symmetry</h2>
<p>This is a new result, not previously available.  In the course of the calculation we found that it is nevertheless
not independent of the double-transverse cross sections: the two are related by a generalized crossing symmetry.
As in ordinary crossing, the transversely polarized incoming parton $b$ of $a^\\uparrow b^\\uparrow\\to c+X$ is moved
to the final state and the observed parton $c$ to the initial state.  The transverse polarization is carried along: the
transversity of the incoming parton $b$ becomes the measured transverse spin of the outgoing $\\bar b$, with a definite
transformation of the spin components and of the color and polarization averages:</p>
$$a^\\uparrow b^\\uparrow\\to c+X\\quad\\longleftrightarrow\\quad a^\\uparrow\\bar c\\to\\bar b^\\uparrow+X,\\qquad
s\\leftrightarrow t,\\quad u\\to u,\\quad v\\to-\\frac{{v}}{{1-v}},\\quad w\\to w.$$
<p>The symmetry is exact: it holds in $D$ dimensions and before any integration, for the squared matrix elements of the
Born, virtual (one-loop) and real-emission contributions alike.  With $\\mathcal C$ the crossing of these unintegrated expressions and
$T^{{\\rm in}}$ the components of the double-transverse squares, the spin components transform as</p>
$$T_{{XX}} = -\\varsigma\\, r\\,\\mathcal C\\big[T_{{XX}}^{{\\rm in}}\\big],\\qquad T_{{NN}} = +\\varsigma\\, r\\,\\mathcal C\\big[T_{{NN}}^{{\\rm in}}\\big],$$
<p>with $(\\varsigma,r)=(1,1)$ when the unpolarized incoming parton $\\bar c$ of the transfer process is a quark or antiquark
and $(-1,[(D-2)C_F]^{{-1}})$ when it is a gluon (its color and polarization average), kept exact in $D$ until the Laurent
expansion.</p>
<p>Here $\\varsigma$ is the fermion crossing sign.  The spin factor of a quark ($\\not k$, $\\not p$, or
$\\not k\\gamma_5\\not s$) is linear in its momentum, which crossing reverses, so each crossed quark brings a factor $-1$;
the polarization sum of a gluon is even in its momentum and brings none.  The crossed parton $b$ is always a quark, since
a gluon has no transversity, so $\\varsigma=(-1)^2=1$ when $\\bar c$ is a quark or antiquark and $\\varsigma=-1$ when it is
a gluon.  The factor $r$ is the change of the color and polarization averages.  The opposite signs of the two
components come from the spin axes: the crossed in-plane axis reverses the azimuth, $\\varphi\\to\\pi-\\varphi$, which gives
$-1$ for $XX$, while the normal enters both spin vectors and its reflection cancels in $NN$.  Consequently the double-transverse projection
$(T^{{\\rm in}}_{{XX}}-T^{{\\rm in}}_{{NN}})/2$ becomes the transfer projection $(T_{{XX}}+T_{{NN}})/2$.  The symmetry also connects a
property of each calculation to the other: the symmetric part $(T^{{\\rm in}}_{{XX}}+T^{{\\rm in}}_{{NN}})/2$ of the
double-transverse squares is evanescent and leaves no finite remainder, which after crossing is the statement that the
transfer components $XX$ and $NN$ have equal finite coefficients, so that the transfer depends on
$\\cos(\\varphi_a-\\varphi_c)$ only.</p>
<p>Step 5 of the transfer workflow (<code>5-Crossing.wl</code>) uses the symmetry as a second, independent derivation of
all seven classes from the double-transverse squares:</p>
<ul>
<li>the crossed Born terms and Born table equal those of the direct calculation identically;</li>
<li>the crossed one-loop correction equals the direct one identically after reduction to master integrals;</li>
<li>the crossed real-emission squares equal the direct ones at three random exact rational points (they are too large
for a symbolic difference);</li>
<li>integrated in the transfer region and combined with the counterterms of the transfer process, the crossed squares
give pole-free finite coefficients equal to those of the direct calculation in exact arithmetic at three rational points
of $s$, $v$, $w$ and $N_f$, with the three scales symbolic.</li>
</ul>
<p>Direct and crossed squares pass through the same scalar reduction and finite-subtraction implementation.  Since the
double-transverse classes agree with the code of Mukherjee, Stratmann and Vogelsang, the symmetry ties the new transfer
results to an independent published calculation.</p>

<h2>Checks</h2>
<div class="wide"><table><tr><th>Class</th><th>Process</th><th>Born</th><th>Crossed from TT class</th>
<th>FeynCalc = FORM = Tracer</th><th>Crossing</th></tr>{rows}</table></div>
<p>All seven classes are reduced and assembled separately in each of the three algebra routes.
The squared matrix elements are compared at two exact-rational points, and the finite coefficients at six sampled points.
Inside the notebooks, poles cancel, the scale derivatives are those of the evolution equations, the
mixed components $XN$, $NX$ vanish, the real emission is gauge invariant, and the Born table equals the
lowest-order Collins hard factors of Yuan (Phys. Rev. Lett. 100 (2008) 032003, Eq. (6)) for the channels given there.</p>

<h2>Files</h2>
<p>Notebooks: <a href="{src('Processes/HadronicTransfer')}">Processes/HadronicTransfer</a> (the crossing is
<code>5-Crossing.wl</code>); conventions: <a href="{src('docs/TRANSVERSE-PP.md')}">docs/TRANSVERSE-PP.md</a>.
Coefficients: <a href="{src('Coefficients/ProtonProtonTransfer.wl')}">Coefficients/ProtonProtonTransfer.wl</a>
(classes 1, 3, 6, 10, 11, 13, 17).</p>
<pre>python3 Coefficients/evaluate.py transfer 1 XX 0.3 0.2 --digits 30</pre>
"""
    process_page("transfer", "Transverse spin transfer in proton–proton collisions", body)


def lepton_pages():
    process_page("leptonquark", "Lepton–quark scattering", f"""
<p class="lead">$\\ell\\,q\\to q+X$ at NLO in QCD, unpolarized and with the helicity difference of the incoming
lepton and quark.</p>
{tiles("leptonquark-", 1)}
<h2>What is computed</h2>{distribution_form(LEPTON_INVARIANTS, LEPTON_SCALES)}
<p>Here $P=\\alpha_{{\\rm em}}^2\\alpha_s C_F e_q^2/(\\pi s v w)$, $\\sigma_{{\\rm Born}}=2\\alpha_{{\\rm em}}^2e_q^2(1+v^2)/[s v(1-v)^2]$
unpolarized and $2\\alpha_{{\\rm em}}^2e_q^2(1+v)/[s v(1-v)]$ for the helicity difference.  The photon-in-lepton
subtraction is included.</p>
<h2>Checks</h2>
<p>Exact symbolic agreement with Hinderer, Schlegel and Vogelsang, Phys. Rev. D 92 (2015) 014001 and
Phys. Rev. D 96 (2017) 014002, in all three algebra routes; pole cancellation; the vector form factor.</p>
<h2>Files</h2>
<p><a href="{src('Processes/LeptonQuark')}">Processes/LeptonQuark</a> ·
<a href="{src('Coefficients/LeptonQuark.wl')}">Coefficients/LeptonQuark.wl</a></p>""")
    process_page("offdiagonal", "Off-diagonal lepton channels", f"""
<p class="lead">$\\ell\\,q\\to g+X$ and $\\ell\\,g\\to q+X$, which start at order $\\alpha_{{\\rm em}}^2\\alpha_s$.</p>
{tiles("leptonoffdiagonal-", 1)}
<p>There is no Born term, so the finite coefficients contain only the regular part $c_{{\\rm reg}}(v,w)$, with
$P=\\alpha_{{\\rm em}}^2\\alpha_s C_F e_q^2/(\\pi s v w)$ for $\\ell\\,q\\to g$ and
$P=\\alpha_{{\\rm em}}^2\\alpha_s T_R e_q^2/(\\pi s v w)$ for $\\ell\\,g\\to q$ (invariants and scales as for
<a href="leptonquark.html">lepton–quark scattering</a>).  The collinear
counterterm is the fragmentation $q\\to g$ of the Born quark or the parton distribution $g\\to q$, and the
photon-in-lepton subtraction uses the photon Born processes.</p>
<h2>Checks</h2>
<p>Exact agreement with Hinderer, Schlegel and Vogelsang (2015, 2017), unpolarized and polarized, in all
three algebra routes.</p>
<h2>Files</h2>
<p><a href="{src('Processes/LeptonOffDiagonal')}">Processes/LeptonOffDiagonal</a> ·
<a href="{src('Coefficients/LeptonOffDiagonal.wl')}">Coefficients/LeptonOffDiagonal.wl</a></p>""")
    process_page("collins", "Transverse spin transfer in lepton–quark scattering", f"""
<p class="lead">The NLO hard part of $\\ell\\,q^\\uparrow\\to q^\\uparrow+X$: a transversely polarized quark
scatters off the lepton and the transverse spin of the outgoing quark is measured.  At the Electron-Ion Collider it enters
measurements of quark transversity through the transverse spin of the observed parton, as for a transversely polarized
hyperon, the spin asymmetry of a single hadron or a dihadron inside a jet, or one- and two-point energy correlators in
the jet.  {NEW}</p>
{tiles("collins-", 1)}
<h2>What is computed</h2>
<p>$T=\\sigma(k^\\uparrow\\to p^\\uparrow)-\\sigma(k^\\uparrow\\to p^\\downarrow)$ from the densities
$\\tfrac12\\not k\\,\\gamma_5\\not s_k$ and $\\not p\\,\\gamma_5\\not s_p$, along the normal to the scattering
plane ($T_{{NN}}$) and in the plane ($T_{{XX}}$); the Collins projection is $T_C=\\tfrac12(T_{{XX}}+T_{{NN}})$, with
Born value $T_C^{{(0)}}=-4e^4e_q^2\\,su/t^2$ (the evanescent terms of $T_{{XX}}$ and $T_{{NN}}$ cancel).  The
transversity kernel $P^T_{{qq}}$ acts on both quark legs.  The resolved-photon contributions are given at Born level
for all parton–parton channels.</p>
{distribution_form(LEPTON_INVARIANTS, COLLINS_SCALES)}
<p>Here $\\sigma_{{\\rm Born}}=4\\alpha_{{\\rm em}}^2e_q^2/[s(1-v)^2]$ and $P=\\alpha_{{\\rm em}}^2e_q^2C_F\\alpha_s/(\\pi s v w)$.</p>
<h2>Checks</h2>
<p>The three algebra routes agree exactly; poles cancel; the $\\mu_I$ and $\\mu_O$ dependence equals the
convolutions with the transversity kernel, the $\\mu_L$ dependence cancels against the photon distribution convolved
with the photon Born term, and there is no $\\mu_R$ dependence at this order; the coefficients are those given in the paper.  The same
files reproduce the unpolarized coefficient of Hinderer, Schlegel and Vogelsang as a check of the pipeline.</p>
{plot("leptonquark-regular-v05", "The regular coefficient $c_{\\rm reg}(v=0.5,w)$ of $\\ell\\,q\\to q+X$: unpolarized, transverse spin transfer ($T_C$ projection) and helicity difference; all scales $\\sqrt s$.", 1)}
{plot("leptonquark-spin-transfer", "Ratios of the spin-dependent to the unpolarized coefficients as functions of $v$.  At threshold every ratio (of the $[1/(1-w)]_+$, $[\\ln(1-w)/(1-w)]_+$ and $\\delta(1-w)$ coefficients) equals the Born ratio exactly, $2v/(1+v^2)$ for the transverse transfer and $(1-v^2)/(1+v^2)$ for the helicity difference: the spin dependence beyond the Born term sits in the regular part.", 1)}
<h2>Files</h2>
<p><a href="{src('Processes/Collins')}">Processes/Collins</a> ·
<a href="{src('Coefficients/Collins.wl')}">Coefficients/Collins.wl</a></p>""")
    process_page("tmd", "TMD beam functions and soft function", f"""
<p class="lead">The NLO matching coefficients $I_{{qq}}, I_{{qg}}, I_{{gq}}, I_{{gg}}$ of the TMD beam functions,
quark and gluon soft functions, and soft-subtracted TMDPDF matching coefficients, with the exponential rapidity regulator.</p>
{tiles("", 1, ["tmd-igg", "tmd-sq"], css_class="tmd")}
<p>Wilson lines and operator insertions are fields of a FeynArts model, so the operator matrix elements
are generated like cross sections.  The Fourier transforms become Euler integrals reduced with Kira and
evaluated with SubTropica.</p>
<h2>Checks</h2>
<p>Exact agreement with the ancillary files of Luo, Yang, Zhu and Zhu, JHEP 06 (2021) 115; Ward identities of
the generated currents; the dependence on the rapidity scale cancels between the beam and soft functions; the
three algebra routes agree.</p>
<h2>Files</h2>
<p><a href="{src('Processes/TMD')}">Processes/TMD</a> · <a href="{src('Coefficients/TMD.wl')}">Coefficients/TMD.wl</a></p>""")
    process_page("born", "The Born table", f"""
<p class="lead">The $D$-dimensional squared matrix elements of all $2\\to2$ partonic channels, spin averaged
and with the helicity difference of the incoming partons (HVBM scheme).</p>
<p>They enter the collinear counterterms of every NLO calculation here, with their order-$\\epsilon$ terms.
Checked against the table of Ellis, Stirling and Webber (<i>QCD and Collider Physics</i>, Table 7.1), helicity
conservation and crossing; the FORM and Tracer tables equal
the FeynCalc one for any $D$ and $N_c$.</p>
<p><a href="{src('Processes/PartonicBorn')}">Processes/PartonicBorn</a> ·
<a href="{src('Coefficients/PartonicBorn.wl')}">Coefficients/PartonicBorn.wl</a></p>""")


def method_page():
    body = r"""
<h1>Method</h1>
<p class="lead">Every process follows the same sequence of steps, one notebook each; the processes differ in the
integral families and counterterms they need.  The squaring is done three times, independently, by <a href="https://feyncalc.github.io">FeynCalc</a>,
<a href="https://form-dev.github.io/form/">FORM</a> and
<a href="https://library.wolfram.com/infocenter/MathSource/2987/">Tracer</a>; from the integrals onward the same reduction and assembly act on each workflow's squares, and the
results are required to agree.</p>

<div class="wide"><table>
<tr><th>Step</th><th>Notebook</th><th>Tools</th><th>Output</th></tr>
<tr><td>1</td><td>Diagrams and amplitudes</td><td>FeynArts</td><td>plain records: coefficient, open fermion chains, closed loops</td></tr>
<tr><td>2</td><td>Spin and color algebra</td><td>FeynCalc &nbsp;|&nbsp; FORM + <code>color.h</code> &nbsp;|&nbsp; Tracer</td><td>Born table, virtual (one-loop) interference, real-emission squares</td></tr>
<tr><td>3</td><td>Integrals</td><td>reverse unitarity; Kira (IBP reduction); SubTropica (master integrals)</td><td>real and virtual contributions as Laurent series in $\epsilon$</td></tr>
<tr><td>4–5</td><td>Assembly</td><td>Mathematica</td><td>distributions, MSbar renormalization and factorization, finite coefficients, comparison with the references</td></tr>
</table></div>

<h2>1. Diagrams</h2>
<p><a href="https://www.feynarts.de/">FeynArts</a> generates the Born, virtual (one-loop) and real-emission diagrams of each process (closed quark loops counted by
$N_f$).  They are written once as plain records, a commuting coefficient times open fermion chains and closed
fermion loops, from which each algebra workflow computes its squares on its own; in the proton–proton processes the
FORM and Tracer notebooks read the FeynCalc squares only at their end, to compare with them.</p>

<h2>2. Three independent algebra workflows</h2>
<div class="wide"><table>
<tr><th>Workflow</th><th>Dirac algebra</th><th>Color</th></tr>
<tr><td>FeynCalc</td><td>traces and contractions in Mathematica, each amplitude simplified once</td><td>SU(N) algebra of FeynCalc, once per color network</td></tr>
<tr><td>FORM</td><td>one FORM job per product $M_i M_j^*$; Mathematica writes the diagrams as color structures times Lorentz parts</td><td>the FORM package <code>color.h</code>, before the Dirac algebra</td></tr>
<tr><td>Tracer</td><td>the Tracer package, in a kernel without FeynCalc</td><td>an explicit SU(N) reduction written in the notebook (traces, Fierz identity, closed loops)</td></tr>
</table></div>
<p>Spin densities, for an incoming quark of momentum $k$ and an outgoing quark of momentum $p$: $\tfrac12\not k$ for
an unpolarized incoming quark and $\not p$ for the spin sum of an outgoing one, $\tfrac12\not k\gamma_5$ for the
helicity difference, $\tfrac12\not k\gamma_5\not s_T$ for a transversely polarized incoming quark, and the projection
$\not p\gamma_5\not s_T'$ on an observed transversely polarized quark (antiquarks with the charge-conjugated
matrices).  Gluons have physical polarization sums with a reference vector $n$,
$-g^{\mu\nu}+(k^\mu n^\nu+n^\mu k^\nu)/(k\cdot n)$, averaged with $1/(D-2)$ for an unpolarized incoming gluon, and
$i\epsilon^{\mu\nu\alpha\beta}k_\alpha n_\beta/(2k\cdot n)$, with $n$ lightlike, for the helicity difference of an
incoming gluon; for selected processes gauge invariance is checked by changing the reference vectors.  $\gamma_5$
follows the HVBM scheme (<a href="https://doi.org/10.1016/0550-3213(72)90279-9">'t Hooft and Veltman</a>; <a href="https://doi.org/10.1007/BF01609069">Breitenlohner and Maison</a>): it anticommutes with the four-dimensional and commutes with the $(D-4)$-dimensional Dirac matrices,
and the hatted components of loop and unobserved momenta are kept.</p>
<p>In the proton–proton processes, at the end of step 2 each workflow's squares are compared with the FeynCalc ones in
exact rational arithmetic at two random points of all variables, with a separate random value assigned to each scalar
integral.  These are sampled comparisons of linear combinations of integral coefficients, not coefficientwise symbolic
identities.  In every process the finite coefficients of the three workflows are compared at the end: symbolically
for the lepton–proton and TMD calculations, at test points for proton–proton collisions.</p>

<h2>3. Integrals</h2>
<p>Two programs do the integration.  <a href="https://kira.hepforge.org/">Kira</a> uses integration-by-parts
identities to reduce both virtual loop integrals and, through reverse unitarity, real-emission phase-space
integrals to a set of master integrals per family.  <a href="https://subtropi.ca/">SubTropica</a>
then evaluates each master analytically from its Feynman-parameter (Euler-integral) representation, with the
integrator HyperIntica, and returns its Laurent series in $\epsilon$ to the order needed.  In the collinear processes
the masters that are Gamma or hypergeometric functions are also compared with their closed forms.</p>
<h3>Real emission: reverse unitarity</h3>
<p>At fixed observed momentum $p$, the real emission is integrated over the unobserved massless pair with total
momentum $P_X$ and invariant mass $P_X^2=s_4=sv(1-w)$, where $P_X=k_a+k_b-p$ for parton–parton and $P_X=l+k-p$ for
lepton–parton scattering:
$$\int d\Phi_2^{(D)}\,F(r)=(2\pi)^{2-D}\int d^Dr\,\delta_+(r^2)\,\delta_+\big((P_X-r)^2\big)\,F(r),\qquad
\delta_+(k^2)=\theta(k^0)\,\delta(k^2).$$
Reverse unitarity writes each on-shell constraint as a cut propagator,
$$\delta(X)=\frac{1}{2\pi i}\Big(\frac{1}{X-i0}-\frac{1}{X+i0}\Big),$$
so that the phase-space integrals become cut integrals, reduced by integration-by-parts identities like loop
integrals; the positive-energy conditions are kept in the definition of the cut masters.  Partial fractions rewrite
the integrand as a sum of terms with at most two linearly independent ordinary denominators, so each integral family
has the two cut propagators and two ordinary denominators (six families in proton–proton collisions, three in
lepton–quark scattering); integrals with a nonpositive power of a cut propagator vanish.  Kira generates the
reduction identities at a symbolic integration dimension $d$, evaluated at $d=D$ for ordinary integrals and at the
shifted dimensions below for integrals with evanescent numerators.  The masters are the two-body phase-space volume
and angular integrals with both ordinary denominators; up to prefactors, the angular masters are values of one
dimension-dependent function at different cross ratios.  Their singular terms take the form
$(1-w)^{-1-\lambda\epsilon}B_\lambda(w,\epsilon)$ with $B_\lambda$ regular at $w=1$, and are expanded in
distributions in step 4.</p>
<h3>Virtual corrections</h3>
<p>The virtual integrals reduce to an off-shell massless bubble in lepton–quark scattering, and to massless bubbles
and boxes in the hadronic processes, written in the three box families (s,t), (s,u), (t,u).  The bubbles are Euler
integrals; the boxes are evaluated with SubTropica from their Feynman-parameter representation in the
Euclidean region and continued with the Feynman prescription to the physical region $s>0$ before the real part of the
interference is taken.  In the hadronic processes the loop normalization $2/(4\pi)^{D/2}$ is derived from the
FeynArts loop measure and checked twice: the coefficient of $1/\epsilon^2$ must equal the universal
$-\frac{\alpha_s}{2\pi}\sum_i C_i$ times the Born term, with $C_i=C_F$ for a quark and $C_A$ for a gluon summed over
the external partons $i$, and the poles must cancel against the real emission.</p>
<h3>HVBM: the dimension shift</h3>
<p>In the HVBM scheme a $D$-dimensional momentum $q$ is split into a four-dimensional part and a
$(D-4)$-dimensional part $\widehat q$; the external momenta are four-dimensional.  After the traces, the squared
matrix elements therefore contain powers of $\widehat q^{\,2}$ for the loop momentum and for the unobserved momenta of
the real emission, in addition to the ordinary propagators.  Such evanescent numerators are removed by raising the
dimension of the integral.  For a scalar integrand $F(q)$ with four-dimensional external momenta, and with
$\mu_q^2=-\widehat q^{\,2}$ in the mostly-minus metric, the normalized one-loop identity is</p>
$$\int\frac{d^Dq}{i\pi^{D/2}}\,(\mu_q^2)^n F(q)=(-\epsilon)_n\int\frac{d^{D+2n}q}{i\pi^{(D+2n)/2}}\,F(q),
\qquad D=4-2\epsilon,$$
<p>where $n$ is a nonnegative integer and $(-\epsilon)_n=\Gamma(n-\epsilon)/\Gamma(-\epsilon)$.  For an insertion
$(\widehat q^{\,2})^n$ the coefficient of the shifted integral is $(-1)^n(-\epsilon)_n$, equal to $+\epsilon$ for
$n=1$.  For the two-body phase-space measure of the real emission above, the corresponding factor is
$(4\pi)^n(-1)^n(-\epsilon)_n$.</p>
<p>Only the integral moves to $D+2n$ dimensions: its IBP coefficients and master integrals are evaluated at
$d=D+2n$, while $D=4-2\epsilon$ is kept in the trace coefficients, spin projectors and spin averages.  Shifting those
as well would change, for instance, the $\delta(1-w)$ coefficient of the helicity difference by a rational term
proportional to the Born.  The evanescent numerator factors are retained because their products with poles can
remain finite.</p>
<h3>Transverse spins</h3>
<p>The momenta $k_a$, $k_b$ of the two incoming partons and $p$ of the observed parton span the scattering plane
(the code calls the incoming momenta $p_a$, $p_b$).  Spin vectors are decomposed on in-plane axes $X$ and on the
normal $N\propto\epsilon^{\mu k_a k_b p}$, with $N^2=-1$; the squares are computed as the components XX, NN, XN, NX,
and the mixed ones must vanish.</p>
<p>The normal components of the unobserved momenta and of the loop momentum enter only through $(N\cdot q)^2$, where
$q$ is the momentum in question.  Write $q=q_\parallel+c_N N+\widehat q$, with $q_\parallel$ in the scattering plane,
$c_N$ a number and $\widehat q$ the $(D-4)$-dimensional part introduced above.  Since $N$ is orthogonal to the plane and
to the $(D-4)$-dimensional directions, $q^2=q_\parallel^2-c_N^2+\widehat q^{\,2}$, so that</p>
$$(N\cdot q)^2=c_N^2=q_\parallel^2-q^2+\widehat q^{\,2}.$$
<p>The in-plane part is fixed by the three scalar products of $q$ with $k_a$, $k_b$ and $p$.  Writing $p_1=k_a$,
$p_2=k_b$, $p_3=p$ and $G_{ij}=p_i\cdot p_j$ for the $3\times3$ table of their scalar products (the Gram matrix),</p>
$$q_\parallel^2=\sum_{i,j=1}^{3}(q\cdot p_i)\,\big(G^{-1}\big)_{ij}\,(q\cdot p_j).$$
<p>For the unobserved partons of the real emission $q^2=0$; for the loop momentum, $q^2$ and the products $q\cdot p_i$
are combinations of the propagators, and $\widehat q^{\,2}$ is handled by the dimension shift above.  Odd powers of
$N\cdot q$ vanish: in the phase-space integral by reflection through the scattering plane, and in the loop because
its tensor integrals can only be built from $k_a$, $k_b$ and $p$.  No new integrals appear.</p>

<h2>4. Distributions and counterterms</h2>
<p>Near threshold the real emission is expanded in distributions,
$$(1-w)^{-1-\lambda\epsilon} = -\frac{\delta(1-w)}{\lambda\epsilon} + \Big[\frac{1}{1-w}\Big]_+ - \lambda\epsilon\Big[\frac{\ln(1-w)}{1-w}\Big]_+ + \dots,$$
with $\lambda$ the exponent of the singular term, the virtual correction is renormalized in MSbar, and the collinear
singularities are removed with the $D$-dimensional Born table (its order-$\epsilon$ terms included) and the
unpolarized, helicity or transversity kernels, e.g. $P^T_{qq}=C_F\big[\frac{2z}{(1-z)_+}+\frac32\delta(1-z)\big]$
($\delta P_{qq}$ in the code).  For the helicity difference the finite term
$-(\alpha_s/2\pi)\,\Delta f_{qq}\otimes\sigma_{LL}^{(0)}$, with $\Delta f_{qq}(z)=-4C_F(1-z)$, is added on each
longitudinally polarized incoming quark or antiquark; it restores the helicity convention of the polarized MSbar
parton distributions.  Transversity needs no such term.  Double and single poles must cancel exactly.</p>

<h2>5. Checks</h2>
<ul>
<li>the three workflows: equal finite coefficients (with separate reduction and finite assembly for each route), and
for proton–proton collisions equal squares (exact arithmetic at random rational points);</li>
<li>for selected processes, gauge invariance of the real-emission squares under a change of the gluon reference
vectors, and Ward identities of the TMD currents;</li>
<li>the universal double pole of the virtual correction, and the cancellation of all poles;</li>
<li>the scale dependence: the derivatives with respect to the factorization scales equal the Born terms convolved with
the leading-order splitting kernels;</li>
<li>helicity conservation in the pure-annihilation channel $q_i\bar q_i\to q_k+X$, $k\ne i$:
$\sigma_{LL}=-\sigma$ in each distribution sector;</li>
<li>in the TMD calculation, the cancellation of the rapidity-scale dependence between the beam and soft functions;</li>
<li>independent results: the Fortran codes of Aversa et al. (UU), Jäger, Schäfer, Stratmann, Vogelsang (LL) and
Mukherjee, Stratmann, Vogelsang (TT), and the published unpolarized and longitudinal lepton–proton and TMD results;</li>
<li>the generalized crossing symmetry between the double-transverse cross sections and the transverse spin transfer.</li>
</ul>
<p>Every check is a <code>Require</code> in the notebook that performs it: a run that does not pass it stops.</p>
"""
    page("method.html", "Method", body)


def reproduce_page():
    body = f"""
<h1>Reproduce</h1>
<h2>Use the coefficients</h2>
<p>Only Mathematica or Python is needed.</p>
<pre>pp = Get["Coefficients/ProtonProton.wl"];
pp["Classes"][13]["LL"] /. {{Nc -> 3, Nf -> 5, LR -> 0, LF -> 0, LFF -> 0, v -> 3/10, w -> 1/5}} // N</pre>
<pre>pip install mpmath                                          # the only requirement
python3 Coefficients/evaluate.py --check                    # built-in 50-digit checks
python3 Coefficients/evaluate.py pp 13 LL 0.3 0.2 --digits 30</pre>
<h2>Run the calculation</h2>
{getting_started("h3", "With an agent (recommended)")}
<p>When the check reports READY, open <code>Start.nb</code> and evaluate it, or ask the agent to run a process.</p>
<h3>Manually</h3>
<p>Install Mathematica with FeynCalc and FeynArts, FORM, Tracer, Kira (with Fermat) and SubTropica (with polymake)
(<a href="{src('docs/INSTALL.md')}">installation notes</a>; on Apple silicon, Kira needs the
<a href="{src('docs/KIRA-CORRECTNESS.md')}">patched build</a>).  Then check the toolchain; the check performs small
exact calculations with each tool, including eight Kira reductions, and reports what to fix:</p>
<pre>python3 feynfacet.py check                # add --all-routes for FORM and Tracer</pre>
<p>When it reports READY, open <code>Start.nb</code> and evaluate it: it loads the package and runs the lepton–quark
process end to end.  Every process can also be run from a terminal:</p>
<pre>python3 Checks/run_process.py Processes/LeptonQuark
python3 Checks/run_process.py Processes/PartonicBorn                         # Born table, read by the proton–proton processes
python3 Checks/run_process.py Processes/Hadronic --set class=13              # FeynCalc route first
python3 Checks/run_process.py Processes/Hadronic --set class=13 --route FORM # compared with the FeynCalc squares</pre>
<p>Each step runs in a fresh kernel; the results of a step are files that the next one reads.  The process READMEs give
the order for each process (the crossing check of the spin transfer, for instance, needs the double-transverse classes).
The headless runner records evaluated cells and SHA-256 hashes of calculation sources and recorded inputs and outputs in
<code>Processes/&lt;process&gt;/Results/&lt;step&gt;[-&lt;choices&gt;].receipt.json</code>. A receipt is marked complete after the last cell succeeds.</p>
{COPY_SCRIPT}
"""
    page("reproduce.html", "Reproduce", body)


EXAMPLE_PROMPT = """1. Learn FeynFacet

The FeynFacet repository is supplied in this working directory.
Read and follow its AGENTS.md, README.md and relevant
documentation, and study the existing processes: their conventions,
shared routines, execution procedure and validation checks.
Run python3 feynfacet.py check --all-routes. If a tool is missing or the
check does not report READY, install or fix it following docs/INSTALL.md
(docs/INSTALL-macOS.md on a Mac), and ask me before anything that needs
an administrator password.
Identify the existing process implementations closest to the task
below and use them as templates, explaining your choice.

2. Specify the process and observable

Derive and execute the analytical LO and NLO QCD calculation for
inclusive prompt-photon production, p + p -> photon + X, for both
unpolarized and double-longitudinally polarized beams.

For this example, compute the q + qbar -> photon + X direct partonic
channel, without photon isolation. Treat X inclusively and include
every contributing NLO final state for this incoming channel.

Derive the result independently: do not search online for, or use,
published results for this process.

3. Choose and apply a workflow

Use one supported FeynFacet algebra workflow: FeynCalc, FORM or Tracer.
Identify the chosen route. Follow FeynFacet's cross-section normalization,
kinematic variables, dimensional and spin conventions.

Derive the coefficients from amplitudes using reverse unitarity,
dimensional shifts where applicable, Kira reduction and SubTropica
integration, following the existing workflow.

4. Renormalization and factorization

Determine and implement the renormalization and factorization terms
required to obtain finite hard coefficients, following FeynFacet's
conventions. Identify and derive any process-specific ingredients
missing from the existing examples.

For this non-isolated photon process, include photon-fragmentation
factorization in the MSbar scheme. Explain the subtraction conventions
and retain the dependence on the relevant scales.

5. Execute, validate and provide results

Work in a new directory, preserving the supplied sources and results.
Develop and execute the calculation. Apply all internal consistency
and validation checks applicable to this calculation and the chosen
workflow, following FeynFacet's existing implementations.

Provide readable, executable sources and notebooks, complete finite
analytic coefficients including their distribution and scale dependence,
the conventions used, and a summary of the checks and their results."""

SETUP_PROMPT = """Set up FeynFacet on this machine. Install the tools it needs, following
docs/INSTALL.md (docs/INSTALL-macOS.md on a Mac; on Apple silicon, build
Kira with tools/install_kira.py as described there), and run
python3 feynfacet.py check --all-routes until it reports READY.
Ask me before anything that needs an administrator password."""

COMPARISON_PROMPT = """Compare these results with Gordon–Vogelsang, Phys. Rev. D 48 (1993) 3136:
unpolarized LO Eq. (12) and polarized LO Eq. (14); NLO Eq. (36) with
Appendices D/E. Translate normalization and factorization conventions as
necessary; for the polarized initial state, Gordon–Vogelsang use the MS_P
scheme of their Eqs. (32)–(33), so report the scheme-transformation term
separately. Report exact symbolic differences. If you cannot open the paper,
stop and ask me for the PDF instead of working from memory."""

COPY_SCRIPT = """<script>
function copyPrompt(button) {
  var text = button.parentElement.querySelector('pre').innerText;
  function done() { button.textContent = 'Copied'; setTimeout(function () { button.textContent = 'Copy'; }, 1500); }
  function select() {
    var range = document.createRange(); range.selectNodeContents(button.parentElement.querySelector('pre'));
    var sel = window.getSelection(); sel.removeAllRanges(); sel.addRange(range);
    button.textContent = 'Selected: press Ctrl/Cmd+C';
  }
  try { navigator.clipboard.writeText(text).then(done, select); } catch (e) { select(); }
}
</script>"""


def prompt_box(text):
    return f'<div class="prompt"><button type="button" onclick="copyPrompt(this)">Copy</button><pre>{html.escape(text)}</pre></div>'


def getting_started(heading="h2", title="Getting started"):
    return f"""
<{heading}>{title}</{heading}>
<p>An AI coding agent installs and checks the tools for you; you do not need to set them up one by one.</p>
<ol class="steps">
<li>Install Mathematica and clone the <a href="{REPO}">repository</a>.</li>
<li>Start the agent in the repository folder, with permission to run Mathematica and the external tools without
approving each command: Codex with <code>codex --approve-for-me</code> (or “Approve for me” under
<code>/permissions</code> at the start of a session); Claude Code in auto mode (Shift+Tab, or
<code>claude --permission-mode auto</code>).  The tools write to their user and temporary directories, which a
restrictive sandbox blocks.</li>
<li>Paste the setup prompt below.  The agent installs FeynCalc, FeynArts, FORM, Tracer, Kira, Fermat and SubTropica
and runs the toolchain check until it reports READY.</li>
</ol>
{prompt_box(SETUP_PROMPT)}
<p>The <a href="{src('docs/INSTALL.md')}">installation notes</a> are what the agent follows.  Mathematica and its
license stay with you, and the agent asks before anything that needs your administrator password; on a Mac this
includes Homebrew and Apple's Command Line Tools for Xcode.  On Apple silicon the agent builds a patched Kira, since the
Kira 3.1 release build can silently return wrong reductions there
(<a href="{src('docs/KIRA-CORRECTNESS.md')}">details</a>); this one-time build takes a few minutes.</p>"""


def example_page():
    body = f"""
<h1>Example: a new process with an AI agent</h1>
<p class="lead">FeynFacet is meant to be extended.  This example uses a coding agent (for instance Codex or Claude Code) to compute a
process that is not in the repository, prompt-photon production in proton–proton collisions, $p+p\\to\\gamma+X$, and to check the
result against a published calculation.</p>

{getting_started("h2", "Before you start")}
<p>If the tools are not yet installed, you can also go straight to the prompt below: its first step installs whatever
is missing.</p>

<h2>The prompt</h2>
<p>We tested this example with Codex using GPT-6-Astra at maximum reasoning effort, and with Claude Code using
Claude Opus 5.5.</p>
<p>Steps 1, 3 and 5 form a template that applies to any new calculation; steps 2 and 4 describe this process and the
renormalization and factorization it needs.  To compute another process, change those two steps.</p>
{prompt_box(EXAMPLE_PROMPT)}

<h2>Optional: compare with the literature</h2>
<p>The annihilation channel was computed at NLO by Gordon and Vogelsang,
<a href="https://doi.org/10.1103/PhysRevD.48.3136"><i>Polarized and unpolarized prompt photon production beyond the
leading order</i>, Phys. Rev. D 48 (1993) 3136</a>, who give the complete coefficients in their appendices.  After the
calculation, a second prompt asks the agent to compare:</p>
{prompt_box(COMPARISON_PROMPT)}
<p>The agent needs to read the paper itself.  If it cannot open the article online, for example without a journal
subscription, download the PDF yourself and save it where the agent can read it, such as the repository folder.  Then
add its location to the prompt, for example: “The paper is at ~/Downloads/GordonVogelsang1993.pdf.”</p>
<p>The unpolarized coefficients should agree exactly.  For the polarized initial state, Gordon and Vogelsang use the
MS<sub>P</sub> scheme of their Eqs. (32)–(33), while FeynFacet uses HVBM with the finite term
$\\Delta f_{{qq}}(z)=-4C_F(1-z)$ that restores the helicity convention of the polarized MSbar parton distributions
(<a href="{src('docs/CONVENTIONS.md')}">conventions</a>); the two differ by a finite factorization term, which the
comparison should identify explicitly.</p>

<h2>What to expect</h2>
<ul>
<li>A new process folder with numbered steps, following the structure of
<a href="{src('Processes/LeptonQuark')}">Processes/LeptonQuark</a>: amplitudes, spin and color algebra, reduction and
master integrals, assembly, and comparison; each step as a <code>.wl</code> file with a generated <code>.nb</code>
notebook. The agent chooses where to create it, which may be outside the supplied repository; ask it for the full
path.</li>
<li>The finite coefficients in FeynFacet's form: a regular term, the distributions
$[1/(1-w)]_+$ and $[\\ln(1-w)/(1-w)]_+$, and $\\delta(1-w)$, with their dependence on the renormalization,
factorization and photon-fragmentation scales.</li>
<li>A summary of the conventions and of the checks: double and single poles cancel, the real emission is independent
of the gluon reference vectors, and the scale dependence follows the evolution equations.</li>
</ul>
<p>The computation itself takes minutes; an agent session that builds and validates the new workflow takes longer.
Agent output varies from run to run: read the notebooks as you would any calculation, and rely on the checks and the
comparison rather than on the agent's summary.</p>
{COPY_SCRIPT}
"""
    page("example.html", "Example", body)


REFERENCES = [
    ("F. Aversa, P. Chiappetta, M. Greco, J.-P. Guillet", "QCD corrections to parton-parton scattering processes",
     "Nucl. Phys. B 327 (1989) 105", "https://doi.org/10.1016/0550-3213(89)90288-5",
     "unpolarized proton–proton channel classes; its Fortran code is the reference of the UU comparison"),
    ("B. Jäger, A. Schäfer, M. Stratmann, W. Vogelsang", "Next-to-leading order QCD corrections to high-pT pion production in longitudinally polarized pp collisions",
     "Phys. Rev. D 67 (2003) 054005", "https://arxiv.org/abs/hep-ph/0211007", "polarized proton–proton classes; reference of the LL comparison"),
    ("A. Mukherjee, M. Stratmann, W. Vogelsang", "Next-to-leading order QCD corrections to single-inclusive hadron production in transversely polarized pp and p̄p collisions",
     "Phys. Rev. D 72 (2005) 034011", "https://arxiv.org/abs/hep-ph/0506315",
     "double-transverse proton–proton classes; its Fortran code is the reference of the TT comparison"),
    ("J. Soffer, M. Stratmann, W. Vogelsang", "Accessing transversity in double-spin asymmetries at the BNL-RHIC",
     "Phys. Rev. D 65 (2002) 114024", "https://arxiv.org/abs/hep-ph/0204058",
     "lowest-order double-transverse cross sections (Born table of the TT classes)"),
    ("P. Hinderer, M. Schlegel, W. Vogelsang", "Single-inclusive production of hadrons and jets in lepton-nucleon scattering at NLO",
     "Phys. Rev. D 92 (2015) 014001", "https://arxiv.org/abs/1505.06415", "unpolarized lepton–quark and off-diagonal coefficients"),
    ("P. Hinderer, M. Schlegel, W. Vogelsang", "Double-longitudinal spin asymmetry in single-inclusive lepton scattering at NLO",
     "Phys. Rev. D 96 (2017) 014002", "https://arxiv.org/abs/1703.10872", "polarized lepton–quark and off-diagonal coefficients"),
    ("M.-X. Luo, T.-Z. Yang, H. X. Zhu, Y. J. Zhu", "Unpolarized quark and gluon TMD PDFs and FFs at N3LO",
     "JHEP 06 (2021) 115", "https://arxiv.org/abs/2012.03256", "NLO TMD matching coefficients and soft function"),
    ("F. Yuan", "Azimuthal asymmetric distribution of hadrons inside a jet at hadron collider",
     "Phys. Rev. Lett. 100 (2008) 032003", "https://arxiv.org/abs/0709.3272",
     "lowest-order Collins hard factors; reference of the transfer Born table"),
    ("R. K. Ellis, W. J. Stirling, B. R. Webber", "QCD and Collider Physics",
     "Cambridge University Press (1996)", "https://doi.org/10.1017/CBO9780511628788",
     "Table 7.1: the four-dimensional spin-averaged Born table"),
    ("L. E. Gordon, W. Vogelsang", "Polarized and unpolarized prompt photon production beyond the leading order",
     "Phys. Rev. D 48 (1993) 3136", "https://doi.org/10.1103/PhysRevD.48.3136",
     "prompt-photon production at NLO; the benchmark of the example calculation"),
    ("A. P. Contogouris, Z. Merebashvili, O. Korakianitis", "Inclusive leptoproduction of a transversely polarized baryon on a transversely polarized proton",
     "Phys. Lett. B 365 (1996) 334", "https://doi.org/10.1016/0370-2693(95)01095-5", "an earlier transverse calculation for a different observable"),
    ("G. 't Hooft, M. Veltman", "Regularization and renormalization of gauge fields",
     "Nucl. Phys. B 44 (1972) 189", "https://doi.org/10.1016/0550-3213(72)90279-9",
     "dimensional regularization and the treatment of $\\gamma_5$ in $D$ dimensions: the HVBM scheme used for the helicity and transverse-spin results"),
    ("P. Breitenlohner, D. Maison", "Dimensional renormalization and the action principle",
     "Commun. Math. Phys. 52 (1977) 11", "https://doi.org/10.1007/BF01609069",
     "the consistent formulation of the HVBM scheme"),
]


def references_page():
    items = "".join(f'<li>{a}, <a href="{url}"><i>{t}</i></a>, {j}.<br><span class="note">{note}</span></li>'
                    for a, t, j, url, note in REFERENCES)
    tools = ('<a href="https://www.feynarts.de/">FeynArts</a>, <a href="https://feyncalc.github.io">FeynCalc</a>, '
             '<a href="https://form-dev.github.io/form/">FORM</a> and its '
             '<a href="https://www.nikhef.nl/~form/maindir/packages/color/color.html">color package</a>, '
             '<a href="https://library.wolfram.com/infocenter/MathSource/2987/">Tracer</a>, '
             '<a href="https://kira.hepforge.org/">Kira</a> with <a href="https://home.bway.net/lewis/">Fermat</a>, '
             '<a href="https://subtropi.ca/">SubTropica</a> with HyperIntica.')
    page("references.html", "References", f"<h1>References</h1><ol class=refs>{items}</ol><h2>Tools</h2><p>{tools}</p>")


def paper_page():
    if PAPER:
        body = f"""
<h1>Paper</h1>
<p class="lead"><i>FeynFacet: AI agents for analytic NLO QCD</i>, Z.-B. Kang, W. Vogelsang, C. Zhang.</p>
<p><a href="{PAPER}">arXiv:{html.escape(PAPER.rsplit("/", 1)[-1])}</a></p>"""
    else:
        body = """
<h1>Paper</h1>
<p class="lead"><i>FeynFacet: AI agents for analytic NLO QCD</i>, Z.-B. Kang, W. Vogelsang, C. Zhang.</p>
<p>In preparation; the link to the paper will appear here.</p>"""
    page("paper.html", "Paper", body)


def copy_assets():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    shutil.copyfile(ROOT / "LICENSE", OUT / "LICENSE.txt")
    figures = SITE / "figures"
    if figures.exists():
        shutil.copytree(figures, OUT / "figures", ignore=shutil.ignore_patterns("*.wls", "*.wl", "*.log"))
    (OUT / ".nojekyll").write_text("")
    if NOINDEX:
        (OUT / "robots.txt").write_text("User-agent: *\nDisallow: /\n")


if __name__ == "__main__":
    copy_assets()
    overview()
    results_overview()
    pp_page()
    tt_page()
    transfer_page()
    lepton_pages()
    method_page()
    reproduce_page()
    example_page()
    references_page()
    paper_page()
    print("written:", sorted(str(p.relative_to(OUT)) for p in OUT.rglob("*.html")))
