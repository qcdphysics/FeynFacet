# Small exact Kira regressions

These fixtures test the reduction tool, independently of any process workflow.
`expectations.json` is shared by the Wolfram toolchain gate and the optional
standalone SymPy verifier. Nothing reads or updates process results.

## Independently derived coefficients

Use the loop measure `d^d k/(i pi^(d/2))` and propagators `k^2+i0`.
Define the dimensionless massless bubble factor

```
G(a,b) = Gamma(a+b-d/2) Gamma(d/2-a) Gamma(d/2-b)
         / (Gamma(a) Gamma(b) Gamma(d-a-b)).
B(a,b;q2) = (-1)^(a+b) (-q2-i0)^(d/2-a-b) G(a,b).
```

This follows by combining the denominators with one Feynman parameter,
shifting the loop momentum and performing the beta integral. It supplies
the expectations without using another Kira run as the reference.

`python3 Checks/kira/derive_expectations.py` verifies the Gamma recurrences
against the shared expected coefficients (optional; requires SymPy).

### Triangle (LeptonQuark configuration, independent scratch directory)

Denominators are `L^2`, `(L+p1)^2`, `(L+p2)^2`, with `p1^2=p2^2=0` and
`(p1-p2)^2=t`. Feynman parameters give

```
T(1,1,1) = -(-t-i0)^(d/2-3)
           Gamma(3-d/2) Gamma(d/2-2)^2 / Gamma(d-3).
T(0,1,1) = B(1,1;t).
T(1,1,1)/T(0,1,1) = -2(d-3)/((d-4)t).
```

The ratios use only `Gamma(z+1)=z Gamma(z)`. The other two targets,
`T(1,1,0)` and `T(1,0,1)`, have only a lightlike external momentum and
vanish as scaleless integrals. Kira sets `t=1` in this fixture.

### Two-loop massless sunset

The first three denominators are `k^2`, `l^2`, `(p-k-l)^2`, with `p^2=q2`.
The last two complete the scalar-product basis and have index zero for
all targets. Integrating the `l` bubble and then the `k` bubble gives

```
J(a,b,c) = (-1)^(a+b+c) (-q2-i0)^(d-a-b-c)
           G(b,c) G(a,b+c-d/2).
```

Relative to `J(1,1,1)`, repeated use of the Gamma recurrence gives:

| Target | Coefficient at q2=1 |
|---|---|
| J(2,1,1) | `-(d-3)(3d-8)/(d-4)` |
| J(3,1,1) | `(d-3)(3d-10)(3d-8)/(2(d-6))` |
| J(2,2,1) | `(d-3)(3d-10)(3d-8)/(d-4)` |
| J(2,2,2) | `-3(d-5)(d-3)(3d-10)(3d-8)/(d-4)` |
| J(1,1,0) | `0` (factorized scaleless tadpoles) |

Restore `(q2)^(3-a-b-c)` for dimensionful coefficients. These are rational
identities in generic `d`, not evaluations at singular integer dimensions.

## Running

Normally run `python3 feynfacet.py check`; its Wolfram implementation needs
no additional Python package. It tests the Kira selected by `ToolPath`.

For release/CI without Wolfram, install SymPy in a separate environment and run:

```sh
python3 -m venv /tmp/feynfacet-kira-check-env
/tmp/feynfacet-kira-check-env/bin/pip install sympy==1.14.0
/tmp/feynfacet-kira-check-env/bin/python Checks/check_kira.py \
  --kira /absolute/path/to/kira --fermat /absolute/path/to/fer64 \
  --workdir /tmp/feynfacet-kira-check-new
```

The work directory must not exist. Each reduction has a 60-second timeout.
Unknown masters, missing rules, wrong coefficients and false scaleless zeros
fail. `report.json` and Kira logs remain in that directory.

The separate `Checks/kira_mod_mul.cpp` harness links against the actual
build's `libpyred.a`. It compares 100,072 edge/random cases against an
independent repeated-addition implementation (no floating point or 128-bit
division in the reference). The installer runs it once before installation.
This source/build diagnostic complements the executable reduction tests;
compiling a fixed standalone function alone would not validate a selected
installed Kira executable.
