# Memory and parallel runs

Memory use depends on the channel, algebra route and calculation stage.
Representative process-tree resident-memory measurements are:

| Calculation | Stage | Approximate sampled peak |
|---|---|---|
| Incoming-TT class 12, q↑ q̄↑ → g + X | Finite assembly (step 4) | 10.9 GiB |
| Transverse spin transfer, class 13, q↑ g → q↑ + X | Crossing calculation (step 5, including finite assembly) | 13.4–14.4 GiB |

These periodic measurements include the Wolfram process and its subprocesses.
They are lower bounds on the true stage peaks, not sufficient-memory
requirements. Leave memory for the operating system and other programs.
Start with one heavy scalar stage at a time and measure the chosen stage
before increasing concurrency.

The `workers` setting controls step 2 algebra concurrency. Its default `auto`
uses at most eight workers. After generating a class's amplitudes, use
`workers=1` to process its algebra jobs sequentially:

```sh
FEYNFACET_SET="route=FeynCalc;class=12;workers=1" \
  wolframscript -file Checks/run_step.wls \
  Processes/HadronicTT/2-Algebra-FeynCalc.wl
```

Steps 3 and 4 use one Wolfram kernel for scalar reduction and finite assembly;
Kira runs as an external process. These steps do not read `workers`, so
reducing algebra concurrency does not limit their memory use. Finite assembly
can require more memory than scalar reduction. Run large classes and crossing
pairs separately until their resource requirements are known.

The separate [independent assembly checker](../Checks/transfer/README.txt)
reads existing `Reduced.wl` results. Its 8 GiB-per-kernel guidance applies to
that checker. Loading the distributed [coefficient tables](../Coefficients/README.md)
does not require rerunning any algebra or scalar reduction.
