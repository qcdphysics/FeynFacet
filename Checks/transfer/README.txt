Independent transverse-spin-transfer assembly check

Requirements: a Wolfram kernel, normally invoked using wolframscript. Allow at
least 8 GiB of available memory per kernel. The scientific programs launch no
parallel workers and do not require Python, Linux /proc, FeynCalc, FORM, TRACER
or the FeynFacet package at execution time.

Inputs. First run the ordinary process steps 1--4 to generate Algebra.wl,
Reduced.wl and Transfer.wl in Processes/HadronicTransfer/Results/<route>/class-n.
Algebra and Reduced provide shared mathematical inputs; Transfer is only the
comparison target. Use a fresh output directory outside the source repository.
From the directory containing the two checkers, one class can be run as:

 INDEPENDENT_SOURCE=/path/to/FeynFacet \
 INDEPENDENT_ROUTE=FeynCalc \
 INDEPENDENT_CT_ONLY=0 \
 INDEPENDENT_CLASS=1 \
 INDEPENDENT_OUTPUT=/path/to/new-direct-checks \
 wolframscript -file IndependentAssembly.wls

The class-1 subdirectory receives independent.wl and summary.json. An unresolved
identity exits nonzero; there is no numerical fallback. If your installation
uses a command named math, replace the final line with
 math -noprompt -script IndependentAssembly.wls
while retaining the same environment variables.

To check all seven direct classes sequentially in a POSIX shell:

 feynfacet_root=/path/to/FeynFacet
 direct_checks=/path/to/new-direct-checks
 mkdir "$direct_checks" || exit 1
 for partonic_class in 1 3 6 10 11 13 17; do
   INDEPENDENT_SOURCE="$feynfacet_root" \
   INDEPENDENT_ROUTE=FeynCalc \
   INDEPENDENT_CT_ONLY=0 \
   INDEPENDENT_CLASS="$partonic_class" \
   INDEPENDENT_OUTPUT="$direct_checks" \
   wolframscript -file IndependentAssembly.wls || exit 1
 done

The same checker supports FORM and Tracer result routes. Set
INDEPENDENT_CT_ONLY=1 only for a short counterterm-only check; this omits the
real/virtual reconstruction and is explicitly identified in its summary.

Crossing. After the ordinary process step 5 has generated the crossed class-nx
Algebra/Reduced inputs, independently assemble all seven crossed classes:

 crossed_checks=/path/to/new-crossed-checks
 mkdir "$crossed_checks" || exit 1
 for partonic_class in 1x 3x 6x 10x 11x 13x 17x; do
   INDEPENDENT_SOURCE="$feynfacet_root" \
   INDEPENDENT_ROUTE=FeynCalc \
   INDEPENDENT_CT_ONLY=0 \
   INDEPENDENT_CLASS="$partonic_class" \
   INDEPENDENT_OUTPUT="$crossed_checks" \
   wolframscript -file IndependentAssembly.wls || exit 1
 done

It uses the physical species of class n with the crossed D-dimensional Born
tensors and scalar integrals. Compare the two newly constructed results using:

 INDEPENDENT_DIRECT_RESULTS="$direct_checks" \
 INDEPENDENT_CROSSED_RESULTS="$crossed_checks" \
 INDEPENDENT_CROSSING_OUTPUT=/path/to/new-crossing-comparison \
 wolframscript -file CompareIndependentCrossing.wls

Alternatively use math -noprompt -script CompareIndependentCrossing.wls with
the same environment variables. The comparison reads independent.wl on both
sides and compares all 56 finite-component slots exactly. The reduced real
integrals, Born tensors and master values already agree between these inputs;
the virtual expressions can differ in form. This checks assembly consistency,
while the nontrivial crossing tests are in production step 5. The comparison
does not consume the original crossing-assembly result.

Coverage. Exact checks retain symbolic kinematics, color and flavour factors,
and separate initial/final factorization and renormalization scales. They test
both XX and NN, every regular/D0/D1/delta coefficient, every finite counterterm,
the real and virtual finite pieces, double/single-pole cancellation, and equality
of the independently assembled XX and NN finite tensors. The SU(N) relation
CF=(CA^2-1)/(2CA) and the physical region s>0, 0<v,w<1 are explicit. Each
summary binds its mathematical inputs and the checker by SHA256 hashes;
independent.wl contains the reconstructed expressions and individual verdicts.

Comparison controls add polynomial perturbations to all 56 coefficient
differences; they test rejection by the comparison routine. Four class-1
physics-term omission tests instead rebuild the actual counterterms without
the evanescent Born terms or the fragmentation epsilon factor, separately for
XX and NN, and require the resulting finite mismatch to be detected.
Logarithm products are expanded only after real nonzero factor signs and a
positive complete product have been certified. Negative, complex and
unknown-sign arguments retain their original logarithms and phases.

Independence. This is a separate distribution-expansion and finite-assembly/
subtraction implementation from shared scalar-integral inputs. It does not
independently perform amplitude generation, spin algebra, IBP reduction or
master integration. DERIVATION.txt explains its Jacobians, dimensional factors
and endpoint-subtraction identity.
