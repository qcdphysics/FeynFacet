* Born.frm: |M|^2 of every massless 2 -> 2 partonic channel in D dimensions,
* summed over the colours of all four partons and over the spins of the
* outgoing ones; for the incoming partons either summed over their spins
* (SPIN=U) or weighted with their helicity densities (SPIN=P, the helicity
* difference).  The notebook 2-Algebra-FORM.wl runs it as
*   tform -w3 -d SPIN=U -d CHANNELS=29 Born.frm      (and again with SPIN=P)
* in the directory that holds amplitudes.inc and declarations.inc, both
* written by the notebook.  The squares are written, in channel order, to
* BornU.m or BornP.m as one Mathematica list.  The spin and colour averages of
* the incoming partons are applied in the notebook.
*
* For each channel i:
*   1. Square = Amp_i * conjugate(Amp_i)        (procedure conjugate, common.h)
*   2. sums over all colour indices with FORM's colour library color.h
*   3. the fermion lines closed into Dirac traces      (procedure sew)
*   4. spin densities, gluon polarization sums         (procedure densities)
*   5. propagators, D-dimensional traces, the epsilon tensors of the HVBM
*      scheme kept four-dimensional, contraction of all indices, invariants.
#-
Off statistics;
#include common.h
#include amplitudes.inc
.sort
* the squares are written in Mathematica syntax
Format Mathematica;
* The 29 channel amplitudes Amp1 ... Amp29 are hidden: they stay in storage
* untouched, and each is copied into the active expressions when its turn comes.
Hide;

*--------------------------------------------------------------------------
* Normalization of the HVBM gamma5 (SPIN=P only).
* common.h writes gamma5 = -(i/24) e_(h1,h2,h3,h4) g_(h1) g_(h2) g_(h3) g_(h4)
* with four-dimensional indices.  With this normalization the trace
* Tr[gamma^a gamma^b gamma^c gamma^d gamma5] must equal -4 i e_(a,b,c,d), the
* convention of the FeynCalc route ($LeviCivitaSign = -1).  G5Check is the
* difference of the two; it is printed to the log, and the notebook requires
* it to be 0.
*--------------------------------------------------------------------------
#if "`SPIN'" == "P"
L G5Check = g_(1,h1,h2,h3,h4)*(-i_/24)*e_(e1,e2,e3,e4)*g_(1,e1,e2,e3,e4) + 4*i_*e_(h1,h2,h3,h4);
tracen,1;
contract;
.sort
Print +s G5Check;
.sort
Drop G5Check;
#endif

* opening bracket of the output list
#write <Born`SPIN'.m> "HoldComplete[{"
#do i = 1, `CHANNELS'

*--------------------------------------------------------------------------
* 1. The square of the amplitude of channel i: a copy Amp, a copy Bar that
*    is complex conjugated (with primed summed indices), and their product.
*--------------------------------------------------------------------------
  L Amp = Amp`i';
  L Bar = Amp`i';
  .sort
  #call conjugate(Bar)
  L Square`i' = Amp*Bar;
  .sort

*--------------------------------------------------------------------------
* 2. Colour.  docolor (color.h) evaluates every colour contraction: the
*    generators T(i,j,a) = (T^a)_ij, the structure constants f(a,b,c) and the
*    deltas d_ of the amplitude times its conjugate, summed over all colour
*    indices, reduce to group invariants.  This happens before any Dirac
*    algebra; the colour factor of each term is a number.
*    color.h works with its own index declarations of dimension NA, so the
*    default dimension is set to NA while it runs and back to D afterwards.
*    The invariants of a general group are then specialized to SU(Nc) with
*    the normalization Tr(T^a T^b) = delta^ab/2 of the FeynCalc route:
*      I2R = 1/2 (the trace normalization), NR = Nc, NA = Nc^2 - 1,
*      cA = Nc (adjoint Casimir), cR = CF = (Nc^2-1)/(2 Nc) (quark Casimir),
*      and color.h's combination [cR-cA/2] = CF - CA/2 = -1/(2 Nc).
*    A generator, structure constant or higher invariant (d33) left over
*    would mean that the colour sum is incomplete; FORM then stops.
*--------------------------------------------------------------------------
  Dimension NA;
  Drop Amp,Bar;
  #call docolor
  .sort
  Dimension D;
  id I2R = 1/2;
  id NA = Nc^2 - 1;
  id NR = Nc;
  id cA = Nc;
  id cR = Nc/2 - 1/2/Nc;
  id [cR-cA/2] = -1/2/Nc;
  if ( count(T,1,f,1,d33,1) > 0 ) exit "colour algebra incomplete";

*--------------------------------------------------------------------------
* 3.-5. Dirac and Lorentz algebra.
*    sew          closes the fermion chains into one or two traces;
*    densities    inserts the quark spin sums or helicity densities and the
*                 gluon polarization sums or circular densities;
*    propagators  writes the propagator denominators as 1/s, 1/t, 1/u;
*    tracen,1/2   the D-dimensional traces of lines 1 and 2 (tracen is the
*                 trace for a general dimension, without 4-dimensional
*                 shortcuts);
*    ProjectEpsilon (declarations.inc) moves every Lorentz index that sits
*                 inside an epsilon tensor to a four-dimensional twin index:
*                 in the HVBM scheme the epsilon tensor is four-dimensional,
*                 so an index summed through it runs over four dimensions;
*    contract     contracts all repeated indices (epsilon tensors included);
*    scalarproducts  replaces the scalar products by s, t, u.
*    Finally the sign MINKOWSKISIGN (common.h) converts FORM's contraction
*    of two epsilon tensors to the Minkowski-space result.
*--------------------------------------------------------------------------
  #call sew
  #call densities
  #call propagators
  tracen,1;
  tracen,2;
  #call ProjectEpsilon
  #call scalarproducts
  .sort
  contract;
  #call scalarproducts
  Multiply `MINKOWSKISIGN';
  .sort

* append the square of channel i to the output list (comma-separated)
  #if `i' < `CHANNELS'
    #write <Born`SPIN'.m> "%E,", Square`i'
  #else
    #write <Born`SPIN'.m> "%E", Square`i'
  #endif
  Drop Square`i';
  .sort
#enddo
#write <Born`SPIN'.m> "}]"
.end
