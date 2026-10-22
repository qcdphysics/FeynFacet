* square.frm: one product M_I M_J^* of the squaring method of step 2.
*
* Run in a process directory prepared by step 2 (see common.h) as
*    form -d SPIN=U -d I=3 -d J=5 -d OUT=name square.frm
* M_I is diagram I (amp<I>.inc), M_J^* the complex conjugate of diagram J
* (bar<J>.inc); for the one-loop interference amp<I>.inc is the one-loop
* graph I and bar1.inc the conjugate of the sum of the Born diagrams.  The
* result, summed over the spins and colours of all partons (with the
* densities of the polarized incoming partons for SPIN = P, and without the
* averages, which step 2 applies), is written to <OUT>.m as the Mathematica
* list {real part, imaginary part}.
*
* Each diagram file defines
*    Amp = sum_x c_x LorentzPartA(x)     (Bar = sum_y c_y^* LorentzPartB(y)
*                                         for the conjugate)
* with the distinct colour structures c_x of the diagram, and the Lorentz parts
*    LorA1, LorA2, ...                   (LorB1, ... ; `NAMP' and `NBAR' of them).
#-
Off statistics;
#include common.h
#include amp`I'.inc
#include bar`J'.inc
.sort

* 1. Transversality of the gluon polarizations, applied once to every Lorentz
*    part before any product: eps(l).l = 0 and eps(l).r = 0, with r the
*    reference vector of the later polarization sum (the procedure is written
*    into process.h by step 2).
#call transversality
.sort

* 2. Colour first: the colour factor of every pair of colour structures
*    c_x c_y^* is evaluated by color.h while the Lorentz parts are hidden
*    (only (number of structures)^2 small terms).  color.h works with the
*    default dimension NA of its adjoint indices, so the dimension is switched
*    to NA for this part.
#do x = 1, `NAMP'
  Hide LorA`x';
#enddo
#do y = 1, `NBAR'
  Hide LorB`y';
#enddo
L Square = Amp*Bar;
.sort
Dimension NA;
Drop Amp, Bar;
#call docolor
#call sunConstants
.sort

* 3. Back to D dimensions.  The Lorentz parts are put in, and the spinor
*    lines of the diagram and of the conjugate diagram are sewn into traces
*    with the spin sums and the densities of the polarized incoming quarks
*    (procedure sew of common.h).  Only Square is kept.
Dimension D;
#do x = 1, `NAMP'
  Unhide LorA`x';
#enddo
#do y = 1, `NBAR'
  Unhide LorB`y';
#enddo
if ( count(T,1,f,1,d33,1,d44,1) > 0 ) exit "square.frm: colour algebra not completed";
#do x = 1, `NAMP'
  id LorentzPartA(`x') = LorA`x';
#enddo
#do y = 1, `NBAR'
  id LorentzPartB(`y') = LorB`y';
#enddo
.sort
Drop;
NDrop Square;
#call sew
.sort

* 4. The traces in D dimensions (lines 1, 2: the lines of the external quarks;
*    3: the quark loop of a one-loop graph).
tracen,1;
tracen,2;
tracen,3;
.sort

* 5. The auxiliary vectors (propagator momenta written as one vector inside
*    the gamma strings) are expanded, and the polarization sums and densities
*    of the gluons are multiplied in (common.h).
#call auxiliary
#call polarizationSums
.sort

* 6. The epsilon tensors (helicity difference only) are contracted in four
*    dimensions (procedure hvbm of common.h), with the Minkowski sign.
#call hvbm
Multiply `MINKOWSKISIGN';
.sort

* 7. Kinematics: the scalar products in terms of s, t, u (and A, B, C, spe3
*    for the real emission, L2, Lpa, Lpb, Lp, Lhat2 for the loop momentum).
#if "`KINEMATICS'" == "Real"
  #call realKinematics
#else
  #call bornKinematics
  #if "`KINEMATICS'" == "Virtual"
    #call loopProducts
  #endif
#endif
.sort
if ( count(pa,1,pb,1,p,1,k1,1,k2,1,loop,1,k1b,1,k2b,1,loopb,1) > 0 ) exit "square.frm: a scalar product was left";

* 8. Real and imaginary parts, written in Mathematica syntax.
Bracket i_;
.sort
L ReSquare = Square[1];
L ImSquare = Square[i_];
.sort
Drop Square;
Format Mathematica;
Format 250;
#write <`OUT'.m> "{(%E),\n(%E)}", ReSquare, ImSquare
.end
