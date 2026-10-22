(* ::Package:: *)

(* Tracer.wl: load the classic Tracer package (M. Jamin and M. E. Lautenbacher,
   Comput. Phys. Commun. 74 (1993) 265) in a kernel of its own, and build Tracer
   gamma strings from amplitude records.  Tracer's own commands are used
   directly in the notebooks:
     G[line, p, {mu}, q, ...]   the ordered string  p-slash gamma^mu q-slash ...
                                on fermion line "line"; G5 inside it is gamma_5
     Spur[line]                 take the trace on that line when G is evaluated
     VectorDimension[D]         work in D dimensions
     AntiCommute[off]           HVBM gamma_5 (needed for helicity differences)
     S[p, q]                    p.q, S[p, {mu}] = p^mu, S[{mu},{nu}] = g^{mu nu}
     T[p], H[p]                 four-dimensional and (D-4)-dimensional parts
     OnShell[on, {p,0}, {p,q,s/2}, ...]   install masses and scalar products *)

BeginPackage["FeynFacet`"];

LoadTracer::usage = "LoadTracer[] loads Tracer.m.  Tracer redefines Dot and \
NonCommutativeMultiply on startup; these are restored so that ordinary \
Mathematica expressions keep working.  Tracer must not share a kernel with \
FeynCalc: run the amplitude-generation step first and this step in a new kernel.";

TracerString::usage = "TracerString[Chain[Ubar[p], ..., Uspinor[k]]] returns the \
list of Tracer arguments ({mu} for GammaMatrix[mu], p for Slash[p], G5) of the \
gamma matrices between the spinors, ready to be inserted into Tracer`G.  The \
spinors themselves are not included: the notebook supplies the density matrices.";

TracerStringReversed::usage = "TracerStringReversed[chain] is the reversed gamma \
string of the chain, i.e. the Dirac adjoint of the matrix product.";

Begin["`Private`"];

LoadTracer[] := Module[{path, contextPath = $ContextPath,
    dot = DownValues[Dot], ncm = DownValues[NonCommutativeMultiply]},
  Require[!MemberQ[$Packages, "FeynCalc`"],
    "Load Tracer in a fresh kernel, not together with FeynCalc."];
  If[MemberQ[$Packages, "Tracer`"], Return["Tracer already loaded"]];
  path = ToolPath["Tracer"];
  Block[{$ContextPath = {"System`"}},
    Quiet[Get[path], {Remove::rmnsm, Syntax::com, General::shdw}]];
  Unprotect[Dot, NonCommutativeMultiply];
  DownValues[Dot] = dot; DownValues[NonCommutativeMultiply] = ncm;
  Protect[Dot, NonCommutativeMultiply];
  $ContextPath = DeleteDuplicates[Append[contextPath, "Tracer`"]];
  "Tracer loaded from " <> path];

item[Global`GammaMatrix[a_]] := {a};
item[Global`Slash[p_]] := p;
item[Global`GammaFive] := Tracer`G5;

TracerString[chain_Global`Chain] := item /@ (List @@ chain)[[2 ;; -2]];
TracerStringReversed[chain_Global`Chain] := Reverse[TracerString[chain]];

End[];
EndPackage[];
