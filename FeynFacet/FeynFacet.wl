(* ::Package:: *)

(* FeynFacet: shared, physics-free mechanics for the NLO notebooks.

   Every module exports a handful of documented functions.  Nothing here
   contains a hard-scattering coefficient: the physics (projectors, traces,
   integral families, master integrals, counterterms) is written out in the
   process notebooks under Processes/.  Load with

       Needs["FeynFacet`"]

   after adding the repository root to $Path, or via Get on this file. *)

BeginPackage["FeynFacet`"];

$FeynFacetDirectory::usage = "$FeynFacetDirectory is the repository root.";

EndPackage[];

$FeynFacetDirectory = DirectoryName[ExpandFileName[$InputFileName], 2];

Scan[Get[FileNameJoin[{$FeynFacetDirectory, "FeynFacet", #}]] &,
  {"Tools.wl", "Preprocessing.wl", "AlgebraReuse.wl", "Amplitudes.wl", "Kira.wl", "Masters.wl", "Distributions.wl",
   "Factorization.wl", "Compare.wl", "FORM.wl", "Tracer.wl"}];
