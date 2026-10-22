(* Only needed if CheckInstallation[] reports a program as NOT FOUND although
   it is installed.  Copy this file to Paths.wl (same directory), uncomment
   the line for that program and give its full path.  Paths.wl is read
   automatically when the package loads; it is private to your computer. *)

(*
$ToolPaths["Fermat"] = "/absolute/path/to/Fermat/Ferl7/fer64";
$ToolPaths["Kira"] = "/absolute/path/to/kira";
$ToolPaths["FORM"] = "/absolute/path/to/form";
$ToolPaths["TFORM"] = "/absolute/path/to/tform";
$ToolPaths["Polymake"] = "/absolute/path/to/polymake";
$ToolPaths["Tracer"] = "/absolute/path/to/Tracer.m";
$ToolPaths["FORMColor"] = "/absolute/path/to/color.h";
*)

(* A Wolfram package (SubTropica, FeynCalc) installed outside the standard
   locations: load it from here, e.g.
PacletDirectoryLoad["/absolute/path/to/SubTropica"];
AppendTo[$Path, "/absolute/path/to/the/folder/containing/FeynCalc"];
*)
