(* ::Title:: *)
(*Start here*)


(* ::Text:: *)
(*This notebook loads the package, checks that every program is installed, and runs the first process.  Evaluate the cells one after the other (Shift+Enter).  The calculations themselves are in the notebooks under Processes/; open them to read and evaluate each step.*)


(* ::Section:: *)
(*1. Load the package*)


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "FeynFacet", "FeynFacet.wl"}]]


(* ::Section:: *)
(*2. Check the installation*)


(* ::Text:: *)
(*Every row must show a location (and a version where the program has one).  A row marked NOT FOUND shows what to install; after installing, evaluate the cell again.  If a program is installed but not found, give its path in Paths.wl (see Paths.example.wl).*)


(* ::Input:: *)
CheckInstallation[]


(* ::Section:: *)
(*3. Run the lepton-quark process*)


(* ::Text:: *)
(*RunSteps evaluates the five step files of Processes/LeptonQuark in order, each in a fresh kernel: FeynArts amplitudes, the FeynCalc algebra, the Kira reduction and master integrals, and the unpolarized and polarized assemblies, which end by comparing with the published coefficients.  About five minutes.  Use "FORM" or "Tracer" as the second argument for the other algebra routes.*)


(* ::Input:: *)
RunSteps["LeptonQuark", "FeynCalc"]


(* ::Input:: *)
Get[FileNameJoin[{NotebookDirectory[], "Processes", "LeptonQuark", "Results", "FeynCalc", "Polarized.wl"}]]["Coefficients"] // TableForm


(* ::Text:: *)
(*The results are now in Processes/LeptonQuark/Results/FeynCalc/.  The step notebooks are in Processes/LeptonQuark/:*)


(* ::Input:: *)
If[$FrontEnd =!= Null, SystemOpen[FileNameJoin[{NotebookDirectory[], "Processes", "LeptonQuark"}]]]


(* ::Section:: *)
(*4. Run the TMD process*)


(* ::Text:: *)
(*The NLO TMD beam-function matching coefficients and the soft function (Processes/TMD), generated from a FeynArts model in which Wilson lines and operators are fields; about one minute.  Regenerated diagrams are written to Processes/TMD/Results/figures/.*)


(* ::Input:: *)
RunSteps["TMD", "FeynCalc"]


(* ::Section:: *)
(*5. The other processes*)


(* ::Text:: *)
(*The off-diagonal lepton channels l q -> g + X and l g -> q + X (Processes/LeptonOffDiagonal, about four minutes, compared with Hinderer, Schlegel and Vogelsang), the transverse spin transfer l q -> q + X for the Collins asymmetry with its resolved-photon Born table (Processes/Collins, about ten minutes, compared with the paper), the table of all 2 -> 2 Born channels (Processes/PartonicBorn, eight minutes) and the 16 proton-proton classes (Processes/Hadronic, one class at a time, minutes to hours; compared with the reference codes) run the same way:*)


(* ::Input:: *)
RunSteps["LeptonOffDiagonal", "FeynCalc"]


(* ::Input:: *)
RunSteps["Collins", "FeynCalc"]


(* ::Input:: *)
RunSteps["PartonicBorn", "FeynCalc"]


(* ::Text:: *)
(*A proton-proton class is selected with the class variable of its step files (the fourth argument of RunSteps overrides it); the Born table above must exist.  Class 2 (q q' -> g + X, no virtual correction) takes about ten minutes, class 1 about twenty:*)


(* ::Input:: *)
RunSteps["Hadronic", "FeynCalc", All, "class=2"]
