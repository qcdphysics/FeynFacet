(* Usage:
   wolframscript -file Checks/check_lepton_scales.wls ROOT OUTPUT [ROUTE]
   ROOT contains FeynFacet/, Coefficients/ and Processes/.
   OUTPUT is a new local directory for this check's report and residue table.
   ROUTE defaults to FeynCalc; use All to check all three algebra routes.
   The process notebooks must already have produced their result files.
   This check performs exact symbolic identities, without phenomenology. *)
$HistoryLength=0;
Get[FileNameJoin[{DirectoryName[$InputFileName], "ScriptArguments.wl"}]];
args=ScriptArguments[];
If[Length[args]<2,Print["Give repository ROOT and a new OUTPUT directory."];Quit[2]];
root=ExpandFileName[args[[1]]]; output=ExpandFileName[args[[2]]];
routeChoice=If[Length[args]>=3,args[[3]],"FeynCalc"];
routes=If[routeChoice==="All",{"FeynCalc","FORM","Tracer"},{routeChoice}];
If[!And@@(MemberQ[{"FeynCalc","FORM","Tracer"},#]&/@routes),Quit[2]];
If[DirectoryQ[output],Print["OUTPUT already exists; choose a new directory."];Quit[2]];
CreateDirectory[output,CreateIntermediateDirectories->True];
module=FileNameJoin[{root,"FeynFacet","LeptonScaleDependence.wl"}];
If[!FileExistsQ[module],Print["Missing ",module];Quit[2]];
Get[module];
require[condition_,message_]:=If[!TrueQ[condition],Print[message];Quit[1]];
sha[file_]:=IntegerString[FileHash[file,"SHA256"],16,64];
inputs=<|module->sha[module]|>;
read[file_]:=(require[FileExistsQ[file],"Missing input: "<>file];inputs[file]=sha[file];Get[file]);
region=s>0 && muR>0 && muPDF>0 && muFF>0 && muPhoton>0 && 0<v<1 && 0<w<1;
zeroQ[x_]:=And@@(TrueQ[FullSimplify[Together[#],region]===0]&/@Flatten[{x}]);
couplings={ee->Sqrt[4 Pi alphaEM],gs->Sqrt[4 Pi alphaS]};
threshold[{r_,p_,d_}]:={r,p,0,d};
rows={}; residueTable=<||>;
Do[
 table=read[FileNameJoin[{root,"Coefficients",process<>".wl"}]];
 Do[
  directory=FileNameJoin[{root,"Processes",process,"Results",route}];
  algebra=read[FileNameJoin[{directory,"Algebra.wl"}]];
  Do[
   result=read[FileNameJoin[{directory,spin<>".wl"}]];
   channels=If[process==="LeptonQuark",{"QQ"},{"QG","GQ"}];
   bornD=Normal[Series[algebra["Born"][spin]/(16 Pi^2 s v)/.D->4-2 eps/.
     couplings/.{t->-s(1-v),u->-s v},{eps,0,1}]];
   y=(1-v)/(1-v w);
   Do[
    photonSquare=If[channel==="QQ",algebra["Photon"][spin],algebra["Photon"][channel][spin]];
    photonD=Normal[Series[photonSquare/(16 Pi^2 y s(1-v w))/.D->4-2 eps/.
      couplings/.{sh->y s,th->-s(1-v),uh->-y s v w},{eps,0,1}]];
    residues=FeynFacet`LeptonScales`LeptonScaleResidues[
      bornD,photonD,spin,channel,{v,w,eps},{alphaS,alphaEM,CF,TR}];
    totalResidue=Total[Values[residues]];
    finite=If[channel==="QQ",result["Finite"],result["Finite"][channel]];
    published=If[channel==="QQ",table["Results"][spin],table["Results"][spin][channel]];
    pref=published["Prefactor"];
    fullTable=MapThread[Times,{published["Coefficients"],{pref,pref/.w->1,pref/.w->1,pref/.w->1}}];
    savedResidue=If[channel==="QQ",-Map[Coefficient[#,EulerGamma]&,result["Counterterms"]],
      threshold[-Map[Coefficient[#,EulerGamma]&,result["Counterterms"]["parton"][channel]+result["Counterterms"]["photon"][channel]]]];
    checks=<|"ExistingHardScaleDerivative"->zeroQ[(muR/2)D[finite,muR]+totalResidue],
      "PublishedTableScaleDerivative"->zeroQ[Map[Coefficient[#,LR]&,fullTable]+totalResidue],
      "SavedCollinearResidue"->zeroQ[savedResidue-totalResidue]|>;
    If[channel=!="QQ",
      checks["SeparatePhotonResidue"]=zeroQ[threshold[-Map[Coefficient[#,EulerGamma]&,result["Counterterms"]["photon"][channel]]]-residues["Photon"]];
      leg=If[channel==="QG","FF","PDF"];
      checks["SeparatePartonResidue"]=zeroQ[threshold[-Map[Coefficient[#,EulerGamma]&,result["Counterterms"]["parton"][channel]]]-residues[leg]]];
    AppendTo[rows,<|"Process"->process,"Route"->route,"Spin"->spin,"Channel"->channel,"Checks"->checks|>];
    require[And@@(TrueQ/@Values[checks]),"Scale identity failed: "<>process<>" "<>spin<>" "<>channel<>" "<>route];
    residueTable[spin<>" "<>channel]=Map[FullSimplify[#,region]&,residues],
    {channel,channels}],{spin,{"Unpolarized","Polarized"}}],{route,routes}],
 {process,{"LeptonQuark","LeptonOffDiagonal"}}];
require[Length[rows]===6 Length[routes],"Incomplete channel coverage"];
require[And@@KeyValueMap[FileExistsQ[#1]&&sha[#1]===#2&,inputs],"An input changed during the check"];
Put[residueTable,FileNameJoin[{output,"LeptonScaleResidues.wl"}]];
Export[FileNameJoin[{output,"scale-check.json"}],
 <|"Passed"->True,"Rows"->rows,"InputSHA256"->inputs,"WolframVersion"->$Version|>,"RawJSON"];
Print[Length[rows]," exact channel checks passed."];
Quit[0];
