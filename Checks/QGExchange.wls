(* Reconstruct q g -> q g g from its 16-amplitude interference cache.
   No production squaring, real-emission sum or integration helper is loaded.
   Exchange comparisons are sampled at saved exact-rational D=4 points. *)
$HistoryLength=0;
sourceRoot=Environment["FEYNFACET_SOURCE"];
outputRoot=Environment["FEYNFACET_BENCHMARK_OUTPUT"];
If[!And@@(StringQ/@{sourceRoot,outputRoot}),Print["Set FEYNFACET_SOURCE and FEYNFACET_BENCHMARK_OUTPUT."];Quit[2]];
If[!DirectoryQ[outputRoot],CreateDirectory[outputRoot,CreateIntermediateDirectories->True]];
cacheRoot=FileNameJoin[{sourceRoot,"Processes","Hadronic","Results","FeynCalc","class-13","pairs"}];
checks={};failures={};results=<||>;inputs=<||>;
record[name_,test_]:=(AppendTo[checks,<|"name"->name,"passed"->TrueQ[test]|>];
 If[!TrueQ[test],AppendTo[failures,name];Print["FAIL ",name]];
 Export[FileNameJoin[{outputRoot,"progress.json"}],<|"checks"->checks,"failures"->failures|>,"RawJSON"];TrueQ[test]);
finiteRational[value_]:=MatchQ[value,_Integer|_Rational];
zeroQ[value_]:=TrueQ[Cancel[Together[value]]===0];
expectedPairs=Flatten[Table[{i,j},{i,16},{j,i,16}],1];
points={
 {s->17,t->-7,u->-6,A->-2,B->-5,C->3,Nc->3,gs->1},
 {s->23,t->-8,u->-9,A->-4,B->-6,C->4,Nc->5,gs->1},
 {s->31,t->-12,u->-11,A->-5,B->-9,C->6,Nc->7,gs->1}};
wrongPoint={s->29,t->-11,u->-10,A->-3,B->-9,C->4,Nc->3,gs->1};
(* k1+k2=pa+pb-p implies this nonthreshold exchange, with s4=s+t+u. *)
exchangeMap={A->-(s+t)-A,B->-(s+u)-B,C->-t-u-C};
record["map/constraint",zeroQ[(A+B+C+s+t+u)/.exchangeMap/.C->-(s+t+u)-A-B]];
record["map/involution",And@@(zeroQ[#]&/@((({A,B,C}/.exchangeMap)/.exchangeMap)-{A,B,C}))];
Do[
 files=FileNames["real-q_g_q_g_g-"<>spin<>"-*.wl",cacheRoot];
 pairs=(ToExpression/@Take[StringSplit[FileBaseName[#],"-"],-2])&/@files;
 If[!record[spin<>"/136-upper-triangle-products",Sort[pairs]===Sort[expectedPairs]],Quit[1]];
 hashes=AssociationThread[files,IntegerString[FileHash[#,"SHA256"],16,64]&/@files];
 inputs=Join[inputs,hashes];
 data=Get/@files;
 If[!record[spin<>"/cache-names",And@@MapThread[#1["Name"]===FileBaseName[#2]&,{data,files}]],Quit[1]];
 values=Lookup[data,"Value"];
 (* Count literal exact zeros in the saved generic expressions, without
    inferring whether their origin is colour orthogonality or Dirac algebra. *)
 zeroPositions=Flatten[Position[values,0,{1}]];
 record[spin<>"/12-exact-zero-products",Length[zeroPositions]===12];
 coefficients=If[#[[1]]===#[[2]],1,2]&/@pairs;
 sampleValues={};probeProducts={};
 Do[
  point=points[[probe]];
  swapped=Join[Thread[{A,B,C}->(({A,B,C}/.exchangeMap)/.point)],point];
  leftValues=Cancel[Together[#/.{D->4,spe3->0}/.point]]&/@values;
  rightValues=Cancel[Together[#/.{D->4,spe3->0}/.swapped]]&/@values;
  good=And@@(finiteRational/@Join[leftValues,rightValues]);
  record[spin<>"/probe-"<>ToString[probe]<>"/finite",
   good && ((A+B+C+s+t+u)/.point)===0 && ((A+B+C+s+t+u)/.swapped)===0];
  (* At real rational kinematics Re removes any imaginary cached interference;
     diagonal products carry weight1, and off-diagonal pairs carry weight2. *)
  left=Total[coefficients Re[leftValues]];right=Total[coefficients Re[rightValues]];
  record[spin<>"/probe-"<>ToString[probe]<>"/exchange",good && left===right];
  AppendTo[probeProducts,leftValues];
  AppendTo[sampleValues,<|"Point"->point,"ExchangedPoint"->swapped,"Left"->left,"Right"->right|>],
  {probe,Length[points]}];
 (* A nonzero value at any saved point proves an expression is not identically
    zero. This supplements the generic literal-zero census. *)
 nonzero=Count[Transpose[probeProducts],row_/;AnyTrue[row,finiteRational[#] && #=!=0&]];
 record[spin<>"/124-proven-nonzero-products",nonzero===124];
 wrongMap={A->-(s+t)-A,B->-(s+u)-B,C->s-C};
 wrongSwapped=Join[Thread[{A,B,C}->(({A,B,C}/.wrongMap)/.wrongPoint)],wrongPoint];
 leftValues=Cancel[Together[#/.{D->4,spe3->0}/.wrongPoint]]&/@values;
 wrongValues=Cancel[Together[#/.{D->4,spe3->0}/.wrongSwapped]]&/@values;
 left=Total[coefficients Re[leftValues]];wrong=Total[coefficients Re[wrongValues]];
 record[spin<>"/control/wrong-exchange-map",And@@(finiteRational/@Join[leftValues,wrongValues]) && left=!=wrong];
 results[spin]=<|"Products"->Length[values],"ExactZeroProducts"->Length[zeroPositions],
  "ProvenNonzeroProducts"->nonzero,"Samples"->sampleValues,"WrongMapControl"->{left,wrong}|>,
 {spin,{"UU","LL"}}];
record["272-inputs",Length[inputs]===272];
record["input-integrity",And@@KeyValueMap[IntegerString[FileHash[#1,"SHA256"],16,64]===#2&,inputs]];
Put[<|"Results"->results,"Inputs"->inputs,"Checks"->checks|>,FileNameJoin[{outputRoot,"exchange.wl"}]];
Export[FileNameJoin[{outputRoot,"summary.json"}],<|"passed"->(failures==={}),"checkCount"->Length[checks],"checks"->checks,"failures"->failures,"inputs"->inputs,
 "dimension"->4,"comparison"->"Three exact-rational points, Nc=3,5,7; not a generic exchange proof",
 "checkerSHA256"->IntegerString[FileHash[$InputFileName,"SHA256"],16,64]|>,"RawJSON"];
Print["QG EXCHANGE: ",Length[checks]," checks, ",Length[failures]," failures"];
Quit[If[failures==={},0,1]];
