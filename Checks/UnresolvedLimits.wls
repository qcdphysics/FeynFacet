(* Exact unresolved-limit checks of selected hadronic algebra files.
   Scope: qq' -> qq'g, class1, UU/LL and transverse transfer XX/NN;
   exchange of the two unobserved gluons in class15 gg -> ggg, UU/LL.
   Four-dimensional unregularized splitting kernels are used. No original
   algebra, integration or finite-assembly helper is loaded. *)
$HistoryLength = 0;
sourceRoot = Environment["FEYNFACET_SOURCE"];
outputRoot = Environment["FEYNFACET_BENCHMARK_OUTPUT"];
If[!And @@ (StringQ /@ {sourceRoot, outputRoot}),
 Print["Set FEYNFACET_SOURCE and FEYNFACET_BENCHMARK_OUTPUT."];Quit[2]];
paths = FileNameJoin[{sourceRoot,"Processes",#[[1]],"Results","FeynCalc","class-"<>#[[2]],"Algebra.wl"}] & /@
 {{"Hadronic","1"},{"HadronicTransfer","1"},{"Hadronic","15"}};
inputHashes=AssociationThread[paths,IntegerString[FileHash[#,"SHA256"],16,64]&/@paths];
{ordinary,transfer,gluon}=Get/@paths;
checks={};failures={};values=<||>;
If[!DirectoryQ[outputRoot],CreateDirectory[outputRoot,CreateIntermediateDirectories->True]];
record[name_,test_]:=(AppendTo[checks,<|"name"->name,"passed"->TrueQ[test]|>];
 If[!TrueQ[test],AppendTo[failures,name];Print["FAIL ",name]];
 Export[FileNameJoin[{outputRoot,"progress.json"}],<|"checks"->checks,"failures"->failures|>,"RawJSON"];TrueQ[test]);
zeroQ[expression_]:=TrueQ[Cancel[Together[expression]]===0];
(* Rational valuation reads the actual leading power of a combined
   numerator/denominator. A stronger pole is returned and rejected; it is
   never discarded when extracting the expected constant. *)
rationalLeading[expression_, variable_] := Module[{rational, numerator, denominator, n, d, order, leading},
 rational=Together[expression];
 If[rational===0,Return[<|"Order"->Infinity,"Value"->0|>]];
 numerator=Numerator[rational];denominator=Denominator[rational];
 If[!PolynomialQ[numerator,variable] || !PolynomialQ[denominator,variable],
  Return[Failure["NonRationalLimit",<||>]]];
 n=Exponent[numerator,variable,Min];d=Exponent[denominator,variable,Min];
 If[!IntegerQ[n] || !IntegerQ[d],Return[Failure["UnsupportedValuation",<||>]]];
 order=n-d;
 leading=Coefficient[numerator,variable,n]/Coefficient[denominator,variable,d];
 <|"Order"->order,"Value"->If[order>0,0,leading variable^order]|>];
supportedLeading[data_] := AssociationQ[data] &&
 TrueQ[data["Order"]===Infinity || IntegerQ[data["Order"]] && data["Order"]>=0];
control=rationalLeading[(2+3x)/(5+x),x];
record["valuation/control-finite",supportedLeading[control] && control["Value"]===2/5];
control=rationalLeading[(1+x)/x^2,x];
record["valuation/control-stronger-pole",AssociationQ[control] && control["Order"]===-2 && !supportedLeading[control]];
control=rationalLeading[x^2/(2+x),x];
record["valuation/control-vanishing",supportedLeading[control] && control["Order"]===2 && control["Value"]===0];
control=rationalLeading[(1+x)/x-1/x,x];
record["valuation/control-cancellation",supportedLeading[control] && control["Order"]===0 && control["Value"]===1];
control=rationalLeading[(x-x)/(1+x),x];
record["valuation/control-zero-polynomial",supportedLeading[control] && control["Value"]===0];

If[!record["input/identical-gluon-channel",gluon["Spec"]["Name"]==="g g -> g"],Quit[1]];
colourRules={CA->Nc,CF->(Nc^2-1)/(2Nc)};
physical4D={D->4,spe3->0};

(* Fierz derivation in colour flows {D=delta31 delta42,X=delta32 delta41}.
   M=(X-D/Nc)/2, and the flow Gram matrix is {{Nc^2,Nc},{Nc,Nc^2}}.
   Both incoming colour charges carry a minus sign; incoming/outgoing
   products have one minus sign. These matrices follow directly from
   T^a_ij T^a_kl=(delta_il delta_kj-delta_ij delta_kl/Nc)/2. *)
flowGram={{Nc^2,Nc},{Nc,Nc^2}};
bornColour={-1/Nc,1}/2;
charge12={{-1/Nc,1},{1,-1/Nc}}/2;
charge13={{-(Nc^2-1)/(2Nc),-1/2},{0,1/(2Nc)}};
norm=Cancel[bornColour.flowGram.bornColour];
r12=Cancel[(bornColour.flowGram.charge12.bornColour)/norm];
r13=Cancel[(bornColour.flowGram.charge13.bornColour)/norm];
r14=Cancel[-(Nc^2-1)/(2Nc)-r12-r13];
record["colour/Born-norm",zeroQ[norm-(Nc^2-1)/4]];
record["colour/T1T2",zeroQ[r12+1/Nc]];
record["colour/T1T3",zeroQ[r13-1/(2Nc)]];
record["colour/conservation",zeroQ[r12+r13+r14+(Nc^2-1)/(2Nc)]];

rho=qa+qb-qc;
softMap={u->-s-t+lambda rho,A->-s-t+lambda qa,
 B->t-lambda rho+lambda qb,C->s-lambda rho-lambda qc};
record["soft/phase-space-constraint",zeroQ[(A+B+C+s+t+u)/.softMap]];
(* qi=2pi.k/lambda. The squared eikonal current is
   -2gs^2 Sum_{i<j}rij*pi.pj/(pi.k pj.k), yielding the factor -4 below. *)
softKernel=-4gs^2(r12 s/(qa qb)+r13(-t)/(qa qc)+r14(s+t)/(qa rho)+
 r14(s+t)/(qb qc)+r13(-t)/(qb rho)+r12 s/(qc rho));

(* The emitted gluon k2 approaches (1-w)pa and delta=-2pa.k2.
   The reduced incoming momentum is w pa; the initial-state splitting
   density has the flux/spin-density factor1/w and denominator-delta.
   Squared helicity-conserving splitting amplitudes give1+w^2; their
   transverse-helicity interference gives2w, both divided by1-w. *)
collinearMap={A->-(s+t)-delta,B->w t,C->s w+delta,u->-(s+t)w};
record["collinear/phase-space-constraint",zeroQ[(A+B+C+s+t+u)/.collinearMap]];
record["collinear/virtuality",zeroQ[((-s-t-A)/.collinearMap)-delta]];

Do[
 data=If[index<=2,ordinary,transfer];spin={"UU","LL","XX","NN"}[[index]];
 real=Expand[(data["Real"][spin]/.physical4D)/.colourRules];
 born=(data["Born"][spin]/.physical4D)/.colourRules;
 Print["SOFT ",spin," ",DateString[]];
 (* Hard-scale homogeneity sets s=1,t=-v. The arbitrary normalization
    of the soft reference momentum is absorbed in lambda, allowing qa=1;
    qb/qa and qc/qa remain generic. No soft direction is sampled. *)
 softUnits={s->1,t->-v,qa->1,qb->softB,qc->softC};
 softData=rationalLeading[lambda^2(real/.softMap)/.softUnits,lambda];
 record[spin<>"/soft-leading-order",supportedLeading[softData]];
 If[!AssociationQ[softData],Quit[1]];softLeading=softData["Value"];
 expected=softKernel born/.softUnits;
 residual=Cancel[Together[softLeading-expected]];
 record[spin<>"/soft",zeroQ[residual]];
 record[spin<>"/control/soft-colour",!zeroQ[softLeading-(expected+(gs^2 born s/(qa qb)/.softUnits))]];
 Print["COLLINEAR ",spin," ",DateString[]];
 collinearData=rationalLeading[delta(real/.collinearMap)/.{s->1,t->-v},delta];
 record[spin<>"/collinear-leading-order",supportedLeading[collinearData]];
 If[!AssociationQ[collinearData],Quit[1]];collinearLeading=collinearData["Value"];
 bornReduced=(born/.{s->w s,t->w t})/.{s->1,t->-v};
 splitting=If[index<=2,(Nc^2-1)/(2Nc)(1+w^2)/(1-w),
  (Nc^2-1)/(2Nc)2w/(1-w)];
 collinearExpected=-2gs^2/w splitting bornReduced;
 collinearResidual=Cancel[Together[collinearLeading-collinearExpected]];
 record[spin<>"/collinear",zeroQ[collinearResidual]];
 record[spin<>"/control/collinear-kernel",!zeroQ[collinearLeading-(collinearExpected+gs^2 bornReduced)]];
 values[spin]=<|"SoftResidual"->residual,"CollinearResidual"->collinearResidual|>,
 {index,4}];

(* For k1+k2=pa+pb-p, C'=-t-u-C=s-s4-C, not s-C. *)
exchangeMap={A->-(s+t)-A,B->-(s+u)-B,C->-t-u-C};
record["exchange/constraint",zeroQ[(A+B+C+s+t+u)/.exchangeMap/.C->-(s+t+u)-A-B]];
record["exchange/involution",And@@(zeroQ[#]&/@((({A,B,C}/.exchangeMap)/.exchangeMap)-{A,B,C}))];
record["exchange/identical-gluon-weight",Cases[gluon["Spec"]["Real"],Rule[{"g","g","g","g","g"},x_]:>x]==={1/2}];
(* These deterministic exact-rational probes preserve the phase-space
   constraint on both sides. They are sampled checks, not symbolic identities. *)
exchangePoints={
 {s->17,t->-7,u->-6,A->-2,B->-5,C->3,Nc->3,gs->1},
 {s->23,t->-8,u->-9,A->-4,B->-6,C->4,Nc->5,gs->1},
 {s->31,t->-12,u->-11,A->-5,B->-9,C->6,Nc->7,gs->1}};
wrongPoint={s->29,t->-11,u->-10,A->-3,B->-9,C->4,Nc->3,gs->1};
finiteRational[value_] := MatchQ[value,_Integer|_Rational];
Do[
 Print["EXCHANGE ",spin," ",DateString[]];
 (* Nf=0 extracts the three-gluon real sector, including its weight1/2. *)
 real=(gluon["Real"][spin]/.physical4D/.Nf->0)/.colourRules;
 sampleValues={};
 Do[
  point=exchangePoints[[probe]];
  swappedPoint=Join[Thread[{A,B,C}->(({A,B,C}/.exchangeMap)/.point)],point];
  left=Cancel[Together[real/.point]];right=Cancel[Together[real/.swappedPoint]];
  record[spin<>"/exchange-probe-"<>ToString[probe]<>"/finite",
   finiteRational[left] && finiteRational[right] &&
    ((A+B+C+s+t+u)/.point)===0 && ((A+B+C+s+t+u)/.swappedPoint)===0];
  record[spin<>"/exchange-probe-"<>ToString[probe],left===right];
  AppendTo[sampleValues,<|"Point"->point,"SwappedPoint"->swappedPoint,"Left"->left,"Right"->right|>],
  {probe,Length[exchangePoints]}];
 wrongMap={A->-(s+t)-A,B->-(s+u)-B,C->s-C};
 wrongSwappedPoint=Join[Thread[{A,B,C}->(({A,B,C}/.wrongMap)/.wrongPoint)],wrongPoint];
 left=Cancel[Together[real/.wrongPoint]];wrong=Cancel[Together[real/.wrongSwappedPoint]];
 record[spin<>"/control/wrong-exchange-map",
  finiteRational[left] && finiteRational[wrong] && left=!=wrong];
 values[spin<>"-exchange"]=<|"Samples"->sampleValues,"WrongMapControl"->{left,wrong}|>,
 {spin,{"UU","LL"}}];
record["input-integrity",inputHashes===AssociationThread[paths,IntegerString[FileHash[#,"SHA256"],16,64]&/@paths]];
If[!DirectoryQ[outputRoot],CreateDirectory[outputRoot,CreateIntermediateDirectories->True]];
Put[<|"Coverage"->"class1 soft/collinear UU LL transferXX NN; class15 gg->ggg sampledexchangeUU LL; D=4", "Values"->values,"Checks"->checks,"Inputs"->inputHashes|>,FileNameJoin[{outputRoot,"limits.wl"}]];
Export[FileNameJoin[{outputRoot,"summary.json"}],<|"passed"->(failures==={}),"checkCount"->Length[checks],"checks"->checks,"failures"->failures,"inputs"->inputHashes,"dimension"->4,"exchangeComparison"->"Three deterministic exact-rational points, Nc=3,5,7","checkerSHA256"->IntegerString[FileHash[$InputFileName,"SHA256"],16,64]|>,"RawJSON"];
Print["UNRESOLVED LIMITS: ",Length[checks]," checks, ",Length[failures]," failures"];
Quit[If[failures==={},0,1]];
