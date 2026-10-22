$HistoryLength = 0;
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "FeynFacet", "AlgebraReuse.wl"}]];
passed = 0;
check[name_, condition_] := If[TrueQ[condition], passed++; Print["PASS ", name], Print["FAIL ", name]; Exit[1]];
grouped = GroupIdenticalAlgebraTerms[a tensor[x] + b tensor[x] + c tensor[y], _tensor];
check["identical structures share a coefficient", Length[grouped] === 2 &&
 Total[(#["Coefficient"] #["Structure"]) & /@ grouped] === (a+b) tensor[x]+c tensor[y]];
check["rational cancellation before evaluation",
 GroupIdenticalAlgebraTerms[(1/(s-t)-(s+t)/(s^2-t^2)) tensor[x], _tensor] === {}];
check["zero yields no work", GroupIdenticalAlgebraTerms[0, _tensor] === {}];
check["noncommutative order is retained", With[{rows =
 GroupIdenticalAlgebraTerms[a word[x,y]+b word[y,x], _word]}, Length[rows] === 2]];
check["failure propagates out of grouping", FailureQ[GroupIdenticalAlgebraTerms[Failure["bad",<||>], _tensor]]];
evaluations = 0;
cache = CreateBoundedAlgebraCache[Function[input, evaluations++; input^2], 2, 1000000];
check["initial cache evaluation", cache["Apply"][{"Real",x},x] === x^2];
second = cache["Apply"][{"Real",x},x];
check["identical key is reused", second === x^2 && evaluations === 1];
composed = 1 + cache["Apply"][{"Real",x},x];
check["cache hit remains a value inside arithmetic", composed === 1+x^2 && FreeQ[composed,_Return]];
multiplied = Expand[(a+b) cache["Apply"][{"Real",x},x]];
check["cache hit remains a value under Expand", multiplied === a x^2+b x^2 && FreeQ[multiplied,_Return]];
check["kinematics separate keys", cache["Apply"][{"Born",x},x] === x^2 && evaluations === 2];
cache["Apply"][{"Real",y},y];
check["entry bound is enforced", cache["Statistics"][]["Entries"] <= 2 && cache["Statistics"][]["Resets"] === 1];
check["failed input is rejected", FailureQ[cache["Apply"]["bad",Failure["bad",<||>]]]];
bad = CreateBoundedAlgebraCache[Function[input, Failure["bad",<||>]]];
check["failed evaluator propagates", FailureQ[bad["Apply"]["x",1]] && bad["Statistics"][]["Entries"] === 0];
check["invalid bounds are rejected", FailureQ[CreateBoundedAlgebraCache[Identity,0,1000]]];
Print[passed, " algebra reuse tests passed."];
Exit[0];
