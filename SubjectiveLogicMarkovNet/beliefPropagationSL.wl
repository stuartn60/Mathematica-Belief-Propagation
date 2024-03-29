(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["beliefPropagationSL`"];
(*\[Copyright] Copyright Stuart Nettleton 2020, 2021*)
(*Disclaimer: This software is provided without any warranty express or implied Stuart Nettleton 2020*)
(*Requires Wolfram Mathematica 12.1.0.0 or higher*)
(*All code files and executable documents are with the license GPL 3.0.For details see http://www.gnu.org/licenses/.*)
(*All documents are with the license Creative Commons Attribution 4.0 International (CC BY 4.0).For details see https://creativecommons.org
licenses/by/4.0/.*)
(*As amended based on original Markov Net schema by Daphne Koller, Stanford University Subject Probabilistic Graphical Models, Copyright (C) Daphne Koller, Stanford University, 2012*)

ind2AssSL::usage ="ind2AssSL[card] calculates an structural assignment for a index of a factor with cardinality card";

ass2IndSL::usage ="ass2IndSL[card,assignments] calculates an index for a structural assignments of a factor with cardinality card";

convertFactorToDirichletSL::usage="convertFactorToDirichletSL[f,uArray:0.15,bArray:1,isLogsInFactor:False] converts probabilities to Dirichlet beliefs:
 ... f is a group of factors in standard Bayes/Markov network order (the first variable varies most frequently & the last variable varies the least frequently)
 ... output retains probabilities in Bayes/Markov order and appends probabilities & beliefs in usual array order, uncertainty, and base rates in usualarray order
 ... uArray is the Uncertainty to apply ... default of 0.15 ... can be a single number to apply to all factors or a vector with length of the number of factors
 ... bArray is the Base Rates to apply ... default of 1 for an equal vector for each factor ... can be a vector with length of the number of factors containing a base rate vector for each factor";

factorProductOrSumSL::usage ="factorProductOrSumSL[f1,f2,isMax:False] calculates the product of two factors where
 ... each factor is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms";

factorMarginalizationSL::usage ="factorMarginalizationSL[a,v,isMax_] removes variables v from a factor a by summing-out ... 
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) which finds the maximum factor value over all possible assignments to the marginalized variables and leaves the output as unnormalised natural logarithms";

(*reorderFactorSL::usage="reorderFactorSL[f1,isMax] ... reoders a factor so variables are in increasing order
 ... f1 is a factor
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms";*)

observeEvidenceSL::usage ="observeEvidenceSL[f,e,emethod] modifies a set of factors (f) for observed evidence (e) where
 ... each factor is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... evidence is a list of {variable, evidence value}";

computeJointDistributionSL::usage ="computeJointDistributionSL[f] calculates the joint distribution of a set of factors (f) by cumulative SumProduct of the fators
 ... each factor is a list containing the list of factor variables, list of variable cardinalities & list of variable values";

computeMarginalSL::usage ="computeMarginalSL[v,f,e,emethod] calculates the marginal over variables (v) in the factors (f) with observed evidence (e)
 ... each factor is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... evidence is a list of {variable, evidence value}";

computeInitialPotentialsSL::usage ="computeInitialPotentialsSL[cliques_] calculates the standard form of a clique tree (cliques) where
 ... the clique tree is a list contining the adjacency matrix of the tree and the factors
 ... assigns each factor to a clique
 ... computes factor products over all factors assigned to each clique";

pruneTreeSL::usage ="pruneTreeSL[cliques] compacts the clique tree (cliques) by removing redundant clique subsets while maintaining the supersets and the running intersection property";

createCliqueTreeSL::usage ="createCliqueTreeSL[f,e,emethod,method:'Kruskal',print:False] calculates the clique tree for the variables in f given the evidence e
 ... f is the factor list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value}
 ... reurned cliques have the same structure as factors
 ... the Maximum Spanning Tree Method may be the default 'Kruskal' for undirected graphs or 'Prim'
 ... if print is True then the node-clique and the clique-Maximum Spanning Tree graphs are printed";

printFactorAndCliqueGraphsSL::usage="printFactorAndCliqueGraphsSL[f,e,emethod,method:'Kruskal'] creates & displays factor & clique graphs
 ... identical to createCliqueTree except only prints the graphs
 ... f is the factor list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value}
 ... returned cliques have the same structure as factors
 ... the Maximum Spanning Tree Method may be the default 'Kruskal' for undirected graphs or 'Prim'";

getNextCliquesSL::usage ="getNextCliquesSL[p,messages] finds a pair of cliques {i,j} where clique i is ready to transmit a message to clique j
 ... the clique tree (p) 
 ... current messages (messages)
 ... a clique i is ready to transmit to its neighbor j when i has received messages from all of its neighbors except from j
 ... each message is passed only once
 ... where more than one message is ready to be transmitted the numerically smallest pair is returned";

cliqueTreeCalibrateSL::usage = "cliqueTreeCalibrateSL[p,isMax] calculates the final potentials for a clique tree (p)
 ... two adjacent cliques i and j are calibrated if they agree on the marginals over their shared variables
... if isMax is False the Sum-Product algorithm is used for message passing in the clique tree and the messages are normalize such that the values in the message sum to 1
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms";

computeMarginalsBPSL::usage="computeMarginalsBPSL[cliques,isMax] is a subfunction for computing marginals
 ... cliques are the input cliques
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms";

computeExactMarginalsBPSL::usage ="computeExactMarginalsBPSL[f,e:{},emethod,isMax:False,method:'Kruskal',print:False] calculates the exact posterior marginal probability distribution for each variable in the initial factors f given the evidence e
 ... f is the clique list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value} 
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms
 ... the Maximum Spanning Tree Method may be the default 'Kruskal' for undirected graphs or 'Prim'
 ... if print is True then the node-clique and the clique-Maximum Spanning Tree graphs are printed";

maxDecodedSL::usage ="maxDecodedSL[marginals] calculates the best assignment and alphabetic equivalent for each Max-Product marginal variable";

naiveGetNextClustersSL::usage ="naiveGetNextClustersSL[p,m] finds a pair of clusters {i,j} where cluster i is ready to transmit the m+1 message to its neighbour cluster j where
 ... the cluster graph (p)
 ... current messages (messages)
 ... the method iterates over the messages (cluster pairs) in increasing order where messages are sorted in ascending ordered by their destination index & ties are broken based on the origin index. If m is 0, [i j] will be the pair of clusters with the lowest j value & (of those pairs over this j) lowest i value as this is the 'first' element in our ordering
 ... where more than one message is ready to be transmitted the numerically smallest pair is returned";

getNextClustersSL::usage ="getNextClustersSL[p,messages,oldMessages,m,useSmart] finds a pair of nodes {i,j} where node i is ready to transmit the m+1 message to node j
 ... the cluster graph (p)
 ... current messages (messages)
 ... oldMessages[[i,j]] contains the value that Messages(i,j) contained immediately before it was updated to its current value
 ... m is the index of the message to be passed
 ... the method iterates over the messages in increasing order where messages are sorted in ascending ordered by their destination index & ties are broken based on the origin index. If m is 0, [i j] will be the pair of clusters with the lowest j value & (of those pairs over this j) lowest i value as this is the 'first' element in our ordering
 ... where more than one message is ready to be transmitted the numerically smallest pair is returned
 ... if useSmartMP is True then an enhanced apporach is used to obtain the next message cluster for processing";

smartGetNextClustersSL::usage ="smartGetNextClustersSL[p,messages,oldMessages,m] finds a pair of nodes {i,j} where node i is ready to transmit the m+1 message to clique j
 ... the cluster graph (p)
 ... current messages (messages)
 ... oldMessages[[i,j]] contains the value that Messages(i,j) contained immediately before it was updated to its current value
 ... m is the index of the message to be passed
 ... the method iterates over the messages in increasing order where messages are sorted in ascending ordered by their destination index & ties are broken based on the origin index. If m is 0, [i j] will be the pair of clusters with the lowest j value & (of those pairs over this j) lowest i value as this is the 'first' element in our ordering
 ... where more than one message is ready to be transmitted the numerically smallest pair is returned
 ... other message passing schedules can further improve the convergence of the algorithm";

createClusterGraphSL::usage ="createClusterGraphSL[f,e,emethod] calculates an adjacency matrix and list of clusters for a Bethe cluster graph with nodes representing single variable clusters or pairwise clusters for the variables in f given the evidence e where:
 ... a Bethe cluster graph is the simplest way to create a cluster graph from a network as a separate cluster is created for each factor and attached to form the cluster graphdd
 ... f is the cluster list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value} 
 ... each cluster in the returned cluster list contains a list of factor variables, list of variable cardinalities & list of variable values initialized to the initial potential";

checkConvergenceSL::usage ="checkConvergenceSL[new_,old_,threshhold:10^-6] determines whether iterative message passing has converged to a threshhold of 10e-6
 ... new and old are message matrices with respective values of the message from cluster i to cluster j
 ... each message in the matrix is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... default threshhold for convergence is 10^-6";

clusterGraphCalibrateSL::usage ="clusterGraphCalibrateSL[p,isMax,useSmartMP,threshhold:10^-6] returns the adjacency matrix for a cluster graph (p) and the initial potentials of the clusers using loopy belief propagation
 ... each cluster is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms
 ... if useSmartMP is True then an enhanced approach is used to obtain the next message cluster for processing
 ... default threshhold for convergence is 10^-6";

computeApproxMarginalsBPSL::usage ="computeApproxMarginalsBP[f,e:{},emethod,isMax:False,useSmartMP:False,threshhold:10^-6] calculates the approximate posterior marginal probability distribution over each variable in f given the evidence e 
 ... f is the clique list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value}
 ... if isMax is False the algorithm is SumProduct with ouput of normalized probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalized natural logarithms
 ... if useSmartMP is True then an enhanced approach is used to obtain the next message cluster for processing
 ... default threshhold for convergence is 10^-6";

combineEvidenceSL::usage="combineEvidenceSL[beliefs] combines beliefs (not probabilities) using Linear Opinion Pool Consensus & Compromise";

toDirichletsSL::usage="toDirichletsSL[beliefs] calculates {\[Alpha]x,rx,axrx,ux}
 ... beliefs are {{beliefs},uncertainty,{base rates}}";

toBeliefsSL::usage="toBeliefsSL[dirichlets] calculates {{beliefs},uncertainty,{base rates}}
 ... dirichlets are {{rx1,axrx1},{rx2,axrx2}, ...}";

convertSamiamNetworkSL::usage ="convertSamiamNetworkSL[file] converts a SAMAIM hugin .net file to factors
 ... return a list of SAMIAM node names and list of factors";

sortFactorVariablesSL::usage = "sortFactorVariablesSL[factors] sorts the numerical factors leaving the left-most the same (which is the node) and the remainder in ascending numerical order"; 

sortFactorVariablesNonSL::usage ="sortFactorVariablesNonSL[factors] sorts the non-SL numerical factors leaving the left-most the same (which is the node) and the remainder in ascending numerical order"; 

renameFactorsSL::usage ="renameFactorsSL[factors,list] renames the numerical factors and sorts the factors into canonical order
 ... list is a list of pairs whose node numbers are to be swapped e.g. {{1,2}},{3,4}}";

normalizeCliqueSL::usage="normalizeCliqueSL[clique,isMax] normalizes the clique or factor to a legal probability distribution
 ... if isMax is False the ouput is normalized probabilities
 ... if isMax is True the input is returned unchanged (i.e. unnormalized)";

adjustBeliefsSL::usage="adjustBeliefsSL[beliefs,cards] normalizes the clique or factor to a legal probability distribution";

Begin["`Private`"];

safeDivideSL[a_,b_]:=Limit[a/(b/.{0|0.->\[Epsilon]}),\[Epsilon]->0,Direction->"FromAbove"];

normalizeCliqueSL[clique_,isMax_]:=Module[{f1,fnNorm},
       fnNorm@f_:=Module[{\[ScriptCapitalD]3,\[ScriptCapitalD]7},If[!isMax,\[ScriptCapitalD]3=f[[3]]/Total@f[[3]];\[ScriptCapitalD]7=Flatten[#/Total@#&/@Partition[f[[7]],f[[2,1]]]],\[ScriptCapitalD]3=f[[3]];\[ScriptCapitalD]7=f[[7]]];
       {f[[1]],f[[2]],\[ScriptCapitalD]3,f[[4]],f[[5]],f[[6]],\[ScriptCapitalD]7}];
 If[Length@Cases[clique,{{__},{__},{__},{__},{__},__,{__}}]==0,fnNorm@clique,fnNorm/@clique]];

ind2AssSL[cards_]:=Reverse/@Flatten[Array[List,Reverse@cards],Length@cards-1];
ass2IndSL[cards_,ass_]:=Flatten[Position[ind2AssSL@cards,#]&/@Flatten[{ass},Length@Dimensions@{ass}-2]];

adjustBeliefsSL[beliefs_,cards_]:=Module[{pos,bx,ux,bx1,bx2},
    If[Length@ Cases[beliefs,{{__},__}]==0,
	bx={First@beliefs};ux={beliefs[[2]]},
	bx=First/@beliefs;ux=beliefs[[All,2]]];
    pos=Position[#<0&/@#&/@bx,True];
    bx1= If[pos=={},bx,ReplacePart[bx,Thread[pos->0]]];
    bx2=MapThread[Flatten[safeDivideSL[(1-#2)Partition[#1,First@#3],Flatten[Total/@Partition[#1,First@#3]]]]&,{bx1,ux,cards}];
    If[Length@ Cases[beliefs,{{__},__}]==0,First@bx2,bx2]];

convertFactorToDirichletSL[f_,uArray_:1,bArray_:1,isLogsInFactor_:False]:=Module[{f1,vars,cards0,cards,vals,uncert0,base0,
bArray1,indx,\[ScriptCapitalD]4,\[ScriptCapitalD]5,\[ScriptCapitalD]6,\[ScriptCapitalD]7,\[ScriptCapitalD]8,\[ScriptCapitalD]9},
  f1=If[Length@Cases[f,{{__},{__},{__}}]==0,{f},f];
  vars=f1[[All,1]];cards=f1[[All,2]];vals=f1[[All,3]];
  (*For each factor {1.vars, 2.cards, 3.vals, 4.normalised probabilities, 5.beliefs, 6.uncertainty, 7.base rates*)
  uncert0=Table[0.15,Length@f1];
  \[ScriptCapitalD]6=Which[
     MemberQ[{{},{1}},Union@Flatten@{uArray}],uncert0,(*if no uncertainty structure*)
     MemberQ[{{0},{0.}},0Union@Flatten@{uArray}],Flatten@Table[Union@Flatten@{uArray},Length@f1],
     Equal[0*Flatten@{uArray},0*uncert0],Flatten@{uArray},(*if full uncertainty structure*)
     True,Print["Error in Uncertainty specification: the Uncertainty vector has ",Length@uArray," elements while the number of factors is " ,Length@f1,". Note: specify a vector, a single uncertainty to apply to all factors, or use the default of 0.15 by either no specification or an empty vector {}."];Abort[]];
  base0=Flatten[#/Total@#&/@#]&/@(Table[1,#]&/@#&/@((Length/@Partition[ind2AssSL@#,#[[1]]])&/@cards));
  bArray1=If[Length@Cases[bArray,{__}]==0,{bArray},bArray];
  \[ScriptCapitalD]7=Which[
      MemberQ[{{},{1},{0.5}},Union@Flatten@{bArray1}],base0,(*if no base rate structures*)
      Equal[0*base0,0* bArray1],bArray1,(*if full base rate structures*)
      True,Print["Error in Base Rate specification: the Base rate vector has dimensions ",Dimensions@bArray," elements while the factor values have dimensions " ,Length@#&/@f1[[All,3]]," ... Note: specify an array of Base Rate vectors, a single Base Rate vector to apply to all factors, or use the default of equal Base Rates by either no specification or an empty vector {}."];Abort[]];
\[ScriptCapitalD]5=adjustBeliefsSL[MapThread[{Flatten[safeDivideSL[Partition[#1,#2],Total/@Partition[#1,#2]]]-#4 #3,#3}&,{If[isLogsInFactor==True,Exp/@vals,vals],cards[[All,1]],\[ScriptCapitalD]6,\[ScriptCapitalD]7}],cards];
    \[ScriptCapitalD]4=\[ScriptCapitalD]5+\[ScriptCapitalD]6 \[ScriptCapitalD]7;
   Transpose@{vars,cards,vals,\[ScriptCapitalD]4,\[ScriptCapitalD]5,\[ScriptCapitalD]6,\[ScriptCapitalD]7}];

multiplicationUncertaintySL[beliefs_] := If[Times@@beliefs[[All,2]]==0,Plus@@beliefs[[All,2]],Mean@beliefs[[All,2]]];(*Module[{b,u,a,uXYMin,uXYMax,px,ax,uXYcandidates},
	b=First/@beliefs;u=beliefs[[All,2]];a=Last/@beliefs;
	uXYMin=Times@@u;uXYMax=Min[1,Plus@@u-uXYMin];
	ax=#/Total@#&/@ReplacePart[a,Position[Last/@beliefs,0]->10^-10];
	px=#[[1]]+#[[2]]#[[3]]&/@beliefs;
	uXYcandidates=Flatten[(Outer[Times,##]&@@px-Outer[Times,##]&@@b)/Outer[Times,##]&@@ax];
	(*Min@@Pick[uXYcandidates,#>=uXYMin&&#<=uXYMax&/@uXYcandidates,True]*)
	Mean@u];*)

factorProductOrSumSL[f1_,f2_,isMax_:False]:=Module[{f,outVars,outCards,vals3,vals7,\[ScriptCapitalD]3,\[ScriptCapitalD]6,\[ScriptCapitalD]7,\[ScriptCapitalD]7a},
  (*For each factor {1.vars, 2.cards, 3.vals, 4.normalised probabilities, 5.beliefs, 6.uncertainty, 7.base rates*)
  f={f1,f2};
  {outVars,outCards}=Transpose@Union@Flatten[Transpose/@f[[All,;;2]],1];
  vals3=Map[ind2AssSL[outCards][[All,Flatten[Position[outVars,#]&/@#[[1]]]]]/.Thread[ind2AssSL@#[[2]]->#[[3]]]&,f];
  vals7=Map[ind2AssSL[outCards][[All,Flatten[Position[outVars,#]&/@#[[1]]]]]/.Thread[ind2AssSL@#[[2]]->#[[7]]]&,f];
  \[ScriptCapitalD]3=If[isMax,Plus,Times]@@@Transpose@vals3;
  \[ScriptCapitalD]6=multiplicationUncertaintySL@f[[All,5;;7]];
  \[ScriptCapitalD]7a=If[isMax,Plus,Times]@@@Transpose@vals7;
  \[ScriptCapitalD]7=If[isMax,\[ScriptCapitalD]7a,Flatten[#/Total@#&/@Partition[\[ScriptCapitalD]7a,outCards[[1]]]]];
  First@convertFactorToDirichletSL[{outVars,outCards,\[ScriptCapitalD]3},\[ScriptCapitalD]6,\[ScriptCapitalD]7,isMax]];

factorMarginalizationSL[f_,mv_:{},isMax_:False]:=Module[{pos0,pos1,outVars,outCards,colOrder,\[ScriptCapitalD]3,\[ScriptCapitalD]6,\[ScriptCapitalD]7,\[ScriptCapitalD]7a},
    If[ Or[First @ f == {},Length@f[[1]]==1],Return@f];
    pos0=Flatten[Position[f[[1]],#]&/@mv,1];
    pos1=Complement[Thread[{Range@Length@f[[1]]}],pos0];
    If[pos1 == {},Print["Error:Resultant factor has empty scope"];Abort[]];
    {outVars,outCards}=Transpose@SortBy[Transpose@{Extract[f[[1]],pos1],Extract[f[[2]],pos1]},First];
   colOrder=Flatten[Position[outVars,#]&/@Select[f[[1]],!MemberQ[mv,#]&]];
    \[ScriptCapitalD]3=ind2AssSL[outCards][[All,colOrder]]/.Normal@GroupBy[Thread[Delete[#,pos0]&/@ind2AssSL@f[[2]]->f[[3]]],First->Last,If[isMax,Max,Plus]@@#&];
    \[ScriptCapitalD]7a=ind2AssSL[outCards][[All,colOrder]]/.Normal@GroupBy[Thread[Delete[#,pos0]&/@ind2AssSL@f[[2]]->f[[7]]],First->Last,If[isMax,Max,Plus]@@#&];
    \[ScriptCapitalD]6=f[[6]];
    \[ScriptCapitalD]7=If[isMax,\[ScriptCapitalD]7a,Flatten[#/Total@#&/@Partition[\[ScriptCapitalD]7a,outCards[[1]]]]];
    First@convertFactorToDirichletSL[{outVars,outCards,\[ScriptCapitalD]3},\[ScriptCapitalD]6,\[ScriptCapitalD]7,isMax]];

observeEvidenceSL[f_,e_,emethod_:1]:=Module[{f1,e1,pos1,pos2},
  If[MemberQ[{{},{{}}},e],Return@f,f1=f];
  e1=If[Length@Cases[e,{{_},{_},{__},{__},{__},_,{__}}]==0,{e},e];
  Switch[emethod,
     1,(*Koller Factor Reduction Method ... General Use ... Evaluates Evidence Nodes*)
              pos1=Position[First/@f1,#]&/@First/@First/@e1;
              pos2=Flatten[MapIndexed[Thread[{#1,First@#2}]&,First/@#&/@pos1],1];
              Scan[(f1[[First@#]]=factorProductOrSumSL[f1[[First@#]],e1[[Last@#]]])&,pos2],

       2,(*Modified Koller Factor Reduction ... Fixes Evidence on Leaf Nodes ... Restricted to Leaf Node Evidence*)
             f1=Select[f1,!MemberQ[e1[[All,1]],#[[1]]]&];
             pos1=Position[First/@f1,#]&/@First/@First/@e1;
             pos2=Flatten[MapIndexed[Thread[{#1,First@#2}]&,First/@#&/@pos1],1];
             Scan[(f1[[First@#]]=factorProductOrSumSL[f1[[First@#]],e1[[Last@#]]])&,pos2];
             Scan[(f1[[First@#]]=factorMarginalizationSL[f1[[First@#]],e1[[Last@#]][[1]]])&,pos2];
            f1=SortBy[
                        If[Length@Dimensions@#==1,#,Flatten[#,1]]&/@(Switch[Length@#,1,#,2,factorProductOrSumSL@@#,_,Fold[factorProductOrSumSL[#1,#2]&,First@#,Rest@#]]&/@GatherBy[Join[f1,e1],First]),
                    #[[1,1]]&]
      ];f1];

computeJointDistributionSL[f_]:=Which[f=={},Print["Error:empty factor list"];{},
	Length@f==1||Length@Dimensions@f==1,sortFactorVariablesSL@f,
	True,Fold[factorProductOrSumSL[#1,#2,False]&,First@f,Rest@f]];

computeMarginalSL[v_,f_,e_,emethod_:1]:=Module[{fe,joint,marginal0},
	If[f=={},Print["Error:empty factor list"];Abort[]];
	fe=If[e=={},f,observeEvidenceSL[f,If[Length@Dimensions@e==1,{e},e],emethod]];
	joint=computeJointDistributionSL@fe;
	marginal0=factorMarginalizationSL[joint,Complement[First@joint,v],False];
	normalizeCliqueSL[#,False]&/@If[Length@Dimensions@marginal0==1,{marginal0},marginal0]];

computeInitialPotentialsSL[cliqueSet_,method_:"Kruskal",print_:False]:=Module[{cliqueList, factorList,newCliques,spanningTree,potentials},
       cliqueList = cliqueSet[[1]];
       factorList = cliqueSet[[-1]];  newCliques=Sort@Normal@GroupBy[Flatten[Thread[#[[1]]->#[[2]]]&/@MapIndexed[{Min@@Intersection@@(Position[cliqueList,#][[All,1]]&/@#1),First@#2}&,factorList[[All,1]]]],First->Last];
       spanningTree=maxSpanningTreeSL[cliqueList[[newCliques[[All,1]]]],method,print];
       potentials =Which[Length@#==1,factorList[[First@#]],True,computeJointDistributionSL@factorList[[#]]]&/@newCliques[[All,2]];
      {spanningTree, potentials}]; 

pruneTreeSL[cliques_]:=Module[{cnodes,cedges,neighboursI,toRemove,toKeep,i,j,nk},
	toRemove={};cnodes=First@cliques;cedges=cliques[[2]];
	Do[If[MemberQ[toRemove,i],Continue[]];
	  neighboursI=Flatten@Position[cedges[[i]],1];
	 Do[j=neighboursI[[neighbour]];
	    If[MemberQ[toRemove,j],Continue[]];
	    If[Length@Select[MemberQ[cnodes[[j]],#]&/@cnodes[[i]],#==True&]==Length@cnodes[[i]],
		Do[If[Length@Intersection[cnodes[[i]],cnodes[[nk]]]==Length@cnodes[[i]],
		cedges[[Complement[neighboursI,{nk}],nk]]=1;
		cedges[[nk,Complement[neighboursI,{nk}]]]=1;
		Break[]],
	  {nk,neighboursI}];
	 cedges[[i,All]]=0;cedges[[All,i]]=0;
	  toRemove=Append[toRemove,i]],{neighbour,Length@neighboursI}],
	{i,Length@cnodes}];
        toKeep=Complement[Range@Length@cnodes,toRemove];
        cedges[[toKeep,toKeep]]];

createCliqueTreeSL[f_,e_,emethod_,method_:"Kruskal",print_:False]:=Module[{f1,nodeEdges, nodeGraph, cliques0,
          cliques, cliqueNodes,c},
         f1=If[e == {},f ,observeEvidenceSL[f, e,emethod]];
	nodeEdges =DeleteDuplicates[UndirectedEdge@@#&/@Sort/@Flatten[Subsets[#,{2}]&/@First /@ f1,1]];
	nodeGraph = Graph[Union @@ First /@ f1, nodeEdges, VertexLabels-> "Name"];
	cliques = Sort@FindClique[nodeGraph, Infinity, All];
	cliqueNodes=Array[c,Length@cliques];
	If[print,Print[Graph[Join[Union @@ First /@ f1, cliqueNodes], Flatten@MapThread[Thread[#1 \[UndirectedEdge] #2]&, {cliqueNodes, cliques}], VertexLabels-> "Name"]]];
         computeInitialPotentialsSL[{cliques,f1(*If[e == {},f ,observeEvidenceSL[f, e,emethod]]*)},method,print]]; 

maxSpanningTreeSL[cliques_,method_:"Kruskal",print_:False]:=Module[{ cliqueNodes,c, cliqueEdges,cliqueGraph, lengthSepSets,spanningTree},
         cliqueNodes=Array[c,Length@cliques];
	cliqueEdges = UndirectedEdge@@#&/@Pick[Subsets[cliqueNodes,{2}],Not@SameQ[#,{}]&/@(Intersection@@#&/@Subsets[cliques,{2}])];
	cliqueGraph = Graph[cliqueNodes, cliqueEdges, VertexLabels ->"Name"];
	lengthSepSets = Length @ (Intersection @@ {Extract[cliques, First @ Position[cliqueNodes, First @ #]],
	    Extract[cliques, First @ Position[cliqueNodes, Last @ #]]})& /@ cliqueEdges;
	spanningTree=FindSpanningTree[cliqueGraph, EdgeWeight -> -1*lengthSepSets, Method ->method];
         If[print,Print[HighlightGraph[cliqueGraph, spanningTree, GraphHighlightStyle-> "Thick "]]];
         pruneTreeSL[{cliques, Normal @ AdjacencyMatrix @ spanningTree}]];

getNextCliquesSL[p_,messages_]:=Module[{edges,requiredEvidence,providedEvidence,lackingEvidence,availableEvidence,reByLe,aeByLe,nextMsg},
	edges=Last@p;
	requiredEvidence=Position[edges,1];
	providedEvidence=Position[Thread[#=={}]&/@#&/@messages,False];
	lackingEvidence=Complement[requiredEvidence,providedEvidence];
	availableEvidence=Position[Thread[#=={}]&/@#&/@messages,False];
	reByLe=Table[Complement[Select[requiredEvidence,Last@#==First@i&],{Reverse@i}],{i,lackingEvidence}];
	aeByLe=Table[Complement[Select[availableEvidence,Last@#==First@i&],{Reverse@i}],{i,lackingEvidence}];
	nextMsg=Pick[lackingEvidence,MapThread[#1==#2&,{aeByLe,reByLe}]];
	If[nextMsg=={},{0,0},First@nextMsg]];

cliqueTreeCalibrateSL[p_,isMax_]:=Module[{i,j,k,m,cliques,cliqueEdges,fp,requiredEvidence,messages,cliqueNodes,cliqueGraph,leafNodes,edges,sum},
	cliqueEdges=First@p;
         cliques=If[isMax==True,Join[#[[1;;2]]&/@Last@p,N@{#}/.{Indeterminate->-\[Infinity]}&/@(Log/@#[[3]]&/@Last@p),#[[4;;]]&/@Last@p,2],Last@p];
	messages=Table[{},Length@cliques,Length@cliques];
	cliqueNodes=Range@Length@cliques;
         cliqueGraph=AdjacencyGraph[cliqueNodes,cliqueEdges,VertexLabels->"Name"];
         leafNodes=Pick[cliqueNodes,#==1&/@VertexDegree@cliqueGraph];
	Scan[(messages[[#,First@First@Position[cliqueEdges[[#]],1]]]=normalizeCliqueSL[cliques[[#]],isMax])&,leafNodes];
	While[True,
	  {i,j}=getNextCliquesSL[{cliques,cliqueEdges},messages];
	   If[{i,j}=={0,0},Break[]];
	   edges=Complement[Flatten@Position[cliqueEdges[[i]],1],{j}];
	   fp=Fold[factorProductOrSumSL[#1,#2,isMax]&,cliques[[i]],messages[[edges,i]]];
	   messages[[i,j]]=normalizeCliqueSL[factorMarginalizationSL[fp,Complement[First@fp,cliques[[i,1]]],isMax],isMax];
         ];
	Do[edges=Complement[Flatten@Position[cliqueEdges[[k]],1],{k}];
	  fp=Fold[factorProductOrSumSL[#1,#2,isMax]&,cliques[[k]],messages[[edges,k]]];
	  cliques[[k]]=factorMarginalizationSL[fp,Complement[First@fp,cliques[[k,1]]],isMax],
	{k,Length@messages}];
	{cliqueEdges,cliques}];

computeMarginalsBPSL[cliques_,isMax_]:=Module[{vars,pos1,pos2},
        vars=Union@Flatten@(First/@cliques);
        pos1=Position[cliques[[All,1]],#][[1,1]]&/@vars;
        pos2=MapThread[Complement,{cliques[[pos1,1]],{#}&/@vars}];
       normalizeCliqueSL[MapThread[If[#2=={},#1,factorMarginalizationSL[#1,#2,isMax]]&,{cliques[[pos1]],pos2}],isMax]];

computeExactMarginalsBPSL[f_,e_:{},emethod_:1,isMax_:False,method_:"Kruskal",print_:False]:=computeMarginalsBPSL[Last@cliqueTreeCalibrateSL[createCliqueTreeSL[f,e,emethod,method,print],isMax],isMax];

maxDecodedSL[marginals_]:=Module[{res},
	res=Fold[Join[#1,First@Position[#2,Max@#2]]&,{},Last/@marginals];
	{res,StringJoin[res/.Thread[Range@26->Alphabet[]]]}];

(*LOOPY BELIEF PROPAGATION*)
naiveGetNextClustersSL[p_,m_]:=Module[{pos,m1},
	pos=Position[First@p,1];
	m1=Mod[m+1,Length@pos];
	If[m1==0,m1=Length@pos];
	Reverse@pos[[m1]]];

getNextClustersSL[p_,messages_,oldMessages_,m_,useSmart_]:=If[useSmart==False,
	naiveGetNextClustersSL[p,m],smartGetNextClustersSL[p,messages,oldMessages,m]];

smartGetNextClustersSL[p_,messages_,oldMessages_,m_]:=Module[{pos,m1,from,to},
	pos=Position[First@p,1];
	m1=Mod[m+1,Length@pos]+1;
	Reverse@pos[[m1]]];

createClusterGraphSL[f_,e_,emethod_]:=Module[{f1,clusterNodes,factorNodes,edges0,edges,ajm},
	f1=observeEvidenceSL[f,e,emethod];
	clusterNodes=Range@Length@f1;
	edges=Flatten@MapThread[Thread[#1\[UndirectedEdge]#2]&,Transpose@Select[Transpose@{clusterNodes,First/@f1},Length@#[[2]]>1&]];
	{Normal@AdjacencyMatrix@Graph[clusterNodes,edges,VertexLabels->"Name"],f1}];

checkConvergenceSL[new_,old_]:=If[Or@@Thread[Flatten@Abs[new-old]<10^-6],True,False];

clusterGraphCalibrateSL[p_,isMax_,useSmartMP_,printIterations_]:=Module[{clusterEdges,clusters,vars,card,mVars,mCards,mVals,mProbabilities,
         mBeliefs,mUncertainty,mBaseRates,mrx,maxrx,messages,tic,iter,i,j,lastMessages,prevMessage,edges,fp,sum,pos,clusterList,m},
	clusterEdges=First@p;clusters=Last@p;pos=Position[clusterEdges,1];
	vars=Union@Flatten@(First/@clusters);
	card=First/@(Extract[clusters[[All,2]],#]&/@Position[First/@clusters,#]&/@vars);
	(*For each factor {1.vars, 2.cards, 3.vals, 4.normalised probabilities, 5.beliefs, 6.uncertainty, 7.base rates*)
	mVars=Table[If[clusterEdges[[i,j]]==1,Intersection[clusters[[i,1]],clusters[[j,1]]],{}],{i,Length@clusters},{j,Length@clusters}];
	mCards=Table[If[clusterEdges[[i,j]]==1,card[[mVars[[i,j]]]],{}],{i,Length@clusters},{j,Length@clusters}];
	mVals=Table[If[clusterEdges[[i,j]]==1,Table[1/Times@@mCards[[i,j]],Times@@mCards[[i,j]]],{0}],{i,Length@clusters},{j,Length@clusters}];
	mProbabilities=Table[If[clusterEdges[[i,j]]==1,Table[1/Times@@mCards[[i,j]],Times@@mCards[[i,j]]],{0}],{i,Length@clusters},{j,Length@clusters}];
	mBeliefs=Table[If[clusterEdges[[i,j]]==1,Table[1/Times@@mCards[[i,j]],Times@@mCards[[i,j]]],{0}],{i,Length@clusters},{j,Length@clusters}];
	mUncertainty=Table[0,{i,Length@clusters},{j,Length@clusters}];
	mBaseRates=Table[If[clusterEdges[[i,j]]==1,Table[1/Times@@mCards[[i,j]],Times@@mCards[[i,j]]],{0}],{i,Length@clusters},{j,Length@clusters}];
	messages=Table[{mVars[[i,j]],mCards[[i,j]],mVals[[i,j]],mProbabilities[[i,j]],mBeliefs[[i,j]],mUncertainty[[i,j]],mBaseRates[[i,j]]},
	{i,Length@clusters},{j,Length@clusters}];
	tic=AbsoluteTime[];iter=0;lastMessages=messages;
	While[True,iter++;
	  {i,j}=getNextClustersSL[p,messages,lastMessages,iter,useSmartMP];
	  prevMessage=messages[[i,j]];
	  edges=Complement[Flatten@Position[clusterEdges[[i]],1],{j}];
	  fp=Fold[factorProductOrSumSL[#1,#2,isMax]&,clusters[[i]],messages[[edges,i]]];
	  messages[[i,j]]=normalizeCliqueSL[factorMarginalizationSL[fp,Complement[First@fp,Intersection[clusters[[i,1]],clusters[[j,1]]]],isMax],isMax];
	  If[useSmartMP==True,lastMessages[[i,j]]=prevMessage];
	    If[Mod[iter,Length@pos]==0,
               If[printIterations==True,Print["LBP Messages Passed: ",iter]];
	      If[checkConvergenceSL[messages,lastMessages]==True,Break[]];
	      If[iter>1500,Print["Reached ",iter," iterations ... aborting"];Abort[]]
            ];
	     If[useSmartMP==False,lastMessages=messages]
             ];
	  If[printIterations==True,Print["Total number of messages passed: ",iter," in ",AbsoluteTime[]-tic," seconds"]];
	  clusterList=clusters;
	  Do[clusterList[[Last@pos[[m]]]]=factorProductOrSumSL[clusterList[[Last@pos[[m]]]],messages[[First@pos[[m]],Last@pos[[m]]]],isMax],{m,Length@pos}];
	  {clusterList,messages}];

computeApproxMarginalsBPSL[f_,e_:{},emethod_:1,isMax_:False,useSmartMP_:False,printIterations_:False]:=computeMarginalsBPSL[First@clusterGraphCalibrateSL[createClusterGraphSL[f,e,emethod],isMax,useSmartMP,printIterations],isMax];

hyperProbabilityTableSL[output_,singleax_]:=Module[{y,vars,subsets,axR,extractBr,axRel,muF,compositeSet,compositeSetAndNotSubset,
pickOutput,bVague,bSharp,probability},
	If[Total@output!=1,Print["ccfSL error: Hyperdimensional Total is ",Total@output]];
	vars=Array[y,Length@singleax];
	subsets=Complement[Subsets[vars],{{},vars}];
	axR=Apply[Plus,subsets/.Thread[vars->singleax],1];
	extractBr[element_]:=Extract[axR,Position[subsets, element]];
	axRel[i_,j_]:=If[Or[Intersection[i,j]=={},extractBr[j]=={0}],{0},extractBr[Intersection[i,j]]/extractBr[j]];
	muF=axR Last@output;
	compositeSet=Complement[subsets,Thread[{vars}]];
	compositeSetAndNotSubset=Cases[#,Except[{}]]&/@Table[If[Not@SameQ[j,i],j,{}],{i,subsets},{j,compositeSet}];
	pickOutput[subset_]:=Part[output,First@Position[subsets,subset]];
	bVague=Flatten@Table[Total[pickOutput[#]axRel[subsets[[i]],#]&/@compositeSetAndNotSubset[[i]]],{i,Length@compositeSetAndNotSubset}];
	bSharp=Flatten@Total[Table[pickOutput[{#}]&/@i,{i, subsets}],{2}]+PadLeft[Take[Most@output,Length@vars-Length@subsets],Length@subsets];
	probability=bSharp+bVague+muF];

ccfSL[bx_,ux_,singleax_]:=Module[{bconsX,y,vars,subsets,axR,extractBr,axRel,bres,bRres,tup,pickbRres,intersectionAllYEqualsX,
unionAllYEqualXAndintersectionAllYNotEqualEmptySet,unionAllYEqualXAndintersectionAllYDoesEqualEmptySet, tot1, tot2,
tot3a,tot3b,tot3,tot4,bcompX,bconsXH,bcompXH,ccfSLOutput,uxpre,output,i,j,k,outputH,probability,nu},
	bconsX=Min/@Transpose@bx;
	vars=Array[y,Length@First@bx];
	subsets=Complement[Subsets[vars],{{},vars}];
	axR=Apply[Plus,subsets/.Thread[vars->singleax],1];
	extractBr[element_]:=Extract[axR,Position[subsets, element]];
	axRel[i_,j_]:=If[Or[Intersection[i,j]=={},extractBr[j]=={0}],{0},extractBr[Intersection[i,j]]/extractBr[j]];
         bres=#-(Min/@Transpose@bx)&/@bx;
	bRres=PadRight[#,Length@subsets]&/@bres;
	tup=Tuples[subsets,Length@bx];
	intersectionAllYEqualsX=Table[Select[tup,i==Intersection@@#&],{i,subsets}];
	unionAllYEqualXAndintersectionAllYNotEqualEmptySet=Table[Complement[Select[tup,And[i==Union@@#,Not[{}==Intersection@@#]]&],Extract[intersectionAllYEqualsX,First@Position[subsets,i]]],{i,subsets}];
	unionAllYEqualXAndintersectionAllYDoesEqualEmptySet=Table[Select[tup,And[i==Union@@#,{}==Intersection@@#]&],{i,subsets}];
	pickbRres[bel_,subset_]:=Part[bRres,bel,First@Position[subsets,subset]];
	tot1=Total[bRres Table[Times@@Drop[ux,{i}],{i,Length@ux}]];
	tot2=Table[Total[Product[pickbRres[j,#[[j]]]axRel[#[[j]],subsets[[i]]],{j,Length@bx}]&/@intersectionAllYEqualsX[[i]]],{i,Length@intersectionAllYEqualsX}];
	tot3a=Table[Product[pickbRres[j,#[[j]]],{j,Length@bx}]&/@unionAllYEqualXAndintersectionAllYNotEqualEmptySet[[i]],{i,Length@unionAllYEqualXAndintersectionAllYNotEqualEmptySet}];
	tot3b=Table[Product[axRel[#[[j]],subsets[[i]]],{j,Length@bx}]&/@unionAllYEqualXAndintersectionAllYNotEqualEmptySet[[i]],{i,Length@unionAllYEqualXAndintersectionAllYNotEqualEmptySet}];
	tot3=Total[(1-tot3b)tot3a,{2}];
	tot4=Table[Total[Product[pickbRres[j,#[[j]]],{j,Length@bx}]&/@unionAllYEqualXAndintersectionAllYDoesEqualEmptySet[[i]]],{i,Length@unionAllYEqualXAndintersectionAllYDoesEqualEmptySet}];
	bcompXH=Total[Flatten/@{tot1,tot2,tot3,tot4}];
	uxpre=Times@@ux;
	ccfSLOutput[bcons_,bconsux_,bcomp_,bcompux_]:=If[Total@bcomp+bcompux==0,Append[bcons,bconsux],
             nu=(1-Total@bcons-bconsux)/(Total@bcomp+bcompux);
	    Append[bcons+nu bcomp,bconsux+nu bcompux]];
	bconsXH=PadRight[bconsX,Length@subsets];
	outputH=ccfSLOutput[bconsXH,uxpre,bcompXH,uxpre];
	probability=hyperProbabilityTableSL[outputH,singleax];
	output={Take[probability,Length@vars]-singleax Last@outputH,Last@outputH};
	If[Total@Flatten@output!=1,Print["Error ... total not 1: ",Total@output]];
	output];

toDirichletsSL[beliefs_]:=Module[{W = 2,bx,ux,ax, rx, sumrx, \[Alpha]x, axbelief, axrx},
	bx=First/@beliefs;ux=beliefs[[All,2]];ax=Last/@beliefs;
	sumrx=safeDivideSL[W,ux]-W;
	rx=bx(W+sumrx);\[Alpha]x=rx+W ax;
	If[Or@@(Thread[#<0]&/@\[Alpha]x),Print["Warning:some \[Alpha]s less than 0"]];
	axbelief=MapThread[#1(1-#2)&,{ax,ux}];
	axrx=axbelief(W+sumrx);
	{\[Alpha]x,rx,axrx,ux}];

toBeliefsSL[dirichlets_]:=Module[{W = 2, rx,axrx,sumrx,bx, ux, axbelief,ax},
        If[Length@Dimensions@dirichlets<3,
	rx={First@dirichlets}; axrx={Last@dirichlets},
	rx=First/@dirichlets; axrx=Last/@dirichlets];
         sumrx = Total/@rx;
         bx =rx/(W + sumrx);ux =W/(W + sumrx);
         axbelief = axrx/(W + Total/@axrx);
         ax=If[Total@#==0,Table[1,Length@#]/Length@#,#/Total@#]&/@axbelief;
         Transpose@{bx,ux,ax}];

combinePrelimiarySL[beliefs_]:=Module[{\[Alpha]x,rx,axrx,ux,avWeighter,certaintyWeights,ccfSLaxrx,axtable,axbelief,ccfSLsingleax,ccfSLcorebelief,res},
	{\[Alpha]x,rx,axrx,ux} = toDirichletsSL[beliefs];
	avWeighter[weights_]:=Plus@@safeDivideSL[MapThread[#1*#2&,{#,weights}],Plus @@ weights]&;
	certaintyWeights=1-ux;
	ccfSLaxrx=(avWeighter@certaintyWeights)@axrx;
	axtable=Table[1,Length@certaintyWeights]/Length@certaintyWeights;
	axbelief=N@Which[
	Plus@@certaintyWeights==0,{Append[axtable*0,1],axtable},
	Times@@certaintyWeights==1,{Append[axtable,0],axtable},
	True,First@toBeliefsSL[{ccfSLaxrx,ccfSLaxrx}]];
	ccfSLsingleax=First@axbelief+axbelief[[2]] Last@axbelief;
	ccfSLcorebelief=ccfSL[beliefs[[All,1]],beliefs[[All,2]],ccfSLsingleax];
	{##,ccfSLsingleax}&@@ccfSLcorebelief];

combineEvidenceSL[beliefs_]:=Module[{vars,b1,bel,prob},
	vars=Union@(First/@beliefs);
	If[Length@vars>1,Print["Error: combining more than one variable: ",vars];Abort[]];
	b1=beliefs[[All,5;;7]];
	bel=Fold[combinePrelimiarySL[{#1,#2}]&,First@b1,Rest@b1];
	prob=bel[[1]]+bel[[2]]bel[[3]];
	{{First@vars,beliefs[[1,2]],prob,prob,bel[[1]],bel[[2]],bel[[3]]}}];

convertSamiamNetworkSL[file_]:= Module[
        {txt, txt1, txt2, val1, val2, subst1, subst2, subst3, out1,out2},
        txt = Import[file, {"Text", "Words"}];
        txt1 = MapThread[Join[{txt[[#1 + 2]]}, txt[[#1 + 4 ;; #2 - 3] ]]&, {Flatten @ Position[txt, "potential"], Flatten @ Position[txt, "data"]}];
        txt2 = MapThread[txt[[#1 + 2 ;; #2]]&, {Flatten @ Position[txt,"data"], Join[Rest @ Flatten @ Position[txt, "potential"], {-1}]}];
        txt2 = StringJoin @@ Riffle[#, ","]& /@ txt2;
        txt2 = StringDelete[#, {";,},potential", ";,}"}]& /@ txt2;
        txt2 = StringReplace[#, {"(," -> "(", ",)" -> ")", ";,}," -> ""}]& /@ txt2;
        txt2 = StringReplace[#, {"(" -> "{", ")" -> "}"}]& /@ txt2;
        val1 = ToExpression /@ txt2;
        val2 = Flatten @ #& /@ val1;
        subst1 = Partition[Flatten @ #, Last @ Dimensions @ #]& /@ (Array[List, #]& /@ Dimensions /@ val1);
        subst2 = Table[Sort @ Thread[{RotateLeft @ txt1[[i]], #}]& /@subst1[[i]], {i, Length @ subst1}];
        subst3 = Flatten @ MapThread[Thread[#1 -> #2]&, {subst2, val2}];
        out1 = ind2AssSL@#&/@RotateRight/@ (Dimensions/@ val1);
        out2 = Table[Sort @ Thread[{txt1[[i]], #}]& /@ out1[[i]], {i,Length @ out1}] /. subst3;
        Transpose[{txt1, RotateRight /@ (Dimensions /@ val1), out2}]]; 

sortFactorVariablesNonSL[factors_]:=Module[{reOrder},
    reOrder@f_:=Module[{pos1,outVars,outCards,colOrder, \[ScriptCapitalD]3},
        If[ Or[First @ f== {},Length@f[[1]]==1],Return@f];
        pos1=Array[List,Length@f[[1]]];
       {outVars,outCards}=Transpose@SortBy[Transpose@{Extract[f[[1]],pos1],Extract[f[[2]],pos1]},First];
        colOrder=Flatten[Position[outVars,#]&/@f[[1]]];
        \[ScriptCapitalD]3=ind2AssSL[outCards][[All,colOrder]]/.Thread[#&/@ind2AssSL@f[[2]]->f[[3]]];
       {outVars,outCards,\[ScriptCapitalD]3}];
       If[Length@Cases[factors,{{__},{__},{__}}]==0,reOrder@factors, reOrder/@factors]]

sortFactorVariablesSL[factors_]:= If[Length@Cases[factors,{{__},{__},{__},{__},{__},_,{__}}]==0,
        convertFactorToDirichletSL[sortFactorVariablesNonSL@{factors[[1]],factors[[2]],factors[[3]]},factors[[6]],sortFactorVariablesNonSL[factors[[1]],factors[[2]],factors[[7]]][[3]]][[1]], 
        convertFactorToDirichletSL[sortFactorVariablesNonSL@Transpose@{factors[[All,1]],factors[[All,2]],factors[[All,3]]},factors[[All,6]],sortFactorVariablesNonSL[Transpose@{factors[[All,1]],factors[[All,2]],factors[[All,7]]}][[All,3]]]]

renameFactorsSL[factors_,list_]:=Module[{vars1,vars2,vars},
	vars1=Union@Join[First/@list,Last/@list];
	vars2="v"<>ToString@#&/@vars1;
	vars=(First/@factors)/.Thread[vars1->vars2];
	SortBy[Join[Thread[{vars}]/.Join[Thread["v"<>ToString@#&/@(First/@list)->Last/@list],
	Thread["v"<>ToString@#&/@(Last/@list)->First/@list]],
	factors[[All,2;;-1]],2],First@First@#&]];

SetAttributes[beliefPropagationSL,{Protected,ReadProtected,Locked}]
End[];
EndPackage[]
