(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["beliefPropagation`" ]; 
(*\[Copyright] Copyright Stuart Nettleton 2020, 2021*)
(*Disclaimer: This software is provided without any warranty express or implied Stuart Nettleton 2020*)
(*Requires Wolfram Mathematica 12.1.0.0 or higher*)
(*All code files and executable documents are with the license GPL 3.0.For details see http://www.gnu.org/licenses/.*)
(*All documents are with the license Creative Commons Attribution 4.0 International (CC BY 4.0).For details see https://creativecommons.org
licenses/by/4.0/.*)
(*As amended based on original Markov Net schema by Daphne Koller, Stanford University Subject Probabilistic Graphical Models, Copyright (C) Daphne Koller, Stanford University, 2012*)

ind2Ass::usage = "ind2Ass[cards,pos:All] calculates a structural assignment with cardinality cards and optional position"; 

ass2Ind::usage = "ass2Ind[cards,ass] calculates an index with cardinality cards for a structural assignment of a factor"; 

factorProductOrSum::usage = "factorProductOrSum[f1,f2,isMax:False] calculates the product of two factors where
 ... each factor is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms"; 

factorMarginalization::usage = "factorMarginalization[a,v,isMax] removes variables v from a factor a by summing-out ... 
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) which finds the maximum factor value over all possible assignments to the marginalized variables and leaves the output as unnormalised natural logarithms"; 

observeEvidence::usage = "observeEvidence[f,e] modifies a set of factors (f) for observed evidence (e) where
 ... each factor is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... evidence is a list of {variable, evidence value}"; 

computeJointDistribution::usage = "computeJointDistribution[f] calculates the joint distribution of a set of factors (f) by cumulative SumProduct of the fators
 ... each factor is a list containing the list of factor variables, list of variable cardinalities & list of variable values"; 

computeMarginal::usage = "computeMarginal[v,f,e] calculates the marginal over variables (v) in the factors (f) with observed evidence (e)
 ... each factor is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... evidence is a list of {variable, evidence value}"; 

computeInitialPotentials::usage = "computeInitialPotentials[cliques_] calculates the standard form of a clique tree (cliques) where
 ... the clique tree is a list contining the adjacency matrix of the tree and the factors
 ... assigns each factor to a clique
 ... computes factor products over all factors assigned to each clique"; 

pruneTree::usage = "pruneTree[cliques] compacts the clique tree (cliques) by removing redundant clique subsets while maintaining the supersets and the running intersection property"; 

createCliqueTree::usage = "createCliqueTree[f,e,method:'Kruskal',print:False] calculates the clique tree for the variables in f given the evidence e
 ... f is the factor list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value}
 ... reurned cliques have the same structure as factors
 ... the Maximum Spanning Tree Method may be the default 'Kruskal' for undirected graphs or 'Prim'
 ... if print is True then the node-clique and the clique-Maximum Spanning Tree graphs are printed"; 

printFactorAndCliqueGraphs::usage = "printFactorAndCliqueGraphs[f,e,method:'Kruskal'] creates & displays factor & clique graphs
 ... identical to createCliqueTree except only prints the graphs
 ... f is the factor list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value}
 ... returned cliques have the same structure as factors
 ... the Maximum Spanning Tree Method may be the default 'Kruskal' for undirected graphs or 'Prim'"; 

getNextCliques::usage = "getNextCliques[p,messages] finds a pair of cliques {i,j} where clique i is ready to transmit a message to clique j
 ... the clique tree (p) 
 ... current messages (messages)
 ... a clique i is ready to transmit to its neighbor j when i has received messages from all of its neighbors except from j
 ... each message is passed only once
 ... where more than one message is ready to be transmitted the numerically smallest pair is returned"; 

cliqueTreeCalibrate::usage = "cliqueTreeCalibrate[p,isMax] calculates the final potentials for a clique tree (p)
 ... two adjacent cliques i and j are calibrated if they agree on the marginals over their shared variables
... if isMax is False the Sum-Product algorithm is used for message passing in the clique tree and the messages are normalize such that the values in the message sum to 1
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms"; 

computeExactMarginalsBP::usage = "computeExactMarginalsBP[f,e:{},isMax:False,method:'Kruskal',print:False] calculates the exact posterior marginal probability distribution for each variable in the initial factors f given the evidence e
 ... f is the clique list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value} 
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms
 ... the Maximum Spanning Tree Method may be the default 'Kruskal' for undirected graphs or 'Prim'
 ... if print is True then the node-clique and the clique-Maximum Spanning Tree graphs are printed"; 

maxDecoded::usage = "maxDecoded[marginals] calculates the best assignment and alphabetic equivalent for each Max-Product marginal variable"; 

naiveGetNextClusters::usage = "naiveGetNextClusters[p,m] finds a pair of clusters {i,j} where cluster i is ready to transmit the m+1 message to its neighbour cluster j where
 ... the cluster graph (p)
 ... current messages (messages)
 ... the method iterates over the messages (cluster pairs) in increasing order where messages are sorted in ascending ordered by their destination index & ties are broken based on the origin index. If m is 0, [i j] will be the pair of clusters with the lowest j value & (of those pairs over this j) lowest i value as this is the 'first' element in our ordering
 ... where more than one message is ready to be transmitted the numerically smallest pair is returned"; 

getNextClusters::usage = "getNextClusters[p,messages,oldMessages,m,useSmart] finds a pair of nodes {i,j} where node i is ready to transmit the m+1 message to node j
 ... the cluster graph (p)
 ... current messages (messages)
 ... oldMessages[[i,j]] contains the value that Messages(i,j) contained immediately before it was updated to its current value
 ... m is the index of the message to be passed
 ... the method iterates over the messages in increasing order where messages are sorted in ascending ordered by their destination index & ties are broken based on the origin index. If m is 0, [i j] will be the pair of clusters with the lowest j value & (of those pairs over this j) lowest i value as this is the 'first' element in our ordering
 ... where more than one message is ready to be transmitted the numerically smallest pair is returned
 ... if useSmartMP is True then an enhanced apporach is used to obtain the next message cluster for processing"; 

smartGetNextClusters::usage = "smartGetNextClusters[p,messages,oldMessages,m] finds a pair of nodes {i,j} where node i is ready to transmit the m+1 message to clique j
 ... the cluster graph (p)
 ... current messages (messages)
 ... oldMessages[[i,j]] contains the value that Messages(i,j) contained immediately before it was updated to its current value
 ... m is the index of the message to be passed
 ... the method iterates over the messages in increasing order where messages are sorted in ascending ordered by their destination index & ties are broken based on the origin index. If m is 0, [i j] will be the pair of clusters with the lowest j value & (of those pairs over this j) lowest i value as this is the 'first' element in our ordering
 ... where more than one message is ready to be transmitted the numerically smallest pair is returned
 ... other message passing schedules can further improve the convergence of the algorithm"; 

createClusterGraph::usage = "createClusterGraph[f,e] calculates an adjacency matrix and list of clusters for a Bethe cluster graph with nodes representing single variable clusters or pairwise clusters for the variables in f given the evidence e where
 ... a Bethe cluster graph is the simplest way to create a cluster graph from a network as a separate cluster is created for each factor and attached to form the cluster graph
 ... f is the cluster list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value} 
 ... each cluster in the returned cluster list contains a list of factor variables, list of variable cardinalities & list of variable values initialized to the initial potential"; 

checkConvergence::usage = "checkConvergence[new_,old_,threshhold:10^-6] determines whether iterative message passing has converged to a threshhold of 10e-6
 ... new and old are message matrices with respective values of the message from cluster i to cluster j
 ... each message in the matrix is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... default threshhold for convergence is \!\(\*SuperscriptBox[\(10\), \(-6\)]\)"; 

clusterGraphCalibrate::usage = "clusterGraphCalibrate[p,isMax,useSmartMP,threshhold:10^-6] returns the adjacency matrix for a cluster graph (p) and the initial potentials of the clusers using loopy belief propagation
 ... each cluster is a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... if isMax is False the algorithm is SumProduct with ouput of normalised probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalised natural logarithms
 ... if useSmartMP is True then an enhanced approach is used to obtain the next message cluster for processing
 ... default threshhold for convergence is \!\(\*SuperscriptBox[\(10\), \(-6\)]\)"; 

computeApproxMarginalsBP::usage = "computeApproxMarginalsBP[f,e:{},isMax:False,useSmartMP:False,threshhold:10^-6] calculates the approximate posterior marginal probability distribution over each variable in f given the evidence e 
 ... f is the clique list with each being a list containing the list of factor variables, list of variable cardinalities & list of variable values
 ... e is the observed evidence being a list of {variable, evidence value}
 ... if isMax is False the algorithm is SumProduct with ouput of normalized probabilities
 ... if isMax is True the algorithm is MaxProduct (Maximum a Posteriori) with output unnormalized natural logarithms
 ... if useSmartMP is True then an enhanced approach is used to obtain the next message cluster for processing
 ... default threshhold for convergence is \!\(\*SuperscriptBox[\(10\), \(-6\)]\)"; 

convertSamiamNetwork::usage = "convertSamiamNetwork[file] converts a SAMAIM hugin .net file to factors
 ... return a list of SAMIAM node names and list of factors"; 

renameFactors::usage = "renameFactors[factors,list] renames the numerical factors and sorts the factors into canonical order
 ... list is a list of pairs whose node numbers are to be swapped e.g. {{1,2}},{3,4}}"; 

sortFactorVariables::usage = "sortFactorVariables[factors] sorts the numerical factors leaving the left-most the same (which is the node) and the remainder in ascending numerical order"; 

normalizeClique::usage = "normalizeClique[clique,isMax] normalizes the clique or factor to a legal probability distribution
 ... if isMax is False the ouput is normalized probabilities
 ... if isMax is True the input is returned unchanged (i.e. unnormalized)"; 

Begin["`Private`"]; 

normalizeClique[clique_, isMax_] :=Module[{f1,fnNorm},
    f1=If[Length@Cases[clique,{{__},{__},{__}}]==0,{clique},clique];
    fnNorm@f_:= {f[[1]],f[[2]],If[!isMax,(*Flatten[#/Total@#&/@Partition[f[[3]],f[[2,1]]]]*)f[[3]]/Total@f[[3]],f[[3]]]};
      First[fnNorm/@f1 ]]

ind2Ass[cards_,pos_:All]:=Part[Reverse@Apply[List,#]&/@Flatten@Array[0,Reverse@cards],pos]
ass2Ind[cards_,ass_]:=Flatten[Position[ind2Ass@cards,#]&/@Flatten[{ass},Length@Dimensions@{ass}-2]]

factorProductOrSum[f1_, f2_, isMax_:False] :=Module[{f,outVars,outCards,vals},
    f={f1,f2};
    {outVars,outCards}=Transpose@Union@Flatten[Transpose/@f[[All,;;2]],1];
     vals=Map[ind2Ass[outCards][[All,Flatten[Position[outVars,#]&/@#[[1]]]]]/.Thread[ind2Ass@#[[2]]->#[[3]]]&,f];
    {outVars,outCards,If[isMax,Plus,Times]@@@Transpose@vals}];

factorMarginalization[f_, mv_:{}, isMax_:False] :=Module[{pos0,pos1,outVars,outCards,dsOut,colOrder,vals},
    If[Or[First @ f == {},Length@f[[1]]==1],Return@f];
    pos0=Flatten[Position[f[[1]],#]&/@mv,1];
    pos1=Complement[Thread[{Range@Length@f[[1]]}],pos0];
    {outVars,outCards}=Transpose@SortBy[Transpose@{Extract[f[[1]],pos1],Extract[f[[2]],pos1]},First];
    If[outVars == {},Print["Error:Resultant factor has empty scope"];Abort[]];
    colOrder=Flatten[Position[outVars,#]&/@Select[f[[1]],!MemberQ[mv,#]&]];
    vals=ind2Ass[outCards][[All,colOrder]]/.Normal@GroupBy[Thread[Delete[#,pos0]&/@ind2Ass@f[[2]]->f[[3]]],First->Last,If[isMax,Max,Plus]@@#&];
    {outVars,outCards,vals}];

observeEvidence[f_, e_] :=Module[{vars, cards,pos1, pos2,e1,f1},
    If[e=={},Return@f,f1=f];
    vars = First /@ e;
    pos1=Position[First/@f1,#]&/@vars;
    cards = First /@ (Extract[f1[[All, 2]], #]& /@ pos1);
    e1=MapThread[{{#1},{#2},ReplacePart[Table[0,#2],#3->1]}&,{vars,cards,Last/@e}];
    pos1=Position[First/@f1,#]&/@First/@First/@e1;
    pos2=Flatten[MapIndexed[Thread[{#1,First@#2}]&,First/@#&/@pos1],1];
    Scan[(f1[[First@#]]=factorProductOrSum[f1[[First@#]],e1[[Last@#]]])&,pos2];
    f1];

computeJointDistribution[f_] := Which[f == {},Print["Error:empty factor list"];{} ,
     Length @ f == 1,factorProductOrSum[First @ f, Join[Most @ First @ f, {1 + 0 Last @ First @ f}], False],
     True,Fold[factorProductOrSum[#1, #2, False]&, First @ f, Rest @ f]]

computeMarginal[v_, f_, e_:{}] :=Module[{factorsEvidence, joint0, joint, marginal},
      If[Length @ f == 0,Print["Error:empty factor list"];Return[{{}, {}, {}}]];
      factorsEvidence = If[e == {},f ,observeEvidence[f,If[Length @ Dimensions @ e == 1,{e},e]]];
      joint = normalizeClique[computeJointDistribution[factorsEvidence], False];
      marginal = factorMarginalization[joint, Complement[First @ joint,v], False];
      normalizeClique[marginal, False]]

computeInitialPotentials[cliqueSet_,method_:"Kruskal",print_:False]:=Module[{cliqueList, factorList,newCliques,spanningTree,potentials},
       cliqueList = cliqueSet[[1]];
       factorList = cliqueSet[[-1]];  newCliques=Sort@Normal@GroupBy[Flatten[Thread[#[[1]]->#[[2]]]&/@MapIndexed[{Min@@Intersection@@(Position[cliqueList,#][[All,1]]&/@#1),First@#2}&,factorList[[All,1]]]],First->Last];
       spanningTree=maxSpanningTree[cliqueList[[newCliques[[All,1]]]],method,print];
       potentials =Which[Length@#==1,factorList[[First@#]],True,computeJointDistribution@factorList[[#]]]&/@newCliques[[All,2]];
      {spanningTree, potentials}]; 

pruneTree[cliques_] := Module[{cnodes, cedges, neighboursI, toRemove, toKeep, i, j, nk},
    toRemove = {};
    cnodes = First @ cliques;
    cedges = cliques[[2]];
    Do[
      If[MemberQ[toRemove, i],Continue[]];
      neighboursI = Flatten @ Position[cedges[[i]], 1];
      Do[
        j = neighboursI[[neighbour]];
        If[MemberQ[toRemove, j],Continue[]];
        If[Length @ Select[MemberQ[cnodes[[j]], #]& /@ cnodes[[i]],# == True&] == Length @ cnodes[[i]],
        Do[
           If[Length @ Intersection[cnodes[[i]], cnodes[[nk]]] ==Length @ cnodes[[i]],
          cedges[[Complement[neighboursI, {nk}], nk]] = 1;
          cedges[[nk, Complement[neighboursI, {nk}]]] = 1;
          Break[]] ,
       {nk, neighboursI}];
     cedges[[i, All]] = 0;
     cedges[[All, i]] = 0;
     toRemove = Append[toRemove, i]] ,{neighbour, Length @ neighboursI}],
     {i, Length @ cnodes}];
    toKeep = Complement[Range @ Length @ cnodes, toRemove];
    cedges[[toKeep, toKeep]]]

createCliqueTree[f_, e_, method_:"Kruskal", print_:False] :=Module[{nodeEdges, nodeGraph, cliques0, cliques, cliqueNodes,
cliqueEdges,cliqueGraph, lengthSepSets, maxSpanningTree},
    nodeEdges =DeleteDuplicates[UndirectedEdge@@#&/@Sort/@Flatten[Subsets[#,{2}]&/@First /@ f,1]];
    nodeGraph = Graph[Union @@ First /@ f, nodeEdges, VertexLabels-> "Name"];
    cliques = Sort@FindClique[nodeGraph, Infinity, All];
    cliqueNodes=Array[c,Length@cliques];
    If[print,Print[Graph[Join[Union @@ First /@ f, cliqueNodes],Flatten@ MapThread[Thread[#1 \[UndirectedEdge] #2]&, {cliqueNodes, cliques}], VertexLabels-> "Name"]]];
    computeInitialPotentials[{cliques,If[e == {},f ,observeEvidence[f, e]]},method,print]]; 

maxSpanningTree[cliques_,method_:"Kruskal",print_:False]:=Module[{ cliqueNodes,c, cliqueEdges,cliqueGraph, lengthSepSets,spanningTree},
         cliqueNodes=Array[c,Length@cliques];
	cliqueEdges = UndirectedEdge@@#&/@Pick[Subsets[cliqueNodes,{2}],Not@SameQ[#,{}]&/@(Intersection@@#&/@Subsets[cliques,{2}])];
	cliqueGraph = Graph[cliqueNodes, cliqueEdges, VertexLabels ->"Name"];
	lengthSepSets = Length @ (Intersection @@ {Extract[cliques, First @ Position[cliqueNodes, First @ #]],
	    Extract[cliques, First @ Position[cliqueNodes, Last @ #]]})& /@ cliqueEdges;
	spanningTree=FindSpanningTree[cliqueGraph, EdgeWeight -> -1*lengthSepSets, Method ->method];
         If[print,Print[HighlightGraph[cliqueGraph, spanningTree, GraphHighlightStyle-> "Thick "]]];
         pruneTree[{cliques, Normal @ AdjacencyMatrix @ spanningTree}]];

printFactorAndCliqueGraphs[f_, e_, method_:"Kruskal"] :=Module[{},createCliqueTree[f, e, "Kruskal", True];{}]; 

getNextCliques[p_, messages_] := Module[{edges, requiredEvidence, providedEvidence, lackingEvidence,
         availableEvidence, reByLe, aeByLe, nextMsg},
        edges = Last @ p;
        requiredEvidence = Position[edges, 1];
        providedEvidence = Position[Thread[# == {{}, {}, {0}}]& /@ #&/@ messages, False];
        lackingEvidence = Complement[requiredEvidence, providedEvidence];
        availableEvidence = Position[Thread[# == {{}, {}, {0}}]& /@ #& /@ messages, False];
        reByLe = Table[Complement[Select[requiredEvidence, Last @ # ==First @ i&], {Reverse @ i}], {i, lackingEvidence}];
        aeByLe = Table[Complement[Select[availableEvidence, Last @ # ==First @ i&], {Reverse @ i}], {i, lackingEvidence}];
        nextMsg = Pick[lackingEvidence, MapThread[#1 == #2&, {aeByLe,reByLe}]];
        If[nextMsg == {},{0, 0},First @ nextMsg ]]

cliqueTreeCalibrate[p_, isMax_] := Module[{i, j, k, m, cliques, cliqueEdges, fp, requiredEvidence, messages,
         cliqueNodes, cliqueGraph, leafNodes, edges, sum},
        cliqueEdges = First @ p;
        cliques =If[isMax == True, Join[Most /@ Last @ p, N @ {#} /. {Indeterminate -> - \[Infinity]}& /@ (Log /@ Last /@ Last @ p), 2],Last @ p];
        messages = Table[{{}, {}, {0}}, Length @ cliques, Length @ cliques];
        cliqueNodes = Range @ Length @ cliques;
        cliqueGraph = AdjacencyGraph[cliqueNodes, cliqueEdges, VertexLabels -> "Name"];
        leafNodes = Pick[cliqueNodes, # == 1& /@ VertexDegree @ cliqueGraph];
        Scan[(messages[[#, First @ First @ Position[cliqueEdges[[#]],1]]] =normalizeClique[cliques[[#]], isMax])&, leafNodes];
        While[ True,
            {i, j} = getNextCliques[{cliques, cliqueEdges}, messages];
            If[{i, j} == {0, 0},Break[]];
            edges = Complement[Flatten @ Position[cliqueEdges[[i]], 1 ], {j}];
            fp = Fold[factorProductOrSum[#1, #2, isMax]&, cliques[[i]],
            messages[[edges, i]] ];
         messages[[i, j]] = normalizeClique[factorMarginalization[fp, Complement[First @ fp, cliques[[i, 1]]], isMax], isMax]];
         Do[edges = Complement[Flatten @ Position[cliqueEdges[[k]],1], {k}];
         fp = Fold[factorProductOrSum[#1, #2, isMax]&, cliques[[k]],
         messages[[edges, k]]];
         cliques[[k]] = factorMarginalization[fp, Complement[First@ fp, cliques[[k, 1]]], isMax],{k, Length @ messages}];
        {cliqueEdges, cliques}]

computeMarginalsBP[cliques_, isMax_] :=Module[{vars, pos1, pos2},
        vars = Union @ Flatten @ (First /@ cliques);
        pos1 = First /@ First /@ (Position[cliques[[All, 1]], #]& /@vars);
        pos2 = MapThread[Complement,{First /@ cliques[[pos1]], Thread[{Range @ Length @ vars}]}];
        MapThread[normalizeClique[factorMarginalization[#1, #2, isMax], isMax]&, {cliques[[pos1]], pos2}]]

computeExactMarginalsBP[f_, e_:{}, isMax_:False, method_:"Kruskal", print_:False]:=
        computeMarginalsBP[Last @ cliqueTreeCalibrate[createCliqueTree[f,e, method, print], isMax], isMax]; 

maxDecoded[marginals_] :=Module[{res},
        res = Fold[Join[#1, First @ Position[#2, Max @ #2]]&, {}, Last/@ marginals];
        {res, StringJoin[res /. Thread[Range @ 26 -> Alphabet[]]]}]

(*LOOPY BELIEF PROPAGATION*)
naiveGetNextClusters[p_, m_] :=Module[{pos, m1},
        pos = Position[First @ p, 1];
        m1 = Mod[m + 1, Length @ pos];
        If[m1 == 0, m1 = Length @ pos];
        Reverse @ pos[[m1]]]

getNextClusters[p_, messages_, oldMessages_, m_, useSmart_] :=
    If[useSmart == False,naiveGetNextClusters[p, m],smartGetNextClusters[p, messages, oldMessages, m]]

smartGetNextClusters[p_, messages_, oldMessages_, m_] :=Module[{pos, m1, from, to},
        pos = Position[First @ p, 1];
        m1 = Mod[m + 1, Length @ pos] + 1;
        Reverse @ pos[[m1]]]

createClusterGraph[f_, e_] :=Module[{f1, clusterNodes, factorNodes, edges0, edges, ajm},
        f1 = observeEvidence[f, e];
        clusterNodes = Range @ Length @ f1;
        edges = Flatten @ MapThread[Thread[#1 \[UndirectedEdge] #2]&, Transpose @ Select[Transpose @ {clusterNodes, First /@ f1}, Length @ #[[2]] > 1&]];
      {Normal @ AdjacencyMatrix @ Graph[clusterNodes, edges, VertexLabels-> "Name"], f1} ]

checkConvergence[new_, old_, threshhold_:10^-6] :=If[Or @@ Thread[Flatten @ Abs[new - old] < threshhold],True,False]

clusterGraphCalibrate[p_, isMax_, useSmartMP_, threshhold_:10^-6] :=Module[{clusterEdges, clusters, vars, card, messageVars,
    messageCards,messageVals, messages, tic, iter, i, j, lastMessages, prevMessage, edges,fp, sum, pos, clusterList, m},
        clusterEdges = First @ p;
        clusters = Last @ p;
        pos = Position[clusterEdges, 1];
        vars = Union @ Flatten @ (First /@ clusters);
        card = First /@ (Extract[clusters[[All, 2]], #]& /@ Position[First /@ clusters, #]& /@ vars);
        messageVars =Table[If[clusterEdges[[i, j]] == 1,Intersection[clusters[[i, 1]], clusters[[j, 1]]],{}],{i, Length @ clusters},{j, Length @ clusters}];
        messageCards =Table[If[clusterEdges[[i, j]] == 1,card[[messageVars[[i, j]]]],{}],{i, Length @ clusters},{j, Length @ clusters}];
        messageVals = Table[If[clusterEdges[[i, j]] == 1,Table[1, Times @@ messageCards[[i, j]]],{0}],{i, Length @ clusters},{j, Length @ clusters}];
        messages = Table[{messageVars[[i, j]], messageCards[[i, j]],  messageVals[[i, j]]}, {i, Length @ clusters}, {j, Length @ clusters}];
        tic = AbsoluteTime[];
        iter = 0;
        lastMessages = messages;
        While[ True,
            iter++;
            {i, j} = getNextClusters[p, messages, lastMessages, iter,useSmartMP];
            prevMessage = messages[[i, j]];
            edges = Complement[Flatten @ Position[clusterEdges[[i]],1], {j}];
            fp = Fold[factorProductOrSum[#1, #2, isMax]&, clusters[[i]],messages[[edges, i]]];
            messages[[i, j]] = normalizeClique[factorMarginalization[fp, Complement[First @ fp, Intersection[clusters[[i, 1]], clusters[[j,1]]]],isMax], isMax];
            If[useSmartMP == True,lastMessages[[i, j]] = prevMessage];
            If[Mod[iter, Length @ pos] == 0,Print["LBP Messages Passed: ", iter];
                If[checkConvergence[messages, lastMessages, threshhold] == True,Break[]];
                If[iter > 1500,Print["Reached ", iter, " iterations ... aborting" ];
                    Abort[]] ];
            If[useSmartMP == False,lastMessages = messages ]
        ];
        Print["Total number of messages passed: ", iter, " in ", AbsoluteTime[]- tic, " seconds"];
        clusterList = clusters;
        Do[clusterList[[Last @ pos[[m]]]] = factorProductOrSum[clusterList[[Last @ pos[[m]]]], messages[[First @ pos[[m]], Last @ pos[[m]]]], isMax ],{m, Length @ pos}];
        {clusterList, messages}    ]

computeApproxMarginalsBP[f_, e_:{}, isMax_:False, useSmartMP_:True, threshhold_:10^-6] :=
computeMarginalsBP[First @ clusterGraphCalibrate[createClusterGraph[f, e], isMax, useSmartMP, threshhold], isMax]

convertSamiamNetwork[file_] :=Module[{txt, txt1, txt2, val1, val2, subst1, subst2, subst3, out1,out2},
        txt = Import[file, {"Text", "Words"}];
        txt1 = MapThread[Join[{txt[[#1 + 2]]}, txt[[#1 + 4 ;; #2 - 3] ]]&, {Flatten @ Position[txt, "potential"], Flatten @ Position[txt, "data" ]}];
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
        out1 = ind2Ass/@RotateRight/@ (Dimensions/@ val1);
        out2 = Table[Sort @ Thread[{txt1[[i]], #}]& /@ out1[[i]], {i,Length @ out1}] /. subst3;
        Transpose[{txt1, RotateRight /@ (Dimensions /@ val1), out2}]]; 

sortFactorVariables[factors_]:=Module[{reOrder},
    reOrder@f_:=Module[{pos1,outVars,outCards,colOrder, \[ScriptCapitalD]3},
        If[ Or[First @ f== {},Length@f[[1]]==1],Return@f];
        pos1=Array[List,Length@f[[1]]];
       {outVars,outCards}=Transpose@SortBy[Transpose@{Extract[f[[1]],pos1],Extract[f[[2]],pos1]},First];
        colOrder=Flatten[Position[outVars,#]&/@f[[1]]];
        \[ScriptCapitalD]3=ind2Ass[outCards][[All,colOrder]]/.Thread[#&/@ind2Ass@f[[2]]->f[[3]]];
       {outVars,outCards,\[ScriptCapitalD]3}];
    reOrder/@factors];

renameFactors[factors_, list_] := Module[{vars1, vars2, vars},
        vars1 = Union @ Join[First /@ list, Last /@ list];
        vars2 = "v" <> ToString @ #& /@ vars1;
        vars = (First /@ factors) /. Thread[vars1 -> vars2];
        SortBy[Join[Thread[{vars}] /. Join[Thread["v" <> ToString @ #& /@ (First /@ list) -> Last /@ list], Thread["v" <> ToString @ #& /@
             (Last /@ list) -> First /@ list]], factors[[All, 2 ;; -1]], 2], First@ First @ #&]]

SetAttributes[beliefPropagation, {Protected, ReadProtected, Locked}]

End[]; 

EndPackage[]
