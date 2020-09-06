# Mathematica-Belief-Propagation
Exact &amp; Loopy Belief Propagation for Undirected (Markov) Networks in Mathematica

Mission statement 

This open source project is for Wolfram Language implementations of decision science algorithms in Subjective Logic, directed Bayesian and undirected Markov networks or graphs that can be used for data analysis, forecast, prediction and recommendation systems.

Organization

The packages are developed & maintained in the latest release of Wolfram Mathematica, which is currently  12.1.0.0. This version or a newer version is required to run the software.

Algorithms implementations are provided in Wolfram package files (“*.wl”). A list of algorithms and usage may be obtained using the Mathematica notebook command:
?beliefPropagationSL`*

Tests and illustrative use are provided in a Mathematica notebook files (“*.nb”), in PDF files, or in Markdown files.

In SubjectiveLogicMarkovNets thre notebooks are provided: beliefPropagationSLtests.nb (unit tests), Gradebook Example.nb and Cancer Example.nb

In Non-Subjective Logic two notebooks are provided: beliefPropagationTests.nb (unit tests) and grade_network.nb (comparative methods for illustration & two procedure for creating factors ... the first procedure is using SAMIAM for automatic creation; and the second procedure is by inspecting Bayesian network directed conditional probabilities). Each notebook requires that the directory path be set to a working directory containing the Wolfram language package with test data in a subdirectiory.

Markov Nets

Markov networks are generally preferred over Bayesian networks for better perspectives on problems where there is a probabilistic interaction between neighbouring variables and in particular when the interaction between variables isn't decisively directed.

Belief Propagation is a message-passing algorithm for performing Bayesian inference in Markov networks. While an exact algorithm (computeMarginals) may be used for small scale marginal distribution inference in undirected networks the computation is NP hard and processing times exponentially increase with the scale of the network. Belief Passing algorithms (computeExactMarginalsBP) are used for improved performance in exact inference in undirected or Markov networks.

Belief propagation is performed with either SumProduct Message Passing (also called Shafer-Shenoy) or MaxProduct Message Passing for Maximum a Posteriori (MAP) estimation. The resulting SumProduct marginal distributions are normalised to sum to 1 while MaxProduct marginal distributions remain unnormalised natural logarithms for accuracy.

A second Belief Passing algorithm called Loopy Belief Propagation (computeApproximateMarginalsBP) can be used for large undirected networks and usually converges although it may not.

The use of undirected Markov networks is not as intuitive as directed Bayesian networks because the factors do not always correspond to probabilities or conditional probabilities. Factors for unirected networks can be readily derived from a corresponding directed Bayesian network either manually or automatically using SAMIAM (both are demonstrated in the grade_network notebook). Factors may be also learned from data using Markov Chain Monte Carlo with Gibbs or Metropolis Hastings sampling and maximum likelihood estimation. This topic is external to this repository as others have contributed Mathematica code (for example StackExchange user Sander "Fast implementation of Metropolis-Hastings update with multiple sets" https://mathematica.stackexchange.com/questions/112677/fast-implementation-of-metropolis-hastings-update-with-multiple-sets ; Burkart, J. "Mathematica Markov Chain Monte Carlo" https://github.com/joshburkart/mathematica-mcmc ; and Gregory, P. "Wolfram Demonstrations Project: Markov Chain Monte Carlo Simulation Using the Metropolis Algorithm" http://demonstrations.wolfram.com/MarkovChainMonteCarloSimulationUsingTheMetropolisAlgorithm/ ).

License matters

These Belief Propagation Algorithms originate with software I auhored as part of Daphne Koller’s Stanford University Subject ‘Probabilistic Graphical Models’, Copyright (c) Daphne Koller, Stanford University, 2012.

All code files and executable documents are with the license GPL 3.0. For details see http://www.gnu.org/licenses/

All documents are with the license Creative Commons Attribution 4.0 International (CC BY 4.0). For details see https://creativecommons.org/licenses/by/4.0/

Disclaimer: This software should be considered as 'alpha stage' & 'not ready for production.' It is provided with no express or implied warranty.
 
Enhancements

The provided algorithms are considered basic and may be progressively enhanced. Suggestions for improvements, extensions and bug-fixes are welcome.

Thank you for experimenting with this software.

Stuart Nettleton 28 July 2020 (initial Non-SubjectiveLogicMarkovNet commit) & 6 September 2020 (initial SubjectiveLogicMarkovNet commit)
