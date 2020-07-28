# Mathematica-Belief-Propagation
Exact &amp; Loopy Belief Propagation in Mathematica

Mission statement 

This open source project is for Wolfram Language implementations of decision science algorithms in Subjective Logic & Bayesian networks that can be used for data analysis, forecast, prediction and recommendation systems.

An exact algorithm (computeMarginals) may be used for small scale marginal distribution inference. However this computation is NP hard and processing times exponentially increase with the scale of the network.

Belief Propagation is a message-passing algorithm for performing Bayesian inference. A Belief Passing algorithm (computeExactMarginalsBP) may be used for improved performance in exact inference.

Belief propagation is performed with either SumProduct Message Passing (also called Shafer-Shenoy) or MaxProduct Message Passing for Maximum a Posteriori (MAP) estimation. The resulting SumProduct marginal distributions are normalised to sum to 1. However MaxProduct marginal distributions remain unnormalised natural logarithms.

A second Belief Passing algorithm called Loopy Belief Propagation (computeApproximateMarginalsBP) can be for large networks and usually converges although it may not.

License matters

The Belief Propagation Algorithms derive from Octave software I developed for Daphne Koller’s Stanford University Subject ‘Probabilistic Graphical Models’, Copyright (c) Daphne Koller, Stanford University, 2012.

All code files and executable documents are with the license GPL 3.0. For details see http://www.gnu.org/licenses/

All documents are with the license Creative Commons Attribution 4.0 International (CC BY 4.0). For details see https://creativecommons.org/licenses/by/4.0/

Disclaimer: This software should be considered as 'alpha stage' & 'not ready for production.' It is provided with no express or implied warranty
 
Organization

The packages are developed & maintained in the latest release of Wolfram Mathematica, which is currently  12.1.0.0. This version or a newer version is required to run the software.

Algorithms implementations are provided in Wolfram package files (“*.wl”). A list of algorithms and usage may be obtained using the Mathematica notebook command:
?beliefPropagation`*

Tests and illustrative use are provided in a Mathematica notebook files (“*.nb”), in PDF files, or in Markdown files. Set the directory path in the Mathematica notebook to your own directory.

Enhancements

The provided algorithms are considered basic and may be progressively enhanced. Suggestions for improvements, extensions and bug-fixes are welcome.

Thank you for experimenting with this software.

Stuart Nettleton 28 July 2020 (initial commit)
