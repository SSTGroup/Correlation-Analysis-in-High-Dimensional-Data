###############
# Description #
###############

This matlab package contains files implementing the method for model selection for multiset canonical correlation analysis described in the paper "Complete model selection in multiset canonical correlation analysis" by Tim Marrinan, Tanuj Hasija, Christian Lameiro and Peter J. Schreier, submitted to Proceedings of the 26th European Signal Processing Conference (EUSIPCO), Rome, Italy, 2018.

This matlab package includes the following files:
01. CorrelationStructureGen.m
02. EUSIPCO_ScenarioSpecification.m
03. hypothesis_testing_bt.m
04. IMCCA_GenEmpiricalDistribution.m
05. IMCCA.m
06. main.m
07. MCCA_CompleteModelSelection.m
08. MCCA_ModelSelectionComparison.m
09. MultisetDataGen_CorrMeans_Script.m
10. MultisetDataGen_CorrMeans.m
11. ProdOfCoherence_WithBootstrap.m
12. readme.txt

############
# Abstract #
############

Traditional model-order selection for canonical correlation analysis infers latent correlations between two sets of noisy data. In this scenario it is enough to count the number of correlated signals, and thus the model order is a scalar. When the problem is generalized to a collection of three or more data sets, signals can demonstrate correlation between all sets or some subset, and one number cannot completely describe the correlation structure. We present a method for estimating multiset correlation structure that combines source extraction in the style of joint blind source separation with pairwise model order selection.  The result is a general technique that describes the complete correlation structure of the collection.

##############
# File Usage #
##############
In the file main.m, the variables 'num_iterations', 'printout', 'plots', and 'queued' can be modified to control which experiments are run for how many iterations, and what kind of information is returned.

To re-run the experiments and generate the plots from the paper:
01. Run main.m
On a laptop with a 2.6 GHz Intel Core i7 processor and 16 GB of RAM, the experiments from the paper run for num_iterations = 1000 takes approximately 40 minutes. With num_iterations = 100, they take approximately 10 minutes. 


To run your own experiments with the model selection method:
01. In the file EUSIPCO_ScenarioSpecification.m, modify the parameters of the 'custom' scenario as desired.
(A correlation structure can be randomly selected for correlations across a specified number of data sets, or the correlation structure can be input explicitly by the user.  See the file for details.)
02. In the file main.m, set the value of the 'queued' variable to 3.
03. In the file main.m, modify the variables 'num_iterations', 'printout', and 'plots' as desired.
04. Run main.m


###########
# Contact #
###########

In case of questions, suggestions, problems etc. please send an email.

Tim Marrinan:
tim.marrinan@sst.upb.de

This matlab package is hosted at:
https://github.com/SSTGroup/Correlation-Analysis-in-High-Dimensional-Data
