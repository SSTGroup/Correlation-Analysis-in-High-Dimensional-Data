###############
# Description #
###############

This matlab package contains files implementing the method for model selection and correlation structure estimation in multiple data sets described in the paper "Determining the Dimension and Structure of the Subspace Correlated Across Multiple Data Sets," by Tanuj Hasija, Christian Lameiro, Tim Marrinan and Peter J. Schreier, Submitted.

This matlab package includes the following files:
01. Experiment12.m
02. Experiment3.m
03. Experiment4.m
04. CorrelationStructureGen.m
05. MultisetDataGen_CorrMeans.m
06. Eval_Evec_Tests_Bootstrap_Multiple_Datasets.m
06. MCCA_CompleteModelSelection.m
07. PCM_Bootstrap_Multiple_Datasets.m
08. Max_Min_Multiple_Datasets.m
09. ITC_MDL_Wu_Multiple_Datasets.m
10. MCCA_KPD.m
11. mcca.m
12. data_edit_index.m
13. README.md

############
# Abstract #
############

Detecting the components common or correlated across multiple data sets is challenging due to a large number of possible correlation structures among the components. Even more challenging is to determine the precise structure of these correlations. Traditional work has focused on determining only the model order, i.e., the dimension of the correlated subspace, a number that depends on how the model-order problem is defined. Moreover, identifying the model order is often not enough to understand the relationship among the components in different data sets. We aim at solving the complete model selection problem, i.e., determining which components are correlated across which data sets. We prove that the eigenvalues and eigenvectors of the normalized covariance matrix of the composite data vector, under certain conditions, completely characterize the underlying correlation structure. We use these results to solve the model-selection problem by employing bootstrap-based hypothesis testing.

##############
# File Usage #
##############
Experiment12.m - Runs experiments 1 and 2 in the paper. The parameters of the two experiments are set automatically in 'scen1' and 'scen2'. To run experiments with different settings modify the parameters of the 'custom' scenario as desired.  
Experiment3.m - Runs experiment 3 in the paper. The parameters of the experiment are set automatically. 
Experiment4.m - Runs experiment 4 in the paper. The parameters of the experiment are set automatically in 'scen1'. To run experiments with different settings modify the parameters of the 'custom' scenario as desired.    

In all the above files, the variables 'num_iter', 'verbose', etc can be modified to control for how many iterations, and what kind of information is printed while the experiments are run. Moreover, in 'custom' scenario, a correlation structure can be automatically generated for correlations across a specified number of data sets.

###########
# Contact #
###########

In case of questions, suggestions, problems etc. please send an email.

Tanuj Hasija:
tanuj.hasija@sst.upb.de

This matlab package is hosted at:
https://github.com/SSTGroup/Correlation-Analysis-in-High-Dimensional-Data
