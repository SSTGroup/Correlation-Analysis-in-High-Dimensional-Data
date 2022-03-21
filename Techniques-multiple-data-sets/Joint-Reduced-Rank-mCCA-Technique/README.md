###############
# Description #
###############

This matlab package contains files implementing the method for estimating the number of correlated components in multiple data sets described in the paper "A GLRT for estimating the number of correlated components in sample-poor mCCA," by Tanuj Hasija and Timothy Marrinan, Submitted.

This matlab package includes the following files:
01. Experiment.m
02. CorrelationStructureGen.m
03. MultisetDataGen_CorrMeans.m
04. Eval_Evec_Tests_Bootstrap_Multiple_Datasets.m
05. Num_PCA_components
06. MCCA_CompleteModelSelection.m
07. mcca.m
08. maxmin.m
09. jRRmCCA.m
10. Threshold_sam_poor_n_sets10M200rmax15pfa0.01MCtrials500.mat
11. README.md

############
# Abstract #
############

In many applications, components correlated across multiple data sets represent meaningful patterns and commonalities. Estimates of these patterns can be improved when the number correlated components is known, but since data exploration often occurs in an unsupervised setting, the number of correlated components is generally not known. In this paper, we derive a generalized likelihood ratio test (GLRT) for estimating the number of components correlated across multiple data sets. In particular, we are concerned with the scenario where the number of available samples is small. As a result of the small sample support, correlation coefficients and other summary statistics are significantly overestimated by traditional methods. The proposed test combines linear dimensionality reduction with a GLRT based on a measure of multiset correlation referred as the generalized variance cost function (mCCA-GENVAR). By jointly estimating the rank of the dimensionality reduction and the
number of correlated components, we are able to provide high-accuracy estimates in the challenging sample-poor setting. These advantages are illustrated in numerical experiments that compare and contrast the proposed method with existing techniques.

##############
# File Usage #
##############

Experiment.m - Runs experiment in the paper. The parameters are set automatically in 'scen1'. To run experiments with different settings and correlation structure modify the parameters of the 'custom' scenario as desired. The variable 'num_iter' can be modified to control for how many iterations the experiments are run. 

###########
# Contact #
###########

In case of questions, suggestions, problems etc. please send an email.

Tanuj Hasija:
tanuj.hasija@sst.upb.de

This matlab package is hosted at:
https://github.com/SSTGroup/Correlation-Analysis-in-High-Dimensional-Data
