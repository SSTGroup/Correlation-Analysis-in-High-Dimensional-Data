#################################################
# Correlation-Analysis-in-High-Dimensional-Data #
#################################################

This project is concerned with identifying coupled effects in high-dimensional data. The goal is to extract only a few modes
that explain much of the joint variability between two or more sets of data. This is a common objective, with numerous applications
in many areas of the natural and social sciences and engineering. In this proposal, we deal with high-dimensional problems with 
extremely low sample support, where the coupling must be identified from very few measurements. In such a small sample scenario, 
only very few, dominant, modes are trustworthy, as the remaining modes are due to the spurious effects of noise or sample variability. 
This requires the right trade-off between bias and variance. Too simple a model is a poor representation of the data, causing large bias. 
Too complicated a model overfits the data, causing large variance. We need to strike the right balance between underfitting and 
overfitting. Getting this balance right is the problem of model-order selection. Model-order selection for single data sets, where a 
single model order needs to be determined and sample support is sufficient, is a well-studied problem. There are also techniques that 
work for small sample support and techniques that work for two or more datasets. However, the combination - two or more datasets with 
small sample support - is still a very challenging open problem.


###########
# Contact #
###########

In case of questions, suggestions, problems etc. please send an email.

Tanuj Hasija:
tanuj.hasija@sst.upb.de

Christian Lameiro:
christian.lameiro@sst.upb.de

##############
# References #
##############

[[1]](Techniques-two-data-sets/) Y. Song, P. J. Schreier, D. Ramirez, and T. Hasija, "Canonical correlation analysis of high-dimensional data with very small sample support," Signal Processing, vol. 128, pp. 449-458, 2016.

[[2]](Techniques-two-data-sets/Cross-Validation/) C. Lameiro, and P. J. Schreier, "Cross-validation techniques for determining the number of correlated components between two data sets when the number of samples is very small," Proc. Asilomar Conf. Signals Syst. Computers, Pacific Grove, CA, USA, November 2016.

[[3]](Techniques-multiple-data-sets/Bootstrap/) T. Hasija, Y. Song, P. J. Schreier and D. Ramirez, "Bootstrap-based Detection of the Number of Signals Correlated across Multiple Data Sets," Proc. Asilomar Conf. Signals Syst. Computers, Pacific Grove, CA, USA, November 2016.

[[4]](Techniques-two-data-sets/Sparse-CCA/) C. Lameiro, and P. J. Schreier, "A sparse CCA algorithm with application to model-order selection for small sample support," Proc. IEEE Int. Conf. Acoustics, Speech and Signal Process., New Orleans, LA, USA, March 2017.

[[5]](Techniques-one-data-set/Improper-Signal-Subpsace-Detection/) T. Hasija,  C. Lameiro and P. J. Schreier, "Determining the Dimension of the Improper Signal Subspace in Complex-Valued Data," IEEE Signal Processing Letters, vol. 24, no. 11, pp. 1606-1610, Nov. 2017.

[[6]](Techniques-multiple-data-sets/Complete-Model-Selection/) T. Marrinan, T. Hasija, C. Lameiro and P. J. Schreier,"Complete model selection in multiset canonical correlation analysis," Proc. 26th European Signal Processing Conference (EUSIPCO), Rome, Italy, 2018.

[[7]](Techniques-multiple-data-sets/Complete-Model-Selection-Eigenvalue-Eigenvector-Test-Technique/) T. Hasija, C. Lameiro, T. Marrinan,  and P. J. Schreier,"Determining the Dimension and Structure of the Subspace Correlated Across Multiple Data Sets," Signal Processing, Volume 176, 2020.

[[8]](Techniques-multiple-data-sets/Joint-Reduced-Rank-mCCA-Technique/) T. Hasija and T. Marrinan,"A GLRT for estimating the number of correlated components in sample-poor mCCA," Submitted.

