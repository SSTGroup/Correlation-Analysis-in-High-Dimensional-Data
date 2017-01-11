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

[1] Y. Song, P. J. Schreier, D. Ramírez, and T. Hasija, “Canonical correlation analysis of high-dimensional data with very small sample support,” Signal Processing, vol. 128, pp. 449–458, 2016.

[[2]](Cross-Validation/) C. Lameiro, and P. J. Schreier, “Cross-validation techniques for determining the number of correlated components between two data sets when the number of samples is very small,” Proc. Asilomar Conf. Signals Syst. Computers, Pacific Grove, CA, USA, November 2016.

[[3]](Bootstrap/) Tanuj Hasija, Yang Song, Peter J. Schreier and David Ramirez,“Bootstrap-based Detection of the Number of Signals Correlated across Multiple Data Sets,” Proc. Asilomar Conf. Signals Syst. Computers, Pacific Grove, CA, USA, November 2016.

[[4]](Sparse-CCA/) C. Lameiro, and P. J. Schreier, “A sparse CCA algorithm with application to model-order selection for small sample support,” Proc. IEEE Int. Conf. Acoustics, Speech and Signal Process., New Orleans, LA, USA, March 2017.


