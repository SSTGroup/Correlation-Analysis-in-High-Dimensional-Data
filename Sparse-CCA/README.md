###############
# Description #
###############

Files implementing the SCCA methods described in the paper "A sparse CCA algorithm with application to model-order selection for small sample support" by Christian Lameiro and Peter J. Schreier,
Proceedings of the International Conference on Accoustics, Speech and Signal Processing (ICASSP), New Orleans, LA, USA, 2017

############
# Abstract #
############

We address the problem of determining the number of signals correlated between two high-dimensional data sets with small sample support. In this setting, conventional techniques based on canonical correlation analysis (CCA) cannot be directly applied since the canonical correlations are significantly overestimated when computed from few samples. To overcome this problem, a principal component analysis (PCA) preprocessing step is usually performed to reduce the dimension of the data. However, PCA reduces the dimension of each data set individually without taking the correlation between the data sets into account. In this paper we propose a sparse CCA (SCCA) algorithm as an alternative to the PCA-CCA approach. This algorithm is based on l1-norm penalization, which optimizes the weight of the $\ell_1$-norm to keep a prescribed number of non-zero components. The number of correlated components is then selected based on an information-theoretic criterion.

###########
# Contact #
###########

In case of questions, suggestions, problems etc. please send an email.

Christian Lameiro:
christian.lameiro@sst.upb.de
