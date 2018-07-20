%   File: MCCA_ModelSelectionComparison.m
%   Copyright (c) <2018> <University of Paderborn>
%   Permission is hereby granted, free of charge, to any person
%   obtaining a copy of this software and associated documentation
%   files (the "Software"), to deal in the Software without restriction,
%   including without limitation the rights to use, copy, modify and
%   merge the Software, subject to the following conditions:
%
%   1.) The Software is used for non-commercial research and
%       education purposes.
%
%   2.) The above copyright notice and this permission notice shall be
%       included in all copies or substantial portions of the Software.
%
%   3.) Publication, Distribution, Sublicensing, and/or Selling of
%       copies or parts of the Software requires special agreements
%       with the University of Paderborn and is in general not permitted.
%
%   4.) Modifications or contributions to the software must be
%       published under this license. The University of Paderborn
%       is granted the non-exclusive right to publish modifications
%       or contributions in future versions of the Software free of charge.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
%   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%   OTHER DEALINGS IN THE SOFTWARE.
%
%   Persons using the Software are encouraged to notify the
%   Signal and System Theory Group at the University of Paderborn
%   about bugs. Please reference the Software in your publications
%   if it was used for them.
% ------------------------------------------------------------------------
% SYNTAX:   
%
% This is a script, not a function.  Parameters are set within.
%
% OVERVIEW: 
%
% This script should compare the accuracy of the complete model selection
% method [1] with the product of coherence matrices method [3] and IMCCA
% [4]. It should also compute precision and recall for the nonzero
% correlations.                 
%
% DEPENDENCIES:
%
% [1] MutlisetDataGen_CorrMeans_Script.m
%
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% [2] CorrelationStructureGen.m
%
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% [3] MultiSetDataGen_CorrMeans.m
%
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% [4] MCCA_CompleteModelSelection.m
% 
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% [5] ProdOfCoherence_WithBootstrap.m
%
%       Not included within.
%       
% REFERENCES:
%
% Complete model selection method:
% [1]   T. Marrinan, T. Hasija, C. Lameiro, and P. Schreier. "Complete 
%       model selection in multiset canonical correlation analysis." 
%       (In preparation.)
%
% The maxmin method (used in [1]):
% [2]   Song, Yang, Peter J. Schreier, David Ramírez, and Tanuj Hasija. 
%       "Canonical correlation analysis of high-dimensional data with very 
%       small sample support." Signal Processing 128 (2016): 449-458.
%
% Product of coherence matrices method:
% [3]   Hasija, Tanuj, Yang Song, Peter J. Schreier, and David Ramírez. 
%       "Bootstrap-based detection of the number of signals correlated 
%       across multiple data sets." In Signals, Systems and Computers, 
%       2016 50th Asilomar Conference on, pp. 610-614. IEEE, 2016.
%
% IMCCA method:
% [4]   Asendorf, Nicholas, and Raj Rao Nadakuditi. "Improving multiset 
%       canonical correlation analysis in high dimensional sample 
%       deficient settings." In Signals, Systems and Computers, 2015 
%       49th Asilomar Conference on, pp. 112-116. IEEE, 2015. Harvard
%
% ------------------------------------------------------------------------
% CREATED:      21/02/2018 by Tim Marrinan
%
% LAST EDITED:  21/02/2018 by Tim Marrinan
%               
% NOTES: 
% ------------------------------------------------------------------------

clc; clear variables;


% This generates data. Modify this file to change scenario parameters.
run('MultisetDataGen_CorrMeans_Script.m')
n_combs = size(x_corrs,1);
CorrTruth = zeros(n_combs, tot_dims);
CorrTruth(:,1:signum) = p>0;

% Other parameters
ITC = 'ht';                 % Which detector to use for maxmin PCA-CCA.
Pfa = 0.001;                % Probability of false alarm.
boot_iters = 500;           % Number of bootstap resamplings to do.
k_hat_percent   = .90;      % Method (3) requires us to keep some  number 
                            % of singular vectors from each subspace. This 
                            % tells us what percentage of the variance to 
                            % retain in the kept singular vectors.
emp_dist_iters  = 1*1e3;    % Samples to generate for the IMCCA empirical 
                            % distribution

%% Run the complete model selection method.
% ----------------------------------------------------------------
[CorrEst,~] = MCCA_CompleteModelSelection(X,ITC,Pfa);

%% Comparions methods
% ----------------------------------------------------------------
% Find orthonormal bases for data.
fullU = cell(n_sets,1);
origSigma = cell(n_sets,1);
for i = 1 : n_sets
    [fullU{i},origSigma{i},~] = svd(X{i}',0);
    fullU{i} = fullU{i}(:,1:subspace_dims(i));
end

% Run product of coherence matrices method [3].
% ----------------------------------------------------------------
d_PoC = ProdOfCoherence_WithBootstrap(n_sets,X,fullU,M,Pfa,boot_iters);                            
                            
% Run IMCCA method [4].
% ----------------------------------------------------------------
% There is no clear explanation of how to choose \hat{k}_j in the
% paper. They seem to know this value apriori. We will take the 
% singular vectors with 90% of the variance for each data set.

% Generate empirical distribution
emp_evals = zeros(emp_dist_iters,n_sets*tot_dims);
for j = 1 : emp_dist_iters
    N = cell(n_sets,1);
    un = cell(n_sets,1);
    emp_meta = zeros(M,tot_dims*n_sets);
    for i = 1 : n_sets
        N{i} = sqrt(sigmaN)*randn(subspace_dims(i),M); % Noise in the ith channel
        [un{i},~,~] = svd(N{i}',0);
        emp_meta(:,((i-1)*tot_dims)+1:i*tot_dims) = un{i};
    end
    emp_evals(j,:) = max(svd(emp_meta,0)-1 ,0);
end

% Find informative dimensions
k_hat = zeros(n_sets,1);
for i = 1 : n_sets
    k_hat(i) = find(cumsum(diag(origSigma{i}))./...
        sum(diag(origSigma{i}))>k_hat_percent,1);
end
IMCCA_META = zeros(M,sum(k_hat));
for i = 1 : n_sets
    IMCCA_META(:,sum(k_hat(1:i-1))+1:sum(k_hat(1:i))) = ...
        fullU{i}(:,1:k_hat(i));
end
lam = svd(IMCCA_META,0);

% % Define the Marcenko-pastur distribution with lambda = M/tot_dims
% lm = M/tot_dims;
% lm = 2;
% s_sq = sigmaN^2;
% lm_plus = s_sq*(1+sqrt(lm))^2;
% lm_minus = s_sq*(1-sqrt(lm))^2;
% mp =@(x)real( sqrt((lm_plus-x).*(x-lm_minus))./(2*pi*s_sq*lm.*x) );
% x = lm_minus:0.1:lm_plus;
% y = mp(x);

% Generate empirical signal-free distribution
imcca_check = zeros(length(lam),1);
for j = 1 : length(lam)
    [h,edges] = histcounts(emp_evals(:,j),1000,...
        'Normalization','probability');
    thresh = edges(sum(cumsum(h)<(1-Pfa))+1);
    imcca_check(j) = max(lam(j)-1,0)>thresh;
end
d_IMCCA = sum(imcca_check);

%% Performance metrics
% ----------------------------------------------------------------
% Precision and recall.
true_positives = sum(sum((CorrEst+CorrTruth) == 2));
false_positives = sum(sum((CorrEst-CorrTruth) == 1));
false_negatives = sum(sum((CorrEst-CorrTruth) == -1));
true_negatives = sum(sum((CorrEst+CorrTruth) == 0));
precision = true_positives./(true_positives+false_positives);
precision(isnan(precision)) = 0;
recall = true_positives./(true_positives+false_negatives);
recall(isnan(recall)) = 0;

CorrTruth
CorrEst
fprintf('\n\tPrecision: %f\n\tRecall: %f\n',...
    precision, recall);

% Product of coherence matrices comparison
answers = sum(sum(CorrTruth,1)==size(x_corrs,1));
d_comp = sum(sum(CorrEst,1) == n_combs,2);
fprintf('\n\tSignals correlated across all sets: %d\n',answers);
fprintf('\tEstimated by the proposed method estimated: %d\n',d_comp);
fprintf('\tEstimated by product of coherence matrices: %d\n',d_PoC);

    
% IMCCA comparison
answers = sum(sum(CorrTruth,1)>0);
d_tot = sum(sum(CorrEst,1) >0,2);
fprintf('\n\tTotal correlated signals: %d\n',answers);
fprintf('\tEstimated by the proposed method estimated: %d\n',d_tot);
fprintf('\tEstimated by IMCCA: %d\n',d_IMCCA);

