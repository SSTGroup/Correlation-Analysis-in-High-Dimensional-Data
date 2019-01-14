%   File: IMCCA.m
%   Copyright (c) <2018> <University of Paderborn>
%   Signal and System Theory Group, Univ. of Paderborn, http://sst.upb.de
%   https://github.com/SSTGroup/Correlation-Analysis-in-High-Dimensional-Data
%
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
% [d_IMCCA] = IMCCA(X,empEvals,k_hat_percent,Pfa)
%
% OVERVIEW: 
%
% MATLAB implementation of the IMCCA method from [2]. See notes.
%       
% INPUTS:
%
% 'X'               Cell array of size n_sets x 1.
%                   The pth cell contains a matrix of size n_p x M.  It is 
%                   the matrix of observations of the pth data channel.
%
% 'k_hat_percent'   Scalar.
%
% 'empEvals'        
%
% 'Pfa'
%
% 'varargin'        Optional inputs.
%
% 'varargin{1}'     Cell array of the same size as 'X'.
%                   If used, contains orthonormal bases for the data sets.
%
% 'varargin{2}'     Cell array of the same size as 'X'.
%                   If used, contains singular values of the data sets.
%                   
% OUTPUTS:  
%
% 'd_IMCCA'
%                   
%
% DEPENDENCIES:
%
% IMCCA_GenEmpiricalDistribution.m
%       Not included within.
%       Written by Tim Marrinan. See file for documentation.
%       
% REFERENCES:
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
%
% There is no clear explanation of how to choose \hat{k}_j in the paper. 
% They seem to know this value apriori. We will take the singular vectors 
% with 'k_hat_percent'% of the variance for each data set. [TM 21/02/2018]
%
% We set a probability of false alarm for detecting the number of total
% signals by creating an empirical distribution for eigenvalues of matrices
% generated from the same correlation structure. [TM 21/02/2018]
% ------------------------------------------------------------------------
function [d_IMCCA] = IMCCA(X,empEvals,k_hat_percent,Pfa,varargin)

n_sets = size(X,1);
M = size(X{1},2);
if isempty(varargin)
    % Find orthonormal bases for the data sets.
    fullU = cell(n_sets,1);
    origSigma = cell(n_sets,1);
    for i = 1 : n_sets
        [fullU{i},origSigma{i},~] = svd(X{i}',0);
    end
else
    fullU = varargin{1};
    origSigma = varargin{2};
end 

% Find informative dimensions
k_hat = zeros(n_sets,1);
for i = 1 : n_sets
    k_hat(i) = find(cumsum(diag(origSigma{i}))./...
        sum(diag(origSigma{i}))>k_hat_percent,1);
end

% Compute singular values of sample augmented coherence matrix
IMCCA_META = zeros(M,sum(k_hat));
for i = 1 : n_sets
    IMCCA_META(:,sum(k_hat(1:i-1))+1:sum(k_hat(1:i))) = ...
        fullU{i}(:,1:k_hat(i));
end
lam = svd(IMCCA_META,0);

% Generate empirical distribution
imcca_check = zeros(length(lam),1);
for j = 1 : length(lam)
    [h,edges] = histcounts(empEvals(:,j),1000,...
        'Normalization','probability');
    thresh = edges(sum(cumsum(h)<(1-Pfa))+1);
    imcca_check(j) = max(lam(j)-1,0)>thresh;
end
d_IMCCA = sum(imcca_check);