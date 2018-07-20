%   File: CorrelationStructureGen.m
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
% [p,sigma_signals,R] = CorrelationStructureGen(n_sets,tot_corr,...
%   corr_std,signum,sigmad,sigmaf,maxIters)
%
% OVERVIEW: This function generates a correlation structure and
%           augmented covariance matrix for multiple data sets. Correlation
%           coefficients for each signal are randomly selected from a
%           normal distribution with a user prescribed mean and standard
%           deviation. This function enforces the 'transitive correlation 
%           condition.'  This condition is not necessary, but makes 
%           bookkeeping easier.
%
% INPUTS:
%
% 'n_sets'          Scalar.
%                   # of data sets to generate.
%
% 'tot_corr'        Vector of size 1 x # of correlated signals.
%                   Elements are the # of datasets correlated with each
%                   signal.
%
% 'corr_means'      Vector of size 1 x # of correlated signals.
%                   The ith element is the mean of the correlation 
%                   coefficients associated with the ith correlated signal.
%
% 'corr_std'        Vector of size 1 x # of correlated signals.
%                   The ith element is the standard deviation of the
%                   correlation coefficients associated with the ith
%                   correlated signal.
%
% 'signum'          Scalar.
%                   # of correlated + independent signals per set.
%
% 'sigmad'          Scalar.
%                   The variance of the correlated components.
%
% 'sigmaf'          Scalar.
%                   The variance of the idependent components.
%
% OUTPUTS:
%
% 'p'               Matrix of size 'n_sets choose two' x signum.
%                   Rows have the same order as x_corrs.  The ith element
%                   of the jth row is the correlation coefficient between
%                   the ith signals in the data sets indicated by the jth
%                   row of x_corrs.
%
% 'sigma_signals'   Matrix of size 'n_sets choose two' x signum.
%                   Rows have the same order as x_corrs.  The ith element
%                   of the jth row is the variance of the ith signals in 
%                   the data sets indicated by the jth row of x_corrs.
%
% 'R'               Matrix of size (n_sets x signum) x (n_sets x signum).
%                   Augmented block correlation matrix of all the data 
%                   sets. Each block is of size signum x signum and the
%                   i-jth block is the correlation matrix between data set
%                   i and data set j.
%
% 'maxIters'        Scalar.
%                   Number of random draws of correlation coefficients
%                   allowed to find a positive definite correlation matrix.
%                   Function throws an error is a PD matrix has not been
%                   found in 'maxIters' attempts.
%
% DEPENDENCIES:
%
% N/A       
%
% REFERENCES:
%
%
% ------------------------------------------------------------------------
% CREATED:      08/01/2018 by Tim Marrinan
%
% LAST EDITED:  17/01/2018 by Tim Marrinan
%               08/01/2018 by Tim Marrinan
%               
%
% NOTES:
%
% See MultisetDataGen_CorrMeans_Script.m for a usage example.
% ------------------------------------------------------------------------
function [p,sigma_signals,R] = CorrelationStructureGen(n_sets,tot_corr,...
    corr_means,corr_std,signum,sigmad,sigmaf,maxIters)

if ~exist('maxIters','var')
    maxIters = 99;
end
corrnum = size(tot_corr,2);
x_corrs = combnk(1:n_sets,2);
n_combi = size(x_corrs,1);

% Check that the number of means and variances matches the number of
% correlated signals requested
if corrnum ~= size(corr_means,2) || corrnum ~= size(corr_std,2)
    error('\n\tA mean correlation and standard deviation need to be\n\tspecified for each correlated signal.\n\n\tThe number of elements in corr_means is:\t%g\n\tThe number of elements in corr_std is:\t\t%g\n\tThe value of totcorr is:\t\t\t%g\n\n\tThe value of totcorr is equal to number of elements in\n\tcorr_across plus the value of fullcorr.\n\n\tNo data has been generated.\n',size(corr_means,2),size(corr_std,2),corrnum)
end

% If the augmented correlation matrix is not positive definite (PD), then 
% the correlation structure is not realizable. Attempt to redraw 
% correlation values a few times if R is not PD, and throw an error if it
% is taking too long.
minEig = -1;
attempts = 0;
while(minEig <= 0)
p = zeros(n_combi,signum);
sigma_signals = zeros(n_combi,signum);
    for j = 1 : corrnum
        % This should implement the method for correlation that requires us
        % to close the loop (transitivity of correlation).
        t = randperm(n_sets,tot_corr(j));
        t = ismember(x_corrs,t);
        t = sum(t,2) == 2;
        temp = corr_means(j) + corr_std(j)*randn(sum(t),1);
        p(t,j) = max(min(temp,1),0);
        sigma_signals(t,j) = sigmad;
        sigma_signals(~t,j) = sigmaf;
    end
   
    if corrnum<signum
        sigma_signals(:,corrnum+1:signum) = ...
            sigmaf.*ones(n_combi,signum-corrnum);
    end
    %% Compute pairwise correlation matrices
    Rxy = cell(nchoosek(n_sets,2),1);
    for i = 1 :size(x_corrs)
        Rxy{i} = sqrt(diag(sigma_signals(i,:))*diag(sigma_signals(i,:)))*diag(p(i,:));
    end

    %% Assemble correlation matrices into augmented block correlation matrix
    for i = 1 : n_sets
        t=sum(x_corrs == i,2);
        temp = max((sigma_signals(find(t),:))==sigmad,[],1);

        R((i-1)*signum+1:i*signum,(i-1)*signum+1:i*signum) = diag(temp*...
            sigmad + ~temp*sigmaf);

        for j = i+1 : n_sets
            a = sum(x_corrs==i,2);
            b = sum(x_corrs==j,2);
            c = find(a.*b);
            R((i-1)*signum+1:i*signum,(j-1)*signum+1:j*signum) = Rxy{c};
            R((j-1)*signum+1:j*signum,(i-1)*signum+1:i*signum) = Rxy{c};    
        end
    end
    attempts = attempts + 1;
    minEig = min(eig(R));
    if attempts > maxIters && minEig <= 0
        error('\n\tA positive definite correlation matrix could\n\tnot be found with the prescribed correlation\n\tstructure after %g attempts. Please adjust\n\tcorr_means and corr_std, and try again.\n\n\tNo data has been generated.\n',attempts)
    end 
end