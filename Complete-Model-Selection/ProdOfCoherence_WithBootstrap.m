%   File: ProdOfCoherence_WithBootstrap.m
%   Copyright (c) <2017> <University of Paderborn>
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
%  [dEst] = ProdOfCoherence_WithBootstrap(X,M,Pfa,boot_iters)
%
% INPUTS:
% 'X'           Cell array of size n_sets x 1. 
%               Each cell contains a data matrix.
%               Number of rows in each data matrix is fixed at M and 
%               correpsonds to the number of samples, but the number of 
%               columns can vary between sets.
%
% 'Pfa'         Scalar.
%               Indicates the probability of false alarm for the hypothesis
%               test.
%
% 'boot_iters'  Scalar.
%               Indicates the number of bootstrap iterations to perform.
%
% 'varargin'        Optional inputs.
%
% 'varargin{1}'     Cell array of the same size as 'X'.
%                   If used, contains orthonormal bases for the data sets.
%
% OUTPUTS:
%
% 'dEst'        Scalar.
%               Estimated number of signals correlated across the entire
%               collection of data sets.
%
% DEPENDENCIES:
%               
% [1] bs_mat_org.m 
%
%       Included within.
%       Computes the ordered pairs of indices necessary for all the
%       coherence matrices. See reference for details.
%
% [2] hypothesis_testing_bt.m
%
%       Not included within.
%       Written by Tanuj Hasija, see file for documentation.
%       Computes the estimated model order for correlation across multiple
%       data sets via traditional hypothesis testing of the rank of the
%       product of coherence matrices against the empirical distribution of
%       this rank computed from boot strap samples.
%       
% REFERENCES:
%
% [1]   Hasija, Tanuj, Yang Song, Peter J. Schreier, and David Ramírez. 
%       "Bootstrap-based detection of the number of signals correlated 
%       across multiple data sets." In Signals, Systems and Computers, 
%       2016 50th Asilomar Conference on, pp. 610-614. IEEE, 2016.
%
% ------------------------------------------------------------------------
% CREATED:      
% 
% 19/10/2017 by Tim Marrinan
%
% LAST EDITED:  
% 
% 26/10/2017 by Tim Marrinan
%               
% 19/10/2017 by Tim Marrinan
% 
% NOTES:
% ------------------------------------------------------------------------
function [dEst] = ProdOfCoherence_WithBootstrap(X,Pfa,boot_iters,varargin)

n_sets = size(X,1);
M = size(X{1},2);
if isempty(varargin)
    % Find orthonormal bases for the data sets.
    fullU = cell(n_sets,1);
    for i = 1 : n_sets
        [fullU{i},~,~] = svd(X{i}',0);
    end
else
    fullU = varargin{1};
end 

[pairs] = bs_mat_org_FixedButSlow(n_sets);
[r,~] = size(pairs);
[~,firstmat] = size(fullU{pairs(1,1)});
prodMat = eye(firstmat);
for i = 1 : r
    prodMat = prodMat*(fullU{pairs(i,1)}'*fullU{pairs(i,2)});
end
K = svd(prodMat);
for b = 1 : boot_iters
    X_star = cell(n_sets,1);
    fullU_star = cell(n_sets,1);

    [X_star{1},I] = datasample(X{1},M,2);
    for i = 1 : n_sets
        X_star{i} = X{i}(:,I);
        X_star{i} = X_star{i} - repmat(mean(X_star{i},2),1,M);
        [fullU_star{i},~,~] = svd(X_star{i}','econ');
    end
    [~,firstmat] = size(fullU_star{pairs(1,1)});
    prodMat = eye(firstmat);
    for i = 1 : r
        prodMat = prodMat*(fullU_star{pairs(i,1)}'*fullU_star{pairs(i,2)});
    end
    K_star_matrix(:,b) = svd(prodMat);
end
dEst = hypothesis_testing_bt(Pfa,K,K_star_matrix,boot_iters);
end

function [pairs] = bs_mat_org(L)
if mod(L,2) == 1
    N = (L-1)/2;
    pairs = zeros(L*N,2);
    for n = 1 : N
        q = 1; 
        p = 1;
        j = 0;
        k = L*(n-1)+p;
        while j ~= 1
            if q + n ~= L
                j = mod(q+n,L);
            else
                j = L;
            end
            pairs(k,:) = [q,j];
            q = j;
            p = p + 1;
            k = L*(n-1)+p;
        end
    end
else
    N = L/2;
    L = L + 1;
    pairs = zeros(L*N,2);
    for n = 1 : N
        q = 1; 
        p = 1;
        j = 0;
        k = L*(n-1)+p;
        while j ~= 1
            if q + n ~= L
                j = mod(q+n,L);
            else
                j = L;
            end
            pairs(k,:) = [q,j];
            q = j;
            p = p + 1;
            k = L*(n-1)+p;
        end
     end
    doubles = find(sum(pairs==L,2));
    pairs(doubles(1:2:end),2) = pairs(doubles(2:2:end),2);
    pairs(doubles(2:2:end),:) = [];
end
end

function [pairs] = bs_mat_org_FixedButSlow(L)
    if mod(L,2) == 1
        N = (L-1)/2;
    else
        N = L/2;
    end
    k = 1;
    %pairs = zeros(L*N,2);
    for n = 1 : N
        p = 1;
        for q = 1 : gcd(n,L)
            i = q; 
            j = 0;
            %k = L*(n-1)+p;
            while j ~= q
                if i + n ~= L
                    j = mod(i+n,L);
                else
                    j = L;
                end
                pairs(k,:) = [i,j];
                i = j;
                p = p + 1;
                %k = L*(n-1)+p;
                k = k + 1;
            end
            if q < gcd(n,L)
                pairs(k,:) = [q,q+1];
                p = p + 1;
                %k = L*(n-1)+p;
                k = k + 1;
            elseif q ~= 1
                pairs(k,:) = [q,1];
                p = p + 1;
                %k = L*(n-1)+p;
                k = k + 1;
            end
        end
    end
end
