%   File: MultisetDataGen_CorrMeans.m
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
% [X,R,A,S] = MultisetDataGen(subspace_dims,signum,x_corrs,mixing,...
%               sigmaN,color,n_sets,p,sigma_signals,M,MAcoeff,ARcoeff)
%
% OVERVIEW: This script should provide an example of how to generate
%           multiple data sets with a prescribed correlation structure.
%           It enforces the transitive correlation condition.  This
%           condition is not necessary, but makes bookkeeping easier.
%
% INPUTS:
%
% 'subspace_dims'   Vector of size 1 x n_sets;
%                   The ith element is the subspace dimension of the ith
%                   data set.
%
% 'signum'          Scalar.
%                   # of correlated + independent signals per set.
%
% 'x_corrs'         Matrix of size 'n_sets choose two' x 2.
%                   Each row is one of all the possible pairs of indices
%                   for pairwise comparisons of n_sets.
%
% 'mixing'          String.
%                   Either 'orth' or 'randn'. Describes the type of mixing
%                   matrix.
%
% 'sigmad'          Scalar.
%                   The variance of the correlated components.
%
% 'sigmaf'          Scalar.
%                   The variance of the idependent components.
%
% 'sigmaN'          Scalar.
%                   The variance of the noise components.
%
% 'color'           String.
%                   Only 'white' is coded at the moment. To implement 
%
% 'n_sets'          Scalar.
%                   # of data sets to generate.
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
% 'M'               Scalar.
%                   # of samples per data set.
%
% 'MAcoeff'         Vector of size 'degree of MA dependency' x 1.
%                   Moving average coefficients for colored noise.
%                   Se 'help filter' for more details.
%
% 'ARcoeff'         Vector of size 'degree of AR dependency' x 1.
%                   Auto-regressive coefficients for colored noise.
%                   Se 'help filter' for more details.
%
% OUTPUTS:
%
% 'X'               Cell array of size n_sets x 1.
%                   The ith cell contains a matrix of size 'ith element of
%                   subspace_dims' by M.  It is the matrix of observations
%                   of the ith data set plus noise.
%
% 'R'               Matrix of size (n_sets x signum) x (n_sets x signum).
%                   Augmented block correlation matrix of all the data 
%                   sets. Each block is of size signum x signum and the
%                   i-jth block is the correlation matrix between data set
%                   i and data set j.
%
% 'A'               Cell array of size n_sets x 1.
%                   The ith cell contains a matrix of size 'ith element of
%                   subspace_dims' by signum.  It is the mixing matrix
%                   associated with the ith data set.
%
% 'S'               Cell array of size n_sets x 1.
%                   The ith cell contains a matrix of size 'ith element of
%                   subspace_dims' by M.  It is the matrix of noise-free
%                   observations of the ith data set.
%
% DEPENDENCIES:
%
% N/A       
%
% REFERENCES:
%
% Adapted from 'TwoChannelModel.m' written by Christian Lamiero 2016.  Not
% included within.
%
% ------------------------------------------------------------------------
% CREATED:      22/08/2017 by Tim Marrinan
%
% LAST EDITED:  17/01/2018 by Tim Marrinan
%               16/01/2018 by Tim Marrinan
%               15/12/2017 by Tim Marrinan
%               
%
% NOTES:
%
% See MultisetDataGen_CorrMeans_Script.m for a usage example.
%
% Does not generate complex data.  That change shouldn't be hard to make.
%
% ------------------------------------------------------------------------
function [X,R,A,S] = MultisetDataGen_CorrMeans(subspace_dims,signum,...
    x_corrs,mixing,sigmad,sigmaf,sigmaN,color,n_sets,p,sigma_signals,...
    M,MAcoeff,ARcoeff)

%% Generate mixing matrices
A = cell(n_sets,1);
switch lower(mixing)
    case 'orth'
        for i = 1 : n_sets
            A{i}=orth(randn(subspace_dims(i),signum));
        end
    case 'randn'
        for i = 1 : n_sets
            A{i}=randn(subspace_dims(i),signum);
        end
    otherwise
        error('Unknown mixing matrix property');
end

%% Compute pairwise correlation matrices
Rxy = cell(nchoosek(n_sets,2),1);
for i = 1 :size(x_corrs)
    Rxy{i} = sqrt(diag(sigma_signals(i,:))*diag(sigma_signals(i,:)))*diag(p(i,:));
end

%% Assemble correlation matrices into augmented block correlation matrix
for i = 1 : n_sets
    t=sum(x_corrs == i,2);
    temp = max((sigma_signals(t>0,:))==sigmad,[],1);
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

%% Generate M realizations of each data set along with Gaussian noise.
S = cell(n_sets,1);
N = cell(n_sets,1);
fullS=sqrtm(R)*randn(n_sets*signum,M);
for i = 1 : n_sets
    S{i}  = fullS((i-1)*signum+1:i*signum,:); % Realizations of S{i}
    N{i} = sqrt(sigmaN)*randn(subspace_dims(i),M); % Noise in the ith channel
end

%% Filter the noise to be colored if specified
switch lower(color)
    case 'white'
        % nothing needs to happen
    case 'colored'
        for i = 1 : n_sets
            N{i} = filter(MAcoeff,ARcoeff,N{i});
        end
    otherwise
        error('Unknown noise color option');
end

%% Compute final noisy observations
X = cell(n_sets,1);
for i = 1 : n_sets
    X{i}=A{i}*S{i}+N{i}; % ith channel observations
end