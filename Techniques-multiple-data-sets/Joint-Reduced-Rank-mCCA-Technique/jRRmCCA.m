% ## ----------------------------------------------------------------------------
% ##
% ##   File: jRRmCCA.m
% ##   Copyright (c) <2022> <University of Paderborn>
% ##   Permission is hereby granted, free of charge, to any person
% ##   obtaining a copy of this software and associated documentation
% ##   files (the "Software"), to deal in the Software without restriction,
% ##   including without limitation the rights to use, copy, modify and
% ##   merge the Software, subject to the following conditions:
% ##
% ##   1.) The Software is used for non-commercial research and
% ##       education purposes.
% ##
% ##   2.) The above copyright notice and this permission notice shall be
% ##       included in all copies or substantial portions of the Software.
% ##
% ##   3.) Publication, Distribution, Sublicensing, and/or Selling of
% ##       copies or parts of the Software requires special agreements
% ##       with the University of Paderborn and is in general not permitted.
% ##
% ##   4.) Modifications or contributions to the software must be
% ##       published under this license. The University of Paderborn
% ##       is granted the non-exclusive right to publish modifications
% ##       or contributions in future versions of the Software free of charge.
% ##
% ##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% ##   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% ##   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% ##   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% ##   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% ##   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% ##   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% ##   OTHER DEALINGS IN THE SOFTWARE.
% ##
% ##   Persons using the Software are encouraged to notify the
% ##   Signal and System Theory Group at the University of Paderborn
% ##   about bugs. Please reference the Software in your publications
% ##   if it was used for them.
% ##
% ##
% ##   Author: Tanuj Hasija
% ##   Created: Tanuj Hasija 01/02/2022
% ##   Edited: Tanuj Hasija 01/02/2022
% ##   Dependencies:
% ##
% ##   [1] mcca.m
%
%       Not included within.
%       Written by Yiou Li (UMBC), see file for documentation.
% ## ----------------------------------------------------------------------------

function [d_est,r_est] = jRRmCCA(X,rmax,tau)

%%%%%%%%%% Input Parametes %%%%%%%%%%%
% X_cell: Data sets in a cell format, each data set is an element of the cell

% rmax: Maximum dimensionality reduction rank

% tau: Matrix of thresholds computed for different reduced ranks and signal
% index, See Experiment script and
% GenerateEmpiricalDistMultipleDsets_RandCorrStruc function for more details

n_sets = size(X,1); % Number of data sets
M = size(X{1},2); % Numer of samples

for r=1:rmax
    X_pca = cell(n_sets,1);
    X_mcca = [];
    % PCA and Using Genvar canonical variables
    for l=1:n_sets
        [U,~,~] = svd(X{l},'econ');
        X_pca{l} = U(:,1:r).'*X{l};
        X_mcca(:,:,l)  = X_pca{l};
    end
    [W,~] = mcca(X_mcca,r,'genvar');
    S_est = cell(n_sets,1);
    for l=1:n_sets
        S_est{l} = W(:,:,l)*X_mcca(:,:,l);
    end
    Scv = cell(r,1);
    for j=1:r
        for l=1:n_sets
            Scv{j}(l,:) = S_est{l}(j,:);
        end
    end
    
    for s=0:r-1
        Eta_s = 0;
        for j=s+1:r
            
            E_est = real(eig(Scv{j}*Scv{j}.'/M));
            sum(E_est);
            Eta_s= Eta_s + -M*log(prod(E_est));
        end
        
        idx = s+1;
        if(Eta_s < tau(r,idx))
            d_hat_r = s;
            break;
        elseif(s== r-1)
            d_hat_r = r;
        end
    end
    d_vec_r(r) = d_hat_r;
end
[d_est, r_est] = max(d_vec_r);
end