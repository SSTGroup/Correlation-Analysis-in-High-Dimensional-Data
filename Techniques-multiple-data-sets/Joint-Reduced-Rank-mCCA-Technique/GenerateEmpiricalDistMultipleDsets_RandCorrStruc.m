% ## ----------------------------------------------------------------------------
% ##
% ##   File: GenerateEmpiricalDistMultipleDsets_RandCorrStruc.m
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

function [tau,f,x,Eta_s_vec] = GenerateEmpiricalDistMultipleDsets_RandCorrStruc(n_sets,tot_dims,M,corr_sig,pfa,emp_dist_iters)

%% signals correlated randomly across the data sets
subspace_dims = repmat(tot_dims,[1,n_sets]);
signum = tot_dims;
full_corr = 0;
x_corrs = combnk(1:n_sets,2);
corr_means = ones(1,corr_sig).*0.99;
corr_std = zeros(1,corr_sig);
sigmad = 1; sigmaf=1;maxIters=99;
sigmaN = 1; mixing = 'orth';
color           = 'white';
MAcoeff         = 1;
ARcoeff         = 1;

for j = 1 : emp_dist_iters
    corr_across = sort(randi([2,n_sets],1,corr_sig),'descend');
    tot_corr = [repmat(n_sets,[1,full_corr]) corr_across];
    
    [p,sigma_signals,~] = CorrelationStructureGen(n_sets,tot_corr,...
        corr_means,corr_std,signum,sigmad,sigmaf,maxIters);

    [X_cell,~,~,~] = MultisetDataGen_CorrMeans(subspace_dims,signum,x_corrs,...
        mixing,sigmad,sigmaf,sigmaN,color,n_sets,p,sigma_signals,M,...
        MAcoeff,ARcoeff);
    
    % Using Genvar canonical variables
    for l=1:n_sets
        X_mcca(:,:,l)  = X_cell{l};
    end
    [W,~] = mcca(X_mcca,tot_dims,'genvar');
    S_est = cell(n_sets,1);
    for l=1:n_sets
        S_est{l} = W(:,:,l)*X_mcca(:,:,l);
    end
    Scv = cell(tot_dims,1);
    for k=1:tot_dims
        for l=1:n_sets
            Scv{k}(l,:) = S_est{l}(k,:);
        end
    end
    
    s = corr_sig;
    
    Eta_s = 0;
    for k=s+1:tot_dims
        
        E_est = real(eig(Scv{k}*Scv{k}.'/M));
        sum(E_est);
        Eta_s= Eta_s + -M*log(prod(E_est));
    end
    Eta_s_vec(j) = Eta_s;
    
end
[f,x] = ecdf(Eta_s_vec);
[a,~] = find(f > (1-pfa));
tau = x(a(1));
end