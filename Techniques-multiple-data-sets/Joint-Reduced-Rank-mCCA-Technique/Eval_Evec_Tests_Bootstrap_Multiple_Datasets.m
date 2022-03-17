% ## ----------------------------------------------------------------------------
% ##
% ##   File: Eval_Evec_Tests_Bootstrap_Multiple_Datasets.m
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
% ##   Created: Tanuj Hasija 07/01/2018
% ##   Edited: Tanuj Hasija 07/01/2018
% ##   Dependencies:
% ##
% ## ----------------------------------------------------------------------------

function [d_cap,U_struc] = Eval_Evec_Tests_Bootstrap_Multiple_Datasets(X_cell, P_fa_eval, P_fa_evec, B)


%%%%%%%%%% Input Parametes %%%%%%%%%%%
% X_cell: Data sets in a cell format, each data set is an element of the cell
% P_fa: Probability of false alarm for hypothesis testing for eigenvalue
% test (P_fa_eval) and eigenvector tests (P_fa_evec)
% B: No. of bootstrap iterations

% display('started');

L = size(X_cell,1); % Number of data sets
M = size(X_cell{1},2); % Numer of samples

% U_cell{L} = []; % Contains the left singular vector matrices of L data sets
% V_cell{L} = []; % Contains the right singular vector matrices of L data sets
tot_dim = 0;
N = size(X_cell{1},1);
X_aug = zeros(N*L,M);
Rxx_mH{L} = [];

for l=1:L
    %     [U_temp,~,V_temp] = svd(X_cell{l},'econ');
    %     U_cell{l} = U_temp;
    %     V_cell{l} = V_temp;
    
%     X_aug = [X_aug;X_cell{l}];
    X_aug((l-1)*N+1:l*N,:) = X_cell{l}; 
    Rxx_mH{l} = sqrtm(inv(X_cell{l}*X_cell{l}.'/M));
    
    m{l} = size(X_cell{l},1); % dimension of each data set
    tot_dim = tot_dim + m{l};
    
end

Rxx_aug = X_aug*X_aug.'/M;

% Generating augmented matrix
Rd_mH = Rxx_mH{1};
aug_dim{L} = []; aug_dim{1} = m{1}; % Augumented dimension required for eigenvector test

for l=2:L
    Rd_mH = blkdiag(Rd_mH,Rxx_mH{l});
    aug_dim{l} = aug_dim{l-1} + m{l};
end

Cxx_aug = Rd_mH*Rxx_aug*Rd_mH;
[U, E] = eig(Cxx_aug);
E = diag(E);
E = abs(E); % Sometimes for very few samples, E might become complex
[E,I] = sort(E,'descend');
U = U(:,I);

X_cell_star{L} = [];
% Bootstrap Operation
for b=1:B
    
    [~,I] = datasample(X_cell{1},M,2); % Bootstrap data matrix
    
%     X_aug_star = []; %
    X_aug_star = zeros(N*L,M);
%     Rxx_mH_star{L} = [];
    Rxx_mH_star = cell(L,1);
    for l=1:L
        X_cell_star{l} = X_cell{l}(:,I);
        X_cell_star{l} = X_cell_star{l} - repmat(mean(X_cell_star{l},2),1,M);
%         X_aug_star = [X_aug_star;X_cell_star{l}];
         X_aug_star((l-1)*N+1:l*N,:) = X_cell_star{l}; 
        Rxx_mH_star{l} = sqrtm(inv(X_cell_star{l}*X_cell_star{l}.'/M));
    end
    Rd_mH_star = Rxx_mH_star{1};
    for l=2:L
        Rd_mH_star = blkdiag(Rd_mH_star,Rxx_mH_star{l});
        
    end
    Rxx_aug_star = X_aug_star*X_aug_star.'/M;
    Cxx_aug_star = Rd_mH_star*Rxx_aug_star*Rd_mH_star;
    [U_star, E_star] = eig(Cxx_aug_star);
    E_star = diag(E_star);
    E_star = abs(E_star);
    [E_star,I] = sort(E_star,'descend');
    U_star = U_star(:,I);
    
    
    E_star_matrix(:,b) = E_star;
    U_star_matrix(:,:,b) = U_star;
    
end

d_cap = hypothesis_testing_bt_eval_aug_coh_matrix(L,tot_dim, P_fa_eval,E,E_star_matrix,B);
U_struc = hypothesis_testing_bt_evec_aug_coh_matrix(L,aug_dim, P_fa_evec,d_cap,U,U_star_matrix,B);
end


function k_cap = hypothesis_testing_bt_eval_aug_coh_matrix(L, tot_dim, P_fa_eval,E,E_star_matrix,B)
% Using the fact that atleast L (no. of datasets) eigenvalues will be equal to 1
% Testing statistic T2 = sum of squared elements = 0;
% L (eigenvalues-1) will be equal to 0.
Lambda = E-1;
Lambda_star_matrix = E_star_matrix-1;

for k=0:sum(tot_dim)-L
    
    T2 = sum(Lambda(k+1:k+L).^2);
    
    for b=1:B
        T2_star(b) =  sum(Lambda_star_matrix(k+1:k+L,b).^2); % calc_statistic(statistic_type,L_star_matrix(:,b),k);
        T2_null(b) = T2_star(b) - T2;
        if(abs(T2) <= abs(T2_null(b)))
            Indicator(b) = 1;
        else
            Indicator(b) = 0;
        end
    end
    p_value = sum(Indicator)/B;
    if(p_value >= P_fa_eval)
        k_cap = k;
        break;
    end
end
end

function U_struc = hypothesis_testing_bt_evec_aug_coh_matrix(L,aug_dim, P_fa_evec,d_cap,U,U_star_matrix,B)

U_struc = zeros(d_cap,L);

for i=1:d_cap
    
    for l=1:L
        present = 1;
        if l==1
            dim1 = 1;
        else
            dim1 = aug_dim{l-1}+1;
        end
        
        dim2 = aug_dim{l};
        T2 = sum(U(dim1:dim2,i).^2);
        
        for b=1:B
            T2_star(b) = sum(U_star_matrix(dim1:dim2,i,b).^2);
            T2_null(b) = T2_star(b) - T2;
            if(abs(T2) <= abs(T2_null(b)))
                Indicator(b) = 1;
            else
                Indicator(b) = 0;
            end
        end
        p_value = sum(Indicator)/B;
        if(p_value >= P_fa_evec)
            present = 0;
        end
        U_struc(i,l) = present;
        
    end
end
end
