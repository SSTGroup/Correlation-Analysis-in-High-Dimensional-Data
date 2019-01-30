%   File: PCM_Bootstrap_Multiple_Datasets.m
%   Copyright (c) <2019> <University of Paderborn>
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
% SYNTAX:   [d_cap] = PCM_Bootstrap_Multiple_Datasets(X_cell, P_fa, B)
%
% INPUTS:
%
% 'X_cell'          Data sets in a cell format, each data set is an element of the cell
%
% 'P_fa'            Probability of false alarm for hypothesis testing 
%
% 'B'               Number of bootstrap iterations
%                   
% OUTPUTS:  
%
% 'd_cap'           Estimated number of correlated components
%            
%
% OVERVIEW: % MATLAB implementation of the technique in [1].
%
%
% DEPENDENCIES: None
%
%
% REFERENCES:
%
% The product of coherence matrices method for multiple data sets:
% [1]   T. Hasija, Y. Song, P. J. Schreier, and D. Ramírez, Bootstrap-based 
% detection of the number of signals correlated across multiple data sets",
% Asilomar Conf on Signals, Systems and Computers, Asilomar, PG California, 
% USA, 2016
%
%
%
% ------------------------------------------------------------------------
% CREATED:      07/01/2018 by Tanuj Hasija
%
% LAST EDITED:  07/01/2018 by Tanuj Hasija
%
% NOTES:
%
% ------------------------------------------------------------------------

function [d_cap] = PCM_Bootstrap_Multiple_Datasets(X_cell, P_fa, B)

L = size(X_cell,1); % Number of data sets
M = size(X_cell{1},2); % Numer of samples

V_cell{L} = []; % Contains the right singular vector matrices of L data sets
for l=1:L
    [~,~,V_temp] = svd(X_cell{l},'econ');
    V_cell{l} = V_temp;
end

[pairs] = bs_mat_org_Fixed_Edited(L);
[p,~] = size(pairs);

PCM_mat = eye(size(V_cell{pairs(1,1)}));
for i=1: p
    PCM_mat = PCM_mat*(V_cell{pairs(i,1)}'*V_cell{pairs(i,2)});
end
K = svd(PCM_mat);

V_cell_star{L} = [];
X_cell_star{L} = [];
% Bootstrap Operation
for b=1:B
    
    [~,I] = datasample(X_cell{1},M,2); % Bootstrap data indices
    
    for l=1:L
        X_cell_star{l} = X_cell{l}(:,I); 
        X_cell_star{l} = X_cell_star{l} - repmat(mean(X_cell_star{l},2),1,M);
        [~,~,V_temp] = svd(X_cell_star{l},'econ');
         V_cell_star{l} = V_temp;
    end
    
    PCM_mat = eye(size(V_cell_star{pairs(1,1)}));
    for i=1: p
        PCM_mat = PCM_mat*(V_cell_star{pairs(i,1)}'*V_cell_star{pairs(i,2)});
    end
    K_star_mat(:,b) = svd(PCM_mat);
    
end

d_cap = hypothesis_testing_bt(P_fa,K,K_star_mat,B);

end

function [pairs] = bs_mat_org_Fixed_Edited(L)
% Returns pairs of data sets whose coherence matrices participate in the
% product
if mod(L,2) == 1
    N = (L-1)/2;
    k = 1;
    for n = 1 : N
        p = 1;
        for q = 1 : gcd(n,L)
            i = q;
            j = 0;
            
            while j ~= q
                if i + n ~= L
                    j = mod(i+n,L);
                else
                    j = L;
                end
                pairs(k,:) = [i,j];
                i = j;
                p = p + 1;
                k = k + 1;
            end
            if q < gcd(n,L)
                pairs(k,:) = [q,q+1];
                p = p + 1;
                k = k + 1;
            elseif q ~= 1
                pairs(k,:) = [q,1];
                p = p + 1;
                k = k + 1;
            end
        end
    end
else
    N = L/2;
    L=L+1;
    k = 1;
    for n = 1 : N
        p = 1;
        for q = 1 : gcd(n,L)
            i = q;
            j = 0;
            while j ~= q
                if i + n ~= L
                    j = mod(i+n,L);
                else
                    j = L;
                end
                pairs(k,:) = [i,j];
                i = j;
                p = p + 1;
                k = k + 1;
            end
            if q < gcd(n,L)
                pairs(k,:) = [q,q+1];
                p = p + 1;
                k = k + 1;
            elseif q ~= 1
                pairs(k,:) = [q,1];
                p = p + 1;
                k = k + 1;
            end
        end
    end
    doubles = find(sum(pairs==L,2));
    pairs(doubles(1:2:end),2) = pairs(doubles(2:2:end),2);
    pairs(doubles(2:2:end),:) = [];
end
end

function d_cap = hypothesis_testing_bt(alpha,K,K_star_matrix,B)

m1 = length(K);
d_cap = m1-1;
for d=0:m1-1
    
    T = calc_statistic(K,d);
    
    for b=1:B
        T_star(b) = calc_statistic(K_star_matrix(:,b),d);
        T_null(b) = T_star(b) - T;
        if(abs(T) <= abs(T_null(b)))
            Indicator(b) = 1;
        else
            Indicator(b) = 0;
        end
    end
    p_value = sum(Indicator)/B;
    if(p_value >= alpha)
        d_cap = d;
        break;
    end
end


end

function T = calc_statistic(K,d)
m1 = length(K);

T = ( (1/(m1-d))*sum(K(d+1:m1)) ) - (prod(K(d+1:m1).^(1/(m1-d))) );

end

