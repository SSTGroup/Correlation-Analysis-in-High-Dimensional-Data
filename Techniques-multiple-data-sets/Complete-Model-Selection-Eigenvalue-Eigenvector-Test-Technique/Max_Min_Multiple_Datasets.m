%   File: Max_Min_Multiple_Datasets.m
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
% SYNTAX:   [d_cap, r_d] = Max_Min_Multiple_Datasets(X_cell, P_fa, rmax, RealComp) 
%
% INPUTS:
%
% 'X_cell'          Data sets in a cell format, each data set is an element of the cell
%
% 'P_fa'            Probability of false alarm for hypothesis testing
%
% 'rmax'            Maximum PCA rank
%
% 'RealComp'        'real' for real-valued data and 'comp' for complex-valued data 
%                   
% OUTPUTS:  
%
% 'd_cap'           Estimated number of correlated components
%
% 'r_d'             Estimated PCA rank required to contain d_cap components
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
% The maxmin method for multiple data sets:
% [1]   T. Hasija, Y. Song, P. J. Schreier, and D. Ramírez, “Detecting the
%       dimension of the subspace correlated across multiple data sets in
%       the sample poor regime,” in Proceedings of the IEEE Workshop on
%       Statistical Signal Processing, 2016.
%
%
%
% ------------------------------------------------------------------------
% CREATED:      14/01/2019 by Tanuj Hasija
%
% LAST EDITED:  14/01/2019 by Tanuj Hasija
%
% NOTES:
%
% ------------------------------------------------------------------------

function [d_cap, r_d] = Max_Min_Multiple_Datasets(X_cell, P_fa, rmax, RealComp) 

% MATLAB Implementation based on the work "Detecting the dimension of the
% subspace correlated across multiple data sets in the sample poor regime",
% T. Hasija, Y. Song, P. J. Schreier, and D. Ramirez,IEEE Workshop on
% Statistical Signal Processing, Palma de Mallorca, Spain, 2016

L = size(X_cell,1); % Number of data sets
M = size(X_cell{1},2); % Numer of samples
V_cell = [];

for l=1:L
    [~,~,V_temp] = svd(X_cell{l},'econ');
    V_cell{l} = V_temp;
%     display(['svd loop 1, number=' num2str(l)]);
end

for i=1:L-2
    X_vec = [];
    for j=i+1:L
        X_vec = [X_vec;X_cell{j}];
    end
    [~,~,V_temp] = svd(X_vec,'econ');
    V_cell{end+1} = V_temp;
%      display(['svd loop 2, number=' num2str(i)]);
end

d_cap_vector = zeros(1,rmax);

for r=1: rmax
%     display(['rank=' num2str(r)]);
    K_cell = [];
    for i=1:L-1
        if (i~= L-1)
            K_cell{i} = svd(V_cell{i}(:,1:r)'*V_cell{L+i}(:,1:r));
        else
            K_cell{i} = svd(V_cell{i}(:,1:r)'*V_cell{i+1}(:,1:r));
        end
    end
    
    for s=0:r
        C_s_r = 0;
        
        for i=1:L-1
            switch lower(RealComp)
                case 'real'
                    C_s_r = C_s_r - (M-(2*r+1)/2)*log(prod((1-K_cell{i}(s+1:r).^2)));
                case 'comp'
                    C_s_r = C_s_r - (2*M-(2*r+1))*log(prod((1-K_cell{i}(s+1:r).^2)));
            end
        end
        
        switch lower(RealComp)
            case 'real'
                dof = (L-1)*((r-s)^2);
            case 'comp'
                dof = 2*(L-1)*((r-s)^2);
        end
        
        T_s_r = chi2inv(1-P_fa,dof);
        if (C_s_r < T_s_r)
            break;
        end
    end
    d_cap_vector(r) = s;
end

[d_cap, r_d] = max(d_cap_vector);

end
