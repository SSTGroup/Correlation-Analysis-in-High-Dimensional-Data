%   File: ITC_MDL_Wu_Multiple_Datasets.m
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
% SYNTAX:   [d_cap] = ITC_MDL_Wu_Multiple_Datasets(X_cell, RealComp)
%
% INPUTS:
%
% 'X_cell'          Data sets in a cell format, each data set is an element of the cell
%
% 'RealComp'        'real' for real-valued data and 'comp' for complex-valued data
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
% The ITC MDL method for multiple data sets:
% [1]   Y.Wu, K.W.Tam, and F.Li, "Determination of number of sources with
% multiple arrays in correlated noise fields", IEEE Transactions on Signal
% Processing, 2002
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

function [d_cap] = ITC_MDL_Wu_Multiple_Datasets(X_cell, RealComp) %detector_type

%%%%%%%%%% Input Parametes %%%%%%%%%%%
% X_cell: Data sets in a cell format, each data set is an element of the cell
% RealComp: 'real' for real-valued data sets and 'complex' for
% complex-valued data sets

L = size(X_cell,1); % Number of data sets
M = size(X_cell{1},2); % Numer of samples
V_cell = [];

for l=1:L
    [~,~,V_temp] = svd(X_cell{l},'econ');
    V_cell{l} = V_temp;
    m(l) = size(X_cell{l},1); % dimension of the data sets
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

d_max = min(m)-1;
d_vector = 0:d_max;
ITC = zeros(size(d_vector));
d_cap_vector = zeros(1,min(m));

K_cell = [];
for i=1:L-1
    if (i~= L-1)
        K_cell{i} = svd(V_cell{i}'*V_cell{L+i});
    else
        K_cell{i} = svd(V_cell{i}'*V_cell{i+1});
    end
end

for i=1:length(d_vector)
    s = d_vector(i);
    L_X = 0;
    V_n = 0;
    
    for r=1:L-1
        p = min(m(r:end)); % represents p_j in paper, which changes with iteration
        L_X = L_X + log(prod(1-K_cell{r}(s+1:p).^2));
        m1 = m(r); m2 = sum(m(r+1:end));
        switch lower(RealComp)
            case 'real'
                V_n =  V_n + (m1*m2 -(m1-s)*(m2-s) + m1*(m1+1)/2 + m2*(m2+1)/2);
            case 'comp'
                V_n =  V_n + 2*s*(m1+ m2-s);
        end
        
    end
    
    switch lower(RealComp)
        case 'real'
            obj(i) = -M/2*L_X + V_n*log(M)/2;
        case 'comp'
            obj(i) = -M*L_X + V_n*log(M)/2;
    end
    
    
end

[~,index] = min(obj);
d_cap = d_vector(index);


end





