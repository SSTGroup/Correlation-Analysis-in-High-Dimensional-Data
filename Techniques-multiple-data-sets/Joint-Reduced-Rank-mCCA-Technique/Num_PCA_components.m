%   File: Num_PCA_components.m
%   Copyright (c) <2022> <University of Paderborn>
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
% OVERVIEW:
%
% This function implements the technique for estimating the number of
% signal components in a signal-plus-noise model from [1] 
%
%
% DEPENDENCIES:
%
% N/A
%
% REFERENCES:
%
% [1]  R. R. Nadakuditi and A. Edelman, “Sample eigenvalue based detection
% of high-dimensional signals in white noise using relatively few samples,”
% IEEE Transactions on Signal Processing, vol. 56, no. 7, pp. 2625–2638,
% 2008.
% ------------------------------------------------------------------------
% CREATED:      01/02/2022 by Tanuj Hasija
%
% LAST EDITED:  01/02/2022 by Tanuj hasija
%
function [Kmin obj] =  Num_PCA_components(L,n,m,beta)

% beta = 2; % For complex case
k_vector = 0:min(n,m)-1;

for i=1:length(k_vector)
    
    k = k_vector(i);
    tk = ((sum(L(k+1:n).^2)/sum(L(k+1:n))^2)*(n-k) - (1+ n/m))*n - (2/beta-1)*(n/m);
    obj(i) = ((beta/4)*((m/n)^2)*tk^2) + 2*(k+1);
    
end

[obj_min,index] = min(obj);
Kmin = k_vector(index);
end