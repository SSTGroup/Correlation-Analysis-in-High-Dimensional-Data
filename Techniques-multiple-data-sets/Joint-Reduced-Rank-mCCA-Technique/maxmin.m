% ## ----------------------------------------------------------------------------
% ##
% ##   File: maxmin.m
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
% ##   Author: Yang Song

function [d, r1Est, r2Est] = maxmin(Vx,Vy,M,Pfa,RealComp,equal_range,rmax,ITC)
% IC_maxmin.m - Hypothesis test based max-min method
%
%
% Input:
%
% Vx:        contains right singular vectors of X
% Vy:        contains right singular vectors of Y
% M:         number of samples
% RealComp: 'real' for real-valued data; 
%           'comp' for complex-valued data
% equal_range:         1 if PCA ranks for both data sets are assumed equal and 0 otherwise
%                        (saves computaion time without siginificant performance loss)
% rmax       {r1,r2} = 1,...,rmax
%             r1: the rank PCA keeps in X
%             r2: the rank PCA keeps in Y
%
% ITC        if 'HT', the threshold for the hypothesis test is determined
%            by a given probability of false alarm, like the traditional 
%            hypothesis test
%
%            if 'MDL', the threshold for the hypothesis test is determined 
%            by MDL
%
%
%
% Output:
%
% d:         number of correlated signal between X and Y
% r1Est:     the optimum rank that PCA should keep for X
% r2Est:     the optimum rank that PCA should keep for Y
%
% Edited 
% 01/22/2022  Tanuj Hasija 

for r1 = 1:rmax
    if (equal_range), r2range  = r1; else, r2range = 1:rmax; end
    for r2 = r2range
        r = min(r1,r2);
        C = zeros(r,1);
        T = zeros(r,1);
        ga = svd(Vx(:,1:r1)'*Vy(:,1:r2));
        for r3 = 0:r-1
            switch lower(RealComp)
                case 'real'
                    switch lower(ITC)
                        case 'ht'
                            C(r3+1) = -(M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))*log(prod(1-ga(r3+1:r).^2));
                            T(r3+1) = chi2inv(1-Pfa,(r1-r3)*(r2-r3));
                        case 'mdl'
                            C(r3+1) = -(M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))*log(prod(1-ga(r3+1:r).^2));
                            T(r3+1) = (M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))/M*log(M)*(r1-r3)*(r2-r3); % MDL threshold
                        case 'aic'
                            C(r3+1) = -(M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))*log(prod(1-ga(r3+1:r).^2));
                            T(r3+1) = (M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))/M*2*(r1-r3)*(r2-r3); % AIC threshold
                    end
                case 'comp'
                    switch lower(ITC)
                        case 'ht'
                            C(r3+1) = -2*(M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))*log(prod(1-ga(r3+1:r).^2));
                            T(r3+1) = chi2inv(1-Pfa,2*(r1-r3)*(r2-r3));
                        case 'mdl'
                            C(r3+1) = -2*(M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))*log(prod(1-ga(r3+1:r).^2));
                            T(r3+1) = (M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))/M*log(M)*2*(r1-r3)*(r2-r3); % MDL threshold
%                             T(r3+1) = (M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))/M*1.5*sqrt(M)/log(M)*2*(r1-r3)*(r2-r3); % MDL threshold
                        case 'aic'
                            C(r3+1) = -2*(M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))*log(prod(1-ga(r3+1:r).^2));
                            T(r3+1) = (M-r3-1/2*(r1+r2+1)+sum(ga(1:r3).^-2))/M*2*2*(r1-r3)*(r2-r3); % AIC threshold
                    end
            end
        end
        IDX = find(C-T<0);
        if isempty(IDX); smin(r1,r2) = r;
        else smin(r1,r2) = IDX(1)-1; end
    end
end
d = max(max(smin));
[~, location] = max(smin(:));
[r1Est,r2Est] = ind2sub(size(smin),location);