% ## ----------------------------------------------------------------------------
% ##
% ##   File: IC_maxmin.m
% ##   Copyright (c) <2016> <University of Paderborn>
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
% ##
% ## ----------------------------------------------------------------------------
%
%
%
% IC_maxmin.m - MDL based max-min method
%
% Thre                'MDL' or 'AIC' ('AIC' is used only for test)
% RealComp            'real' or 'comp'
% rmax                {r1,r2} = 1,...,rmax
%
%
%
% Input:
%
% Vx:        contains right singular vectors of X
% Vy:        contains right singular vectors of Y
% M:         number of samples
% Thre:     'MDL' ('AIC' is used only for test)
% RealComp: 'real' for real-valued data;
%           'comp' for complex-valued data
% f:         f(1) is the number of independent signals in X
%            f(2) is the number of independent signals in Y
% rmax       {r1,r2} = 1,...,rmax
%             r1: the rank PCA keeps in X
%             r2: the rank PCA keeps in Y
%
% Output:
%
% d:         number of correlated signal between X and Y
% r1Est:     the optimum rank that PCA should keep for X
% r2Est:     the optimum rank that PCA should keep for Y

function [d, r1Est, r2Est] = IC_maxmin(Vx,Vy,M,Thre,RealComp,f,rmax)

for r1 = 1:rmax
    % If f1=f2, we let r1=r2.
    if f(1)==f(2), r2range  = r1; else r2range = 1:rmax; end
    for r2 = r2range
        
        ga = sort(svd(Vx(:,1:r1)'*Vy(:,1:r2)),'descend');
        for r3 = 0:min(r1,r2)-1
            switch lower(RealComp)
                case 'real'
                    dof = r1*r2 - (r1-r3)*(r2-r3);
                    Loglike = M/2*log(prod(1-ga(1:r3).^2));
                case 'comp'
                    dof = 2*r3*(r1+r2-r3);
                    Loglike = M*log(prod(1-ga(1:r3).^2));
            end
            switch lower(Thre)
                case 'mdl'
                    pen = 1/2*log(M)*dof;
                case 'mdlbywu'
                    pen = 0.65*sqrt(M)/log(M)*dof;
                case 'aic'
                    pen = dof;
            end
            IC(r3+1) = Loglike + pen;
        end
        if min(IC)==-inf,
            idx(r1,r2)=0;
        else
            idx(r1,r2) = find(IC==min(IC))-1;
        end
        clear IC
    end
end
d = max(max(idx));

[~, location] = max(idx(:));
[r1Est,r2Est] = ind2sub(size(idx),location);
