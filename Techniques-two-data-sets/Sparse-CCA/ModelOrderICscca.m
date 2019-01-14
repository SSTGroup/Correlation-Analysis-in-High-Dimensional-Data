% ## ----------------------------------------------------------------------------
% ##
% ##   File: ModelOrderICscca.m
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
% ##   Author: Christian Lameiro
% ##
% ## ----------------------------------------------------------------------------
%
%
%
%   Matlab script to estimate, from canonical correlations obtained through a SCCA algorithm, the number of 
%   correlated components by MDL criterion.
%
%   Call:
%       dhat=ModelOrderICscca(k,M);
%   Input:
%       k:    matrix of canonical correlations. The ith column contains the
%             canonical correlations, in decreasing order, for i non-zero
%             components in the SCCA algorithm
%       M:    number of observations
%   Output:
%       dhat: estimated number of correlated components
%
%

function dhat=ModelOrderICscca(k,M)

rmax=length(k);
MaxScore=zeros(rmax+1,1);
dopt=zeros(rmax+1,1);

for r=1:rmax
    aux=zeros(r,1);
    for d=1:r
        aux(d)=-M/2*log(prod(1-k(1:d,r).^2))-1/2*log(M)*(4*r-2*d+1).*d;
    end
    [MaxScore(r+1),dopt(r+1)]=max(aux);
    if MaxScore(r+1)<0
        MaxScore(r+1)=0;
        dopt(r+1)=0;
    end    
end

[~,ind]=max(MaxScore);
dhat=dopt(ind);
