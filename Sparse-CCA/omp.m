% ## ----------------------------------------------------------------------------
% ##
% ##   File: omp.m
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
%   Orthogonal matching pursuit algorithm to find an sparse vector b that
%   minimizes ||y-Xb||^2.
%
%   Call:
%       b=omp(y,X,k);
%   Input:
%       y: data vector
%       X: measurement matrix
%       k: number of non-zero elements of b
%   Output:
%       b: sparse solution
%
%

function b=omp(y,X,k)

for i=1:size(X,2)
    if sum(abs(X(:,i)))>1e-3
        X(:,i)=X(:,i)/norm(X(:,i)); % column normalization
    end
end
r=y; % Initial residual
c=[]; % Initial set of selected variables
for i=1:k
    [m,t]=max(abs(X'*r));
    c=[c,t]; % Update selected variables
    P=X(:,c)*(X(:,c)'*X(:,c))^-1*X(:,c)';
    r=(eye(length(P))-P)*y; % Update residual
end
Xt=X(:,c);
b=zeros(size(X,2),1);
bnz=(Xt'*Xt)^-1*Xt'*y;
b(c)=bnz;
 