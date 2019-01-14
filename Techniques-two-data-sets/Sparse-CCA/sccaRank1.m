% ## ----------------------------------------------------------------------------
% ##
% ##   File: sccaRank1.m
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
%   Matlab implementation of the sparse CCA algorithm described in
%
%   Sparse canonical correlation analysis based on rank-1 matrix approximation and its application for fMRI signals 
%   A. Aïssa-El-Bey and A. K. Seghouane
%   Proc. IEEE Int. Conf. Acoustics, Speech and Signal Process., Shanghai, China, March 2016, pp. 4678? 4682
%
%   Call:
%       [k,S,T]=sccaRank1(X,Y,r);
%   Input:
%       X: X-channel data
%       Y: Y-channel data
%       r: number of non-zero components of each column of S and T
%   Output:
%       k: vector of estimated canonical correlations
%       S: canonical loadings for channel X
%       T: canonical loadings for channel Y
%
%

function [k,S,T]=sccaRank1(X,Y,r)

%% General parameters
M=size(X,2); % Number of samples
n=size(X,1); % Dimension of channel x
m=size(Y,1); % Dimension of channel y

%% Sample covariance matrices (ML estimates)
Rxx=1/M*X*X';
Ryy=1/M*Y*Y';
Rxy=1/M*X*Y';

[Ux,D,~]=svd(Rxx);
d=diag(D);
d(1:min(M,n))=1./d(1:min(M,n));
piRxx=Ux*diag(d)*Ux'; % Pseudoinverse of Rxx

[Uy,D,~]=svd(Ryy);
d=diag(D);
d(1:min(M,m))=1./d(1:min(M,m));
piRyy=Uy*diag(d)*Uy'; % Pseudoinverse of Ryy

%% Beginning of the algorithm
Px=1/M*X'*piRxx*X;
Py=1/M*Y'*piRyy*Y;
Kxy=Px*Py;
S=zeros(n,r);
T=zeros(m,r);
for i=1:r
    [U,D,V]=svd(Kxy);
    ut=U(:,1);
    vt=V(:,1);
    for it=1:20
        S(:,i)=omp(Kxy*vt,X',r); % Orthogonal matching pursuit algorithm
        ut=X'*S(:,i)/norm(X'*S(:,i));
        T(:,i)=omp(Kxy'*ut,Y',r); % Orthogonal matching pursuit algorithm
        vt=Y'*T(:,i)/norm(Y'*T(:,i));
    end
    Kxy=Kxy-(ut'*Kxy*vt)*ut*vt';
end
for i=1:r
    S(:,i)=S(:,i)/sqrt(S(:,i)'*Rxx*S(:,i));
    T(:,i)=T(:,i)/sqrt(T(:,i)'*Ryy*T(:,i));
end
k=sort(abs(diag(S'*Rxy*T)),'descend');
        