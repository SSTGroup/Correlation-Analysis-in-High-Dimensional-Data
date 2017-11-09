% This is a MATLAB implementation supporting the paper
%
% "Determining the Dimension of the Improper Signal
% Subspace in Complex-Valued Data" by Tanuj Hasija, Christian Lameiro and
% Peter J. Schreier, IEEE Signal Processing Letters, vol. 24, no. 11, pp. 1606-1610, Nov. 2017.
% 
% and implementation of a detector of the paper - 
%
% "Noncircular principal component analysis and its application to model
% selection", by X.-L. Li, T. Adali, and M. Anderson, IEEE Transactions on
% Signal Processing, vol. 59, no. 10, pp. 4516-4528, 2011
%
% ## ----------------------------------------------------------------------------
% ##
% ##   File: ncPCA_Detector.m
% ##   Copyright (c) <2017> Signal and System Theory Group, 
% ##                        Univ. of Paderborn, http://sst.upb.de
% ##                        https://github.com/SSTGroup/Correlation-Analysis-in-High-Dimensional-Data
% ##
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
% ##   3.) Publication, distribution, sublicensing, and/or selling of
% ##       copies or parts of the Software requires special agreements
% ##       with the Signal and System Theory Group, University of Paderborn,
% ##       and is in general not permitted.
% ##
% ##   4.) Modifications or contributions to the Software must be
% ##       published under this license. 
% ##   
% ##   5.) Any publication that was created with the help of this Software  
% ##       or parts of this Software must include a citation of the papers 
% ##       referenced above.
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
% ##   about bugs. 
% ##
% ##
% ##   Author: Tanuj Hasija, <tanuj.hasija@sst.upb>
% ##   Date: 18.08.2017
% ##   Version: v1 
% ## ----------------------------------------------------------------------------

function [n_cap, p_cap] = ncPCA_Detector(X)
% ncPCA_Detector.m - ncPCA Detector based on [9] for number of improper signals
%
%
% Input:
%
% X:         data matrix
%
% Output:
%
% p_cap:     estimated number of improper signals
% n_cap:     estimated total number of signals

T=size(X,2); % Number of samples
N=size(X,1); % Number of sensors
Chat=1/T*X*X'; % Sample covariance matrix
Phat=1/T*X*X.'; % Sample complementary covariance matrix

ITC = Inf+zeros(N);

for M=0:N-1
    
    idx1 = M+1;
    
    for p=0:M
        idx2 = p+1;
        
        [U,Lambda2,~]=eig(Chat);
        [~,I] = sort(abs(diag(Lambda2)),'descend');
        Uc=U(:,I(1:M));
        
        Lambda2c = real(diag(Uc'*Chat*Uc));
        Lambdac = sqrt(Lambda2c);
        sigma2c = real(trace(Chat) - trace(Uc'*Chat*Uc))/(N-M);
        
        Q = diag(1./Lambdac)*Uc'*Phat*conj(Uc)*diag(1./Lambdac);
        [V,~] = eig(Q*conj(Q));
        K = diag(V'*Q*conj(V));
        [~,I] = sort(abs(K),'descend');
        Kc = K(I(1:p));
        Vc = V(:,I(1:p));
        
        L(idx1,idx2) = ( (N-M)*log(sigma2c) + sum(log(Lambda2c)) + 0.5*sum(log(1-abs(Kc).^2)) );
        r(idx1,idx2) = M*(2*N - M) + p*(2*M -p+1);
        
        ITC(idx1,idx2) = L(idx1,idx2) + log(T)*r(idx1,idx2)/(2*T);
    end
    
end

ITC = real(ITC);
[I1, idx2] = min(ITC, [], 2);
[~, idx1] = min(I1, [], 1);
idx2=idx2(idx1);
n_cap = idx1-1;
p_cap = idx2-1;

end
