%   File: MCCA_KPD.m
%
% OVERVIEW: % MATLAB implementation of MCCA-KPD Technique in "Estimation of common
%       subspace order across multiple datasets: Application to multi-subject
%       fMRI data," [1].
%
%
% DEPENDENCIES:
%
% [1] mcca.m
%
%       Not included within.
%       Written by Yiou Li (UMBC), see file for documentation.
%
%
% REFERENCES:
%
% MCCA-Knee point detection method:
% [1]   S. Bhinge, Y. Levin-Schwartz, and T. Adalı, “Estimation of common
%       subspace order across multiple datasets: Application to multi-subject
%       fMRI data,” in Proceedings of the 51st Annual Conference on Information
%       Sciences and Systems, 2017.
%
% ------------------------------------------------------------------------
% CREATED:      Suchita Bhinge
%               

function order = MCCA_KPD(X)
% X : Input data, N X V X K, K is the number of datasets and V is the number
% of samples

W = mcca(X,size(X,1),'genvar'); % perform MCCA 
%% Obtain SCVs
for k = 1 : size(W,3)
    Y(:,:,k) = W(:,:,k)*X(:,:,k);
    for n = 1 : size(W,1)
        SCV(k,:,n) = Y(n,:,k);
    end
end
clear Y
%% Obtain correlation of SCV
for n = 1 : size(W,1)
    CC_SCV(:,:,n) = corr(SCV(:,:,n)');
    j = 1;
    for k = 1 : size(W,3)
        for kk = 1 : size(W,3)
            if(k < kk)
                CC_SCV_vec(j,n) = abs(CC_SCV(k,kk,n));
                j = j + 1;
            end
        end
    end
end
clear SCV CC_SCV k kk n W
[tmp,~] = min(CC_SCV_vec);
tmp = sort(tmp,'descend');
%         figure();plot(tmp);
order = SORTE(tmp); 
end

function order = SORTE(CC)

M = length(CC);
CC = sort(CC,'descend');
for i = 1 : M-1
    diff_CC(i) = CC(i)-CC(i+1);
end

for n = 1 : M-3
    if (var(diff_CC(n:M-1)) > 0)
        sorte(n) = var(diff_CC(n+1:M-1))/var(diff_CC(n:M-1));
    elseif (var(diff_CC(n:M-1)) == 0)
        sorte(n) = +inf;
    end
end
[~,order] = min(sorte(1:M-3));
end