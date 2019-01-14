% ## ----------------------------------------------------------------------------
% ##
% ##   File: CrossValidation.m
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
%   Matlab implementation of the cross-validation techniques for
%   model-order selection described in
%
%   Cross-validation techniques for determining the number of correlated components between two data sets when the number of samples is very small
%   Christian Lameiro and Peter J. Schreier
%   Proc. Asilomar Conf. Signals Syst. Computers, Pacific Grove, CA, USA, November 2016
%
%   Call:
%       [dhat,rhat]=CrossValidation(X,Y,k,rmax,mode)
%   Input:
%       X:    Channel one data
%       Y:    Channel two data
%       k:    Number of subsets (folds)
%       rmax: Maximum number of correlated components
%       mode: Rank-reduction: 'standard' (CV), 'PCA' (PCA-CV), 'XCOV'
%             (XCOV-CV), '3Set' (3Set-CV)
%   Output:
%       dhat: Estimated number of correlated components
%       rhat: Optimal rank of the rank-reduction step
%
%

function [dhat,rhat]=CrossValidation(X,Y,k,rmax,mode)

%% Initialitation
M=size(X,2); % Number of observations
n=size(X,1); % Data dimension of channel X
m=size(Y,1); % Data dimension of channel Y
SetLength=fix(M/k); % Length of the subsample sets
rmax=min(rmax,M-SetLength); % Correction of maximum number of correlated components
MaxIter=20; % Maximum number of iterations of alternating optimization (for 3set-CV)

if strcmpi(mode,'standard') % No rank reduction is performed
    rmin=rmax;
elseif strcmpi(mode,'PCA') || strcmpi(mode,'3set') || strcmpi(mode,'XCOV')
    rmin=1;
else
    error('Unknown rank-reduction option');
end
e=Inf+zeros(rmax,rmax); % Matrix containing the estimation errors
eI=zeros(rmax,rmax); % Auxiliary matrix
eP=zeros(rmax+1,1); % Auxiliary matrix

%% Cross-validation starts here
for r=rmin:rmax % We evaluate the error for all number of kept components
    if ~strcmpi(mode,'3set') && ~strcmpi(mode,'XCOV')
        eP=zeros(r+1,1);
    end
    for iSet=1:k % Model is estimated in all subsets except the kth one and validated in the kth one
        if iSet<k
            Xest=X(:,[1:(iSet-1)*SetLength iSet*SetLength+1:end]); % estimation samples
            Yest=Y(:,[1:(iSet-1)*SetLength iSet*SetLength+1:end]); % estimation samples
            Xval=X(:,(iSet-1)*SetLength+1:iSet*SetLength); % validation sample
            Yval=Y(:,(iSet-1)*SetLength+1:iSet*SetLength); % validation sample
        else
            Xest=X(:,1:(k-1)*SetLength); % model-estimation samples
            Yest=Y(:,1:(k-1)*SetLength); % model-estimation samples
            Xval=X(:,(k-1)*SetLength+1:end); % validation sample
            Yval=Y(:,(k-1)*SetLength+1:end); % validation sample
        end
        Mest=size(Xest,2); % Length of estimation set
        % Sample covariance matrices
        Ryy=1/Mest*Yest*Yest';
        Rxy=1/Mest*Xest*Yest';
        Rxx=1/Mest*Xest*Xest';
        [Ux,~,Uy]=svd(Rxy);
        
        if strcmpi(mode,'3set')           
            
            % Alternating optimization to find a local optimum of Fx and Fy            
            Fx=Ux(:,1:r); % Initial point for Fx
            Fy=Uy(:,1:r); % Initial point for Fy
            
            % Auxiliary matrices for fast computation of inner LOOCV
            Aux=1/Mest*Yest'*Fy*(Fy'*Ryy*Fy)^-1*Fy'*Yest;
            D1=diag(diag(Aux));
            D2=diag(diag(eye(Mest)-Aux));
            Ax=(Rxy*Fy*(Fy'*Ryy*Fy)^-1*Fy'*Yest-Xest*D1)*D2^-1;
            Zx=(2*Xest-Ax)*Ax';
            
            diffFx=Inf;
            diffFy=Inf;
            iter=1;
            
            while (diffFx>1e-3 || diffFy>1e-3) && iter<=MaxIter && r>0 % Convergence in terms of number of iterations and update of Fx and Fy matrices
                iter=iter+1;
                [U,E]=eig(-Zx);
                FxOld=Fx;
                [~,index]=sort(real(diag(E)));
                Fx=U(:,index(1:r)); % Fx update
                diffFx=abs(1-min(svd(FxOld'*Fx))); % Difference between new and old Fx in terms of spanned subspace
                
                % Auxiliary matrices for fast computation of inner LOOCV
                Aux=1/Mest*Xest'*Fx*(Fx'*Rxx*Fx)^-1*Fx'*Xest;
                D1=diag(diag(Aux));
                D2=diag(diag(eye(Mest)-Aux));
                Ay=(Rxy'*Fx*(Fx'*Rxx*Fx)^-1*Fx'*Xest-Yest*D1)*D2^-1;
                Zy=(2*Yest-Ay)*Ay';
                
                [U,E]=eig(-Zy);
                FyOld=Fy;
                [~,index]=sort(real(diag(E)));
                Fy=U(:,index(1:r));  % Fy update
                diffFy=abs(1-min(svd(FyOld'*Fy))); % Difference between new and old Fy in terms of spanned subspace
                
                % Auxiliary matrices for fast computation of inner LOOCV
                Aux=1/Mest*Yest'*Fy*(Fy'*Ryy*Fy)^-1*Fy'*Yest;
                D1=diag(diag(Aux));
                D2=diag(diag(eye(Mest)-Aux));
                Ax=(Rxy*Fy*(Fy'*Ryy*Fy)^-1*Fy'*Yest-Xest*D1)*D2^-1;
                Zx=(2*Xest-Ax)*Ax';
            end
            
            Exy=Xval-Fx*Fx'*Rxy*Fy*(Fy'*Ryy*Fy)^-1*Fy'*Yval; % Error estimating X from Y
            Eyx=Yval-Fy*Fy'*Rxy'*Fx*(Fx'*Rxx*Fx)^-1*Fx'*Xval; % Error estimating Y from X
            eP(r+1)=eP(r+1)+1/k*(1/n*trace(Exy*Exy')+1/m*trace(Eyx*Eyx'));
            if r==1
                Exy=Xval;
                Eyx=Yval;
                eP(r)=eP(r)+1/k*(1/n*trace(Exy*Exy')+1/m*trace(Eyx*Eyx'));
            end
            
        elseif strcmpi(mode,'XCOV')
            Fx=Ux(:,1:r);
            Fy=Uy(:,1:r);
            Exy=Xval-Fx*Fx'*Rxy*Fy*(Fy'*Ryy*Fy)^-1*Fy'*Yval; % Error estimating X from Y
            Eyx=Yval-Fy*Fy'*Rxy'*Fx*(Fx'*Rxx*Fx)^-1*Fx'*Xval; % Error estimating Y from X
            eP(r+1)=eP(r+1)+1/k*(1/n*trace(Exy*Exy')+1/m*trace(Eyx*Eyx'));
            if r==1
                Exy=Xval;
                Eyx=Yval;
                eP(r)=eP(r)+1/k*(1/n*trace(Exy*Exy')+1/m*trace(Eyx*Eyx'));
            end
        else
            [Ux,~,~]=svd(Xest);
            [Uy,~,~]=svd(Yest);            
            if strcmpi(mode,'standard')
                Fx=Ux;
                Fy=Uy;
            else % PCA-CV
                % PCA projectors
                Fx=Ux(:,1:r); 
                Fy=Uy(:,1:r);
                % Projected covariance matrices
                Rxx=Fx*Fx'*Rxx*Fx*Fx';
                Ryy=Fy*Fy'*Ryy*Fy*Fy';
            end            
            [U,D,V]=svd(Fx*Fx'*Rxy*Fy*Fy');
            for rxy=0:r % We evaluate the error for each guess on the number of correlated components
                Exy=Xval-U(:,1:rxy)*D(1:rxy,1:rxy)*V(:,1:rxy)'*pinv(Ryy)*Yval;
                Eyx=Yval-V(:,1:rxy)*D(1:rxy,1:rxy)*U(:,1:rxy)'*pinv(Rxx)*Xval;
                eP(rxy+1)=eP(rxy+1)+1/k*(1/n*trace(Exy*Exy')+1/m*trace(Eyx*Eyx'));
            end
        end
    end
    if ~strcmpi(mode,'3set') && ~strcmpi(mode,'XCOV')
        [e(r,r),eI(r,r)]=min(eP);
    end
end
% Estimation of the number of correlated components
if strcmpi(mode,'3set') || strcmpi(mode,'XCOV')
    [~,eI]=min(eP);
    dhat=eI-1;
    rhat=dhat;
else
    [minr,indr]=min(e,[],1);
    [~,indc]=min(minr);
    dhat=eI(indr(indc),indc)-1;
    rhat=indr(indc);
end

