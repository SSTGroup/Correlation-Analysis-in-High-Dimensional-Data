% ## ----------------------------------------------------------------------------
% ##
% ##   File: scca.m
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
%   Matlab implementation of the sparse CCA techniques for
%   model-order selection described in
%
%   A sparse CCA algorithm with application to model-order selection for small sample support
%   Christian Lameiro, and Peter J. Schreier
%   Proc. IEEE Int. Conf. Acoustics, Speech and Signal Process., New Orleans, LA, USA, March 2017
%
%   Call:
%       [k,S,T]=scca(X,Y,r);
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

function [k,S,T]=scca(X,Y,r)

%% Initialization

% General parameters
M=size(X,2); % Number of samples
n=size(X,1); % Dimension of channel x
m=size(Y,1); % Dimension of channel y
maxIter=20; % Maximum number of iterations

% Sample covariance matrices (ML estimates)
Rxx=1/M*X*X';
Ryy=1/M*Y*Y';
Rxy=1/M*X*Y';

% Auxiliar matrices
Is=eye(n);
It=eye(m);

[Ux,Dx,~]=svd(Rxx);
dx=diag(Dx);
dx(1:min(M,n))=1./dx(1:min(M,n));
piRxx=Ux*diag(dx)*Ux';

[Uy,Dy,~]=svd(Ryy);
dy=diag(Dy);
dy(1:min(M,m))=1./dy(1:min(M,m));
piRyy=Uy*diag(dy)*Uy';
if M<n
    Nxx=Ux(:,M+1:end);
else
    Nxx=zeros(n,0);
end
if M<m
    Nyy=Uy(:,M+1:end);
else
    Nyy=zeros(m,0);
end

%% Beginning of the algorithm

% Initial point
t=[1;zeros(m-1,1)];
t=t/sqrt(t'*Ryy*t);

% First canonical variable
conv=0;
costO=0;
iter=1;
while conv==0 && iter<=maxIter % Alternating optimization
    index=[]; % Vector containing indexes of non-zero components
    indexP=[]; % Vector containing indexes of positive components
    indexN=[]; % Vector containing indexes of negative components
    for iiS=1:r+1  % Subsequent identification of non-zero components for the first column of S  
        ival=setdiff(1:n,index);
        tIs=Is(ival,:);
        
        % Auxiliar matrices to compute closed-form solution
        M1=tIs*[piRxx , Nxx];
        Ms=pinv(M1);
        Mo=null(M1);
        Mm=Ms(1:n,:);
        m1=-tIs*piRxx*Rxy*t;
        Mo(1:n,:)=(eye(n)-Nxx*Nxx')*Mo(1:n,:);
        [U,~,~]=svd(Mo);
        Mo=U(:,1:length(index));
        Mmo=U(1:n,1:length(index));
        Mmoa=[Mmo(indexP,:) ; Mmo(indexN,:)];
        Mma=[Mm(indexP,:) ; Mm(indexN,:)];
        muTilde=(Mm-Mmo*Mmoa^-1*Mma)*m1;
        vec1=ones(n,1)-Mmo*Mmoa^-1*[ones(length(indexP),1);-ones(length(indexN),1)];
        vec1(index)=Inf;
        vec2=ones(n,1)+Mmo*Mmoa^-1*[ones(length(indexP),1);-ones(length(indexN),1)];
        vec2(index)=Inf;
        
        % Obtain minimum lambdax to keep iiS-1 non-zero components and
        % index of the next iiS-th non-zero component
        [lambdax,newIndex]=max([muTilde.*(muTilde>0)./vec1 ; -muTilde.*(muTilde<0)./vec2]);
        c=Mmoa^-1*(lambdax*[ones(length(indexP),1);-ones(length(indexN),1)]-Mma*m1);
        if newIndex>n
            indexN=[indexN ; newIndex-n];
        else
            indexP=[indexP ; newIndex];
        end
        index=[indexP ; indexN];        
    end    
    sol=Ms*m1+Mo*c; % Solution to the optimization problem
    muS=sol(1:n);
    phi=sol(n+1:end);
    saux=piRxx*(Rxy*t+muS)+Nxx*phi;
    s=saux/sqrt(saux'*Rxx*saux); % Sparse vector of canonical loadings


    index=[];
    indexP=[];
    indexN=[];
    for iiT=1:r+1 % Subsequent identification of non-zero components for the first column of T    
        ival=setdiff(1:m,index);
        tIt=It(ival,:);
        
        % Auxiliar matrices to compute closed-form solution
        M1=tIt*[piRyy , Nyy];
        Ms=pinv(M1);
        Mo=null(M1);
        Mm=Ms(1:m,:);
        m1=-tIt*piRyy*Rxy'*s;
        Mo(1:m,:)=(eye(m)-Nyy*Nyy')*Mo(1:m,:);
        [U,~,~]=svd(Mo);
        Mo=U(:,1:length(index));
        Mmo=U(1:n,1:length(index));
        Mmoa=[Mmo(indexP,:) ; Mmo(indexN,:)];
        Mma=[Mm(indexP,:) ; Mm(indexN,:)];
        muTilde=(Mm-Mmo*Mmoa^-1*Mma)*m1;
        vec1=ones(m,1)-Mmo*Mmoa^-1*[ones(length(indexP),1);-ones(length(indexN),1)];
        vec1(index)=Inf;
        vec2=ones(m,1)+Mmo*Mmoa^-1*[ones(length(indexP),1);-ones(length(indexN),1)];
        vec2(index)=Inf;
        
        % Obtain minimum lambday to keep iiT-1 non-zero components and
        % index of the next iiT-th non-zero component
        [lambday,newIndex]=max([muTilde.*(muTilde>0)./vec1 ; -muTilde.*(muTilde<0)./vec2]);
        c=Mmoa^-1*(lambday*[ones(length(indexP),1);-ones(length(indexN),1)]-Mma*m1);
        if newIndex>m
            indexN=[indexN ; newIndex-m];
        else
            indexP=[indexP ; newIndex];
        end
        index=[indexP ; indexN];        
    end    
    sol=Ms*m1+Mo*c; % Solution to the optimization problem
    muT=sol(1:m);
    phi=sol(m+1:end);
    taux=piRyy*(Rxy'*s+muT)+Nyy*phi;
    t=taux/sqrt(taux'*Ryy*taux); % Sparse vector of canonical loadings
    
    % Check convergence of alternating optimization
    costN=s'*Rxy*t;
    if abs(costN-costO)<1e-3 
        conv=1;
    else
        costO=costN;
    end
    iter=iter+1;
end
S=s;
T=t;

for ii=1:r-1 % Obtain next r-1 columns of S and T
    % Initial point
    As=Rxx*S;
    At=Ryy*T;
    aux=At*pinv(At'*At)^0.5;
    t=(eye(m)-aux*aux')*ones(m,1);
    t=real(t/sqrt(t'*Ryy*t));
    if isnan(t)
        t=randn(m,1);
        t=real(t/sqrt(t'*Ryy*t));
    end
    
    conv=0;
    costO=0;
    iter=1;
    
    % Auxiliar matrices
    [U,D,V]=svd(As);
    Ns=U(:,ii+1:end);
    
    [Ux,Dx,~]=svd(Ns'*Rxx*Ns);
    dx=diag(Dx);
    dx(1:min(M,n)-ii)=1./dx(1:min(M,n)-ii);
    piRxx=Ux*diag(dx)*Ux';
    Nxx=Ns*Ux(:,min(M,n)-ii+1:end);
    
    [U,D,V]=svd(At);
    Nt=U(:,ii+1:end);
    
    [Uy,Dy,~]=svd(Nt'*Ryy*Nt);
    dy=diag(Dy);
    dy(1:min(M,m)-ii)=1./dy(1:min(M,m)-ii);
    piRyy=Uy*diag(dy)*Uy';

    Nyy=Nt*Uy(:,min(M,m)-ii+1:end);
    
    while conv==0 && iter<=maxIter % Alternating optimization
        index=[];
        indexP=[];
        indexN=[];
        for iiS=1:r+1 % Subsequent identification of non-zero components for the (r+1)-th column of S 
            ival=setdiff(1:n,index);
            tIs=Is(ival,:);

            % Auxiliar matrices to compute closed-form solution
            M1=tIs*[Ns*piRxx*Ns' , Nxx];
            Ms=pinv(M1);
            Mo=null(M1);
            Mm=Ms(1:n,:);
            m1=-tIs*Ns*piRxx*Ns'*Rxy*t;
            Mo(1:n,:)=(eye(n)-Nxx*Nxx')*Mo(1:n,:);
            [U,~,~]=svd(Mo);
            Mo=U(:,1:length(index)); 
            Mmo=U(1:n,1:length(index));
            Mmoa=[Mmo(indexP,:) ; Mmo(indexN,:)];
            Mma=[Mm(indexP,:) ; Mm(indexN,:)];
            muTilde=(Mm-Mmo*Mmoa^-1*Mma)*m1;
            vec1=ones(n,1)-Mmo*Mmoa^-1*[ones(length(indexP),1);-ones(length(indexN),1)];
            vec1(index)=Inf;
            vec2=ones(n,1)+Mmo*Mmoa^-1*[ones(length(indexP),1);-ones(length(indexN),1)];
            vec2(index)=Inf;

        
            % Obtain minimum lambdax to keep iiS-1 non-zero components and
            % index of the next iiS-th non-zero component
            [lambdax,newIndex]=max([muTilde.*(muTilde>0)./vec1 ; -muTilde.*(muTilde<0)./vec2]);
            c=Mmoa^-1*(lambdax*[ones(length(indexP),1);-ones(length(indexN),1)]-Mma*m1);
            if newIndex>n
                indexN=[indexN ; newIndex-n];
            else
                indexP=[indexP ; newIndex];
            end
            index=[indexP ; indexN]; 

        end
        sol=Ms*m1+Mo*c; % Solution to the optimization problem
        muS=sol(1:n);
        phi=sol(n+1:end);
        saux=Ns*piRxx*Ns'*(Rxy*t+muS)+Nxx*phi;
        s=saux/sqrt(saux'*Rxx*saux); % Sparse vector of canonical loadings

        
        index=[];
        indexP=[];
        indexN=[];
        for iiT=1:r+1 % Subsequent identification of non-zero components for the (r+1)-th column of T
            ival=setdiff(1:m,index);
            tIt=Is(ival,:);

            % Auxiliar matrices to compute closed-form solution
            M1=tIt*[Nt*piRyy*Nt' , Nyy];
            Ms=pinv(M1);
            Mo=null(M1);
            Mm=Ms(1:m,:);
            m1=-tIt*Nt*piRyy*Nt'*Rxy'*s;
            Mo(1:m,:)=(eye(m)-Nyy*Nyy')*Mo(1:m,:);
            [U,~,~]=svd(Mo);
            Mo=U(:,1:length(index)); 
            Mmo=U(1:m,1:length(index));
            Mmoa=[Mmo(indexP,:) ; Mmo(indexN,:)];
            Mma=[Mm(indexP,:) ; Mm(indexN,:)];
            muTilde=(Mm-Mmo*Mmoa^-1*Mma)*m1;
            vec1=ones(m,1)-Mmo*Mmoa^-1*[ones(length(indexP),1);-ones(length(indexN),1)];
            vec1(index)=Inf;
            vec2=ones(m,1)+Mmo*Mmoa^-1*[ones(length(indexP),1);-ones(length(indexN),1)];
            vec2(index)=Inf;

            % Obtain minimum lambdax to keep iiT-1 non-zero components and
            % index of the next iiT-th non-zero component
            [lambday,newIndex]=max([muTilde.*(muTilde>0)./vec1 ; -muTilde.*(muTilde<0)./vec2]);
            c=Mmoa^-1*(lambday*[ones(length(indexP),1);-ones(length(indexN),1)]-Mma*m1);
            if newIndex>m
                indexN=[indexN ; newIndex-n];
            else
                indexP=[indexP ; newIndex];
            end
            index=[indexP ; indexN]; 

        end
        sol=Ms*m1+Mo*c; % Solution to the optimization problem
        muT=sol(1:m);
        phi=sol(m+1:end);
        taux=Nt*piRyy*Nt'*(Rxy'*s+muT)+Nyy*phi;
        t=taux/sqrt(taux'*Ryy*taux); % Sparse vector of canonical loadings

        % Check convergence of alternating optimization
        costN=s'*Rxy*t;
        if abs(costN-costO)<1e-3
            conv=1;
        else
            costO=costN;
        end
        iter=iter+1;
    end
    if costN==0 % If optimal s or t is an all-zero vector we stop
        break;
    end
    S=[S s];
    T=[T t];
end
k=diag(S'*Rxy*T); % Estimated canonical correlations
