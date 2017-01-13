% ## ----------------------------------------------------------------------------
% ##
% ##   File: Main.m
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
%   Main file to evaluate the cross-validation techniques for
%   model-order selection described in
%
%   Cross-validation techniques for determining the number of correlated components between two data sets when the number of samples is very small
%   Christian Lameiro and Peter J. Schreier
%   Proc. Asilomar Conf. Signals Syst. Computers, Pacific Grove, CA, USA, November 2016
%
%   Different scenarios can be provided using an specific format that is
%   described next. The variable parameter is detected automatically
%   (example, the number M of observations, or the variance of independent components).
%
%   The names of the parameters match those in the paper:
%
%       n:        data dimension of channel X
%       m:        data dimension of channel Y
%       M:        number of observations
%       d:        number of correlated components
%       fx:       number of indendent signals in channel X
%       fy:       number of independent signals in channel Y
%       p:        correlation coefficients
%       sigmadx2: variance of correlated signals in channel X
%       sigmady2: variance of correlated signals in channel Y
%       sigmafx2: variance of independent signals in channel X
%       sigmafy2: variance of independent signals in channel Y
%       sigma2:   noise variance
%       mixing:   type of mixing matrix ('orth' for random unitary matrix and 'randn' for normal entries)
%       noise:    'white' for white noise and 'colored' for colored noise
%       MAcoeff:  coefficients of MA model for colored noise
%       ARcoeff:  coefficients of AR model for colored noise
%
%   For multiple correlated and/or independent signals, their corresponding
%   variance and/or correlation coefficients must have the corresponding
%   dimension and be provided as column vectors. Similarly, MAcoeff and
%   ARcoeff must be column vectors if they are multidimensional. If there
%   is a variable parameter, their different values must be provided along
%   the second dimension (e.g., M=[10 20 30];). There can be up to one variable parameter. There are also some special
%   cases to consider:
%
%       - Use mn or nm as a variable when n=m is a variable
%       parameter (to evaluate the performance versus the dimension).
%       - Use f as a variable when fx=fy is a variable paramenter (to
%       evaluate the performance versus the number of independent
%       components in both channels).
%       - Use sigmad2 as a variable when sigmadx2=sigmady2 is a variable
%       parameter.
%       - Use sigmaf2 as a variable when sigmafx2=sigmafy2 is a variable
%       parameter.
%       - When fx, fy or d are variable parameters, provide the
%       corresponding variance and/or correlation coefficient for the
%       maximum fx, fy or d.
%
%   Please see the provided examples. 'scen1', 'scen2' and 'scen3' generate
%   Fig. 3, Fig. 4 and Fig. 5 of the paper, respectively.
%

%% Initialization
clc;clear;close all;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));
nSim=100; % Number of Monte-Carlo simulations

%% Scenario specification
scen='scen1';
switch lower(scen)
    case 'custom'
        n=20; % Data dimension in channel X
        m=20; % Data dimension in channel Y
        M=10; % Number of observations
        d=1; % Number of correlated components
        fx=4; % Number of independent signals in channel X
        fy=4; % Number of independent signals in channel Y
        p=0.8*ones(max(d),1); % Correlation coefficients
        
        sigmadx2=4*ones(max(d),1); % Variance of correlated signals in channel X
        sigmafx2=1*ones(max(fx),1); % Variance of independent signals in channel X
        sigmady2=4*ones(max(d),1); % Variance of correlated signals in channel Y
        sigmafy2=1*ones(max(fy),1); % Variance of independent signals in channel Y
        sigma2=0.1; % Noise variance
        mixing='orth'; % Type of mixing matrix ('orth' for random unitary matrix and 'randn' for normal entries)
        noise='white'; % 'white' for white noise and 'colored' for colored noise
        MAcoeff=1; % Coefficients of MA model for colored noise
        ARcoeff=1; % Coefficients of AR model for colored noise
    case 'scen1'
        n=20; % Data dimension in channel X
        m=20; % Data dimension in channel Y
        M=10:5:30; % Number of observations
        d=1; % Number of correlated components
        fx=4; % Number of independent signals in channel X
        fy=4; % Number of independent signals in channel Y
        p=0.7*ones(max(d),1); % Correlation coefficients
        
        sigmadx2=4*ones(max(d),1); % Variance of correlated signals in channel X
        sigmafx2=1*ones(max(fx),1); % Variance of independent signals in channel X
        sigmady2=4*ones(max(d),1); % Variance of correlated signals in channel Y
        sigmafy2=1*ones(max(fy),1); % Variance of independent signals in channel Y
        sigma2=0.1; % Noise variance
        mixing='orth'; % Type of mixing matrix ('orth' for random unitary matrix and 'randn' for normal entries)
        noise='white'; % 'white' for white noise and 'colored' for colored noise
        MAcoeff=1; % Coefficients of MA model for colored noise
        ARcoeff=1; % Coefficients of AR model for colored noise
    case 'scen2'
        n=40; % Data dimension in channel X
        m=40; % Data dimension in channel Y
        M=40:10:100; % Number of observations
        d=2; % Number of correlated components
        fx=20; % Number of independent signals in channel X
        fy=20; % Number of independent signals in channel Y
        p=[0.92;0.9]; % Correlation coefficients
        
        sigmadx2=10*ones(max(d),1); % Variance of correlated signals in channel X
        sigmafx2=8*ones(max(fx),1); % Variance of independent signals in channel X
        sigmady2=10*ones(max(d),1); % Variance of correlated signals in channel Y
        sigmafy2=8*ones(max(fy),1); % Variance of independent signals in channel Y
        sigma2=1; % Noise variance
        mixing='orth'; % Type of mixing matrix ('orth' for random unitary matrix and 'randn' for normal entries)
        noise='white'; % 'white' for white noise and 'colored' for colored noise
        MAcoeff=1; % Coefficients of MA model for colored noise
        ARcoeff=1; % Coefficients of AR model for colored noise
    case 'scen3'
        n=20; % Data dimension in channel X
        m=20; % Data dimension in channel Y
        M=20; % Number of observations
        d=1; % Number of correlated components
        fx=2; % Number of independent signals in channel X
        fy=2; % Number of independent signals in channel Y
        p=0.8; % Correlation coefficients
        
        sigmadx2=5*ones(max(d),1); % Variance of correlated signals in channel X
        sigmaf2=1:10; % Variance of independent signals in channel X and Y
        sigmady2=5*ones(max(d),1); % Variance of correlated signals in channel Y
        sigma2=1; % Noise variance
        mixing='orth'; % Type of mixing matrix ('orth' for random unitary matrix and 'randn' for normal entries)
        noise='white'; % 'white' for white noise and 'colored' for colored noise
        MAcoeff=1; % Coefficients of MA model for colored noise
        ARcoeff=1; % Coefficients of AR model for colored noise
    otherwise
        error('Unknown scenario');
end
wspace=whos;
variableIndex=cellfun(@strcmpi,{wspace.class},mat2cell(repmat('double',1,length(wspace)),1,6*ones(length(wspace),1)));
wspace=wspace(variableIndex);

for ii=1:length(wspace)
    switch wspace(ii).name
        case 'n'
            wspace(ii).label='data dimension of channel X (n)';
        case 'm'
            wspace(ii).label='data dimension of channel Y (m)';
        case 'd'
            wspace(ii).label='number of correlated components (d)';
        case 'fx'
            wspace(ii).label='number of independent components in channel X (f_x)';
        case 'fy'
            wspace(ii).label='number of independent components in channel Y (f_y)';
        case 'M'
            wspace(ii).label='number of observations (M)';
        case 'sigmadx2'
            wspace(ii).label='variance of correlated components in channel X (\sigma_x^2)';
        case 'sigmady2'
            wspace(ii).label='variance of correlated components in channel Y (\sigma_y^2)';
        case 'sigmafx2'
            wspace(ii).label='variance of independent components in channel X';
        case 'sigmafy2'
            wspace(ii).label='variance of independent components in channel Y';
        case 'sigma2'
            wspace(ii).label='noise variance (\sigma^2)';
        case 'p'
            wspace(ii).label='correlation coefficient (\rho)';
        case 'nm' 
            wspace(ii).label='data dimension (n,m)';
        case 'mn' 
            wspace(ii).label='data dimension (n,m)';
        case 'sigmad2'
            wspace(ii).label='variance of correlated components (\sigma_x^2,\sigma_y^2)';
        case 'sigmaf2'
            wspace(ii).label='variance of independent components';
    end
end

aux=cell2mat({wspace.size}');
[sweepL,sweepI]=max(aux(:,2));

eval(['sweepVar=',wspace(sweepI).name,';']);

dEstCVstd=zeros(sweepL,nSim);
dEstCVpca=zeros(sweepL,nSim);
dEstCV3set=zeros(sweepL,nSim);
dEstCVxcov=zeros(sweepL,nSim);
dEstDet2=zeros(sweepL,nSim);
dEstDet3=zeros(sweepL,nSim);

%% Monte-Carlo simulations for the specified scenario
for sim=1:nSim
    s=sprintf('Simulation %i/%i',sim,nSim);
    fprintf(s);
    
    for ii=1:sweepL
        eval([wspace(sweepI).name,'=sweepVar(ii);']);
        if strcmpi(wspace(sweepI).name,'f')
            fx=f;
            fy=f;
        elseif strcmpi(wspace(sweepI).name,'sigmad2')
            sigmadx2=repmat(sigmad2,[d 1]);
            sigmady2=repmat(sigmad2,[d 1]);
        elseif strcmpi(wspace(sweepI).name,'sigmaf2')
            sigmafx2=repmat(sigmaf2,[fx 1]);
            sigmafy2=repmat(sigmaf2,[fy 1]);
        elseif strcmpi(wspace(sweepI).name,'nm') || strcmpi(wspace(sweepI).name,'mn')
            eval(['n=',wspace(sweepI).name,';']);
            eval(['m=',wspace(sweepI).name,';']);
        end
            
        [X,Y]=TwoChannelModel(n,m,d,[fx,fy],p(1:d),sigmadx2(1:d),sigmady2(1:d),...
            sigmafx2(1:fx),sigmafy2(1:fy),sigma2,M,mixing,noise,MAcoeff,ARcoeff); % Generate observations
        
        K=M; % Number of folds (M for LOOCV)
        rmax=min([12,floor(0.3*M),min(m,n)]); % Maximum number of correlated components (keep it low to reduce computation time)
        
        [~,~,Vx]=svd(X);
        [~,~,Vy]=svd(Y);
        
        dEstCVstd(ii,sim)=CrossValidation(X,Y,K,rmax,'standard');
        dEstCVpca(ii,sim)=CrossValidation(X,Y,K,rmax,'pca');
        dEstCV3set(ii,sim)=CrossValidation(X,Y,K,rmax,'3set');
        dEstCVxcov(ii,sim)=CrossValidation(X,Y,K,rmax,'xcov');
        dEstDet2(ii,sim)=maxmin(Vx,Vy,M,0,'real',[fx,fy],rmax,'mdl');
        dEstDet3(ii,sim)=IC_maxmin(Vx,Vy,M,'mdl','real',[fx,fy],rmax);

    end
    fprintf(repmat('\b',[1 length(s)]));
end
fprintf('Simulation %i/%i\n',nSim,nSim);

if sweepL>1
    figure, plot(sweepVar,[mean(dEstCVxcov==d,2),mean(dEstCV3set==d,2),...
        mean(dEstCVpca==d,2),mean(dEstCVstd==d,2),mean(dEstDet2==d,2),mean(dEstDet3==d,2)]);
    xlabel(wspace(sweepI).label); ylabel('probability of detection');
    grid on, box on;
    legend('XCOV-CV','3Set-CV','PCA-CV','CV','Detector 2 from [6]','Detector 3 from [6]','location','best');
else
    fprintf('Probability of detection for XCOV-CV: %1.2f\n',mean(dEstCVxcov==d,2));
    fprintf('Probability of detection for 3Set-CV: %1.2f\n',mean(dEstCV3set==d,2));
    fprintf('Probability of detection for PCA-CV: %1.2f\n',mean(dEstCVpca==d,2));
    fprintf('Probability of detection for CV: %1.2f\n',mean(dEstCVstd==d,2));
    fprintf('Probability of detection for Detector 2 from [6]: %1.2f\n',mean(dEstDet2==d,2));
    fprintf('Probability of detection for Detector 3 from [6]: %1.2f\n',mean(dEstDet3==d,2));
end
    
