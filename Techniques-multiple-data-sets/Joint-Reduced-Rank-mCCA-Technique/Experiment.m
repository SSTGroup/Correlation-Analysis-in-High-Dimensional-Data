%   File: Experiment.m
%   Copyright (c) <2022> <University of Paderborn>
%   Signal and System Theory Group, Univ. of Paderborn, http://sst.upb.de
%   https://github.com/SSTGroup/Correlation-Analysis-in-High-Dimensional-Data
%
%   Permission is hereby granted, free of charge, to any person
%   obtaining a copy of this software and associated documentation
%   files (the "Software"), to deal in the Software without restriction,
%   including without limitation the rights to use, copy, modify and
%   merge the Software, subject to the following conditions:
%
%   1.) The Software is used for non-commercial research and
%       education purposes.
%
%   2.) The above copyright notice and this permission notice shall be
%       included in all copies or substantial portions of the Software.
%
%   3.) Publication, Distribution, Sublicensing, and/or Selling of
%       copies or parts of the Software requires special agreements
%       with the University of Paderborn and is in general not permitted.
%
%   4.) Modifications or contributions to the software must be
%       published under this license. The University of Paderborn
%       is granted the non-exclusive right to publish modifications
%       or contributions in future versions of the Software free of charge.
%
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
%   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%   OTHER DEALINGS IN THE SOFTWARE.
%
%   Persons using the Software are encouraged to notify the
%   Signal and System Theory Group at the University of Paderborn
%   about bugs. Please reference the Software in your publications
%   if it was used for them.
% ------------------------------------------------------------------------
% SYNTAX:   This is a script, not a function. Parameters are set within.
%
% OVERVIEW:  MATLAB implementation of Numerical example in "A GLRT for
% estimating the number of correlated components in sample-poor mCCA," [1].
% The experiments evaluate the performance of proposed and competing
% techniques for model-order selection. The performance plots
% show mean accuracy and average mean value of D (number of components
% correlated across multiple data sets) as a function of the SNR.
%
%
% DEPENDENCIES:
%
% [1] CorrelationStructureGen.m
%
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% [2] MultiSetDataGen_CorrMeans.m
%
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% [3] Eval_Evec_Tests_Bootstrap_Multiple_Datasets.m
%
%       Not included within.
%       Written by Tanuj Hasija, see file for documentation.
%
% [4] PairwisemCCAHT.m
%
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% [5] jRRmCCA.m
%
%       Not included within.
%       Written by Tanuj Hasija, see file for documentation.
%
% [6] Num_PCA_components.m
%
%       Not included within.
%       Written by Tanuj Hasija, see file for documentation.
%
%
% REFERENCES:
% [1]   T. Hasija and T. Marrinan ,"A GLRT for estimating the number of correlated
% components in sample-poor mCCA," Submitted.
%
% PairwisemCCA-HT:
% [2]   T. Marrinan, T. Hasija, C. Lameiro, and P. Schreier. "Complete
%       model selection in multiset canonical correlation analysis."
%       Proc. 26th European Signal Processing Conference (EUSIPCO), Rome, Italy, 2018.
%
% HybridjointEVD
% [3]   T. Hasija, T. Marrinan, C. Lameiro, and P. J. Schreier, “Determining the
% dimension and structure of the subspace correlated across multiple data
% sets,” Signal Processing, p. 107613, 2020.
% [4]   N. Asendorf and R. R. Nadakuditi, “Improving multiset canonical
% correlation analysis in high dimensional sample deficient settings,” in
% 2015 49th Asilomar Conference on Signals, Systems and Computers.
% IEEE, 2015
% [5] R. R. Nadakuditi and A. Edelman, “Sample eigenvalue based detection
% of high-dimensional signals in white noise using relatively few samples,”
% IEEE Transactions on Signal Processing, vol. 56, no. 7, pp. 2625–2638,
% 2008.
%
%
% ------------------------------------------------------------------------
% CREATED:      01/02/2022 by Tanuj Hasija
%
% LAST EDITED:  17/03/2022 by Tanuj Hasija
%
% NOTES:
%
% ------------------------------------------------------------------------

% Scenario specification

clc; clear variables; close all;
rng('default');
scen='scen1';
switch lower(scen)
    case 'scen1'
        n_sets          = 10;        % # of data sets
        signum          = 7;        % # of correlated + independent signals per set
        tot_dims        = 50;        % # of sensors per set
        M               = 200;     % # of samples per set
        num_iter = 1e1; % Originally for 5e2 iterations in the paper
        SNR_vec = [-6:3:18]; %-5:5:20 -6:3:18
        full_corr       = 1;        % # of signals correlated across all data sets (scalar)
        corr_across     = [ 9 8 7 6 ]; %9 7 5      % across how many data sets should each
        % additional signal be correlated?
        corr_means      = [.85 .8 .75 .7 .6]; % mean of the correlation
        % coefficients of each signal for all
        % data sets
        corr_std        = [.01 .01 .01 .01 .01  ]*0; % std of the correlation coefficients
        % of each signal for all data sets
        RealComp        = 'real';   % real/complex data
        % (only real is coded for now)
        sigmad          = 1;       % variance of correlated signals
        sigmaf          = 2;        % variance of independent signals
        add_weak_uncorr_sig = 1;      % 1 if weak uncorrelated signals are added after signum components, 0 otherwise
        sigmafweak      = sigmad/2;  % variance of weak uncorrelated signals
        num_signumweak  = 5;         % number of weak uncorrelated signals
        mixing          = 'orth';  % mixing matrix type ('orth'/'randn')
        color           = 'white';  % noise type ('white'/'colored')
        MAcoeff         = [1];        % moving average coefficients for
        % colored noise
        ARcoeff         = [1];        % auto-regressive coefficients for
        % colored noise
        maxIters        = 99;       % maximum # of random draws allowed to
        % find a positive definite covariance
        % matrix
    case 'custom'
        % Use this scene to play around
        n_sets          = 10;
        signum          = 10;
        tot_dims        = 30;
        M               = 200;
        SNR_vec         = [-10:5:10];
        num_iter        = 1e0;
        full_corr       = 0;
        corr_across     = [ 9 8 6];
        corr_means      = [.9 .8 .6];
        corr_std        = [.01 .01 .01];
        RealComp        = 'real';
        sigmad          = 2;
        sigmaf          = 5;
        add_weak_uncorr_sig= 0;
        sigmafweak      = sigmad/2;  
        num_signumweak  = 1;       
        mixing          = 'randn';
        color           = 'white';
        MAcoeff         = 1;
        ARcoeff         = 1;
        maxIters        = 99;
    otherwise   
        error('Unknown scenario');
end

% Initialize some variables
x_corrs = combnk(1:n_sets,2);       % This is all the possible pairs of
% indices
subspace_dims = repmat(tot_dims,[1,n_sets]); % If you want the data sets
% to have different numbers of
% dimensions, you will have to change
% this snippet. Other things probably
% have to change too.
tot_corr = [repmat(n_sets,[1,full_corr]) corr_across];

% Define correlation structure
% This implements the method for correlation that requires correlation to
% be transitive. It is not a necessary assumption, it just makes
% bookkeeping easier. There are no diagonal correlations included here.
%
% p is a matrix containing the pairwise correlations. The rows are indexed
% in the order is the same as the order of indices in x_corrs, and the
% columns are indexed by the signal they represent.

[p,sigma_signals,~] = CorrelationStructureGen(n_sets,tot_corr,...
    corr_means,corr_std,signum,sigmad,sigmaf,maxIters);

% Call the data generation function with this correlation structure
n_combs = size(x_corrs,1);

%% Generate thresholds for empirical distribution for jRR-mCCA
pfa = 0.01;
emp_dist_iters = 5*1e2;
rmax = 15; %should be smaller than M/n_sets

% Check if results already exist
file_name = ['Threshold_sam_poor_','n_sets',num2str(n_sets),'M',num2str(M),'rmax',...
    num2str(rmax),'pfa',num2str(pfa),'MCtrials',num2str(emp_dist_iters),'.mat'];
if isfile(file_name)
    disp("Thresholds for distribution of statistic found");
    load(file_name,'tau_rc');
else
    disp("Thresholds for distribution of statistic NOT found, generating them");
    tau_rc = zeros(rmax);
    for r=1:rmax
        display(['r=' num2str(r)]);
        parfor s=0:r-1 % parfor makes the computations faster
            %             display(['s=' num2str(s)]);
            idx=s+1;
            [tau_rc(r,s+1),~,~,~] =  GenerateEmpiricalDistMultipleDsets_RandCorrStruc(n_sets,r,M,s,pfa,emp_dist_iters) ;
        end
    end
    save(filename);
end

%% MC trials for different SNR values
for snr =1:length(SNR_vec)
    sigmaN =sigmad/(10^(SNR_vec(snr)/10));
    display(['snr =' num2str(snr)]);
    for iter=1:num_iter
        %display(['iteration =' num2str(iter)]);
        % This generates data. Modify this file to change scenario parameters.
        [X,R,A,S] = MultisetDataGen_CorrMeans(subspace_dims,signum,x_corrs,...
            mixing,sigmad,sigmaf,sigmaN,color,n_sets,p,sigma_signals,M,...
            MAcoeff,ARcoeff,add_weak_uncorr_sig,sigmafweak,num_signumweak);
        
        % Other parameters
        ITC = 'ht';                 % Which detector to use for maxmin PCA-CCA used in pairwise mCCA-HT.
        pfa_maxmin = 0.005;                % Probability of false alarm for maxmin PCA-CCA
        
        % Run the pairwise mCCA-HT
        % ----------------------------------------------------------------
        [CorrEst,~] = MCCA_CompleteModelSelection(X,ITC,pfa_maxmin);
        d_est_pmht(iter) = sum(sum(CorrEst,1) >0,2);
        
        % Run the hybrid jointEVD
        % ----------------------------------------------------------------
        % PCA using [5] 
        for l=1:n_sets
            Rxx = X{l}*X{l}.'/M;
            L = eig(Rxx);
            L = sort(L,'descend');
            [numcomp(l),~] = Num_PCA_components(L,subspace_dims(l),M,1);
        end
        r_hjevd(iter) = round(mean(numcomp));
        
        if(r_hjevd(iter))
            X_pca = cell(n_sets,1);
            for l=1:n_sets
                [U,~,~] = svd(X{l},'econ');
                X_pca{l} = U(:,1:r_hjevd(iter)).'*X{l};
            end
            % jointEVD 
            B = 500; % Number of bootstrap resamples
            [d_est_hjevd(iter),~] = Eval_Evec_Tests_Bootstrap_Multiple_Datasets(X_pca,pfa,pfa,B);
        else
            d_est_hjevd(iter) = 0;
        end
        
        % Run the jRR-mCCA 
        % ----------------------------------------------------------------
        [d_est_jm(iter),r_est_jm(iter)] = jRRmCCA(X,rmax,tau_rc);
    end
    d_true = full_corr+ length(corr_across);
    
    p_d_pmht(snr) = length(find(d_est_pmht == d_true))/num_iter;
    p_d_hjevd(snr) = length(find(d_est_hjevd == d_true))/num_iter;
    p_d_jm(snr) = length(find(d_est_jm == d_true))/num_iter;
    
    d_cap_func_hjevd(snr) = mean(d_est_hjevd);
    d_cap_func_pmht(snr) = mean(d_est_pmht);
    d_cap_func_jm(snr) = mean(d_est_jm);

    r_cap_func_jm(snr) = mean(r_est_jm);
    r_cap_func_hjevd(snr) = mean(r_hjevd);

    
end

%% Figures

figure();
plot(SNR_vec,d_cap_func_pmht,'m^:','markersize',12,'Linewidth',2); hold on;
plot(SNR_vec,d_cap_func_hjevd,'gs-','markersize',12,'Linewidth',2);
plot(SNR_vec,d_cap_func_jm,'ko-','markersize',12,'Linewidth',2);
plot(SNR_vec,ones(1,length(SNR_vec))*d_true,'r-','markersize',12,'Linewidth',2);
a = xlabel('SNR (dB)','fontsize',20,'FontName','Times New Roman');
b = ylabel('Mean $\hat{D}$','fontsize',20,'FontName','Times New Roman');
c = legend('Pairwise mCCA-HT','Hybrid jointEVD','Proposed jRR-mCCA','True value');
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(c,'interpreter','latex');
xlim([-6 18]); %ylim([0 8]);
set(gca,'FontName','Times New Roman','FontSize',24);
% saveas(gcf,'mean_acc_d'); saveas(gcf,'mean_acc_d','epsc');

figure();
plot(SNR_vec,p_d_pmht,'m^:','markersize',12,'Linewidth',2); hold on;
plot(SNR_vec,p_d_hjevd,'gs-','markersize',12,'Linewidth',2);
plot(SNR_vec,p_d_jm,'ko-','markersize',12,'Linewidth',2);
a = xlabel('SNR (dB)','fontsize',20,'FontName','Times New Roman');
b = ylabel('Mean accuracy of $\hat{D}$','fontsize',20,'FontName','Times New Roman');
c = legend('Pairwise mCCA-HT','Hybrid jointEVD','Proposed jRR-mCCA');
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(c,'interpreter','latex');
xlim([-6 18]);
set(gca,'FontName','Times New Roman','FontSize',24);


figure(); hold on;
plot(SNR_vec,r_cap_func_hjevd,'gs-','markersize',12,'Linewidth',2);
plot(SNR_vec,r_cap_func_jm,'ko-','markersize',12,'Linewidth',2);
a = xlabel('SNR (dB)','fontsize',20,'FontName','Times New Roman');
b = ylabel('Mean PCA rank','fontsize',20,'FontName','Times New Roman');
c = legend('Hybrid jointEVD','Proposed jRR-mCCA');
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(c,'interpreter','latex');
xlim([-6 18]);
set(gca,'FontName','Times New Roman','FontSize',24);
