%   File: Experiment12.m
%   Copyright (c) <2019> <University of Paderborn>
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
% OVERVIEW:  MATLAB implementation of Experiments 1 and 2 in "Determining the Dimension and
% Structure of the Subspace Correlated Across Multiple Data Sets," [1].
% The experiments evaluate the performance of proposed and competing
% techniques for model-order selection with special (experiment 1) and 
% arbitrary (experiment 2) correlation structures. The performance plots
% show mean accuracy and average mean value of d_all (number of components correlated across all 
% pairs of data sets) as a function of the SNR.
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
% [4] MCCA_CompleteModelSelection.m
%
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% [5] PCM_Bootstrap_Multiple_Datasets.m
%
%       Not included within.
%       Written by Tanuj Hasija, see file for documentation.
%
% [6] Max_Min_Multiple_Datasets.m
%
%       Not included within.
%       Written by Tanuj Hasija, see file for documentation.
%
% [7] ITC_MDL_Wu_Multiple_Datasets.m
%
%       Not included within.
%       Written by Tanuj Hasija, see file for documentation.
%
% [8] MCCA_KPD.m
%
%       Not included within.
%       Written by Suchita Bhinge (UMBC), see file for documentation.
%
%
% REFERENCES:
% Complete model selection using eigenvalue and eigenvectors of coherence
% matrix:
% [1]   T. Hasija, C. Lameiro, T. Marrinan and P. J. Schreier,"Determining
%       the Dimension and Structure of the Subspace Correlated Across
%       Multiple Data Sets," Submitted.
%
% MCCA-HT Complete model selection method:
% [2]   T. Marrinan, T. Hasija, C. Lameiro, and P. Schreier. "Complete
%       model selection in multiset canonical correlation analysis."
%       Proc. 26th European Signal Processing Conference (EUSIPCO), Rome, Italy, 2018.
%
% The maxmin method (used in [2]):
% [3]   Song, Yang, Peter J. Schreier, David Ramirez, and Tanuj Hasija.
%       "Canonical correlation analysis of high-dimensional data with very
%       small sample support." Signal Processing 128 (2016): 449-458.
%
% Product of coherence matrices (PCM) method:
% [4]   Hasija, Tanuj, Yang Song, Peter J. Schreier, and David Ramirez.
%       "Bootstrap-based detection of the number of signals correlated
%       across multiple data sets." In Signals, Systems and Computers,
%       2016 50th Asilomar Conference on, pp. 610-614. IEEE, 2016.
%
% The maxmin method for multiple data sets:
% [5]   T. Hasija, Y. Song, P. J. Schreier, and D. Ramirez,"Detecting the
%       dimension of the subspace correlated across multiple data sets in
%       the sample poor regime," in Proceedings of the IEEE Workshop on
%       Statistical Signal Processing, 2016.
%
% The ITC MDL method for multiple data sets:
% [6]   Y. Wu, K. W. Tam, and F. Li,"Determination of number of sources with
%       multiple arrays in correlated noise fields," IEEE Transactions on Signal
%       Processing, vol. 50, no. 6, pp. 1257–1260, 2002.
%
% MCCA-Knee point detection method:
% [7]   S. Bhinge, Y. Levin-Schwartz, and T. Adali,"Estimation of common
%       subspace order across multiple datasets: Application to multi-subject
%       fMRI data," in Proceedings of the 51st Annual Conference on Information
%       Sciences and Systems, 2017.
%
% ------------------------------------------------------------------------
% CREATED:      14/01/2019 by Tanuj Hasija
%
% LAST EDITED:  11/10/2019 by Tanuj Hasija
%
% NOTES:
%
% ------------------------------------------------------------------------
clc; clear variables; close all;
%% Scenario specification
scen='scen1'; % Select 'scen1' for Experiment 1, 'scen2' for Experiment 2 and 'custom' 
% for trying different scenarios 
switch lower(scen)
    case 'scen1'
        n_sets          = 4;        % # of data sets
        signum          = 7;        % # of correlated + independent signals
        % per set
        tot_dims        = 7;        % # of sensors per set
        M               = 350;      % # of samples per set
        num_iter = 1*1e1;           % # of trials for each data point
        SNR_vec = [-10:3:15];       % SNR vector ranging from -10 to 15dB
        full_corr       = 3;        % # of signals correlated across all data sets
        corr_across     = [];       % across how many data sets should each
        % additional signal be correlated?
        RealComp        = 'real';   % real/complex data
        % (only real is coded)
        Distr = 'gaussian';        % gaussian or laplacian sources
        sigmad          = 1;       % variance of correlated signals
        sigmaf          = 1;        % variance of independent signals
        mixing          = 'orth';  % mixing matrix type ('orth'/'randn')
        color           = 'white';  % noise type ('white'/'colored')
        MAcoeff         = [1];        % moving average coefficients for
        % colored noise
        ARcoeff         = [1];        % auto-regressive coefficients for
        % colored noise
        maxIters        = 99;       % maximum # of random draws allowed to
        % find a positive definite covariance matrix
        verbose         = 1;        % variable to print iteration counter 
        % for trial and snr in the command window, choose 0 to not print 
    case 'scen2'
        n_sets          = 4;        % # of data sets
        signum          = 7;        % # of correlated + independent signals
        % per set
        tot_dims        = 7;        % # of sensors per set
        M               = 350;      % # of samples per set
        num_iter = 1*1e1;           % # of trials for each data point
        SNR_vec = [-10:3:15];       % SNR vector ranging from -10 to 15dB
        full_corr       = 1;        % # of signals correlated across all data sets
        corr_across     = [3 2];       % across how many data sets should each
        % additional signal be correlated?
        RealComp        = 'real';   % real/complex data
        % (only real is coded)
        Distr = 'gaussian';        % gaussian or laplacian sources
        sigmad          = 1;       % variance of correlated signals
        sigmaf          = 1;        % variance of independent signals
        mixing          = 'orth';  % mixing matrix type ('orth'/'randn')
        color           = 'white';  % noise type ('white'/'colored')
        MAcoeff         = [1];        % moving average coefficients for
        % colored noise
        ARcoeff         = [1];        % auto-regressive coefficients for
        % colored noise
        maxIters        = 99;       % maximum # of random draws allowed to
        % find a positive definite covariance matrix
        verbose         = 1;        % variable to print iteration counter 
        % for trial and snr in the command window, choose 0 to not print 
    case 'custom'
        % Use this scene to play around
        n_sets          = 6;
        signum          = 6;
        tot_dims        = 6;
        M               = 1000;
        num_iter = 1*1e1;
        full_corr       = 1;
        corr_across     = [5 3];    % across how many data sets should each
        % additional signal be correlated?
        % (vector)
        corr_means      = [.8 .7 .6]; % mean of the correlation .7 .5 .6
        % coefficients of each signal for all
        % data sets
        corr_std        = [.1 .1 .1 ]; % std of the correlation coefficients
        % of each signal for all data sets
        RealComp        = 'real';
        Distr = 'laplacian';        % gaussian or laplacian sources
        sigmad          = 10;
        sigmaf          = 3;
        %sigmaN          = 10;
        SNR_vec = [-10:3:15];       % SNR vector ranging from -10 to 15dB
        mixing          = 'randn';
        color           = 'white';
        MAcoeff         = 1;
        ARcoeff         = 1;
        maxIters        = 99;
        verbose         = 1;        
    otherwise
        error('Unknown scenario');
end

%% Initialize some variables
x_corrs = combnk(1:n_sets,2);       % This is all the possible pairs of
% indices
subspace_dims = repmat(tot_dims,[1,n_sets]); % If you want the data sets
% to have different numbers of
% dimensions, you will have to change
% this snippet. Other things probably
% have to change too.
tot_corr = [repmat(n_sets,[1,full_corr]) corr_across];

%% Define correlation structure

% p is a matrix containing the pairwise correlations. The rows are indexed
% in the order is the same as the order of indices in x_corrs, and the
% columns are indexed by the signal they represent.

switch lower(scen)
    case 'scen1'
        p_c = [                     % Correlation coefficients of correlated components
            0.9186    0.9104    0.6234
            0.6442    0.8203    0.7156
            0.8199    0.7156    0.5701
            0.6939    0.7475    0.7206
            0.7890    0.6793    0.8181
            0.6328    0.6209    0.8490];
        p = [p_c,zeros(size(x_corrs,1),signum-length(tot_corr))];
        sigma_signals = sigmaf*ones(size(p));
        sigma_signals(p>0) = sigmad;
        
    case 'scen2'
        p_c = [                     % Correlation coefficients of correlated components
            0.9186    0.9104    0
            0.6442    0.8203    0.7156
            0.8199    0.7156    0
            0.6939    0    0
            0.7890    0    0
            0.6328    0    0];
        p = [p_c,zeros(size(x_corrs,1),signum-length(tot_corr))];
        sigma_signals = sigmaf*ones(size(p));
        sigma_signals(p>0) = sigmad;
        
    case 'custom'
        % This implements the method for generating the correlation coefficients
        
        [p,sigma_signals,~] = CorrelationStructureGen(n_sets,tot_corr,...
            corr_means,corr_std,signum,sigmad,sigmaf,maxIters);
end

%% Call the data generation function with this correlation structure
n_combs = size(x_corrs,1);
for snr =1:length(SNR_vec)
    
    sigmaN =sigmad/(10^(SNR_vec(snr)/10));
    if(verbose)
    display(['snr =' num2str(snr)]);
    end
    
    for iter=1:num_iter
        
        if(verbose)
        display(['iteration =' num2str(iter)]);
        end
        % This generates data.
        [X,R,A,S] = MultisetDataGen_CorrMeans(subspace_dims,signum,x_corrs,...
            mixing,sigmad,sigmaf,sigmaN,color,n_sets,p,sigma_signals,M,...
            MAcoeff,ARcoeff,Distr);
        
        CorrTruth = zeros(n_combs, tot_dims);
        CorrTruth(:,1:signum) = p>0;
        
        % Proposed technique
        Pfa_eval = 0.05; Pfa_evec= 0.05; % Probability of false alarm.
        B = 1000; % Number of bootstap resamplings to do.
        [d_cap,CorrEst_bt] = Eval_Evec_Tests_Bootstrap_Multiple_Datasets(X,Pfa_eval,Pfa_evec,B);
        d_cap_bt(iter) = d_cap;
        
        
        % Performance metrics
        % ----------------------------------------------------------------
        % Precision and recall.
        CorrEst_bt = CorrEst_bt.';
        true_positives_bt = sum(sum((CorrEst_bt+CorrTruth) == 2));
        false_positives_bt = sum(sum((CorrEst_bt-CorrTruth) == 1));
        false_negatives_bt = sum(sum((CorrEst_bt-CorrTruth) == -1));
        true_negatives_bt = sum(sum((CorrEst_bt+CorrTruth) == 0));
        precision_bt = true_positives_bt./(true_positives_bt+false_positives_bt);
        precision_bt(isnan(precision_bt)) = 0;
        recall_bt = true_positives_bt./(true_positives_bt+false_negatives_bt);
        recall_bt(isnan(recall_bt)) = 0;
        
        prec_vec_bt(iter) = precision_bt;
        recall_vec_bt(iter) = recall_bt;
        d_all_bt(iter)= sum(sum(CorrEst_bt,1) == n_combs,2);
        
        %% Comparions methods
        
        % Run the MCCA-HT complete model selection method [2].
        % ----------------------------------------------------------------
        ITC = 'ht';                 % Which detector to use for maxmin PCA-CCA.
        Pfa = 0.001;                % Probability of false alarm.
        [CorrEst_mht,~] = MCCA_CompleteModelSelection(X,ITC,Pfa);
        
        true_positives_mht = sum(sum((CorrEst_mht+CorrTruth) == 2));
        false_positives_mht = sum(sum((CorrEst_mht-CorrTruth) == 1));
        false_negatives_mht = sum(sum((CorrEst_mht-CorrTruth) == -1));
        true_negatives_mht = sum(sum((CorrEst_mht+CorrTruth) == 0));
        precision_mht = true_positives_mht./(true_positives_mht+false_positives_mht);
        precision_mht(isnan(precision_mht)) = 0;
        recall_mht = true_positives_mht./(true_positives_mht+false_negatives_mht);
        recall_mht(isnan(recall_mht)) = 0;
        
        prec_vec_mht(iter) = precision_mht;
        recall_vec_mht(iter) = recall_mht;
        
        d_cap_mht(iter) = sum(sum(CorrEst_mht,1) >0,2);
        d_all_mht(iter) = sum(sum(CorrEst_mht,1) == n_combs,2);
        
        % Run product of coherence matrices method [4].
        % ----------------------------------------------------------------
        Pfa_bt = 0.05;
        d_pcm(iter) = PCM_Bootstrap_Multiple_Datasets(X,Pfa_bt,B);
        
        % Run max-min multiple data sets [5].
        % ----------------------------------------------------------------
        rmax = min([floor(2*M/(3*n_sets)),tot_dims]);
        d_maxmin_mult(iter) = Max_Min_Multiple_Datasets(X,Pfa,rmax,'real');
        
        % Run ITC-MDL multiple data sets [6].
        % ----------------------------------------------------------------
        d_itc_mdl(iter) = ITC_MDL_Wu_Multiple_Datasets(X, 'real');
        
        % Run MCCA-KPD Method [7].
        % ----------------------------------------------------------------
        for i=1:n_sets
            X_kpd(:,:,i) = X{i};
        end
        d_kpd(iter) = MCCA_KPD(X_kpd);
        
    end
    %% Functions with respect to SNR
    % ----------------------------------------------------------------
    % Mean value of d_all
    d_all_func_itc_mdl(snr) = mean(d_itc_mdl);
    d_all_func_maxmin_mult(snr) = mean(d_maxmin_mult);
    d_all_func_pcm(snr) = mean(d_pcm);
    d_all_func_mht(snr) = mean(d_all_mht);
    d_all_func_bt(snr) = mean(d_all_bt);
    d_all_func_kpd(snr) = mean(d_kpd);
    
    % Mean Accuracy of d_all
    p_d_all_itc_mdl(snr) = length(find(d_itc_mdl == full_corr))/num_iter;
    p_d_all_maxmin_mult(snr) = length(find(d_maxmin_mult == full_corr))/num_iter;
    p_d_all_pcm(snr) = length(find(d_pcm == full_corr))/num_iter;
    p_d_all_mht(snr) = length(find(d_all_mht == full_corr))/num_iter;
    p_d_all_bt(snr) = length(find(d_all_bt == full_corr))/num_iter;
    p_d_all_kpd(snr) = length(find(d_kpd == full_corr))/num_iter;
    
    d_true = full_corr+ length(corr_across);
    
    p_d_mht(snr) = length(find(d_cap_mht == d_true))/num_iter;
    p_d_bt(snr) = length(find(d_cap_bt == d_true))/num_iter;
    
    d_func_bt(snr) = mean(d_cap_bt);
    d_func_mht(snr) = mean(d_cap_mht);
    
    prec_func_mht(snr) = mean(prec_vec_mht);
    recall_func_mht(snr) = mean(recall_vec_mht);
    
    prec_func_bt(snr) = mean(prec_vec_bt);
    recall_func_bt(snr) = mean(recall_vec_bt);
    
end

%% Plot performance

figure();
plot(SNR_vec,d_all_func_itc_mdl,'ro--','markersize',12,'Linewidth',2); hold on;
plot(SNR_vec,d_all_func_maxmin_mult,'bd--','markersize',12,'Linewidth',2);
plot(SNR_vec,d_all_func_pcm,'k*-','markersize',12,'Linewidth',2);
plot(SNR_vec,d_all_func_kpd,'c+-','markersize',12,'Linewidth',2);
plot(SNR_vec,d_all_func_mht,'m.--','markersize',12,'Linewidth',2); hold on;
plot(SNR_vec,d_all_func_bt,'gs-','markersize',12,'Linewidth',2);
a = xlabel('SNR (dB)','fontsize',20,'FontName','Times New Roman');
b = ylabel('Mean value of $\hat{d}_\textrm{all}$','fontsize',20,'FontName','Times New Roman');
c = legend('ITC MDL [16]','Max-Min [17]','[19]','mCCA-KPD [20]','mCCA-HT [23]','Proposed');
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(c,'interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',20);
title('Mean value vs SNR');

figure();
plot(SNR_vec,p_d_all_itc_mdl,'ro--','markersize',12,'Linewidth',2); hold on;
plot(SNR_vec,p_d_all_maxmin_mult,'bd--','markersize',12,'Linewidth',2);
plot(SNR_vec,p_d_all_pcm,'k*-','markersize',12,'Linewidth',2);
plot(SNR_vec,p_d_all_kpd,'c+-','markersize',12,'Linewidth',2);
plot(SNR_vec,p_d_all_mht,'m^-.','markersize',12,'Linewidth',2); hold on;
plot(SNR_vec,p_d_all_bt,'gs-','markersize',12,'Linewidth',2);
a = xlabel('SNR (dB)','fontsize',20,'FontName','Times New Roman');
b = ylabel('Mean accuracy of $\hat{d}_\textrm{all}$','fontsize',20,'FontName','Times New Roman');
c = legend('ITC MDL [16]','Max-Min [17]','[19]','mCCA-KPD [20]','mCCA-HT [23]','Proposed');
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(c,'interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',20);
title('Mean accuracy vs SNR');
