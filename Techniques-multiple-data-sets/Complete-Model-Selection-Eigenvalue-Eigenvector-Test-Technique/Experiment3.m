%   File: Experiment3.m
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
% OVERVIEW: % MATLAB implementation of Experiment 3 in "Determining the Dimension and
% Structure of the Subspace Correlated Across Multiple Data Sets," [1].
% The experiment evaluate the performance of proposed 
% technique for model-order d (the number of correlated components)
% when the pairwise correlation coefficients do not
% exceed the threshold used for the proof of Theorem 1. The performance plot
% shows the mean accuracy of d as a function of the correlation coefficient \rho 
% (corresponding to the components in first data set) for 3 different SNRs.
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

n_sets          = 5;        % # of data sets
signum          = 7;        % # of correlated + independent signals
% per set
tot_dims        = 7;        % # of sensors per set
M               = 350;      % # of samples per set
num_iter = 1*1e1;           % # of trials for each data point
SNR_vec = [-2.5,0,2.5];     % SNR vector ranging from -10 to 15dB
rho_fix = 0.75;              % Correlation coefficients for all data
% sets except first
rho_vec = [.88,.8:-0.1:0.1]; % Correlation coefficient vector for
% first data set
full_corr       = 2;        % # of signals correlated across all data sets
corr_across     = [];       % across how many data sets should each
% additional signal be correlated?
RealComp        = 'real';   % real/complex data
% (only real is coded)
Distr = 'laplacian';        % gaussian or laplacian sources
sigmad          = 1;       % variance of correlated signals
sigmaf          = 1;        % variance of independent signals
mixing          = 'orth';  % mixing matrix type ('orth'/'randn')
color           = 'white';  % noise type ('white'/'colored')
MAcoeff         = [1];        % moving average coefficients for
% colored noise
ARcoeff         = [1];        % auto-regressive coefficients for
% colored noise
maxIters        = 99;       % maximum # of random draws allowed to
% find a positive definite covariance
% matrix
verbose         = 1;        % variable to print iteration counter 
% for trial and snr in the command window, choose 0 to not print 
        
%% Initialize some variables
x_corrs = combnk(1:n_sets,2);       % This is all the possible pairs of
% indices
subspace_dims = repmat(tot_dims,[1,n_sets]); % If you want the data sets
% to have different numbers of
% dimensions, you will have to change
% this snippet. Other things probably
% have to change too.
tot_corr = [repmat(n_sets,[1,full_corr]) corr_across];


%% Call the data generation function with this correlation structure
n_combs = size(x_corrs,1);

for snr =1:length(SNR_vec)
    
    sigmaN =sigmad/(10^(SNR_vec(snr)/10));
    if(verbose)
    display(['snr =' num2str(snr)]);
    end
    
    for r = 1:length(rho_vec)
        %% Define correlation structure
        
        % p is a matrix containing the pairwise correlations. The rows are indexed
        % in the order is the same as the order of indices in x_corrs, and the
        % columns are indexed by the signal they represent.
        
        rho = rho_vec(r);
        if(verbose)
        display(['rho =' num2str(rho)]);
        end
        p = zeros(n_combs,signum);
        p_rho = [rho_fix*ones(n_combs-n_sets+1,full_corr);rho*ones(n_sets-1,full_corr)]; % Varying pairwise corr coeffs for first data set
        p(:,1:full_corr) = p_rho;
        sigma_signals = sigmaf*ones(size(p));
        sigma_signals(p>0) = sigmad;
        
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
            
        end
        %% Functions with respect to rho
        % ----------------------------------------------------------------
        % Mean value of d_all
        d_all_func_mht(r) = mean(d_all_mht);
        d_all_func_bt(r) = mean(d_all_bt);
        
        % Mean Accuracy of d_all
        p_d_all_mht(r) = length(find(d_all_mht == full_corr))/num_iter;
        p_d_all_bt(r) = length(find(d_all_bt == full_corr))/num_iter;
        
        d_true = full_corr+ length(corr_across);
        
        % Mean Accuracy of d
        p_d_mht(r) = length(find(d_cap_mht == d_true))/num_iter;
        p_d_bt(r) = length(find(d_cap_bt == d_true))/num_iter;
        
        % Mean value of d
        d_func_bt(r) = mean(d_cap_bt);
        d_func_mht(r) = mean(d_cap_mht);
        
        % Precision and recall functions
        prec_func_mht(r) = mean(prec_vec_mht);
        recall_func_mht(r) = mean(recall_vec_mht);
        
        prec_func_bt(r) = mean(prec_vec_bt);
        recall_func_bt(r) = mean(recall_vec_bt);
        
    end
    
    %% Functions with respect to SNR
    
    p_d_mht_snr(:,snr) = p_d_mht;
    p_d_bt_snr(:,snr) = p_d_bt;
    mean_d_bt_snr(:,snr) =  d_func_bt;
    mean_d_all_bt_snr(:,snr) =  d_all_func_bt;
    mean_d_mht_snr(:,snr) =  d_func_mht;
    mean_d_all_mht_snr(:,snr) =  d_all_func_mht;
    p_d_all_mht_snr(:,snr) = p_d_all_mht;
    p_d_all_bt_snr(:,snr) = p_d_all_bt;
    
end

%% Plot performance

figure(); plot(rho_vec,p_d_bt_snr(:,1),'go:','markersize',12,'Linewidth',2); hold on;
plot(rho_vec,p_d_bt_snr(:,2),'m^-.','markersize',12,'Linewidth',2);
plot(rho_vec,p_d_bt_snr(:,3),'kd--','markersize',12,'Linewidth',2);
a = xlabel('Pairwise correlation coefficient $\rho$','fontsize',20,'FontName','Times New Roman');
b = ylabel('Mean accuracy of $\hat{d}$','fontsize',20,'FontName','Times New Roman');
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',20);
set(gca,'xdir','reverse');
xlim([0.1,1]);
x = [.64 .64];
y = [0 1];
line(x,y,'Color','red','LineStyle','--','Linewidth',1.5);
c = legend('SNR = -2.5dB','SNR = 0dB','SNR = 2.5dB','$\epsilon$ = 0.64');
set(c,'interpreter','latex');
title('Mean accuracy vs pairwise correlation coefficients');