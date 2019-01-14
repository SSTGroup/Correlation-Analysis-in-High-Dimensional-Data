%   File: EUSIPCO_ScenarioSpecification.m
%   Copyright (c) <2018> <University of Paderborn>
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
% OVERVIEW: 
%
% Scenario specifications associated with the experiments from [1], as well
% as a custom scenario for running user-defined experiments.
%
% INPUTS:
%
% OUTPUTS:    
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
%
% REFERENCES:
%
% [1]   T. Marrinan, T. Hasija, C. Lameiro, and P. Schreier. "Complete 
%       model selection in multiset canonical correlation analysis." 
%       (In preparation.)
% ------------------------------------------------------------------------
% CREATED:      21/02/2018 by Tim Marrinan
%
% LAST EDITED:  21/02/2018 by Tim Marrinan
%
% NOTES: 
%
% ------------------------------------------------------------------------


%% Scenario specification
% ------------------------------------------------------------------------
switch lower(scen)
    case 'corrstruct1' 
        n_sets          = 4;        % # of data sets
        signum          = 5;        % # of signals per set
        tot_dims        = 5;        % # of sensors per set
        M               = 500;      % # of samples per set
        RealComp        = 'real';   % 'real'/'complex' (only real is coded)
        sigmad          = 10;       % variance of correlated signals
        sigmaf          = 3;        % variance of independent signals
        sigmaN          = [0.3162 1 3.1623 10 31.6228 100];  % variance of the noise
        mixing          = 'orth';   % mixing matrix type ('orth'/'randn')
        color           = 'white';  % noise type ('white'/'colored')
        MAcoeff         = 0;        % moving average coefficients for 
                                    % colored noise
        ARcoeff         = 0;        % auto-regressive coefficients for 
                                    % colored noise
        Pfa             = 1e-3;     % Probability of false alarm.
        ITC             = 'ht';     % Detector used for maxmin PCA-CCA [2],
                                    % ('ht'/'mdl'/'aic').
        boot_iters      = 500;      % Number of bootstap resamplings to do
                                    % for comparison method (2).
        k_hat_percent   = .90;      % Method (3) requires us to keep some
                                    % number of singular vectors from each
                                    % subspace. This tells us what
                                    % percentage of the variance to retain
                                    % in the kept singular vectors.
        emp_dist_iters  = 1*1e3;    % Samples to generate for the IMCCA
                                    % empirical distribution
        % p is analogous to Table 1a in [1], except it is transposed.                            
        p               = [ 0.9 0.9 0.9 0.0 0.0 0.0;
                            0.0 0.8 0.0 0.8 0.0 0.8;
                            0.0 0.7 0.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 0.0]';
    case 'corrstruct2'
        n_sets          = 4;        % # of data sets
        signum          = 5;        % # of signals per set
        tot_dims        = 5;        % # of sensors per set
        M               = 500;      % # of samples per set
        RealComp        = 'real';   % 'real'/'complex' (only real is coded)
        sigmad          = 10;       % variance of correlated signals
        sigmaf          = 3;        % variance of independent signals
        sigmaN          = [0.3162 1 3.1623 10 31.6228 100];  % variance of the noise
        mixing          = 'orth';   % mixing matrix type ('orth'/'randn')
        color           = 'white';  % noise type ('white'/'colored')
        MAcoeff         = 0;        % moving average coefficients for 
                                    % colored noise
        ARcoeff         = 0;        % auto-regressive coefficients for 
                                    % colored noise
        Pfa             = 1e-3;     % Probability of false alarm.
        ITC             = 'ht';     % Detector used for maxmin PCA-CCA [2],
                                    % ('ht'/'mdl'/'aic').
        boot_iters      = 500;      % Number of bootstap resamplings to do
                                    % for comparison method (2).
        k_hat_percent   = .90;      % Method (3) requires us to keep some
                                    % number of singular vectors from each
                                    % subspace. This tells us what
                                    % percentage of the variance to retain
                                    % in the kept singular vectors.
        emp_dist_iters  = 1*1e3;    % Samples to generate for the IMCCA
                                    % empirical distribution
        % p is analogous to Table 1b in [1], except it is transposed.
        p               = [ 0.9 0.9 0.9 0.9 0.9 0.9;
                            0.0 0.0 0.0 0.0 0.8 0.0;
                            0.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 0.0]';
                                   
    case 'custom'
        n_sets          = 4;        % # of data sets
        signum          = 3;        % # of correlated + independent signals
                                    % per set
        tot_dims        = 10;        % # of sensors per set
        M               = 500;      % # of samples per set
        RealComp        = 'real';   % real/complex data 
                                    % (only real is coded)
        sigmad          = 10;       % variance of correlated signals
        sigmaf          = 3;        % variance of independent signals
        sigmaN          = [0.3162 1 3.1623 10 31.6228 100]; % variance of the noise
        mixing          = 'orth';   % mixing matrix type ('orth'/'randn')
        color           = 'white';  % noise type ('white'/'colored')
        MAcoeff         = 0;        % moving average coefficients for 
                                    % colored noise
        ARcoeff         = 0;        % auto-regressive coefficients for 
                                    % colored noise
        maxIters        = 99;       % maximum # of random draws allowed to
                                    % find a positive definite covariance
                                    % matrix
        Pfa             = 1e-3;     % Probability of false alarm.
        ITC             = 'ht';     % Detector used for maxmin PCA-CCA [2],
                                    % ('ht'/'mdl'/'aic').
        boot_iters      = 500;      % Number of bootstap resamplings to do
                                    % for comparison method (2).
        k_hat_percent   = .90;      % Method (3) requires us to keep some
                                    % number of singular vectors from each
                                    % subspace. This tells us what
                                    % percentage of the variance to retain
                                    % in the kept singular vectors.
        emp_dist_iters  = 1*1e3;    % Samples to generate for the IMCCA
                                    % empirical distribution
        % p is a randomly generated table of correlations analagous to 
        % those in Table 1 of [1]. It requires the following additional
        % parameters to be set.
        full_corr       = 1;        % # of signals correlated across all
                                    % data sets (scalar)
        corr_across     = [3,2];    % across how many data sets should each 
                                    % additional signal be correlated?
                                    % (vector)
        corr_means      = [.9 .8 .7];  % mean of the correlation
                                    % coefficients of each signal for all
                                    % data sets
        corr_std        = [0 0 0];  % std of the correlation coefficients
                                    % of each signal for all data sets
        tot_corr = [repmat(n_sets,[1,full_corr]) corr_across];
        [p,~,~] = CorrelationStructureGen(n_sets,tot_corr,corr_means,...
            corr_std,signum,sigmad,sigmaf,maxIters);                            
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
sigma_signals = sigmaf*ones(size(p));
sigma_signals(p>0) = sigmad;
