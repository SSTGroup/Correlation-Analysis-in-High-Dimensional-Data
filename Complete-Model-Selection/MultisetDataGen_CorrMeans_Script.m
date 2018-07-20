%   File: MultisetDataGen_CorrMeans_Script.m
%   Copyright (c) <2017> <University of Paderborn>
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
% OVERVIEW: This script should provide an example of how to generate
%           multiple data sets with a prescribed correlation structure.
%           It enforces the transitive correlation condition.  This
%           condition is not necessary, but makes bookkeeping easier.
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
% ------------------------------------------------------------------------
% CREATED:      14/12/2017 by Tim Marrinan
%
% LAST EDITED:  17/01/2018 by Tim Marrinan
%               14/12/2017 by Tim Marrinan
%               15/12/2017 by Tim Marrinan
%
% NOTES: 
%
% ------------------------------------------------------------------------
clc; clear variables; %close all;
%% Scenario specification
scen='scen1';
switch lower(scen)
    case 'scen1' 
        n_sets          = 4;        % # of data sets
        signum          = 3;        % # of correlated + independent signals
                                    % per set
        tot_dims        = 4;        % # of sensors per set
        M               = 1000;     % # of samples per set
        full_corr       = 1;        % # of signals correlated across all
                                    % data sets (scalar)
        corr_across     = [3];    % across how many data sets should each 
                                    % additional signal be correlated?
                                    % (vector)
        corr_means      = [.9 .6]; % mean of the correlation
                                      % coefficients of each signal for all
                                      % data sets
        corr_std        = [.01 .01]; % std of the correlation coefficients
                                      % of each signal for all data sets
        RealComp        = 'real';   % real/complex data 
                                    % (only real is coded)
        sigmad          = 10;       % variance of correlated signals
        sigmaf          = 3;        % variance of independent signals
        sigmaN          = 1;       % variance of the noise
        mixing          = 'randn';  % mixing matrix type ('orth'/'randn')
        color           = 'white';  % noise type ('white'/'colored')
        MAcoeff         = 1;        % moving average coefficients for 
                                    % colored noise
        ARcoeff         = 1;        % auto-regressive coefficients for 
                                    % colored noise
        maxIters        = 99;       % maximum # of random draws allowed to
                                    % find a positive definite covariance
                                    % matrix
    case 'custom'
        % Use this scene to play around
        n_sets          = 4;
        signum          = 3;
        tot_dims        = 4;
        M               = 1000;
        full_corr       = 3;
        corr_across     = [];
        corr_means      = [.9 .8 .6];
        corr_std        = [.05 .05 .05];
        RealComp        = 'real';
        sigmad          = 10;
        sigmaf          = 3;       
        sigmaN          = 10;       
        mixing          = 'randn';
        color           = 'white';
        MAcoeff         = 1;
        ARcoeff         = 1;
        maxIters        = 99;
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
% This implements the method for correlation that requires correlation to 
% be transitive. It is not a necessary assumption, it just makes 
% bookkeeping easier. There are no diagonal correlations included here.
%
% p is a matrix containing the pairwise correlations. The rows are indexed
% in the order is the same as the order of indices in x_corrs, and the
% columns are indexed by the signal they represent.

[p,sigma_signals,~] = CorrelationStructureGen(n_sets,tot_corr,...
    corr_means,corr_std,signum,sigmad,sigmaf,maxIters);


%% Call the data generation function with this correlation structure
[X,R,A,S] = MultisetDataGen_CorrMeans(subspace_dims,signum,x_corrs,...
    mixing,sigmad,sigmaf,sigmaN,color,n_sets,p,sigma_signals,M,...
    MAcoeff,ARcoeff);
    
