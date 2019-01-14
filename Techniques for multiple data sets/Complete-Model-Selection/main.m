%   File: main.m
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
% SYNTAX:   
%
% This is a script, not a function.  Parameters are set within.
%
% OVERVIEW: 
%
% MATLAB implementation of "Complete model selection in multiset canonical 
% correlation analysis" [1].                 
%
% DEPENDENCIES:
%
% EUSIPCO_ScenarioSpecification.m
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% MutlisetDataGen_CorrMeans_Script.m
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% CorrelationStructureGen.m
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% MultiSetDataGen_CorrMeans.m
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% MCCA_CompleteModelSelection.m 
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% ProdOfCoherence_WithBootstrap.m
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% IMCCA.m
%       Not included within.
%       Written by Tim Marrinan, see file for documentation.
%
% REFERENCES:
%
% Complete model selection method:
% [1]   T. Marrinan, T. Hasija, C. Lameiro, and P. Schreier. "Complete 
%       model selection in multiset canonical correlation analysis." 
%       (In preparation.)
%
% The maxmin method (used in [1]):
% [2]   Song, Yang, Peter J. Schreier, David Ramírez, and Tanuj Hasija. 
%       "Canonical correlation analysis of high-dimensional data with very 
%       small sample support." Signal Processing 128 (2016): 449-458.
%
% Product of coherence matrices method:
% [3]   Hasija, Tanuj, Yang Song, Peter J. Schreier, and David Ramírez. 
%       "Bootstrap-based detection of the number of signals correlated 
%       across multiple data sets." In Signals, Systems and Computers, 
%       2016 50th Asilomar Conference on, pp. 610-614. IEEE, 2016.
%
% IMCCA method:
% [4]   Asendorf, Nicholas, and Raj Rao Nadakuditi. "Improving multiset 
%       canonical correlation analysis in high dimensional sample 
%       deficient settings." In Signals, Systems and Computers, 2015 
%       49th Asilomar Conference on, pp. 112-116. IEEE, 2015. Harvard
%
% ------------------------------------------------------------------------
% CREATED:      21/02/2018 by Tim Marrinan
%
% LAST EDITED:  21/02/2018 by Tim Marrinan
%               23/02/2018 by Tim Marrinan
%
% NOTES:
%
% 
% ------------------------------------------------------------------------

clear variables;close all;clc;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)))

scenario{1}     = 'corrstruct1';
scenario{2}     = 'corrstruct2';
scenario{3}     = 'custom';
results         = cell(size(scenario,2),1);
num_iterations  = 1*1e3;            % number of monte carlo trials
printout        = 'verbose';        % 'verbose'/'terse'/'none'
plots           = 'on';             % 'on'/'off'
queued          = 1:2;              % the indices of each scenario to run
finished        = [];

for current = queued

    % Scenario specification.
    % Modify the 'custom' scenario to exeriment.
    % Change 'printout' to show or hide output.
    % ----------------------------------------------------------------
    scen = scenario{current};
    run('EUSIPCO_ScenarioSpecification.m') 
    
    switch lower(printout)
        case 'verbose'
            fprintf('\t--------------------------------------------------------'); 
                fprintf('\n\tSimulating "%s" scenario:\n',scenario{current});
                fprintf('\t--------------------------------------------------------');
        case 'terse'
        otherwise
    end
    
    results{current}.precision  = zeros(num_iterations,size(sigmaN,2));
    results{current}.recall     = zeros(num_iterations,size(sigmaN,2));
    results{current}.corrEst    = cell(num_iterations,size(sigmaN,2));
    results{current}.compEst    = zeros(num_iterations,size(sigmaN,2));
    results{current}.totEst     = zeros(num_iterations,size(sigmaN,2));
    results{current}.pocEst     = zeros(num_iterations,size(sigmaN,2));
    results{current}.imccaEst   = zeros(num_iterations,size(sigmaN,2));
    
    % Create ground truth
    n_combs                 = size(x_corrs,1);
    corrTruth               = zeros(n_combs, tot_dims);
    corrTruth(:,1:signum)   = p>0;
    compTruth               = sum(sum(corrTruth,1)==size(x_corrs,1));
    totTruth                = sum(sum(corrTruth,1)>0);
    
    for iVar = 1 : size(sigmaN,2)   % Test with different levels of noise
        
        switch lower(printout)
            case 'verbose'
                    fprintf('\n\tExperiment:\t\t\t\t\t%d/%d\n',...
                        iVar,size(sigmaN,2));
                fprintf('\tNoise variance:\t\t\t\t\t%.4f\n',sigmaN(iVar));
            case 'terse'
            otherwise
        end
        
        % Generate empirical distribution for IMCCA [4].
        % ----------------------------------------------------------------
        [empEvals] =  IMCCA_GenEmpiricalDistribution(n_sets,M,tot_dims,...
            sigmaN(iVar),emp_dist_iters);
    
        for iter = 1 : num_iterations       % Run 'num_iterations' trials
        
            % Generate data
            % ------------------------------------------------------------
            [X,R,A,S] = MultisetDataGen_CorrMeans(subspace_dims,signum,...
                x_corrs,mixing,sigmad,sigmaf,sigmaN(iVar),color,n_sets,...
                p,sigma_signals,M,MAcoeff,ARcoeff);

            % All methods start by computing an orthonormal basis for the 
            % data. Precompute this to save time.
            % ------------------------------------------------------------
            fullU = cell(n_sets,1);
            origSigma = cell(n_sets,1);
            for i = 1 : n_sets
                [fullU{i},origSigma{i},~] = svd(X{i}',0);
            end 

            % Run model (order) selection
            % ------------------------------------------------------------

            % Proposed method [1].
            [corrEst,~] = MCCA_CompleteModelSelection(X,ITC,Pfa,fullU);

            % Product of coherence matrices method [3].
            pocEst = ProdOfCoherence_WithBootstrap(X,Pfa,boot_iters,fullU); 

            % IMCCA method [4].
            imccaEst = IMCCA(X,empEvals,k_hat_percent,Pfa,fullU,...
                origSigma);

            % Performance metrics
            % ------------------------------------------------------------

            % Precision and recall.
            true_positives = sum(sum((corrEst+corrTruth) == 2));
            false_positives = sum(sum((corrEst-corrTruth) == 1));
            false_negatives = sum(sum((corrEst-corrTruth) == -1));
            true_negatives = sum(sum((corrEst+corrTruth) == 0));
            precision = true_positives./(true_positives+false_positives);
            precision(isnan(precision)) = 0;
            recall = true_positives./(true_positives+false_negatives);
            recall(isnan(recall)) = 0;

            % Product of coherence matrices comparison
            compEst = sum(sum(corrEst,1) == n_combs,2);

            % IMCCA comparison
            totEst = sum(sum(corrEst,1) >0,2);
            
            % Store results
            results{current}.name                   = scenario{current};   
            results{current}.corrTruth              = corrTruth;
            results{current}.compTruth              = compTruth;
            results{current}.totTruth               = totTruth;
            results{current}.corrEst{iter,iVar}     = corrEst;              
            results{current}.precision(iter,iVar)   = precision;
            results{current}.recall(iter,iVar)      = recall;
            results{current}.compEst(iter,iVar)     = compEst;
            results{current}.totEst(iter,iVar)      = totEst;
            results{current}.pocEst(iter,iVar)      = pocEst;
            results{current}.imccaEst(iter,iVar)    = imccaEst;
            results{current}.SNR                    = 10*log10(...
                sigmad./sigmaN);
            
            
            % Print iteration output
            switch lower(printout)
                case 'verbose'
                    %fprintf('\n\tSimulation of "%s" scenario (%d/%d)\n',...
                    %scenario{current},iVar,size(sigmaN,2));
                    %fprintf('\t-----------------------------------------------');    
                    %fprintf('\n\tNoise variance: %f\n',sigmaN(iVar));
                    %fprintf('\tIteration:\t\t\t\t\t%d/%d\n',iter,num_iterations);
                case 'terse'
                    clc
                    fprintf('\n\tScenario:\t"%s"\n\tExperiment:\t%d/%d\n\tIteration:\t%d/%d\n',...
                        scenario{current},iVar,size(sigmaN,2),iter,num_iterations);
                    otherwise
            end
        end
        
        % Print average performance
        switch lower(printout)
            case 'verbose'
                fprintf('\n\tAvg. Precision:\t\t\t\t\t%.3f\n\tAvg. Recall:\t\t\t\t\t%.3f\n',...
                    mean(results{current}.precision(:,iVar)),...
                    mean(results{current}.recall(:,iVar)));

%                 fprintf('\n\tSignals correlated across all sets: %d\n',...
%                     results{current}.compTruth);
                fprintf('\n\tAvg. accuracy of the proposed method:\t\t%.3f\n',...
                    mean(results{current}.compEst(:,iVar)==...
                    results{current}.compTruth));
                fprintf('\tAvg. accuracy of product of coherence matrices:\t%.3f\n',...
                    mean(results{current}.pocEst(:,iVar)==...
                    results{current}.compTruth));

%                 fprintf('\n\tTotal correlated signals: %d\n',...
%                     results{current}.totTruth);
                fprintf('\n\tAvg. accuracy of the proposed method:\t\t%.3f\n',...
                    mean(results{current}.totEst(:,iVar)==...
                    results{current}.totTruth));
                fprintf('\tAvg. accuracy of IMCCA:\t\t\t\t%.3f\n',...
                    mean(results{current}.imccaEst(:,iVar)==...
                    results{current}.totTruth));
                fprintf('\t--------------------------------------------------------'); 
            case 'terse'

            otherwise
        end
    end
    switch lower(printout)
        case 'verbose'
            fprintf('\n\t%d iterations run. "%s" scenario complete.\n',...
                num_iterations,scenario{current});
            fprintf('\t--------------------------------------------------------\n\n'); 
        case 'terse'
            fprintf('\tscenario complete.\n');
        otherwise
    end
end

switch lower(printout)
    case 'verbose'
        fprintf('\n\t%d scenarios simulated.\n\tAll done.\n', size(queued,2))
    case 'terse'
        clc
        fprintf('\t%d scenarios simulated.\n\tAll done.\n', size(queued,2))
    otherwise
end
                



%% Generate plots
if strcmp(lower(plots),'on')
    % Color and line-style denote method
    MethodColors    = cell(6,1);
    MethodColors{1} = [152,78,163]/255;
    MethodColors{2} = [179,222,105]/255;
    MethodColors{3} = [251,128,114]/255;
    MethodColors{4} = [253,180,98]/255;
    MethodColors{5} = [141,211,199]/255;
    MethodColors{6} = [228,26,28]/255;
    
    % Markers denote scenario
    SceneMarkers = cell(3,1);
    SceneMarkers{1} = 'o';
    SceneMarkers{2} = 'x';
    SceneMarkers{3} = 's';
    
    f_size = 16;
    m_size = 10;
    l_width = 3;
    fig_width = 12;
    fig_height = 12;
    
    f       = cell(3,1);
    f{1}    = figure('units','centimeters','Position',...
                [10,10,fig_width,fig_height]);
    f{2}    = figure('units','centimeters','Position',...
                [10,10,fig_width,fig_height]);
    f{3}    = figure('units','centimeters','Position',...
                [10,10,fig_width,fig_height]);    
    g       = cell(3,4,size(queued,2));
    
    for current = queued
        compMean    = zeros(size(sigmaN,2),1);
        PoCMean     = zeros(size(sigmaN,2),1);
        totMean     = zeros(size(sigmaN,2),1);
        imccaMean   = zeros(size(sigmaN,2),1);
        precMean    = zeros(size(sigmaN,2),1);
        recMean     = zeros(size(sigmaN,2),1);
        
        for iVar = 1 : size(sigmaN,2)
            compMean(iVar) = mean(results{current}.compEst(:,iVar)==...
                results{current}.compTruth);
            PoCMean(iVar) = mean(results{current}.pocEst(:,iVar)==...
                    results{current}.compTruth);
                
            totMean(iVar) = mean(results{current}.totEst(:,iVar)==...
                results{current}.totTruth);
            imccaMean(iVar) = mean(results{current}.imccaEst(:,iVar)==...
                    results{current}.totTruth);
            
            precMean(iVar) = mean(results{current}.precision(:,iVar));
            recMean(iVar) = mean(results{current}.recall(:,iVar));    
        end
        
        % CVIP vs. Product of Coherence Matrices
        set(0,'CurrentFigure',f{1})
        g{1,1,current} = plot(results{current}.SNR,compMean,'-',...
            'Color',MethodColors{1},'LineWidth',l_width);
        hold on
        g{1,2,current} = plot(results{current}.SNR,PoCMean,'--',...
            'Color', MethodColors{2},'LineWidth',l_width);
        g{1,3,current} = plot(results{current}.SNR,compMean,...
            SceneMarkers{current},'Color','black','MarkerSize',m_size);
        g{1,4,current} = plot(results{current}.SNR,PoCMean,...
            SceneMarkers{current},'Color','black','MarkerSize',m_size);
        axis([-10 15 0 1])
        
        
        % CVIP vs. IMCCA
        set(0,'CurrentFigure',f{2})
        g{2,1,current} = plot(results{current}.SNR,totMean,'-',...
            'Color',MethodColors{1},'LineWidth',l_width);
        hold on
        g{2,2,current} = plot(results{current}.SNR,imccaMean,'--',...
            'Color',MethodColors{3},'LineWidth',l_width);
        g{2,3,current} = plot(results{current}.SNR,totMean,...
            SceneMarkers{current},'Color','black','MarkerSize',m_size);
        g{2,4,current} = plot(results{current}.SNR,imccaMean,...
            SceneMarkers{current},'Color','black','MarkerSize',m_size);
        axis([-10 15 0 1])
        
        % Precision/Recall
        set(0,'CurrentFigure',f{3})
        g{3,1,current} = plot(results{current}.SNR,precMean,'-',...
            'Color',MethodColors{5},'LineWidth',l_width);
        hold on
        g{3,2,current} = plot(results{current}.SNR, recMean,'--',...
            'Color',MethodColors{6},'LineWidth',l_width);
        g{3,3,current} = plot(results{current}.SNR,precMean,...
            SceneMarkers{current},'Color','black','MarkerSize',m_size);
        g{3,4,current} = plot(results{current}.SNR,recMean,...
            SceneMarkers{current},'Color','black','MarkerSize',m_size);
        axis([-10 15 0 1])
    end

    % Currently this plotting stuff only seems to work if the SNR vectors
    % are the same for everything, but I don't know why
    t = min(queued);
    L{1,1} = [g{1,1,t}, g{1,2,t}];
    L{1,2} = [{'Proposed method'}, {'Prod.of Coherence'}];
    L{2,1} = [g{2,1,t}, g{2,2,t}];
    L{2,2} = [{'Proposed method'}, {'IMCCA'}];
    L{3,1} = [g{3,1,t}, g{3,2,t}];
    L{3,2} = [{'Precision'}, {'Recall'}];
    for i = 1 : 3
        for current = queued
            L{i,1} = [L{i,1}, g{i,3,current}];
            L{i,2} = [L{i,2}, {results{current}.name}];
        end
        set(0,'CurrentFigure',f{i})
        legend(L{i,1}, L{i,2},'Location','SouthEast');
        xlabel('SNR per component (dB)');
        ylabel('Mean accuracy');
    end
    
end