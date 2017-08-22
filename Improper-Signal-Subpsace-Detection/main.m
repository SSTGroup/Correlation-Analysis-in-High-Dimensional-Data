% This is a MATLAB implementation supporting the paper
%
% "Determining the Dimension of the Improper Signal
% Subspace in Complex-Valued Data" by Tanuj Hasija, Christian Lameiro and
% Peter J. Schreier, Submitted in IEEE signal proceesing letters, July 2017
% 
% ## ----------------------------------------------------------------------------
% ##
% ##   File: main.m
% ##   Copyright (c) <2017> Signal and System Theory Group, 
% ##                        Univ. of Paderborn, http://sst.upb.de
% ##                        https://github.com/SSTGroup/Correlation-Analysis-in-High-Dimensional-Data
% ##
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
% ##   3.) Publication, distribution, sublicensing, and/or selling of
% ##       copies or parts of the Software requires special agreements
% ##       with the Signal and System Theory Group, University of Paderborn,
% ##       and is in general not permitted.
% ##
% ##   4.) Modifications or contributions to the Software must be
% ##       published under this license. 
% ##   
% ##   5.) Any publication that was created with the help of this Software  
% ##       or parts of this Software must include a citation of the paper 
% ##       referenced above.
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
% ##   about bugs. 
% ##
% ##
% ##   Author: Tanuj Hasija, <tanuj.hasija@sst.upb>
% ##   Date: 18.08.2017
% ##   Version: v1 
% ## ----------------------------------------------------------------------------

clc;
clear variables;
close all;

m=60;  % total number of sensors in data set
d = 4; % total number of improper signals
f=0;  % total number of proper signals

RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));

var_d = 5*ones(1,d); % variance of improper signals
var_f = 2*ones(1,f); % variance of proper signals
K = [1,0.9,0.8,0.6]; % circularity coefficients of d improper signals

colorness = 'white'; % color / white noise
var_n = 1; % variance of noise
AR= [ 1/2 sqrt(7)/4 1/2 1/4 ]; % colored noise filter coefficient

num_iterations =  5*1e2; % number of iterations
Samples_Vector = [10:10:100,120,150,170,200:50:300];

for sample = 1:length(Samples_Vector)
    
    M = Samples_Vector(sample);
    display(['Samples=' num2str(M)]);
    
    parfor i = 1:num_iterations
        display(['iteration=' num2str(i)]);
        
        X = GenerateData_Array_Complex(m,d,f,var_d,var_f,K,var_n,colorness,AR,M);
        
        [~,~,Vx1] = svd(X,'econ');
        [~,~,Vx2] = svd(conj(X),'econ');
        
        rmax = min(floor(M/3),m);
        Pfa1 = 0.001; Pfa2 = 0.005;
        d_cap_max_min_ITC(i) = Max_Min_ITC_Impropriety(M,Vx1,Vx2,rmax); % detector (14)
        d_cap_max_min_pfa1(i) = Max_Min_HT_Impropriety(M,Pfa1,Vx1,Vx2,rmax); % detector (20)
        d_cap_max_min_pfa2(i) = Max_Min_HT_Impropriety(M,Pfa2,Vx1,Vx2,rmax); % detector (20)
        [~, d_cap] = ncPCA_Detector(X); % ncPCA detector [9]
        d_cap_ncpca(i) = d_cap;
        
        warning('off','all');
        
    end
    
    P_d_max_min_pfa1(sample) = length(find(d_cap_max_min_pfa1 == d))/num_iterations;
    P_d_max_min_pfa2(sample) = length(find(d_cap_max_min_pfa2 == d))/num_iterations;
    P_d_max_min_ITC(sample) = length(find(d_cap_max_min_ITC == d))/num_iterations;
    P_d_ncpca(sample) = length(find(d_cap_ncpca == d))/num_iterations;
    
end

figure();
hold on;
plot(Samples_Vector,P_d_max_min_ITC,'bo-','markersize',10);
plot(Samples_Vector,P_d_max_min_pfa1,'rs-','markersize',10);
plot(Samples_Vector,P_d_max_min_pfa2,'r+-','markersize',10);
plot(Samples_Vector,P_d_ncpca,'kd-','markersize',10);
a = xlabel('Number of Samples','fontsize',20,'FontName','Times New Roman');
b = ylabel('Probability of Detection','fontsize',20,'FontName','Times New Roman');
c = legend('Detector (14)','Detector (20) $P_\textrm{fa} = 0.001$','Detector (20) $P_\textrm{fa} = 0.005$','ncPCA detector [9]');
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(c,'interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',20);

