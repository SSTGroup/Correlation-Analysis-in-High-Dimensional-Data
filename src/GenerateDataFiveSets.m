## ----------------------------------------------------------------------------
##
##   File: GenerateDataFiveSets.m
##   Copyright (c) <2016> <University of Paderborn>
##   Permission is hereby granted, free of charge, to any person
##   obtaining a copy of this software and associated documentation
##   files (the "Software"), to deal in the Software without restriction,
##   including without limitation the rights to use, copy, modify and
##   merge the Software, subject to the following conditions:
##
##   1.) The Software is used for non-commercial research and
##       education purposes.
##
##   2.) The above copyright notice and this permission notice shall be
##       included in all copies or substantial portions of the Software.
##
##   3.) Publication, Distribution, Sublicensing, and/or Selling of
##       copies or parts of the Software requires special agreements
##       with the University of Paderborn and is in general not permitted.
##
##   4.) Modifications or contributions to the software must be
##       published under this license. The University of Paderborn
##       is granted the non-exclusive right to publish modifications
##       or contributions in future versions of the Software free of charge.
##
##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
##   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
##   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
##   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
##   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
##   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
##   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
##   OTHER DEALINGS IN THE SOFTWARE.
##
##   Persons using the Software are encouraged to notify the
##   Department of Communications Engineering at the University of Paderborn
##   about bugs. Please reference the Software in your publications
##   if it was used for them.
##
##
##   Author: Tanuj Hasija
##
## ----------------------------------------------------------------------------

function [X1,X2,X3,X4,X5] =  GenerateDataFiveSets(m,Q,M,var_1,var_2,var_3,var_4,var_5,R12,R31,R14,R15,R23,R24,R25,R34,R35,R45,var_n,A)
% Generate 5 Data Sets with linear signal-plus-noise model
% X = A*S + N;

m1=m(1); m2=m(2); m3=m(3); m4=m(4); m5=m(5);
Q1=Q(1); Q2=Q(2); Q3=Q(3); Q4=Q(4); Q5=Q(5);

% Source Generation
% Generating Source Matrices according to the variance and correlation coefficient specified above

idx1= Q1; idx2=Q2+Q1; idx3=Q1+Q2+Q3; idx4=Q1+Q2+Q3+Q4; idx5=Q1+Q2+Q3+Q4+Q5;

E = zeros(idx5);

E(1:idx1,1:idx1) = diag(var_1/2);
E(idx1+1:idx2,idx1+1:idx2) = diag(var_2/2);
E(idx2+1:idx3,idx2+1:idx3) = diag(var_3/2);
E(idx3+1:idx4,idx3+1:idx4) = diag(var_4/2);
E(idx4+1:idx5,idx4+1:idx5) = diag(var_5/2);

E(1:idx1,idx1+1:idx2) = R12; E(idx1+1:idx2,1:idx1) = R12';
E(idx1+1:idx2,idx2+1:idx3) = R23; E(idx2+1:idx3,idx1+1:idx2) = R23';
E(idx2+1:idx3,1:idx1) = R31; E(1:idx1,idx2+1:idx3) = R31';
E(1:idx1,idx3+1:idx4) = R14; E(idx1+1:idx2,idx3+1:idx4) = R24; E(idx2+1:idx3,idx3+1:idx4) = R34;
E(idx3+1:idx4,1:idx1) = R14'; E(idx3+1:idx4,idx1+1:idx2) = R24'; E(idx3+1:idx4,idx2+1:idx3) = R34';
E(1:idx1,idx4+1:idx5) = R15; E(idx1+1:idx2,idx4+1:idx5) = R25; E(idx2+1:idx3,idx4+1:idx5) = R35; E(idx3+1:idx4,idx4+1:idx5) = R45;
E(idx4+1:idx5,1:idx1) = R15'; E(idx4+1:idx5,idx1+1:idx2) = R25'; E(idx4+1:idx5,idx2+1:idx3) = R35'; E(idx4+1:idx5,idx3+1:idx4) = R45';

mu = zeros(idx5,1);
F = mvnrnd(mu,E,M);

S1_r = F(:,1:idx1)'; S2_r = F(:,idx1+1:idx2)';  S3_r= F(:,idx2+1:idx3)'; S4_r= F(:,idx3+1:idx4)'; S5_r= F(:,idx4+1:idx5)';

F = mvnrnd(mu,E,M);
S1_c = F(:,1:idx1)'; S2_c = F(:,idx1+1:idx2)';  S3_c= F(:,idx2+1:idx3)'; S4_c= F(:,idx3+1:idx4)'; S5_c= F(:,idx4+1:idx5)';

S1 = (S1_r + 1i*S1_c); S2 = (S2_r + 1i*S2_c); S3 = (S3_r + 1i*S3_c); S4 = (S4_r + 1i*S4_c); S5 = (S5_r + 1i*S5_c);

% Mixing Matrices Generation

A1 = orth(rand(m1,Q1)+ 1i*rand(m1,Q1));
A2 = orth(rand(m2,Q2)+ 1i*rand(m2,Q2));
A3 = orth(rand(m3,Q3)+ 1i*rand(m3,Q3));
A4 = orth(rand(m4,Q4)+ 1i*rand(m4,Q4));
A5 = orth(rand(m5,Q5)+ 1i*rand(m5,Q5));

% Noise Matrices Generation

var_n_1 = var_n; var_n_2 = var_n; var_n_3 = var_n; var_n_4 = var_n; var_n_5 = var_n; % Variance of noise

% White Noise
N1 = sqrt(var_n_1/2)*(randn(m1,M)+ 1i*randn(m1,M));
N2 = sqrt(var_n_2/2)*(randn(m2,M)+ 1i*randn(m2,M));
N3 = sqrt(var_n_3/2)*(randn(m3,M)+ 1i*randn(m3,M));
N4 = sqrt(var_n_4/2)*(randn(m4,M)+ 1i*randn(m4,M));
N5 = sqrt(var_n_5/2)*(randn(m5,M)+ 1i*randn(m5,M));

%  Colored Noise

V_N1 = filter(1,A,N1);
V_N2 = filter(1,A,N2);
V_N3 =  filter(1,A,N3);
V_N4 =  filter(1,A,N4);
V_N5 =  filter(1,A,N5);

% Data sets generation

X1 = A1*S1 + V_N1;
X2 = A2*S2 + V_N2;
X3 = A3*S3 + V_N3;
X4 = A4*S4 + V_N4;
X5 = A5*S5 + V_N5;


end
