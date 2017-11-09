% This is a MATLAB implementation supporting the paper
%
% "Determining the Dimension of the Improper Signal
% Subspace in Complex-Valued Data" by Tanuj Hasija, Christian Lameiro and
% Peter J. Schreier, IEEE Signal Processing Letters, vol. 24, no. 11, pp. 1606-1610, Nov. 2017.
% 
% ## ----------------------------------------------------------------------------
% ##
% ##   File: GenerateData_Array_Complex.m
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

function X =  GenerateData_Array_Complex(m,d,f,var_d,var_f,K,var_n,colorness,AR,M)
% GenerateData_Array_Complex.m - generate complex-valued ULA data.
%
% X = A*S + N
%
%
% A is the steering matrix
% S is the source matrix
% N is the noise matrix
%
% Input:
%
% m:         dimension of data / number of sensors in X
% d:         number of improper signals in S
% f:         number of proper signals in S
% var_d:     variance of improper signals
% var_f:     variance of proper signals
% K:         circularity coefficients of d improper signals
% var_n:     variance of noise
% colorness  'color' / 'white' noise;
% AR         Colored autoregressive noise filter coeffcient
% M:         number of samples
%
% Output:
%
% X:         data matrix m-by-M
%

Q = d+f; % total number of signals
fs=2000*1e4; %sampling frequency for electromagnetic waves
Ts=1/fs;
delta = 5; % degree separation between sources
degreesVector= [10*ones(1,Q)] + [0:delta:delta*(Q-1)];
c = 3e8; % speed of propogation
fc=[ones(1,Q)]*fs; % frequency of the signal
Lambda= c./fc; % wavelength of the signal
d =min(Lambda)/2;  % distance between adjacent sensors
w = 2*pi*fc;
kk=0 : m-1;
k=(1 :M)*Ts;

var_s = [var_d,var_f];  % variance of sources in the data set

if (d ~= 0)
    Rsrsi = diag([K,zeros(1,f)].*var_s/2); % source complementary covariance matrix
else
    Rsrsi = diag(zeros(1,n));
end

% Generating source matrix according to the variance and circularity coefficients specified above

E = zeros(2*Q);
E(1:2*Q,1:2*Q) = diag([var_s/2,var_s/2]);
E(1:Q,Q+1:2*Q) = Rsrsi; E(Q+1:2*Q,1:Q) = Rsrsi';

mu = zeros(2*Q,1);
Data = mvnrnd(mu,E,M);

S_r = Data(:,1:Q)'; S_i = Data(:,Q+1:2*Q)';

S = S_r + 1i*S_i;
A = zeros(m,Q); S_exp = zeros(Q,M);

for sig=1 : Q
    A(:,sig)=exp(1i.*kk*2*pi/Lambda(sig)*d*cos(degreesVector(sig)*pi/180)).';  % correponds to the steering vector , a(theta)
    S_exp(sig,:) = exp(1i*(w(sig)*k+2*pi*rand));
end

S = S.* S_exp;

switch lower(colorness)
    case 'white'
        N = sqrt(var_n/2)*(randn(m,M) + 1i*randn(m,M));
    case 'color'
        N = sqrt(var_n/2)*(randn(m,M) + 1i*randn(m,M));
        N = filter(1,AR,N);
end

X = A*S + N;

end