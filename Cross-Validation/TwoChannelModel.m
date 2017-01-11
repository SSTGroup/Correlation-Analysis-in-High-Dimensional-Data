% ## ----------------------------------------------------------------------------
% ##
% ##   File: TwoChannelModel.m
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
%   Sample generator for the two-channel model.
%
%   Call:
%       [X,Y]=TwoChannelModel(n,m,d,f,p,sigmadx2,sigmady2,sigmafx2,sigmafy2,sigma2,M,mixing,color,MAcoeff)
%   Input:
%       n:       Dimension of x-channel observations
%       m:       Dimension of y-channel observations
%       d:       Number of correlated components
%       f:       Number of uncorrelated components ([fx,fy])
%       p:       Correlation coefficients
%       sigmad2: Vector of variances of correlated components
%       sigmaf2: Vector of variances of uncorrelated components
%       sigma2:  Noise variance
%       M:       Number of observations
%       mixing:  'orth' if mixing matrices are unitary, 'randn' if each
%                entry follows a standard normal distribution
%       color:   'white' for white noise, 'colored' for colored noise
%       MAcoeff: moving average coefficients for colored noise
%       ARcoeff: auto-regressive coefficients for colored noise
%   Output:
%       X:       X-channel observations
%       Y:       Y-channel observations
%
%

function [X,Y]=TwoChannelModel(n,m,d,f,p,sigmadx2,sigmady2,sigmafx2,sigmafy2,sigma2,M,mixing,color,MAcoeff,ARcoeff)

fx=f(1);
fy=f(2);

switch lower(mixing)
    case 'orth'
        Ax=orth(randn(n,d+fx));
        Ay=orth(randn(m,d+fy));
    case 'randn'
        Ax=randn(n,d+fx);
        Ay=randn(m,d+fy);
    otherwise
        error('Unknown mixing matrix property');
end

Rsxsx=blkdiag(diag(sigmadx2),diag(sigmafx2)); % covariance matrix of sx
Rsysy=blkdiag(diag(sigmady2),diag(sigmafy2)); % covariance matrix of sy
Rsxsy=blkdiag(sqrt(diag(sigmadx2)*diag(sigmady2))*diag(p),zeros(fx,fy)); % cross-covariance matrix of sx and sy
Rs=[Rsxsx  Rsxsy ; Rsxsy'  Rsysy]; % covariance matrix of the jointly Gaussian distribution

S=sqrtm(Rs)*randn(2*d+fx+fy,M);
Sx=S(1:d+fx,:); % realizations of sx
Sy=S(d+fx+1:end,:); % realizations of sy
Nx=sqrt(sigma2)*randn(n,M); % noise samples in the x channel
Ny=sqrt(sigma2)*randn(m,M); % noise samples in the y channel

switch lower(color)
    case 'white'
        Nx=Nx;
        Ny=Ny;
    case 'colored'
        Nx=filter(MAcoeff,ARcoeff,Nx);
        Ny=filter(MAcoeff,ARcoeff,Ny);
    otherwise
        error('Unknown noise color option');
end

X=Ax*Sx+Nx; % x channel observations
Y=Ay*Sy+Ny; % y channel observations
