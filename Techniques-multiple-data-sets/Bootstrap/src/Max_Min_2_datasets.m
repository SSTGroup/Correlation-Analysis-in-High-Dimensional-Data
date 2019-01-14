% ## ----------------------------------------------------------------------------
% ##
% ##   File: Max_Min_2_datasets.m
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
% ##   Author: Tanuj Hasija
% ##
% ## ----------------------------------------------------------------------------

function [d_cap, r_d_1, r_d_2] = Max_Min_2_datasets(M,P_fa, Vx, Vy, rmax)
% MATLAB implementation based on "Determining the number of correlated signals
% between two data sets using PCA-CCA when sample support is extremely small",
% Y. Song, P. J. Schreier, and N. J. Roseveare, ICASSP 2015

d_cap_vector = zeros(1,rmax);

for r1=1: rmax
    
    for r2=1: rmax
        
        K_r = svd(Vx(:,1:r1)'*Vy(:,1:r2));
        
        r=min(r1,r2);
        for s=0:r-1
            
            F_K_r = 1 - K_r(s+1:r).^2;
            
            C_s_r = -(2*M- (r1+r2+1))* log(prod(F_K_r));
            
            dof = 2*(r1-s)*(r2-s);
            T_s_r = chi2inv(1-P_fa,dof);
            if (C_s_r < T_s_r)
                break;
            end
            if(s == r-1)
                s=r;
            end
        end
        
        d_cap_vector(r1,r2) = s;
    end
end

[L1, idx2] = max(d_cap_vector, [], 2);
[d_cap, idx1] = max(L1, [], 1);
idx2=idx2(idx1);
r_d_1 = idx1; r_d_2 = idx2;

end
