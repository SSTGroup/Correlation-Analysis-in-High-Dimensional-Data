## ----------------------------------------------------------------------------
##
##   File: ITC_MDL_Wu.m
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

function [Kmin, obj] = ITC_MDL_Wu(M,m,L,K,Cn,beta)
% MATLAB implementation based on ?Determination of number of sources with
% multiple arrays in correlated noise fields?, Y.Wu,K.W.Tam,andF.Li,
% Transactions on Signal Processing, 2002

% M - No. of Samples
% m - Vector of Dimension of Data sets
% L - No. of Data sets
% K - Sample Canonical Correlation Matrix (p-1 rows)
% Cn - Parameter for penalty term
% beta = 1 for Real Case, 2 for Complex

eta_max = min(m)-1;
eta_vector = 0:eta_max;
E = K.^2;
obj = zeros(size(eta_vector));

for i=1:length(eta_vector)
    eta = eta_vector(i);
    L_X = 0;
    V_n = 0;
    
    for r=1:L-1
        p = min(m(r:end)); % represents p_j in paper, which changes with iteration
        L_X = L_X + log(prod(1-E(r,eta+1:p)));
        
        if(beta==2)
            V_n =  V_n + beta*eta*(m(r)+ sum(m(r+1:end))-eta);
            
        else if(beta==1)
                m1 = m(r); m2 = sum(m(r+1:end));
                V_n =  V_n + (m1*m2 -(m1-eta)*(m2-eta) + m1*(m1+1)/2 + m2*(m2+1)/2);
                
            end
        end
        
    end
    
    if(beta==2)
        obj(i) = -M*L_X + V_n*Cn;
        
    else if(beta==1)
            obj(i) = -M/2*L_X + V_n*Cn; 
        end
    end
    
end

[~,index] = min(obj);
Kmin = eta_vector(index);

end
