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
