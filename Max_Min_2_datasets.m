function [d_cap, r_d_1, r_d_2] = Max_Min_2_datasets(M,P_fa, Vx, Vy, rmax)
% MATLAB implementation based on ?Determining the number of correlated signals
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
