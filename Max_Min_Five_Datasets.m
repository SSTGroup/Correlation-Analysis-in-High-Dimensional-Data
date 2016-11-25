function [d_cap, r_d] = Max_Min_Five_Datasets(M, P_fa, Vx1,Vx2x3x4x5,Vx2,Vx3x4x5...
    ,Vx3, Vx4x5, Vx4, Vx5, rmax)

% MATLAB Implementation based on the work "Detecting the dimension of the
% subspace correlated across multiple data sets in the sample poor regime",
% T. Hasija, Y. Song, P. J. Schreier, and D. Ram?rez,IEEE Workshop on
% Statistical Signal Processing, Palma de Mallorca, Spain, 2016

d_cap_vector = zeros(1,rmax);

parfor r=1: rmax
    
    Kx1_x2x3x4x5 = svd(Vx1(:,1:r)'*Vx2x3x4x5(:,1:4*r));
    Kx2_x3x4x5 = svd(Vx2(:,1:r)'*Vx3x4x5(:,1:3*r));
    Kx3_x4x5 = svd(Vx3(:,1:r)'*Vx4x5(:,1:2*r));
    Kx4_x5 = svd(Vx4(:,1:r)'*Vx5(:,1:r));
    
    for s=0:r-1
        
        F_Kx1_x2x3x4x5 = 1 - Kx1_x2x3x4x5(s+1:r).^2;
        F_Kx2_x3x4x5 = 1 - Kx2_x3x4x5(s+1:r).^2;
        F_Kx3_x4x5 = 1 - Kx3_x4x5(s+1:r).^2;
        F_Kx4_x5 = 1 - Kx4_x5(s+1:r).^2;
        
        C_s_r_5 = -(2*M-(5*r+1))*log(prod(F_Kx1_x2x3x4x5)) -(2*M-(4*r+1))*log(prod(F_Kx2_x3x4x5)) -(2*M-(3*r+1))*log(prod(F_Kx3_x4x5))...
            -(2*M-(2*r+1))*log(prod(F_Kx4_x5));
        
        dof = 2*(r-s)*(4*r-s) + 2*(r-s)*(3*r-s) + 2*(r-s)*(2*r-s)+ 2*(r-s)*(r-s);
        T_s_r = chi2inv(1-P_fa,dof);
        if (C_s_r_5 < T_s_r)
            break;
        end
    end
    d_cap_vector(r) = s;
end

[d_cap, r_d] = max(d_cap_vector);

end
