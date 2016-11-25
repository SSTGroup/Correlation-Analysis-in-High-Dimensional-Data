% MATLAB Implementation of the work - 
% "Bootstrap-based Detection of the Number of Signals Correlated across
% Multiple Data Sets" by Tanuj Hasija, Yang Song, Peter J. Schreier and
% David Ramirez, Proceedings of the 50th Asilomar Conference on Signals,
% Systems and Computers 

clc; clear variables; close all;

Q1=6; Q2= 6; Q3=6; Q4=6; Q5=6; % Total no of signals in data sets
m1=14;m2=18;m3=20; m4=22; m5=24;% Total no of sensors in data sets
Q = [Q1,Q2,Q3,Q4,Q5];
m=[m1,m2,m3,m4,m5];

var_1 = [1 1 1 1 1 1];  % Variance of sources in data sets
var_2 = [1 1 1 1 1 1];
var_3 = [1 1 1 1 1 1];
var_4 = [1 1 1 1 1 1];
var_5 = [1 1 1 1 1 1];

M = 400; % Number of samples

% Cross Correlation matrices satisfying the assumption in [6],[12] and [13]
% i.e. d_ij = d for all (i,j) in {1,2,..L} and i != j
R12 = diag([0.8,  0.8, 0,   0    0   0 ].*sqrt(var_1/2).*sqrt(var_2/2));
R23 = diag([0.9,  0.9, 0  , 0    0   0 ].*sqrt(var_2/2).*sqrt(var_3/2));
R31 = diag([0.8,  0.8, 0,   0    0   0 ].*sqrt(var_3/2).*sqrt(var_1/2));
R14 = diag([0.85, 0.8, 0,   0,   0   0 ].*sqrt(var_1/2).*sqrt(var_4/2));
R24 = diag([0.9,  0.8, 0,   0    0   0 ].*sqrt(var_2/2).*sqrt(var_4/2));
R34 = diag([0.9,  0.8, 0,   0    0   0 ].*sqrt(var_3/2).*sqrt(var_4/2));
R15 = diag([0.85, 0.8, 0,   0    0   0 ].*sqrt(var_1/2).*sqrt(var_5/2));
R25 = diag([0.9,  0.75, 0,   0    0   0 ].*sqrt(var_2/2).*sqrt(var_5/2));
R35 = diag([0.9,  0.8, 0,   0    0   0 ].*sqrt(var_3/2).*sqrt(var_5/2));
R45 = diag([0.85, 0.8, 0,   0    0   0 ].*sqrt(var_4/2).*sqrt(var_5/2));

% Cross Correlation matrices not satisfying the assumption in [6],[12] and [13]
% R12 = diag([0.8,  0.8, 0.7,   0    0.7   0 ].*sqrt(var_1/2).*sqrt(var_2/2));
% R23 = diag([0.9,  0.9, 0.5  , 0    0   0.6 ].*sqrt(var_2/2).*sqrt(var_3/2));
% R31 = diag([0.8,  0.8, 0.7,   0    0   0.5 ].*sqrt(var_3/2).*sqrt(var_1/2));
% R14 = diag([0.85, 0.8, 0.7,   0.7,   0   0 ].*sqrt(var_1/2).*sqrt(var_4/2));
% R24 = diag([0.9,  0.8, 0.6,   0.6    0   0.6 ].*sqrt(var_2/2).*sqrt(var_4/2));
% R34 = diag([0.9,  0.8, 0,   0    0.5   0.6 ].*sqrt(var_3/2).*sqrt(var_4/2));
% R15 = diag([0.85, 0.8, 0,   0.7    0.5   0 ].*sqrt(var_1/2).*sqrt(var_5/2));
% R25 = diag([0.9,  0.75, 0,   0    0.5   0.6 ].*sqrt(var_2/2).*sqrt(var_5/2));
% R35 = diag([0.9,  0.8, 0,   0.5    0   0.5 ].*sqrt(var_3/2).*sqrt(var_5/2));
% R45 = diag([0.85, 0.8, 0,   0.6    0.6   0.6 ].*sqrt(var_4/2).*sqrt(var_5/2));

d = 2; % Number of signals correlated across all data sets

B = 500; % Number of Bootstrap iterations
alpha = 0.05; % Probability of false alarm

num_iterations = 0.05*1e2; % Number of Monte Carlo trials
SNR_Vector = [ -10:2.5:-5, -4 -3 -2.5 -1, 0:5:15 ];
A= [sqrt(8) sqrt(1) sqrt(1) ]; % AR Noise Filter coefficient vector

for snr = 1:length(SNR_Vector)
    
    var_n = var_1(1)/(10^(SNR_Vector(snr)/10));
    
    for iter=1:num_iterations
        
        display(['iteration=' num2str(iter)]);
        
        [X1,X2,X3,X4,X5] =  GenerateDataFiveSets(m,Q,M,var_1,var_2,var_3,var_4,var_5,...
                                    R12,R31,R14,R15,R23,R24,R25,R34,R35,R45,var_n,A);
        
        % Estimating the Singular Vectors of the data sets
        
        [~,~,Vx1] = svd(X1,'econ');
        [~,~,Vx2] = svd(X2,'econ');
        [~,~,Vx3] = svd(X3,'econ');
        [~,~,Vx4] = svd(X4,'econ');
        [~,~,Vx5] = svd(X5,'econ');
        
        % Singular values of the product of coherence matrices
        % svd((Vx1'*Vx2)) is same as svd(Cx1x2), where Cx1x2 is the
        % coherence matrix of X1 and X2. Refer to "Empirical canonical
        % correlation analysis in subspaces", Pezeshki et al, 2004 for more
        % details
        K = svd((Vx1'*Vx2)*(Vx2'*Vx3)*(Vx3'*Vx4)*(Vx4'*Vx5)*(Vx5'*Vx3)*(Vx3'*Vx1)*(Vx1'*Vx4)*(Vx4'*Vx2)*(Vx2'*Vx5)*(Vx5'*Vx1));
        
        % Bootstrap Operation
        parfor b=1:B
            
            [X1_star,I] = datasample(X1,M,2); % Bootstrap data matrix
            X2_star = X2(:,I);
            X3_star = X3(:,I);
            X4_star = X4(:,I);
            X5_star = X5(:,I);
            
            
            % Centering around the mean
            X1_star = X1_star - repmat(mean(X1_star,2),1,M);
            X2_star = X2_star - repmat(mean(X2_star,2),1,M);
            X3_star = X3_star - repmat(mean(X3_star,2),1,M);
            X4_star = X4_star - repmat(mean(X4_star,2),1,M);
            X5_star = X5_star - repmat(mean(X5_star,2),1,M);
            
            [~,~,Vy1] = svd(X1_star,'econ');
            [~,~,Vy2] = svd(X2_star,'econ');
            [~,~,Vy3] = svd(X3_star,'econ');
            [~,~,Vy4] = svd(X4_star,'econ');
            [~,~,Vy5] = svd(X5_star,'econ');
            K_star = svd((Vy1'*Vy2)*(Vy2'*Vy3)*(Vy3'*Vy4)*(Vy4'*Vy5)*(Vy5'*Vy3)*(Vy3'*Vy1)*(Vy1'*Vy4)*(Vy4'*Vy2)*(Vy2'*Vy5)*(Vy5'*Vy1));
            
            K_star_matrix(:,b) = K_star;
            
        end
        
        % Hypothesis Testing using Bootstrap
 
        d_cap_bt(iter) = hypothesis_testing_bt(alpha,K,K_star_matrix,B);
        
         % ITC MDL Detector in [6]
        
        [~,~, Vx2x3x4x5] = svd([X2;X3;X4;X5],'econ');
        [~,~, Vx3x4x5] = svd([X3;X4;X5],'econ');
        [~,~, Vx4x5] = svd([X4;X5],'econ');
        
        Kx1_x2x3x4x5 = svd(Vx1'*Vx2x3x4x5);
        Kx2_x3x4x5 = svd(Vx2'*Vx3x4x5);
        Kx3_x4x5 = svd(Vx3'*Vx4x5);
        Kx4_x5 = svd(Vx4'*Vx5);
        
        L=5; m = [m1,m2,m3,m4,m5];   
        K_ITC_MDL = zeros(L-1,max(m));
        K_ITC_MDL(1,1:m1) = Kx1_x2x3x4x5'; K_ITC_MDL(2,1:m2) = Kx2_x3x4x5'; K_ITC_MDL(3,1:m3) = Kx3_x4x5'; K_ITC_MDL(4,1:m4) = Kx4_x5';
        Det_ITC_MDL(iter) = ITC_MDL_Wu(M,m,L,K_ITC_MDL,log(M)/2,2);
        
        % Max-Min Multiple data sets [12]
        
        P_fa = 0.001;
        rmax = min([floor(2*M/(3*5)),m1]);
        [Det_Max_Min_5(iter), ~] = Max_Min_Five_Datasets(M,P_fa, Vx1,Vx2x3x4x5,Vx2, Vx3x4x5, Vx3, Vx4x5, Vx4, Vx5, rmax);
        
        % Detector [13]
        Det_13(iter) = Det_13_Song(X1,X2,X3,X4,X5,M,P_fa);
        
        warning('off','all')
        
    end
    
    Prob_detection_bt(snr) = length(find(d_cap_bt == d))/num_iterations;
    Prob_detection_ITC_MDL(snr) = length(find(Det_ITC_MDL == d))/num_iterations;  
    Prob_detection_Max_Min_5(snr) = length(find(Det_Max_Min_5 == d))/num_iterations;
    Prob_detection_Det_13(snr) = length(find(Det_13 == d))/num_iterations;
     
    Mean_bt(snr) = mean(d_cap_bt)
    Var_bt(snr) = var(d_cap_bt);
    
    Mean_ITC_MDL(snr) = mean(Det_ITC_MDL)
    Var_ITC_MDL(snr) = var(Det_ITC_MDL);
    Mean_Max_Min_5(snr) = mean(Det_Max_Min_5)
    Var_Max_Min_5(snr) = var(Det_Max_Min_5);
    Mean_Det_13(snr) = mean(Det_13)
    Var_Det_13(snr) = var(Det_13);
    
end


figure();
h3 = plot(SNR_Vector,Mean_bt,'ro--','markersize',12,'Linewidth',2); hold on;
h4 = plot(SNR_Vector,Mean_ITC_MDL,'ks-.','markersize',12,'Linewidth',2);
h5 = plot(SNR_Vector,Mean_Max_Min_5,'b*--','markersize',12,'Linewidth',2);
h6 = plot(SNR_Vector,Mean_Det_13,'mx-.','markersize',12,'Linewidth',2);
h7 = plot(SNR_Vector,d*ones(size(SNR_Vector)),'g-','Linewidth',2);
h = [h3,h4,h5,h6,h7];
a = xlabel('SNR (dB)','fontsize',20,'FontName','Times New Roman');
b = ylabel('Mean value of $\hat{d}$','fontsize',20,'FontName','Times New Roman');
c = legend('Proposed detector','ITC MDL [6]','Max-Min detector [12]','Detector based on [13]','True number','FontName','Times New Roman','fontsize',20);
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(c,'interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',20);

