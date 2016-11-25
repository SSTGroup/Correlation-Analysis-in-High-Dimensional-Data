function dcap = Det_13_Song(X1,X2,X3,X4,X5,M,P_fa)
% MATLAB implementation based on ?Determining the number of signals correlated 
% across multiple data sets for small sample support", 
% Y. Song, T. Hasija, P. J. Schreier, and D. Ram?rez, EUSIPCO 2016 


% Max-Min for two data sets to estimate r1,r2,r3,r4,r5

[Ux1,Lx1,Vx1] = svd(X1,'econ');
[Ux2,Lx2,Vx2] = svd(X2,'econ');
[Ux3,Lx3,Vx3] = svd(X3,'econ');
[Ux4,Lx4,Vx4] = svd(X4,'econ');
[Ux5,Lx5,Vx5] = svd(X5,'econ');

% X1 and other data sets

m1 = size(X1,1); m2 = size(X2,1); m3 = size(X3,1); m4 = size(X4,1); m5 = size(X5,1);
rmax = min([floor(M/3),m1]);  % Assuming m1<=m2<=..<=m5
[d12, r_12_1, r_12_2] = Max_Min_2_datasets(M,P_fa, Vx1, Vx2, rmax);
[d13, r_13_1, r_13_3] = Max_Min_2_datasets(M,P_fa, Vx1, Vx3, rmax);
[d14, r_14_1, r_14_4] = Max_Min_2_datasets(M,P_fa, Vx1, Vx4, rmax);
[d15, r_15_1, r_15_5] = Max_Min_2_datasets(M,P_fa, Vx1, Vx5, rmax);

% X2 and other data sets
rmax = min([floor(M/3),m2]);  % Assuming m1<=m2<=..<=m5
[d23, r_23_2, r_23_3] = Max_Min_2_datasets(M,P_fa, Vx2, Vx3, rmax);
[d24, r_24_2, r_24_4] = Max_Min_2_datasets(M,P_fa, Vx2, Vx4, rmax);
[d25, r_25_2, r_25_5] = Max_Min_2_datasets(M,P_fa, Vx2, Vx5, rmax);

% X3 and other data sets
rmax = min([floor(M/3),m3]);  % Assuming m1<=m2<=..<=m5
[d34, r_34_3, r_34_4] = Max_Min_2_datasets(M,P_fa, Vx3, Vx4, rmax);
[d35, r_35_3, r_35_5] = Max_Min_2_datasets(M,P_fa, Vx3, Vx5, rmax);

% X4 and other data sets
rmax = min([floor(M/3),m4]); % Assuming m1<=m2<=..<=m5
[d45, r_45_4, r_45_5] = Max_Min_2_datasets(M,P_fa, Vx4, Vx5, rmax);

r1 = max([r_12_1, r_13_1, r_14_1, r_15_1]);
r2 = max([r_12_2, r_23_2, r_24_2, r_25_2]);
r3 = max([r_13_3, r_23_3, r_34_3, r_35_3]);
r4 = max([r_14_4, r_24_4, r_34_4, r_45_4]);
r5 = max([r_15_5, r_25_5, r_35_5, r_45_5]);

X1r = Ux1(:,1:r1)*Lx1(1:r1,1:r1)*Vx1(:,1:r1)';
X2r = Ux2(:,1:r2)*Lx2(1:r2,1:r2)*Vx2(:,1:r2)';
X3r = Ux3(:,1:r3)*Lx3(1:r3,1:r3)*Vx3(:,1:r3)';
X4r = Ux4(:,1:r4)*Lx4(1:r4,1:r4)*Vx4(:,1:r4)';
X5r = Ux5(:,1:r5)*Lx5(1:r5,1:r5)*Vx5(:,1:r5)';

R12r = 1/M*X1r*(X2r');
R13r = 1/M*X1r*(X3r');
R14r = 1/M*X1r*(X4r');
R15r = 1/M*X1r*(X5r');
R23r = 1/M*X2r*(X3r');
R24r = 1/M*X2r*(X4r');
R25r = 1/M*X2r*(X5r');
R34r = 1/M*X3r*(X4r');
R35r = 1/M*X3r*(X5r');
R45r = 1/M*X4r*(X5r');

R11r_sqrtinv = sqrt(M)*Ux1(:,1:r1)/(Lx1(1:r1,1:r1))*Ux1(:,1:r1)';
R22r_sqrtinv = sqrt(M)*Ux2(:,1:r2)/(Lx2(1:r2,1:r2))*Ux2(:,1:r2)';
R33r_sqrtinv = sqrt(M)*Ux3(:,1:r3)/(Lx3(1:r3,1:r3))*Ux3(:,1:r3)';
R44r_sqrtinv = sqrt(M)*Ux4(:,1:r4)/(Lx4(1:r4,1:r4))*Ux4(:,1:r4)';
R55r_sqrtinv = sqrt(M)*Ux5(:,1:r5)/(Lx5(1:r5,1:r5))*Ux5(:,1:r5)';

C12 = R11r_sqrtinv * R12r * R22r_sqrtinv;
C13 = R11r_sqrtinv * R13r * R33r_sqrtinv;
C14 = R11r_sqrtinv * R14r * R44r_sqrtinv;
C15 = R11r_sqrtinv * R15r * R55r_sqrtinv;
C23 = R22r_sqrtinv * R23r * R33r_sqrtinv;
C24 = R22r_sqrtinv * R24r * R44r_sqrtinv;
C25 = R22r_sqrtinv * R25r * R55r_sqrtinv;
C34 = R33r_sqrtinv * R34r * R44r_sqrtinv;
C35 = R33r_sqrtinv * R35r * R55r_sqrtinv;
C45 = R44r_sqrtinv * R45r * R55r_sqrtinv;

[F12,K12,G12] = svd(C12);
[F13,K13,G13] = svd(C13);
[F14,K14,G14] = svd(C14);
[F15,K15,G15] = svd(C15);
[F23,K23,G23] = svd(C23);
[F24,K24,G24] = svd(C24);
[F25,K25,G25] = svd(C25);
[F34,K34,G34] = svd(C34);
[F35,K35,G35] = svd(C35);
[F45,K45,G45] = svd(C45);

X1new  = F12(:,1:d12)'*R11r_sqrtinv*X1r;
X2new  = F13(:,1:d13)'*R11r_sqrtinv*X1r;
X3new  = F14(:,1:d14)'*R11r_sqrtinv*X1r;
X4new  = F15(:,1:d15)'*R11r_sqrtinv*X1r;
X5new  = F23(:,1:d23)'*R22r_sqrtinv*X2r;
X6new  = F24(:,1:d24)'*R22r_sqrtinv*X2r;
X7new  = F25(:,1:d25)'*R22r_sqrtinv*X2r;
X8new  = F34(:,1:d34)'*R33r_sqrtinv*X3r;
X9new  = F35(:,1:d35)'*R33r_sqrtinv*X3r;
X10new = F45(:,1:d45)'*R44r_sqrtinv*X4r;

Y1new  = G12(:,1:d12)'*R22r_sqrtinv*X2r;
Y2new  = G13(:,1:d13)'*R33r_sqrtinv*X3r;
Y3new  = G14(:,1:d14)'*R44r_sqrtinv*X4r;
Y4new  = G15(:,1:d15)'*R55r_sqrtinv*X5r;
Y5new  = G23(:,1:d23)'*R33r_sqrtinv*X3r;
Y6new  = G24(:,1:d24)'*R44r_sqrtinv*X4r;
Y7new  = G25(:,1:d25)'*R55r_sqrtinv*X5r;
Y8new  = G34(:,1:d34)'*R44r_sqrtinv*X4r;
Y9new  = G35(:,1:d35)'*R55r_sqrtinv*X5r;
Y10new = G45(:,1:d45)'*R55r_sqrtinv*X5r;

Xnew = [X1new;X2new;X3new;X4new;X5new;X6new;X7new;X8new;X9new;X10new];
Ynew = [Y1new;Y2new;Y3new;Y4new;Y5new;Y6new;Y7new;Y8new;Y9new;Y10new];
dim = [d12,d13,d14,d15,d23,d24,d25,d34,d35,d45];
r = 0;
for i=1:10 % L(L-1)/2 
    for j=i+1:10
        r = r + 1;
        Z1 = Xnew(sum(dim(1:i))-dim(i)+1:sum(dim(1:i)),:);
        Z2 = Ynew(sum(dim(1:j))-dim(j)+1:sum(dim(1:j)),:);
        Rz11 = 1/M*Z1*(Z1');
        Rz12 = 1/M*Z1*(Z2');
        Rz22 = 1/M*Z2*(Z2');
        C = sqrtm(inv(Rz11)) * Rz12 * sqrtm(inv(Rz22));
        ga(r,1:min(dim(i),dim(j))) = sort(svd(C),'descend');
    end
end

for d=0:min([d12,d13,d14,d15,d23,d24,d25,d34,d35,d45])
    f = @(n,ga,rx,ry,r)n*log(prod(1-ga(1:r).^2))+1/2*log(n)*2*r*(rx+ry-r);
    ICsum = 0;
    r = 0;
    for i=1:10
        for j=i+1:10
            r = r+1;
            ICsum = ICsum + f(M,ga(r,:),dim(i),dim(j),d);
        end
    end
    IC(:,d+1) = ICsum;
end
dcap = find(IC==min(IC))-1;

end
