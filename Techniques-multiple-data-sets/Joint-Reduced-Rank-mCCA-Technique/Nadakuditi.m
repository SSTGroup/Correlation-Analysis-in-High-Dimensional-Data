function [Kmin obj] =  Nadakuditi(L,n,m,beta)

% beta = 2; % For complex case
k_vector = 0:min(n,m)-1;

for i=1:length(k_vector)
    
    k = k_vector(i);
    tk = ((sum(L(k+1:n).^2)/sum(L(k+1:n))^2)*(n-k) - (1+ n/m))*n - (2/beta-1)*(n/m);
    obj(i) = ((beta/4)*((m/n)^2)*tk^2) + 2*(k+1);
    
end

[obj_min,index] = min(obj);
Kmin = k_vector(index);
end