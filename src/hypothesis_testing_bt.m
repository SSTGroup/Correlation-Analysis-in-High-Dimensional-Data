function d_cap = hypothesis_testing_bt(alpha,K,K_star_matrix,B)


m1 = length(K);
d_cap = m1-1;
for d=0:m1-1
    
    T3 = calc_statistic(K,d);
    
    parfor b=1:B
        T3_star(b) = calc_statistic(K_star_matrix(:,b),d);
        T3_null(b) = T3_star(b) - T3;
        if(abs(T3) <= abs(T3_null(b)))
            Indicator(b) = 1;
        else
            Indicator(b) = 0;
        end
    end
    p_value = sum(Indicator)/B;
    if(p_value >= alpha)
        d_cap = d;
        break;
    end
end


end

function T = calc_statistic(K,d)
m1 = length(K);

T = ( (1/(m1-d))*sum(K(d+1:m1)) ) - (prod(K(d+1:m1).^(1/(m1-d))) );


end
