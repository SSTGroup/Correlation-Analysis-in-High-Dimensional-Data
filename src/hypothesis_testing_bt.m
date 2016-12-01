% ## ----------------------------------------------------------------------------
% ##
% ##   File: hypothesis_testing_bt.m
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


function d_cap = hypothesis_testing_bt(alpha,K,K_star_matrix,B)


m1 = length(K);
d_cap = m1-1;
for d=0:m1-1
    
    T3 = calc_statistic(K,d);
    
    for b=1:B
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
