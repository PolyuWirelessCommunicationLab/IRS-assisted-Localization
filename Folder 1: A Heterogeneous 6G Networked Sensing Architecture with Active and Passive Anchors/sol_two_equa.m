function [x_pssbl,y_pssbl] = sol_two_equa(x_bs,y_bs,r1,r2)
    a1 = x_bs(1);
    b1 = y_bs(1);
    a2 = x_bs(2);
    b2 = y_bs(2);
    [x_pssbl,y_pssbl] = circcirc(a1,b1,r1,a2,b2,r2);
end