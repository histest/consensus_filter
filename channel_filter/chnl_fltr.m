function [x_CF, P_CF] = chnl_fltr(x_p, P_p, x_1, P_1, x_2, P_2)

P_CF = inv(inv(P_1) + inv(P_2) - inv(P_p) );
x_CF = P_CF * (inv(P_1) * x_1 + inv(P_2) * x_2 - inv(P_p) * x_p);