function [x_upd, P_upd] = KF_upd(z, x_pred, P_pred, H, R)

IM = H * x_pred;
IS = R + H * P_pred * H';
K = P_pred * H' / IS;
x_upd = x_pred + K * (z - IM);
P_upd = P_pred - K * IS * K';