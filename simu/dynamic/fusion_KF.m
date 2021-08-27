function [c_new, p_new] = fusion_KF(c_old, p_old, acc1, acc2, k)

fs = 1/60;
Ts = 60;
% C_0 = [8.7e-4, 6.5e-3, 3.1e-2, 0.12]'; %need to be modified
H_1 =[(k^5 - (k - 1)^5) / 5 * Ts^5, (k^4 - (k - 1)^4) / 4 * Ts^4,...
      (k^3 - (k - 1)^3) / 3 * Ts^3, (k^2 - (k - 1)^2) / 2 * Ts^2] ;
delta_x = acc1 - acc2; %low - high
Q1 = diag([1e1, 1e1, 1e1, 1e1]); %need to be modified
r1 = 1.01e-8; %need to be modified
p_minus1 = p_old + Q1;
gain1 = p_minus1 * H_1' * inv(r1 + H_1 * p_minus1 *H_1');
c_new = c_old + gain1 * (delta_x - H_1 * c_old); %out
p_new = (eye(4) - gain1 * H_1) * p_minus1; %out
drift = H_1 * c_new;

end