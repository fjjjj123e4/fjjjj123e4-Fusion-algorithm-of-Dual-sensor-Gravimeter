function [acc_new, p_new] = noise_KF(acc_old, p_old, z)
Q1 = 1e-6; %need to be modified
r1 = 1.01e-2; %need to be modified
p_minus1 = p_old + Q1;
gain1 = p_minus1 * inv(r1 + p_minus1);
acc_new = acc_old + gain1 * (z - acc_old); %out
p_new = (eye(1) - gain1) * p_minus1; %out
end