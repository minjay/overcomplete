load('svd.mat')
load('theta_phi.mat')

r = double(r);
% remove the obs when theta=0
% remove the constant obs
r = r(:, 361:(end-360));
phi = phi(:, 2:(end-1));
theta = theta(:, 2:(end-1));

T = size(r, 1);

L = 3;
M = 3;
multi = 4;

times = 1:10:T;
emp_std_vec = std(r(times, :), 0, 1);

resid_all = zeros(length(times), size(r, 2));
for t = 1:length(times)
    [coef, resid, y_hat, ~] = SCHA_regr(reshape(r(times(t), :)./emp_std_vec, size(phi)), theta, phi, L, M, multi);
    resid_all(t, :) = resid'.*emp_std_vec;
end