% 1. EOF residuals
% 2. get residuals by spherical harmonic regression

load('svd.mat')
load('theta_phi.mat')

r = double(r);
% remove the obs when theta=0
r = r(:, 361:end);
phi = phi(:, 2:end);
theta = theta(:, 2:end);

T = size(r, 1);

L = 3;
M = 3;

emp_std_vec = std(r, 0, 1);

times = 1:10:T;
resid_all = zeros(length(times), size(r, 2));
for t = 1:length(times)
    [coef, resid, y_hat, ~] = SCHA_regr(reshape(r(times(t), :)./emp_std_vec, size(phi)), theta, phi, L, M);
    resid_all(t, :) = resid'.*emp_std_vec;
end

index = 1;
plot_pot(reshape(resid_all(index, :), size(phi)), phi, theta, 1000, max(abs(resid_all(index, :))));

save('data_EOF_regr_new.mat', 'resid_all', 'theta', 'phi', '-v7.3')