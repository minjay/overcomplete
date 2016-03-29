% 1. EOF residuals
% 2. get residuals by spherical harmonic regression

load('svd.mat')
load('theta_phi.mat')

r = double(r);
r = r(:, 1:16920);
phi = phi(:, 1:end-1);
theta = theta(:, 1:end-1);

emp_std_vec = std(r, 0, 1);

T = size(r, 1);

L = 3;
M = 3;

resid_all = zeros(size(r));
for t = 1:100 
    [coef, resid, y_hat, ~] = SCHA_regr(reshape(r(t, :)./emp_std_vec, size(phi)), theta, phi, L, M);
    resid_all(t, :) = resid.*emp_std_vec';
end

index = 1;
plot_pot(reshape(resid_all(index, :), size(phi)), phi, theta, 1000, max(abs(resid_all(index, :))));

resid_all=resid_all(1:100,:);
save('data_EOF_regr_new.mat', 'resid_all', 'theta', 'phi', '-v7.3')