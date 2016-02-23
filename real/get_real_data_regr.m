% get residuals by spherical harmonic regression

load('svd.mat')
load('theta_phi.mat')

r = double(r);

T = size(r, 1);

L = 3;
M = 3;

resid_all = zeros(size(r));
for t = 1:T
    [coef, resid, y_hat, ~] = SCHA_regr(reshape(r(t, :), size(phi)), theta, phi, L, M);
    resid_all(t, :) = resid;
end

index = 1;
plot_pot(reshape(resid_all(index, :), size(phi)), phi, theta, 1000, max(abs(resid_all(index, :))));

save('data_EOF_regr.mat', 'resid_all', 'theta', 'phi', '-v7.3')