% get residuals by spherical harmonic regression

load('WHI_quad.mat')

N = length(all_Pot_N);

L = 18;
M = 3;

resid_all = zeros(N, numel(all_Pot_N{1}));
for t = 1:N
    [coef, resid, y_hat, ~] = SCHA_regr(all_Pot_N{t}, theta, phi, L, M);
    resid_all(t, :) = resid;
end

save('resid_regr_all.mat', 'resid_all', 'theta', 'phi')