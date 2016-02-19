% get residuals by spherical harmonic regression

load('WHI_quad.mat')

all_Pot_N = all_Pot_N(1:720);

L = 18;
M = 3;

resid_all = zeros(720, numel(all_Pot_N{1}));
for t = 1:720
    [coef, resid, y_hat, ~] = SCHA_regr(all_Pot_N{t}, theta, phi, L, M);
    resid_all(t, :) = resid;
end

index = 1;
plot_pot(reshape(resid_all(index, :), size(phi)), phi, theta, 1000, max(abs(resid_all(index, :))));

resid = resid_all(index, :);
save('data_regr.mat', 'resid', 'theta', 'phi')