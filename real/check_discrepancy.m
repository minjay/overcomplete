load('svd.mat')
load('theta_phi.mat')

r = double(r);

L = 2;
M = 2;
t = 1;
[coef, resid, y_hat, ~] = SCHA_regr(reshape(r(t, :), size(theta)), theta, phi, L, M);
max_caxis = max(abs(r(t, :)));
subplot(1, 3, 1)
plot_pot_lite(reshape(r(t, :), size(phi)), phi, theta, 1000, max_caxis)
title('Resid')
subplot(1, 3, 2)
plot_pot_lite(reshape(y_hat, size(phi)), phi, theta, 1000, max_caxis)
title('Fitted by SPH')
subplot(1, 3, 3)
plot_pot_lite(reshape(resid, size(phi)), phi, theta, 1000, max_caxis)
title('Resid-SPH')

colorbar('Position', [0.925 0.05 0.025 0.9])
colormap(b2r(-max_caxis, max_caxis))