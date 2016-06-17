clear

load('post_samples_2dot5_1.mat')

theta_vec = 0:0.01:pi;
% non-stationary variance function
knots = [0 0 0 0 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec);

b_mat(:, 1) = 1;

eta_est_2dot5 = mean(post_samples.eta, 2);

load('post_samples_3_1.mat')

eta_est_3 = mean(post_samples.eta, 2);

load('post_samples_4_1.mat')

eta_est_4 = mean(post_samples.eta, 2);

plot(theta_vec, exp(b_mat*eta_est_2dot5), 'LineWidth', 2)
hold on
plot(theta_vec, exp(b_mat*eta_est_3), 'LineWidth', 2)
plot(theta_vec, exp(b_mat*eta_est_4), 'LineWidth', 2)
plot(theta_vec, exp(-(theta_vec-pi/2).^2/2/(pi/4)^2), 'LineWidth', 2)
lg = legend('2.5', '3', '4', 'true');
lg.FontSize = 12;
set(gca, 'FontSize', 12)
