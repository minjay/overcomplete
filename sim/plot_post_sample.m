load('post_samples_rand_sim.mat')

figure
for i = 1:5
    subplot(3, 2, i)
    plot(post_samples.eta(i, :))
    title(['eta', num2str(i)])
    axis tight
end

subplot(3, 2, 6)
plot(post_samples.sigma_j_sq(2, :))
title('sigma2sq')
axis tight

theta = 0:0.01:pi;

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta);

b_mat(:, 1) = 1;

std_vec = exp(b_mat*eta);

figure
plot(theta, std_vec, 'b', 'LineWidth', 2)
hold on
eta_est = mean(post_samples.eta(:, 1001:end), 2);
std_vec_est = exp(b_mat*eta_est);
plot(theta, std_vec_est, 'r', 'LineWidth', 2)
legend('True', 'Fitted')
