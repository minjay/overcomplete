% summarize sim results
eta_est_all = [];
sigma_j_est_all = [];
tau_est_all = [];
for i = 1:10
    load(['sim_rep', num2str(i), '_2.mat'])
    eta_est_all = [eta_est_all eta_est];
    sigma_j_est_all = [sigma_j_est_all sigma_j_est(2, :)];
    tau_est_all = [tau_est_all tau_est];
end

theta = 0:0.01:pi;

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta);

b_mat(:, 1) = 1;

subplot('position', [0.1 0.575 0.6 0.375])
boxplot(eta_est_all', 'Labels', {'0', '1', '2', '3', '4'})
yl = ylim;
ylim([min(eta)-0.1 yl(2)])
for i = 1:5
    line([i-0.4 i+0.4], [eta(i) eta(i)], 'LineWidth', 2)
end
title('\eta')
set(gca, 'FontSize', 12)

subplot('position', [0.1 0.1 0.6 0.375])
plot(theta, exp(b_mat*eta_est_all))
hold on
plot(theta, exp(b_mat*eta), 'k', 'LineWidth', 2)
axis tight
xlabel('Co-latitude')
ylabel('Standard deviation')
set(gca, 'FontSize', 12)

subplot('position', [0.8 0.575 0.15 0.375])
boxplot(sigma_j_est_all', 'widths', 0.4)
line([1-0.4 1+0.4], [sigma_j(2) sigma_j(2)], 'LineWidth', 2)
title('\sigma_3')
set(gca, 'FontSize', 12)

subplot('position', [0.8 0.1 0.15 0.375])
boxplot(tau_est_all', 'widths', 0.4)
line([1-0.4 1+0.4], [tau tau], 'LineWidth', 2)
title('\tau')
set(gca, 'FontSize', 12)