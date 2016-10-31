clear

% summarize sim results
eta_est_all = [];
sigma_j_est_all = [];
tau_est_all = [];
for i = 1:10
    load(['sim_rep', num2str(i), '.mat'])
    eta_est_all = [eta_est_all eta_est];
    sigma_j_est_all = [sigma_j_est_all sigma_j_est(2, :)];
    tau_est_all = [tau_est_all tau_est];
end

theta = 0:0.01:pi;

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta);

b_mat(:, 1) = 1;

subplot('position', [0.1 0.1 0.7 0.85])
fitted = exp(b_mat*eta_est_all);
med = median(fitted, 2);
lb = quantile(fitted, 0.05, 2);
ub = quantile(fitted, 0.95, 2);
plot(theta, exp(b_mat*eta), 'k', 'LineWidth', 2)
hold on
plot(theta, med, 'b--', 'LineWidth', 2)
plot(theta, lb, 'r-.', 'LineWidth', 2)
plot(theta, ub, 'r-.', 'LineWidth', 2)
axis tight
xlabel('Co-latitude (rad)')
ylabel('Standard deviation')
legend('True', 'Pointwise median', '5% pointwise quantile',...
    '95% pointwise quantile', 'Location', 'Best')
set(gca, 'FontSize', 12)

subplot('position', [0.875 0.575 0.1 0.375])
boxplot(sigma_j_est_all', 'widths', 0.4)
line([1-0.4 1+0.4], [sigma_j(2) sigma_j(2)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
title('\sigma_j, j=3')
set(gca, 'FontSize', 12)

subplot('position', [0.875 0.1 0.1 0.375])
boxplot(tau_est_all', 'widths', 0.4)
line([1-0.4 1+0.4], [tau tau], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
title('\tau')
set(gca, 'FontSize', 12)
