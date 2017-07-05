clear

load('sim_rep2.mat')

theta = 0:0.01:pi;

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta);

b_mat(:, 1) = 1;

subplot('position', [0.1 0.575 0.5 0.375])
bh = boxplot(eta_est_all(2:end, :)', 'widths', 0.4, 'Labels', {'1', '2', '3', '4'});
set(bh(6, :), 'LineWidth', 1.5)
yl = ylim;
ylim([min(eta)-0.1 yl(2)])
for i = 1:4
    line([i-0.4 i+0.4], [eta(i+1) eta(i+1)], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k')
end
title('\eta')
set(gca, 'FontSize', 12)

subplot('position', [0.1 0.1 0.6 0.375])
fitted = exp(b_mat*eta_est_all);
med = median(fitted, 2);
lb = quantile(fitted, 0.05, 2);
ub = quantile(fitted, 0.95, 2);
plot(theta, med, 'r', 'LineWidth', 2)
hold on
plot(theta, exp(b_mat*eta), 'k--', 'LineWidth', 2)
plot(theta, lb, 'b-.', 'LineWidth', 2)
plot(theta, ub, 'b-.', 'LineWidth', 2)
axis tight
xlabel('Co-latitude (rad)')
ylabel('Standard deviation')
legend('Pointwise median', 'True', '5% pointwise quantile',...
    '95% pointwise quantile', 'Location', 'Best')
set(gca, 'FontSize', 12)

subplot('position', [0.7 0.575 0.25 0.375])
bh = boxplot(sigma_j_est_all', 'widths', 0.4, 'Labels', {'2', '3'});
set(bh(6, :), 'linewidth', 1.5)
for i = 1:2
    line([i-0.4 i+0.4], [sigma_j(i) sigma_j(i)], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k')
end
title('\sigma')
set(gca, 'FontSize', 12)

subplot('position', [0.8 0.1 0.15 0.375])
bh = boxplot(tau_est_all', 'widths', 0.4, 'Labels', {''});
set(bh(6, :), 'linewidth', 1.5)
line([1-0.4 1+0.4], [tau tau], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k')
title('\tau')
set(gca, 'FontSize', 12)

suptitle('(a)')