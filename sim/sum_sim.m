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

n_line = size(eta_est_all, 2);
color_mat = colormap(jet(n_line));
subplot('position', [0.1 0.1 0.7 0.85])
fitted = exp(b_mat*eta_est_all);
hold on
for i = 1:n_line
    h = plot(theta, fitted(:, i));
    % 0.3 represents the transparency
    h.Color = [color_mat(i, :) 0.5];
end
plot(theta, exp(b_mat*eta), 'k', 'LineWidth', 2)
axis tight
xlabel('Co-latitude (rad)')
ylabel('Standard deviation')
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
