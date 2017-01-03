clear

load('ns.mat')
load('ns_deriv.mat')
load('theta_phi_R.mat')

grey = [0.4, 0.4, 0.4];

subplot('position', [0.075 0.1 0.4 0.8])
plot(theta(1, :)*4/pi*180, b_mat, 'LineWidth', 2)
y_range = get(gca, 'ylim');
hold on
plot([10 10]*4, y_range, '--', 'Color', grey, 'LineWidth', 1)
plot([20 20]*4, y_range, '--', 'Color', grey, 'LineWidth', 1)
plot([30 30]*4, y_range, '--', 'Color', grey, 'LineWidth', 1)
axis tight
xlabel('\theta^\prime (degree)')
ylabel('b(\theta^\prime)')
title('(a)')
set(gca, 'FontSize', 12)

subplot('position', [0.575 0.1 0.4 0.8])
plot(theta(1, :)*4/pi*180, b_mat_deriv, 'LineWidth', 2)
y_range = get(gca, 'ylim');
hold on
plot([10 10]*4, y_range, '--', 'Color', grey, 'LineWidth', 1)
plot([20 20]*4, y_range, '--', 'Color', grey, 'LineWidth', 1)
plot([30 30]*4, y_range, '--', 'Color', grey, 'LineWidth', 1)
axis tight
xlabel('\theta^\prime (degree)')
ylabel('db/d\theta^\prime')
title('(b)')
set(gca, 'FontSize', 12)
