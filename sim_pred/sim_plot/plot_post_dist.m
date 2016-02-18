clear
load('data_sim_2dot5')
load('Y_pred_needlet_2dot5_1_long.mat')
load('Y_pred_Matern_2dot5_1.mat')

phi_samples = phi_vec(index);
theta_samples = theta_vec(index);
Y_samples = Y(index);

scatter(phi_samples, theta_samples, [], Y_samples, '.')
axis equal
axis tight
colorbar
caxis([-max(abs(Y_samples)) max(abs(Y_samples))])
hold on
outlier = find(std_Y_pred_Matern<std(Y_pred_all(index_pred, :), 0, 2));
plot(phi_vec(index_pred(outlier)), theta_vec(index_pred(outlier)), 'ko')

figure
scatter(phi_samples, theta_samples, [], Y_samples, '.')
axis equal
axis tight
colorbar
caxis([-max(abs(Y_samples)) max(abs(Y_samples))])
hold on
outlier = find(abs(Y_pred_Matern(index_pred))<abs(Y_pred_needlet(index_pred)));
plot(phi_vec(index_pred(outlier)), theta_vec(index_pred(outlier)), 'ko')

Y_err_needlet = Y(index_pred)-Y_pred_needlet(index_pred);
Y_err_Matern = Y(index_pred)-Y_pred_Matern(index_pred);
figure
scatter(phi_samples, theta_samples, [], Y_samples, '.')
axis equal
axis tight
colorbar
caxis([-max(abs(Y_samples)) max(abs(Y_samples))])
hold on
outlier = find(abs(Y_err_Matern)<abs(Y_err_needlet));
plot(phi_vec(index_pred(outlier)), theta_vec(index_pred(outlier)), 'ko')

loc1 = 900;
[f, xi] = ksdensity(Y_pred_all(index_pred(loc1), :));
figure
plot(xi, f, 'LineWidth', 1.5)
hold on
ff = normpdf(xi, Y_pred_Matern(index_pred(loc1)), std_Y_pred_Matern(loc1));
plot(xi, ff, 'r', 'LineWidth', 1.5)
 
loc2 = 2000;
mean_Matern = Y_pred_Matern(index_pred(loc2));
std_Matern = std_Y_pred_Matern(loc2);
[f, xi] = ksdensity(Y_pred_all(index_pred(loc2), :), (mean_Matern-3*std_Matern):0.1:(mean_Matern+4*std_Matern));
figure
plot(xi, f, 'LineWidth', 1.5)
hold on
ff = normpdf(xi, mean_Matern, std_Matern);
plot(xi, ff, 'r', 'LineWidth', 1.5)

loc3 = 2010;
[f, xi] = ksdensity(Y_pred_all(index_pred(loc3), :));
figure
plot(xi, f, 'LineWidth', 1.5)
hold on
ff = normpdf(xi, Y_pred_Matern(index_pred(loc3)), std_Y_pred_Matern(loc3));
plot(xi, ff, 'r', 'LineWidth', 1.5)

loc4 = 2330;
[f, xi] = ksdensity(Y_pred_all(index_pred(loc4), :));
figure
plot(xi, f, 'LineWidth', 1.5)
hold on
ff = normpdf(xi, Y_pred_Matern(index_pred(loc4)), std_Y_pred_Matern(loc4));
plot(xi, ff, 'r', 'LineWidth', 1.5)

figure
scatter(phi_samples, theta_samples, [], Y_samples, '.')
axis equal
axis tight
colorbar
caxis([-max(abs(Y_samples)) max(abs(Y_samples))])
hold on
text(phi_vec(index_pred(loc1)), theta_vec(index_pred(loc1)), '1')
text(phi_vec(index_pred(loc2)), theta_vec(index_pred(loc2)), '2')
text(phi_vec(index_pred(loc3)), theta_vec(index_pred(loc3)), '3')
text(phi_vec(index_pred(loc4)), theta_vec(index_pred(loc4)), '4')