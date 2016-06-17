function cmp_pred(seed, name, extra)

% load data
load(['data_sim_', name, '.mat'])
load(['Y_pred_needlet_', name, '_', num2str(seed), extra, '.mat'])
load(['Y_pred_Matern_', name, '_', num2str(seed), extra, '.mat'])

err_needlet = Y-Y_pred_needlet;
err_Matern = Y-Y_pred_Matern;

theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);

subplot(2, 2, 1)
pcolor(phi_mat, theta_mat, reshape(Y, size(phi_mat)))
shading flat
axis image
hold on
plot([pi-pi/18 pi-pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
plot([pi+pi/18 pi+pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
scatter(phi_samples, theta_samples, 'k.')
caxis([-max(abs(Y)) max(abs(Y))])
colorbar
xlabel('Longitude (rad)')
ylabel('Co-latitude (rad)')
title('Data')
colormap(jet)

cmax = max(max(abs(err_needlet)), max(abs(err_Matern)));
subplot(2, 2, 2)
pcolor(phi_mat, theta_mat, reshape(err_needlet, size(phi_mat)))
shading flat
axis image
hold on
plot([pi-pi/18 pi-pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
plot([pi+pi/18 pi+pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
scatter(phi_samples, theta_samples, 'k.')
caxis([-cmax, cmax])
colorbar
xlabel('Longitude (rad)')
ylabel('Co-latitude (rad)')
title('Error (nonGau-need)')
colormap(jet)

subplot(2, 2, 3)
pcolor(phi_mat, theta_mat, reshape(err_Matern, size(phi_mat)))
shading flat
axis image
hold on
plot([pi-pi/18 pi-pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
plot([pi+pi/18 pi+pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
scatter(phi_samples, theta_samples, 'k.')
caxis([-cmax cmax])
colorbar
xlabel('Longitude (rad)')
ylabel('Co-latitude (rad)')
title('Error (Gau-Matern)')
colormap(jet)

subplot(2, 2, 4)
err_diff = abs(err_Matern)-abs(err_needlet);
pcolor(phi_mat, theta_mat, reshape(err_diff, size(phi_mat)))
shading flat
axis image
hold on
plot([pi-pi/18 pi-pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
plot([pi+pi/18 pi+pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
scatter(phi_samples, theta_samples, 'k.')
caxis([-max(abs(err_diff)), max(abs(err_diff))])
colorbar
xlabel('Longitude (rad)')
ylabel('Co-latitude (rad)')
title('Error diff')
colormap(jet)

end
