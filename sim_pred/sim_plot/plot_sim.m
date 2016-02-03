function plot_sim(seed, model, name)
% e.g.: seed = 1; model = 'Matern'; name = '3';

% load data
load(['data_sim_', name, '.mat'])
load(['Y_pred_', model, '_', name, '_', num2str(seed), '.mat'])

% determine Matern or needlet
if strcmp(model, 'Matern')
    Y_pred = Y_pred_Matern;
else
    Y_pred = Y_pred_needlet;
end

theta_samples = theta_vec(index);
phi_samples = phi_vec(index);
cmin = min(min(Y), min(Y_pred));
cmax = max(max(Y), max(Y_pred));

% plot1
figure
subplot('position', [0 0.55 0.9 0.35])
pcolor(phi_mat, theta_mat, reshape(Y, size(phi_mat)))
shading flat
axis image
hold on
plot([pi-pi/18 pi-pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
plot([pi+pi/18 pi+pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
scatter(phi_samples, theta_samples, 'k.')
caxis([cmin cmax])
xlabel('Longitude')
ylabel('Colatitdue')
title('Data')

subplot('position', [0 0.1 0.9 0.35])
pcolor(phi_mat, theta_mat, reshape(Y_pred, size(phi_mat)))
shading flat
axis image
hold on
plot([pi-pi/18 pi-pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
plot([pi+pi/18 pi+pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
caxis([cmin cmax])
xlabel('Longitude')
ylabel('Colatitdue')
if strcmp(model, 'Matern')
    title('Predictions using non-stationary Matern model')
else
    title('Predictions using needlet model')
end

% colorbar
h = colorbar;
set(h, 'Position', [.85 .1 .05 .8]);

% plot2
figure
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[0.1 0.2 0.85 0.7]);
err = Y-Y_pred;
pcolor(phi_mat, theta_mat, reshape(err, size(phi_mat)))
shading flat
axis image
hold on
plot([pi-pi/18 pi-pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
plot([pi+pi/18 pi+pi/18], [0 pi], 'k--', 'LineWidth', 1.5)
scatter(phi_samples, theta_samples, 'k.')
caxis([-max(abs(err)) max(abs(err))])
colorbar
xlabel('Longitude')
ylabel('Colatitude')

% add text
axes(ax1)
text(0.6, 0.1, ['max-err: ', num2str(max(err)), '    min-err: ', num2str(min(err))])

end