clear

load('Y_pred_need.mat')
load('Y_pred_Gau_need.mat')
load('Y_pred_Matern.mat')
load('theta_phi_R.mat')

% get theta_samples and phi_samples
load('data_EOF_regr_new.mat')
resid = resid_all(10, :);

rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
[pot_samples, index] = datasample(resid', 4000, 'Replace', false,...
    'Weights', w);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.15]);

phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);

phi_samples_rot = phi_samples+pi/2;
[x_samples, y_samples] = pol2cart(phi_samples_rot, theta_samples/pi*180);

Y_err_need = Y_err_need/1e3;
Y_err_Gau_need = Y_err_Gau_need/1e3;
Y_err_Matern = Y_err_Matern/1e3;

cmax = max([max(abs(Y_err_need)) max(abs(Y_err_Gau_need)) max(abs(Y_err_Matern))]);

subplot(1, 3, 1)
cf = reshape(Y_err_need, size(phi));
vmin = min(cf(:));
vmax = max(cf(:));
vmag = linspace(vmin, vmax, 10);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
caxis([-cmax cmax])
title('nonGau-need')
xlabel(sprintf('Min %6.3f  Max %5.3f [kV]',vmin,vmax),'FontName','times','Fontsize',10)
hold on
scatter(x_samples(:), y_samples(:), 5, '.')

subplot(1, 3, 2)
cf = reshape(Y_err_Gau_need, size(phi));
vmin = min(cf(:));
vmax = max(cf(:));
vmag = linspace(vmin, vmax, 10);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
caxis([-cmax cmax])
title('Gau-need')
xlabel(sprintf('Min %6.3f  Max %5.3f [kV]',vmin,vmax),'FontName','times','Fontsize',10)
hold on
scatter(x_samples(:), y_samples(:), 5, '.')

subplot(1, 3, 3)
cf = reshape(Y_err_Matern, size(phi));
vmin = min(cf(:));
vmax = max(cf(:));
vmag = linspace(vmin, vmax, 10);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
caxis([-cmax cmax])
title('Gau-Matern')
xlabel(sprintf('Min %6.3f  Max %5.3f [kV]',vmin,vmax),'FontName','times','Fontsize',10)
hold on
scatter(x_samples(:), y_samples(:), 5, '.')

h = colorbar;
set(h, 'Position', [.9 .05 .025 .9]);