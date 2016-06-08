clear

load('Y_pred_need.mat')
load('Y_pred_Gau_need.mat')
load('Y_pred_Matern.mat')
load('theta_phi_R.mat')

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.15]);

phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);

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
xlabel(sprintf('Min %6.1f  Max %5.1f [V]',vmin,vmax),'FontName','times','Fontsize',10)

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
xlabel(sprintf('Min %6.1f  Max %5.1f [V]',vmin,vmax),'FontName','times','Fontsize',10)

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
xlabel(sprintf('Min %6.1f  Max %5.1f [V]',vmin,vmax),'FontName','times','Fontsize',10)

h = colorbar;
set(h, 'Position', [.9 .05 .025 .9]);