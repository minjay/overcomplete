clear

load('WHI_quad.mat')
% reload theta, phi
load('data_EOF_regr_new.mat')

whole_field = all_Pot_N{1}(:);
whole_field = whole_field(361:(end-360));
resid = resid_all(1, :)';
large_scale = whole_field-resid;

phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.1], [0.2 0.05], [0.05 0.05]);

cmax = max(max(abs(whole_field/1e3)), max(abs(large_scale/1e3)));

subplot(1, 3, 1)
cf = reshape(whole_field/1e3, size(phi));
vmag = linspace(min(cf(:)), max(cf(:)), 10);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
caxis([-cmax cmax])
colormap(jet)
title('Sum')
text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)

subplot(1, 3, 2)
cf = reshape(large_scale/1e3, size(phi));
vmag = linspace(min(cf(:)), max(cf(:)), 10);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
caxis([-cmax cmax])
colormap(jet)
title('Large-scale')
text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)

h = colorbar('southoutside');
set(h, 'Position', [.05 .075 .575 .075]);

cmax = max(abs(resid/1e3));
subplot(1, 3, 3)
cf = reshape(resid/1e3, size(phi));
vmag = linspace(min(cf(:)), max(cf(:)), 10);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
caxis([-cmax cmax])
colormap(jet)
title('Resid')
text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)

h = colorbar('southoutside');
set(h, 'Position', [.7125 .075 .2375 .075]);
