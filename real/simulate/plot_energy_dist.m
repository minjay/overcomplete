load('energy_both.mat')
load('energy_large_scale.mat')
load('theta_phi_R.mat')

cf = reshape(energy_large_scale, size(phi));
vmag = linspace(min(energy_large_scale), max(energy_large_scale), 10);
phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
h = colorbar;
set(h, 'Position', [.85 .1 .05 .8]);

% set up subplot
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.2]);

cmax = max(max(energy(:, 1:9)));
for i = 1:9
    subplot(3, 3, i)
    cf = reshape(energy(:, i), size(phi));
    vmag = linspace(min(cf(:)), max(cf(:)), 10);
    h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
    delete(h)
    shading flat
    caxis([0 cmax])
end
h = colorbar;
set(h, 'Position', [.85 .05 .05 .9]);

cmax = max(max(energy_Gau_need(:, 1:9)));
for i = 1:9
    subplot(3, 3, i)
    cf = reshape(energy_Gau_need(:, i), size(phi));
    vmag = linspace(min(cf(:)), max(cf(:)), 10);
    h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
    delete(h)
    shading flat
    caxis([0 cmax])
end
h = colorbar;
set(h, 'Position', [.85 .05 .05 .9]);
