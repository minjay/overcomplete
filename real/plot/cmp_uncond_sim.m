clear

% load data 
load('Y_sim_need.mat')
load('Y_sim_Matern.mat')
load('theta_phi_R.mat')

% rescale
Y_sim_need = Y_sim_need(5:8, :)/1e3;
Y_sim_Gau_need = Y_sim_Gau_need(5:8, :)/1e3;
Y_sim_Matern = Y_sim_Matern(5:8, :)/1e3;

figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.05], [0.05 0.05], [0.05 0.1]);

phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);

% plot nonGau-need
cmax = max(abs(Y_sim_need(:)));
for t = 1:4
    subplot(3, 4, t)
    cf = reshape(Y_sim_need(t, :), size(phi));
    vmag = linspace(min(cf(:)), max(cf(:)), 10);
    h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
    delete(h)
    shading flat
    caxis([-cmax cmax])
    text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
    text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)
    if t==1
        ylabel('AXING-need')
    end
end
h = colorbar;
colormap(jet)
set(h, 'Position', [.925 1-0.275 .025 .2]);

% plot Gau-need
cmax = max(abs(Y_sim_Gau_need(:)));
for t = 1:4
    subplot(3, 4, t+4)
    cf = reshape(Y_sim_Gau_need(t, :), size(phi));
    vmag = linspace(min(cf(:)), max(cf(:)), 10);
    h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
    delete(h)
    shading flat
    caxis([-cmax cmax])
    text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
    text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)
    if t==1
        ylabel('Gau-need')
    end
end
h = colorbar;
colormap(jet)
set(h, 'Position', [.925 1-0.3*2 .025 .2]);

% plot Gau-Matern
cmax = max(abs(Y_sim_Matern(:)));
for t = 1:4
    subplot(3, 4, t+8)
    cf = reshape(Y_sim_Matern(t, :), size(phi));
    vmag = linspace(min(cf(:)), max(cf(:)), 10);
    h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
    delete(h)
    shading flat
    caxis([-cmax cmax])
    text(-50, -50, sprintf('Min\n%2.1f',min(cf(:))),'FontName','times','Fontsize',10)
    text(30, -50, sprintf('Max\n%2.1f [kV]',max(cf(:))),'FontName','times','Fontsize',10)
    if t==1
        ylabel('Gau-Matern')
    end
end
h = colorbar;
colormap(jet)
set(h, 'Position', [.925 0.075 .025 .2]);

print -painters -depsc cmp_uncond_sim.eps