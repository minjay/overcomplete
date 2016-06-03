load('sim_energy.mat')
load('sim_energy_Gau_need.mat')

load('WHI_quad_cond.mat')
P_cond = all_Cond_N{1}(:);

plot_pot(reshape(P_cond, size(phi)), phi, theta, 1000, max(abs(P_cond)))
P_cond = P_cond(361:(end-360));
P_cond = P_cond(1:5:end);

R = 6.5*1e6;
energy = relative_energy/R^2;
for i = 1:size(energy, 2)
    energy(:, i) = energy(:, i).*P_cond;
end

energy_Gau_need = relative_energy_Gau_need/R^2;
for i = 1:size(energy_Gau_need, 2)
    energy_Gau_need(:, i) = energy_Gau_need(:, i).*P_cond;
end

% set up colormap
map = brewermap(2, 'Set1');

% figure
% for i = 1
%     t = 5000;
%     energy_norm = energy(t, :)/std(energy(t, :));
%     [f, xi] = ksdensity(energy_norm);
%     plot(xi, f)
%     hold on
%     energy_Gau_need_norm = energy_Gau_need(t, :)/std(energy_Gau_need(t, :));
%     [f, xi] = ksdensity(energy_Gau_need_norm);
%     plot(xi, f, 'r')
%     axis tight
% end

% set up subplot
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.05], [0.05 0.05], [0.05 0.01]);

figure
for i = 1:9
    subplot(3, 3, i)
    t = i*300;
    % normalize the energy
    energy_norm = energy(t, :)/std(energy(t, :));
    % set up bins
    xs = linspace(0, max(energy_norm), 200);
    histf(energy_norm, xs,...
        'facecolor', map(1, :), 'facealpha', .5, 'edgecolor', 'none')
    hold on
    % normalize the energy
    energy_Gau_need_norm = energy_Gau_need(t, :)/std(energy_Gau_need(t, :));
    histf(energy_Gau_need_norm, xs,...
        'facecolor', map(2, :), 'facealpha', .5, 'edgecolor', 'none')
    axis tight
    % only show the left half of the plot
    xlim([0 max(energy_norm)/2])
    % enlarge font size
    set(gca, 'FontSize', 12)
    box off
    % add legend
    if i==9 
        legalpha('nonGau-need','Gau-need','location', 'northeast')
        legend boxoff
    end
    title(['Location ', num2str(i)]) 
end

% save figure
% 2.5 times resolution
export_fig cmp_energy_dist.png -m2.5

% 
%  [f, xi] = ksdensity(neg_Y_theta(10000, :)/std(neg_Y_theta(10000, :)));
%  plot(xi, f)
%     hold on
%  [f, xi] = ksdensity(neg_Y_theta_Gau_need(10000, :)/std(neg_Y_theta_Gau_need(10000, :)), xi);
%  plot(xi, f, 'r')
