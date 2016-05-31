load('sim_energy.mat')

load('WHI_quad_cond.mat')
P_cond = all_Cond_N{1}(:);

plot_pot(reshape(P_cond, size(phi)), phi, theta, 1000, max(abs(P_cond)))
P_cond = P_cond(361:(end-360));

R = 6.5*1e6;
energy = relative_energy/R^2;
for i = 1:size(energy, 2)
    energy(:, i) = energy(:, i).*P_cond;
end

figure
for i = 1:25
    subplot(5, 5, i)
    ksdensity(energy(i*500, :), 'support', 'positive')
end

figure
for i = 1:25
    subplot(5, 5, i)
    hist(energy(i*500, :), 100)
end
