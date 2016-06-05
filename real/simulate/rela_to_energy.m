load('sim_energy.mat')
load('sim_energy_Gau_need.mat')

load('WHI_quad_cond.mat')
P_cond = all_Cond_N{1}(:);

P_cond = P_cond(361:(end-360));

R = 6.5*1e6;
energy = relative_energy/R^2;
for i = 1:size(energy, 2)
    energy(:, i) = energy(:, i).*P_cond;
end

energy_Gau_need = relative_energy_Gau_need/R^2;
for i = 1:size(energy_Gau_need, 2)
    energy_Gau_need(:, i) = energy_Gau_need(:, i).*P_cond;
end

save('energy_both.mat', 'energy', 'energy_Gau_need', '-v7.3')