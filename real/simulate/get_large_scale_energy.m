load('WHI_quad.mat')
% reload theta, phi
load('data_EOF_regr_new.mat')

whole_field = all_Pot_N{1}(:);
whole_field = whole_field(361:(end-360));
resid = resid_all(1, :)';
large_scale = whole_field-resid;
large_scale = reshape(large_scale, size(phi));

plot_pot(large_scale, phi, theta, 1000, max(abs(large_scale(:))))

HX = theta(1, :);
HY = phi(:, 1);
[FX, FY] = gradient(large_scale, HX, HY);
R = 6.5*1e6;
FX = -FX/R;
FY = -FY./(R*sin(theta));

figure
plot_pot(FX, phi, theta, 1000, max(abs(FX(:))))
figure
plot_pot(FY, phi, theta, 1000, max(abs(FY(:))))

energy_large_scale = FX.^2+FY.^2;

load('WHI_quad_cond.mat')
P_cond = all_Cond_N{1}(:);

P_cond = P_cond(361:(end-360));
energy_large_scale = energy_large_scale(:).*P_cond;

save('energy_large_scale.mat', 'energy_large_scale')
