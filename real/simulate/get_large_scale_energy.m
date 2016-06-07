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
E_theta_large_scale = -FX/R;
E_phi_large_scale = -FY./(R*sin(theta));

save('energy_large_scale.mat', 'E_theta_large_scale', 'E_phi_large_scale')
