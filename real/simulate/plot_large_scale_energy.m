load('WHI_quad.mat')
load('WHI_quad_cond.mat')
% reload theta, phi
% note the inconsistency between WHI_quad and WHI_quad_cond
load('data_EOF_regr_new.mat')

for t = 1:720
    t_new = (t-1)*10+1;
    whole_field = all_Pot_N{t_new}(:);
    whole_field = whole_field(361:(end-360));
    resid = resid_all(t, :)';
    large_scale = whole_field-resid;
    large_scale = reshape(large_scale, size(phi));

    HX = theta(1, :);
    HY = phi(:, 1);
    [FX, FY] = gradient(large_scale, HX, HY);
    R = 6.5*1e6;
    FX = -FX/R;
    FY = -FY./(R*sin(theta));

    energy_large_scale = FX.^2+FY.^2;

    P_cond = all_Cond_N{t_new}(:);

    P_cond = P_cond(361:(end-360));
    energy_large_scale = energy_large_scale(:).*P_cond;

    plot_pot(reshape(energy_large_scale, size(phi)), phi, theta, 1000, max(abs(energy_large_scale)))
    print(['./plots/', 'fig', num2str(t)], '-dpng')
    close
end