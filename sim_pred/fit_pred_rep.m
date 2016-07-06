% 20 times
T = 20;

% width of the region
width = pi/9;

MSPE_needlet_out_region_all = zeros(T, 1);
MSPE_needlet_region_all = zeros(T, 1);
MAE_needlet_out_region_all = zeros(T, 1);
MAE_needlet_region_all = zeros(T, 1);
for t = 1:T
    t
    % no save
    [post_samples, ~, index, index_region] = fit_sim_needlet(t, 0, name, width);
    [MSPE_needlet_out_region, MSPE_needlet_region, MAE_needlet_out_region, MAE_needlet_region] = ...
        pred_sim_needlet(name, post_samples, index, index_region);
    MSPE_needlet_out_region_all(t) = MSPE_needlet_out_region;
    MSPE_needlet_region_all(t) = MSPE_needlet_region;
    MAE_needlet_out_region_all(t) = MAE_needlet_out_region;
    MAE_needlet_region_all(t) = MAE_needlet_region;
end

save(['sim_pred_err_rep_', name, '.mat'], 'MSPE_needlet_out_region_all', 'MSPE_needlet_region_all',...
    'MAE_needlet_out_region_all', 'MAE_needlet_region_all')