name = '4';

% 20 times
T = 20;

% width of the region
width = pi/9;

MSPE_Matern_out_region_all = zeros(T, 1);
MSPE_Matern_region_all = zeros(T, 1);
MAE_Matern_out_region_all = zeros(T, 1);
MAE_Matern_region_all = zeros(T, 1);
for t = 1:T
    t
    % no save
    [beta_hat, index, index_region] = fit_sim_nonsta_Matern(t, 0, name, width);
    [MSPE_Matern_out_region, MSPE_Matern_region, MAE_Matern_out_region, MAE_Matern_region] = ...
        pred_sim_nonsta_Matern(name, index, index_region);
    MSPE_Matern_out_region_all(t) = MSPE_Matern_out_region;
    MSPE_Matern_region_all(t) = MSPE_Matern_region;
    MAE_Matern_out_region_all(t) = MAE_Matern_out_region;
    MAE_Matern_region_all(t) = MAE_needlet_region;
end

save(['sim_pred_err_Matern_rep_', name, '.mat'], 'MSPE_Matern_out_region_all', 'MSPE_Matern_region_all',...
    'MAE_Matern_out_region_all', 'MAE_Matern_region_all')