function pred_err_Matern(name)
% Evaluates prediction errors of nonstationary Matern model
% Criteria: MSPE, MAE, CRPS

T = 20;

MSPE_out_all = zeros(T, 1);
MSPE_in_all = zeros(T, 1);
MAE_out_all = zeros(T, 1);
MAE_in_all = zeros(T, 1);
CRPS_out_all = zeros(T, 1);
CRPS_in_all = zeros(T, 1);
len_90_out_all = zeros(T, 1);
len_90_in_all = zeros(T, 1);
len_50_out_all = zeros(T, 1);
len_50_in_all = zeros(T, 1);
cp_90_out_all = zeros(T, 1);
cp_90_in_all = zeros(T, 1);
cp_50_out_all = zeros(T, 1);
cp_50_in_all = zeros(T, 1);

for t = 1:T
    t
    load(['beta_hat_', name, '_', num2str(t), '.mat'])
    [MSPE_out, MSPE_in, MAE_out, MAE_in, CRPS_out, CRPS_in] = ...
        pred_sim_nonsta_Matern(name, beta_hat, index, index_region);
    MSPE_out_all(t) = MSPE_out;
    MSPE_in_all(t) = MSPE_in;
    MAE_out_all(t) = MAE_out;
    MAE_in_all(t) = MAE_in;
    CRPS_out_all(t) = CRPS_out;
    CRPS_in_all(t) = CRPS_in;
    len_90_out_all(t) = avg_len_90_out;
    len_90_in_all(t) = avg_len_90_in;
    len_50_out_all(t) = avg_len_50_out;
    len_50_in_all(t) = avg_len_50_in;
    cp_90_out_all(t) = cp_90_out;
    cp_90_in_all(t) = cp_90_in;
    cp_50_out_all(t) = cp_50_out;
    cp_50_in_all(t) = cp_50_in;
end

save(['sim_pred_err_Matern_', name, '.mat'], 'MSPE_out_all', 'MSPE_in_all',...
    'MAE_out_all', 'MAE_in_all', 'CRPS_out_all', 'CRPS_in_all',...
    'len_90_out_all', 'len_90_in_all', 'len_50_out_all', 'len_50_in_all',...
    'cp_90_out_all', 'cp_90_in_all', 'cp_50_out_all', 'cp_50_in_all')

end