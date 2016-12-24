function [MAE_in, MAE_out, MSPE_in, MSPE_out, CRPS_in, CRPS_out] = combine_results(name)
% combine results from the outputs of sim_pred

MAE_in = [];
MAE_out = [];
MSPE_in = [];
MSPE_out = [];
CRPS_in = [];
CRPS_out = [];

methods = {'nonGau_needlet', 'Gau_needlet', 'Matern'};
    
for i = 1:3
    load(['sim_pred_err_', methods{i}, '_', name, '.mat'])
    MAE_in = [MAE_in MAE_in_all];
    MAE_out = [MAE_out MAE_out_all];
    MSPE_in = [MSPE_in MSPE_in_all];
    MSPE_out = [MSPE_out MSPE_out_all];
    CRPS_in = [CRPS_in CRPS_in_all];
    CRPS_out = [CRPS_out CRPS_out_all];
end

end