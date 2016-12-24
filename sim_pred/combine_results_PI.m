function [len_50_in, len_50_out, len_90_in, len_90_out,...
    cp_50_in, cp_50_out, cp_90_in, cp_90_out] = combine_results_PI(name)
% combine results from the outputs of sim_pred
% PI

len_50_in = [];
len_50_out = [];
len_90_in = [];
len_90_out = [];
cp_50_in = [];
cp_50_out = [];
cp_90_in = [];
cp_90_out = [];

methods = {'nonGau_needlet', 'Gau_needlet', 'Matern'};
    
for i = 1:3
    load(['sim_pred_err_', methods{i}, '_', name, '.mat'])
    len_50_in = [len_50_in len_50_in_all];
    len_50_out = [len_50_out len_50_out_all];
    len_90_in = [len_90_in len_90_in_all];
    len_90_out = [len_90_out len_90_out_all];
    cp_50_in = [cp_50_in cp_50_in_all];
    cp_50_out = [cp_50_out cp_50_out_all];
    cp_90_in = [cp_90_in cp_90_in_all];
    cp_90_out = [cp_90_out cp_90_out_all];
end

end