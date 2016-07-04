function [MSPE_Matern_out_region, MSPE_Matern_region, MAE_Matern_out_region, MAE_Matern_region] = ...
    pred_sim_nonsta_Matern(name, index, index_region)

load(['data_sim_', name, '.mat'])

% sampling
N = length(Y);
n = length(index);
Y_samples = Y(index);

[x, y, z] = trans_coord(theta_vec, phi_vec);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance function
knots = [0 0 0 0 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec);

b_mat(:, 1) = 1;

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat = get_cov_nonsta_Matern(beta, r, b_mat)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(Y_samples, n, 1);

index_pred = setdiff(1:N, index);
index_pred_out_region = setdiff(index_pred, index_region);
index_pred_region = index_region;
SigmaP0 = cov_mat(:, index);
Y_pred_Matern = SigmaP0*tmp;

Y_err_Matern_out_region = Y(index_pred_out_region)-Y_pred_Matern(index_pred_out_region);
Y_err_Matern_region = Y(index_pred_region)-Y_pred_Matern(index_pred_region);

% MSPE
MSPE_Matern_out_region = mean(Y_err_Matern_out_region.^2);
MSPE_Matern_region = mean(Y_err_Matern_region.^2);

% MAE
MAE_Matern_out_region = mean(abs(Y_err_Matern_out_region));
MAE_Matern_region = mean(abs(Y_err_Matern_region));

end
