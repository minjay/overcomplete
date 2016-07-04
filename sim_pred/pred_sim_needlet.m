function [MSPE_needlet_out_region, MSPE_needlet_region, MAE_needlet_out_region, MAE_needlet_region] = ...
    pred_sim_needlet(name, post_samples, index, index_region)

load(['data_sim_', name, '.mat'])

% specify parameters
B = 2;
j_min = 2;
j_max = 3;

% design matrix A
[~, ~, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec);
N = size(A, 1);

% non-stationary variance function
knots = [0 0 0 0 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec);

b_mat(:, 1) = 1;

% predict
T = size(post_samples.eta, 2);
eta_mean = mean(post_samples.eta, 2);
tau_mean = mean(1./sqrt(post_samples.tau_sq_inv));
std_vec = exp(b_mat*eta_mean);
Ac = A*post_samples.c;
Y_pred_all = zeros(N, T);
for t = 1:T
    Y_pred_all(:, t) = std_vec.*Ac(:, t);
end
Y_pred_all = Y_pred_all+tau_mean*randn(N, T);

Y_pred_needlet = mean(Y_pred_all, 2);

index_pred = setdiff(1:N, index);
index_pred_out_region = setdiff(index_pred, index_region);
index_pred_region = index_region;

Y_err_needlet_out_region = Y(index_pred_out_region)-Y_pred_needlet(index_pred_out_region);
Y_err_needlet_region = Y(index_pred_region)-Y_pred_needlet(index_pred_region);

% MSPE
MSPE_needlet_out_region = mean(Y_err_needlet_out_region.^2);
MSPE_needlet_region = mean(Y_err_needlet_region.^2);

% MAE
MAE_needlet_out_region = mean(abs(Y_err_needlet_out_region));
MAE_needlet_region = mean(abs(Y_err_needlet_region));

end
