function [MSPE_out, MSPE_in, MAE_out, MAE_in, CRPS_out, CRPS_in] = ...
    pred_sim_nonGau_needlet(name, post_samples, index, index_region)

load(['data_sim_', name, '.mat'])

post_samples.c = squeeze(post_samples.c);

% specify parameters
B = 2;
j_min = 2;
j_max = 3;

% design matrix A
[~, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);
N = size(A, 1);

% predict
T = size(post_samples.eta, 2);
Ac = A*post_samples.c;
Y_pred_all = zeros(N, T);
std_vec_all = exp(b_mat*post_samples.eta);
err_all = randn(1, T).*(1./sqrt(post_samples.tau_sq_inv));
for t = 1:T
    Y_pred_all(:, t) = std_vec_all(:, t).*Ac(:, t)+err_all(t);
end

Y_pred_needlet = mean(Y_pred_all, 2);

index_pred = setdiff(1:N, index);
index_pred_out = setdiff(index_pred, index_region);
index_pred_in = index_region;

Y_err_out = Y(index_pred_out)-Y_pred_needlet(index_pred_out);
Y_err_in = Y(index_pred_in)-Y_pred_needlet(index_pred_in);

% MSPE
MSPE_out = mean(Y_err_out.^2);
MSPE_in = mean(Y_err_in.^2);

% MAE
MAE_out = mean(abs(Y_err_out));
MAE_in = mean(abs(Y_err_in));

% CRPS
CRPS_out_MC = zeros(length(index_pred_out), T);
CRPS_in_MC = zeros(length(index_pred_in), T);
for t = 1:T
    CRPS_out_MC(:, t) = CRPS(Y(index_pred_out),...
        std_vec_all(index_pred_out, t).*Ac(index_pred_out, t), 1/post_samples.tau_sq_inv(t));
    CRPS_in_MC(:, t) = CRPS(Y(index_pred_in),...
        std_vec_all(index_pred_in, t).*Ac(index_pred_in, t), 1/post_samples.tau_sq_inv(t));
end
CRPS_out = mean(mean(CRPS_out_MC));
CRPS_in = mean(mean(CRPS_in_MC));

end
