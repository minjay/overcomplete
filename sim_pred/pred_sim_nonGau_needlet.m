function [MSPE_out, MSPE_in, MAE_out, MAE_in, CRPS_out, CRPS_in,...
    avg_len_90_out, avg_len_90_in, avg_len_50_out, avg_len_50_in,...
    cp_90_out, cp_90_in, cp_50_out, cp_50_in, QS_95_out, QS_95_in,...
    QS_05_out, QS_05_in] = ...
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

% get quantiles
Y_lb_90 = quantile(Y_pred_all, 0.05, 2);
Y_ub_90 = quantile(Y_pred_all, 0.95, 2);
Y_lb_50 = quantile(Y_pred_all, 0.25, 2);
Y_ub_50 = quantile(Y_pred_all, 0.75, 2);

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
mu_all = std_vec_all(index_pred, :).*Ac(index_pred, :);
sigma_sq_all = 1./post_samples.tau_sq_inv;
n_index_pred = length(index_pred);
% create the full array
CRPS_all = zeros(N, 1);

for i = 1:n_index_pred
    term2 = CRPS_term2(mu_all(i, :), sigma_sq_all);
    term1 = mean(fun_A(Y(index_pred(i))-mu_all(i, :), sigma_sq_all));
    CRPS_all(index_pred(i)) = term1-term2;
end

CRPS_out = mean(CRPS_all(index_pred_out));
CRPS_in = mean(CRPS_all(index_pred_in));

% PI
avg_len_90 = Y_ub_90-Y_lb_90;
avg_len_90_out = mean(avg_len_90(index_pred_out));
avg_len_90_in = mean(avg_len_90(index_pred_in));
avg_len_50 = Y_ub_50-Y_lb_50;
avg_len_50_out = mean(avg_len_50(index_pred_out));
avg_len_50_in = mean(avg_len_50(index_pred_in));
cp_90 = Y>=Y_lb_90 & Y<=Y_ub_90;
cp_90_out = mean(cp_90(index_pred_out));
cp_90_in = mean(cp_90(index_pred_in));
cp_50 = Y>=Y_lb_50 & Y<=Y_ub_50;
cp_50_out = mean(cp_50(index_pred_out));
cp_50_in = mean(cp_50(index_pred_in));

% quantiles
QS_95_all = zeros(N, 1);
QS_05_all = zeros(N, 1);
for i = 1:n_index_pred
    QS_95_all(index_pred(i)) = QS(0.95, Y(index_pred(i)), Y_ub_90(index_pred(i)));
    QS_05_all(index_pred(i)) = QS(0.05, Y(index_pred(i)), Y_lb_90(index_pred(i)));
end
QS_95_out = mean(QS_95_all(index_pred_out));
QS_95_in = mean(QS_95_all(index_pred_in));
QS_05_out = mean(QS_05_all(index_pred_out));
QS_05_in = mean(QS_05_all(index_pred_in));

end
