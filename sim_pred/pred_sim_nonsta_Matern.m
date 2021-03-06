function [MSPE_out, MSPE_in, MAE_out, MAE_in, CRPS_out, CRPS_in,...
    avg_len_90_out, avg_len_90_in, avg_len_50_out, avg_len_50_in,...
    cp_90_out, cp_90_in, cp_50_out, cp_50_in, QS_95_out, QS_95_in,...
    QS_05_out, QS_05_in] = ...
    pred_sim_nonsta_Matern(name, beta_hat, index, index_region)

load(['data_sim_', name, '.mat'])

% sampling
N = length(Y);
n = length(index);
Y_samples = Y(index);

[x, y, z] = trans_coord(theta, phi);

% get distance matrix
r_dist = get_chordal_dist(x, y, z);

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat = get_cov_nonsta_Matern(beta, r_dist, b_mat)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(Y_samples, n, 1);

index_pred = setdiff(1:N, index);
index_pred_out = setdiff(index_pred, index_region);
index_pred_in = index_region;
SigmaP0 = cov_mat(:, index);
Y_pred_Matern = SigmaP0*tmp;

SigmaPP = cov_mat;
Sigma0P = SigmaP0';
Var_Y_pred_Matern = diag(SigmaPP-SigmaP0*(Sigma00\Sigma0P));

Y_err_out = Y(index_pred_out)-Y_pred_Matern(index_pred_out);
Y_err_in = Y(index_pred_in)-Y_pred_Matern(index_pred_in);

% get quantiles
Y_lb_90 = Y_pred_Matern-norminv(0.95)*sqrt(Var_Y_pred_Matern);
Y_ub_90 = Y_pred_Matern+norminv(0.95)*sqrt(Var_Y_pred_Matern);
Y_lb_50 =  Y_pred_Matern-norminv(0.75)*sqrt(Var_Y_pred_Matern);
Y_ub_50 = Y_pred_Matern+norminv(0.75)*sqrt(Var_Y_pred_Matern);

% MSPE
MSPE_out = mean(Y_err_out.^2);
MSPE_in = mean(Y_err_in.^2);

% MAE
MAE_out = mean(abs(Y_err_out));
MAE_in = mean(abs(Y_err_in));

% CRPS
CRPS_out = mean(CRPS(Y(index_pred_out), Y_pred_Matern(index_pred_out),...
    Var_Y_pred_Matern(index_pred_out)));
CRPS_in = mean(CRPS(Y(index_pred_in), Y_pred_Matern(index_pred_in),...
    Var_Y_pred_Matern(index_pred_in)));

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
n_index_pred = length(index_pred);
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
