% assess predictive performance

function [MSPE, MAE, CRPS, QS_95, QS_05] = post_assess_eta0(name)

% load fitted result
load(['post_samples_real_reparam_nu', name, '.mat'])

% load data
load('data_EOF_regr_new.mat')
Y = resid_all(1, :)';

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

% non-stationary variance function
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(N, 1) b_mat];

range = 2001:3000;
post_samples.eta = post_samples.eta(:, range);
post_samples.c = post_samples.c(:, range);
post_samples.tau_sq_inv = post_samples.tau_sq_inv(range);

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

Y_err = Y(index_pred)-Y_pred_needlet(index_pred);

% MSPE
MSPE = mean(Y_err.^2);

% MAE
MAE = mean(abs(Y_err));

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

CRPS = mean(CRPS_all(index_pred));

% quantiles
QS_95_all = zeros(N, 1);
QS_05_all = zeros(N, 1);
for i = 1:n_index_pred
    QS_95_all(index_pred(i)) = QS(0.95, Y(index_pred(i)), Y_ub_90(index_pred(i)));
    QS_05_all(index_pred(i)) = QS(0.05, Y(index_pred(i)), Y_lb_90(index_pred(i)));
end
QS_95 = mean(QS_95_all(index_pred));
QS_05 = mean(QS_05_all(index_pred));

end
