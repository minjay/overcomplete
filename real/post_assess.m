% assess predictive performance

clear

% load fitted result
load('post_samples_real_exp3_nu2dot5.mat')

% load data
load('data_EOF_regr_new.mat')
resid = resid_all(1, :);

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

theta_vec = theta(:);

% non-stationary variance function
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(N, 1) b_mat];

st = 3001;
en = 5000;
T = (en-st+1)/2;
Ac = A*post_samples.c(:, st:2:en);
Y_pred_all = zeros(N, T);
std_vec_all = exp(b_mat*post_samples.eta(:, st:2:en));
rng(1)
err_all = randn(1, T).*(1./sqrt(post_samples.tau_sq_inv(st:2:en)));

for t = 1:T
    Y_pred_all(:, t) = (std_vec_all(:, t).*Ac(:, t)+err_all(t))*1e3;
end

Y_pred_need = mean(Y_pred_all, 2);

% get quantiles
index_pred = setdiff(1:size(A, 1), index);
Y_lb_90 = quantile(Y_pred_all(index_pred, :), 0.05, 2);
Y_ub_90 = quantile(Y_pred_all(index_pred, :), 0.95, 2);
Y = resid(index_pred)';
n_index_pred = length(index_pred);
% quantiles
QS_95_all = zeros(n_index_pred, 1);
QS_05_all = zeros(n_index_pred, 1);
for i = 1:n_index_pred
    QS_95_all(i) = QS(0.95, Y(i), Y_ub_90(i));
    QS_05_all(i) = QS(0.05, Y(i), Y_lb_90(i));
end

mean((Y-Y_pred_need(index_pred)).^2)
mean(abs(Y-Y_pred_need(index_pred)))
max(Y-Y_pred_need(index_pred))
mean(QS_95_all)
mean(QS_05_all)

figure
Y_err_need = resid'-Y_pred_need;
plot_pot_with_obs(reshape(Y_err_need, size(phi)), phi, theta, phi_samples, theta_samples, 1000)

save('Y_pred_need.mat', 'Y_pred_need', 'Y_err_need')
