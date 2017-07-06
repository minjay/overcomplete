addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

clear

load('post_samples_real_reparam_nu4.mat')

% load data
load('data_EOF_regr_new.mat')
Y = resid_all(1, :)';
n = length(index);
Y_samples = Y(index);

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

% non-stationary variance function
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(N, 1) b_mat];

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(Y_samples, n, 1);

index_pred = setdiff(1:N, index);
SigmaP0 = cov_mat(:, index);
Y_pred_needlet = SigmaP0*tmp;

SigmaPP = cov_mat;
Sigma0P = SigmaP0';
Var_Y_pred_needlet = diag(SigmaPP-SigmaP0*(Sigma00\Sigma0P));

Y_err = Y(index_pred)-Y_pred_needlet(index_pred);

% get quantiles
Y_lb_90 = Y_pred_needlet-norminv(0.95)*sqrt(Var_Y_pred_needlet);
Y_ub_90 = Y_pred_needlet+norminv(0.95)*sqrt(Var_Y_pred_needlet);
Y_lb_50 =  Y_pred_needlet-norminv(0.75)*sqrt(Var_Y_pred_needlet);
Y_ub_50 = Y_pred_needlet+norminv(0.75)*sqrt(Var_Y_pred_needlet);

% MSPE
MSPE = mean(Y_err.^2);

% MAE
MAE = mean(abs(Y_err));

% CRPS
CRPS = mean(CRPS(Y(index_pred), Y_pred_needlet(index_pred),...
    Var_Y_pred_needlet(index_pred)));

% quantiles
n_index_pred = length(index_pred);
QS_95_all = zeros(N, 1);
QS_05_all = zeros(N, 1);
for i = 1:n_index_pred
    QS_95_all(index_pred(i)) = QS(0.95, Y(index_pred(i)), Y_ub_90(index_pred(i)));
    QS_05_all(index_pred(i)) = QS(0.05, Y(index_pred(i)), Y_lb_90(index_pred(i)));
end
QS_95 = mean(QS_95_all(index_pred));
QS_05 = mean(QS_05_all(index_pred));

save('post_assess_Gau.mat', 'MSPE', 'MAE', 'CRPS', 'QS_95', 'QS_05')
