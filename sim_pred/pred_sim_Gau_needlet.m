function [MSPE_out, MSPE_in, MAE_out, MAE_in, CRPS_out, CRPS_in] = ...
    pred_sim_Gau_needlet(name, beta_hat, index, index_region)

load(['data_sim_', name, '.mat'])

% sampling
N = length(Y);
n = length(index);
Y_samples = Y(index);

B = 2;
j_min = 2;
j_max = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(Y_samples, n, 1);

index_pred = setdiff(1:N, index);
index_pred_out = setdiff(index_pred, index_region);
index_pred_in = index_region;
SigmaP0 = cov_mat(:, index);
Y_pred_needlet = SigmaP0*tmp;

SigmaPP = cov_mat;
Sigma0P = SigmaP0';
Var_Y_pred_needlet = diag(SigmaPP-SigmaP0*(Sigma00\Sigma0P));

Y_err_out = Y(index_pred_out)-Y_pred_needlet(index_pred_out);
Y_err_in = Y(index_pred_in)-Y_pred_needlet(index_pred_in);

% MSPE
MSPE_out = mean(Y_err_out.^2);
MSPE_in = mean(Y_err_in.^2);

% MAE
MAE_out = mean(abs(Y_err_out));
MAE_in = mean(abs(Y_err_in));

% CRPS
CRPS_out = mean(CRPS(Y(index_pred_out), Y_pred_needlet(index_pred_out),...
    Var_Y_pred_needlet(index_pred_out)));
CRPS_in = mean(CRPS(Y(index_pred_in), Y_pred_needlet(index_pred_in),...
    Var_Y_pred_needlet(index_pred_in)));

end
