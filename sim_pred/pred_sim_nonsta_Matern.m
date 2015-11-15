function pred_sim_nonsta_Matern(seed, flag, name)

load(['data_sim_', name, '.mat'])
load(['beta_hat_', name, '_', num2str(seed), '.mat']) 

% sampling
N = length(Y);
n = length(index);
Y_samples = Y(index);

[x, y, z] = trans_coord(theta_vec, phi_vec);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance function
m = 4;
lambda_inv = 2.5;
b_mat = get_nonsta_var(m, lambda_inv, theta_vec);

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat = get_cov_nonsta_Matern(beta, r, b_mat)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(Y_samples, n, 1);

index_pred = setdiff(1:N, index);
SigmaP0 = cov_mat(:, index);
Y_pred_Matern = SigmaP0*tmp;

Y_err_Matern = Y(index_pred)-Y_pred_Matern(index_pred);
Y_err_Matern_in = Y(index)-Y_pred_Matern(index);

% MSPE
MSPE_Matern = mean(Y_err_Matern.^2);
MSPE_Matern_in = mean(Y_err_Matern_in.^2);
fprintf('Out-of-sample MSPE of Matern is %5f', MSPE_Matern)
fprintf('In-sample MSPE of Matern is %5f', MSPE_Matern_in)

% MAE
MAE_Matern = mean(abs(Y_err_Matern));
MAE_Matern_in = mean(abs(Y_err_Matern_in));
fprintf('Out-of-sample MAE of Matern is %5f', MAE_Matern)
fprintf('In-sample MAE of Matern is %5f', MAE_Matern_in)

if flag
    save(['Y_pred_Matern_', name, '_', num2str(seed), '.mat'],...
        'Y_pred_Matern', 'Y', 'index', 'index_pred', 'MSPE_Matern',...
        'MAE_Matern')
end

end
