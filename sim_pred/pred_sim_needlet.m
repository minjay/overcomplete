function pred_sim_needlet(seed, flag, name)

load(['data_sim_', name, '.mat'])
load(['post_samples_', name, '_', num2str(seed), '.mat']) 

% specify parameters
B = 2;
j_min = 2;
j_max = 4;

% design matrix A
[~, ~, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec);
N = size(A, 1);

% non-stationary variance function
m = 4;
lambda_inv = 2.5;
b_mat = get_nonsta_var(m, lambda_inv, theta_vec);

% predict
T = 1e3;
Ac = A*post_samples.c;
Y_pred_all = zeros(N, T);
for t = 1:T
    eta = post_samples.eta(:, t);
    tau = 1/sqrt(post_samples.tau_sq_inv(t));
    std_vec = exp(b_mat*eta);
    for i = 1:N
        Y_pred_all(i, :) = std_vec(i)*Ac(i, :);
    end
    Y_pred_all(:, t) = Y_pred_all(:, t)+tau^2*randn(N, 1);
end

Y_pred_needlet = mean(Y_pred_all, 2);

Y_err_needlet = Y(index_pred)-Y_pred_needlet(index_pred);
Y_err_needlet_in = Y(index)-Y_pred_needlet(index);

% MSPE
MSPE_needlet = mean(Y_err_needlet.^2);
MSPE_needlet_in = mean(Y_err_needlet_in.^2);
fprintf('Out-of-sample MSPE of Matern is %5f\n', MSPE_needlet)
fprintf('In-sample MSPE of Matern is %5f\n', MSPE_needlet_in)

% MAE
MAE_needlet = mean(abs(Y_err_needlet));
MAE_needlet_in = mean(abs(Y_err_needlet_in));
fprintf('Out-of-sample MAE of Matern is %5f\n', MAE_needlet)
fprintf('In-sample MAE of Matern is %5f\n', MAE_needlet_in)

if flag
    save(['Y_pred_needlet_', name, '_', num2str(seed), '.mat'],...
        'Y_pred_needlet', 'Y_pred_all', 'Y', 'index', 'index_pred',...
        'MSPE_needlet', 'MAE_needlet')
end

end
