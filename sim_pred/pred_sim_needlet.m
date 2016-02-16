function pred_sim_needlet(seed, flag, name, extra)

load(['data_sim_', name, '.mat'])
load(['post_samples_', name, '_', num2str(seed), extra, '.mat']) 

% specify parameters
B = 2;
j_min = 2;
j_max = 3;

% design matrix A
[~, ~, A] = get_A_ss(B, j_min, j_max, theta_vec, phi_vec);
N = size(A, 1);

% non-stationary variance function
m = 4;
lambda_inv = 2.5;
b_mat = get_nonsta_var(m, lambda_inv, theta_vec);

% predict
T = size(post_samples.eta, 2);
eta_mean = mean(post_samples.eta, 2);
std_vec = exp(b_mat*eta_mean);
Ac = A*post_samples.c;
Y_pred_all = zeros(N, T);
for t = 1:T
    Y_pred_all(:, t) = std_vec.*Ac(:, t);
end

Y_pred_needlet = mean(Y_pred_all, 2);

index_pred = setdiff(1:N, index);

Y_err_needlet = Y(index_pred)-Y_pred_needlet(index_pred);
Y_err_needlet_in = Y(index)-Y_pred_needlet(index);

% MSPE
MSPE_needlet = mean(Y_err_needlet.^2);
MSPE_needlet_in = mean(Y_err_needlet_in.^2);
fprintf('Out-of-sample MSPE of needlet is %5f\n', MSPE_needlet)
fprintf('In-sample MSPE of needlet is %5f\n', MSPE_needlet_in)

% MAE
MAE_needlet = mean(abs(Y_err_needlet));
MAE_needlet_in = mean(abs(Y_err_needlet_in));
fprintf('Out-of-sample MAE of needlet is %5f\n', MAE_needlet)
fprintf('In-sample MAE of needlet is %5f\n', MAE_needlet_in)

if flag
    save(['Y_pred_needlet_', name, '_', num2str(seed), extra, '.mat'],...
        'Y_pred_needlet', 'Y_pred_all', 'Y', 'index', 'index_pred',...
        'MSPE_needlet', 'MAE_needlet')
end

end
