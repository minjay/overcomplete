% assess predictive performance

clear

% load fitted result
load('post_samples_real_exp3.mat')

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

st = 2001;
en = 3000;
T = en-st+1;
Ac = A*post_samples.c(:, st:en);
Y_pred_all = zeros(N, T);
std_vec_all = exp(b_mat*post_samples.eta(:, st:en));
err_all = randn(1, T).*(1./sqrt(post_samples.tau_sq_inv(st:en)));

for t = 1:T
    Y_pred_all(:, t) = (std_vec_all(:, t).*Ac(:, t)+err_all(t))*1e3;
end

Y_pred_need = mean(Y_pred_all, 2);

figure
Y_err_need = resid'-Y_pred_need;
plot_pot_with_obs(reshape(Y_err_need, size(phi)), phi, theta, phi_samples, theta_samples, 1000)

save('Y_pred_need.mat', 'Y_pred_need', 'Y_err_need')
