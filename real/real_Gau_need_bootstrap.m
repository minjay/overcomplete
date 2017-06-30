clear
clc

addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

load('data_EOF_regr_new.mat')
resid = resid_all(1, :);
load('post_samples_exp2.mat')

rng(1)

% sampling
n = 4000;
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
[pot_samples, index] = datasample(resid', n, 'Replace', false,...
    'Weights', w);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

% fit
% parameter specification
B = 2;
j_min = 2;
j_max = 4;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples*4, phi_samples);
M = size(A, 2);

% non-stationary variance function
% the first column is all ones
load('ns.mat')
b_mat_full = kron(b_mat, ones(size(theta, 1), 1));
b_mat = b_mat_full(index, :);
b_mat = [ones(n, 1) b_mat];

r = size(b_mat, 2)-1;

beta = beta_hat(1:end-1);
tau = beta_hat(end);
cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A)+tau^2*eye(n);
R = 200;
Y = mvnrnd(zeros(1, n), cov_mat, R);

lb = [-2*ones(1, r) 0 0 0 1e-3];
ub = [2*ones(1, r) 1 1 1 Inf];

beta_fit_all = zeros(R, length(beta_hat));

parfor i = 1:R
    negloglik1 = @(beta_all) negloglik_Gaussian_needlet(beta_all, b_mat, Y(i, :), Npix, A);
    [beta_fit, f_min] = Gaussian_needlet_fit(negloglik1, beta_hat, lb, ub, false);
    beta_fit_all(i, :) = beta_fit;
end

save('Gau_need_bootstrap.mat', 'beta_fit_all', 'beta_hat')
