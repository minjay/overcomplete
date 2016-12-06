addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

load('data_EOF_regr_new.mat')
% work on the data at the first time point
resid = resid_all(1, :);

rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec);
[pot_samples, index] = datasample(resid', 4000, 'Replace', false,...
    'Weights', w);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

% plot
% plot_samples(theta_vec, index, phi_samples, pot_samples)

% fit
% parameter specification
B = 2;
j_min = 2;
j_max = 4;
nu = 4;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples, phi_samples);
M = size(A, 2);

% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 1 1 1 1]*pi/4;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_samples);

% replace the first B-spline basis function by intercept
b_mat(:, 1) = 1;

r = size(b_mat, 2)-1;

% rescale the observations
Y = pot_samples/1e3;

% get init values for MCMC
beta_init = [zeros(1, r+1) 0.1^2 0.01^2 1e-2];
negloglik1 = @(beta_all) negloglik_Gaussian_needlet(beta_all, b_mat, Y, Npix, A);

lb = [-10*ones(1, r+1) 0 0 1e-3];
ub = [10*ones(1, r+1) 1 1 Inf];

[beta_hat, f_min] = Gaussian_needlet_fit(negloglik1, beta_init, lb, ub, true);

save('post_samples_real_exp.mat', 'beta_hat')