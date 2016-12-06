addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

load('data_exp.mat')

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

delete(gcp)

TT = size(Y, 2);

% init
% c
c_init = zeros(M, TT);
% V
V_inv_init = ones(M, 1); 
% sigma_j_sq
sigma_j_sq_init = beta_hat(end-2:end-1)';
% eta
eta_init = beta_hat(1:r+1)';
% pri_sig of eta_0
tau_sigma_sq = 1e2;
% pri_sig of eta
tau_eta_sq = 1e2;
% tau
tau_init = beta_hat(end);
tau_sq_inv_init = 1/tau_init^2;
% tuning parameters
mu_init = zeros(r+1, 1);
Sigma_init = eye(r+1);
lambda = 0.001;
% the number of MCMC iterations
T = 6e5+1e4;
% the length of the burn-in period
burn_in = 0;
% the length of the thinning interval
thin = 200;
% the length of the interval to report progress
n_report = 100;

model = struct('A', A, 'b_mat', b_mat, 'nu', nu);

data = struct('Y', Y, 'Npix', Npix);

params = struct('c', c_init, 'V', V_inv_init, 'sigma_j_sq', sigma_j_sq_init,...
    'eta', eta_init, 'tau_sigma_sq', tau_sigma_sq, 'tau_eta_sq', tau_eta_sq,...
    'tau', tau_sq_inv_init);

tuning = struct('mu', mu_init, 'Sigma', Sigma_init, 'lambda', lambda);

options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report, 'save', true);

post_samples = Gibbs_sampler_AM_rep_inter(model, data, params, tuning, options);

save('post_samples_real_exp.mat', 'post_samples', 'Npix', 'index', 'theta_samples', 'phi_samples', 'beta_hat')