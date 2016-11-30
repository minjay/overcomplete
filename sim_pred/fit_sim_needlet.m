function [post_samples, Npix, index, index_region] = fit_sim_needlet(seed, flag, name, width)
% seed is the random seed used to sample data
% flag determines whether to save file
% name is one of 2dot5, 3 or 4.
% width is the width of the region

load(['data_sim_', name, '.mat'])
rng(seed)

% sampling
[index, index_region] = rand_sampler(phi, width);
theta_samples = theta(index);
phi_samples = phi(index);
Y = Y(index);

B = 2;
j_min = 2;
j_max = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples, phi_samples);
M = size(A, 2); 

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_samples);

b_mat(:, 1) = 1;

r = size(b_mat, 2)-1;

% get init values for MCMC
beta_init = [zeros(1, r+1) 0.1^2 1e-2];
negloglik1 = @(beta_all) negloglik_Gaussian_needlet(beta_all, b_mat, Y, Npix, A);

lb = [-10*ones(1, r+1) 0 1e-3];
ub = [10*ones(1, r+1) 1 Inf];

[beta_hat, f_min] = Gaussian_needlet_fit(negloglik1, beta_init, lb, ub, true);

% init
% c
c_init = zeros(M, 1);
% V
V_inv_init = ones(M, 1); 
% sigma_j_sq
sigma_j_sq_init = beta_hat(end-1);
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
lambda = 0.05;
% the number of MCMC iterations
T = 3e5;
% the length of the burn-in period
burn_in = 1e5;
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

options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report);

post_samples = Gibbs_sampler_AM_rep_inter(model, data, params, tuning, options);

if flag
    save(['post_samples_', name, '_', num2str(seed), '.mat'], 'post_samples', 'Npix', 'index', 'index_region')
end

end