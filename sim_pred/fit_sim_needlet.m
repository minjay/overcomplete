function [post_samples, Npix, index] = fit_sim_needlet(seed, flag, name)

load(['data_sim_', name, '.mat'])
rng(seed)

% sampling
index = rand_sampler(phi_vec);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);
Y = Y(index);

B = 2;
j_min = 2;
j_max = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples, phi_samples);
M = size(A, 2); 

% non-stationary variance function
m = 4;
lambda_inv = 2.5;
b_mat = get_nonsta_var(m, lambda_inv, theta_samples);

% init
% c
c_init = zeros(M, 1);
% V
V_inv_init = ones(M, 1); 
% eta
eta_init = zeros(m+1, 1);
% pri_sig of eta_0
tau_sigma_sq = 1e4;
% pri_sig of eta
tau_eta_sq = 0.25^2;
% tau
tau_init = 0.01;
tau_sq_inv_init = 1/tau_init^2;
% tuning parameters
mu_init = zeros(m+1, 1);
Sigma_init = eye(m+1);
lambda = 0.002;
% the number of MCMC iterations
T = 4e5;
% the length of the burn-in period
burn_in = 2e5;
% the length of the thinning interval
thin = 200;
% the length of the interval to report progress
n_report = 100;

model = struct('A', A, 'fj_sq', fj_sq, 'b_mat', b_mat, 'nu', nu);

data = struct('Y', Y, 'Npix', Npix);

params = struct('c', c_init, 'V', V_inv_init, 'eta', {eta_init, tau_sigma_sq, tau_eta_sq},...
    'tau', tau_sq_inv_init);

tuning = struct('mu', mu_init, 'Sigma', Sigma_init, 'lambda', lambda);

options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report);

post_samples = Gibbs_sampler_AM(model, data, params, tuning, options);

if flag
    save(['post_samples_', name, '_', num2str(seed), '.mat'], 'post_samples', 'Npix', 'index')
end

end