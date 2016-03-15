load('data_EOF_regr.mat')
resid = resid_all(1, :);
emp_std_vec = std(resid_all, 0, 1);
resid = resid./emp_std_vec;

rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
index = rand_sampler_real(theta_vec*4);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);
pot_samples = resid(:, index)';

% plot
% plot_samples(theta_vec, index, phi_samples, pot_samples)

% fit
% parameter specification
B = 2;
j_min = 2;
j_max = 4;
nu = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples*4, phi_samples);
M = size(A, 2);

% rescale the observations
Y = pot_samples/1e3;

TT = size(Y, 2);

% init
% c
c_init = zeros(M, TT);
% V
V_inv_init = ones(M, 1); 
% sigma_j_sq
sigma_j_sq_init = ones(j_max-j_min+1, 1);
% tau
tau_init = 0.01;
tau_sq_inv_init = 1/tau_init^2;
% the number of MCMC iterations
T = 1e5;
% the length of the burn-in period
burn_in = 0;
% the length of the thinning interval
thin = 100;
% the length of the interval to report progress
n_report = 100;

model = struct('A', A, 'nu', nu);

data = struct('Y', Y, 'Npix', Npix);

params = struct('c', c_init, 'V', V_inv_init, 'sigma_j_sq', sigma_j_sq_init,...
    'tau', tau_sq_inv_init);

options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report);

post_samples = Gibbs_sampler_AM10(model, data, params, options);

save('post_samples_real.mat', 'post_samples', 'Npix', 'index', 'theta_samples', 'phi_samples')