load('data_EOF_regr_new.mat')
resid = resid_all(1, :);

rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
[pot_samples, index] = datasample(resid', 2000, 'Replace', false,...
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
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples*4, phi_samples);
M = size(A, 2);

% non-stationary variance function
knots = [0 0 0 0 1/3 2/3 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_samples*4);

r = size(b_mat, 2);

% rescale the observations
Y = pot_samples/1e3;

TT = size(Y, 2);

% init
% c
c_init = zeros(M, TT);
% V
V_inv_init = ones(M, 1); 
% sigma_j_sq
sigma_j_sq_init = ones(j_max-j_min, 1);
% eta
eta_init = zeros(r, 1);
% pri_sig of eta
tau_eta_sq = 1e2;
% tau
tau_init = 0.01;
tau_sq_inv_init = 1/tau_init^2;
% tuning parameters
mu_init = zeros(r, 1);
Sigma_init = eye(r);
lambda = 0.001;
% the number of MCMC iterations
T = 6e5;
% the length of the burn-in period
burn_in = 0;
% the length of the thinning interval
thin = 250;
% the length of the interval to report progress
n_report = 100;

model = struct('A', A, 'b_mat', b_mat, 'nu', nu);

data = struct('Y', Y, 'Npix', Npix);

params = struct('c', c_init, 'V', V_inv_init, 'sigma_j_sq', sigma_j_sq_init,...
    'eta', eta_init, 'tau_eta_sq', tau_eta_sq, 'tau', tau_sq_inv_init);

tuning = struct('mu', mu_init, 'Sigma', Sigma_init, 'lambda', lambda);

options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report);

post_samples = Gibbs_sampler_AM4(model, data, params, tuning, options);

save('post_samples_real.mat', 'post_samples', 'Npix', 'index', 'theta_samples', 'phi_samples')