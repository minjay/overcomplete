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
n = 4000;
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
[pot_samples, index] = datasample(resid', n, 'Replace', false,...
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
% the first column is all ones
load('ns.mat')
b_mat_full = kron(b_mat, ones(size(theta, 1), 1));
b_mat = b_mat_full(index, :);
b_mat = [ones(n, 1) b_mat];

r = size(b_mat, 2)-1;

% rescale the observations
Y = pot_samples/1e3;

load('post_samples_exp2.mat')

TT = size(Y, 2);

% init
% c
c_init = zeros(M, TT);
% V
V_inv_init = ones(M, 1); 
% sigma_j_sq
sigma_j_sq_init = [1 beta_hat(end-2:end-1)'];
% eta
eta = beta_hat(1:r+1)';
std_vec = exp(b_mat*eta);
DA = zeros(N, M);
for i = 1:N
    DA(i, :) = std_vec(i)*A(i, :);
end
DATDA = DA'*DA;
DATY = DA'*T;
% tau
tau_init = beta_hat(end);
tau_sq_inv_init = 1/tau_init^2;
% the number of MCMC iterations
T = 5e6;
% the length of the burn-in period
burn_in = 0;
% the length of the thinning interval
thin = 1000;
% the length of the interval to report progress
n_report = 100;

model = struct('A', DA, 'DATDA', DATDA, 'DATY', DATY, 'b_mat', b_mat, 'nu', nu);

data = struct('Y', Y, 'Npix', Npix);

params = struct('c', c_init, 'V', V_inv_init, 'sigma_j_sq', sigma_j_sq_init,...
    'eta', eta_init, 'tau_sigma_sq', tau_sigma_sq, 'tau', tau_sq_inv_init);

options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report, 'save', true);

maxNumCompThreads(32);

post_samples = Gibbs_sampler(model, data, params, options);

save('post_samples_real_no_eta.mat', 'post_samples', 'Npix', 'index', 'theta_samples', 'phi_samples', 'beta_hat')