load('data.mat')

rng(1)

% sampling
N = 1e3;
[pot_samples, theta_samples, phi_samples, index] = sampling_data(resid,...
    theta, phi, N, 0);

index = find(theta_samples<=25/180*pi);
N = length(index);
pot_samples = pot_samples(index);
theta_samples = theta_samples(index);
phi_samples = phi_samples(index);

% fit
% parameter specification
B = 2;
j_min = 2;
j_max = 2;
nu = 4;
alpha = 4;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples, phi_samples);
M = size(A, 2);

% non-stationary variance funcion
r = 3;
theta_range = 100/180*pi;
mu = theta_range/(r+1)*(1:r);
lambda = theta_range/(r+1)*2.5/2;
b_mat = zeros(N, r+1);
b_mat(:, 1) = 1;
for i = 2:r+1
    b_mat(:, i) = exp(-(theta_samples*4-mu(i-1)).^2/2/lambda^2);
end

fj_sq = zeros(M, 1);
st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    fj_sq(range) = B^(-alpha*j)*ones(Npix(index_j), 1);
    st = st+Npix(index_j);
end

% rescale the observations
Y = pot_samples'/1e3;

% init
% c
c_init = zeros(M, 1);
% V
V_inv_init = ones(M, 1); 
% eta
eta_init = zeros(r+1, 1);
% pri_sig of eta_0
tau_sigma_sq = 1e4;
% pri_sig of eta
tau_eta_sq = 0.25^2;
% tau
tau_init = 0.01;
tau_sq_inv_init = 1/tau_init^2;
% tuning parameters
mu_init = zeros(r+1, 1);
Sigma_init = eye(r+1);
lambda = 0.01;
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

save('post_samples_real.mat', 'post_samples', 'Npix')