clear
rng(1)

nu = 4;
alpha = 3;

tau = 0.1;

% the grid
B = 2;
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta = zeros(N, 1);
phi = zeros(N, 1);
for i = 1:N
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end

j_min = 2;
j_max = 3;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);
M = size(A, 2); 

sigma_j = B.^(-alpha/2*(j_min:j_max));
sigma_j = sigma_j/sigma_j(1);

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta);

b_mat(:, 1) = 1;

r = size(b_mat, 2)-1;

rng(2)
eta = [1.5; randn(r, 1)];
std_vec = exp(b_mat*eta);
DA = zeros(N, M);
for i = 1:N
    DA(i, :) = std_vec(i)*A(i, :);
end

c = zeros(M, 1);
st = 1;
for j = j_min:j_max
    index_j = j-j_min+1;
    range = st:st+Npix(index_j)-1;
    c(range) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
    st = st+Npix(index_j);
end

Y = DA*c+randn(N, 1)*tau;

% init
% c
c_init = zeros(M, 1);
% V
V_inv_init = ones(M, 1); 
% sigma_j_sq
sigma_j_sq_init = 0.1^2;
% eta
eta_init = zeros(r+1, 1);
% pri_sig of eta_0
tau_sigma_sq = 1e2;
% pri_sig of eta
tau_eta_sq = 0.25^2;
% tau
tau_init = 0.01;
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

save('post_samples.mat', 'post_samples', 'sigma_j', 'eta', 'c', 'tau', 'theta', 'phi')