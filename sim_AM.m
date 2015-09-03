clear

nu = 3;

alpha = 4;
tau = 0.1;

% the grid
B = 2;
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
n = length(tp);
theta = zeros(n, 1);
phi = zeros(n, 1);
for i = 1:n
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end

rng(1)
theta = theta+randn(n, 1)*pi/10;
theta(theta<0) = theta(theta<0)+pi;
theta(theta>pi) = theta(theta>pi)-pi;
phi = phi+randn(n, 1)*2*pi/10;
phi(phi<0) = phi(phi<0)+2*pi;
phi(phi>2*pi) = phi(phi>2*pi)-2*pi;

j_min = 2;
j_max = 4;
j_len = j_max-j_min+1;
n_dist = 1e3;

[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi, n_dist);
[N, M] = size(A);

sigma_j = B.^(-alpha/2*(j_min:j_max));
fj_sq = zeros(M, 1);
c = zeros(M, 1);
for j = j_min:j_max
    index = j-j_min+1;
    range = sum(Npix(1:index))-Npix(index)+1:sum(Npix(1:index));
    fj_sq(range) = B^(-alpha*j)*ones(Npix(index), 1);
    c(range) = sigma_j(index)*trnd(nu, Npix(index), 1);
end

r = 4;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(r+1, N);
b_mat(1, :) = 1;
for i = 2:r+1
    b_mat(i, :) = normpdf(theta, mu(i-1), lambda);
end

eta = [1.5; randn(r, 1)];

Y = diag(exp(b_mat'*eta))*A*c+randn(N, 1)*tau;

% init
% V
V_inv_init = ones(M, 1); 
% eta
eta_init = zeros(r+1, 1);
% pri_sig of eta_0
tau_sigma_sq = 1;
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
T = 150000;
% the length of the burn-in period
burn_in = 50000;
% the length of the thinning interval
thin = 100;
% the length of the interval to report progress
n_report = 100;

model = struct('A', A, 'fj_sq', fj_sq, 'b_mat', b_mat, 'nu', nu);

data = Y;

params = struct('V', V_inv_init, 'eta', {eta_init, tau_sigma_sq, tau_eta_sq},...
    'tau', tau_sq_inv_init);

tuning = struct('mu', mu_init, 'Sigma', Sigma_init, 'lambda', lambda);

options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report);

post_samples = Gibbs_sampler_AM(model, data, params, tuning, options);
