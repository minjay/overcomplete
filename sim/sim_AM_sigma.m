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

% perturbation
theta = theta+randn(N, 1)*pi/10;
theta(theta<0) = theta(theta<0)+pi;
theta(theta>pi) = theta(theta>pi)-pi;
phi = phi+randn(N, 1)*2*pi/10;
phi(phi<0) = phi(phi<0)+2*pi;
phi(phi>2*pi) = phi(phi>2*pi)-2*pi;

j_min = 2;
j_max = 4;

% design matrix A
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta, phi);
M = size(A, 2); 

sigma_j = B.^(-alpha/2*(j_min:j_max));

% non-stationary variance funcion
r = 4;
mu = pi/(r+1)*(1:r);
lambda = pi/(r+1)*2.5/2;
b_mat = zeros(N, r+1);
b_mat(:, 1) = 1;
for i = 2:r+1
    b_mat(:, i) = exp(-(theta-mu(i-1)).^2/2/lambda^2);
end

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
sigma_j_sq_init = ones(j_max-j_min, 1);
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

post_samples = Gibbs_sampler_AM2(model, data, params, tuning, options);

save('post_samples.mat', 'post_samples', 'c', 'eta', 'theta', 'phi', 'sigma_j')