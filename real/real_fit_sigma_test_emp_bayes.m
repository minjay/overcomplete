load('data_EOF_regr.mat')
resid = resid_all(1, :);

resid_norm = resid/1e3;
cf = reshape(resid_norm, size(phi));
vmin = min(resid_norm);
vmax = max(resid_norm);
vmag = vmin:0.1:vmax;
phi_rot = phi+pi/2;
[x, y] = pol2cart(phi_rot, theta/pi*180);
h = mypolar([0 2*pi], [0 max(theta(:))/pi*180], x, y, cf, vmag);
delete(h)
shading flat
title('Electric Potential','FontName','times','Fontsize',10)
xlabel(sprintf('Min %6.1f  Max %5.1f [kV]',vmin,vmax),'FontName','times','Fontsize',10)
% plot_pot(reshape(resid, size(phi)), phi, theta, 1000, max(abs(resid)));


rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
index = rand_sampler_real(theta_vec*4);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);
pot_samples = resid(index)';

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

% non-stationary variance function
r = 4;
lambda_inv = 2/2.5;
b_mat = get_nonsta_var(r, lambda_inv, theta_samples*4);

% rescale the observations
Y = pot_samples/1e3;

% non-stationary variance function
r = 4;
lambda_inv = 2/2.5;
XX = get_nonsta_var(r, lambda_inv, theta_vec*4);

% get emp prior
emp_std_vec = std(resid_all, 0, 1);
yy = log(emp_std_vec');
eta_est_regr = (XX'*XX)\(XX'*yy);
eta_hat = eta_est_regr(2:r+1);

% init
% c
c_init = zeros(M, 1);
% V
V_inv_init = ones(M, 1); 
% sigma_j_sq
sigma_j_sq_init = ones(j_max-j_min, 1);
% eta
eta_init = zeros(r+1, 1);
eta_init(2:r+1) = eta_hat;
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
lambda = 0.001;
% the number of MCMC iterations
T = 5e5;
% the length of the burn-in period
burn_in = 25*1e4;
% the length of the thinning interval
thin = 250;
% the length of the interval to report progress
n_report = 100;

model = struct('A', A, 'b_mat', b_mat, 'nu', nu);

data = struct('Y', Y, 'Npix', Npix);

params = struct('c', c_init, 'V', V_inv_init, 'sigma_j_sq', sigma_j_sq_init,...
    'eta', eta_init, 'tau_sigma_sq', tau_sigma_sq, 'tau_eta_sq', tau_eta_sq,...
    'tau', tau_sq_inv_init, 'eta_hat', eta_hat);

tuning = struct('mu', mu_init, 'Sigma', Sigma_init, 'lambda', lambda);

options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report);

post_samples = Gibbs_sampler_AM3(model, data, params, tuning, options);

save('post_samples_real.mat', 'post_samples', 'Npix', 'index', 'theta_samples', 'phi_samples')