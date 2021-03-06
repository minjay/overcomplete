% batch file for sim_AM_sigma_rep.m

addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

parpool(10)

eta_seed = 3;
R = 10;
n_seed = 10;

eta_est_cell = cell(1, n_seed);
sigma_j_est_cell = cell(1, n_seed);
tau_est_cell = cell(1, n_seed);

tau = 0.1;

B = 2;
j_min = 2;
j_max = 3;
alpha = 3;

sigma_j = B.^(-alpha/2*(j_min:j_max))*10;

rng(eta_seed)
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
r = 4;
eta = [0; randn(r, 1)];

eta

parfor seed = 1:n_seed
    maxNumCompThreads(4);
    [eta_est, sigma_j_est, tau_est] = sim_AM_sigma_rep(R, seed, eta, sigma_j, tau);
    eta_est_cell{seed} = eta_est;
    sigma_j_est_cell{seed} = sigma_j_est;
    tau_est_cell{seed} = tau_est;
end

eta_est_all = cell2mat(eta_est_cell);
sigma_j_est_all = cell2mat(sigma_j_est_cell);
tau_est_all = cell2mat(tau_est_cell);

filename = ['sim_rep', num2str(eta_seed), '.mat'];
save(filename, 'eta_est_all', 'sigma_j_est_all', 'tau_est_all', 'eta', 'sigma_j', 'tau')
