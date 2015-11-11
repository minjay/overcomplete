parpool(8)

load('data_sim.mat')

rng(1)

% sampling
N = length(Y);
n = 1e3;
index = randsample(N, n);
theta_samples = theta(index);
phi_samples = phi(index);

[x, y, z] = trans_coord(theta_samples, phi_samples);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance funcion
m = 4;
lambda_inv = 2/2.5;
b_mat = get_nonsta_var(m, lambda_inv, theta_samples);

beta_init = [zeros(1, m+1) 0.5 1 0.1];
negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r, b_mat, pot_samples);

lb = [-Inf -Inf -Inf -Inf -Inf 0 0 0];
ub = [Inf Inf Inf Inf Inf 5 Inf Inf];

[beta_hat, f_min] = nonsta_Matern_fit(negloglik1, beta_init, lb, ub, true);

save('beta_hat.mat', 'beta_hat')

delete(gcp)