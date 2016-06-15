function [beta_hat, index] = fit_sim_nonsta_Matern(seed, flag, name)

load(['data_sim_', name, '.mat'])

rng(seed)

% sampling
index = rand_sampler(theta_vec, phi_vec);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);
Y = Y(index);

[x, y, z] = trans_coord(theta_samples, phi_samples);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance function
knots = [0 0 0 0 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_samples);

b_mat(:, 1) = 1;

m = size(b_mat, 2)-1;

beta_init = [zeros(1, m+1) 2 10 0.1];
negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r, b_mat, Y);

lb = [-Inf(1, m+1) 0 0 1e-3];
ub = [Inf(1, m+1) 10 Inf Inf];

[beta_hat, f_min] = nonsta_Matern_fit(negloglik1, beta_init, lb, ub, true);

if flag
    save(['beta_hat_', name, '_', num2str(seed), '.mat'], 'beta_hat', 'index')
end

end