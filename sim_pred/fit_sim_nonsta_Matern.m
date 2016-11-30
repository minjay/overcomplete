function [beta_hat, index, index_region] = fit_sim_nonsta_Matern(seed, flag, name, width)

load(['data_sim_', name, '.mat'])

rng(seed)

% sampling
[index, index_region] = rand_sampler(phi, width);
theta_samples = theta(index);
phi_samples = phi(index);
Y = Y(index);

[x, y, z] = trans_coord(theta_samples, phi_samples);

% get distance matrix
r_dist = get_chordal_dist(x, y, z);

% non-stationary variance function
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_samples);

b_mat(:, 1) = 1;

r = size(b_mat, 2)-1;

beta_init = [zeros(1, r+1) 2 10 1e-2];
negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r_dist, b_mat, Y);

lb = [-10*ones(1, r+1) 0 0 1e-3];
ub = [10*ones(1, r+1) 10 Inf Inf];

[beta_hat, f_min] = nonsta_Matern_fit(negloglik1, beta_init, lb, ub, true);

if flag
    save(['beta_hat_', name, '_', num2str(seed), '.mat'], 'beta_hat', 'index', 'index_region')
end

end