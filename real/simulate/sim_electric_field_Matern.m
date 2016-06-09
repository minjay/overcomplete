load('data_EOF_regr_new.mat')
load('beta_hat.mat')

% set seed
rng(1)

theta_vec = theta(:);
phi_vec = phi(:);
% stretching
[x, y, z] = trans_coord(theta_vec*4, phi_vec);

N = length(x);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec*4);

b_mat(:, 1) = 1;

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
% measurement error is not included given that it is not differentiable
cov_mat_Matern = get_cov_nonsta_Matern(beta, r, b_mat);

% figure
T = 5e3;
Y_sim_Matern = mvnrnd(zeros(T, N), cov_mat_Matern)*1000;

HX = theta(1, :);
HY = phi(:, 1);
R = 6.5*1e6;

E_theta_Matern = zeros(N, T);
E_phi_Matern = zeros(N, T);

for t = 1:T
    [FX, FY] = gradient(reshape(Y_sim_Matern(t, :), size(phi)), HX, HY);
    tmp = -FX/R;
    E_theta_Matern(:, t) = tmp(:);
    tmp = -FY./(R*sin(theta));
    E_phi_Matern(:, t) = tmp(:);
end

save('sim_energy_Matern.mat', 'E_theta_Matern', 'E_phi_Matern')