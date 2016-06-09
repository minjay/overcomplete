clear

% load data
load('data_EOF_regr_new.mat')
load('post_samples_real_new_resid_j24_pen_B_spline_cubic_2knots_aurora_n4000_sigma_j_sq_0.01_0.0002_nu4_long_init0.mat')

% get theta and phi
theta_vec = theta(1, :)';
phi_vec = phi(1, :)';

N = length(theta_vec);

post_samples_eta = post_samples.eta(:, 1501:end);

eta_est = mean(post_samples_eta, 2);

% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec*4);

b_mat(:, 1) = 1;

std_vec_est = exp(b_mat*eta_est);

std_vec_post = exp(b_mat*post_samples_eta);

nu = 4;
sigma_j = sqrt([1 0.01 0.0002]);
variance = plot_corr_fun_flex(2, sigma_j, 2, 4, nu, 1000, 'r');
std_vec_est = std_vec_est*sqrt(variance);
std_vec_post = std_vec_post*sqrt(variance);

CI = quantile(std_vec_post', [0.025 0.975]);

% plot empirical std function
figure
emp_std_vec = std(resid_all, 0, 1);
h1 = plot(theta(:)/pi*180, emp_std_vec/1e3, '.');

hold on

%%% plot fitted std function by nonGau-need
h2 = plot(theta_vec/pi*180, std_vec_est, 'b-', 'LineWidth', 2);
hold on
h3 = plot(theta_vec/pi*180, CI(1, :), 'r--', 'LineWidth', 2);
plot(theta_vec/pi*180, CI(2, :), 'r--', 'LineWidth', 2)

%%% plot fitted std function by Gau-need
load('beta_hat_needlet.mat')
beta = beta_hat(1:end-1);
B = 2;
j_min = 2;
j_max = 4;
[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_vec*4, phi_vec);
cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A);
h4 = plot(theta_vec/pi*180, sqrt(diag(cov_mat)), 'g', 'LineWidth', 2);

%%% plot fitted std function by Matern
load('beta_hat.mat')
beta = beta_hat(1:end-1);
% stretching
[x, y, z] = trans_coord(theta_vec*4, phi_vec);
% get distance matrix
r = get_chordal_dist(x, y, z);
cov_mat = get_cov_nonsta_Matern(beta, r, b_mat);
h5 = plot(theta_vec/pi*180, sqrt(diag(cov_mat)), 'c', 'LineWidth', 2);

lg = legend([h1 h2 h3 h4 h5], {'Empirical std function', 'Est. std function (nonGau-need)',...
    '95% CI endpoints', 'Est. std function (Gau-need)', 'Est. std function (Gau-Matern)'});
set(lg, 'FontSize', 12);
xlabel('Co-latitude (\circ)')
ylabel('Standard deviation (kV)')
set(gca, 'FontSize', 12);
axis tight
